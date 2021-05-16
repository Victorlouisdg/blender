#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"
#include "BLI_float3x3.hh"
#include "BLI_sparse_matrix.hh"

#include "BLI_index_range.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"

#include <tuple>

// #include "bmesh.h"

using blender::Array;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using blender::Span;
using blender::SparseMatrix;

/* Utilities */

struct int3 {
  int i, j, k;

  int3() = default;

  int3(const int *ptr) : i{ptr[0]}, j{ptr[1]}, k{ptr[2]}
  {
  }

  int3(int i, int j, int k) : i{i}, j{j}, k{k}
  {
  }

  operator int *()
  {
    return &i;
  }
};

struct int4 {
  int i, j, k, l;

  int4() = default;

  int4(const int *ptr) : i{ptr[0]}, j{ptr[1]}, k{ptr[2]}, l{ptr[3]}
  {
  }

  int4(int i, int j, int k, int l) : i{i}, j{j}, k{k}, l{l}
  {
  }

  operator int *()
  {
    return &i;
  }
};

/* Copied from point_distirbute. */
static Span<MLoopTri> get_mesh_looptris(const Mesh &mesh)
{
  /* This only updates a cache and can be considered to be logically const. */
  const MLoopTri *looptris = BKE_mesh_runtime_looptri_ensure(const_cast<Mesh *>(&mesh));
  const int looptris_len = BKE_mesh_runtime_looptri_len(&mesh);
  return {looptris, looptris_len};
}

/* Cloth Simulator based off of the popular paper "Large steps in cloth simulation."
 * by Baraff and Witkin, hence the name. In comments below Ithis paper is referred to
 * as BW98.
 *
 * Each timestep, the simulator builds a linear system of equations that's filled with
 * several conditions that encode how we want the cloth to behave. (e.g. how it responds to
 * stretch, or to collisions.) sidenote: collisions response can be done in several ways.
 */
class ClothSimulatorBaraffWitkin {
 public:
  int n_vertices;
  int n_triangles;

  /* State */
  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  /* Indices */
  Array<int3> triangle_vertex_indices;
  Array<int4> bending_vertex_indices;

  /* Parameters / Configuration */

  /* Note about the "uv" in the names below: these letters are what's used in the original paper.
   * These uv coordinates is used to represent the 2D rest state of the cloth, and should thus not
   * necessarily be shared with the uv coordinates for texturing.
   *
   * A more descriptive name for vertex_positions_uv would have been vertex_rest_positions_2D.
   */
  Array<float2>
      vertex_positions_uv;  // This is only needed for initialization so maybe should be removed.
  Array<float> vertex_masses;
  Array<float> triangle_stretch_stiffness_u;
  Array<float> triangle_stretch_stiffness_v;
  Array<float> triangle_areas_uv;
  Array<float3> edges_normals;

  float density = 0.2f;            /* in kg/m^2*/
  float standard_gravity = -9.81f; /* in m/s^2*/

  /* Datastructures for intermediate computations. */
  SparseMatrix vertex_force_derivatives;  // <float3x3> elements
  Array<float3> vertex_forces;
  Array<float3> triangle_normals;
  /* The inverted mass matrix is stored as an Array because it's a diagonal matrix.
   * It's components are float3 instead of float so we can use the "mass modification"
   * technique on it's components to make vertices immovable in certain directions.
   */
  Array<float3> vertex_mass_matrix_inverted;

  /* Precomputed quantities. */
  Array<float2x2> triangle_inverted_delta_u_matrices;
  Array<float3> triangle_wu_derivatives;
  Array<float3> triangle_wv_derivatives;

  /* Currently the simulator has its own local copy of all the necessary
   * mesh attributes e.g. vertex positions. This was easier to implement and I believe could make
   * the simulator faster due to better data locality. I'll make setters() for certain attributes
   * later that do the required recalculations.
   */
  void initialize(const Mesh &mesh)
  {
    std::cout << "Cloth simulation initialisation" << std::endl;

    n_vertices = mesh.totvert;

    vertex_force_derivatives = SparseMatrix(n_vertices);

    initialize_vertex_attributes(mesh);
    initialize_triangle_attributes(mesh);
  };

  void initialize_vertex_attributes(const Mesh &mesh)
  {
    vertex_positions = Array<float3>(n_vertices);
    vertex_velocities = Array<float3>(n_vertices, float3(0.0f));
    vertex_forces = Array<float3>(n_vertices, float3(0.0f));
    vertex_positions_uv = Array<float2>(n_vertices);

    for (const int i : IndexRange(mesh.totvert)) {
      MVert vertex = mesh.mvert[i];
      vertex_positions[i] = vertex.co;
    }

    /* TODO: rethink this UV position initialization in the context of seams/multiple cloth pieces.
     */
    for (const int i : IndexRange(mesh.totloop)) {
      const MLoop &loop = mesh.mloop[i];
      const float2 uv = mesh.mloopuv[i].uv;
      vertex_positions_uv[loop.v] = uv;
    }
  }

  void initialize_triangle_attributes(const Mesh &mesh)
  {
    /* While technically a vertex attribute, the vertex_masses get initialized here because
     * vertex mass is set based on the area of it's adjacent triangles.
     */
    vertex_masses = Array<float>(n_vertices, 0.0f);

    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    n_triangles = looptris.size();

    triangle_vertex_indices = Array<int3>(n_triangles);
    triangle_areas_uv = Array<float>(n_triangles);
    triangle_normals = Array<float3>(n_triangles);
    triangle_inverted_delta_u_matrices = Array<float2x2>(n_triangles);
    triangle_wu_derivatives = Array<float3>(n_triangles);
    triangle_wv_derivatives = Array<float3>(n_triangles);

    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop = looptri.tri[0];
      const int v1_loop = looptri.tri[1];
      const int v2_loop = looptri.tri[2];
      const int v0_index = mesh.mloop[v0_loop].v;
      const int v1_index = mesh.mloop[v1_loop].v;
      const int v2_index = mesh.mloop[v2_loop].v;
      const float3 v0_pos = vertex_positions[v0_index];
      const float3 v1_pos = vertex_positions[v1_index];
      const float3 v2_pos = vertex_positions[v2_index];

      triangle_vertex_indices[looptri_index] = {v0_index, v1_index, v2_index};

      const float3 normal = float3::cross(v1_pos - v0_pos, v2_pos - v0_pos);
      /* Storing this normal is not that useful because it will probablty need to be recalculated.
       */
      triangle_normals[looptri_index] = normal;

      /* The length of the cross product vector is equal to the area of the parallelogram spanned
       * by the two input vectors, the triangle area is conveniently half of this area. */
      const float triangle_area_uv = normal.length() / 2;
      triangle_areas_uv[looptri_index] = triangle_area_uv;

      /* This particular way of initializing vertex mass is described in BW98 section 2.2.
       * Different methods could definitely be tried, when Simulation nodes comes, vertex_mass
       * could even just be an attribute provided by the user.  */
      float vertex_mass = density * triangle_area_uv / 3.0f;
      vertex_masses[v0_index] += vertex_mass;
      vertex_masses[v1_index] += vertex_mass;
      vertex_masses[v2_index] += vertex_mass;

      float2 uv0 = vertex_positions_uv[v0_index];
      float2 uv1 = vertex_positions_uv[v1_index];
      float2 uv2 = vertex_positions_uv[v2_index];

      float u0 = uv0[0];
      float u1 = uv1[0];
      float u2 = uv2[0];
      float v0 = uv0[1];
      float v1 = uv1[1];
      float v2 = uv2[1];

      float delta_u1 = u1 - u0;
      float delta_u2 = u2 - u0;
      float delta_v1 = v1 - v0;
      float delta_v2 = v2 - v0;

      /* If anyone knows how to do this in one line, let me know. */
      float array[2][2] = {{delta_u1, delta_u2}, {delta_v1, delta_v2}};
      float2x2 delta_u = float2x2(array);
      triangle_inverted_delta_u_matrices[looptri_index] = delta_u.inverted();

      float dw_denominator = 1.0f / (delta_u1 * delta_v2 - delta_u2 * delta_v1);

      float3 dwu_dx;
      dwu_dx[0] = (delta_v1 - delta_v2) * dw_denominator;
      dwu_dx[1] = delta_v2 * dw_denominator;
      dwu_dx[2] = -delta_v1 * dw_denominator;

      float3 dwv_dx;
      dwv_dx[0] = (delta_u2 - delta_u1) * dw_denominator;
      dwv_dx[1] = -delta_u2 * dw_denominator;
      dwv_dx[2] = delta_u1 * dw_denominator;

      triangle_wu_derivatives[looptri_index] = dwu_dx;
      triangle_wv_derivatives[looptri_index] = dwv_dx;
    }
  };

  void step()
  {
    std::cout << "Cloth simulation step" << std::endl;
    reset_forces_and_derivatives();
    calculate_forces_and_derivatives();

    for (int i : IndexRange(n_vertices)) {
      vertex_mass_matrix_inverted[i] = float3(vertex_masses[i]);
      /* TODO later set mass of pinned vertices to zero etc. */
    }
  };

  void reset_forces_and_derivatives()
  {
    vertex_forces.fill(float3(0.0f));
  }

  void calculate_forces_and_derivatives()
  {
    for (int vertex_index : IndexRange(n_vertices)) {
      calculate_gravity(vertex_index);
    }

    for (int triangle_index : IndexRange(n_triangles)) {
      auto [wu, wv] = calculate_w_uv(triangle_index);
      calculate_stretch(triangle_index, wu, wv);
    }
  }

  std::tuple<float3, float3> calculate_w_uv(int triangle_index)
  {
    auto [v0_index, v1_index, v2_index] = triangle_vertex_indices[triangle_index];
    float3 x0 = vertex_positions[v0_index];
    float3 x1 = vertex_positions[v1_index];
    float3 x2 = vertex_positions[v2_index];

    float3 delta_x1 = x1 - x0;
    float3 delta_x2 = x2 - x0;

    float2x2 delta_u = triangle_inverted_delta_u_matrices[triangle_index];

    /* This is basically a hard coded 3x2 @ 2x2 matrix multiplication.
     * I didn't want to create a float3x2 class specifically for this.
     */
    float wu0 = delta_x1[0] * delta_u.values[0][0] + delta_x2[0] * delta_u.values[1][0];
    float wu1 = delta_x1[1] * delta_u.values[0][0] + delta_x2[1] * delta_u.values[1][0];
    float wu2 = delta_x1[2] * delta_u.values[0][0] + delta_x2[2] * delta_u.values[1][0];
    float wv0 = delta_x1[0] * delta_u.values[0][1] + delta_x2[0] * delta_u.values[1][1];
    float wv1 = delta_x1[1] * delta_u.values[0][1] + delta_x2[1] * delta_u.values[1][1];
    float wv2 = delta_x1[2] * delta_u.values[0][1] + delta_x2[2] * delta_u.values[1][1];

    float3 wu = float3(wu0, wu1, wu2);
    float3 wv = float3(wv0, wv1, wv2);

    return {wu, wv};
  }

  void calculate_gravity(int vertex_index)
  {
    float3 gravity = float3(0.0f);
    gravity.z = vertex_masses[vertex_index] * standard_gravity;  // F_grav = m * g
    vertex_forces[vertex_index] += gravity;
  }

  void calculate_stretch(int triangle_index, float3 wu, float3 wv)
  {
    float area_uv = triangle_areas_uv[triangle_index];

    float wu_norm = wu.normalize_and_get_length();
    float wv_norm = wv.normalize_and_get_length();

    float Cu = area_uv * (wu_norm - 1.0);
    float Cv = area_uv * (wv_norm - 1.0);

    int3 vertex_indices = triangle_vertex_indices[triangle_index];
    float3 dwu_dx = triangle_wu_derivatives[triangle_index];
    float3 dwv_dx = triangle_wv_derivatives[triangle_index];

    for (int m : IndexRange(3)) {
      int vertex_index = vertex_indices[m];
    }
  }

  ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin constructed" << std::endl;
  }

  ~ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin destructed" << std::endl;
  }
};