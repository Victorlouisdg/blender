#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"
#include "BLI_float3x3.hh"
#include "BLI_int3.hh"
#include "BLI_int4.hh"
#include "BLI_sparse_matrix.hh"

#include "BLI_index_range.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"

#include <tuple>

/* Note about variable names in the simulator: in general I've tried to give the member
 * variables/attributes long and clear names. I deliberately read these out into local variable
 * with much shorter (often 2 letter) names. This is because long variable names would obfuscate
 * all the structure of the math operations.
 *
 * For the loop index variables in this file I've mantained this convention:
 * i, j: index the Arrays of all vertices/triangles
 * m, n: index the vertices of a single triangle (or 2 triangles for the bending condition)
 * s, t: index the components (x, y, z) of a vector.
 */

using blender::Array;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using blender::int3;
using blender::int4;
using blender::solve_gauss_seidel;
using blender::solve_pcg_filtered;
using blender::Span;
using blender::SparseMatrix;

/* Utilities */

/* Copied from point_distirbute. */
static Span<MLoopTri> get_mesh_looptris(const Mesh &mesh)
{
  /* This only updates a cache and can be considered to be logically const. */
  const MLoopTri *looptris = BKE_mesh_runtime_looptri_ensure(const_cast<Mesh *>(&mesh));
  const int looptris_len = BKE_mesh_runtime_looptri_len(&mesh);
  return {looptris, looptris_len};
}

/* Cloth Simulator based off of the popular paper "Large steps in cloth simulation."
 * by Baraff and Witkin, hence the name. In comments below this paper is referred to
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

  int n_substeps;
  float step_time; /* Currently assuming 30 fps. */
  float substep_time;

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
  Array<float2> vertex_positions_uv;  // Only needed for initialization so maybe should be removed.
  Array<float> vertex_masses;
  Array<float> triangle_stretch_stiffness_u;
  Array<float> triangle_stretch_stiffness_v;
  Array<float> triangle_areas_uv;
  Array<float3> edges_normals;

  float density = 0.2f;            /* in kg/m^2 */
  float standard_gravity = -9.81f; /* in m/s^2 */

  /* Datastructures for intermediate computations. */
  SparseMatrix vertex_force_derivatives;  // <float3x3> elements
  SparseMatrix *vertex_force_derivatives_pointer;

  Array<float3> vertex_forces;
  Array<float3> triangle_normals;
  /* The inverted mass matrix is stored as an Array because it's a diagonal matrix.
   * It's components are float3 instead of float so we can use the "mass modification"
   * technique on it's components to make vertices immovable in certain directions.
   */
  Array<float3x3> vertex_mass_matrix_inverted;

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

    n_substeps = 10;
    step_time = 1.0f / 30.0f; /* Currently assuming 30 fps. */
    substep_time = step_time / n_substeps;

    n_vertices = mesh.totvert;

    initialize_vertex_attributes(mesh);
    initialize_triangle_attributes(mesh);
    vertex_force_derivatives_pointer = new SparseMatrix(n_vertices);
    vertex_force_derivatives = *vertex_force_derivatives_pointer;

    // verify_w_derivatives(); /* For debugging, should become a test. */
  };

  void step()
  {
    for (int substep : IndexRange(n_substeps)) {
      UNUSED_VARS(substep);
      reset_forces_and_derivatives();
      calculate_forces_and_derivatives();
      fill_inverted_mass_matrix();

      // integrate_explicit_forward_euler();
      // integrate_implicit_backward_euler_pcg();
      integrate_implicit_backward_euler_pcg_filtered();
    }
  };

  void verify_w_derivatives()
  {
    /* This method does a finite difference verification of the w_derivatives.
        Recall that the derivative of a function is defined as: [f(x + h) - f(x)] / h
    */
    for (int i : IndexRange(n_triangles)) {

      std::cout << std::endl << "Triangle: " << i << std::endl;
      std::cout << "Analytic wu_derivative" << triangle_wu_derivatives[i] << std::endl;

      float h = 0.00001f;

      int3 vertex_indices = triangle_vertex_indices[i];

      for (int m : IndexRange(3)) {    // vertex0, vertex1, vertex2
        for (int s : IndexRange(3)) {  // x, y, z
          auto [wu, wv] = calculate_w_uv(i);
          vertex_positions[vertex_indices[m]][s] += h;
          auto [wu_h, wv_h] = calculate_w_uv(i);
          std::cout << (wu_h - wu) / h << std::endl;
        }
      }
    }
  };

  void initialize_vertex_attributes(const Mesh &mesh)
  {
    vertex_positions = Array<float3>(n_vertices);
    vertex_velocities = Array<float3>(n_vertices, float3(0.0f));
    vertex_forces = Array<float3>(n_vertices, float3(0.0f));
    vertex_positions_uv = Array<float2>(n_vertices);

    /* vertex_mass_matrix_inverted is only allocated here, it's values are refilled each timestep.
     */
    vertex_mass_matrix_inverted = Array<float3x3>(n_vertices);

    for (const int i : IndexRange(n_vertices)) {
      MVert vertex = mesh.mvert[i];
      vertex_positions[i] = vertex.co;
    }
  }

  void initialize_triangle_attributes(const Mesh &mesh)
  {
    /* TODO set vertex mass based on Voronoi area of vertex. I read this in a paper somewhere but
     * can't find it. */
    float fixed_vertex_mass = 1.0f / n_vertices;
    vertex_masses = Array<float>(n_vertices, fixed_vertex_mass);

    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    n_triangles = looptris.size();

    triangle_vertex_indices = Array<int3>(n_triangles);
    triangle_areas_uv = Array<float>(n_triangles);
    triangle_normals = Array<float3>(n_triangles);
    triangle_inverted_delta_u_matrices = Array<float2x2>(n_triangles);
    triangle_wu_derivatives = Array<float3>(n_triangles);
    triangle_wv_derivatives = Array<float3>(n_triangles);

    float stretch_stiffness = 10000.0f;
    triangle_stretch_stiffness_u = Array<float>(n_triangles, stretch_stiffness);
    triangle_stretch_stiffness_v = Array<float>(n_triangles, stretch_stiffness);

    /* Think about disjoint pieces of cloth and seams etc. */
    for (const int i : IndexRange(mesh.totloop)) {
      const MLoop &loop = mesh.mloop[i];
      const float2 uv = mesh.mloopuv[i].uv;
      vertex_positions_uv[loop.v] = uv;  // TODO: fix UVs being normalized.
    }

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
       * could even just be an attribute provided by the user.
       *
       * I've commented this method out for now become it gave weird mass distributions. E.g. think
       * of a plane consisting of 2 triangles, the vertices on the shared edge would have double
       * the mass of the non-shared vertices.
       */
      // float vertex_mass = density * triangle_area_uv / 3.0f;
      // vertex_masses[v0_index] += vertex_mass;
      // vertex_masses[v1_index] += vertex_mass;
      // vertex_masses[v2_index] += vertex_mass;

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

  void reset_forces_and_derivatives()
  {
    vertex_forces.fill(float3(0.0f));
    vertex_force_derivatives.clear();
  }

  void calculate_forces_and_derivatives()
  {
    for (int vertex_index : IndexRange(n_vertices)) {
      calculate_gravity(vertex_index);
    }

    for (int triangle_index : IndexRange(n_triangles)) {

      // int i = triangle_index;
      // int3 vertex_indices = triangle_vertex_indices[i];

      auto [wu, wv] = calculate_w_uv(triangle_index);
      calculate_stretch(triangle_index, wu, wv);

      // std::cout << vertex_force_derivatives.get(vertex_indices[0], vertex_indices[1]]);
    }
  }

  void fill_inverted_mass_matrix()
  {
    for (int i : IndexRange(n_vertices)) {
      float3x3 matrix = float3x3();
      float inverted_mass = 1.0f / vertex_masses[i];
      matrix.values[0][0] = inverted_mass;
      matrix.values[1][1] = inverted_mass;
      matrix.values[2][2] = inverted_mass;
      vertex_mass_matrix_inverted[i] = matrix;

      /* Temporary hardcoded pin of vertices 0 and 1. Later pinning should be governed by an
       * attribute. */
      if (i == 0 || i == 1) {
        vertex_mass_matrix_inverted[i] = float3x3(0.0f);
      }
    }
  }

  void integrate_explicit_forward_euler()
  {
    /* vertex_accelerations are only needed for debugging with explicit integration. */
    Array<float3> vertex_accelerations = Array<float3>(n_vertices);

    /* TODO: look into doing these Array-based operations with std::transform() etc. */
    for (int i : IndexRange(n_vertices)) {
      vertex_accelerations[i] = vertex_mass_matrix_inverted[i] * vertex_forces[i];
      vertex_velocities[i] += vertex_accelerations[i] * substep_time;
      vertex_positions[i] += vertex_velocities[i] * substep_time;
    }
  }

  void integrate_implicit_backward_euler_pcg()
  {

    /* Assemble linear system */
    float h = substep_time;
    SparseMatrix &dfdx = vertex_force_derivatives;
    Array<float3> &v0 = vertex_velocities;
    Array<float3> &f0 = vertex_forces;
    Array<float3x3> &M_inv = vertex_mass_matrix_inverted;

    std::cout << "dfdx symmetric: " << dfdx.is_symmetric() << std::endl;

    /* Making b */
    Array<float3> dfdx_v0 = Array<float3>(n_vertices);

    dfdx.multiply(v0, dfdx_v0);

    blender::multiply_float_inplace(h, dfdx_v0);
    blender::add(dfdx_v0, f0);

    for (int i : IndexRange(n_vertices)) {
      dfdx_v0[i] = M_inv[i] * dfdx_v0[i];
    }

    blender::multiply_float_inplace(h, dfdx_v0);
    Array<float3> &b = dfdx_v0;

    /* Making A */
    dfdx.multiply_block_diagonal_matrix(M_inv, h * h);  // h^2 * M_inv * dfdx
    SparseMatrix &A = dfdx;

    for (int i : IndexRange(A.n_rows)) {
      A.insert(i, i, float3x3::identity() - A.get(i, i));
    }

    /* Solving the system. */
    Array<float3> delta_v = Array<float3>(n_vertices, float3(0.0f));
    solve_pcg(A, b, delta_v);

    std::cout << "A symmetric: " << A.is_symmetric() << std::endl;
    // solve_gauss_seidel(A, b, delta_v);

    for (int i : IndexRange(n_vertices)) {
      vertex_velocities[i] += delta_v[i];
      vertex_positions[i] += vertex_velocities[i] * substep_time;
    }
  }

  void integrate_implicit_backward_euler_pcg_filtered()
  {
    /* Assemble linear system */
    float h = substep_time;
    SparseMatrix &dfdx = vertex_force_derivatives;
    Array<float3> &v0 = vertex_velocities;
    Array<float3> &f0 = vertex_forces;

    /* Making b */
    Array<float3> dfdx_v0 = Array<float3>(n_vertices);
    dfdx.multiply(v0, dfdx_v0);
    blender::multiply_float_inplace(h, dfdx_v0);
    blender::add(dfdx_v0, f0);
    blender::multiply_float_inplace(h, dfdx_v0);
    Array<float3> &b = dfdx_v0;

    /* Making A */
    SparseMatrix &A = dfdx;
    A.multiply_float(h * h);
    for (int i : IndexRange(A.n_rows)) {
      float3x3 mass_matrix = float3x3(0.0f);
      mass_matrix.values[0][0] = vertex_masses[i];
      mass_matrix.values[1][1] = vertex_masses[i];
      mass_matrix.values[2][2] = vertex_masses[i];
      A.insert(i, i, mass_matrix - A.get(i, i));
    }

    /* Solving the system. */
    Array<float3> delta_v = Array<float3>(n_vertices, float3(0.0f));
    solve_pcg_filtered(A, b, delta_v);

    std::cout << "A symmetric: " << A.is_symmetric() << std::endl;
    std::cout << "delta_v: " << delta_v[0] << " " << delta_v[1] << std::endl;

    for (int i : IndexRange(n_vertices)) {
      vertex_velocities[i] += delta_v[i];
      vertex_positions[i] += vertex_velocities[i] * substep_time;
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
    int i = triangle_index;

    float area_uv = triangle_areas_uv[i];
    float wu_norm = wu.normalize_and_get_length();
    float wv_norm = wv.normalize_and_get_length();

    /* Stretch condition: Equation (10) in BW98. */
    float Cu = area_uv * (wu_norm - 1.0);
    float Cv = area_uv * (wv_norm - 1.0);

    int3 vertex_indices = triangle_vertex_indices[i];
    float3 dwu_dx = triangle_wu_derivatives[i];
    float3 dwv_dx = triangle_wv_derivatives[i];

    float ku = triangle_stretch_stiffness_u[i];
    float kv = triangle_stretch_stiffness_v[i];

    Array<float3> dCu_dx = Array<float3>(3);
    Array<float3> dCv_dx = Array<float3>(3);

    /* Forces */
    for (int m : IndexRange(3)) {
      int vertex_index = vertex_indices[m];

      dCu_dx[m] = area_uv * dwu_dx[m] * wu;
      dCv_dx[m] = area_uv * dwv_dx[m] * wv;

      float3 force_m_u = -ku * Cu * dCu_dx[m];
      float3 force_m_v = -kv * Cv * dCv_dx[m];

      vertex_forces[vertex_index] += force_m_u + force_m_v;
    }

    /* Force derivatives */
    float3x3 I_wu_wuT = float3x3::identity() - float3x3::outer(wu, wu);
    float3x3 I_wv_wvT = float3x3::identity() - float3x3::outer(wv, wv);

    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dCu_dx_mn = area_uv / wu_norm * dwu_dx[m] * dwu_dx[n] * I_wu_wuT;
        float3x3 dfu_dx_mn = -ku * float3x3::outer(dCu_dx[m], dCu_dx[n]) + Cu * dCu_dx_mn;

        float3x3 dCv_dx_mn = area_uv / wv_norm * dwv_dx[m] * dwv_dx[n] * I_wv_wvT;
        float3x3 dfv_dx_mn = -kv * float3x3::outer(dCv_dx[m], dCv_dx[n]) + Cv * dCv_dx_mn;

        int i = vertex_indices[m];
        int j = vertex_indices[n];
        vertex_force_derivatives.add(i, j, dfu_dx_mn + dfv_dx_mn);
      }
    }
  }

  ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin constructed" << std::endl;
  }

  ~ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin destructed" << std::endl;
    delete vertex_force_derivatives_pointer;
  }
};