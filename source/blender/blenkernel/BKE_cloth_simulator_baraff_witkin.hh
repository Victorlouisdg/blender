#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"
#include "BLI_float3x3.hh"
#include "BLI_int3.hh"
#include "BLI_int4.hh"
#include "BLI_sparse_matrix.hh"
#include "BLI_vector.hh"

#include "BLI_index_range.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_cloth_simulator_constants.hh"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"

#include <tuple>

extern "C" {
#include "draw_debug.h"
}

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
using blender::ConjugateGradientSolver;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using blender::int3;
using blender::int4;
using blender::Span;
using blender::SparseMatrix;
using blender::Vector;

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

  ConjugateGradientSolver solver;

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
  Array<float> triangle_shear_stiffness;

  Array<float> triangle_areas_uv;
  Array<float3> edges_normals;

  float density = 0.2f;            /* in kg/m^2 */
  float standard_gravity = -9.81f; /* in m/s^2 */

  /* Datastructures for intermediate computations. */
  SparseMatrix vertex_force_derivatives;  // <float3x3> elements

  Array<float3> vertex_forces;
  Array<float3> triangle_normals;

  /* Precomputed quantities. */
  Array<float2x2> triangle_inverted_delta_u_matrices;
  Array<float3> triangle_wu_derivatives;
  Array<float3> triangle_wv_derivatives;

  /* Kinematic Collision stuff */
  bool kinematic_collisions_enabled = false;
  int n_collision_vertices;
  int n_collision_triangles;

  Array<int3> collision_triangle_vertex_indices;
  Array<float3> collision_vertex_positions;
  Array<float3> collision_triangle_normals;
  Array<float3> collision_edge_normals;

  /* Maybe turn this into a class LoopTriEdge that overwrite the hash() function? */
  /* Currently storing looptri edges as a stored pair of vertex indices. */
  std::map<std::pair<int, int>, float3> looptri_edge_normals;

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

    vertex_force_derivatives = SparseMatrix(n_vertices);
    solver = ConjugateGradientSolver(n_vertices);

    // verify_w_derivatives(); /* For debugging, could become a test. */
  };

  void pin(int vertex_index)
  {
    solver.setConstraint(vertex_index, float3(0.0f), float3x3(0.0f));
  }

  void set_collision_mesh(const Mesh &collision_mesh)
  {
    kinematic_collisions_enabled = true;
    n_collision_vertices = collision_mesh.totvert;

    collision_vertex_positions = Array<float3>(n_collision_vertices);

    for (const int i : IndexRange(n_collision_vertices)) {
      MVert vertex = collision_mesh.mvert[i];
      collision_vertex_positions[i] = vertex.co;
    }

    /* Triangle attributes */
    Span<MLoopTri> looptris = get_mesh_looptris(collision_mesh);
    n_collision_triangles = looptris.size();

    int n_collision_edges = collision_mesh.totedge;

    collision_triangle_vertex_indices = Array<int3>(n_collision_triangles);
    collision_triangle_normals = Array<float3>(n_collision_triangles);
    collision_edge_normals = Array<float3>(n_collision_edges, float3(0.0f));

    looptri_edge_normals = std::map<std::pair<int, int>, float3>();

    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop_index = looptri.tri[0];
      const int v1_loop_index = looptri.tri[1];
      const int v2_loop_index = looptri.tri[2];
      const MLoop v0_loop = collision_mesh.mloop[v0_loop_index];
      const MLoop v1_loop = collision_mesh.mloop[v1_loop_index];
      const MLoop v2_loop = collision_mesh.mloop[v2_loop_index];
      const int v0_index = v0_loop.v;
      const int v1_index = v1_loop.v;
      const int v2_index = v2_loop.v;
      const float3 v0_pos = collision_vertex_positions[v0_index];
      const float3 v1_pos = collision_vertex_positions[v1_index];
      const float3 v2_pos = collision_vertex_positions[v2_index];

      collision_triangle_vertex_indices[looptri_index] = {v0_index, v1_index, v2_index};

      float3 normal = float3::cross(v1_pos - v0_pos, v2_pos - v0_pos).normalized();

      collision_triangle_normals[looptri_index] = normal;

      looptri_edge_normals[std::minmax(v0_index, v1_index)] += normal;
      looptri_edge_normals[std::minmax(v0_index, v2_index)] += normal;
      looptri_edge_normals[std::minmax(v1_index, v2_index)] += normal;
    }

    std::cout << "looptri edge normals" << std::endl;
    for (std::pair<std::pair<int, int>, float3> p : looptri_edge_normals) {
      std::cout << p.first.first << ", " << p.first.second << ": " << p.second << std::endl;
    }

    std::cout << "Collision mesh set." << std::endl;
  }

  void step()
  {
    for (int substep : IndexRange(n_substeps)) {
      UNUSED_VARS(substep);
      reset_forces_and_derivatives();

      calculate_kinematic_collisions();
      calculate_forces_and_derivatives();

      // integrate_explicit_forward_euler();
      // integrate_implicit_backward_euler_pcg_filtered();
    }

    // calculate_kinematic_collisions(); /* Placed here for debugging. */
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
    float shear_stiffness = 100.0f;

    triangle_stretch_stiffness_u = Array<float>(n_triangles, stretch_stiffness);
    triangle_stretch_stiffness_v = Array<float>(n_triangles, stretch_stiffness);
    triangle_shear_stiffness = Array<float>(n_triangles, shear_stiffness);

    /* Maps from the 2 vertices of a looptri edge to an oppossing vertex. */
    std::multimap<std::pair<int, int>, int> looptri_edge_opposing_vertices =
        std::multimap<std::pair<int, int>, int>();

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

      /* Bending edges. */
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v0_index, v1_index), v2_index));
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v0_index, v2_index), v1_index));
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v1_index, v2_index), v0_index));
    }

    Vector<int4> bending_indices_vector = Vector<int4>();

    /* Here we need to check which looptri edges have 2 opposing vertices, these edges will
     * become bending edges. */
    for (std::multimap<std::pair<int, int>, int>::iterator it =
             looptri_edge_opposing_vertices.begin();
         it != looptri_edge_opposing_vertices.end();
         it++) {
      std::pair<int, int> edge = it->first;

      /* Note that indices of the shared edge are stored in position 1 and 2, this is the
       * convention used in the "Implementing Baraff Wikin" by David Pritchard. */
      int v0 = it->second;
      int v1 = edge.first;
      int v2 = edge.second;

      std::multimap<std::pair<int, int>, int>::iterator it_next = std::next(it, 1);

      if (it_next != looptri_edge_opposing_vertices.end()) {
        std::pair<int, int> edge_next = it_next->first;
        if (edge == edge_next) { /* This means two looptris share this edge. */
          int v3 = it_next->second;
          bending_indices_vector.append(int4(v0, v1, v2, v3));
        }
      }
    }

    bending_vertex_indices = Array<int4>(bending_indices_vector.as_span());
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
      auto [wu, wv] = calculate_w_uv(triangle_index);
      calculate_stretch(triangle_index, wu, wv);
      calculate_shear(triangle_index, wu, wv);
    }

    for (int bending_index : IndexRange(bending_vertex_indices.size())) {
      calculate_bend(bending_index);
    }
  }

  void integrate_explicit_forward_euler()
  {
    /* vertex_accelerations are only needed for debugging with explicit integration. */
    Array<float3> vertex_accelerations = Array<float3>(n_vertices);

    /* TODO: look into doing these Array-based operations with std::transform() etc. */
    for (int i : IndexRange(n_vertices)) {
      vertex_accelerations[i] = 1.0f / vertex_masses[i] * vertex_forces[i];
      vertex_velocities[i] += vertex_accelerations[i] * substep_time;
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
    Array<float3> delta_v = Array<float3>(n_vertices);
    solver.solve(A, b, delta_v);

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

  void calculate_shear(int triangle_index, float3 wu, float3 wv)
  {
    int i = triangle_index;
    float area_uv = triangle_areas_uv[i];

    /* Stretch condition: section 4.3 in [BW98]. */
    float C = area_uv * float3::dot(wu, wv);

    int3 vertex_indices = triangle_vertex_indices[i];
    float3 dwu_dx = triangle_wu_derivatives[i];
    float3 dwv_dx = triangle_wv_derivatives[i];

    float k = triangle_shear_stiffness[i];

    Array<float3> dC_dx = Array<float3>(3);

    /* Forces */
    for (int m : IndexRange(3)) {
      int vertex_index = vertex_indices[m];

      dC_dx[m] = area_uv * (dwu_dx[m] * wv + dwv_dx[m] * wu);
      float3 force_m = -k * C * dC_dx[m];
      vertex_forces[vertex_index] += force_m;
    }

    /* Force derivatives */
    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dC_dx_mn = area_uv * (dwu_dx[m] * dwv_dx[n] + dwu_dx[n] * dwv_dx[m]) *
                            float3x3::identity();
        float3x3 df_dx_mn = -k * (float3x3::outer(dC_dx[m], dC_dx[n]) + C * dC_dx_mn);

        int i = vertex_indices[m];
        int j = vertex_indices[n];
        vertex_force_derivatives.add(i, j, df_dx_mn);
      }
    }
  }

  void calculate_bend(int bending_index)
  {
    auto [v0_index, v1_index, v2_index, v3_index] = bending_vertex_indices[bending_index];
    float3 x0 = vertex_positions[v0_index];
    float3 x1 = vertex_positions[v1_index];
    float3 x2 = vertex_positions[v2_index];
    float3 x3 = vertex_positions[v3_index];

    const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};
    DRW_debug_line_v3v3(x0, x3, green);
  }

  void calculate_kinematic_collisions()
  {
    Array<float> collision_distances = Array<float>(n_vertices, FLT_MAX);
    Array<float3> collision_closest_points = Array<float3>(n_vertices);

    /* For each cloth vertex find the closest point on the collision geometry. */
    for (int i : IndexRange(n_vertices)) {
      float3 x = vertex_positions[i];

      /* Currently we search for collision by looping over all triangles, later this should be
       * accelerated with an AABB tree or similar. */
      for (int t : IndexRange(n_collision_triangles)) {
        int3 vertex_indices = collision_triangle_vertex_indices[t];

        float3 v0 = collision_vertex_positions[vertex_indices[0]];
        float3 v1 = collision_vertex_positions[vertex_indices[1]];
        float3 v2 = collision_vertex_positions[vertex_indices[2]];

        float3 closest_point = float3();
        closest_on_tri_to_point_v3(closest_point, x, v0, v1, v2);

        float distance = float3::distance(closest_point, x);

        if (distance < collision_distances[i]) {
          collision_distances[i] = distance;
          collision_closest_points[i] = closest_point;
        }
      }
    }

    const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};
    const float red[4] = {1.0f, 0.0f, 0.0f, 1.0f};

    float offset = 0.1f;

    for (int i : IndexRange(n_vertices)) {
      float3 x = vertex_positions[i];
      float3 closest_point = collision_closest_points[i];
      float distance = collision_distances[i];

      if (distance > offset) {
        // DRW_debug_line_v3v3(x, closest_point, green);
      }
      else {
        // DRW_debug_line_v3v3(x, closest_point, red);

        // Stick particle to surface -> fully constraint delta_v to 0
        solver.setConstraint(i,
                             -vertex_velocities[i],
                             float3x3(0.0f)); /* Fully constrain particle velocity to zero. */
      }
    }
  }
};