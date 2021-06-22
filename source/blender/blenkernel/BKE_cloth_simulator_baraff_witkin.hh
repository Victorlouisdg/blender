#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"
#include "BLI_float3x3.hh"
#include "BLI_int2.hh"
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
#include <unordered_set>

extern "C" {
#include "draw_debug.h"
}

// #include <eigen3/Eigen/Dense>

/* Note about variable names in the simulator: in general I've tried to give the member
 * variables/attributes long and clear names. I deliberately read these out into local variable
 * with much shorter (often 2 letter) names. This is because long variable names would obfuscate
 * all the structure of the math operations.
 *
 * For the loop index variables in this file I've mantained this convention:
 * i, j: index the Arrays of all vertices
 * ti: indexes the Array of all triangles
 * m, n: index the vertices of a single triangle (or 2 triangles for the bending condition)
 * s, t: index the components (x, y, z) of a vector
 */

using blender::Array;
using blender::ConjugateGradientSolver;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using blender::int2;
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
  bool enable_shear;
  bool enable_bending;

  bool damp_stretch;
  bool damp_shear;
  bool damp_bending;
  bool damp_springs;

  bool use_explicit_integration;

  float stretch_damping_factor;
  float shear_damping_factor;
  float bending_damping_factor;
  float spring_damping_factor;

  int n_vertices;
  int n_triangles;
  int n_bending_edges;

  int n_substeps;
  float step_time; /* Currently assuming 30 fps. */
  float substep_time;

  int current_substep; /* Mostly used for debug drawing at last substep. */

  ConjugateGradientSolver solver;

  /* State */
  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  /* Indices */
  Array<int2> spring_vertex_indices;
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
  Array<float> bending_stiffness;
  Array<float> bending_rest_lengths;

  Array<float> triangle_areas_uv_square_root;
  Array<float3> edges_normals;

  int n_springs;
  Array<float> spring_rest_lengths;
  Array<float> spring_stiffness;

  float density = 0.2f;            /* in kg/m^2 */
  float standard_gravity = -9.81f; /* in m/s^2 */

  Vector<int> pinned_vertices;

  /* Datastructures for intermediate computations. */
  SparseMatrix vertex_force_derivatives;           // <float3x3> elements
  SparseMatrix vertex_force_velocity_derivatives;  // <float3x3> elements

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

  Array<float3> vertex_position_alterations;

  std::map<int, float3> collision_constrained_vertices_directions;

  /* Maybe turn this into a class LoopTriEdge that overwrite the hash() function? */
  /* Currently storing looptri edges as a stored pair of vertex indices. */
  std::map<std::pair<int, int>, float3> looptri_edge_normals;

  /* Currently the simulator has its own local copy of all the necessary
   * mesh attributes e.g. vertex positions. This was easier to implement and I believe could make
   * the simulator faster due to better data locality. I'll make setters() for certain attributes
   * later that do the required recalculations.
   */
  void initialize(const Mesh &mesh,
                  int n_substeps,
                  float stretch_stiffness_value,
                  float shear_stiffness_value,
                  float bending_stiffness_value,
                  float spring_stiffness_value,
                  float stretch_damping_factor,
                  float shear_damping_factor,
                  float bending_damping_factor,
                  float spring_damping_factor,
                  bool enable_shear,
                  bool enable_bending,
                  bool damp_stretch,
                  bool damp_shear,
                  bool damp_bending,
                  bool damp_springs,
                  bool use_explicit_integration)
  {
    std::cout << "Cloth simulation initialisation" << std::endl;

    this->n_substeps = n_substeps;
    step_time = 1.0f / 30.0f; /* Currently assuming 30 fps. */
    substep_time = step_time / n_substeps;

    n_vertices = mesh.totvert;

    this->stretch_damping_factor = stretch_damping_factor;
    this->shear_damping_factor = shear_damping_factor;
    this->bending_damping_factor = bending_damping_factor;
    this->spring_damping_factor = spring_damping_factor;

    this->enable_shear = enable_shear;
    this->enable_bending = enable_bending;
    this->damp_stretch = damp_stretch;
    this->damp_shear = damp_shear;
    this->damp_bending = damp_bending;
    this->damp_springs = damp_springs;
    this->use_explicit_integration = use_explicit_integration;

    /* TODO move this to "cloth_attribute_initialisation.cc" or something */
    initialize_vertex_attributes(mesh);
    initialize_triangle_attributes(
        mesh, stretch_stiffness_value, shear_stiffness_value, bending_stiffness_value);

    vertex_force_derivatives = SparseMatrix(n_vertices);
    vertex_force_velocity_derivatives = SparseMatrix(n_vertices);
    solver = ConjugateGradientSolver(n_vertices);

    collision_constrained_vertices_directions = std::map<int, float3>();
    pinned_vertices = Vector<int>();

    // verify_w_derivatives(); /* For debugging, could become a test. */

    /* Finding all edges that are not an edge of a face. */
    std::unordered_set<int> non_face_vertices = std::unordered_set<int>();

    for (int i : IndexRange(n_vertices)) {
      non_face_vertices.insert(i);
    }

    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop_index = looptri.tri[0];
      const int v1_loop_index = looptri.tri[1];
      const int v2_loop_index = looptri.tri[2];
      const MLoop v0_loop = mesh.mloop[v0_loop_index];
      const MLoop v1_loop = mesh.mloop[v1_loop_index];
      const MLoop v2_loop = mesh.mloop[v2_loop_index];
      const int v0_index = v0_loop.v;
      const int v1_index = v1_loop.v;
      const int v2_index = v2_loop.v;

      /* TODO simplify. */
      auto search0 = non_face_vertices.find(v0_index);
      if (search0 != non_face_vertices.end()) {
        non_face_vertices.erase(search0);
      }
      auto search1 = non_face_vertices.find(v1_index);
      if (search1 != non_face_vertices.end()) {
        non_face_vertices.erase(search1);
      }
      auto search2 = non_face_vertices.find(v2_index);
      if (search2 != non_face_vertices.end()) {
        non_face_vertices.erase(search2);
      }
    }

    std::cout << "Non face vertices" << std::endl;
    for (const int vertex_index : non_face_vertices) {
      std::cout << vertex_index << std::endl;
    }

    Vector<int2> edges_not_in_any_face = Vector<int2>();

    for (int i : IndexRange(mesh.totedge)) {
      const MEdge &edge = mesh.medge[i];

      if (non_face_vertices.find(edge.v1) != non_face_vertices.end() ||
          non_face_vertices.find(edge.v2) != non_face_vertices.end()) {
        edges_not_in_any_face.append(int2(edge.v1, edge.v2));
      }
    }

    spring_vertex_indices = Array<int2>(edges_not_in_any_face.as_span());
    n_springs = spring_vertex_indices.size();
    spring_rest_lengths = Array<float>(n_springs);
    spring_stiffness = Array<float>(n_springs);

    for (int i : IndexRange(n_springs)) {
      int2 spring = spring_vertex_indices[i];
      float3 x0 = vertex_positions[spring.i];
      float3 x1 = vertex_positions[spring.j];

      spring_rest_lengths[i] = (x0 - x1).length();
      spring_stiffness[i] = spring_stiffness_value;
    }
  };

  void pin(int vertex_index)
  {
    pinned_vertices.append(vertex_index);
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
  }

  void step()
  {

    std::cout << "step" << std::endl;
    for (int substep : IndexRange(n_substeps)) {
      current_substep = substep;
      reset_forces_and_derivatives();

      calculate_kinematic_collisions();
      calculate_forces_and_derivatives();

      if (use_explicit_integration) {
        integrate_explicit_forward_euler();
      }
      else {
        integrate_implicit_backward_euler_pcg_filtered();
      }
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

    vertex_position_alterations = Array<float3>(n_vertices, float3(0.0f));

    for (const int i : IndexRange(n_vertices)) {
      MVert vertex = mesh.mvert[i];
      vertex_positions[i] = vertex.co;
    }
  }

  void initialize_triangle_attributes(const Mesh &mesh,
                                      float stretch_stiffness_value,
                                      float shear_stiffness_value,
                                      float bending_stiffness_value)
  {
    /* TODO set vertex mass based on Voronoi area of vertex. I read this in a paper somewhere but
     * can't find it. */
    float fixed_vertex_mass = 1.0f / n_vertices;
    vertex_masses = Array<float>(n_vertices, fixed_vertex_mass);

    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    n_triangles = looptris.size();

    triangle_vertex_indices = Array<int3>(n_triangles);
    triangle_areas_uv_square_root = Array<float>(n_triangles);
    triangle_normals = Array<float3>(n_triangles);
    triangle_inverted_delta_u_matrices = Array<float2x2>(n_triangles);
    triangle_wu_derivatives = Array<float3>(n_triangles);
    triangle_wv_derivatives = Array<float3>(n_triangles);

    triangle_stretch_stiffness_u = Array<float>(n_triangles, stretch_stiffness_value);
    triangle_stretch_stiffness_v = Array<float>(n_triangles, stretch_stiffness_value);
    triangle_shear_stiffness = Array<float>(n_triangles, shear_stiffness_value);

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
      triangle_areas_uv_square_root[looptri_index] = sqrtf(
          triangle_area_uv); /* Sqrt such that energy is linear in area (see Bhat and Pritchard)*/

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
    Vector<float> bending_rest_lengths_vector = Vector<float>();

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
          float3 x0 = vertex_positions[v0];
          float3 x3 = vertex_positions[v3];
          float rest_length = (x3 - x0).length();
          bending_rest_lengths_vector.append(rest_length);
        }
      }
    }

    bending_vertex_indices = Array<int4>(bending_indices_vector.as_span());
    bending_rest_lengths = Array<float>(bending_rest_lengths_vector.as_span());
    n_bending_edges = bending_vertex_indices.size();

    bending_stiffness = Array<float>(n_bending_edges, bending_stiffness_value);
  };

  void calculate_forces_and_derivatives()
  {
    for (int vertex_index : IndexRange(n_vertices)) {
      calculate_gravity(vertex_index);
    }

    for (int spring_index : IndexRange(n_springs)) {
      calculate_spring(spring_index);
    }

    for (int triangle_index : IndexRange(n_triangles)) {
      auto [wu, wv] = calculate_w_uv(triangle_index);
      calculate_stretch(triangle_index, wu, wv);
      if (enable_shear) {
        calculate_shear(triangle_index, wu, wv);
      }
    }

    if (enable_bending) {
      for (int bending_index : IndexRange(bending_vertex_indices.size())) {
        calculate_bend(bending_index);
      }
    }
  }

  void integrate_explicit_forward_euler()
  {
    /* vertex_accelerations are only needed for debugging with explicit integration. */
    Array<float3> vertex_accelerations = Array<float3>(n_vertices);

    /* TODO: look into doing these Array-based operations with std::transform() etc. */
    for (int i : IndexRange(n_vertices)) {
      vertex_accelerations[i] = 1.0f / vertex_masses[i] * vertex_forces[i];
    }

    for (int i : pinned_vertices) {
      vertex_accelerations[i] = float3(0.0f);
    }

    for (int i : IndexRange(n_vertices)) {
      vertex_velocities[i] += vertex_accelerations[i] * substep_time;
      vertex_positions[i] += vertex_velocities[i] * substep_time;
    }
  }

  void integrate_implicit_backward_euler_pcg_filtered()
  {
    /* This function builds and solve the system from equation (18) in BW98. */

    /* Assemble linear system */
    float h = substep_time;
    SparseMatrix &dfdx = vertex_force_derivatives;
    SparseMatrix &dfdv = vertex_force_velocity_derivatives;

    Array<float3> &v0 = vertex_velocities;
    Array<float3> &f0 = vertex_forces;
    Array<float3> &y = vertex_position_alterations;

    /* Making b */
    Array<float3> dfdx_v0 = Array<float3>(n_vertices);
    Array<float3> dfdx_y = Array<float3>(n_vertices);

    dfdx.multiply(y, dfdx_y);

    // std::cout << "dfdx_y" << std::endl;
    // for (int i : IndexRange(n_vertices)) {
    //   std::cout << dfdx_y[i] << std::endl;
    // }

    dfdx.multiply(v0, dfdx_v0);
    blender::multiply_float_inplace(h, dfdx_v0);
    blender::add(dfdx_v0, f0);
    blender::add(dfdx_v0, dfdx_y);
    blender::multiply_float_inplace(h, dfdx_v0);
    Array<float3> &b = dfdx_v0;

    /* Making A */
    SparseMatrix &A = dfdx;
    A.multiply_float(h * h);
    dfdv.multiply_float(h);
    A.add_matrix(dfdv);

    /* A = M - A */
    A.multiply_float(-1.0f);
    for (int i : IndexRange(A.n_rows)) {
      float3x3 mass_matrix = float3x3::diagonal(vertex_masses[i]);
      A.insert(i, i, mass_matrix + A.get(i, i));
    }

    /* Solving the system. */
    Array<float3> delta_v = Array<float3>(n_vertices);
    solver.solve(A, b, delta_v);

    /* Forces a resulting from the enforcement of constraints. */
    Array<float3> e = Array<float3>(n_vertices);
    A.multiply(delta_v, e); /* e = Ax */
    subtract_from(e, b);    /* e = b - Ax */

    /* Releasing constraints that pull towards the surface. */
    std::map<int, float3> &m = collision_constrained_vertices_directions;
    for (auto it = m.cbegin(); it != m.cend(); /* no increment */) {

      int i = it->first;
      float3 constrained_direction = it->second;
      float3 constraining_force = e[i];

      bool force_pulling_towards_surface = float3::dot(constrained_direction, constraining_force) <
                                           0.0f;

      if (force_pulling_towards_surface) {
        solver.releaseConstraint(i);
        m.erase(it++);
      }
      else {
        ++it;
      }
    }

    for (int i : IndexRange(n_vertices)) {
      vertex_velocities[i] += delta_v[i] + y[i];
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

  float3 spring_force(float3 x0, float3 x1, float k, float rest_length)
  {
    float3 spring_direction = x1 - x0; /* Points towards x1 */
    float length = spring_direction.normalize_and_get_length();

    float3 force = k * (length - rest_length) * spring_direction;
    return force;
  }

  float3x3 spring_force_jacobian(float3 x0, float3 x1, float k, float rest_length)
  {
    float3 spring_direction = x1 - x0; /* Points towards x1 */
    float length = spring_direction.normalize_and_get_length();
    float3x3 spring_outer = 1.0f / float3::dot(spring_direction, spring_direction) *
                            float3x3::outer(spring_direction, spring_direction);

    /* TODO investigate setting jacobian to zero when in compression */
    float3x3 force_derivative = -k * spring_outer - k * (1 - rest_length / length) *
                                                        (float3x3::identity() - spring_outer);
    return force_derivative;
  }

  void calculate_spring(int spring_index)
  {
    /* 1D spring experiment. */
    int2 spring = spring_vertex_indices[spring_index];
    int i = spring.i;
    int j = spring.j;
    float3 x0 = vertex_positions[i];
    float3 x1 = vertex_positions[j];

    float rest_length = spring_rest_lengths[spring_index];
    float k = spring_stiffness[spring_index];

    float3 spring_direction = x1 - x0; /* Points towards x1 */
    float length = spring_direction.normalize_and_get_length();
    bool spring_is_in_compression = length < rest_length;

    float3 force = spring_force(x0, x1, k, rest_length);
    float3x3 force_derivative = spring_force_jacobian(x0, x1, k, rest_length);

    vertex_forces[i] += force;
    vertex_forces[j] -= force;

    vertex_force_derivatives.add(i, i, force_derivative);
    vertex_force_derivatives.add(j, j, force_derivative);
    vertex_force_derivatives.add(i, j, -1.0f * force_derivative);
    vertex_force_derivatives.add(j, i, -1.0f * force_derivative);

    if (damp_springs) {
      float3 v0 = vertex_velocities[i];
      float3 v1 = vertex_velocities[j];
      float3 velocity_difference = v1 - v0;

      float kd = spring_damping_factor * k;
      float3 spring_direction = x1 - x0; /* Points towards x1 */
      float length = spring_direction.normalize_and_get_length();

      float3 damping_force = kd * float3::dot(velocity_difference, spring_direction) *
                             spring_direction;

      float3x3 damping_derivative = -kd * float3x3::outer(spring_direction, spring_direction);

      vertex_forces[i] += damping_force;
      vertex_forces[j] -= damping_force;

      vertex_force_velocity_derivatives.add(i, i, damping_derivative);
      vertex_force_velocity_derivatives.add(j, j, damping_derivative);
      vertex_force_velocity_derivatives.add(i, j, -1.0f * damping_derivative);
      vertex_force_velocity_derivatives.add(j, i, -1.0f * damping_derivative);
    }
  }

  void calculate_stretch(int triangle_index, float3 wu, float3 wv)
  {
    int ti = triangle_index;

    float area_uv = triangle_areas_uv_square_root[ti];
    float wu_norm = wu.normalize_and_get_length();
    float wv_norm = wv.normalize_and_get_length();

    /* Stretch condition: Equation (10) in BW98. */
    float Cu = area_uv * (wu_norm - 1.0);
    float Cv = area_uv * (wv_norm - 1.0);

    int3 vertex_indices = triangle_vertex_indices[ti];
    float3 dwu_dx = triangle_wu_derivatives[ti];
    float3 dwv_dx = triangle_wv_derivatives[ti];

    float ku = triangle_stretch_stiffness_u[ti];
    float kv = triangle_stretch_stiffness_v[ti];

    /* Temporary values. */
    float kdu = ku * stretch_damping_factor;
    float kdv = kv * stretch_damping_factor;

    Array<float3> dCu_dx = Array<float3>(3);
    Array<float3> dCv_dx = Array<float3>(3);

    float Cu_dot = 0.0f;
    float Cv_dot = 0.0f;

    /* Forces */
    for (int m : IndexRange(3)) {
      dCu_dx[m] = area_uv * dwu_dx[m] * wu;
      dCv_dx[m] = area_uv * dwv_dx[m] * wv;

      if (damp_stretch) {
        int i = vertex_indices[m];
        Cu_dot += float3::dot(dCu_dx[m], vertex_velocities[i]);
        Cv_dot += float3::dot(dCv_dx[m], vertex_velocities[i]);
      }
    }

    for (int m : IndexRange(3)) {
      float3 force_m_u = -ku * Cu * dCu_dx[m];
      float3 force_m_v = -kv * Cv * dCv_dx[m];

      int i = vertex_indices[m];
      vertex_forces[i] += force_m_u + force_m_v;

      if (damp_stretch) {
        float3 damping_force_m_u = -kdu * Cu_dot * dCu_dx[m];
        float3 damping_force_m_v = -kdv * Cv_dot * dCv_dx[m];
        vertex_forces[i] += damping_force_m_u + damping_force_m_v;
      }
    }

    /* Force derivatives */
    float3x3 I_wu_wuT = float3x3::identity() - float3x3::outer(wu, wu);
    float3x3 I_wv_wvT = float3x3::identity() - float3x3::outer(wv, wv);

    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dCu_dx_mn = area_uv / wu_norm * dwu_dx[m] * dwu_dx[n] * I_wu_wuT;
        float3x3 dCu_outer = float3x3::outer(dCu_dx[m], dCu_dx[n]);
        float3x3 dfu_dx_mn = -ku * (dCu_outer + Cu * dCu_dx_mn);

        float3x3 dCv_dx_mn = area_uv / wv_norm * dwv_dx[m] * dwv_dx[n] * I_wv_wvT;
        float3x3 dCv_outer = float3x3::outer(dCv_dx[m], dCv_dx[n]);
        float3x3 dfv_dx_mn = -kv * (dCv_outer + Cv * dCv_dx_mn);

        int i = vertex_indices[m];
        int j = vertex_indices[n];
        vertex_force_derivatives.add(i, j, dfu_dx_mn + dfv_dx_mn);

        // std::cout << std::endl << "Stretch eigenvalues: " << std::endl;
        // std::cout << dfu_dx_mn.eigenvalues() << std::endl;
        // std::cout << dfv_dx_mn.eigenvalues() << std::endl;

        if (damp_stretch) {
          float3x3 ddu_dx_mn = -kdu * Cu_dot * dCu_dx_mn;
          float3x3 ddv_dx_mn = -kdv * Cv_dot * dCv_dx_mn;
          vertex_force_derivatives.add(i, j, ddu_dx_mn + ddv_dx_mn);

          float3x3 ddu_dv_nm = -kdu * dCu_outer;
          float3x3 ddv_dv_nm = -kdv * dCv_outer;
          vertex_force_velocity_derivatives.add(i, j, ddu_dv_nm + ddv_dv_nm);
        }
      }
    }
  }

  void calculate_shear(int triangle_index, float3 wu, float3 wv)
  {
    int ti = triangle_index;
    float area_uv = triangle_areas_uv_square_root[ti];

    /* Shear condition: section 4.3 in [BW98]. */
    float C = area_uv * float3::dot(wu, wv);

    int3 vertex_indices = triangle_vertex_indices[ti];
    float3 dwu_dx = triangle_wu_derivatives[ti];
    float3 dwv_dx = triangle_wv_derivatives[ti];

    float k = triangle_shear_stiffness[ti];
    float kd = k * shear_damping_factor;

    Array<float3> dC_dx = Array<float3>(3);
    float C_dot = 0.0f;

    /* Forces */
    for (int m : IndexRange(3)) {
      dC_dx[m] = area_uv * (dwu_dx[m] * wv + dwv_dx[m] * wu);

      if (damp_shear) {
        int i = vertex_indices[m];
        C_dot += float3::dot(dC_dx[m], vertex_velocities[i]);
      }
    }

    for (int m : IndexRange(3)) {
      float3 force_m = -k * C * dC_dx[m];

      int i = vertex_indices[m];
      vertex_forces[i] += force_m;

      if (damp_shear) {
        float3 damping_force_m = -kd * C_dot * dC_dx[m];
        vertex_forces[i] += damping_force_m;
      }
    }

    /* Force derivatives */
    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dC_dx_mn = area_uv * (dwu_dx[m] * dwv_dx[n] + dwu_dx[n] * dwv_dx[m]) *
                            float3x3::identity();

        float3x3 dC_outer = float3x3::outer(dC_dx[m], dC_dx[n]);
        float3x3 df_dx_mn = -k * (dC_outer + C * dC_dx_mn);

        int i = vertex_indices[m];
        int j = vertex_indices[n];
        vertex_force_derivatives.add(i, j, df_dx_mn);

        if (damp_shear) {
          float3x3 dd_dx_mn = -kd * C_dot * dC_dx_mn;
          vertex_force_derivatives.add(i, j, dd_dx_mn);

          float3x3 dd_dv_nm = -kd * dC_outer;
          vertex_force_velocity_derivatives.add(i, j, dd_dv_nm);
        }
      }
    }
  }

  void calculate_bend_blender(int bending_index)
  {
    float k = bending_stiffness[bending_index];

    int4 vertex_indices = bending_vertex_indices[bending_index];
    float3 x0 = vertex_positions[vertex_indices[0]];
    float3 x1 = vertex_positions[vertex_indices[1]];
    float3 x2 = vertex_positions[vertex_indices[2]];
    float3 x3 = vertex_positions[vertex_indices[3]];

    float3 e = x1 - x2;
    float3 nA = float3::cross(x2 - x0, x1 - x0);
    float3 nB = float3::cross(x1 - x3, x2 - x3);

    float e_norm = e.normalize_and_get_length();
    float nA_norm = nA.normalize_and_get_length();
    float nB_norm = nB.normalize_and_get_length();

    float cos = float3::dot(nA, nB);
    float sin = float3::dot(float3::cross(nA, nB), e);

    float angle = atan2f(sin, cos);

    float force_magnitude = k * angle;

    float3 forceA = force_magnitude * nA;
    float3 forceB = force_magnitude * nB;

    vertex_forces[vertex_indices[0]] += forceA;
    vertex_forces[vertex_indices[3]] += forceB;

    float3 force_avg = 0.5 * (forceA + forceB);

    vertex_forces[vertex_indices[1]] -= force_avg;
    vertex_forces[vertex_indices[2]] -= force_avg;

    // const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};

    // if (current_substep == n_substeps - 1) {
    //   DRW_debug_line_v3v3(x0, x0 + forceA, green);
    //   DRW_debug_line_v3v3(x3, x3 + forceB, green);
    //   DRW_debug_line_v3v3(x1, x1 - force_avg, green);
    //   DRW_debug_line_v3v3(x2, x2 - force_avg, green);
    // }
  }

  void calculate_bend(int bending_index)
  {
    float k = bending_stiffness[bending_index];
    float kd = k * bending_damping_factor;

    int4 vertex_indices = bending_vertex_indices[bending_index];
    float3 x0 = vertex_positions[vertex_indices[0]];
    float3 x1 = vertex_positions[vertex_indices[1]];
    float3 x2 = vertex_positions[vertex_indices[2]];
    float3 x3 = vertex_positions[vertex_indices[3]];

    float3 e = x1 - x2;
    float3 nA = float3::cross(x2 - x0, x1 - x0);
    float3 nB = float3::cross(x1 - x3, x2 - x3);

    float e_norm = e.normalize_and_get_length();
    float nA_norm = nA.normalize_and_get_length();
    float nB_norm = nB.normalize_and_get_length();

    float cos = float3::dot(nA, nB);
    float sin = float3::dot(float3::cross(nA, nB), e);

    float C = atan2f(sin, cos);

    float3 qA[4] = {x2 - x1, x0 - x2, x1 - x0, float3(0.0f)};
    float3 qB[4] = {float3(0.0f), x2 - x3, x3 - x1, x1 - x2};
    float3x3 qIe[4] = {float3x3(0.0f),
                       1.0f / e_norm * float3x3::diagonal(1.0f),
                       1.0f / e_norm * float3x3::diagonal(-1.0f),
                       float3x3(0.0f)};

    float3 dC_dx[4];
    float3 de_dx[4][3];
    float3 dnA_dx[4][3];
    float3 dnB_dx[4][3];
    float dcos_dx[4][3];
    float dsin_dx[4][3];

    float C_dot = 0.0f;

    /* Forces */
    for (int m : IndexRange(4)) {
      float3x3 SA = 1.0f / nA_norm * float3x3::skew(qA[m]);
      float3x3 SB = 1.0f / nB_norm * float3x3::skew(qB[m]);

      for (int s : IndexRange(3)) {
        de_dx[m][s] = qIe[m].row(s);
        dnA_dx[m][s] = SA.row(s);
        dnB_dx[m][s] = SB.row(s);
        dcos_dx[m][s] = float3::dot(dnA_dx[m][s], nB) + float3::dot(nA, dnB_dx[m][s]);
        dsin_dx[m][s] = float3::dot(
                            float3::cross(dnA_dx[m][s], nB) + float3::cross(nA, dnB_dx[m][s]), e) +
                        float3::dot(float3::cross(nA, nB), de_dx[m][s]);
        dC_dx[m][s] = cos * dsin_dx[m][s] - sin * dcos_dx[m][s];
      }

      if (damp_bending) {
        int i = vertex_indices[m];
        C_dot += float3::dot(dC_dx[m], vertex_velocities[i]);
      }
    }

    for (int m : IndexRange(4)) {
      float3 force_m = -k * C * dC_dx[m];

      int i = vertex_indices[m];
      vertex_forces[i] += force_m;

      if (damp_bending) {
        float3 damping_force_m = -kd * C_dot * dC_dx[m];
        vertex_forces[i] += damping_force_m;
      }
    }

    /* Force derivatives */

    float3x3 dC_dx_mn[4][4];

    for (int m : IndexRange(4)) {
      for (int n : IndexRange(4)) {
        for (int s : IndexRange(3)) {
          for (int t : IndexRange(3)) {

            /* TODO think over whether s and t might need to be swapped. */

            float3 dnA_dx_mnst = normalA_second_derivatives[m][n][s][t] / nA_norm;
            float3 dnB_dx_mnst = normalA_second_derivatives[m][n][s][t] / nB_norm;

            float3 dnA_dx_ms = dnA_dx[m][s];
            float3 dnA_dx_nt = dnA_dx[n][t];
            float3 dnB_dx_ms = dnB_dx[m][s];
            float3 dnB_dx_nt = dnB_dx[n][t];

            float dcos_dx_mnst = float3::dot(dnA_dx_mnst, nB) + float3::dot(dnB_dx_nt, dnA_dx_ms) +
                                 float3::dot(dnA_dx_nt, dnB_dx_ms) + float3::dot(nA, dnB_dx_mnst);

            float dsin_dx_mnst =
                float3::dot(float3::cross(dnA_dx_mnst, nB) + float3::cross(dnA_dx_ms, dnB_dx_nt) +
                                float3::cross(dnA_dx_nt, dnB_dx_ms) +
                                float3::cross(nA, dnB_dx_mnst),
                            e) +
                float3::dot(float3::cross(dnA_dx_ms, nB) + float3::cross(nA, dnB_dx_ms),
                            de_dx[n][t]) +
                float3::dot(float3::cross(dnA_dx_nt, nB) + float3::cross(nA, dnB_dx_nt),
                            de_dx[m][s]);

            /* t ans s are swapped here due to the column wise storage of small matrices in
             * blender. */
            dC_dx_mn[m][n].values[t][s] = dcos_dx[n][t] * dsin_dx[m][s] + cos * dsin_dx_mnst -
                                          dsin_dx[n][t] * dcos_dx[m][s] - sin * dcos_dx_mnst;
          }
        }

        float3x3 dC_outer = float3x3::outer(dC_dx[m], dC_dx[n]);
        float3x3 df_dx_mn = -k * (dC_outer + dC_dx_mn[m][n] * C);

        int i = vertex_indices[m];
        int j = vertex_indices[n];
        vertex_force_derivatives.add(i, j, df_dx_mn);

        if (damp_bending) {
          float3x3 dd_dx_mn = -kd * C_dot * dC_dx_mn[m][n];
          vertex_force_derivatives.add(i, j, dd_dx_mn);

          float3x3 dd_dv_nm = -kd * dC_outer;
          vertex_force_velocity_derivatives.add(i, j, dd_dv_nm);
        }
      }
    }
  }

  void calculate_kinematic_collisions()
  {
    Array<float> collision_distances = Array<float>(n_vertices, FLT_MAX);
    Array<float3> collision_closest_points = Array<float3>(n_vertices);
    Array<int> collision_closest_triangle_indices = Array<int>(n_vertices);

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
          collision_closest_triangle_indices[i] = t;
        }
      }
    }

    const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};
    const float red[4] = {1.0f, 0.0f, 0.0f, 1.0f};

    float offset = 0.02f;

    for (int i : IndexRange(n_vertices)) {
      float3 x = vertex_positions[i];
      float3 closest_point = collision_closest_points[i];
      float distance = collision_distances[i];

      if (distance > offset) {

        if (current_substep == n_substeps - 1) {
          // DRW_debug_line_v3v3(x, closest_point, green);
        }
      }
      else {
        if (current_substep == n_substeps - 1) {
          // DRW_debug_line_v3v3(x, closest_point, red);
        }
        int closest_triangle_index = collision_closest_triangle_indices[i];
        float3 normal = collision_triangle_normals[closest_triangle_index];

        float3x3 S = float3x3::identity() - float3x3::outer(normal, normal);
        /* TODO maybe multiply constraint by h as Pritchard suggests? */
        solver.setConstraint(i, -vertex_velocities[i], S);

        collision_constrained_vertices_directions[i] = normal;

        float3 desired_position = collision_closest_points[i] + normal * offset;
        float3 required_position_alteration = desired_position - vertex_positions[i];
        vertex_position_alterations[i] = required_position_alteration;

        std::cout << "Positions of " << i << " should be altered by "
                  << vertex_position_alterations[i] << std::endl;
      }
    }
  }

  void reset_forces_and_derivatives()
  {
    vertex_forces.fill(float3(0.0f));
    vertex_force_derivatives.clear();
    vertex_force_velocity_derivatives.clear();

    vertex_position_alterations.fill(float3(0.0f));
  }
};