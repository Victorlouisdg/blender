#include <iostream>

#include "BKE_cloth_simulator_bw_attribute_provider.hh"
#include "BKE_cloth_simulator_bw_force_elements.hh"

#include "BLI_index_range.hh"
#include "BLI_span.hh"
#include "BLI_timeit.hh"
#include "BLI_vector.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

extern "C" {
#include "draw_debug.h"
}

using blender::IndexRange;
using blender::MutableSpan;
using blender::Vector;

typedef Eigen::Triplet<float> Triplet;

using Eigen::ConjugateGradient;
using Eigen::Lower;
using Eigen::SparseMatrix;
using Eigen::Upper;
using Eigen::VectorXf;

class ClothSimulatorBW {
 public:
  int substeps;
  float step_time; /* Currently assuming 30 fps. */
  float substep_time;
  int current_substep;

  bool use_explicit_integration;

  int amount_of_vertices;
  int amount_of_triangles;
  int amount_of_bend_edges;
  int amount_of_springs;
  int system_size;

  float standard_gravity = -9.81f; /* in m/s^2 */

  ClothBWAttributeProvider ap;

  /* Indices */
  Span<int2> spring_vertex_indices;
  Span<int3> triangle_vertex_indices;
  Span<int4> bend_vertex_indices;

  MutableSpan<float3> vertex_positions;
  MutableSpan<float3> vertex_velocities;

  /* Read-Only Attributes, eventually provided by e.g. Geometry Nodes */
  Span<float> vertex_masses;

  Span<float> triangle_stretch_stiffness_u;
  Span<float> triangle_stretch_stiffness_v;
  Span<float> triangle_shear_stiffness;
  Span<float> bend_stiffness;
  Span<float> spring_stiffness;

  Span<float> triangle_stretch_damping_u;
  Span<float> triangle_stretch_damping_v;
  Span<float> triangle_shear_damping;
  Span<float> bend_damping;
  Span<float> spring_damping;

  Span<float3> triangle_normals;
  Span<float> triangle_area_factors;
  Span<float2x2> triangle_inverted_delta_u_matrices;
  Span<float3> triangle_wu_derivatives;
  Span<float3> triangle_wv_derivatives;

  Span<float> bend_rest_lengths;
  Span<float> spring_rest_lengths;

  Vector<int> pinned_vertices;

  /* State */
  VectorXf vertex_positions_eigen;
  VectorXf vertex_velocities_eigen;

  /* Local forces */
  Array<DeformationGradient> deformation_gradients;
  Array<GravityForceElement> gravity_force_elements;
  Array<StretchForceElementBW> stretch_force_elements;
  Array<ShearForceElementBW> shear_force_elements;
  Array<BendForceElementBW> bend_force_elements;
  Array<SpringForceElement> spring_force_elements;

  /* Global system (Eigen) */
  VectorXf vertex_forces;
  SparseMatrix<float> vertex_force_derivatives;
  SparseMatrix<float> vertex_force_velocity_derivatives;
  SparseMatrix<float> mass_matrix;
  SparseMatrix<float> identity_matrix;

  /* These vectors of Triplets are used to initialize the corresponding SparseMatrices. */
  std::vector<Triplet> force_derivative_triplets;
  std::vector<Triplet> force_velocity_derivative_triplets;

  void initialize(const Mesh &mesh,
                  const ClothBWModifierData &modifier_data,
                  const Object &cloth_object)
  {
    substeps = modifier_data.n_substeps;
    step_time = 1.0f / 24.0f; /* Currently assuming 24 fps. */
    substep_time = step_time / substeps;

    use_explicit_integration = modifier_data.use_explicit_integration;

    ap = ClothBWAttributeProvider(mesh, modifier_data, cloth_object);
    amount_of_vertices = ap.get_amount_of_vertices();
    amount_of_triangles = ap.get_amount_of_triangles();
    amount_of_bend_edges = ap.get_amount_of_bend_edges();
    amount_of_springs = ap.get_amount_of_springs();

    triangle_vertex_indices = ap.get_triangle_vertex_indices();
    bend_vertex_indices = ap.get_bend_vertex_indices();
    spring_vertex_indices = ap.get_spring_vertex_indices();

    pinned_vertices = ap.get_pinned_vertices();

    vertex_positions = ap.get_vertex_positions();
    vertex_velocities = ap.get_vertex_velocities();

    /* Setting the read-only attributes. */
    vertex_masses = ap.get_vertex_masses();

    triangle_stretch_stiffness_u = ap.get_triangle_stretch_stiffness_u();
    triangle_stretch_stiffness_v = ap.get_triangle_stretch_stiffness_v();
    triangle_shear_stiffness = ap.get_triangle_shear_stiffness();
    bend_stiffness = ap.get_bend_stiffness();
    spring_stiffness = ap.get_spring_stiffness();

    triangle_stretch_damping_u = ap.get_triangle_stretch_damping_u();
    triangle_stretch_damping_v = ap.get_triangle_stretch_damping_v();
    triangle_shear_damping = ap.get_triangle_shear_damping();
    bend_damping = ap.get_bend_damping();
    spring_damping = ap.get_spring_damping();

    triangle_area_factors = ap.get_triangle_area_factors();
    triangle_inverted_delta_u_matrices = ap.get_triangle_inverted_delta_u_matrices();
    triangle_wu_derivatives = ap.get_triangle_wu_derivatives();
    triangle_wv_derivatives = ap.get_triangle_wv_derivatives();

    bend_rest_lengths = ap.get_bend_rest_lengths();
    spring_rest_lengths = ap.get_spring_rest_lengths();

    /* Creating the local storage for the forces. */
    deformation_gradients = Array<DeformationGradient>(amount_of_triangles);
    gravity_force_elements = Array<GravityForceElement>(amount_of_vertices);
    stretch_force_elements = Array<StretchForceElementBW>(amount_of_triangles);
    shear_force_elements = Array<ShearForceElementBW>(amount_of_triangles);
    bend_force_elements = Array<BendForceElementBW>(amount_of_bend_edges);
    spring_force_elements = Array<SpringForceElement>(amount_of_springs);

    /* Creating the global matrixs and vectors for the Eigen solver. */
    system_size = 3 * amount_of_vertices;

    vertex_positions_eigen = VectorXf(system_size);
    vertex_velocities_eigen = VectorXf(system_size);
    vertex_forces = VectorXf(system_size);
    vertex_force_derivatives = SparseMatrix<float>(system_size, system_size);
    vertex_force_velocity_derivatives = SparseMatrix<float>(system_size, system_size);
    identity_matrix = SparseMatrix<float>(system_size, system_size);
    identity_matrix.setIdentity();

    initialize_mass_matrix();

    /* Currently only initial positions and velocities are copied. Later maybe these need to be
     * copied each timestep to account for outside influence. */
    for (int i : IndexRange(amount_of_vertices)) {
      for (int s : IndexRange(3)) {
        vertex_positions_eigen[3 * i + s] = vertex_positions[i][s];
        vertex_velocities_eigen[3 * i + s] = vertex_velocities[i][s];
      }
    }

    int amount_of_force_derivative_elements = 9 * 9 * stretch_force_elements.size() +
                                              9 * 9 * shear_force_elements.size() +
                                              12 * 12 * bend_force_elements.size();

    force_derivative_triplets = std::vector<Triplet>();
    force_velocity_derivative_triplets = std::vector<Triplet>();
    force_derivative_triplets.reserve(amount_of_force_derivative_elements);
    force_velocity_derivative_triplets.reserve(amount_of_force_derivative_elements);
  }

  void step()
  {
    std::cout << "step eigen" << std::endl;

    for (int substep : IndexRange(substeps)) {
      current_substep = substep;
      {
        SCOPED_TIMER("calculate_forces_and_derivatives");
        calculate_local_forces_and_derivatives();
      }

      {
        SCOPED_TIMER("accumulate_local_forces_to_global");
        accumulate_local_forces_to_global();
      }
      // if (use_explicit_integration) {
      // integrate_explicit_forward_euler();
      // }
      // else {
      {
        SCOPED_TIMER("integrate_implicit_backward_euler");
        integrate_implicit_backward_euler();
      }
      // }

      std::cout << std::endl;
    }
  }

  void initialize_mass_matrix()
  {
    mass_matrix = SparseMatrix<float>(system_size, system_size);

    Array<Triplet> triplets = Array<Triplet>(system_size);
    for (int i : IndexRange(amount_of_vertices)) {
      for (int s : IndexRange(3)) {
        int index = 3 * i + s;
        triplets[index] = Triplet(index, index, vertex_masses[i]);
      }
    }
    mass_matrix.setFromTriplets(triplets.begin(), triplets.end());
  }

  void calculate_local_forces_and_derivatives()
  {
    for (int i : IndexRange(amount_of_vertices)) {
      float vertex_mass = vertex_masses[i];
      gravity_force_elements[i].calculate(vertex_mass, standard_gravity);
    }

    for (int ti : IndexRange(amount_of_triangles)) {
      int3 vertex_indices = triangle_vertex_indices[ti];
      float3 x0 = vertex_positions[vertex_indices[0]];
      float3 x1 = vertex_positions[vertex_indices[1]];
      float3 x2 = vertex_positions[vertex_indices[2]];
      float2x2 delta_u = triangle_inverted_delta_u_matrices[ti];
      deformation_gradients[ti].calculate(x0, x1, x2, delta_u);
      DeformationGradient F = deformation_gradients[ti];

      float area_factor = triangle_area_factors[ti];
      float ku = triangle_stretch_stiffness_u[ti];
      float kv = triangle_stretch_stiffness_v[ti];
      float3 dwu_dx = triangle_wu_derivatives[ti];
      float3 dwv_dx = triangle_wv_derivatives[ti];
      float kdu = triangle_stretch_damping_u[ti];
      float kdv = triangle_stretch_damping_v[ti];

      float3 v0 = vertex_velocities[vertex_indices[0]];
      float3 v1 = vertex_velocities[vertex_indices[1]];
      float3 v2 = vertex_velocities[vertex_indices[2]];
      array<float3, 3> velocities = {v0, v1, v2};

      stretch_force_elements[ti].calculate(
          ku, kv, area_factor, F, dwu_dx, dwv_dx, kdu, kdv, velocities);

      float k_shear = triangle_shear_stiffness[ti];
      float kd_shear = triangle_shear_damping[ti];
      shear_force_elements[ti].calculate(
          k_shear, area_factor, F, dwu_dx, dwv_dx, kd_shear, velocities);
    }

    for (int bi : IndexRange(amount_of_bend_edges)) {
      int4 vertex_indices = bend_vertex_indices[bi];
      float3 x0 = vertex_positions[vertex_indices[0]];
      float3 x1 = vertex_positions[vertex_indices[1]];
      float3 x2 = vertex_positions[vertex_indices[2]];
      float3 x3 = vertex_positions[vertex_indices[3]];
      float k = bend_stiffness[bi];
      float kd = bend_damping[bi];

      float3 v0 = vertex_velocities[vertex_indices[0]];
      float3 v1 = vertex_velocities[vertex_indices[1]];
      float3 v2 = vertex_velocities[vertex_indices[2]];
      float3 v3 = vertex_velocities[vertex_indices[3]];
      array<float3, 4> velocities = {v0, v1, v2, v3};

      bend_force_elements[bi].calculate(k, x0, x1, x2, x3, kd, velocities);
    }

    for (int si : IndexRange(amount_of_springs)) {
      int2 vertex_indices = spring_vertex_indices[si];
      float3 x0 = vertex_positions[vertex_indices[0]];
      float3 x1 = vertex_positions[vertex_indices[1]];
      float k = spring_stiffness[si];
      float rest_length = spring_rest_lengths[si];

      float kd = spring_damping[si];
      float3 v0 = vertex_velocities[vertex_indices[0]];
      float3 v1 = vertex_velocities[vertex_indices[1]];

      spring_force_elements[si].calculate(k, rest_length, x0, x1, kd, v0, v1);
    }
  }

  void integrate_explicit_forward_euler()
  {
    VectorXf vertex_accelerations(3 * amount_of_vertices);

    for (int i : IndexRange(amount_of_vertices)) {
      vertex_accelerations.segment(3 * i, 3) = 1.0 / vertex_masses[i] *
                                               vertex_forces.segment(3 * i, 3);
    }

    for (int i : pinned_vertices) {
      vertex_accelerations.segment(3 * i, 3).setZero();
    }

    vertex_velocities_eigen += substep_time * vertex_accelerations;
    vertex_positions_eigen += substep_time * vertex_velocities_eigen;

    copy_eigen_to_blender();
  }

  SparseMatrix<float> create_S_matrix()
  {
    SparseMatrix<float> S = SparseMatrix<float>(system_size, system_size);
    Array<Triplet> S_triplets = Array<Triplet>(system_size);
    for (int i : IndexRange(amount_of_vertices)) {
      for (int s : IndexRange(3)) {
        int index = 3 * i + s;
        float value = 1.0f;
        if (std::find(pinned_vertices.begin(), pinned_vertices.end(), i) !=
            pinned_vertices.end()) {
          value = 0.0f;
        }
        S_triplets[index] = Triplet(index, index, value);
      }
    }
    S.setFromTriplets(S_triplets.begin(), S_triplets.end());
    return S;
  }

  void integrate_implicit_backward_euler()
  {
    /* First we build the A and b of the linear system Ax = b. Then we use eigen to solve it. */
    float h = substep_time;
    SparseMatrix<float> &dfdx = vertex_force_derivatives;
    SparseMatrix<float> &dfdv = vertex_force_velocity_derivatives;

    SparseMatrix<float> &M = mass_matrix;
    SparseMatrix<float> &I = identity_matrix;

    VectorXf &v0 = vertex_velocities_eigen;
    VectorXf &f0 = vertex_forces;

    /* Creating the b vector, the right hand side of the linear system. */
    VectorXf b(system_size);  // TODO maybe preallocate this memory if needed.
    b = h * (f0 + h * (dfdx * v0));

    VectorXf z(system_size);
    z.fill(0.0f);

    /* Creating the A matrix the left hand side of the linear systems. */
    SparseMatrix<float> A = SparseMatrix<float>(system_size, system_size);
    A = M - h * dfdv - (h * h) * dfdx;

    /* PPCG */
    SparseMatrix<float> S = create_S_matrix();
    SparseMatrix<float> ST = SparseMatrix<float>(S.transpose());
    SparseMatrix<float> LHS = (S * A * ST) + I - S;
    VectorXf c = b - A * z;
    VectorXf rhs = S * c;

    ConjugateGradient<SparseMatrix<float>, Lower | Upper> solver;

    // Solver that doesn't require SPD-ness, could be used as backup in case of non-convergance.
    // Eigen::BiCGSTAB<SparseMatrix<float>> solver;

    VectorXf y;
    {
      SCOPED_TIMER("solve");
      solver.setTolerance(0.01);
      solver.setMaxIterations(50);
      solver.analyzePattern(LHS);
      y = solver.compute(LHS).solve(rhs);
    }

    VectorXf x = y + z;

    vertex_velocities_eigen += x;
    vertex_positions_eigen += h * vertex_velocities_eigen;
    copy_eigen_to_blender();
  }

  // TODO: make a cleaner interface to get the vertex positions & velocities.
  void copy_eigen_to_blender()
  {
    for (int i : IndexRange(amount_of_vertices)) {
      vertex_positions[i][0] = vertex_positions_eigen[3 * i];
      vertex_positions[i][1] = vertex_positions_eigen[3 * i + 1];
      vertex_positions[i][2] = vertex_positions_eigen[3 * i + 2];

      vertex_velocities[i][0] = vertex_velocities_eigen[3 * i];
      vertex_velocities[i][1] = vertex_velocities_eigen[3 * i + 1];
      vertex_velocities[i][2] = vertex_velocities_eigen[3 * i + 2];
    }
  }

  /* Code for accumulating the locally calculated force elements into the global Eigen system.*/
  void insert_float3x3_as_triplets(const float3x3 &matrix,
                                   std::vector<Triplet> &triplets,
                                   int i,
                                   int j)
  {
    for (int s : IndexRange(3)) {
      for (int t : IndexRange(3)) {
        Triplet triplet = Triplet(3 * i + s, 3 * j + t, matrix.values[t][s]);
        triplets.push_back(triplet);
      }
    }
  }

  void accumulate_particle_based_force(float3 force, int vertex_index)
  {
    int i = vertex_index;
    vertex_forces[3 * i] += force.x;
    vertex_forces[3 * i + 1] += force.y;
    vertex_forces[3 * i + 2] += force.z;
  }

  template<int n_vertices>
  void accumulate_element(
      const array<float3, n_vertices> &forces,
      const array<array<float3x3, n_vertices>, n_vertices> &force_derivatives,
      const array<array<float3x3, n_vertices>, n_vertices> &force_velocity_derivatives,
      const array<int, n_vertices> &vertex_indices)
  {
    for (int m : IndexRange(n_vertices)) {
      int i = vertex_indices[m];
      float3 force = forces[m];
      vertex_forces[3 * i] += force.x;
      vertex_forces[3 * i + 1] += force.y;
      vertex_forces[3 * i + 2] += force.z;

      for (int n : IndexRange(n_vertices)) {
        int j = vertex_indices[n];
        insert_float3x3_as_triplets(force_derivatives[m][n], force_derivative_triplets, i, j);
        insert_float3x3_as_triplets(
            force_velocity_derivatives[m][n], force_velocity_derivative_triplets, i, j);
      }
    }
  }

  void accumulate_local_forces_to_global()
  {
    vertex_forces.fill(0.0f);

    for (int i : IndexRange(gravity_force_elements.size())) {
      GravityForceElement element = gravity_force_elements[i];
      accumulate_particle_based_force(element.force, i);
    }

    force_derivative_triplets.clear();
    force_velocity_derivative_triplets.clear();

    for (int ti : IndexRange(stretch_force_elements.size())) {
      StretchForceElementBW element = stretch_force_elements[ti];
      accumulate_element<3>(element.forces,
                            element.force_derivatives,
                            element.force_velocity_derivatives,
                            triangle_vertex_indices[ti]);
    }

    for (int ti : IndexRange(shear_force_elements.size())) {
      ShearForceElementBW element = shear_force_elements[ti];
      accumulate_element<3>(element.forces,
                            element.force_derivatives,
                            element.force_velocity_derivatives,
                            triangle_vertex_indices[ti]);
    }

    for (int bi : IndexRange(bend_force_elements.size())) {
      BendForceElementBW element = bend_force_elements[bi];
      accumulate_element<4>(element.forces,
                            element.force_derivatives,
                            element.force_velocity_derivatives,
                            bend_vertex_indices[bi]);
    }

    for (int si : IndexRange(spring_force_elements.size())) {
      SpringForceElement element = spring_force_elements[si];
      accumulate_element<2>(element.forces,
                            element.force_derivatives,
                            element.force_velocity_derivatives,
                            spring_vertex_indices[si]);
    }

    vertex_force_derivatives.setFromTriplets(force_derivative_triplets.begin(),
                                             force_derivative_triplets.end());

    vertex_force_velocity_derivatives.setFromTriplets(force_velocity_derivative_triplets.begin(),
                                                      force_velocity_derivative_triplets.end());
  }
};