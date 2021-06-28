#include <iostream>

#include "BKE_cloth_simulator_bw_attribute_provider.hh"
#include "BKE_cloth_simulator_bw_force_elements.hh"

#include "BLI_index_range.hh"
#include "BLI_span.hh"
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

typedef Eigen::Triplet<double> Triplet;

using Eigen::COLAMDOrdering;
using Eigen::SparseLU;
using Eigen::SparseMatrix;
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
  int system_size;

  float standard_gravity = -9.81f; /* in m/s^2 */

  ClothBWAttributeProvider ap;

  /* Indices */
  Span<int2> spring_vertex_indices;
  Span<int3> triangle_vertex_indices;
  Span<int4> bending_vertex_indices;

  MutableSpan<float3> vertex_positions;
  MutableSpan<float3> vertex_velocities;
  /* Read-Only Attributes, eventually provided by e.g. Geometry Nodes */

  Span<float> vertex_masses;
  Span<float> triangle_stretch_stiffness_u;
  Span<float> triangle_stretch_stiffness_v;
  Span<float> triangle_shear_stiffness;
  Span<float> bending_stiffness;
  Span<float> bending_rest_lengths;
  Span<float3> triangle_normals;
  Span<float> triangle_area_factors;
  Span<float2x2> triangle_inverted_delta_u_matrices;
  Span<float3> triangle_wu_derivatives;
  Span<float3> triangle_wv_derivatives;

  Vector<int> pinned_vertices;

  /* State */
  VectorXf vertex_positions_eigen;
  VectorXf vertex_velocities_eigen;

  /* Local forces */
  Array<DeformationGradient> deformation_gradients;
  Array<GravityForceElement> gravity_force_elements;
  Array<StretchForceElementBW> stretch_force_elements;

  /* Global system (Eigen) */
  VectorXf vertex_forces;
  SparseMatrix<double> vertex_force_derivatives;

  void initialize(const Mesh &mesh,
                  const ClothBWModifierData &modifier_data,
                  const Object &cloth_object)
  {
    substeps = modifier_data.n_substeps;
    step_time = 1.0f / 30.0f; /* Currently assuming 30 fps. */
    substep_time = step_time / substeps;

    use_explicit_integration = modifier_data.use_explicit_integration;

    ap = ClothBWAttributeProvider(mesh, modifier_data, cloth_object);
    amount_of_vertices = ap.get_amount_of_vertices();
    amount_of_triangles = ap.get_amount_of_triangles();

    triangle_vertex_indices = ap.get_triangle_vertex_indices();

    pinned_vertices = ap.get_pinned_vertices();

    /* Setting the read-only attributes. */
    vertex_positions = ap.get_vertex_positions();
    vertex_velocities = ap.get_vertex_velocities();
    vertex_masses = ap.get_vertex_masses();
    triangle_stretch_stiffness_u = ap.get_triangle_stretch_stiffness_u();
    triangle_stretch_stiffness_v = ap.get_triangle_stretch_stiffness_v();
    triangle_area_factors = ap.get_triangle_area_factors();
    triangle_inverted_delta_u_matrices = ap.get_triangle_inverted_delta_u_matrices();
    triangle_wu_derivatives = ap.get_triangle_wu_derivatives();
    triangle_wv_derivatives = ap.get_triangle_wv_derivatives();

    /* Creating the local storage for the forces. */
    deformation_gradients = Array<DeformationGradient>(amount_of_triangles);
    gravity_force_elements = Array<GravityForceElement>(amount_of_vertices);
    stretch_force_elements = Array<StretchForceElementBW>(amount_of_triangles);

    /* Creating the global matrixs and vectors for the Eigen solver. */
    system_size = 3 * amount_of_vertices;

    vertex_positions_eigen = VectorXf(system_size);
    vertex_velocities_eigen = VectorXf(system_size);
    vertex_forces = VectorXf(system_size);
    vertex_force_derivatives = SparseMatrix<double>(system_size, system_size);

    /* Currently only initial positions and velocities are copied. Later maybe these need to be
     * copied each timestep to account for outside influence. */
    for (int i : IndexRange(amount_of_vertices)) {
      for (int s : IndexRange(3)) {
        vertex_positions_eigen[3 * i + s] = vertex_positions[i][s];
        vertex_velocities_eigen[3 * i + s] = vertex_velocities[i][s];
      }
    }
  }

  void step()
  {
    std::cout << "step eigen" << std::endl;
    for (int substep : IndexRange(substeps)) {
      current_substep = substep;
      calculate_forces_and_derivatives();
      aggregate_forces_local_to_global();

      // if (use_explicit_integration) {
      integrate_explicit_forward_euler();
      // }
      // else {
      //   integrate_implicit_backward_euler();
      // }
    }
  }

  void calculate_forces_and_derivatives()
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

      stretch_force_elements[ti].calculate(ku, kv, area_factor, F, dwu_dx, dwv_dx);

      // std::array<float3, 3> &forces = stretch_force_elements[ti].forces;
      // const float red[4] = {1.0f, 0.0f, 0.0f, 1.0f};
      // const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};
      // const float blue[4] = {0.0f, 0.0f, 1.0f, 1.0f};
      // if (current_substep == substeps - 1) {
      //   DRW_debug_line_v3v3(x0, x0 + forces[0], red);    // zero?
      //   DRW_debug_line_v3v3(x1, x1 + forces[1], green);  // flipped?
      //   DRW_debug_line_v3v3(x2, x2 + forces[2], blue);
      // }
    }
  }

  void insert_triplets_from_float3x3(
      MutableSpan<Triplet> triplets, int &triplet_index, float3x3 force_derivative, int i, int j)
  {
    for (int s : IndexRange(3)) {
      for (int t : IndexRange(3)) {
        /* Note the inversion of indices s and t, this is because blender
         * matrices are stored in column major format. */
        Triplet triplet = Triplet(3 * i + s, 3 * j + t, force_derivative.values[t][s]);
        triplets[triplet_index] = triplet;
        triplet_index++;
      }
    }
  }

  void aggregate_forces_local_to_global()
  {
    vertex_forces.fill(0.0f);

    for (int i : IndexRange(amount_of_vertices)) {
      float3 force = gravity_force_elements[i].force;
      /* x and y are always zero so maybe even just remove this? */
      vertex_forces[3 * i] += force.x;
      vertex_forces[3 * i + 1] += force.y;
      vertex_forces[3 * i + 2] += force.z;
    }

    /* This copying from ForceElements into the Eigen datastructures is a bit awkward for two
     * reason: 1) The float3x3 values are stored in column major format. 2) The eigen sparse matrix
     * is not 3x3 blocked.
     */
    int total_force_derivatives = 9 * 9 * stretch_force_elements.size();
    Array<Triplet> triplets = Array<Triplet>(total_force_derivatives);
    int triplet_index = 0;

    for (int ti : IndexRange(stretch_force_elements.size())) {
      int3 vertex_indices = triangle_vertex_indices[ti];

      // std::cout << "Forces on the triangle" << std::endl;

      for (int m : IndexRange(3)) {
        int i = vertex_indices[m];
        float3 force = stretch_force_elements[ti].forces[m];
        vertex_forces[3 * i] += force.x;
        vertex_forces[3 * i + 1] += force.y;
        vertex_forces[3 * i + 2] += force.z;

        for (int n : IndexRange(3)) {
          int j = vertex_indices[n];

          float3x3 force_derivative = stretch_force_elements[ti].force_derivatives[m][n];
          insert_triplets_from_float3x3(triplets, triplet_index, force_derivative, i, j);
        }
      }
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

    for (int i : IndexRange(amount_of_vertices)) {
      vertex_velocities_eigen.segment(3 * i, 3) += vertex_accelerations.segment(3 * i, 3) *
                                                   substep_time;
      vertex_positions_eigen.segment(3 * i, 3) += vertex_velocities_eigen.segment(3 * i, 3) *
                                                  substep_time;
    }

    for (int i : IndexRange(amount_of_vertices)) {
      // Make a cleaner interface to get the vertex positions.
      vertex_positions[i][0] = vertex_positions_eigen[3 * i];
      vertex_positions[i][1] = vertex_positions_eigen[3 * i + 1];
      vertex_positions[i][2] = vertex_positions_eigen[3 * i + 2];

      vertex_velocities[i][0] = vertex_velocities_eigen[3 * i];
      vertex_velocities[i][1] = vertex_velocities_eigen[3 * i + 1];
      vertex_velocities[i][2] = vertex_velocities_eigen[3 * i + 2];
    }
  }

  // void integrate_implicit_backward_euler()
  // {
  //   /*First we build the A and b of the linear system Ax = b. Then we use eigen to solve it. */
  //   double h = substep_time;
  //   SparseMatrix<double> &dfdx = vertex_force_derivatives;

  //   VectorXd &v0 = vertex_velocities;
  //   VectorXd &f0 = vertex_forces;

  //   /* Creating the b vector, the right hand side of the linear system. */
  //   VectorXd b(3 * n_vertices);  // TODO maybe preallocate this memory if needed.
  //   b = dfdx * v0;
  //   b *= h;
  //   b += f0;
  //   b *= h;

  //   /* Creating the A matrix the left hand side of the linear systems. */
  //   SparseMatrix<double> &A = dfdx;
  //   A *= (h * h);

  //   // TODO A = M (diagonal) - A;
  //   for (int i : IndexRange(n_vertices)) {
  //     A.coeffRef(i, i) = vertex_masses[i] - A.coeffRef(i, i);
  //   }
  // }
};