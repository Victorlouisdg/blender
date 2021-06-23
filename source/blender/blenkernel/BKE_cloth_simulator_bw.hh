#include <iostream>

#include "BKE_cloth_simulator_bw_attribute_provider.hh"
#include "BKE_cloth_simulator_bw_force_elements.hh"

#include "BLI_index_range.hh"
#include "BLI_span.hh"
#include "BLI_vector.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

using blender::IndexRange;
using blender::MutableSpan;
using blender::Vector;

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

typedef Eigen::Triplet<double> Triplet;

using Eigen::COLAMDOrdering;
using Eigen::SparseLU;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

class ClothSimulatorBW {
 public:
  int substeps;
  double step_time; /* Currently assuming 30 fps. */
  double substep_time;

  bool use_explicit_integration;

  int amount_of_vertices;
  int amount_of_triangles;
  int system_size;

  double standard_gravity = -9.81; /* in m/s^2 */

  /* Indices */
  Span<int2> spring_vertex_indices;
  Span<int3> triangle_vertex_indices;
  Span<int4> bending_vertex_indices;

  /* Read-Only Attributes, eventually provided by e.g. Geometry Nodes */
  Span<float3> vertex_positions;
  Span<float3> vertex_velocities;

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
  VectorXd vertex_positions_eigen;
  VectorXd vertex_velocities_eigen;

  /* Local forces */
  Array<StretchForceElementBW> stretch_force_elements;

  /* Global system (Eigen) */
  VectorXd vertex_forces;
  SparseMatrix<double> vertex_force_derivatives;

  void initialize(const Mesh &mesh, const ClothBWModifierData &modifier_data)
  {
    substeps = modifier_data.n_substeps;
    step_time = 1.0 / 30.0; /* Currently assuming 30 fps. */
    substep_time = step_time / substeps;

    use_explicit_integration = modifier_data.use_explicit_integration;

    ClothBWAttributeProvider ap = ClothBWAttributeProvider(mesh, modifier_data);
    amount_of_vertices = ap.get_amount_of_vertices();
    amount_of_triangles = ap.get_amount_of_triangles();

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
    stretch_force_elements = Array<StretchForceElementBW>(amount_of_triangles);

    /* Creating the global matrixs and vectors for the Eigen solver. */
    system_size = 3 * amount_of_vertices;

    vertex_positions_eigen = VectorXd(system_size);
    vertex_velocities_eigen = VectorXd(system_size);
    vertex_forces = VectorXd(system_size);
    vertex_force_derivatives = SparseMatrix<double>(system_size, system_size);
  }

  void step()
  {
    std::cout << "step eigen" << std::endl;
    for (int substep : IndexRange(substeps)) {
      calculate_forces_and_derivatives();
      aggregate_forces_local_to_global();

      // if (use_explicit_integration) {
      //   integrate_explicit_forward_euler();
      // }
      // else {
      //   integrate_implicit_backward_euler();
      // }
    }
  }

  // void pin(int vertex_index)
  // {
  //   pinned_vertices.append(vertex_index);
  //   // solver.setConstraint(vertex_index, float3(0.0f), float3x3(0.0f));
  // }

  void calculate_forces_and_derivatives()
  {
    for (int vertex_index : IndexRange(amount_of_vertices)) {
      // calculate_gravity(vertex_index);
    }

    for (int ti : IndexRange(amount_of_triangles)) {
      auto [wu, wv] = calculate_w_uv(ti);
      float area_factor = triangle_area_factors[ti];
      float ku = triangle_stretch_stiffness_u[ti];
      float kv = triangle_stretch_stiffness_v[ti];
      float3 dwu_dx = triangle_wu_derivatives[ti];
      float3 dwv_dx = triangle_wu_derivatives[ti];

      stretch_force_elements[ti].calculate(ku, kv, area_factor, wu, wv, dwu_dx, dwv_dx);
    }
  }

  /* TODO think where this belongs, maybe a DeformationGradient3x2 class? */
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
    /* This copying from ForceElements into the Eigen datastructures is a bit awkward for two
     * reason: 1) The float3x3 values are stored in column major format. 2) The eigen sparse matrix
     * is not 3x3 blocked.
     */
    int total_force_elements = 9 * 9 * stretch_force_elements.size();
    Array<Triplet> triplets = Array<Triplet>(total_force_elements);
    int triplet_index = 0;

    for (int ti : IndexRange(stretch_force_elements.size())) {
      int3 vertex_indices = triangle_vertex_indices[ti];

      for (int m : IndexRange(3)) {
        int i = vertex_indices[m];
        float3 force = stretch_force_elements[ti].forces[m];
        vertex_forces[3 * i] = force.x;
        vertex_forces[3 * i + 1] = force.y;
        vertex_forces[3 * i + 2] = force.z;

        for (int n : IndexRange(3)) {
          int j = vertex_indices[n];

          float3x3 force_derivative = stretch_force_elements[ti].force_derivatives[m][n];
          insert_triplets_from_float3x3(triplets, triplet_index, force_derivative, i, j);
        }
      }
    }
  }

  // void calculate_gravity(int vertex_index)
  // {
  //   Vector3 gravity = Vector3();
  //   gravity.setZero();
  //   gravity[2] = vertex_masses[vertex_index] * standard_gravity;  // F_grav = m * g
  //   vertex_forces.segment(3 * vertex_index, 3) += gravity;
  // }

  // void integrate_explicit_forward_euler()
  // {
  //   VectorXd vertex_accelerations(3 * n_vertices);

  //   for (int i : IndexRange(n_vertices)) {
  //     vertex_accelerations.segment(3 * i, 3) = 1.0 / vertex_masses[i] *
  //                                              vertex_forces.segment(3 * i, 3);

  //     std::cout << vertex_accelerations.segment(3 * i, 3) << std::endl;
  //   }

  //   for (int i : pinned_vertices) {
  //     vertex_accelerations.segment(3 * i, 3).setZero();
  //   }

  //   for (int i : IndexRange(n_vertices)) {
  //     vertex_velocities.segment(3 * i, 3) += vertex_accelerations.segment(3 * i, 3) *
  //     substep_time; vertex_positions.segment(3 * i, 3) += vertex_velocities.segment(3 * i, 3) *
  //     substep_time;
  //   }
  // }

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

  // void stretchHessian(const Matrix3x2 &F, Matrix6x6 &H) const
  // {
  //   H.setZero();
  //   const Vector2 u(1.0, 0.0);
  //   const Vector2 v(0.0, 1.0);
  //   const double I5u = (F * u).transpose() * (F * u);
  //   const double I5v = (F * v).transpose() * (F * v);
  //   const double invSqrtI5u = 1.0 / sqrt(I5u);
  //   const double invSqrtI5v = 1.0 / sqrt(I5v);

  //   // set the block diagonals, build the rank-three
  //   // subspace with all-(1 / invSqrtI5) eigenvalues
  //   H(0, 0) = H(1, 1) = H(2, 2) = std::max((1.0 - invSqrtI5u), 0.0);
  //   H(3, 3) = H(4, 4) = H(5, 5) = std::max((1.0 - invSqrtI5v), 0.0);

  //   // modify the upper block diagonal, bump the single
  //   // outer-product eigenvalue back to just 1, unless it
  //   // was clamped, then just set it directly to 1
  //   const Vector3 fu = F.col(0).normalized();
  //   const double uCoeff = (1.0 - invSqrtI5u >= 0.0) ? invSqrtI5u : 1.0;
  //   H.block<3, 3>(0, 0) += uCoeff * (fu * fu.transpose());

  //   // modify the lower block diagonal similarly
  //   const Vector3 fv = F.col(1).normalized();
  //   const double vCoeff = (1.0 - invSqrtI5v >= 0.0) ? invSqrtI5v : 1.0;
  //   H.block<3, 3>(3, 3) += vCoeff * (fv * fv.transpose());

  //   // the leading 2 is absorbed by the mu / 2 coefficient
  //   H *= _mu;
  // }
};