#include <iostream>

#include "BKE_cloth_bw_general.hh"
#include "BKE_cloth_simulator_bw_attribute_provider.hh"

#include "BLI_index_range.hh"
#include "BLI_vector.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

using blender::IndexRange;
using blender::Vector;

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

using Eigen::COLAMDOrdering;
using Eigen::SparseLU;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

class ClothSimulatorBWEigen {
 public:
  int substeps;
  double step_time; /* Currently assuming 30 fps. */
  double substep_time;

  bool use_explicit_integration;

  int amount_of_vertices;
  int amount_of_triangles;
  int system_size;

  double standard_gravity = -9.81; /* in m/s^2 */

  /* State */
  VectorXd &vertex_positions;
  VectorXd &vertex_velocities;

  /* Attributes owned by the simulator. */
  VectorXd vertex_forces;
  SparseMatrix<double> vertex_force_derivatives;

  /* Read-Only Attributes, eventually provided by e.g. Geometry Nodes */
  VectorXd &vertex_masses;

  Vector<int> pinned_vertices;

  void initialize(const Mesh &mesh, int substeps, bool use_explicit_integration)
  {
    this->substeps = substeps;
    step_time = 1.0 / 30.0; /* Currently assuming 30 fps. */
    substep_time = step_time / substeps;

    this->use_explicit_integration = use_explicit_integration;

    ClothBWAttributeProvider AP = ClothBWAttributeProvider(mesh);
    amount_of_vertices = AP.get_amount_of_vertices();
    system_size = 3 * amount_of_vertices;

    vertex_forces = VectorXd(system_size);
    vertex_force_derivatives = SparseMatrix<double>(system_size, system_size);

    // vertex_positions = AP.get_vertex_positions();
    // vertex_velocities = AP.get_vertex_velocities();
    // vertex_masses = AP.get_vertex_masses();
  }

  // void step()
  // {
  //   std::cout << "step eigen" << std::endl;
  //   for (int substep : IndexRange(n_substeps)) {
  //     reset_forces_and_derivatives();
  //     calculate_forces_and_derivatives();

  //     if (use_explicit_integration) {
  //       integrate_explicit_forward_euler();
  //     }
  //     else {
  //       integrate_implicit_backward_euler();
  //     }
  //   }
  // }

  // void pin(int vertex_index)
  // {
  //   pinned_vertices.append(vertex_index);
  //   // solver.setConstraint(vertex_index, float3(0.0f), float3x3(0.0f));
  // }

  // void calculate_forces_and_derivatives()
  // {
  //   for (int vertex_index : IndexRange(n_vertices)) {
  //     calculate_gravity(vertex_index);
  //   }

  //   // for (int triangle_index : IndexRange(n_triangles)) {
  //   //   auto [wu, wv] = calculate_w_uv(triangle_index);
  //   //   calculate_stretch(triangle_index, wu, wv);
  //   // }
  // }

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

  // void reset_forces_and_derivatives()
  // {
  //   vertex_forces.fill(Vector3(0.0));
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