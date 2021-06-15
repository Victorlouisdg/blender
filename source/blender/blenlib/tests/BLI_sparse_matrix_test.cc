/* Apache License, Version 2.0 */

#include "testing/testing.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

namespace blender::tests {

using Eigen::COLAMDOrdering;
using Eigen::SparseLU;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

TEST(sparse_matrix, SolveLU)
{
  int n = 100;
  SparseMatrix<double> A(n, n);
  VectorXd x(n), b(n);

  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  x = solver.solve(b);
}

}  // namespace blender::tests
