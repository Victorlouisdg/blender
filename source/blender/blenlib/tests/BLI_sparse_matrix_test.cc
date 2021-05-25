/* Apache License, Version 2.0 */

#include "BLI_float3x3.hh"
#include "BLI_sparse_matrix.hh"
#include "testing/testing.h"
#include <numeric>

namespace blender::tests {

TEST(sparse_matrix, Insert)
{
  int n = 5;
  SparseMatrix matrix = SparseMatrix(n);
  Array<float3> vector = Array<float3>(n, float3(1.0f));
  float array[3][3] = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  float3x3 block = float3x3(array);
  matrix.insert(2, 2, block);
}

TEST(sparse_matrix, MultiplyIdentity)
{
  int n = 5;
  SparseMatrix matrix = SparseMatrix(n);
  for (int i : IndexRange(n)) {
    matrix.insert(i, i, float3x3::identity());
  }

  Array<float3> vector = Array<float3>(n);
  for (int i : IndexRange(n)) {
    float v = 3.0f * i;
    vector[i] = float3({v, v + 1, v + 2});
  }

  Array<float3> result = Array<float3>(n);
  matrix.multiply(vector, result);

  for (int i : IndexRange(n)) {
    for (int s : IndexRange(3)) {
      EXPECT_FLOAT_EQ(vector[i][s], result[i][s]);
    }
  }
}

TEST(sparse_matrix, DistrubutiveProperty)
{
  srand(0);
  int n = 5;
  SparseMatrix A = SparseMatrix(n);
  for (int i : IndexRange(n)) {

    float3x3 matrix = float3x3();
    for (int s : IndexRange(3)) {
      for (int t : IndexRange(3)) {
        matrix.values[s][t] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      }
    }

    A.insert(i, i, matrix);
  }

  Array<float3> b = Array<float3>(n);
  for (int i : IndexRange(n)) {
    float v = 3.0f * i;
    b[i] = float3({v, v + 1, v + 2});
  }

  Array<float3> c = Array<float3>(n, float3(2.0f));

  Array<float3> Ab = Array<float3>(n);
  Array<float3> Ac = Array<float3>(n);
  Array<float3> AbAc = Array<float3>(n);
  Array<float3> Abc = Array<float3>(n);
  Array<float3> bc = Array<float3>(n);

  A.multiply(b, Ab);
  A.multiply(c, Ac);

  for (int i : IndexRange(n)) {
    bc[i] = b[i] + c[i];
    AbAc[i] = Ab[i] + Ac[i];
  }

  A.multiply(bc, Abc);

  for (int i : IndexRange(n)) {
    for (int s : IndexRange(3)) {
      EXPECT_FLOAT_EQ(Abc[i][s], AbAc[i][s]);
    }
  }
}

TEST(sparse_matrix, SolveGaussSeidelSimple3x3)
{
  /* Simple 3x3 system. */
  SparseMatrix A = SparseMatrix(1);
  float array[3][3] = {{2, -2, 2}, {4, -3, 2}, {2, 1, -3}};
  A.insert(0, 0, float3x3(array));
  Array<float3> b = Array<float3>(1, float3(16, -5, -3));
  Array<float3> x = Array<float3>(1, float3(0.0f));

  solve_gauss_seidel(A, b, x);

  EXPECT_FLOAT_EQ(x[0].x, 1.0f);
  EXPECT_FLOAT_EQ(x[0].y, 2.0f);
  EXPECT_FLOAT_EQ(x[0].z, 3.0f);
}

TEST(sparse_matrix, SolvePCGSimple3x3)
{
  /* Simple definite 3x3 system. */
  SparseMatrix A = SparseMatrix(1);
  float array[3][3] = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
  A.insert(0, 0, float3x3(array));
  Array<float3> b = Array<float3>(1, float3(1, 1, 1));
  Array<float3> x = Array<float3>(1, float3(0.0f));

  solve_filtered_pcg(A, b, x);

  EXPECT_FLOAT_EQ(x[0].x, 1.5f);
  EXPECT_FLOAT_EQ(x[0].y, 2.0f);
  EXPECT_FLOAT_EQ(x[0].z, 1.5f);
}

}  // namespace blender::tests
