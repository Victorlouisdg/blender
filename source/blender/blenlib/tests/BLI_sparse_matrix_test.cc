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

}  // namespace blender::tests
