/* Apache License, Version 2.0 */

#include "BLI_float3x3.hh"
#include "BLI_sparse_matrix.hh"
#include "testing/testing.h"

namespace blender::tests {

TEST(sparse_matrix, Insertion)
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

  Array<float3> result = Array<float3>(n);

  matrix.multiply(vector, result);

  for (int i : IndexRange(n)) {
    std::cout << result[i] << std::endl;
  }
}

}  // namespace blender::tests
