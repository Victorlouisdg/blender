#include "BLI_array.hh"
#include "BLI_float3x3.hh"
#include <map>

/* This sparse matrix class was specifically designed for the Baraff-Witkin Cloth Simulator.
 * The matrices used there have some properties we can optimize for:
 *
 * 1) Each row represents a cloth particle. Each particle participates in a few conditions on its
 * neighboring edges/faces. Thus each row of the Sparse matrix will contain at least some entries.
 * For this reason we use a (dense) Array for the rows of the matrix.
 *
 * 2) Besides the conditions between a particle and its neighbors, particles can also interact with
 * non-neighboring particles through collisions. However, in general a particle will collide only
 * with one (or max a few) faces, so the amount of entries per row will be relatively low. This
 * might make if possible to statically preallocate the storage for each row, instead of
 * dynamiacally allocating for each added 3x3 Block.
 *
 * 3) For the Cloth simulator is it not necessary to be able to remove previously inserted blocks,
 * just insertion of new blocks and clearing the entire matrix each timestep are sufficient.
 *
 * 4) The main operations that have to be fast are insertion and matrix-vector multiplication
 * (which required iteration.) So for now I've used an std::map for the storage of the column
 * elements within each row. A potential problem with std::map is that the storage of the values is
 * not local/contiguous, possibly making iteration slow due to cache misses. This could be fixed by
 * a "copy/compress" pass after the matrix assembly is done and before the solving starts.
 *
 * Currently I've written the code specifically for float3x3 elements, however this could be
 * templated later when anyone needs it.
 */

namespace blender {

class SparseMatrix {
 public:
  int n_rows;
  Array<std::map<int, float3x3>> *rows;

  SparseMatrix() = default;

  SparseMatrix(int n_rows)
  {
    std::cout << "SparseMatrix " << std::endl;
    this->n_rows = n_rows;
    rows = new Array<std::map<int, float3x3>>(n_rows);
  }

  ~SparseMatrix()
  {
    std::cout << "~SparseMatrix " << std::endl;
    delete rows;
  }

  void insert(int row, int col, float3x3 value)
  {
    (*rows)[row][col] = value;
  }

  void multiply(const Span<float3> &vector, const MutableSpan<float3> &result)
  {
    std::cout << n_rows << std::endl;
    for (int i : IndexRange(n_rows)) {
      std::cout << "i " << i << std::endl;
      result[i] = float3(0.0f);
      for (std::pair<int, float3x3> element : (*rows)[i]) {
        int j = element.first;
        std::cout << "j " << j << std::endl;
        float3x3 value = element.second;
        result[i] += value * vector[j];
      }
    }
  }
};

}  // namespace blender