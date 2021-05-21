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
    for (int i : IndexRange(n_rows)) {
      result[i] = float3(0.0f);
      for (std::pair<int, float3x3> element : (*rows)[i]) {
        int j = element.first;
        float3x3 value = element.second;
        result[i] += value * vector[j];
      }
    }
  }
};

float dot(const Span<float3> a, const Span<float3> b)
{
  // TODO: Assert a.length == b.length
  float result = 0.0f;
  for (int i : a.index_range()) {
    result += float3::dot(a[i], b[i]);
  }
};

/* Solvers */
void solve_filtered_pcg(const SparseMatrix &A, const Span<float3> &b, const MutableSpan<float3> &x)
{

  Array<float3> d = Array<float3>(A.n_rows);
  Array<float3> Ad = Array<float3>(A.n_rows);
  Array<float3> r = Array<float3>(A.n_rows);

  float delta_new = dot(r, d);

  // r = A * x - b;
  // p = -r;

  // while (!converged) {
  //   float rTr = r.squared_length();

  //   A.multiply(p, Ap);

  //   float alpha = rTr / (p.dot(Ap));
  //   x += alpha * p;
  //   r += alpha * Ap;

  //   float rTr_new = r.squared_length();
  //   float beta = rTr_new / rTr;

  //   p = -r + beta
  // }
};

}  // namespace blender