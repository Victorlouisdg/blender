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
  Array<std::map<int, float3x3>> *rows_pointer;
  Array<std::map<int, float3x3>> rows;

  SparseMatrix() = default;

  SparseMatrix(int n_rows)
  {
    std::cout << "SparseMatrix " << std::endl;
    this->n_rows = n_rows;
    rows_pointer = new Array<std::map<int, float3x3>>(n_rows);
    rows = *rows_pointer;
  }

  ~SparseMatrix()
  {
    std::cout << "~SparseMatrix " << std::endl;
    delete rows_pointer;
  }

  void insert(int row, int col, float3x3 value)
  {
    rows[row][col] = value;
  }

  float3x3 get(int row, int col)
  {
    return rows[row][col];
  }

  void multiply(const Span<float3> &vector, const MutableSpan<float3> &result)
  {
    for (int i : IndexRange(n_rows)) {
      result[i] = float3(0.0f);
      for (std::pair<int, float3x3> element : rows[i]) {
        int j = element.first;
        float3x3 value = element.second;
        result[i] += value * vector[j];
      }
    }
  }

  void multiply_block_diagonal_matrix(const Span<float3x3> &diagonal, const float f)
  {
    /* For each row of A: multiply each element in the row with the diagonal matrix' diagonal block on that row. */
    for (int i : IndexRange(n_rows)) {
      for (std::pair<int, float3x3> element : rows[i]) {
        int j = element.first;
        float3x3 value = element.second;
        rows[i][j] = f * diagonal[i] * value;
      }
    }
  }

  bool is_symmetric()
  {
    for (int i : IndexRange(n_rows)) {
      for (std::pair<int, float3x3> element : rows[i]) {
        int j = element.first;
        float3x3 value_ij = element.second;
        float3x3 value_ji = get(j, i);

        // std::cout << std::endl << value_ij << std::endl << value_ji << std::endl;

        for (int s : IndexRange(3)) {
          for (int t : IndexRange(3)) {
            if (value_ij.values[s][t] != value_ji.values[t][s]) {
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  void clear()
  {
    for (int i : IndexRange(n_rows)) {
      rows[i].clear();
    }
  }
};

/* These operations could probably be done through stl parallel algorithms. */
static float dot(const Span<float3> a, const Span<float3> b)
{
  // TODO: Assert a.length == b.length
  float result = 0.0f;
  for (int i : a.index_range()) {
    result += float3::dot(a[i], b[i]);
  }
  return result;
};

static void add(const MutableSpan<float3> &a, const Span<float3> &b)
{
  for (int i : a.index_range()) {
    a[i] += b[i];
  }
}

static void subtract(const MutableSpan<float3> &a, const Span<float3> &b)
{
  for (int i : a.index_range()) {
    a[i] -= b[i];
  }
}

/* TODO think of a better name. */
static void subtract_from(const MutableSpan<float3> &a, const Span<float3> &b)
{
  for (int i : a.index_range()) {
    a[i] = b[i] - a[i];
  }
}

static MutableSpan<float3> multiply_float_inplace(float f, const MutableSpan<float3> &a)
{
  for (int i : a.index_range()) {
    a[i] = f * a[i];
  }
  return a;
}

static MutableSpan<float3> multiply_float(float f,
                                          const Span<float3> &a,
                                          const MutableSpan<float3> &result)
{
  for (int i : a.index_range()) {
    result[i] = f * a[i];
  }
  return result;
}

static void apply_jacobi_preconditioning(SparseMatrix &A,
                                         const Span<float3> &r,
                                         const MutableSpan<float3> &d)
{
  for (int i : r.index_range()) {
    d[i] = float3::safe_divide(r[i], A.get(i, i).diagonal());
  }
}

/* Solvers */

/* TODO maybe make this into a PCG class that has member variables for the intermediate
 * calculations e.g. for d, Ad, r. */
static void solve_filtered_pcg(SparseMatrix &A,
                               const Span<float3> &b,
                               const MutableSpan<float3> &x)
{
  /* Filtered Preconditioned Conjugate Gradient solver for a system of linear equations.
   * This implementation is partially based off of "Large Steps in Cloth Simultion" by Baraff &
   * Witkin (BW98). BW98 added the filtering procedure, which is used to enforce cloth-object
   * collision constraint exactly. For the rest of the PCG solver they used the algorithm described
   * in "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain" by Shewchuk.
   * A large part of the paper builds up intuition, the PCG algorithm is presented in its entirety
   * on page 51.
   *
   * The parameters of this function represent the linear system:
   * A @ x = b
   * Note that x is allowed to be set to an initial guess for the solution of the system.
   */

  // Array<float3> d_temp = Array<float3>(A.n_rows);
  // Array<float3> Ad_temp = Array<float3>(A.n_rows);

  Array<float3> d = Array<float3>(A.n_rows);
  Array<float3> Ad = Array<float3>(A.n_rows);
  Array<float3> alpha_d = Array<float3>(A.n_rows);
  Array<float3> beta_d = Array<float3>(A.n_rows);
  Array<float3> r = Array<float3>(A.n_rows);
  Array<float3> s = Array<float3>(A.n_rows);

  /* r = b - Ax */
  A.multiply(x, r);  // r = Ax

  // for (int i : IndexRange(A.n_rows)) {
  //   std::cout << "r = " << r[i] << " x = " << x[i] << std::endl;
  // }

  subtract_from(r, b);  // r = b - Ax

  // c = precondition (r) -> c == d
  // apply_jacobi_preconditioning(A, r, d);

  /* d = r */
  for (int i : IndexRange(A.n_rows)) {
    d[i] = r[i];
  }

  float delta_old;
  float delta_new = dot(r, d);  // delta_new = rTr
  float delta_0 = delta_new;

  /* Arbitrarily chosen tolerance. */
  float tolerance = 0.00001f;

  for (int i : IndexRange(100)) {
    UNUSED_VARS(i);

    /* Early stopping when error has decreased enough. */
    if (delta_new < tolerance) {
      break;
    }

    A.multiply(d, Ad);
    float alpha = delta_new / dot(d, Ad);

    multiply_float(alpha, d, alpha_d);
    add(x, alpha_d);

    subtract(r, multiply_float_inplace(alpha, Ad));
    // apply_jacobi_preconditioning(A, r, s);

    delta_old = delta_new;
    // delta_new = dot(r, s);
    delta_new = dot(r, r);

    float beta = delta_new / delta_old;

    multiply_float_inplace(beta, d);
    // add(d, s);
    add(d, r);

    // std::cout << i << " " << delta_new << std::endl;

    // std::cout << "d_temp       d" << std::endl;
    // for (int i : IndexRange(A.n_rows)) {
    //   std::cout << d_temp[i] << "         " << d[i] << std::endl;
    // }

    // A.multiply(d_temp, Ad_temp);
    // std::cout << "conjugacy: " << dot(d, Ad_temp) << std::endl;
  }

  float3 sum = float3(0.0f);
  A.multiply(x, r);
  subtract_from(r, b);
  for (int i : IndexRange(A.n_rows)) {
    sum += r[i];
  }
  std::cout << "(PCG) Sum of r = " << sum << std::endl;
};

static void solve_gauss_seidel(SparseMatrix &A,
                               const Span<float3> &b,
                               const MutableSpan<float3> &x)
{

  Array<float3> r = Array<float3>(A.n_rows);
  for (int n : IndexRange(50)) {

    for (int i : IndexRange(A.n_rows)) {
      x[i] = b[i];

      float3x3 a_ii;
      for (std::pair<int, float3x3> element : A.rows[i]) {
        int j = element.first;

        float3x3 a_ij = element.second;
        if (i == j) {
          a_ii = a_ij;
        }
        else {
          x[i] -= a_ij * x[j];
        }
      }
      x[i] = a_ii.inverted() * x[i];
    }
  }

  float3 sum = float3(0.0f);
  A.multiply(x, r);
  subtract_from(r, b);
  for (int i : IndexRange(A.n_rows)) {
    sum += r[i];
  }
  std::cout << "(GS) Sum of r = " << sum << std::endl;
};

}  // namespace blender