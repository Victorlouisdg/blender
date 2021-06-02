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
  Array<std::map<int, float3x3>> rows;

  SparseMatrix() = default;

  SparseMatrix(int n_rows)
  {
    this->n_rows = n_rows;
    rows = Array<std::map<int, float3x3>>(n_rows);
  }

  void insert(int row, int col, float3x3 value)
  {
    rows[row][col] = value;
  }

  void add(int row, int col, float3x3 value)
  {
    rows[row][col] = value + rows[row][col];
  }

  float3x3 get(int row, int col)
  {
    return rows[row][col];
  }

  void multiply(const Span<float3> &vector, const MutableSpan<float3> &result)
  {
    for (int i : IndexRange(n_rows)) {
      result[i] = float3(0.0f);
      for (std::pair<int, float3x3> p : rows[i]) {
        int j = p.first;
        float3x3 value = p.second;
        result[i] += value * vector[j];
      }
    }
  }

  void multiply_float(const float f)
  {
    for (int i : IndexRange(n_rows)) {
      for (std::pair<int, float3x3> p : rows[i]) {
        int j = p.first;
        float3x3 value = p.second;
        rows[i][j] = f * value;
      }
    }
  }

  void add_matrix(SparseMatrix &B)
  {
    /* This function currently only works for matrices with the same sparsity pattern. */
    for (int i : IndexRange(n_rows)) {
      for (std::pair<int, float3x3> p : rows[i]) {
        int j = p.first;
        float3x3 value = p.second;
        rows[i][j] = value + B.get(i, j);
      }
    }
  }

  void clear()
  {
    for (int i : IndexRange(n_rows)) {
      rows[i].clear();
    }
  }
};

/* Span-based operations that should go elsewhere. */

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

/* Simple Jacobi Preconditioner. */
static void precondition(SparseMatrix &A, const Span<float3> &a, const MutableSpan<float3> &result)
{
  for (int i : a.index_range()) {
    result[i] = float3::safe_divide(a[i], A.get(i, i).diagonal());
  }
}

class ConjugateGradientSolver {
  /* Filtered Preconditioned Conjugate Gradient solver for a system of linear equations.
   *
   * This class is an almost direct implementation of the modified-pcg algorithm in
   * "Large Steps in Cloth Simultion" by Baraff & Witkin [BW98]. BW98 added the filtering
   * procedure, which is used to enforce cloth-object collision constraints exactly. For the rest
   * of the PCG solver they used the algorithm described in "An Introduction to the Conjugate
   * Gradient Method Without the Agonizing Pain" by Shewchuk.
   */

 public:
  int n;
  int max_iterations = 100;   /* If this is reached, something is probably wrong. */
  float tolerance = 0.000001f; /* Not sure what tolerance is good enough. */

  Array<float3> r;
  Array<float3> c;
  Array<float3> q;
  Array<float3> s;
  Array<float3> alpha_c;
  Array<float3> beta_c;

  std::map<int, float3> constraint_values;    /* Holds the zi's from the paper. */
  std::map<int, float3x3> constraint_filters; /* Holds the Si's from the paper.*/

  ConjugateGradientSolver() = default;

  ConjugateGradientSolver(int linear_system_size)
  {
    n = linear_system_size;

    r = Array<float3>(n);
    c = Array<float3>(n);
    q = Array<float3>(n);
    s = Array<float3>(n);
    alpha_c = Array<float3>(n);
    beta_c = Array<float3>(n);

    constraint_values = std::map<int, float3>();
    constraint_filters = std::map<int, float3x3>();
  }

  void setConstraint(int i, float3 value, float3x3 filter)
  {
    constraint_values[i] = value;
    constraint_filters[i] = filter;
  }

  void clearConstraints()
  {
    constraint_values.clear();
    constraint_filters.clear();
  }

  void filter(const MutableSpan<float3> &a)
  {
    for (std::pair<int, float3x3> p : constraint_filters) {
      int i = p.first;
      float3x3 Si = p.second;
      a[i] = Si * a[i];
    }
  }

  void initialize_solution(const MutableSpan<float3> &x)
  {
    // Not sure why "x.fill(float3(0.0f));" doesn't work.
    for (int i : x.index_range()) {
      x[i] = float3(0.0f);
    }

    for (std::pair<int, float3> p : constraint_values) {
      int i = p.first;
      float3 zi = p.second;
      x[i] = zi;
    }
  }

  /* The parameters of this function represent the linear system:
   * A @ x = b
   *
   * The memory for the solution x should be allocated by the caller, but this function will
   * initialize its values.
   */
  void solve(SparseMatrix &A, const Span<float3> &b, const MutableSpan<float3> &x)
  {
    initialize_solution(x); /* delta_v = z */

    A.multiply(x, r);    /* r = Ax */
    subtract_from(r, b); /* r = b - Ax */
    filter(r);

    precondition(A, r, c); /* c = P-1 r */
    filter(c);

    float delta_new = dot(r, c);

    for (int iteration : IndexRange(max_iterations)) {
      UNUSED_VARS(iteration);
      if (delta_new < tolerance) {
        // std::cout << "Early stopping at " << iteration << std::endl;
        break;
      }

      A.multiply(c, q); /* q = Ac*/
      filter(q);

      float alpha = delta_new / dot(c, q);

      multiply_float(alpha, c, alpha_c);
      add(x, alpha_c);

      multiply_float_inplace(alpha, q); /* q = alpha * q */
      subtract(r, q);

      /* TODO reuse q for s to save memory */
      precondition(A, r, s);
      float delta_old = delta_new;
      delta_new = dot(r, s);
      float beta = delta_new / delta_old;
      multiply_float_inplace(beta, c);  // c = beta * c
      add(c, s);                        // c = s + beta * c
      filter(c);
    }
  }
};

}  // namespace blender