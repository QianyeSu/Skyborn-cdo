// This source code is copied from PINGO version 1.5

/* ********************************** */
/* HEADER FOR PARALLEL EIGEN SOLUTION */
/*  -->SEE END OF ROUTINE             */
/* ********************************** */

#include <atomic>
#include <cstdlib>

#include "cdo_options.h"
#include "cdo_output.h"
#include "eigen_solution.h"
#include "cdo_omp.h"

constexpr double FNORM_PRECISION = 1.e-12;
constexpr int MAX_JACOBI_ITER = 12;

// global variables to handle environment settings
static double fnorm_precision;
static int max_jacobi_iter;

static std::atomic<size_t> n_finished;

static void
heap_sort(Varray<double> &eigenValues, Varray2D<double> &a, long n)
{
  long j_next;

  // First part of heap sort:
  for (long i = n / 2 - 1; i >= 0; i--)
  {
    for (long j = i; 2 * j + 1 < n; j = j_next)
    {
      auto k1 = 2 * j + 1;
      auto k2 = 2 * j + 2;
      j_next = j;
      if (eigenValues[k1] < eigenValues[j_next]) j_next = k1;
      if (k2 < n && eigenValues[k2] < eigenValues[j_next]) j_next = k2;
      if (j == j_next) break;
      std::swap(eigenValues[j], eigenValues[j_next]);
      std::swap(a[j], a[j_next]);
    }
  }
  // Second part of head sort:
  for (long i = n - 1; i > 0; i--)
  {
    std::swap(eigenValues[0], eigenValues[i]);
    std::swap(a[0], a[i]);
    for (long j = 0; 2 * j + 1 < i; j = j_next)
    {
      auto k1 = 2 * j + 1;
      auto k2 = 2 * j + 2;
      j_next = j;
      if (eigenValues[k1] < eigenValues[j_next]) j_next = k1;
      if (k2 < i && eigenValues[k2] < eigenValues[j_next]) j_next = k2;
      if (j == j_next) break;
      std::swap(eigenValues[j], eigenValues[j_next]);
      std::swap(a[j], a[j_next]);
    }
  }
}

static void
make_symmetric_matrix_triangular(Varray2D<double> &a, long n, Varray<double> &d, Varray<double> &e)
{
  double f, g, h, hh, scale;

  for (long i = n - 1; i >= 1; i--)
  {
    h = scale = 0;
    if (i > 1)
    {
      for (long k = 0; k < i; ++k) scale += std::fabs(a[i][k]);
      if (fp_is_equal(scale, 0.))
        e[i] = a[i][i - 1];
      else
      {
        for (long k = 0; k < i; ++k)
        {
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f = a[i][i - 1];
        g = (f >= 0) ? -sqrt(h) : std::sqrt(h);
        e[i] = scale * g;
        h -= f * g;
        a[i][i - 1] = f - g;
        f = 0;
        for (long j = 0; j < i; ++j)
        {
          a[j][i] = a[i][j] / h;
          g = 0;
          for (long k = 0; k <= j; ++k) g += a[j][k] * a[i][k];
          for (long k = j + 1; k < i; ++k) g += a[k][j] * a[i][k];
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        hh = f / (2 * h);
        for (long j = 0; j < i; ++j)
        {
          f = a[i][j];
          e[j] = g = e[j] - hh * f;
          for (long k = 0; k <= j; ++k) a[j][k] -= f * e[k] + g * a[i][k];
        }
      }
    }
    else
      e[i] = a[i][i - 1];

    d[i] = h;
  }

  d[0] = e[0] = 0;
  for (long i = 0; i < n; ++i)
  {
    if (std::fabs(d[i]) > 0)
    {
      for (long j = 0; j < i; ++j)
      {
        g = 0;
        for (long k = 0; k < i; ++k) g += a[i][k] * a[k][j];
        for (long k = 0; k < i; ++k) a[k][j] -= g * a[k][i];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1;
    for (long j = 0; j < i; ++j) a[j][i] = a[i][j] = 0;
  }
}

static double
pythagoras(double a, double b)
{
  auto abs_a = std::fabs(a);
  auto abs_b = std::fabs(b);
  if (abs_a > abs_b)
  {
    auto sqr = abs_b / abs_a;
    sqr *= sqr;
    return abs_a * std::sqrt(1.0 + sqr);
  }
  else if (abs_b > abs_a)
  {
    auto sqr = abs_a / abs_b;
    sqr *= sqr;
    return abs_b * std::sqrt(1.0 + sqr);
  }
  else
    return M_SQRT2 * abs_a;
}

static void
eigen_solution_of_triangular_matrix(Varray<double> &d, Varray<double> &e, long n, Varray2D<double> &a, const char *prompt)
{
  constexpr double eps = 1.e-6;
  constexpr long MAX_ITER = 1000;

  for (long i = 1; i < n; ++i) e[i - 1] = e[i];

  e[n - 1] = 0.0;
  for (long l = 0; l < n; ++l)
  {
    long iter = 0;
    while (1)
    {
      long m;
      for (m = l; m < n - 1; ++m)
        if (std::fabs(e[m]) <= eps * (std::fabs(d[m]) + std::fabs(d[m + 1]))) break;
      if (m == l)
      {
        // printf("found solution after %i Iteration\n", iter++);
        break;
      }
      iter++;
      if (iter == MAX_ITER)
      {
        std::fprintf(stderr, "%s: ERROR! Too many iterations while determining the eigensolution!\n", prompt);
        std::exit(1);
      }
      double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
      double r = pythagoras(g, 1.0);
      g = d[m] - d[l] + e[l] / (g + ((std::fabs(g) > 0) ? ((g >= 0) ? std::fabs(r) : -std::fabs(r)) : r));
      double s = 1.0;
      double c = 1.0;
      double p = 0.0;
      long i;
      for (i = m - 1; i >= l; i--)
      {
        double f = s * e[i];
        double b = c * e[i];
        e[i + 1] = r = pythagoras(f, g);

        if (fp_is_equal(r, 0.0))
        {
          d[i + 1] -= p;
          e[m] = 0.0;
          break;
        }

        s = f / r;
        c = g / r;
        g = d[i + 1] - p;
        r = (d[i] - g) * s + 2.0 * c * b;
        p = s * r;
        d[i + 1] = g + p;
        g = c * r - b;
        for (long k = 0; k < n; ++k)
        {
          f = a[k][i + 1];
          a[k][i + 1] = s * a[k][i] + c * f;
          a[k][i] = c * a[k][i] - s * f;
        }
      }

      if (fp_is_equal(r, 0.0) && i >= l) continue;

      d[l] -= p;
      e[l] = g;
      e[m] = 0.0;
    }
  }
}

void
eigen_solution_of_symmetric_matrix(Varray2D<double> &a, Varray<double> &eigenValues, size_t n, const char *prompt)
// After return the rows (!!!) of a are the eigenvectors
{
  {
    Varray<double> e(n);
    make_symmetric_matrix_triangular(a, n, eigenValues, e);
    eigen_solution_of_triangular_matrix(eigenValues, e, n, a, prompt);
  }

  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) std::swap(a[i][j], a[j][i]);

  heap_sort(eigenValues, a, n);
}
/*
static int
lu_decomposition(double **a, long n, int *index, int *sign)
{
  int imax = 0;
  double big, sum, temp;

  std::vector<double> v(n);
  *sign = 1;
  for (long i = 0; i < n; ++i)
    {
      big = 0;
      for (long j = 0; j < n; ++j)
        if ((temp = std::fabs(a[i][j])) > big) big = temp;

      if (fp_is_equal(big, 0.)) return 0;

      v[i] = 1 / big;
    }
  for (long j = 0; j < n; ++j)
    {
      for (long i = 0; i < j; ++i)
        {
          sum = a[i][j];
          for (long k = 0; k < i; ++k) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0;
      for (long i = j; i < n; ++i)
        {
          sum = a[i][j];
          for (long k = 0; k < j; ++k) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((temp = v[i] * std::fabs(sum)) >= big)
            {
              big = temp;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (long k = 0; k < n; ++k)
            {
              temp = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = temp;
            }
          *sign = -*sign;
          v[imax] = v[j];
        }
      index[j] = imax;

      if (fp_is_equal(a[j][j], 0.)) return 0;

      if (j != n)
        {
          temp = 1 / a[j][j];
          for (long i = j + 1; i < n; ++i) a[i][j] *= temp;
        }
    }

  return 1;
}

static void
lu_backsubstitution(double **a, long n, int *index, double *b)
{
  long ii = 0;
  for (long i = 0; i < n; ++i)
    {
      int ip = index[i];
      double sum = b[ip];
      b[ip] = b[i];

      if (ii)
        for (long j = ii; j < i; ++j) sum -= a[i][j] * b[j];
      else if (std::fabs(sum) > 0)
        ii = i;

      b[i] = sum;
    }
  for (long i = n - 1; i >= 0; i--)
    {
      double sum = b[i];
      for (long j = i + 1; j < n; ++j) sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}

static int
solution_of_linear_equation(double **a, double *b, int n)
{
  std::vector<int> index(n);

  int sign;
  int not_singular = lu_decomposition(a, n, index.data(), &sign);

  if (not_singular) lu_backsubstitution(a, n, index.data(), b);

  return not_singular;
}

static int
inverse_of_matrix(double **a, double **b, long n)
{
  int sign;

  std::vector<int> index(n);
  std::vector<double> col(n);

  int not_singular = lu_decomposition(a, n, index.data(), &sign);

  if (not_singular)
    {
      for (long i = 0; i < n; ++i)
        {
          for (long j = 0; j < n; ++j) col[j] = 0;
          col[i] = 1;
          lu_backsubstitution(a, n, index.data(), col.data());
          for (long j = 0; j < n; ++j) b[j][i] = col[j];
        }
    }

  return not_singular;
}
*/

/* ******************************************************************************** */
/* This routine rotates columns/rows i and j of a symmetric Matrix M in a fashion,  */
/* thus that the dot product of columns i and j 0 afterwards                        */
/*                                                                                  */
/* As this is done by a right-multiplication with a rotation matrix, which only     */
/* changes columns i and j, this can be carried out for n/2 pairs of columns at     */
/* the same time.                                                                   */
/* ******************************************************************************** */
static void
annihilate_1side(Varray2D<double> &M, size_t i, size_t j, size_t n)
{
  i--;
  j--;

  if (j < i) std::swap(i, j);

  auto &Mi = M[i];
  auto &Mj = M[j];

  double alpha{ 0.0 }, beta{ 0.0 }, gamma{ 0.0 };
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : alpha) reduction(+ : beta) reduction(+ : gamma)
#endif
  for (size_t r = 0; r < n; ++r)
  {
    alpha += Mj[r] * Mj[r];
    beta += Mi[r] * Mi[r];
    gamma += Mi[r] * Mj[r];
  }

  // 2011-08-15 Cedrick Ansorge: bug fix
  // auto tmp = std::fabs(gamma / std::sqrt(alpha / beta));
  auto tmp = std::fabs(gamma / std::sqrt(alpha * beta));
  if (tmp < fnorm_precision)
  {
    n_finished++;
    return;
  }

  auto zeta = (beta - alpha) / (2.0 * gamma);  // tan(2*theta)
  auto tk = 1.0 / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
  tk = (zeta > 0) ? tk : -tk;               // = cot(2*theta)
  auto ck = 1.0 / std::sqrt(1. + tk * tk);  // = cos(theta)
  auto sk = ck * tk;                        // = sin(theta)

  // calculate a_i,j - tilde
  for (size_t r = 0; r < n; ++r)
  {
    auto mi = Mi[r];
    auto mj = Mj[r];
    Mi[r] = ck * mi + sk * mj;
    Mj[r] = -sk * mi + ck * mj;
  }
}

static int
jacobi_1side(Varray2D<double> &M, Varray<double> &A, size_t n)
{
  Varray<size_t> annihilations1(n * n), annihilations2(n * n);

  size_t count = 0;
  for (size_t k = 1; k <= n; ++k)
  {
    if (k < n)
    {
      {
        auto nmax = (size_t) std::ceil(0.5 * (n - k));
        for (size_t i = 1; i <= nmax; ++i)
        {
          auto j = n - k + 2 - i;
          annihilations1[count] = i;
          annihilations2[count] = j;
          count++;
        }
      }
      if (k > 2)
      {
        auto nmax = n - (size_t) std::floor(0.5 * k);
        for (size_t i = n - k + 2; i <= nmax; ++i)
        {
          auto j = 2 * n - k + 2 - i;
          annihilations1[count] = i;
          annihilations2[count] = j;
          count++;
        }
      }
    }
    else if (k == n)
    {
      auto nmax = (size_t) std::ceil(0.5 * n);
      for (size_t i = 2; i <= nmax; ++i)
      {
        auto j = n + 2 - i;
        annihilations1[count] = i;
        annihilations2[count] = j;
        count++;
      }
    }
  }

  // std::fprintf(stderr, "%d annihilations per sweep\n", count);

  n_finished = 0;

  // override global openmp settings works
  // omp_set_num_threads(2);

  int n_iter = 0;
  while (n_iter < max_jacobi_iter && n_finished < count)
  {
    n_finished = 0;
    if (n % 2 == 1)
    {
      for (size_t m = 0; m < n; ++m)
      {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for (size_t i = 0; i < (n / 2); ++i)
        {
          auto idx = m * (n / 2) + i;
          auto i_ann = annihilations1[idx];
          auto j_ann = annihilations2[idx];
          if (i_ann && j_ann && i_ann != j_ann) annihilate_1side(M, i_ann, j_ann, n);
        }
      }
    }
    else
    {  // n%2 == 0
      for (size_t m = 0; m < n; ++m)
      {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for (size_t i = 0; i < (n / 2 - (m % 2)); ++i)
        {
          auto idx = m / 2 * (n / 2 + n / 2 - 1);
          if (m % 2) idx += n / 2;
          auto i_ann = annihilations1[idx + i];
          auto j_ann = annihilations2[idx + i];
          if (i_ann && j_ann && i_ann != j_ann) annihilate_1side(M, i_ann, j_ann, n);
        }
      }
    }
    n_iter++;
  }

  if (Options::cdoVerbose) cdo_print("Finished one-sided jacobi scheme for eigenvalue computation after %i iterations", n_iter);

  // std::fprintf(stderr,"finished after %i sweeps (n_finished %i)\n",n_iter,n_finished);

  if (n_iter == max_jacobi_iter && n_finished < count)
  {
    std::fprintf(stderr,
            "jacobi_1side (Warning): Eigenvalue computation with one-sided jacobi scheme did not converge properly.\n"
            "                        %zu of %zu pairs of columns did not achieve requested orthogonality of %g\n",
            count - n_finished, count, fnorm_precision);

    if (n_finished == 0)
    {
      // Do not overwrite results in case of insufficient convergence
      cdo_warning("Setting Matrix and Eigenvalues to 0 before return");
      for (size_t i = 0; i < n; ++i) memset(M[i].data(), 0, n * sizeof(double));
      memset(A.data(), 0, n * sizeof(double));
      return -1;
    }
  }

  // calculate  eigen values as std::sqrt(||m_i||)
  for (size_t i = 0; i < n; ++i)
  {
    A[i] = 0.0;
    for (size_t r = 0; r < n; ++r) A[i] += M[i][r] * M[i][r];
    A[i] = std::sqrt(A[i]);
    for (size_t r = 0; r < n; ++r) M[i][r] /= A[i];
  }

  heap_sort(A, M, n);

  return n_iter;
}

/* ******************************************************************************** */
/*                                                                                  */
/*   P A R A L L E L   S O L U T I O N   O F   T H E   E I G E N   P R O B L E M    */
/*                     WITH ONE SIDED JACOBI ALGORITHM                              */
/*                                                                                  */
/* ******************************************************************************** */

void
parallel_eigen_solution_of_symmetric_matrix(Varray2D<double> &M, Varray<double> &A, size_t n, const char func[])
{
  (void) (func);  // unused

  // Get Environment variables if set
  static bool initEnv = true;
  if (initEnv)
  {
    initEnv = false;
    auto envstr = getenv("MAX_JACOBI_ITER");
    max_jacobi_iter = envstr ? atoi(envstr) : MAX_JACOBI_ITER;
    if (Options::cdoVerbose) cdo_print("Using MAX_JACOBI_ITER %i from %s", max_jacobi_iter, envstr ? "Environment" : "default");

    envstr = getenv("FNORM_PRECISION");
    fnorm_precision = envstr ? std::strtod(envstr, nullptr) : FNORM_PRECISION;
    if (Options::cdoVerbose) cdo_print("Using FNORM_PRECISION %g from %s", fnorm_precision, envstr ? "Environment" : "default");
  }

  // eigen_solution_of_symmetric_matrix(M, A, n, func);
  jacobi_1side(M, A, n);

  return;
}
