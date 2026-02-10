#include <math.h>
#include <assert.h>
#include "compute_weights.h"
#include "yac_lapack_interface.h"
#include "geometry.h"

static void
inverse(double *A, size_t n)
{
// LAPACKE_dsytrf_work and LAPACKE_dsytri might not be available (see yac_lapack_interface.h).
#ifdef YAC_LAPACK_NO_DSYTR
  lapack_int ipiv[n + 1];
  double work[n * n];

  for (size_t i = 0; i < n + 1; ++i) ipiv[i] = 0;
  for (size_t i = 0; i < n * n; ++i) work[i] = 0;

  assert(!LAPACKE_dgetrf(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) n, A, (lapack_int) n, ipiv) || "internal ERROR: dgetrf");

  assert(!LAPACKE_dgetri_work(LAPACK_COL_MAJOR, (lapack_int) n, A, (lapack_int) n, ipiv, work, (lapack_int) (n * n))
         || "internal ERROR: dgetri");
#else
  lapack_int ipiv[n];
  double work[n];

  assert(!LAPACKE_dsytrf_work(LAPACK_COL_MAJOR, 'L', (lapack_int) n, A, (lapack_int) n, ipiv, work, (lapack_int) n)
         || "internal ERROR: dsytrf");

  assert(!LAPACKE_dsytri_work(LAPACK_COL_MAJOR, 'L', (lapack_int) n, A, (lapack_int) n, ipiv, work) || "internal ERROR: dsytri");

  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) A[j * n + i] = A[i * n + j];
#endif
}

void
yac_compute_weights_rbf(double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n, double *weights,
                        double const rbf_scale)
{
  double A[n][n];
  double a[n];

  double sum_d = 0.0, scale_d = 1.0;

  // compute distance matrix for all found source points
  for (size_t i = 0; i < n; ++i) A[i][i] = 0.0;
  for (size_t i = 0; i < n - 1; ++i)
    {
      for (size_t j = i + 1; j < n; ++j)
        {
          double d = get_vector_angle(src_coords[i], src_coords[j]);
          A[i][j] = d;
          A[j][i] = d;
          sum_d += d;
        }
    }

  // compute and apply scale factor for distance matrix
  if (sum_d > 0.0) scale_d = ((double) ((n - 1) * n)) / (2.0 * sum_d);
  scale_d /= rbf_scale;

  double sq_scale_d = scale_d * scale_d;

  // compute matrix A[n][n]
  // with A = rbf(A)
  //      rbf(a) => Radial basis function
  for (size_t i = 0; i < n; ++i)
    {
      for (size_t j = 0; j < n; ++j)
        {
          double d = A[i][j];
          // gauß rbf kernel
          A[i][j] = exp(-1.0 * d * d * sq_scale_d);
        }
    }

  // compute inverse of A
  inverse(&A[0][0], n);

  // compute a[NUM_NEIGH]
  // with a_i = rbf(d(x_i, y))
  //      x => vector containing the coordinates of the
  //           n nearest neighbours of the current
  //           target point
  //      y => coordinates of target point
  //      d(a, b) => great circle distance between point a and b
  //      rbf(a) => Radial basis function
  for (size_t i = 0; i < n; ++i)
    {
      double d = get_vector_angle(tgt_coord, src_coords[i]);
      // gauß rbf kernel
      a[i] = exp(-1.0 * d * d * sq_scale_d);
    }

  // compute weights
  // with w_i = SUM(A_inv[i][j]*a[j]) // j = 0..n-1
  for (size_t i = 0; i < n; ++i)
    {
      weights[i] = 0.0;
      for (size_t j = 0; j < n; ++j) weights[i] += A[i][j] * a[j];
    }
}
