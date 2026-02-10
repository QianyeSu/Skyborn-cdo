/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "varray.h"
#include "cdo_options.h"
#include "cdo_omp.h"

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifdef SX
constexpr long NFFT = 1024;
#else
constexpr long NFFT = 8;
#endif

constexpr double QUA = 0.25;
constexpr double QT5 = 0.559016994374947;

constexpr double S36 = 0.587785252292473;
constexpr double S60 = 0.866025403784437;
constexpr double S72 = 0.951056516295154;

constexpr double SQ2 = 0.707106781186547524401;

constexpr double D60 = S60 + S60;

int
fft_set(double *trigs, long *ifax, long n)
{
  if (n == 0 || trigs == nullptr || ifax == nullptr)
  {
    fprintf(stderr, "fft_set: parameter not initialized!\n");
    return -2;
  }

  const long len = n;
  const long nhl = n / 2;
  const double del = 4.0 * std::asin(1.0) / n;

  for (long k = 0; k < nhl; ++k)
  {
    const double angle = k * del;
    trigs[2 * k] = std::cos(angle);
    trigs[2 * k + 1] = std::sin(angle);
  }

  for (long k = 0; k < 9; ++k) ifax[k] = 0;

  ifax[9] = n;

  long nfax = 0;
  // clang-format off
  if    (n % 8 == 0) { ifax[++nfax] = 8; n /= 8; }
  while (n % 6 == 0) { ifax[++nfax] = 6; n /= 6; }
  while (n % 5 == 0) { ifax[++nfax] = 5; n /= 5; }
  while (n % 4 == 0) { ifax[++nfax] = 4; n /= 4; }
  while (n % 3 == 0) { ifax[++nfax] = 3; n /= 3; }
  if    (n % 2 == 0) { ifax[++nfax] = 2; n /= 2; }
  // clang-format on

  ifax[0] = nfax;

#ifdef CRAY
#pragma _CRI novector
#endif
#ifdef SX
#pragma vdir novector
#endif
#ifdef __uxp__
#pragma loop scalar
#endif
  for (long k = 0; k < nfax / 2; ++k)
  {
    const long j = ifax[k + 1];
    ifax[k + 1] = ifax[nfax - k];
    ifax[nfax - k] = j;
  }
  /*
  long prod = ifax[1];
  for (long k = 2; k <= nfax; ++k) prod *= ifax[k];
  printf("nfax = %ld  len %ld  prod %ld  n %ld\n", nfax, len, prod, n);
  for (long k = 1; k <= nfax; ++k) printf(" %ld", ifax[k]);
  printf("\n");
  */
  if (n > 8 || n == 7)
  {
    fprintf(stderr, "FFT does not work with len=%ld (n=%ld)!\n", len, n);
    return -1;
  }

  return 0;
}

static int
rpassc_2(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * jink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long j0 = 0;
  long i1 = i0 + inc1 * (m + m - la);
  long j1 = j0 + jink;

  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i] + a[i1 + i];
        c[j1 + j] = a[i0 + i] - a[i1 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    i0 += iink;
    iink += iink;
    i1 -= iink;
    ibase = 0;
    jbase += jump;
    jump += jump + jink;

    if (i0 != i1)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        ibase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a0m1 = a[i0 + i] - a[i1 + i];
            const double b0p1 = b[i0 + i] + b[i1 + i];

            c[j0 + j] = a[i0 + i] + a[i1 + i];
            d[j0 + j] = b[i0 + i] - b[i1 + i];
            c[j1 + j] = c1 * a0m1 - s1 * b0p1;
            d[j1 + j] = s1 * a0m1 + c1 * b0p1;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        i0 += iink;
        i1 -= iink;
        jbase += jump;
      } /* End FORK */
      if (i0 > i1) return 0;
    } /* End (i0 != i1) */
    ibase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i];
        c[j1 + j] = -b[i0 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i]);
        c[j1 + j] = 2.0 * (a[i0 + i] - a[i1 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
rpassc_3(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * jink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long j0 = 0;
  long i1 = i0 + inc1 * (m + m - la);
  long i2 = i1;
  long j1 = j0 + jink;
  long j2 = j1 + jink;

  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double afa1 = a[i0 + i] - 0.5 * a[i1 + i];
        const double bfa1 = S60 * b[i1 + i];

        c[j0 + j] = a[i0 + i] + a[i1 + i];
        c[j1 + j] = afa1 - bfa1;
        c[j2 + j] = afa1 + bfa1;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    i0 += iink;
    iink += iink;
    i1 += iink;
    i2 -= iink;
    jbase += jump;
    jump += jump + jink;

    if (i0 != i2)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        ibase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a1p2 = a[i0 + i] - 0.5 * (a[i1 + i] + a[i2 + i]);
            const double b1m2 = b[i0 + i] - 0.5 * (b[i1 + i] - b[i2 + i]);
            const double a1m2 = S60 * (a[i1 + i] - a[i2 + i]);
            const double b1p2 = S60 * (b[i1 + i] + b[i2 + i]);

            const double a0mm = a1p2 - b1p2;
            const double a0mp = a1p2 + b1p2;
            const double b0mm = b1m2 - a1m2;
            const double b0mp = b1m2 + a1m2;

            c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
            d[j0 + j] = b[i0 + i] + b[i1 + i] - b[i2 + i];
            c[j1 + j] = c1 * a0mm - s1 * b0mp;
            d[j1 + j] = s1 * a0mm + c1 * b0mp;
            c[j2 + j] = c2 * a0mp - s2 * b0mm;
            d[j2 + j] = s2 * a0mp + c2 * b0mm;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        i0 += iink;
        i1 += iink;
        i2 -= iink;
        jbase += jump;
      } /* End FORK */
      if (i0 > i2) return 0;
    } /* End (i0 != i2) */
    ibase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0mp = 0.5 * a[i0 + i];
        const double b0mp = S60 * b[i0 + i];

        c[j0 + j] = a[i0 + i] + a[i1 + i];
        c[j1 + j] = a0mp - a[i1 + i] - b0mp;
        c[j2 + j] = a[i1 + i] - a0mp - b0mp;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0mp = 2.0 * a[i0 + i] - a[i1 + i];
        const double b0mp = D60 * b[i1 + i];

        c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i]);
        c[j1 + j] = a0mp - b0mp;
        c[j2 + j] = a0mp + b0mp;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
rpassc_4(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * jink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long j0 = 0;
  long i1 = i0 + inc1 * (m + m - la);
  long i3 = i1;
  long i2 = i1 + inc1 * (m + m);
  long j1 = j0 + jink;
  long j2 = j1 + jink;
  long j3 = j2 + jink;

  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p2 = a[i0 + i] + a[i2 + i];
        const double a0m2 = a[i0 + i] - a[i2 + i];

        c[j0 + j] = a0p2 + a[i1 + i];
        c[j1 + j] = a0m2 - b[i1 + i];
        c[j2 + j] = a0p2 - a[i1 + i];
        c[j3 + j] = a0m2 + b[i1 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    i0 += iink;
    iink += iink;
    i1 += iink;
    i2 -= iink;
    i3 -= iink;
    jbase += jump;
    jump += jump + jink;

    if (i1 != i2)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        ibase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a0p2 = a[i0 + i] + a[i2 + i];
            const double a0m2 = a[i0 + i] - a[i2 + i];
            const double a1p3 = a[i1 + i] + a[i3 + i];
            const double a1m3 = a[i1 + i] - a[i3 + i];
            const double b0p2 = b[i0 + i] + b[i2 + i];
            const double b0m2 = b[i0 + i] - b[i2 + i];
            const double b1p3 = b[i1 + i] + b[i3 + i];
            const double b1m3 = b[i1 + i] - b[i3 + i];

            const double a0p2ma1p3 = a0p2 - a1p3;
            const double a0m2pb1p3 = a0m2 + b1p3;
            const double a0m2mb1p3 = a0m2 - b1p3;
            const double b0p2pa1m3 = b0p2 + a1m3;
            const double b0p2ma1m3 = b0p2 - a1m3;
            const double b0m2mb1m3 = b0m2 - b1m3;

            c[j0 + j] = a0p2 + a1p3;
            d[j0 + j] = b0m2 + b1m3;
            c[j2 + j] = c2 * a0p2ma1p3 - s2 * b0m2mb1m3;
            d[j2 + j] = s2 * a0p2ma1p3 + c2 * b0m2mb1m3;
            c[j1 + j] = c1 * a0m2mb1p3 - s1 * b0p2pa1m3;
            d[j1 + j] = s1 * a0m2mb1p3 + c1 * b0p2pa1m3;
            c[j3 + j] = c3 * a0m2pb1p3 - s3 * b0p2ma1m3;
            d[j3 + j] = s3 * a0m2pb1p3 + c3 * b0p2ma1m3;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        i0 += iink;
        i1 += iink;
        i2 -= iink;
        i3 -= iink;
        jbase += jump;
      } /* End FORK */
      if (i1 > i2) return 0;
    } /* End (i1 != i2) */
    ibase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0m1 = a[i0 + i] - a[i1 + i];
        const double b0p1 = b[i0 + i] + b[i1 + i];

        c[j0 + j] = a[i0 + i] + a[i1 + i];
        c[j2 + j] = b[i1 + i] - b[i0 + i];

        c[j1 + j] = SQ2 * (a0m1 - b0p1);
        c[j3 + j] = -SQ2 * (a0m1 + b0p1);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p2 = a[i0 + i] + a[i2 + i];
        const double a0m2 = a[i0 + i] - a[i2 + i];

        c[j0 + j] = 2.0 * (a0p2 + a[i1 + i]);
        c[j1 + j] = 2.0 * (a0m2 - b[i1 + i]);
        c[j2 + j] = 2.0 * (a0p2 - a[i1 + i]);
        c[j3 + j] = 2.0 * (a0m2 + b[i1 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
rpassc_5(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * jink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long j0 = 0;
  long i1 = i0 + inc1 * (m + m - la);
  long i4 = i1;
  long i2 = i1 + inc1 * (m + m);
  long i3 = i2;
  long j1 = j0 + jink;
  long j2 = j1 + jink;
  long j3 = j2 + jink;
  long j4 = j3 + jink;

  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a1p2 = QUA * (a[i1 + i] + a[i2 + i]);
        const double a1m2 = QT5 * (a[i1 + i] - a[i2 + i]);

        const double a0mp = a[i0 + i] - a1p2 + a1m2;
        const double a0mm = a[i0 + i] - a1p2 - a1m2;

        const double b136 = b[i1 + i] * S36;
        const double b172 = b[i1 + i] * S72;
        const double b236 = b[i2 + i] * S36;
        const double b272 = b[i2 + i] * S72;

        c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
        c[j1 + j] = a0mp - b172 - b236;
        c[j2 + j] = a0mm - b136 + b272;
        c[j3 + j] = a0mm + b136 - b272;
        c[j4 + j] = a0mp + b172 + b236;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    i0 += iink;
    iink += iink;
    i1 += iink;
    i2 += iink;
    i3 -= iink;
    i4 -= iink;
    jbase += jump;
    jump += jump + jink;

    if (i1 != i3)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const long ke = kd + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        const double c4 = trigs[ke];
        const double s4 = trigs[ke + 1];
        ibase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a10 = (a[i0 + i] - 0.25 * ((a[i1 + i] + a[i4 + i]) + (a[i2 + i] + a[i3 + i])))
                               + QT5 * ((a[i1 + i] + a[i4 + i]) - (a[i2 + i] + a[i3 + i]));
            const double a20 = (a[i0 + i] - 0.25 * ((a[i1 + i] + a[i4 + i]) + (a[i2 + i] + a[i3 + i])))
                               - QT5 * ((a[i1 + i] + a[i4 + i]) - (a[i2 + i] + a[i3 + i]));
            const double b10 = (b[i0 + i] - 0.25 * ((b[i1 + i] - b[i4 + i]) + (b[i2 + i] - b[i3 + i])))
                               + QT5 * ((b[i1 + i] - b[i4 + i]) - (b[i2 + i] - b[i3 + i]));
            const double b20 = (b[i0 + i] - 0.25 * ((b[i1 + i] - b[i4 + i]) + (b[i2 + i] - b[i3 + i])))
                               - QT5 * ((b[i1 + i] - b[i4 + i]) - (b[i2 + i] - b[i3 + i]));

            const double a11 = S72 * (b[i1 + i] + b[i4 + i]) + S36 * (b[i2 + i] + b[i3 + i]);
            const double a21 = S36 * (b[i1 + i] + b[i4 + i]) - S72 * (b[i2 + i] + b[i3 + i]);
            const double b11 = S72 * (a[i1 + i] - a[i4 + i]) + S36 * (a[i2 + i] - a[i3 + i]);
            const double b21 = S36 * (a[i1 + i] - a[i4 + i]) - S72 * (a[i2 + i] - a[i3 + i]);

            c[j0 + j] = a[i0 + i] + ((a[i1 + i] + a[i4 + i]) + (a[i2 + i] + a[i3 + i]));
            d[j0 + j] = b[i0 + i] + ((b[i1 + i] - b[i4 + i]) + (b[i2 + i] - b[i3 + i]));
            c[j1 + j] = c1 * (a10 - a11) - s1 * (b10 + b11);
            d[j1 + j] = s1 * (a10 - a11) + c1 * (b10 + b11);
            c[j4 + j] = c4 * (a10 + a11) - s4 * (b10 - b11);
            d[j4 + j] = s4 * (a10 + a11) + c4 * (b10 - b11);
            c[j2 + j] = c2 * (a20 - a21) - s2 * (b20 + b21);
            d[j2 + j] = s2 * (a20 - a21) + c2 * (b20 + b21);
            c[j3 + j] = c3 * (a20 + a21) - s3 * (b20 - b21);
            d[j3 + j] = s3 * (a20 + a21) + c3 * (b20 - b21);
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        i0 += iink;
        i1 += iink;
        i2 += iink;
        i3 -= iink;
        i4 -= iink;
        jbase += jump;
      } /* End FORK */
      if (i1 > i3) return 0;
    } /* End (i1 != i3) */
    ibase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
        c[j1 + j]
            = (QT5 * (a[i0 + i] - a[i1 + i]) + (0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) - (S36 * b[i0 + i] + S72 * b[i1 + i]);
        c[j4 + j]
            = -(QT5 * (a[i0 + i] - a[i1 + i]) + (0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) - (S36 * b[i0 + i] + S72 * b[i1 + i]);
        c[j2 + j]
            = (QT5 * (a[i0 + i] - a[i1 + i]) - (0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) - (S72 * b[i0 + i] - S36 * b[i1 + i]);
        c[j3 + j]
            = -(QT5 * (a[i0 + i] - a[i1 + i]) - (0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) - (S72 * b[i0 + i] - S36 * b[i1 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else
  {
    const double qqrt5 = 2.0 * QT5;
    const double ssin36 = 2.0 * S36;
    const double ssin72 = 2.0 * S72;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i] + a[i2 + i]);
        c[j1 + j] = (2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) + qqrt5 * (a[i1 + i] - a[i2 + i]))
                    - (ssin72 * b[i1 + i] + ssin36 * b[i2 + i]);
        c[j2 + j] = (2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) - qqrt5 * (a[i1 + i] - a[i2 + i]))
                    - (ssin36 * b[i1 + i] - ssin72 * b[i2 + i]);
        c[j3 + j] = (2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) - qqrt5 * (a[i1 + i] - a[i2 + i]))
                    + (ssin36 * b[i1 + i] - ssin72 * b[i2 + i]);
        c[j4 + j] = (2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) + qqrt5 * (a[i1 + i] - a[i2 + i]))
                    + (ssin72 * b[i1 + i] + ssin36 * b[i2 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
rpassc_6(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * jink;

  long ibase = 0;
  long jbase = 0;

  long ia = 0;
  long ib = ia + (2 * m - la) * inc1;
  long ic = ib + 2 * m * inc1;
  long id = ic + 2 * m * inc1;
  long ie = ic;
  long iF = ib;
  long ja = 0;
  long jb = ja + jink;
  long jc = jb + jink;
  long jd = jc + jink;
  long je = jd + jink;
  long jf = je + jink;

  if (la != m) /* go to 690 */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = (a[ia + i] + a[id + i]) + (a[ib + i] + a[ic + i]);
        c[jd + j] = (a[ia + i] - a[id + i]) - (a[ib + i] - a[ic + i]);
        c[jb + j] = ((a[ia + i] - a[id + i]) + 0.5 * (a[ib + i] - a[ic + i])) - S60 * (b[ib + i] + b[ic + i]);
        c[jf + j] = ((a[ia + i] - a[id + i]) + 0.5 * (a[ib + i] - a[ic + i])) + S60 * (b[ib + i] + b[ic + i]);
        c[jc + j] = ((a[ia + i] + a[id + i]) - 0.5 * (a[ib + i] + a[ic + i])) - S60 * (b[ib + i] - b[ic + i]);
        c[je + j] = ((a[ia + i] + a[id + i]) - 0.5 * (a[ib + i] + a[ic + i])) + S60 * (b[ib + i] - b[ic + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    ia += iink;
    iink += iink;
    ib += iink;
    ic += iink;
    id -= iink;
    ie -= iink;
    iF -= iink;
    jbase += jump;
    jump += jump + jink;

    if (ic != id) /* go to 660 */
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const long ke = kd + kb;
        const long kf = ke + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        const double c4 = trigs[ke];
        const double s4 = trigs[ke + 1];
        const double c5 = trigs[kf];
        const double s5 = trigs[kf + 1];
        ibase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            double a11 = a[ie + i] + a[ib + i] + a[ic + i] + a[iF + i];
            double a20 = a[ia + i] + a[id + i] - 0.5 * a11;
            double a21 = S60 * ((a[ie + i] + a[ib + i]) - (a[ic + i] + a[iF + i]));
            double b11 = b[ib + i] - b[ie + i] + b[ic + i] - b[iF + i];
            double b20 = b[ia + i] - b[id + i] - 0.5 * b11;
            double b21 = S60 * ((b[ib + i] - b[ie + i]) - (b[ic + i] - b[iF + i]));

            c[ja + j] = a[ia + i] + a[id + i] + a11;
            d[ja + j] = b[ia + i] - b[id + i] + b11;
            c[jc + j] = c2 * (a20 - b21) - s2 * (b20 + a21);
            d[jc + j] = s2 * (a20 - b21) + c2 * (b20 + a21);
            c[je + j] = c4 * (a20 + b21) - s4 * (b20 - a21);
            d[je + j] = s4 * (a20 + b21) + c4 * (b20 - a21);

            a11 = (a[ie + i] - a[ib + i]) + (a[ic + i] - a[iF + i]);
            b11 = (b[ie + i] + b[ib + i]) - (b[ic + i] + b[iF + i]);
            a20 = (a[ia + i] - a[id + i]) - 0.5 * a11;
            a21 = S60 * ((a[ie + i] - a[ib + i]) - (a[ic + i] - a[iF + i]));
            b20 = (b[ia + i] + b[id + i]) + 0.5 * b11;
            b21 = S60 * ((b[ie + i] + b[ib + i]) + (b[ic + i] + b[iF + i]));

            c[jd + j] = c3 * (a[ia + i] - a[id + i] + a11) - s3 * (b[ia + i] + b[id + i] - b11);
            d[jd + j] = s3 * (a[ia + i] - a[id + i] + a11) + c3 * (b[ia + i] + b[id + i] - b11);
            c[jb + j] = c1 * (a20 - b21) - s1 * (b20 - a21);
            d[jb + j] = s1 * (a20 - b21) + c1 * (b20 - a21);
            c[jf + j] = c5 * (a20 + b21) - s5 * (b20 + a21);
            d[jf + j] = s5 * (a20 + b21) + c5 * (b20 + a21);
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        ia += iink;
        ib += iink;
        ic += iink;
        id -= iink;
        ie -= iink;
        iF -= iink;
        jbase += jump;
      }
      if (ic > id) return 0;
    }
    ibase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = a[ib + i] + (a[ia + i] + a[ic + i]);
        c[jd + j] = b[ib + i] - (b[ia + i] + b[ic + i]);
        c[jb + j] = (S60 * (a[ia + i] - a[ic + i])) - (0.5 * (b[ia + i] + b[ic + i]) + b[ib + i]);
        c[jf + j] = -(S60 * (a[ia + i] - a[ic + i])) - (0.5 * (b[ia + i] + b[ic + i]) + b[ib + i]);
        c[jc + j] = S60 * (b[ic + i] - b[ia + i]) + (0.5 * (a[ia + i] + a[ic + i]) - a[ib + i]);
        c[je + j] = S60 * (b[ic + i] - b[ia + i]) - (0.5 * (a[ia + i] + a[ic + i]) - a[ib + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = (2.0 * (a[ia + i] + a[id + i])) + (2.0 * (a[ib + i] + a[ic + i]));
        c[jd + j] = (2.0 * (a[ia + i] - a[id + i])) - (2.0 * (a[ib + i] - a[ic + i]));
        c[jb + j] = (2.0 * (a[ia + i] - a[id + i]) + (a[ib + i] - a[ic + i])) - (D60 * (b[ib + i] + b[ic + i]));
        c[jf + j] = (2.0 * (a[ia + i] - a[id + i]) + (a[ib + i] - a[ic + i])) + (D60 * (b[ib + i] + b[ic + i]));
        c[jc + j] = (2.0 * (a[ia + i] + a[id + i]) - (a[ib + i] + a[ic + i])) - (D60 * (b[ib + i] - b[ic + i]));
        c[je + j] = (2.0 * (a[ia + i] + a[id + i]) - (a[ib + i] + a[ic + i])) + (D60 * (b[ib + i] - b[ic + i]));
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
rpassc_8(const double *a, double *c, long inc1, long inc2, long inc3, long inc4, long lot, long n, long ifac, long la)
{
  const long m = n / ifac;
  long iink = la * inc1;
  long jink = la * inc2;

  long ibase = 0;
  long jbase = 0;

  if (la != m) return 3;

  long i0 = 0;
  long i1 = i0 + iink;
  long i2 = i1 + iink;
  long i3 = i2 + iink;
  long i4 = i3 + iink;
  long i5 = i4 + iink;
  long i6 = i5 + iink;
  long i7 = i6 + iink;
  long j0 = 0;
  long j1 = j0 + jink;
  long j2 = j1 + jink;
  long j3 = j2 + jink;
  long j4 = j3 + jink;
  long j5 = j4 + jink;
  long j6 = j5 + jink;
  long j7 = j6 + jink;

  for (long l = 0; l < la; ++l)
  {
    long i = ibase;
    long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (long ijk = 0; ijk < lot; ++ijk)
    {
      const double a0p7 = a[i0 + i] + a[i7 + i];
      const double a0m7 = a[i0 + i] - a[i7 + i];
      const double a1p5 = a[i1 + i] + a[i5 + i];
      const double a1m5 = a[i1 + i] - a[i5 + i];
      const double a2p6 = a[i2 + i] + a[i6 + i];
      const double a2m6 = a[i2 + i] - a[i6 + i];

      const double p073 = a0p7 + a[i3 + i];
      const double m073 = a0p7 - a[i3 + i];

      const double p074 = 2.0 * (a0m7 + a[i4 + i]);
      const double m074 = 2.0 * (a0m7 - a[i4 + i]);

      const double p152 = M_SQRT2 * (a1m5 + a2p6);
      const double m152 = M_SQRT2 * (a1m5 - a2p6);

      c[j0 + j] = 2.0 * (p073 + a1p5);
      c[j4 + j] = 2.0 * (p073 - a1p5);
      c[j2 + j] = 2.0 * (m073 - a2m6);
      c[j6 + j] = 2.0 * (m073 + a2m6);

      c[j1 + j] = m074 + m152;
      c[j5 + j] = m074 - m152;
      c[j3 + j] = p074 - p152;
      c[j7 + j] = p074 + p152;
      i += inc3;
      j += inc4;
    }
    ibase += inc1;
    jbase += inc2;
  }

  return 0;
}

static int
rpassc(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
       long lot, long n, long ifac, long la)
{
  /*
     rpassc - performs one pass through data as part of multiple real fft (fourier synthesis) routine

     a      is first real input vector
     b      is equivalent to a + la * inc1
     c      is first real output vector
     d      is equivalent to c + ifac * la * inc2
     trigs  is a precalculated list of sines & cosines
     inc1   is the addressing increment for a
     inc2   is the addressing increment for c
     inc3   is the increment between input vectors a
     inc4   is the increment between output vectors c
     lot    is the number of vectors
     n      is the length of the vectors
     ifac   is the current factor of n
     la     is the product of previous factors
     ierr   is an error indicator:
               0 - pass completed without error
               2 - ifac not catered for
               3 - ifac only catered for if la=n/ifac
   */

  switch (ifac)
  {
    case 2: return rpassc_2(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 3: return rpassc_3(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 4: return rpassc_4(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 5: return rpassc_5(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 6: return rpassc_6(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 8: return rpassc_8(a, c, inc1, inc2, inc3, inc4, lot, n, ifac, la);
  }

  return 0;
}

static int
qpassc_2(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * iink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long j0 = 0;
  long i1 = i0 + iink;
  long j1 = j0 + inc2 * (m + m - la);
  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i] + a[i1 + i];
        c[j1 + j] = a[i0 + i] - a[i1 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    j0 += jink;
    jink += jink;
    j1 -= jink;
    ibase += jump;
    jump += jump + iink;

    if (j0 != j1)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        jbase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            c[j0 + j] = a[i0 + i] + c1 * a[i1 + i] + s1 * b[i1 + i];
            c[j1 + j] = a[i0 + i] - c1 * a[i1 + i] - s1 * b[i1 + i];
            d[j0 + j] = c1 * b[i1 + i] - s1 * a[i1 + i] + b[i0 + i];
            d[j1 + j] = c1 * b[i1 + i] - s1 * a[i1 + i] - b[i0 + i];
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        j0 += jink;
        j1 -= jink;
        ibase += jump;
      } /* End FORK */
      if (j0 > j1) return 0;
    } /* End (i0 != i1) */
    jbase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i];
        d[j1 + j] = -a[i1 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    const double z = 1.0 / n;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = z * (a[i0 + i] + a[i1 + i]);
        c[j1 + j] = z * (a[i0 + i] - a[i1 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
qpassc_3(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * iink;

  long ibase = 0;
  long jbase = 0;

  long ia = 0;
  long ib = ia + iink;
  long ic = ib + iink;

  long ja = 0;
  long jb = ja + inc2 * (m + m - la);
  long jc = jb;

  if (la != m) /* else 390 */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = a[ia + i] + a[ib + i] + a[ic + i];
        c[jb + j] = a[ia + i] - 0.5 * (a[ib + i] + a[ic + i]);
        d[jb + j] = S60 * (a[ic + i] - a[ib + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    ja += jink;
    jink += jink;
    jb += jink;
    jc -= jink;
    ibase += jump;
    jump += jump + iink;

    if (ja != jc) /* else  360 */
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        jbase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a1 = c1 * a[ib + i] + s1 * b[ib + i] + c2 * a[ic + i] + s2 * b[ic + i];
            const double b1 = c1 * b[ib + i] - s1 * a[ib + i] + c2 * b[ic + i] - s2 * a[ic + i];
            const double a2 = a[ia + i] - 0.5 * a1;
            const double b2 = b[ia + i] - 0.5 * b1;
            const double a3 = S60 * (c1 * a[ib + i] + s1 * b[ib + i] - c2 * a[ic + i] - s2 * b[ic + i]);
            const double b3 = S60 * (c1 * b[ib + i] - s1 * a[ib + i] - c2 * b[ic + i] + s2 * a[ic + i]);

            c[ja + j] = a[ia + i] + a1;
            d[ja + j] = b[ia + i] + b1;
            c[jb + j] = a2 + b3;
            d[jb + j] = b2 - a3;
            c[jc + j] = a2 - b3;
            d[jc + j] = -b2 - a3;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        ja += jink;
        jb += jink;
        jc -= jink;
        ibase += jump;
      } /* End FORK */
      if (ja > jc) return 0;
    } /* End (ia != ic) */
    jbase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = a[ia + i] + 0.5 * (a[ib + i] - a[ic + i]);
        d[ja + j] = -S60 * (a[ib + i] + a[ic + i]);
        c[jb + j] = a[ia + i] - a[ib + i] + a[ic + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    const double z = 1.0 / n;
    const double y = S60 / n;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[ja + j] = z * (a[ia + i] + a[ib + i] + a[ic + i]);
        c[jb + j] = z * (a[ia + i] - 0.5 * (a[ib + i] + a[ic + i]));
        d[jb + j] = y * (a[ic + i] - a[ib + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
qpassc_4(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * iink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long i1 = i0 + iink;
  long i2 = i1 + iink;
  long i3 = i2 + iink;
  long j0 = 0;
  long j1 = j0 + inc2 * (m + m - la);
  long j2 = j1 + inc2 * (m + m);
  long j3 = j1;

  if (la != m) /*else go to 490 */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p2 = a[i0 + i] + a[i2 + i];
        const double a1p3 = a[i1 + i] + a[i3 + i];

        c[j0 + j] = a0p2 + a1p3;
        c[j2 + j] = a0p2 - a1p3;

        c[j1 + j] = a[i0 + i] - a[i2 + i];
        d[j1 + j] = a[i3 + i] - a[i1 + i];
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    j0 += jink;
    jink += jink;
    j1 += jink;
    j2 -= jink;
    j3 -= jink;
    ibase += jump;
    jump += jump + iink;

    if (j1 != j2) /* else go to 460; */
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        jbase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a0 = a[i0 + i] + c2 * a[i2 + i] + s2 * b[i2 + i];
            const double a2 = a[i0 + i] - c2 * a[i2 + i] - s2 * b[i2 + i];
            const double b0 = b[i0 + i] + c2 * b[i2 + i] - s2 * a[i2 + i];
            const double b2 = b[i0 + i] - c2 * b[i2 + i] + s2 * a[i2 + i];

            const double a1 = c1 * a[i1 + i] + s1 * b[i1 + i] + c3 * a[i3 + i] + s3 * b[i3 + i];
            const double a3 = c1 * a[i1 + i] + s1 * b[i1 + i] - c3 * a[i3 + i] - s3 * b[i3 + i];
            const double b1 = c1 * b[i1 + i] - s1 * a[i1 + i] + c3 * b[i3 + i] - s3 * a[i3 + i];
            const double b3 = c1 * b[i1 + i] - s1 * a[i1 + i] - c3 * b[i3 + i] + s3 * a[i3 + i];

            c[j0 + j] = a0 + a1;
            c[j2 + j] = a0 - a1;
            d[j0 + j] = b0 + b1;
            d[j2 + j] = b1 - b0;
            c[j1 + j] = a2 + b3;
            c[j3 + j] = a2 - b3;
            d[j1 + j] = b2 - a3;
            d[j3 + j] = -b2 - a3;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        j0 += jink;
        j1 += jink;
        j2 -= jink;
        j3 -= jink;
        ibase += jump;
      } /* End FORK */
      if (j1 > j2) return 0;
    } /* End (i1 != i2) */
    jbase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        c[j0 + j] = a[i0 + i] + SQ2 * (a[i1 + i] - a[i3 + i]);
        c[j1 + j] = a[i0 + i] - SQ2 * (a[i1 + i] - a[i3 + i]);
        d[j0 + j] = -a[i2 + i] - SQ2 * (a[i1 + i] + a[i3 + i]);
        d[j1 + j] = a[i2 + i] - SQ2 * (a[i1 + i] + a[i3 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else /* (la != m) */
  {
    const double z = 1.0 / n;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p2 = a[i0 + i] + a[i2 + i];
        const double a1p3 = a[i1 + i] + a[i3 + i];

        c[j0 + j] = z * (a0p2 + a1p3);
        c[j2 + j] = z * (a0p2 - a1p3);
        c[j1 + j] = z * (a[i0 + i] - a[i2 + i]);
        d[j1 + j] = z * (a[i3 + i] - a[i1 + i]);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
qpassc_5(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * iink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long i1 = i0 + iink;
  long i2 = i1 + iink;
  long i3 = i2 + iink;
  long i4 = i3 + iink;
  long j0 = 0;
  long j1 = j0 + inc2 * (m + m - la);
  long j2 = j1 + inc2 * (m + m);
  long j3 = j2;
  long j4 = j1;

  if (la != m) /* else go to 590; */
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a1p4 = a[i1 + i] + a[i4 + i];
        const double a1m4 = a[i1 + i] - a[i4 + i];
        const double a2p3 = a[i2 + i] + a[i3 + i];
        const double a2m3 = a[i2 + i] - a[i3 + i];

        const double a025 = a[i0 + i] - 0.25 * (a1p4 + a2p3);
        const double aqrt = QT5 * (a1p4 - a2p3);

        c[j0 + j] = a[i0 + i] + a1p4 + a2p3;
        c[j1 + j] = a025 + aqrt;
        c[j2 + j] = a025 - aqrt;
        d[j1 + j] = -S72 * a1m4 - S36 * a2m3;
        d[j2 + j] = -S36 * a1m4 + S72 * a2m3;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    j0 += jink;
    jink += jink;
    j1 += jink;
    j2 += jink;
    j3 -= jink;
    j4 -= jink;
    ibase += jump;
    jump += jump + iink;

    if (j1 != j3)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const long ke = kd + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        const double c4 = trigs[ke];
        const double s4 = trigs[ke + 1];
        jbase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double a1p4 = c1 * a[i1 + i] + s1 * b[i1 + i] + c4 * a[i4 + i] + s4 * b[i4 + i];
            const double a1m4 = c1 * a[i1 + i] + s1 * b[i1 + i] - c4 * a[i4 + i] - s4 * b[i4 + i];
            const double a2p3 = c2 * a[i2 + i] + s2 * b[i2 + i] + c3 * a[i3 + i] + s3 * b[i3 + i];
            const double a2m3 = c2 * a[i2 + i] + s2 * b[i2 + i] - c3 * a[i3 + i] - s3 * b[i3 + i];
            const double b1p4 = c1 * b[i1 + i] - s1 * a[i1 + i] + c4 * b[i4 + i] - s4 * a[i4 + i];
            const double b1m4 = c1 * b[i1 + i] - s1 * a[i1 + i] - c4 * b[i4 + i] + s4 * a[i4 + i];
            const double b2p3 = c2 * b[i2 + i] - s2 * a[i2 + i] + c3 * b[i3 + i] - s3 * a[i3 + i];
            const double b2m3 = c2 * b[i2 + i] - s2 * a[i2 + i] - c3 * b[i3 + i] + s3 * a[i3 + i];

            const double a025 = a[i0 + i] - 0.25 * (a1p4 + a2p3);
            const double aqrt = QT5 * (a1p4 - a2p3);
            const double b025 = b[i0 + i] - 0.25 * (b1p4 + b2p3);
            const double bqrt = QT5 * (b1p4 - b2p3);

            const double a0pq = a025 + aqrt;
            const double a0mq = a025 - aqrt;
            const double b0pq = b025 + bqrt;
            const double b0mq = b025 - bqrt;

            const double asps = S72 * a1m4 + S36 * a2m3;
            const double asms = S36 * a1m4 - S72 * a2m3;
            const double bsps = S72 * b1m4 + S36 * b2m3;
            const double bsms = S36 * b1m4 - S72 * b2m3;

            c[j0 + j] = a[i0 + i] + a1p4 + a2p3;
            c[j1 + j] = a0pq + bsps;
            c[j2 + j] = a0mq + bsms;
            c[j3 + j] = a0mq - bsms;
            c[j4 + j] = a0pq - bsps;
            d[j0 + j] = b[i0 + i] + b1p4 + b2p3;
            d[j1 + j] = b0pq - asps;
            d[j2 + j] = b0mq - asms;
            d[j3 + j] = -b0mq - asms;
            d[j4 + j] = -b0pq - asps;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        j0 += jink;
        j1 += jink;
        j2 += jink;
        j3 -= jink;
        j4 -= jink;
        ibase += jump;
      } /* End FORK */
      if (j1 > j3) return 0;
    } /* End (jb != jd) */
    jbase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a1p4 = a[i1 + i] + a[i4 + i];
        const double a1m4 = a[i1 + i] - a[i4 + i];
        const double a2p3 = a[i2 + i] + a[i3 + i];
        const double a2m3 = a[i2 + i] - a[i3 + i];

        const double a025 = a[i0 + i] + 0.25 * (a1m4 - a2m3);
        const double aqrt = QT5 * (a1m4 + a2m3);

        c[j0 + j] = a025 + aqrt;
        c[j1 + j] = a025 - aqrt;
        c[j2 + j] = a[i0 + i] - a1m4 + a2m3;
        d[j0 + j] = -S36 * a1p4 - S72 * a2p3;
        d[j1 + j] = -S72 * a1p4 + S36 * a2p3;

        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else
  {
    const double z = 1.0 / n;
    const double y = QT5 / n;
    const double x = S36 / n;
    const double w = S72 / n;

    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a1p4 = a[i1 + i] + a[i4 + i];
        const double a1m4 = a[i1 + i] - a[i4 + i];
        const double a2p3 = a[i2 + i] + a[i3 + i];
        const double a2m3 = a[i2 + i] - a[i3 + i];

        const double a025 = z * (a[i0 + i] - 0.25 * (a1p4 + a2p3));
        const double aqrt = y * (a1p4 - a2p3);

        c[j0 + j] = z * (a[i0 + i] + a1p4 + a2p3);
        c[j1 + j] = a025 + aqrt;
        c[j2 + j] = a025 - aqrt;
        d[j1 + j] = -w * a1m4 - x * a2m3;
        d[j2 + j] = w * a2m3 - x * a1m4;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
qpassc_6(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
         long lot, long n, long ifac, long la)
{
  const long kstop = (n - ifac) / (2 * ifac);
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;
  long jump = (ifac - 1) * iink;

  long ibase = 0;
  long jbase = 0;

  long i0 = 0;
  long i1 = i0 + iink;
  long i2 = i1 + iink;
  long i3 = i2 + iink;
  long i4 = i3 + iink;
  long i5 = i4 + iink;
  long j0 = 0;
  long j1 = j0 + inc2 * (m + m - la);
  long j2 = j1 + inc2 * (m + m);
  long j3 = j2 + inc2 * (m + m);
  long j4 = j2;
  long j5 = j1;

  if (la != m)
  {
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p3 = a[i0 + i] + a[i3 + i];
        const double a0m3 = a[i0 + i] - a[i3 + i];
        const double a1p4 = a[i1 + i] + a[i4 + i];
        const double a1m4 = a[i1 + i] - a[i4 + i];
        const double a2p5 = a[i2 + i] + a[i5 + i];
        const double a2m5 = a[i2 + i] - a[i5 + i];

        c[j0 + j] = a0p3 + a1p4 + a2p5;
        c[j3 + j] = a0m3 + a2m5 - a1m4;

        c[j1 + j] = a0m3 - 0.5 * (a2m5 - a1m4);
        c[j2 + j] = a0p3 - 0.5 * (a1p4 + a2p5);

        d[j1 + j] = S60 * (-a2m5 - a1m4);
        d[j2 + j] = S60 * (a2p5 - a1p4);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
    j0 += jink;
    jink += jink;
    j1 += jink;
    j2 += jink;
    j3 -= jink;
    j4 -= jink;
    j5 -= jink;
    ibase += jump;
    jump += jump + iink;

    if (j2 != j3)
    {
      for (long k = la; k <= kstop; k += la)
      {
        const long kb = k + k;
        const long kc = kb + kb;
        const long kd = kc + kb;
        const long ke = kd + kb;
        const long kf = ke + kb;
        const double c1 = trigs[kb];
        const double s1 = trigs[kb + 1];
        const double c2 = trigs[kc];
        const double s2 = trigs[kc + 1];
        const double c3 = trigs[kd];
        const double s3 = trigs[kd + 1];
        const double c4 = trigs[ke];
        const double s4 = trigs[ke + 1];
        const double c5 = trigs[kf];
        const double s5 = trigs[kf + 1];
        jbase = 0;
        for (long l = 0; l < la; ++l)
        {
          long i = ibase;
          long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
          for (long ijk = 0; ijk < lot; ++ijk)
          {
            const double ab1a = c1 * a[i1 + i] + s1 * b[i1 + i];
            const double ab1b = c1 * b[i1 + i] - s1 * a[i1 + i];
            const double ab2a = c2 * a[i2 + i] + s2 * b[i2 + i];
            const double ab2b = c2 * b[i2 + i] - s2 * a[i2 + i];
            const double ab3a = c3 * a[i3 + i] + s3 * b[i3 + i];
            const double ab3b = c3 * b[i3 + i] - s3 * a[i3 + i];
            const double ab4a = c4 * a[i4 + i] + s4 * b[i4 + i];
            const double ab4b = c4 * b[i4 + i] - s4 * a[i4 + i];
            const double ab5a = c5 * a[i5 + i] + s5 * b[i5 + i];
            const double ab5b = c5 * b[i5 + i] - s5 * a[i5 + i];

            const double a1p4 = ab1a + ab4a;
            const double a1m4 = ab1a - ab4a;
            const double a2p5 = ab2a + ab5a;
            const double a2m5 = ab2a - ab5a;

            const double b1p4 = ab1b + ab4b;
            const double b1m4 = ab1b - ab4b;
            const double b2p5 = ab2b + ab5b;
            const double b2m5 = ab2b - ab5b;

            const double ap05 = a[i0 + i] + ab3a - 0.5 * (a1p4 + a2p5);
            const double bp05 = b[i0 + i] + ab3b - 0.5 * (b1p4 + b2p5);
            const double am05 = a[i0 + i] - ab3a - 0.5 * (a2m5 - a1m4);
            const double bm05 = -b[i0 + i] + ab3b - 0.5 * (b1m4 - b2m5);

            const double ap60 = S60 * (a2p5 - a1p4);
            const double bp60 = S60 * (b2p5 - b1p4);
            const double am60 = S60 * (-a2m5 - a1m4);
            const double bm60 = S60 * (-b2m5 - b1m4);

            c[j0 + j] = a[i0 + i] + ab3a + a1p4 + a2p5;
            d[j0 + j] = b[i0 + i] + ab3b + b1p4 + b2p5;
            c[j1 + j] = am05 - bm60;
            d[j1 + j] = am60 - bm05;
            c[j2 + j] = ap05 - bp60;
            d[j2 + j] = ap60 + bp05;
            c[j3 + j] = a[i0 + i] - ab3a - a1m4 + a2m5;
            d[j3 + j] = -b[i0 + i] + ab3b + b1m4 - b2m5;
            c[j4 + j] = ap05 + bp60;
            d[j4 + j] = ap60 - bp05;
            c[j5 + j] = am05 + bm60;
            d[j5 + j] = am60 + bm05;
            i += inc3;
            j += inc4;
          }
          ibase += inc1;
          jbase += inc2;
        }
        j0 += jink;
        j1 += jink;
        j2 += jink;
        j3 -= jink;
        j4 -= jink;
        j5 -= jink;
        ibase += jump;
      }
      if (j2 > j3) return 0;
    }
    jbase = 0;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a1p5 = a[i1 + i] + a[i5 + i];
        const double a1m5 = a[i1 + i] - a[i5 + i];
        const double a2p4 = a[i2 + i] + a[i4 + i];
        const double a2m4 = a[i2 + i] - a[i4 + i];

        c[j0 + j] = a[i0 + i] + 0.5 * a2m4 + S60 * a1m5;
        d[j0 + j] = -a[i3 + i] - 0.5 * a1p5 - S60 * a2p4;
        c[j1 + j] = a[i0 + i] - a2m4;
        d[j1 + j] = a[i3 + i] - a1p5;
        c[j2 + j] = a[i0 + i] + 0.5 * a2m4 - S60 * a1m5;
        d[j2 + j] = -a[i3 + i] - 0.5 * a1p5 + S60 * a2p4;
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }
  else
  {
    const double z = 1.0 / n;
    const double y = S60 / n;
    for (long l = 0; l < la; ++l)
    {
      long i = ibase;
      long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
      for (long ijk = 0; ijk < lot; ++ijk)
      {
        const double a0p3 = a[i0 + i] + a[i3 + i];
        const double a0m3 = a[i0 + i] - a[i3 + i];
        const double a1p4 = a[i1 + i] + a[i4 + i];
        const double a1m4 = a[i1 + i] - a[i4 + i];
        const double a2p5 = a[i2 + i] + a[i5 + i];
        const double a2m5 = a[i2 + i] - a[i5 + i];

        c[j0 + j] = z * (a0p3 + a1p4 + a2p5);
        c[j3 + j] = z * (a0m3 + a2m5 - a1m4);

        c[j1 + j] = z * (a0m3 - 0.5 * (a2m5 - a1m4));
        c[j2 + j] = z * (a0p3 - 0.5 * (a1p4 + a2p5));

        d[j1 + j] = y * (-a2m5 - a1m4);
        d[j2 + j] = y * (a2p5 - a1p4);
        i += inc3;
        j += inc4;
      }
      ibase += inc1;
      jbase += inc2;
    }
  }

  return 0;
}

static int
qpassc_8(const double *a, double *c, long inc1, long inc2, long inc3, long inc4, long lot, long n, long ifac, long la)
{
  const long m = n / ifac;
  const long iink = la * inc1;
  long jink = la * inc2;

  long ibase = 0;
  long jbase = 0;

  if (la != m) return 3;

  long i0 = 0;
  long i1 = i0 + iink;
  long i2 = i1 + iink;
  long i3 = i2 + iink;
  long i4 = i3 + iink;
  long i5 = i4 + iink;
  long i6 = i5 + iink;
  long i7 = i6 + iink;
  long j0 = 0;
  long j1 = j0 + jink;
  long j2 = j1 + jink;
  long j3 = j2 + jink;
  long j4 = j3 + jink;
  long j5 = j4 + jink;
  long j6 = j5 + jink;
  long j7 = j6 + jink;
  const double z = 1.0 / n;
  const double y = SQ2 / n;

  for (long l = 0; l < la; ++l)
  {
    long i = ibase;
    long j = jbase;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (long ijk = 0; ijk < lot; ++ijk)
    {
      const double a0p4 = a[i0 + i] + a[i4 + i];
      const double a0m4 = a[i0 + i] - a[i4 + i];
      const double a1p5 = a[i1 + i] + a[i5 + i];
      const double a1m5 = a[i1 + i] - a[i5 + i];
      const double a2p6 = a[i2 + i] + a[i6 + i];
      const double a2m6 = a[i2 + i] - a[i6 + i];
      const double a3p7 = a[i3 + i] + a[i7 + i];
      const double a3m7 = a[i3 + i] - a[i7 + i];

      c[j0 + j] = z * (a0p4 + a1p5 + a2p6 + a3p7);
      c[j7 + j] = z * (a0p4 - a1p5 + a2p6 - a3p7);

      c[j3 + j] = z * (a0p4 - a2p6);
      c[j4 + j] = z * (a3p7 - a1p5);

      c[j1 + j] = z * a0m4 + y * (a1m5 - a3m7);
      c[j5 + j] = z * a0m4 - y * (a1m5 - a3m7);
      c[j2 + j] = -z * a2m6 - y * (a1m5 + a3m7);
      c[j6 + j] = z * a2m6 - y * (a1m5 + a3m7);
      i += inc3;
      j += inc4;
    }
    ibase += inc1;
    jbase += inc2;
  }

  return 0;
}

static int
qpassc(const double *a, const double *b, double *c, double *d, const double *trigs, long inc1, long inc2, long inc3, long inc4,
       long lot, long n, long ifac, long la)
{
  /*
     qpassc - performs one pass through data as part of multiple real fft (fourier analysis) routine

     a      is first real input vector
     b      is equivalent to a + ifac * la * inc1
     c      is first real output vector;
     d      is equivalent to c + la * inc2
     trigs  is a precalculated list of sines & cosines
     inc1   is the addressing increment for a
     inc2   is the addressing increment for c
     inc3   is the increment between input vectors a
     inc4   is the increment between output vectors c
     lot    is the number of vectors
     n      is the length of the vectors
     ifac   is the current factor of n
     la     is the product of previous factors
   */

  switch (ifac)
  {
    case 2: return qpassc_2(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 3: return qpassc_3(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 4: return qpassc_4(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 5: return qpassc_5(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 6: return qpassc_6(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot, n, ifac, la);
    case 8: return qpassc_8(a, c, inc1, inc2, inc3, inc4, lot, n, ifac, la);
  }

  return 0;
}

/* ====================== */
/* Fast Fourier Transform */
/* ====================== */
void
fc2gp(const double *restrict trig, const long *restrict ifax, const double *restrict fc, double *restrict gp, long nlat, long nlon,
      long nlev, long nfc)
{
  /* fc2gp performs fourier to gridpoint transforms using           */
  /* multiple fast fourier transform of length nlon                 */
  /*                                                                */
  /* fc   - real array of fourier coefficients fc[nlev][nfc][nlat]  */
  /* gp   - real array of gridpoints           gp[nlev][nlat][nlon] */
  /* nlat - Number of latitudes                                     */
  /* nlon - Number of longitudes                                    */
  /* nlev - Number of levels                                        */
  /* nfc  - Number of fourier coefficients on 1 latitude            */

  /* x(j) = sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/nlon))             */
  /*        where c(k) = a(k) + i*b(k) and c(n-k) = a(k)-i*b(k)     */

  if (ifax[9] != nlon) fprintf(stderr, "fc2gp: wrong initialization!\n");

  const long nfax = ifax[0];

  const long jump = (nlon + 2);
  const long lot = nlev * nlat;

  const long nx = (nlon % 2 == 1) ? nlon : nlon + 1;
  const long nblox = 1 + (lot - 1) / NFFT;
  long nvexx = lot - (nblox - 1) * NFFT;
  const long nvex0 = nvexx;

  const long nthmax = (nblox < Threading::ompNumMaxThreads) ? nblox : Threading::ompNumMaxThreads;
  const long nvals = lot * jump;

  Varray<double> wfc(nvals);
  Varray2D<double> wgp2d(nthmax, Varray<double>(nvals));

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long lat = 0; lat < nlat; ++lat)
  {
    for (long lev = 0; lev < nlev; ++lev)
    {
      double *restrict wfcx = wfc.data() + jump * (lat + lev * nlat);
      const double *restrict fcx = fc + (lat + lev * nlat * nfc);
      for (long fou = 0; fou < nfc; ++fou) wfcx[fou] = fcx[fou * nlat];
      for (long fou = nfc; fou < jump; ++fou) wfcx[fou] = 0.0;
    }
  }

  std::vector<long> istartv(nblox);

  long istart0 = 0;
  for (long nb = 0; nb < nblox; nb++)
  {
    istartv[nb] = istart0;
    istart0 = istart0 + nvexx * jump;
    nvexx = NFFT;
  }

// printf("nblox %ld nvex0 %ld  lot %ld\n",nblox, nvex0, lot);
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long nb = 0; nb < nblox; nb++)
  {
    auto ompthID = cdo_omp_get_thread_num();
    double *restrict wgp = &wgp2d[ompthID][0];

    long istart = istartv[nb];
    long nvex = (nb == 0) ? nvex0 : NFFT;

    long i = istart;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (long j = 0; j < nvex; ++j)
    {
      wfc[i + 1] = 0.5 * wfc[i];
      i += jump;
    }
    if (nlon % 2 != 1)
    {
      i = istart + nlon;
      for (long j = 0; j < nvex; ++j)
      {
        wfc[i] = 0.5 * wfc[i];
        i += jump;
      }
    }

    long ia = istart + 1;
    long la = 1;
    auto wfcx = wfc.data();
    for (long k = 0; k < nfax; ++k)
    {
      const long ifac = ifax[k + 1];

      if (k & 1)
        rpassc(wgp, wgp + la, wfcx + ia, wfcx + ia + ifac * la, trig, 1, 1, nx, jump, nvex, nlon, ifac, la);
      else
        rpassc(wfcx + ia, wfcx + ia + la, wgp, wgp + ifac * la, trig, 1, 1, jump, nx, nvex, nlon, ifac, la);

      la *= ifac;
      ia = istart;
    }

    // If necessary, copy results back to a

    if (nfax % 2 != 0)
    {
      long ibase = 0;
      long jbase = ia;
      for (long jj = 0; jj < nvex; ++jj)
      {
        for (long ii = 0; ii < nlon; ++ii) wfc[jbase + ii] = wgp[ibase + ii];
        ibase = ibase + nx;
        jbase = jbase + jump;
      }
    }

    // Fill in zeros at end

    long ix = istart + nlon;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (long j = 0; j < nvex; ++j)
    {
      wfc[ix] = 0.0;
      wfc[ix + 1] = 0.0;
      ix += jump;
    }
  }

  auto const &wpt = wfc;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long j = 0; j < lot; ++j)
    for (long lon = 0; lon < nlon; ++lon) gp[lon + j * nlon] = wpt[lon + j * jump];
}

void
gp2fc(const double *trig, const long *ifax, const double *restrict gp, double *restrict fc, long nlat, long nlon, long nlev,
      long nfc)
{
  /* gp2fc performs gridpoint to fourier transforms using           */
  /* multiple fast fourier transform of length nlon                 */
  /*                                                                */
  /* gp   - real array of gridpoints           gp[nlev][nlat][nlon] */
  /* fc   - real array of fourier coefficients fc[nlev][nfc][nlat]  */
  /* nlat - Number of latitudes                                     */
  /* nlon - Number of longitudes                                    */
  /* nlev - Number of levels                                        */
  /* nfc  - Number of fourier coefficients on 1 latitude            */

  /* a(k) =  (1/n) * sum(j=0,...,n-1)(x(j) * cos(2*j*k*pi/n))       */
  /* b(k) = -(1/n) * sum(j=0,...,n-1)(x(j) * sin(2*j*k*pi/n))       */

  if (ifax[9] != nlon) fprintf(stderr, "gp2fc: wrong initialization!\n");

  const long nfax = ifax[0];

  const long jump = (nlon + 2);
  const long lot = nlev * nlat;

  Varray<double> wfc(lot * jump), wgp(lot * jump);

  long rix = 0;
  long wix = 0;
  for (long j = 0; j < lot; ++j)
  {
    for (long lon = 0; lon < nlon; ++lon) wgp[wix + lon] = gp[rix + lon];
    wgp[wix + nlon] = 0.0;
    wgp[wix + nlon + 1] = 0.0;
    rix += nlon;
    wix += jump;
  }

  const long nx = (nlon % 2 == 1) ? nlon : nlon + 1;
  const long nblox = 1 + (lot - 1) / NFFT;
  long nvex = lot - (nblox - 1) * NFFT;

  long istart = 0;
  auto wfcx = wfc.data();
  auto wgpx = wgp.data();
  for (long nb = 0; nb < nblox; nb++)
  {
    long ia = istart;
    long la = nlon;
    for (long k = 0; k < nfax; ++k)
    {
      long ifac = ifax[nfax - k];
      la /= ifac;
      if (k & 1)
        qpassc(wfcx, wfcx + ifac * la, wgpx + ia, wgpx + ia + la, trig, 1, 1, nx, jump, nvex, nlon, ifac, la);
      else
        qpassc(wgpx + ia, wgpx + ia + ifac * la, wfcx, wfcx + la, trig, 1, 1, jump, nx, nvex, nlon, ifac, la);
      ia = istart + 1;
    }

    // If necessary, copy results back to a

    if (nfax % 2 != 0)
    {
      long ibase = 0;
      long jbase = ia;
      for (long jj = 0; jj < nvex; ++jj)
      {
        for (long ii = 0; ii < nlon; ++ii) wgp[jbase + ii] = wfc[ibase + ii];
        ibase = ibase + nx;
        jbase = jbase + jump;
      }
    }

    // Shift a(0) & fill in zero imag parts

    long ix = istart;
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (long j = 0; j < nvex; ++j)
    {
      wgp[ix] = wgp[ix + 1];
      wgp[ix + 1] = 0.0;
      ix += jump;
    }

    if (nlon % 2 != 1)
    {
      long iz = istart + (nlon + 1);
      for (long j = 0; j < nvex; ++j)
      {
        wgp[iz] = 0.0;
        iz += jump;
      }
    }

    istart = istart + nvex * jump;
    nvex = NFFT;
  }

  for (long lev = 0; lev < nlev; ++lev)
  {
    for (long lat = 0; lat < nlat; ++lat)
    {
      rix = jump * (lat + lev * nlat);
      wix = lat + lev * nlat * nfc;
      const auto wpt = wgp.data() + rix;
      auto fct = fc + wix;
      fct[0] = wpt[0];
      fct[nlat] = 0.0;
      for (long fou = 2; fou < nfc; ++fou) fct[fou * nlat] = wpt[fou];
    }
  }
}
