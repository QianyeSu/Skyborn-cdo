/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cmath>
#include <cstring>

#include "cdo_options.h"
#include "cdo_omp.h"
#include "constants.h"
#include "varray.h"
#include "gaussian_latitudes.h"

// phcs - Compute values of Legendre polynomials and their meridional derivatives
static void
phcs(bool needHnm, double *pnm, double *hnm, long waves, double pmu, double *ztemp1, double *ztemp2)
{
  long jnmjk;
  double zcospar, zsinpar, zcosfak, zsinfak;
  double zq, zwm2, zw, zwq, zq2m1, zwm2q2, z2q2, zcnm, zdnm, zenm;

  auto twowaves = waves << 1;

  auto zcos2 = std::sqrt(1.0 - pmu * pmu);
  auto lat = std::acos(pmu);
  auto zan = 1.0;

  ztemp1[0] = 0.5;

  for (long jn = 1; jn < twowaves; jn++)
  {
    auto zsqp = 1.0 / std::sqrt((double) (jn + jn * jn));
    zan *= std::sqrt(1.0 - 1.0 / (4 * jn * jn));

    zcospar = std::cos(lat * jn);
    zsinpar = std::sin(lat * jn) * jn * zsqp;
    zcosfak = 1.0;

    for (long jk = 2; jk < jn; jk += 2)
    {
      jnmjk = jn - jk;
      zcosfak *= (jk - 1.0) * (jn + jnmjk + 2.0) / (jk * (jn + jnmjk + 1.0));
      zsinfak = zcosfak * zsqp * (double) jnmjk;
      zcospar += zcosfak * std::cos(lat * jnmjk);
      zsinpar += zsinfak * std::sin(lat * jnmjk);
    }

    // code for jk == jn

    if ((jn & 1) == 0)
    {
      zcosfak *= (double) ((jn - 1) * (jn + 2)) / (double) (jn * (jn + 1));
      zcospar += zcosfak * 0.5;
    }
    ztemp1[jn] = zan * zcospar;
    ztemp2[jn - 1] = zan * zsinpar;
  }

  std::memcpy(pnm, ztemp1, waves * sizeof(double));
  pnm += waves;
  std::memcpy(pnm, ztemp2, waves * sizeof(double));
  pnm += waves;

  if (needHnm)
  {
    hnm[0] = 0.0;
    for (long jn = 1; jn < waves; jn++)
      hnm[jn] = jn * (pmu * ztemp1[jn] - std::sqrt((jn + jn + 1.0) / (jn + jn - 1.0)) * ztemp1[jn - 1]);

    hnm += waves;

    hnm[0] = pmu * ztemp2[0];

    for (long jn = 1; jn < waves; jn++)
      hnm[jn] = (jn + 1) * pmu * ztemp2[jn]
                - std::sqrt(((jn + jn + 3.0) * ((jn + 1) * (jn + 1) - 1.0)) / (jn + jn + 1.0)) * ztemp2[jn - 1];

    hnm += waves;
  }

  for (long jm = 2; jm < waves; jm++)
  {
    pnm[0] = std::sqrt(1.0 + 1.0 / (jm + jm)) * zcos2 * ztemp2[0];
    if (needHnm) hnm[0] = jm * pmu * pnm[0];
#if defined(CRAY)
#pragma _CRI novector
#endif
#if defined(__uxp__)
#pragma loop scalar
#endif
    for (long jn = 1; jn < (twowaves - jm); jn++)
    {
      zq = jm + jm + jn - 1;
      zwm2 = zq + jn;
      zw = zwm2 + 2.0;
      zwq = zw * zq;
      zq2m1 = zq * zq - 1.0;
      zwm2q2 = zwm2 * zq2m1;
      z2q2 = zq2m1 * 2;
      zcnm = std::sqrt((zwq * (zq - 2.0)) / (zwm2q2 - z2q2));
      zdnm = std::sqrt((zwq * (jn + 1.0)) / zwm2q2);
      zenm = std::sqrt(zw * jn / ((zq + 1.0) * zwm2));
      pnm[jn] = zcnm * ztemp1[jn] - pmu * (zdnm * ztemp1[jn + 1] - zenm * pnm[jn - 1]);
      if (needHnm) hnm[jn] = (jm + jn) * pmu * pnm[jn] - std::sqrt(zw * jn * (zq + 1) / zwm2) * pnm[jn - 1];
    }
    std::memcpy(ztemp1, ztemp2, twowaves * sizeof(double));
    std::memcpy(ztemp2, pnm, twowaves * sizeof(double));
    pnm += waves;
    if (needHnm) hnm += waves;
  }
}

void
after_legini_full(long ntr, long nlat, double *restrict poli, double *restrict pold, double *restrict pdev, double *restrict pol2,
                  double *restrict pol3, double *restrict coslat)
{
  auto dimsp = (ntr + 1) * (ntr + 2);
  auto waves = ntr + 1;
  auto twowaves = waves << 1;

  Varray<double> gmu(nlat), gwt(nlat);
  gaussian_latitudes(nlat, gmu.data(), gwt.data());

  auto needHnm = (pdev != nullptr) || (pol2 != nullptr);

#ifdef _OPENMP
  Varray2D<double> pnm_2(Threading::ompNumMaxThreads);
  Varray2D<double> hnm_2(Threading::ompNumMaxThreads);
  Varray2D<double> ztemp1_2(Threading::ompNumMaxThreads);
  Varray2D<double> ztemp2_2(Threading::ompNumMaxThreads);
  for (long i = 0; i < Threading::ompNumMaxThreads; ++i) pnm_2[i].resize(dimsp);
  if (needHnm)
    for (long i = 0; i < Threading::ompNumMaxThreads; ++i) hnm_2[i].resize(dimsp);
  for (long i = 0; i < Threading::ompNumMaxThreads; ++i) ztemp1_2[i].resize(twowaves);
  for (long i = 0; i < Threading::ompNumMaxThreads; ++i) ztemp2_2[i].resize(twowaves);
#else
  Varray<double> pnm(dimsp), hnm, ztemp1(twowaves), ztemp2(twowaves);
  if (needHnm) hnm.resize(dimsp);
#endif

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long jgl = 0; jgl < nlat; ++jgl)
  {
#ifdef _OPENMP
    auto ompthID = cdo_omp_get_thread_num();
    auto &pnm = pnm_2[ompthID];
    auto &hnm = hnm_2[ompthID];
    auto &ztemp1 = ztemp1_2[ompthID];
    auto &ztemp2 = ztemp2_2[ompthID];
#endif
    auto gmusq = 1.0 - gmu[jgl] * gmu[jgl];
    coslat[jgl] = std::sqrt(gmusq);

    phcs(needHnm, pnm.data(), hnm.data(), waves, gmu[jgl], ztemp1.data(), ztemp2.data());

    auto zgwt = gwt[jgl];
    auto zrafgmusqr = 1.0 / (PlanetRadiusDefault * gmusq);
    auto zradsqrtgmusqr = 1.0 / (-PlanetRadiusDefault * std::sqrt(gmusq));

    auto jsp = jgl;
    for (long jm = 0; jm < waves; ++jm)
      for (long jn = 0; jn < (waves - jm); ++jn)
      {
        if (poli) poli[jsp] = pnm[jm * waves + jn] * 2.0;
        if (pold) pold[jsp] = pnm[jm * waves + jn] * zgwt;
        if (pdev) pdev[jsp] = hnm[jm * waves + jn] * 2.0 * zradsqrtgmusqr;
        if (pol2) pol2[jsp] = hnm[jm * waves + jn] * zgwt * zrafgmusqr;
        if (pol3) pol3[jsp] = pnm[jm * waves + jn] * zgwt * jm * zrafgmusqr;
        jsp += nlat;
      }
  }
}

static inline void
sp2fc_kernel(long nlat, const double *pol, const double *sal, double *restrict far, double *restrict fai)
{
  auto sar = sal[0];
  auto sai = sal[1];
  for (long lat = 0; lat < nlat; lat++)
  {
    far[lat] += pol[lat] * sar;
    fai[lat] += pol[lat] * sai;
  }
}

void
sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long ntp1 = nt + 1;
  long nsp2 = (nt + 1) * (nt + 2);

  std::vector<long> cumindex(ntp1);
  cumindex[0] = 0;
  for (long jmm = 1; jmm < ntp1; jmm++) cumindex[jmm] = cumindex[jmm - 1] + (ntp1 - jmm + 1);

  if (nlev >= Threading::ompNumMaxThreads)
  {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (long lev = 0; lev < nlev; lev++)
    {
      auto sal = sa + lev * nsp2;
      auto fal = fa + lev * nfc * nlat;
      memset(fal, 0, nfc * nlat * sizeof(double));

      for (long jmm = 0; jmm < ntp1; jmm++)
      {
        auto polt = poli + cumindex[jmm] * nlat;
        auto salt = sal + cumindex[jmm] * 2;
        auto far = fal + jmm * 2 * nlat;
        auto fai = far + nlat;
        for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { sp2fc_kernel(nlat, polt + jfc * nlat, salt + jfc * 2, far, fai); }
      }
    }
  }
  else
  {
    for (long lev = 0; lev < nlev; lev++)
    {
      auto sal = sa + lev * nsp2;
      auto fal = fa + lev * nfc * nlat;
      memset(fal, 0, nfc * nlat * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
      for (long jmm = 0; jmm < ntp1; jmm++)
      {
        auto polt = poli + cumindex[jmm] * nlat;
        auto salt = sal + cumindex[jmm] * 2;
        auto far = fal + jmm * 2 * nlat;
        auto fai = far + nlat;
        for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { sp2fc_kernel(nlat, polt + jfc * nlat, salt + jfc * 2, far, fai); }
      }
    }
  }
}

static inline void
fc2sp_kernel(long nlat, const double *pol, const double *far, const double *fai, double *sal)
{
  double sar = 0.0;
  double sai = 0.0;
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sar) reduction(+ : sai)
#endif
  for (long lat = 0; lat < nlat; lat++)
  {
    sar += pol[lat] * far[lat];
    sai += pol[lat] * fai[lat];
  }
  sal[0] = sar;
  sal[1] = sai;
}

void
fc2sp(const double *fa, double *sa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long ntp1 = nt + 1;
  long nsp2 = (nt + 1) * (nt + 2);

  std::vector<long> cumindex(ntp1);
  cumindex[0] = 0;
  for (long jmm = 1; jmm < ntp1; jmm++) cumindex[jmm] = cumindex[jmm - 1] + (ntp1 - jmm + 1);

  if (nlev >= Threading::ompNumMaxThreads)
  {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (long lev = 0; lev < nlev; lev++)
    {
      auto fal = fa + lev * nfc * nlat;
      auto sal = sa + lev * nsp2;

      for (long jmm = 0; jmm < ntp1; jmm++)
      {
        auto polt = poli + cumindex[jmm] * nlat;
        auto salt = sal + cumindex[jmm] * 2;
        auto far = fal + jmm * 2 * nlat;
        auto fai = far + nlat;
        for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { fc2sp_kernel(nlat, polt + jfc * nlat, far, fai, salt + jfc * 2); }
      }
    }
  }
  else
  {
    for (long lev = 0; lev < nlev; lev++)
    {
      auto fal = fa + lev * nfc * nlat;
      auto sal = sa + lev * nsp2;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
      for (long jmm = 0; jmm < ntp1; jmm++)
      {
        auto polt = poli + cumindex[jmm] * nlat;
        auto salt = sal + cumindex[jmm] * 2;
        auto far = fal + jmm * 2 * nlat;
        auto fai = far + nlat;
        for (long jfc = 0; jfc < (ntp1 - jmm); jfc++) { fc2sp_kernel(nlat, polt + jfc * nlat, far, fai, salt + jfc * 2); }
      }
    }
  }
}

/* ======================================== */
/* Convert Spectral Array to new truncation */
/* ======================================== */

void
sp2sp(const double *arrayIn, long truncIn, double *arrayOut, long truncOut)
{
  if (truncOut <= truncIn)
  {
    for (long n = 0; n <= truncOut; ++n)
    {
      for (long m = n; m <= truncOut; ++m)
      {
        *arrayOut++ = *arrayIn++;
        *arrayOut++ = *arrayIn++;
      }
      arrayIn += 2 * (truncIn - truncOut);
    }
  }
  else
  {
    for (long n = 0; n <= truncIn; ++n)
    {
      for (long m = n; m <= truncIn; ++m)
      {
        *arrayOut++ = *arrayIn++;
        *arrayOut++ = *arrayIn++;
      }
      for (long m = truncIn + 1; m <= truncOut; ++m)
      {
        *arrayOut++ = 0.0;
        *arrayOut++ = 0.0;
      }
    }
    for (long n = truncIn + 1; n <= truncOut; ++n)
      for (long m = n; m <= truncOut; ++m)
      {
        *arrayOut++ = 0.0;
        *arrayOut++ = 0.0;
      }
  }
}

/* ======================================== */
/* Cut spectral wave numbers                */
/* ======================================== */

void
spcut(const double *arrayIn, double *arrayOut, long trunc, const int *waves)
{
  for (long n = 0; n <= trunc; ++n)
  {
    for (long m = n; m <= trunc; ++m)
    {
      if (waves[m])
      {
        *arrayOut++ = *arrayIn++;
        *arrayOut++ = *arrayIn++;
      }
      else
      {
        *arrayOut++ = 0.0;
        *arrayOut++ = 0.0;
        arrayIn++;
        arrayIn++;
      }
    }
  }
}
