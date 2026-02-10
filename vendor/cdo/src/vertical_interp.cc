/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"  // restrict
#endif

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "constants.h"
#include "vertical_interp.h"

constexpr auto SCALEHEIGHT = -7000.0;
constexpr auto SCALESLP = 101325.0;

void
height_to_pressure(const double *heightLevels, double *pressureLevels, int numLevels)
{
  for (int k = 0; k < numLevels; ++k)
  {
    // heightLevels: height in meters
    pressureLevels[k] = SCALESLP * std::exp(heightLevels[k] / SCALEHEIGHT);
  }
}
/*
void
pressure_to_height(const double *pressureLevels, double *heightLevels, int numLevels)
{
  for (int k = 0; k < numLevels; ++k) { heightLevels[k] = std::log(pressureLevels[k] / SCALESLP) * SCALEHEIGHT; }
}
*/
template <typename T>
void
vct_to_hybrid_pressure(T *restrict fullPress, T *halfPress, Varray<double> const &vct, const T *restrict ps, long nhlev, long ngp)
{
  assert(ps != nullptr);

  auto halfpresure = halfPress;
  for (long lh = 0; lh < nhlev; lh++)
  {
    auto zp = vct[lh];
    auto ze = vct[lh + nhlev + 1];
    for (long i = 0; i < ngp; ++i) halfpresure[i] = zp + ze * ps[i];
    halfpresure += ngp;
  }
  array_copy(ngp, ps, halfpresure);

  if (fullPress)
  {
    halfpresure = halfPress;
    for (long i = 0; i < ngp * nhlev; ++i) fullPress[i] = 0.5 * (halfpresure[i] + halfpresure[i + ngp]);
  }
}

template void vct_to_hybrid_pressure(float *fullPress, float *halfPress, Varray<double> const &vct, const float *ps, long nhlev,
                                     long ngp);
template void vct_to_hybrid_pressure(double *fullPress, double *halfPress, Varray<double> const &vct, const double *ps, long nhlev,
                                     long ngp);

void
extrapolate_P(double *restrict slp, const double *restrict halfPress, const double *restrict fullPress, const double *restrict geop,
              const double *restrict temp, long ngp)
{
  constexpr auto zlapse = 0.0065;
  auto zrg = 1.0 / PlanetGrav;

  for (long j = 0; j < ngp; ++j)
  {
    if (geop[j] < 0.0001 && geop[j] > -0.0001) { slp[j] = halfPress[j]; }
    else
    {
      auto alpha = PlanetRD * zlapse * zrg;
      auto tstar = (1.0 + alpha * (halfPress[j] / fullPress[j] - 1.0)) * temp[j];

      if (tstar < 255.0) { tstar = 0.5 * (255.0 + tstar); }

      auto tmsl = tstar + zlapse * zrg * geop[j];
      if (tmsl > 290.5 && tstar > 290.5)
      {
        tstar = 0.5 * (290.5 + tstar);
        tmsl = tstar;
      }

      if (tmsl - tstar < 0.000001 && tstar - tmsl < 0.000001) { alpha = 0.0; }
      else if (geop[j] > 0.0001 || geop[j] < -0.0001) { alpha = PlanetRD * (tmsl - tstar) / geop[j]; }

      auto zprt = geop[j] / (PlanetRD * tstar);
      auto zprtal = zprt * alpha;
      slp[j] = halfPress[j] * std::exp(zprt * (1.0 - zprtal * (0.5 - zprtal / 3.0)));
    }
  }
}

static inline double
extrapolate_T(double pres, double halfPress, double fullPress, double geop, double temp)
{
  auto peval = 0.0;
  constexpr auto zlapse = 0.0065;
  auto zrg = 1.0 / PlanetGrav;
  auto tstar = (1.0 + zlapse * PlanetRD * zrg * (halfPress / fullPress - 1.0)) * temp;
  auto ztsz = tstar;
  auto z1 = tstar + zlapse * zrg * geop;

  if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);

  auto ztmsl = tstar + zlapse * zrg * geop;

  if (ztmsl > 290.5 && tstar > 290.5)
  {
    tstar = 0.5 * (290.5 + tstar);
    // ztmsl = tstar;
  }

  // if (ztmsl > 290.5 && tstar <= 290.5) ztmsl = 290.5;

  if (pres <= halfPress) { peval = ((halfPress - pres) * temp + (pres - fullPress) * tstar) / (halfPress - fullPress); }
  else
  {
    ztmsl = z1;
    tstar = ztsz;
    auto zhts = geop * zrg;

    if (zhts > 2000. && z1 > 298.0)
    {
      ztmsl = 298.0;
      if (zhts < 2500.0) ztmsl = 0.002 * ((2500.0 - zhts) * z1 + (zhts - 2000.0) * ztmsl);
    }

    double zalph;
    if ((ztmsl - tstar) < 0.000001)
      zalph = 0.;
    else if (geop > 0.0001 || geop < -0.0001)
      zalph = PlanetRD * (ztmsl - tstar) / geop;
    else
      zalph = PlanetRD * zlapse * zrg;

    auto zalp = zalph * std::log(pres / halfPress);
    peval = tstar * (1.0 + zalp * (1.0 + zalp * (0.5 + 0.16666666667 * zalp)));
  }

  return peval;
}

static inline double
extrapolate_Z(double pres, double halfPress, double fullPress, double geop, double temp)
{
  constexpr auto zlapse = 0.0065;
  constexpr auto ztlim = 290.5;
  auto zrg = 1.0 / PlanetGrav;
  auto alpha = PlanetRD * zlapse * zrg;
  auto tstar = (1.0 + alpha * (halfPress / fullPress - 1.0)) * temp;

  if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);

  auto tmsl = tstar + zlapse * zrg * geop;

  if (tmsl > ztlim && tstar > ztlim)
  {
    tstar = 0.5 * (ztlim + tstar);
    tmsl = tstar;
  }

  if (tmsl > ztlim && tstar <= ztlim) tmsl = ztlim;

  if ((tmsl - tstar) < 0.000001 && (tstar - tmsl) < 0.000001)
    alpha = 0.0;
  else if (geop > 0.0001 || geop < -0.0001)
    alpha = PlanetRD * (tmsl - tstar) / geop;

  auto zalp = std::log(pres / halfPress);
  auto zalpal = zalp * alpha;

  return (geop - PlanetRD * tstar * zalp * (1.0 + zalpal * (0.5 + zalpal / 6.0))) * zrg;
}

template <typename T>
void
vertical_interp_T(const T *restrict geop, const T *restrict gt, T *pt, const T *restrict fullPress, const T *restrict halfPress,
                  const int *vertIndex, const double *restrict plev, long nplev, long ngp, long nhlev, double missval)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long lp = 0; lp < nplev; lp++)
  {
    auto pres = plev[lp];
    const int *restrict vertIndexLev = vertIndex + lp * ngp;
    auto ptl = pt + lp * ngp;

    for (long i = 0; i < ngp; ++i)
    {
      long nl = vertIndexLev[i];
      if (nl < 0) { ptl[i] = missval; }
      else
      {
        if (nl > nhlev - 2)
        {
          ptl[i] = extrapolate_T(pres, halfPress[nhlev * ngp + i], fullPress[(nhlev - 1) * ngp + i], geop[i],
                                 gt[(nhlev - 1) * ngp + i]);
        }
        else
        {
          auto nh = nl + 1;
          ptl[i] = gt[nl * ngp + i]
                   + (pres - fullPress[nl * ngp + i]) * (gt[nh * ngp + i] - gt[nl * ngp + i])
                         / (fullPress[nh * ngp + i] - fullPress[nl * ngp + i]);
        }
      }
    }
  }
}

// Explicit instantiation
template void vertical_interp_T(const float *restrict geop, const float *restrict gt, float *pt, const float *restrict fullPress,
                                const float *restrict halfPress, const int *vertIndex, const double *restrict plev, long nplev,
                                long ngp, long nhlev, double missval);
template void vertical_interp_T(const double *restrict geop, const double *restrict gt, double *pt,
                                const double *restrict fullPress, const double *restrict halfPress, const int *vertIndex,
                                const double *restrict plev, long nplev, long ngp, long nhlev, double missval);

template <typename T>
void
vertical_interp_Z(const T *restrict geop, const T *restrict gz, T *pz, const T *restrict fullPress, const T *restrict halfPress,
                  const int *vertIndex, const T *restrict gt, const double *restrict plev, long nplev, long ngp, long nhlev,
                  double missval)
{
  assert(geop != nullptr);
  assert(gz != nullptr);
  assert(pz != nullptr);
  assert(fullPress != nullptr);
  assert(halfPress != nullptr);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long lp = 0; lp < nplev; lp++)
  {
    auto pres = plev[lp];
    const int *restrict vertIndexLev = vertIndex + lp * ngp;
    auto pzl = pz + lp * ngp;

    for (long i = 0; i < ngp; ++i)
    {
      long nl = vertIndexLev[i];
      if (nl < 0) { pzl[i] = missval; }
      else
      {
        if (pres > halfPress[(nl + 1) * ngp + i]) nl++;

        if (nl > nhlev - 1)
        {
          pzl[i] = extrapolate_Z(pres, halfPress[nhlev * ngp + i], fullPress[(nhlev - 1) * ngp + i], geop[i],
                                 gt[(nhlev - 1) * ngp + i]);
        }
        else
        {
          auto nh = nl + 1;
          pzl[i] = gz[nl * ngp + i]
                   + (pres - halfPress[nl * ngp + i]) * (gz[nh * ngp + i] - gz[nl * ngp + i])
                         / (halfPress[nh * ngp + i] - halfPress[nl * ngp + i]);
        }
      }
    }
  }
}

// Explicit instantiation
template void vertical_interp_Z(const float *restrict geop, const float *restrict gz, float *pz, const float *restrict fullPress,
                                const float *restrict halfPress, const int *vertIndex, const float *restrict gt,
                                const double *restrict plev, long nplev, long ngp, long nhlev, double missval);
template void vertical_interp_Z(const double *restrict geop, const double *restrict gz, double *pz,
                                const double *restrict fullPress, const double *restrict halfPress, const int *vertIndex,
                                const double *restrict gt, const double *restrict plev, long nplev, long ngp, long nhlev,
                                double missval);

template <typename T>
static inline double
vertical_interp_X_kernel(const T *arrayIn, const T *levels3D, long nl, double level, long ngp, long nhlev, double missval)
{
  auto nh = nl + ngp;
  return (nl < 0) ? missval
                  : ((nh >= ngp * nhlev)
                         ? arrayIn[nl]
                         : arrayIn[nl] + (level - levels3D[nl]) * (arrayIn[nh] - arrayIn[nl]) / (levels3D[nh] - levels3D[nl]));
}

template <typename T>
void
vertical_interp_X(const T *restrict arrayIn3D, T *restrict arrayOut3D, const T *levels3D, const int *vertIndex3D,
                  const double *restrict levels, long numLevels, long numGP, long nhlev, double missval)
{
  if (numLevels > 3)
  {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (long levIndex = 0; levIndex < numLevels; levIndex++)
    {
      auto level = levels[levIndex];
      const int *restrict vertIndex = vertIndex3D + levIndex * numGP;
      auto *restrict arrayOut2D = arrayOut3D + levIndex * numGP;
      for (long i = 0; i < numGP; ++i)
      {
        arrayOut2D[i] = vertical_interp_X_kernel(arrayIn3D, levels3D, vertIndex[i] * numGP + i, level, numGP, nhlev, missval);
      }
    }
  }
  else
  {
    for (long levIndex = 0; levIndex < numLevels; levIndex++)
    {
      auto level = levels[levIndex];
      const int *restrict vertIndex = vertIndex3D + levIndex * numGP;
      auto *restrict arrayOut2D = arrayOut3D + levIndex * numGP;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (long i = 0; i < numGP; ++i)
      {
        arrayOut2D[i] = vertical_interp_X_kernel(arrayIn3D, levels3D, vertIndex[i] * numGP + i, level, numGP, nhlev, missval);
      }
    }
  }
}

// Explicit instantiation
template void vertical_interp_X(const float *restrict arrayIn3D, float *arrayOut3D, const float *levels3D, const int *vertIndex3D,
                                const double *levels, long numLevels, long numGP, long nhlev, double missval);
template void vertical_interp_X(const double *restrict arrayIn3D, double *arrayOut3D, const double *levels3D,
                                const int *vertIndex3D, const double *levels, long numLevels, long numGP, long nhlev,
                                double missval);

template <typename T>
void
gen_vert_index(int *vertIndex, const double *restrict plev, const T *restrict levels3D, long ngp, int nplev, int nhlev,
               bool lreverse)
{
  std::ranges::fill_n(vertIndex, ngp * nplev, 0);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (int lp = 0; lp < nplev; lp++)
  {
    T pres = plev[lp];
    auto *restrict vertIndexLev = vertIndex + lp * ngp;
    for (int lh = 0; lh < nhlev; lh++)
    {
      const auto *restrict fullpx = levels3D + lh * ngp;
      if (lreverse)
      {
        for (long i = 0; i < ngp; ++i)
        {
          if (pres < fullpx[i]) vertIndexLev[i] = lh;
        }
      }
      else
      {
        for (long i = 0; i < ngp; ++i)
        {
          if (pres > fullpx[i]) vertIndexLev[i] = lh;
        }
      }
    }
  }
}

// Explicit instantiation
template void gen_vert_index(int *vertIndex, const double *plev, const float *levels3D, long ngp, int nplev, int nhlev,
                             bool lreverse);
template void gen_vert_index(int *vertIndex, const double *plev, const double *levels3D, long ngp, int nplev, int nhlev,
                             bool lreverse);

template <typename T>
void
gen_vert_index_mv(int *vertIndex, const double *restrict plev, long ngp, int nplev, const T *restrict psProg,
                  size_t *restrict pnumMissVals, bool lreverse)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (int lp = 0; lp < nplev; lp++)
  {
    pnumMissVals[lp] = 0;
    T pres = plev[lp];
    auto *restrict vertIndexLev = vertIndex + lp * ngp;

    if (lreverse)
    {
      for (long i = 0; i < ngp; ++i)
      {
        if (pres < psProg[i])
        {
          vertIndexLev[i] = -1;
          pnumMissVals[lp]++;
        }
      }
    }
    else
    {
      for (long i = 0; i < ngp; ++i)
      {
        if (pres > psProg[i])
        {
          vertIndexLev[i] = -1;
          pnumMissVals[lp]++;
        }
      }
    }
  }
}

// Explicit instantiation
template void gen_vert_index_mv(int *vertIndex, const double *plev, long ngp, int nplev, const float *psProg, size_t *pnumMissVals,
                                bool lreverse);
template void gen_vert_index_mv(int *vertIndex, const double *plev, long ngp, int nplev, const double *psProg, size_t *pnumMissVals,
                                bool lreverse);
