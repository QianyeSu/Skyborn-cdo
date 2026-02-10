/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>

#include "constants.h"

#define SQUARE_RADIUS (-PlanetRadiusDefault * PlanetRadiusDefault)

void
dv2ps(const double *restrict div, double *restrict pot, long nlev, long ntr)
{
  long l, m, n;
  double fact;

  for (l = 0; l < nlev; ++l)
  {
    /* m == 0 */
    *pot++ = 0.0;
    *pot++ = 0.0;
    div += 2;

    for (n = 1; n <= ntr; ++n)
    {
      fact = SQUARE_RADIUS / (n * n + n);
      *pot++ = *div++ * fact;
      *pot++ = *div++ * fact;
    }

    /* m >= 0 */
    for (m = 1; m <= ntr; ++m)
      for (n = m; n <= ntr; ++n)
      {
        fact = SQUARE_RADIUS / (n * n + n);
        *pot++ = *div++ * fact;
        *pot++ = *div++ * fact;
      }
  }
}

static void
dv2uv_kernel(const double *d, const double *o, double *restrict u, double *restrict v, const double *f, const double *g, long nt)
{
  // d(nsp,nlev), o(nsp,nlev)     ! divergence, vorticity
  // u(nsp,nlev), v(nsp,nlev)     ! zonal wind, meridional wind
  // f(nsp/2)   , g(nsp/2)        ! factor tables

  long i = 0;
  for (long m = 0; m < nt - 1; ++m)
  {
    // n = m
    if (m == 0)
    {
      *u++ = -g[i + 1] * o[2 * (i + 1)];
      *u++ = -g[i + 1] * o[2 * (i + 1) + 1];
      *v++ = g[i + 1] * d[2 * (i + 1)];
      *v++ = g[i + 1] * d[2 * (i + 1) + 1];
    }
    else
    {
      *u++ = -f[i] * d[2 * i + 1] - g[i + 1] * o[2 * (i + 1)];
      *u++ = f[i] * d[2 * i] - g[i + 1] * o[2 * (i + 1) + 1];
      *v++ = -f[i] * o[2 * i + 1] + g[i + 1] * d[2 * (i + 1)];
      *v++ = f[i] * o[2 * i] + g[i + 1] * d[2 * (i + 1) + 1];
    }
    ++i;

    // m < n < nt-1
    for (long n = m + 1; n < nt - 1; ++n)
    {
      *u++ = g[i] * o[2 * (i - 1)] - f[i] * d[2 * i + 1] - g[i + 1] * o[2 * (i + 1)];
      *u++ = g[i] * o[2 * (i - 1) + 1] + f[i] * d[2 * i] - g[i + 1] * o[2 * (i + 1) + 1];
      *v++ = -g[i] * d[2 * (i - 1)] - f[i] * o[2 * i + 1] + g[i + 1] * d[2 * (i + 1)];
      *v++ = -g[i] * d[2 * (i - 1) + 1] + f[i] * o[2 * i] + g[i + 1] * d[2 * (i + 1) + 1];
      ++i;
    }

    // n = nt-1
    *u++ = g[i] * o[2 * (i - 1)] - f[i] * d[2 * i + 1];
    *u++ = g[i] * o[2 * (i - 1) + 1] + f[i] * d[2 * i];
    *v++ = -g[i] * d[2 * (i - 1)] - f[i] * o[2 * i + 1];
    *v++ = -g[i] * d[2 * (i - 1) + 1] + f[i] * o[2 * i];
    ++i;

    // n = nt
    *u++ = g[i] * o[2 * (i - 1)];
    *u++ = g[i] * o[2 * (i - 1) + 1];
    *v++ = -g[i] * d[2 * (i - 1)];
    *v++ = -g[i] * d[2 * (i - 1) + 1];
    ++i;
  }

  // m = nt-1  and  n = nt-1
  *u++ = -f[i] * d[2 * i + 1];
  *u++ = f[i] * d[2 * i];
  *v++ = -f[i] * o[2 * i + 1];
  *v++ = f[i] * o[2 * i];
  ++i;

  // m = nt-1  and  n = nt
  *u++ = g[i] * o[2 * (i - 1)];
  *u++ = g[i] * o[2 * (i - 1) + 1];
  *v++ = -g[i] * d[2 * (i - 1)];
  *v++ = -g[i] * d[2 * (i - 1) + 1];

  // m = nt  and  n = nt
  *u++ = 0.0;
  *u++ = 0.0;
  *v++ = 0.0;
  *v++ = 0.0;
}

void
dv2uv(const double *d, const double *o, double *u, double *v, const double *f, const double *g, long nt, long nsp, long nlev)
{
  // d(nsp,nlev), o(nsp,nlev)     ! divergence, vorticity
  // u(nsp,nlev), v(nsp,nlev)     ! zonal wind, meridional wind
  // f(nsp/2)   , g(nsp/2)        ! factor tables

#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
  for (long l = 0; l < nlev; ++l) { dv2uv_kernel(d + l * nsp, o + l * nsp, u + l * nsp, v + l * nsp, f, g, nt); }
}

/*
void scaluv(double *fu, const double *rclat, long nlat, long lot)
{
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (long l = 0; l < lot; l++ )
    {
      double *ful = fu + l*nlat;
      for (long lat = 0; lat < nlat; lat++ )
        {
          ful[lat] = rclat[lat];
        }
    }
}
*/

void
scaluv(double *fu, const double *rclat, long nlat, long lot)
{
  for (long l = 0; l < lot; ++l)
    for (long lat = 0; lat < nlat; lat++)
    {
      *fu *= rclat[lat];
      fu++;
    }
}

void
uv2dv(const double *fu, const double *fv, double *sd, double *sv, const double *pol2, const double *pol3, long klev, long nlat,
      long nt)
{
  long nsp2 = (nt + 1) * (nt + 2);
  long nfc = (nt + 1) * 2;

#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
  for (long lev = 0; lev < klev; lev++)
  {
    auto po2 = pol2;
    auto po3 = pol3;
    auto ful = fu + lev * nfc * nlat;
    auto fvl = fv + lev * nfc * nlat;
    auto sdl = sd + lev * nsp2;
    auto svl = sv + lev * nsp2;
    for (long jmm = 0; jmm <= nt; jmm++)
    {
      for (long jfc = jmm; jfc <= nt; jfc++)
      {
        auto ufr = ful;
        auto ufi = ful + nlat;
        auto vfr = fvl;
        auto vfi = fvl + nlat;
        double dir = 0.0;
        double dii = 0.0;
        double vor = 0.0;
        double voi = 0.0;
        for (long lat = 0; lat < nlat; lat++)
        {
          dir += vfr[lat] * po2[lat] - ufi[lat] * po3[lat];
          dii += vfi[lat] * po2[lat] + ufr[lat] * po3[lat];
          vor -= ufr[lat] * po2[lat] + vfi[lat] * po3[lat];
          voi -= ufi[lat] * po2[lat] - vfr[lat] * po3[lat];
        }
        *sdl++ = dir;
        *sdl++ = dii;
        *svl++ = vor;
        *svl++ = voi;
        po2 += nlat;
        po3 += nlat;
      }
      ful += 2 * nlat;
      fvl += 2 * nlat;
    }
  }
}

void
geninx(long ntr, double *f, double *g)
{
  for (long m = 0; m <= ntr; ++m)
  {
    long m2 = m * m;
    for (long n = m; n <= ntr; ++n)
    {
      long n2 = n * n;
      if (n)
      {
        *g++ = -PlanetRadiusDefault / n * std::sqrt((double) (n2 - m2) / (double) (4 * n2 - 1));
        *f++ = -PlanetRadiusDefault * m / (double) (n2 + n);
      }
      else
      {
        *g++ = 0.0;
        *f++ = 0.0;
      }
    }
  }
}
