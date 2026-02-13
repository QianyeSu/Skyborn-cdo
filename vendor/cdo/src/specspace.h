/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef SPECSPACE_H
#define SPECSPACE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_options.h"
#include "cdo_output.h"
#include "transform.h"
#include "varray.h"

enum struct PolFlag
{
  UNDEF,
  SP2FC,
  FC2SP,
  UV2DV,
};

// clang-format off
class  // FC_Transformation
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
FC_Transformation
// clang-format on
{
public:
  bool use_fftw = false;
  long nlon = 0;
  long nlat = 0;
  long ntr = 0;
  long nlev = 0;
  long ifax[10] = {};
  Varray<double> vtrig;

  FC_Transformation()
  {
    if (Options::Use_FFTW)
      {
#ifdef HAVE_LIBFFTW3
        if (Options::cdoVerbose) cdo_print("Using fftw3 lib");
        use_fftw = true;
#else
        if (Options::cdoVerbose) cdo_print("LIBFFTW3 support not compiled in!");
#endif
      }
  }

  void
  init(long _nlon, long _nlat, long _ntr, long _nlev = 0)
  {
    if (_nlon <= 0 || _nlat <= 0 || _ntr <= 0)
      {
        std::fprintf(stderr, "SP_Transformation.init(): parameter not initialized\n");
        return;
      }

    nlon = _nlon;
    nlat = _nlat;
    ntr = _ntr;
    nlev = _nlev;

    // if (nlev > 1) use_fftw = false;

    if (use_fftw == false)
      {
        vtrig.resize(nlon);
        auto status = fft_set(vtrig.data(), ifax, nlon);
        if (status < 0)
          {
            cdo_print("Retry it with the fftw3 library!");
#ifdef HAVE_LIBFFTW3
            cdo_warning("Using fftw3 lib is disabled!");
#else
            cdo_warning("LIBFFTW3 support not compiled in!");
#endif
            cdo_abort("FFT error!");
          }
      }
  }
};

// clang-format off
class  // SP_Transformation
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
SP_Transformation
// clang-format on
{
public:
  FC_Transformation fcTrans;
  long nlat = 0;
  long ntr = 0;
  Varray<double> poli;
  Varray<double> pold;
  Varray<double> pol2;     // only for uv2dv
  Varray<double> pol3;     // only for uv2dv
  Varray<double> coslat;   // only for scaluv with uv2dv
  Varray<double> rcoslat;  // only for scaluv with dv2uv

  SP_Transformation() {}

  void
  init(long _nlon, long _nlat, long _ntr, PolFlag polFlag, long _nlev = 0)
  {
    if (_nlon <= 0 || _nlat <= 0 || _ntr <= 0)
      {
        std::fprintf(stderr, "SP_Transformation.init(): parameter not initialized\n");
        return;
      }

    fcTrans.init(_nlon, _nlat, _ntr, _nlev);

    nlat = _nlat;
    ntr = _ntr;

    long nsp = (ntr + 1) * (ntr + 2);
    long poldim = (nsp / 2) * nlat;

    auto numArrays = (polFlag == PolFlag::UV2DV) ? 2 : 1;
    cdo_print("Allocate %d array%s to hold the Legendre polynoms: size=%ld Bytes", numArrays, numArrays > 1 ? "s" : "", poldim * 8);

    if (polFlag == PolFlag::SP2FC) varrayResize(poli, poldim);
    if (polFlag == PolFlag::FC2SP) varrayResize(pold, poldim);

    if (polFlag == PolFlag::UV2DV) varrayResize(pol2, poldim);
    if (polFlag == PolFlag::UV2DV) varrayResize(pol3, poldim);

    coslat.resize(nlat);
    rcoslat.resize(nlat);

    after_legini_full(ntr, nlat, poli.data(), pold.data(), nullptr, pol2.data(), pol3.data(), coslat.data());

    for (long jgl = 0; jgl < nlat; ++jgl) rcoslat[jgl] = 1.0 / coslat[jgl];
  }
};

// clang-format off
class  // DV_Transformation
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
DV_Transformation
// clang-format on
{
public:
  long ntr = 0;
  long fdim = 0;
  Varray<double> f1;
  Varray<double> f2;

  DV_Transformation() {}

  void
  init(long _ntr)
  {
    ntr = _ntr;

    long dimsp = (ntr + 1) * (ntr + 2);
    fdim = dimsp / 2;

    f1.resize(fdim);
    f2.resize(fdim);

    geninx(ntr, f1.data(), f2.data());
  }
};

void dv2ps(const double *div, double *pot, long nlev, long ntr);

void trans_uv2dv(const SP_Transformation &spTrans, long nlev, int gridID1, Varray<double> const &gu, Varray<double> const &gv,
                 int gridID2, Varray<double> &sd, Varray<double> &svo);

void trans_dv2uv(const SP_Transformation &spTrans, const DV_Transformation &dvTrans, long nlev, int gridID1,
                 Varray<double> const &sd, Varray<double> const &svo, int gridID2, Varray<double> &gu, Varray<double> &gv);

void grid2spec(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);
void spec2grid(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);
void four2spec(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);
void spec2four(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);
void four2grid(const FC_Transformation &fcTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);
void grid2four(const FC_Transformation &fcTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout,
               Varray<double> &arrayOut);

void spec2spec(int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut);
void speccut(int gridIDin, Varray<double> const &arrayIn, Varray<double> &arrayOut, Varray<int> const &waves);

void spcut(const double *arrayIn, double *arrayOut, long ntr, const int *waves);

#endif
