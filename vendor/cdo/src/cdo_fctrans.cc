#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#endif

#include <vector>

#ifndef HAVE_LIBFFTW3
#include "cdo_output.h"
#endif
#include "cdo_options.h"
#include "cdo_fctrans.h"
#include "cdo_omp.h"

#ifdef HAVE_LIBFFTW3
#include <mutex>
static std::mutex fftwMutex;
#endif

void
fc2gp(const double *fc, double *gp, long nlat, long nlon, long nlev, long nfc)
{
#ifdef HAVE_LIBFFTW3
  struct FourierMemory
  {
    fftw_complex *in_fft;
    double *out_fft;
    fftw_plan plan;
  };

  std::vector<FourierMemory> ompmem(Threading::ompNumMaxThreads);

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    ompmem[i].in_fft = fftw_alloc_complex(nlon / 2 + 1);
    ompmem[i].out_fft = (double *) fftw_malloc(nlon * sizeof(double));
    std::scoped_lock lock(fftwMutex);
    ompmem[i].plan = fftw_plan_dft_c2r_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
  }

  static int nprint = 0;
  if (Options::cdoVerbose && nprint++ < 4) fftw_print_plan(ompmem[0].plan);

  for (long ilev = 0; ilev < nlev; ++ilev)
  {
    auto gpx = gp + ilev * nlon * nlat;
    auto fcx = fc + ilev * nfc * nlat;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (long ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      auto out_fft = ompmem[ompthID].out_fft;

      for (long ifc = 0; ifc < nfc / 2; ++ifc)
      {
        in_fft[ifc][0] = fcx[2 * ifc * nlat + ilat];
        in_fft[ifc][1] = fcx[(2 * ifc + 1) * nlat + ilat];
      }
      for (long ifc = nfc / 2; ifc < (nlon / 2 + 1); ++ifc)
      {
        in_fft[ifc][0] = 0.0;
        in_fft[ifc][1] = 0.0;
      }

      fftw_execute(ompmem[ompthID].plan);

      for (long ilon = 0; ilon < nlon; ++ilon) gpx[ilat * nlon + ilon] = out_fft[ilon];
    }
  }

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    fftw_free(ompmem[i].in_fft);
    fftw_free(ompmem[i].out_fft);
    std::scoped_lock lock(fftwMutex);
    fftw_destroy_plan(ompmem[i].plan);
  }
#else
  cdo_abort("FFTW support not compiled in!");
#endif
}

void
gp2fc(const double *gp, double *fc, long nlat, long nlon, long nlev, long nfc)
{
#ifdef HAVE_LIBFFTW3
  double norm = 1. / nlon;
  struct FourierMemory
  {
    double *in_fft;
    fftw_complex *out_fft;
    fftw_plan plan;
  };

  std::vector<FourierMemory> ompmem(Threading::ompNumMaxThreads);

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    ompmem[i].in_fft = (double *) fftw_malloc(nlon * sizeof(double));
    ompmem[i].out_fft = fftw_alloc_complex(nlon / 2 + 1);
    std::scoped_lock lock(fftwMutex);
    ompmem[i].plan = fftw_plan_dft_r2c_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
  }

  static int nprint = 0;
  if (Options::cdoVerbose && nprint++ < 4) fftw_print_plan(ompmem[0].plan);

  for (long ilev = 0; ilev < nlev; ++ilev)
  {
    auto gpx = gp + ilev * nlon * nlat;
    auto fcx = fc + ilev * nfc * nlat;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (long ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      auto out_fft = ompmem[ompthID].out_fft;

      for (long ilon = 0; ilon < nlon; ++ilon) in_fft[ilon] = gpx[ilat * nlon + ilon];

      fftw_execute(ompmem[ompthID].plan);

      for (long ifc = 0; ifc < nfc / 2; ++ifc)
      {
        fcx[2 * ifc * nlat + ilat] = norm * out_fft[ifc][0];
        fcx[(2 * ifc + 1) * nlat + ilat] = norm * out_fft[ifc][1];
      }
    }
  }

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    fftw_free(ompmem[i].in_fft);
    fftw_free(ompmem[i].out_fft);
    std::scoped_lock lock(fftwMutex);
    fftw_destroy_plan(ompmem[i].plan);
  }
#else
  cdo_abort("FFTW support not compiled in!");
#endif
}
