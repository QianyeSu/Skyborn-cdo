#include "cdo_fftw3.h"

#include "cdo_options.h"
#include "cdo_omp.h"
#include <cdi.h>

#ifdef HAVE_LIBFFTW3
#include <mutex>
static std::mutex fftwMutex;
#endif

#ifdef HAVE_LIBFFTW3
void
fourier2grid(int gridID1, Varray<double> const &array1, Varray<double> &array2)
{
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

  struct FourierMemory
  {
    fftw_complex *in_fft;
    double *out_fft;
    fftw_plan plan;
  };

  std::vector<FourierMemory> ompmem(Threading::ompNumMaxThreads);

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
    {
      ompmem[i].in_fft = fftw_alloc_complex(nlon);
      ompmem[i].out_fft = (double *) fftw_malloc(nlon * sizeof(double));
      std::scoped_lock lock(fftwMutex);
      ompmem[i].plan = fftw_plan_dft_c2r_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
    }

  if (Options::cdoVerbose) fftw_print_plan(ompmem[0].plan);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      const auto out_fft = ompmem[ompthID].out_fft;

      for (size_t ifc = 0; ifc < nlon; ++ifc)
        {
          in_fft[ifc][0] = array1[2 * (ilat * nlon + ifc)];
          in_fft[ifc][1] = array1[2 * (ilat * nlon + ifc) + 1];
        }

      fftw_execute(ompmem[ompthID].plan);

      for (size_t ilon = 0; ilon < nlon; ++ilon) array2[ilat * nlon + ilon] = out_fft[ilon];
    }

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
    {
      fftw_free(ompmem[i].in_fft);
      fftw_free(ompmem[i].out_fft);
      std::scoped_lock lock(fftwMutex);
      fftw_destroy_plan(ompmem[i].plan);
    }
}

void
grid2fourier(int gridID1, Varray<double> const &array1, int gridID2, Varray<double> &array2)
{
  (void) gridID2;
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

  double norm = 1.0 / nlon;
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
      ompmem[i].out_fft = fftw_alloc_complex(nlon);
      std::scoped_lock lock(fftwMutex);
      ompmem[i].plan = fftw_plan_dft_r2c_1d(nlon, ompmem[i].in_fft, ompmem[i].out_fft, FFTW_ESTIMATE);
    }

  if (Options::cdoVerbose) fftw_print_plan(ompmem[0].plan);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t ilat = 0; ilat < nlat; ++ilat)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto in_fft = ompmem[ompthID].in_fft;
      auto out_fft = ompmem[ompthID].out_fft;

      for (size_t ilon = 0; ilon < nlon; ++ilon) in_fft[ilon] = array1[ilat * nlon + ilon];

      fftw_execute(ompmem[ompthID].plan);

      for (size_t ifc = 0; ifc < nlon; ++ifc)
        {
          array2[2 * (ilat * nlon + ifc)] = norm * out_fft[ifc][0];
          array2[2 * (ilat * nlon + ifc) + 1] = norm * out_fft[ifc][1];
        }
    }

  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
    {
      fftw_free(ompmem[i].in_fft);
      fftw_free(ompmem[i].out_fft);
      std::scoped_lock lock(fftwMutex);
      fftw_destroy_plan(ompmem[i].plan);
    }
}

void
filter_fftw(int nts, std::vector<int> const &fmasc, fftw_complex *fft_out, fftw_plan *p_T2S, fftw_plan *p_S2T)
{
  fftw_execute(*p_T2S);

  for (int i = 0; i < nts; ++i)
    if (!fmasc[i])
      {
        fft_out[i][0] = 0.0;
        fft_out[i][1] = 0.0;
      }

  fftw_execute(*p_S2T);

  return;
}

#else

#include "cdo_output.h"
void
fourier2grid(int gridID1, Varray<double> const &array1, Varray<double> &array2)
{
  (void) gridID1;
  (void) array1;
  (void) array2;
  cdo_abort("FFTW support not compiled in!");
}

void
grid2fourier(int gridID1, Varray<double> const &array1, int gridID2, Varray<double> &array2)
{
  (void) gridID1;
  (void) gridID2;
  (void) array1;
  (void) array2;
  cdo_abort("FFTW support not compiled in!");
}

#endif
