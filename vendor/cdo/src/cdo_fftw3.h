#ifndef CDO_FFTW3_H
#define CDO_FFTW3_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "varray.h"

void fourier2grid(int gridID1, Varray<double> const &array1, Varray<double> &array2);

void grid2fourier(int gridID1, Varray<double> const &array1, int gridID2, Varray<double> &array2);

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
void filter_fftw(int nts, std::vector<int> const &fmasc, fftw_complex *fft_out, fftw_plan *p_T2S, fftw_plan *p_S2T);
#endif

#endif
