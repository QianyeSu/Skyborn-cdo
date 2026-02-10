#ifndef CDO_FFT_H
#define CDO_FFT_H

namespace cdo
{
void fft(double *real, double *imag, int n, int sign);
void ft_r(double *real, double *imag, int n, int sign, double *work_r, double *work_i);
}  // namespace cdo

#endif
