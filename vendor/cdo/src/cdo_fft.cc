// This source code is copied from PINGO version 1.5

#include <cmath>
#include <algorithm>  // std::swap()

namespace cdo
{

void
fft(double *real, double *imag, int n, int sign)
{
  // n must be a power of 2
  // sign should be 1 (FT) or -1 (reverse FT)

  // Bit reversal part
  for (int i = 0, j = 0, bit; i < n; ++i)  // The bit pattern of i and j are reverse
    {
      if (i > j) std::swap(real[i], real[j]);
      if (i > j) std::swap(imag[i], imag[j]);

      for (bit = n >> 1; j & bit; bit >>= 1) j ^= bit;
      j |= bit;
    }

  // Danielson-Lanczos Part
  for (int step = 1; step < n; step <<= 1)
    {
      auto w_r = std::cos(M_PI / step);
      auto w_i = std::sin(M_PI / step) * sign;
      double ww_r = 1.0;
      double ww_i = 0.0;
      for (int i = 0; i < step; ++i)
        {
          for (int j1 = i, j2 = i + step; j2 < n; j1 += 2 * step, j2 += 2 * step)
            {
              auto temp_r = ww_r * real[j2] - ww_i * imag[j2];
              auto temp_i = ww_r * imag[j2] + ww_i * real[j2];
              real[j2] = real[j1] - temp_r;
              imag[j2] = imag[j1] - temp_i;
              real[j1] += temp_r;
              imag[j1] += temp_i;
            }
          auto temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
    }

  auto norm = 1.0 / std::sqrt(n);
  for (int i = 0; i < n; ++i) real[i] *= norm;
  for (int i = 0; i < n; ++i) imag[i] *= norm;
}
/* not used
void
ft(double *real, double *imag, int n, int sign, double *work_r, double *work_i)
{
  // sign should be 1 (FT) or -1 (reverse FT)
  double temp_r;

  for (int k = 0; k < n; ++k)
    {
      auto w_r = std::cos(2 * M_PI * k / n);
      auto w_i = std::sin(2 * M_PI * k / n) * sign;
      double ww_r = 1.0;
      double ww_i = 0.0;
      double sum_r = 0.0;
      double sum_i = 0.0;
      for (int j = 0; j < n; ++j)
        {
          sum_r += real[j] * ww_r - imag[j] * ww_i;
          sum_i += real[j] * ww_i + imag[j] * ww_r;
          temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
      work_r[k] = sum_r;
      work_i[k] = sum_i;
    }

  auto norm = 1. / std::sqrt(n);
  for (int k = 0; k < n; ++k) real[k] = work_r[k] * norm;
  for (int k = 0; k < n; ++k) imag[k] = work_i[k] * norm;
}
*/
// reentrant version of ft
void
ft_r(double *real, double *imag, int n, int sign, double *work_r, double *work_i)
{
  // sign should be 1 (FT) or -1 (reverse FT)
  double temp_r;

  for (int k = 0; k < n; ++k)
    {
      auto w_r = std::cos(2 * M_PI * k / n);
      auto w_i = std::sin(2 * M_PI * k / n) * sign;
      double ww_r = 1.0;
      double ww_i = 0.0;
      double sum_r = 0.0;
      double sum_i = 0.0;
      for (int j = 0; j < n; ++j)
        {
          sum_r += real[j] * ww_r - imag[j] * ww_i;
          sum_i += real[j] * ww_i + imag[j] * ww_r;
          temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
      work_r[k] = sum_r;
      work_i[k] = sum_i;
    }

  auto norm = 1.0 / std::sqrt(n);
  for (int k = 0; k < n; ++k) real[k] = work_r[k] * norm;
  for (int k = 0; k < n; ++k) imag[k] = work_i[k] * norm;
}

}  // namespace cdo
