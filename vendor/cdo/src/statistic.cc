// This source code is copied from PINGO version 1.5

#include <cfloat>
#include <cmath>
#include <cstdio>

#include "compare.h"
#include "statistic.h"

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390 /* 2/std::sqrt(pi) */
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880 /* std::sqrt(2) */
#endif

namespace cdo
{

// same result as std::lgamma()
static double
lngamma(double x)
{
  constexpr double cof[6] = { 76.18009172947146,  -86.50532032941677,    24.01409824083091,
                              -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };

  double b, a = b = x;
  double temp = a + 5.5;
  temp -= (a + 0.5) * std::log(temp);
  double ser = 1.000000000190015;
  for (int j = 0; j <= 5; ++j) ser += cof[j] / ++b;
  return -temp + std::log(2.5066282746310005 * ser / a);
}

static double
gamma_help_1(double a, double x)
{
  constexpr double eps = 1.e-20;

  double gln = lngamma(a);
  double ap = a;
  double sum, del = sum = 1.0 / a;

  for (int i = 1; i <= 100; ++i)  // 100 iterations
  {
    ap++;
    del *= x / ap;
    sum += del;
    if (std::fabs(del) < std::fabs(sum) * eps) return sum * std::exp(-x + a * std::log(x) - (gln));
  }

  fprintf(stderr, "%s: internal error, too many iterations!\n", __func__);
  exit(1);

  return 0;
}

static double
gamma_help_2(double a, double x)
{
  constexpr double eps = 1.e-20;
  constexpr double very_small = 1000.0 * DBL_MIN;

  double gln = lngamma(a);
  double b = x + 1.0 - a;
  double c = 1.0 / very_small;
  double d = 1.0 / b;
  double h = d;

  for (int i = 1; i <= 100; ++i)  // 100 iterations
  {
    double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < very_small) d = very_small;
    c = b + an / c;
    if (std::fabs(c) < very_small) c = very_small;
    d = 1 / d;
    double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < eps) return std::exp(-x + a * std::log(x) - gln) * h;
  }

  fprintf(stderr, "%s: internal error, too many iterations!\n", __func__);
  exit(1);

  return -1;
}

static double
incomplete_gamma(double a, double x)
{
  if (x < 0.0 || a <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return (x < (a + 1.0)) ? gamma_help_1(a, x) : 1.0 - gamma_help_2(a, x);
}

double
beta(double a, double b)
{
  if (a <= 0.0 || b <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return std::exp(lngamma(a) + lngamma(b) - lngamma(a + b));
}

static double
beta_help(double a, double b, double x)
{
  constexpr double very_small = 1000.0 * DBL_MIN;
  constexpr double eps = 3.e-07;

  auto qab = a + b;
  auto qap = a + 1;
  auto qam = a - 1;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < very_small) d = very_small;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= 100; ++m)
  {
    int m2 = 2 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < very_small) d = very_small;
    c = 1.0 + aa / c;
    if (std::fabs(c) < very_small) c = very_small;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < very_small) d = very_small;
    c = 1.0 + aa / c;
    if (std::fabs(c) < very_small) c = very_small;
    d = 1.0 / d;
    double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < eps) return h;
  }

  fprintf(stderr, "%s: ERROR! Too many iterations in routine!\n", __func__);
  exit(1);

  return -1;
}

static double
incomplete_beta(double a, double b, double x)
{
  if (a <= 0.0 || b <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (x < 0.0 || x > 1.0)
  {
    fprintf(stderr, "%s: Value out of range (0-1)!\n", __func__);
    exit(4);
  }

  double c = (fp_is_equal(x, 0.) || fp_is_equal(x, 1.))
                 ? 0.0
                 : std::exp(lngamma(a + b) - lngamma(a) - lngamma(b) + a * std::log(x) + b * std::log(1 - x));

  if (x < (a + 1) / (a + b + 2.0))
    return c * beta_help(a, b, x) / a;
  else
    return 1.0 - c * beta_help(b, a, 1.0 - x) / b;
}

double
normal_density(double x)
{
  return M_2_SQRTPI / 2.0 / M_SQRT2 * std::exp(-x * x / 2.0);
}

double
normal(double x)
{
  return (x > 0.0)   ? 0.5 * (1.0 + incomplete_gamma(0.5, x * x / 2.0))
         : (x < 0.0) ? 0.5 * (1.0 - incomplete_gamma(0.5, x * x / 2.0))
                     : 0.5;
}

double
normal_inv(double p)
{
  constexpr double eps = 1.e-10;
  static double last_p = 0.5, last_x = 0.0;

  if (p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(p, last_p)) return last_x;

  if (p < 0.5) return -normal_inv(1 - p);
  if (p > 0.5)
  {
    double x = 0.0;
    while (true)
    {
      double xx = x - (normal(x) - p) / normal_density(x);
      if (std::fabs(xx - x) < x * eps) break;
      x = xx;
    }
    last_p = p;
    last_x = x;
    return x;
  }

  return 0;
}

double
student_t_density(double n, double x)
{
  if (n <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return std::exp(lngamma((n + 1) / 2) - lngamma(n / 2)) / std::sqrt(n / 2) * std::pow((1 + x * x / n), -(n + 1) / 2);
}

double
student_t(double n, double x)
{
  if (n <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (x > 0)
    return 1.0 - 0.5 * incomplete_beta(n / 2.0, 0.5, n / (n + x * x));
  else if (x < 0)
    return 0.5 * incomplete_beta(n / 2.0, 0.5, n / (n + x * x));
  else
    return 0.5;
}

double
student_t_inv(double n, double p)
{
  constexpr double eps = 1.e-10;
  static double last_n = 1.0, last_p = 0.5, last_x = 0.0;

  if (n <= 0.0 || p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(n, last_n) && fp_is_equal(p, last_p)) return last_x;

  if (fp_is_equal(p, 0.5))
    return 0.0;
  else if (p < 0.5)
    return -student_t_inv(n, 1.0 - p);
  else
  {
    double x = 0.0;
    while (true)
    {
      double xx = x - (student_t(n, x) - p) / student_t_density(n, x);
      if (std::fabs(xx - x) < x * eps) break;
      x = xx;
    }
    last_n = n;
    last_p = p;
    last_x = x;
    return x;
  }
}

double
chi_square_density(double n, double x)
{
  if (n <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return (x <= 0) ? 0 : std::pow(2, -n / 2) * std::pow(x, n / 2 - 1) * std::exp(-x / 2 - lngamma(n / 2));
}

double
chi_square(double n, double x)
{
  if (n <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return (x <= 0.0) ? 0.0 : incomplete_gamma(n / 2.0, x / 2.0);
}

double
chi_square_inv(double n, double p)
{
  constexpr double eps = 1.e-10;
  static double last_n = -1.0, last_p = -1.0, last_x = -1.0;
  static double last_last_n = -1.0, last_last_p = -1.0, last_last_x = -1.0;

  if (n <= 0.0 || p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(n, last_n) && fp_is_equal(p, last_p)) return last_x;

  if (fp_is_equal(n, last_last_n) && fp_is_equal(p, last_last_p)) return last_last_x;

  double x = n;
  while (true)
  {
    double xx = x - (chi_square(n, x) - p) / chi_square_density(n, x);
    if (std::fabs(xx - x) < x * eps) break;
    if (xx < 0)
      x /= 2.0;
    else
      x = xx;
  }

  last_last_n = last_n;
  last_last_p = last_p;
  last_last_x = last_x;
  last_n = n;
  last_p = p;
  last_x = x;

  return x;
}

void
chi_square_constants(double n, double p, double *c1, double *c2)
{
  constexpr double eps = 1.e-10;
  static double last_n, last_p, last_c1, last_c2;

  if (n <= 0.0 || p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(n, last_n) && fp_is_equal(p, last_p))
  {
    *c1 = last_c1;
    *c2 = last_c2;
    return;
  }

  *c1 = n;
  *c2 = n;

  while (true)
  {
    auto a11 = -chi_square_density(n, *c1);
    auto a12 = chi_square_density(n, *c2);
    auto a21 = -chi_square_density(n + 2, *c1);
    auto a22 = chi_square_density(n + 2, *c2);
    auto b1 = p + chi_square(n, *c1) - chi_square(n, *c2);
    auto b2 = p + chi_square(n + 2, *c1) - chi_square(n + 2, *c2);
    // Solve ((a11,a12),(a21,a22))*(delta_c1,delta_c2)==(b1,b2)
    auto det = a11 * a22 - a12 * a21;
    auto delta_c1 = (b1 * a22 - b2 * a12) / det;
    auto delta_c2 = (b2 * a11 - b1 * a21) / det;
    if (std::fabs(delta_c1) < *c1 * eps && std::fabs(delta_c2) < *c2 * eps) break;
    if (*c1 + delta_c1 >= n)
      *c1 = (n + *c1) / 2.0;
    else if (*c1 + delta_c1 <= 0)
      *c1 /= 2.0;
    else
      *c1 += delta_c1;
    if (*c2 + delta_c2 <= n)
      *c2 = (n + *c2) / 2.0;
    else
      *c2 += delta_c2;
  }

  last_n = n;
  last_p = p;
  last_c1 = *c1;
  last_c2 = *c2;
}

double
beta_distr_density(double a, double b, double x)
{
  if (a <= 0.0 || b <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return (x <= 0) ? 0 : (x >= 1) ? 1 : std::pow(x, a - 1) * std::pow(1 - x, b - 1) / beta(a, b);
}

double
beta_distr(double a, double b, double x)
{
  return incomplete_beta(a, b, x);
}

double
beta_distr_inv(double a, double b, double p)
{
  constexpr double eps = 1.e-10;
  static double last_a = -1.0, last_b, last_p = -1.0, last_x = -1.0;
  static double last_last_a = -1.0, last_last_b = -1.0, last_last_p = -1.0, last_last_x = -1.0;

  if (a <= 0.0 || b <= 0.0 || p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(a, last_a) && fp_is_equal(b, last_b) && fp_is_equal(p, last_p)) return last_x;

  if (fp_is_equal(a, last_last_a) && fp_is_equal(b, last_last_b) && fp_is_equal(p, last_last_p)) return last_last_x;

  double x = a / (a + b);
  while (true)
  {
    double xx = x - (beta_distr(a, b, x) - p) / beta_distr_density(a, b, x);
    if (std::fabs(xx - x) < x * eps) break;
    if (xx <= 0.0)
      x /= 2;
    else if (xx >= 1.0)
      x = (1.0 + x) / 2.0;
    else
      x = xx;
  }
#if 0
  for (x_l = 0, x_r = 1; std::fabs(x_l - x_r) > eps;
       x = (x_l+x_r) / 2.0, (beta_distr(a, b, x) < p) ? (x_l=x):(x_r=x));
#endif

  last_last_a = last_a;
  last_last_b = last_b;
  last_last_p = last_p;
  last_last_x = last_x;
  last_a = a;
  last_b = b;
  last_p = p;
  last_x = x;

  return x;
}

void
beta_distr_constants(double a, double b, double p, double *c1, double *c2)
{
  constexpr double eps = 1.e-10;
  static double last_a, last_b, last_p, last_c1, last_c2;

  if (a <= 0.0 || b <= 0.0 || p <= 0.0 || p >= 1.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  if (fp_is_equal(a, last_a) && fp_is_equal(b, last_b) && fp_is_equal(p, last_p))
  {
    *c1 = last_c1;
    *c2 = last_c2;
    return;
  }

#if 0
  *c1 = a / (a + b);
  *c2 = a / (a + b);
#endif
  *c1 = beta_distr_inv(a, b, p / 2.0);
  *c2 = beta_distr_inv(a, b, 1.0 - p / 2.0);

  while (true)
  {
    auto a11 = -beta_distr_density(a, b, *c1);
    auto a12 = beta_distr_density(a, b, *c2);
    auto a21 = -beta_distr_density(a + 1, b, *c1);
    auto a22 = beta_distr_density(a + 1, b, *c2);
    auto b1 = p + beta_distr(a, b, *c1) - beta_distr(a, b, *c2);
    auto b2 = p + beta_distr(a + 1, b, *c1) - beta_distr(a + 1, b, *c2);
    // Solve ((a11,a12),(a21,a22))*(delta_c1,delta_c2)==(b1,b2)
    auto det = a11 * a22 - a12 * a21;
    auto delta_c1 = (b1 * a22 - b2 * a12) / det;
    auto delta_c2 = (b2 * a11 - b1 * a21) / det;
    if (std::fabs(delta_c1) < *c1 * eps && std::fabs(delta_c2) < *c2 * eps) break;
    if (*c1 + delta_c1 >= a / (a + b))
      *c1 = (a / (a + b) + *c1) / 2.0;
    else if (*c1 + delta_c1 <= 0)
      *c1 /= 2.0;
    else
      *c1 += delta_c1;
    if (*c2 + delta_c2 >= 1.0)
      *c2 = (1.0 + *c2) / 2.0;
    else if (*c2 + delta_c2 <= a / (a + b))
      *c2 = (a / (a + b) + *c2) / 2.0;
    else
      *c2 += delta_c2;
  }

  last_a = a;
  last_b = b;
  last_p = p;
  last_c1 = *c1;
  last_c2 = *c2;
}

double
fisher(double m, double n, double x)
{
  if (m <= 0.0 || n <= 0.0)
  {
    fprintf(stderr, "%s: IMPLEMENTATION ERROR! (Invalid argument)\n", __func__);
    exit(4);
  }

  return incomplete_beta(m / 2.0, n / 2.0, n / (n + m * x));
}

}  // namespace cdo
