/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cmath>
#include <cstring>

#include <algorithm>
#include <unordered_map>

#include "percentiles.h"
#include "util_string.h"
#include "cdo_output.h"
#include "cdo_options.h"

enum struct PercentileMethod
{
  NRANK = 1,
  NIST,
  NUMPY,
  NR8
};

enum struct NumpyMethod
{
  // NumPy
  linear = 1,
  lower,
  higher,
  nearest,
  midpoint,
  // discontinuous sample
  inverted_cdf,
  averaged_inverted_cdf,
  closest_observation,
  // continuous sample
  interpolated_inverted_cdf,
  hazen,
  weibull,
  median_unbiased,
  normal_unbiased,
};

static PercentileMethod percentileMethod = PercentileMethod::NRANK;
static NumpyMethod numpyMethod = NumpyMethod::linear;

template <typename T>
static double
get_nth_element(T *array, size_t n, size_t idx)
{
  std::nth_element(array, array + idx, array + n);
  return array[idx];
}

template <typename T>
static double
percentile_nrank(T *array, size_t n, double quantile)
{
  auto irank = (size_t) std::ceil(n * quantile);
  irank = std::clamp(irank, static_cast<size_t>(1), n);
  return get_nth_element(array, n, irank - 1);
}

template <typename T>
static double
percentile_nist(T *array, size_t n, double quantile)
{
  double rank = (n + 1) * quantile;
  size_t k = (size_t) rank;

  double percentil = 0.0;
  if (k == 0) { percentil = get_nth_element(array, n, 0); }
  else if (k >= n) { percentil = get_nth_element(array, n, n - 1); }
  else
  {
    auto vk = get_nth_element(array, n, k - 1);
    auto vk2 = get_nth_element(array, n, k);
    double d = rank - k;
    percentil = vk + d * (vk2 - vk);
  }

  return percentil;
}

template <typename T>
static double
percentile_numpy(T *array, size_t n, double quantile)
{
  // R code: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/quantile.R
  // Python: https://github.com/numpy/numpy/blob/main/numpy/lib/function_base.py
  double rank = 1.0 + (n - 1) * quantile;
  size_t k = (size_t) rank;

  double percentil = 0.0;
  if (k == 1) { percentil = get_nth_element(array, n, 0); }
  else if (k >= n) { percentil = get_nth_element(array, n, n - 1); }
  else
  {
    if (numpyMethod == NumpyMethod::linear)
    {
      size_t lo = std::floor(rank);
      size_t hi = std::ceil(rank);
      double h = rank - lo;  // > 0	by construction
      percentil = (1.0 - h) * get_nth_element(array, n, lo - 1) + h * get_nth_element(array, n, hi - 1);
    }
    else if (numpyMethod == NumpyMethod::lower)
    {
      size_t lo = std::floor(rank);
      percentil = get_nth_element(array, n, std::clamp(lo, static_cast<size_t>(1), n) - 1);
    }
    else if (numpyMethod == NumpyMethod::higher)
    {
      size_t hi = std::ceil(rank);
      percentil = get_nth_element(array, n, std::clamp(hi, static_cast<size_t>(1), n) - 1);
    }
    else if (numpyMethod == NumpyMethod::nearest)  // numpy is using around(), with rounds to the nearest even value
    {
      size_t j = std::lround(rank);
      percentil = get_nth_element(array, n, std::clamp(j, static_cast<size_t>(1), n) - 1);
    }
    else if (numpyMethod == NumpyMethod::midpoint)
    {
      size_t lo = std::floor(rank);
      size_t hi = std::ceil(rank);
      constexpr double h = 0.5;
      percentil = h * get_nth_element(array, n, lo - 1) + h * get_nth_element(array, n, hi - 1);
    }
    else
    {
      if (numpyMethod >= NumpyMethod::interpolated_inverted_cdf)
      {
        // clang-format off
        std::unordered_map<NumpyMethod, std::pair<double, double>> abMap = {
          { NumpyMethod::interpolated_inverted_cdf, {0.0, 1.0} },
          { NumpyMethod::hazen,                     {0.5, 0.5} },
          { NumpyMethod::weibull,                   {0.0, 0.0} },
          { NumpyMethod::median_unbiased,           {1.0 / 3.0, 1.0 / 3.0} },
          { NumpyMethod::normal_unbiased,           {3.0 / 8.0, 3.0 / 8.0} },
        };
        // clang-format on

        auto [a, b] = abMap[numpyMethod];
        double nppn = a + quantile * (n + 1.0 - a - b);  // n*probs + m
        constexpr double fuzz = 4.0 * std::numeric_limits<double>::epsilon();
        size_t j = std::floor(nppn + fuzz);  // m = a + probs*(1 - a - b)
        double h = nppn - j;
        if (std::fabs(h) < fuzz) h = 0.0;
        if (h > 0.0 && h < 1.0)
          percentil = (1.0 - h) * get_nth_element(array, n, j - 1) + h * get_nth_element(array, n, j);
        else
          percentil = get_nth_element(array, n, (h >= 1.0) ? j : j - 1);
      }
      else
      {
        // discontinuous sample
        double nppm = (numpyMethod == NumpyMethod::closest_observation) ? n * quantile - 0.5 : n * quantile;
        size_t j = std::floor(nppm);
        double h = (numpyMethod == NumpyMethod::inverted_cdf)            ? (nppm > j)
                   : (numpyMethod == NumpyMethod::averaged_inverted_cdf) ? ((nppm > j) + 1.0) / 2.0
                                                                         : ((std::fabs(nppm - j) > 0.0) | ((j % 2) == 1));
        if (h > 0.0 && h < 1.0)
          percentil = (1.0 - h) * get_nth_element(array, n, j - 1) + h * get_nth_element(array, n, j);
        else
          percentil = get_nth_element(array, n, (h >= 1.0) ? j : j - 1);
      }
    }
  }

  return percentil;
}

template <typename T>
static double
percentile_Rtype8(T *array, size_t len, double quantile)
{
  double rank = 1.0 / 3.0 + (len + 1.0 / 3.0) * quantile;
  size_t k = (size_t) rank;

  double percentil = 0.0;
  if (k == 0) { percentil = get_nth_element(array, len, 0); }
  else if (k >= len) { percentil = get_nth_element(array, len, len - 1); }
  else
  {
    auto vk = get_nth_element(array, len, k - 1);
    auto vk2 = get_nth_element(array, len, k);
    double d = rank - k;
    percentil = vk + d * (vk2 - vk);
  }

  return percentil;
}

static void
percentile_check_number(double pn)
{
  if (pn < 0 || pn > 100) cdo_abort("Percentile number %g out of range! Percentiles must be in the range [0,100].", pn);
}

static void
print_percentile_method(size_t len)
{
  const char *method = "unknown";
  // clang-format off
  if      (percentileMethod == PercentileMethod::NR8)    method = "NR8 (R’s type=8)";
  else if (percentileMethod == PercentileMethod::NRANK)  method = "NRANK (Nearest Rank)";
  else if (percentileMethod == PercentileMethod::NIST)   method = "NIST (recommended by NIST)";
  else if (percentileMethod == PercentileMethod::NUMPY)
  {
    if      (numpyMethod == NumpyMethod::linear)                    method = "NumPy linear";
    else if (numpyMethod == NumpyMethod::lower)                     method = "NumPy lower";
    else if (numpyMethod == NumpyMethod::higher)                    method = "NumPy higher";
    else if (numpyMethod == NumpyMethod::nearest)                   method = "NumPy nearest";
    else if (numpyMethod == NumpyMethod::midpoint)                  method = "NumPy midpoint";
    else if (numpyMethod == NumpyMethod::inverted_cdf)              method = "NumPy inverted_cdf";
    else if (numpyMethod == NumpyMethod::averaged_inverted_cdf)     method = "NumPy averaged_inverted_cdf";
    else if (numpyMethod == NumpyMethod::closest_observation)       method = "NumPy closest_observation";
    else if (numpyMethod == NumpyMethod::interpolated_inverted_cdf) method = "NumPy interpolated_inverted_cdf";
    else if (numpyMethod == NumpyMethod::hazen)                     method = "NumPy hazen";
    else if (numpyMethod == NumpyMethod::weibull)                   method = "NumPy weibull";
    else if (numpyMethod == NumpyMethod::median_unbiased)           method = "NumPy median_unbiased";
    else if (numpyMethod == NumpyMethod::normal_unbiased)           method = "NumPy normal_unbiased";
  }
  // clang-format on

  cdo_print("Using percentile method: %s with %zu values", method, len);
}

template <typename T>
double
percentile(T *array, size_t len, double pn)
{
  static auto printMethod = true;
  if (printMethod && Options::cdoVerbose)
  {
    printMethod = false;
    print_percentile_method(len);
  }

  percentile_check_number(pn);

  double quantile = pn / 100.0;
  double percentil = 0.0;

  // clang-format off
  if      (percentileMethod == PercentileMethod::NR8)    percentil = percentile_Rtype8(array, len, quantile);
  else if (percentileMethod == PercentileMethod::NRANK)  percentil = percentile_nrank(array, len, quantile);
  else if (percentileMethod == PercentileMethod::NIST)   percentil = percentile_nist(array, len, quantile);
  else if (percentileMethod == PercentileMethod::NUMPY)  percentil = percentile_numpy(array, len, quantile);
  else cdo_abort("Internal error: percentile method %d not implemented!", (int)percentileMethod);
  // clang-format on

  return percentil;
}

// Explicit instantiation
template double percentile(float *array, size_t len, double pn);
template double percentile(double *array, size_t len, double pn);

static void
set_numpy_method(NumpyMethod npMethod)
{
  percentileMethod = PercentileMethod::NUMPY;
  numpyMethod = npMethod;
}

void
percentile_set_method(std::string const &methodStr)
{
  auto methodName = string_to_lower(methodStr);

  // clang-format off
  if      ("nrank"   == methodName) percentileMethod = PercentileMethod::NRANK;
  else if ("nist"    == methodName) percentileMethod = PercentileMethod::NIST;
  else if ("rtype8"  == methodName) percentileMethod = PercentileMethod::NR8;
  else if ("numpy"   == methodName) percentileMethod = PercentileMethod::NUMPY;
  else if ("linear"  == methodName || "numpy_linear"  == methodName) set_numpy_method(NumpyMethod::linear);
  else if ("lower"   == methodName || "numpy_lower"   == methodName) set_numpy_method(NumpyMethod::lower);
  else if ("higher"  == methodName || "numpy_higher"  == methodName) set_numpy_method(NumpyMethod::higher);
  else if ("nearest" == methodName || "numpy_nearest" == methodName) set_numpy_method(NumpyMethod::nearest);
  else if ("midpoint"                  == methodName) set_numpy_method(NumpyMethod::midpoint);
  else if ("inverted_cdf"              == methodName) set_numpy_method(NumpyMethod::inverted_cdf);
  else if ("averaged_inverted_cdf"     == methodName) set_numpy_method(NumpyMethod::averaged_inverted_cdf);
  else if ("closest_observation"       == methodName) set_numpy_method(NumpyMethod::closest_observation);
  else if ("interpolated_inverted_cdf" == methodName) set_numpy_method(NumpyMethod::interpolated_inverted_cdf);
  else if ("hazen"                     == methodName) set_numpy_method(NumpyMethod::hazen);
  else if ("weibull"                   == methodName) set_numpy_method(NumpyMethod::weibull);
  else if ("median_unbiased"           == methodName) set_numpy_method(NumpyMethod::median_unbiased);
  else if ("normal_unbiased"           == methodName) set_numpy_method(NumpyMethod::normal_unbiased);
  else cdo_abort("Percentile method %s not available!", methodStr);
  // clang-format on
}

/*
  CDO percentile

#/bin/sh
#
cdo -f nc input,r5x1 testfile <<EOF
 15 20 35 40 50
EOF
cdo -f nc input,r6x1 testfile <<EOF
 15 20 35 40 50 55
EOF
#
cat > per_cdo.sh << EOR
#/bin/sh
CDO=cdo
PERS="30 40 50 75 100"
METS="nrank nist rtype8 numpy numpy_lower numpy_higher numpy_nearest"
METS="linear lower higher nearest midpoint inverted_cdf averaged_inverted_cdf closest_observation interpolated_inverted_cdf hazen
weibull median_unbiased normal_unbiased"
#
for MET in \$METS; do
    echo "\${MET}"
    for PER in \$PERS; do
        $CDO -s --percentile \$MET outputf,%g -fldpctl,\${PER} testfile
    done
done
EOR
sh ./per_cdo.sh > per_cdo_result
*/
/*
  numpy percentile

cat > per_numpy.py << EOR
#python with numpy 1.9.0
import numpy as np
np.version.version
a=np.array([15.0, 20.0, 35.0, 40.0, 50.0, 55.0])
for m in ['linear', 'lower', 'higher', 'nearest', 'midpoint', 'inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
'interpolated_inverted_cdf', 'hazen', 'weibull', 'median_unbiased', 'normal_unbiased'] : print(m) for p in [30, 40, 50, 75, 100] :
print(np.percentile(a, p, method=m)) EOR python per_numpy.py > per_numpy_result
*/
/*
  R percentile

x = c(15.0, 20.0, 35.0, 40.0, 50.0, 55.0)
quantile(x, probs = c(.30, .40, .50, .75, 1.00), type=9)
*/
