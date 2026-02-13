/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Quast

*/

#include <cstdlib>
#include <climits>
#include <algorithm>

#include <cdi.h>

#include "cdo_output.h"
#include "percentiles.h"
#include "percentiles_hist.h"

#define FLT_CAPACITY(n, s) ((int) (((n) * (s)) / sizeof(float)))
#define FLT_PTR(p) ((float *) (p))
#define INT_PTR(p) ((unsigned int *) (p))
#define SHR_PTR(p) ((unsigned short *) (p))

static int
histGetEnvNumBins()
{
  constexpr int NBINS_DEFAULT = 101;
  constexpr int NBINS_MINIMUM = 11;

  auto const *str = getenv("CDO_PCTL_NBINS");

  return (str != nullptr) ? std::max(atoi(str), NBINS_MINIMUM) : NBINS_DEFAULT;
}

template <typename T>
static void
hist_init_bins(int numBins, T *bins)
{
  for (int i = 0; i < numBins; ++i) bins[i] = 0;
}

static void
histDefBounds(HistogramEntry &hist, float a, float b)
{
  assert(hist.numBins > 0);

  hist.nsamp = 0;
  hist.min = std::min(a, b);
  hist.max = std::max(a, b);
  hist.step = (hist.max - hist.min) / hist.numBins;

  if (hist.isUint32)
    hist_init_bins(hist.numBins, INT_PTR(hist.ptr));
  else
    hist_init_bins(hist.numBins, SHR_PTR(hist.ptr));
}

static inline int
calc_bin(int numBins, float histMin, float histStep, float value)
{
  return (histStep > 0.0f) ? std::min((int) ((value - histMin) / histStep), numBins - 1) : 0;
}

template <typename T>
static void
histBinAddValue(HistogramEntry &hist, T *bins, float value)
{
  auto bin = calc_bin(hist.numBins, hist.min, hist.step, value);
  if (bin >= 0 && bin < hist.numBins) bins[bin]++;
}

template <typename T>
static void
histBinSubValue(HistogramEntry &hist, T *bins, float value)
{
  auto bin = calc_bin(hist.numBins, hist.min, hist.step, value);
  if (bin >= 0 && bin < hist.numBins && bins[bin] > 0) bins[bin]--;
}

static void
histBin(HistogramEntry &hist)
{
  assert(hist.nsamp == hist.capacity);

  std::vector<float> values(hist.nsamp);

  auto const *fltptr = FLT_PTR(hist.ptr);
  for (int i = 0; i < hist.nsamp; ++i) values[i] = fltptr[i];

  if (hist.isUint32)
    hist_init_bins(hist.numBins, INT_PTR(hist.ptr));
  else
    hist_init_bins(hist.numBins, SHR_PTR(hist.ptr));

  if (hist.isUint32)
    for (int i = 0; i < hist.nsamp; ++i) histBinAddValue(hist, INT_PTR(hist.ptr), values[i]);
  else
    for (int i = 0; i < hist.nsamp; ++i) histBinAddValue(hist, SHR_PTR(hist.ptr), values[i]);
}
/* unused
static int
histReset(HistogramEntry &hist)
{
  assert(hist.numBins > 0);

  if (hist.nsamp < hist.capacity)
    {
      for (int i = 0; i < hist.nsamp; ++i) FLT_PTR(hist.ptr)[i] = 0.;
    }
  else
    {
      if (hist.isUint32)
        for (int i = 0; i < hist.numBins; ++i) INT_PTR(hist.ptr)[i] = 0;
      else
        for (int i = 0; i < hist.numBins; ++i) SHR_PTR(hist.ptr)[i] = 0;
    }

  hist.nsamp = 0;

  return 0;
}
*/

static void
histCheckValue(const HistogramEntry &hist, float &value)
{
  // 2011-08-01 Uwe Schulzweida: added check for rounding errors
  if (value < hist.min && (hist.min - value) < 1.e5f) value = hist.min;
  if (value > hist.max && (value - hist.max) < 1.e5f) value = hist.max;
}

static int
histAddValue(HistogramEntry &hist, float value)
{
  assert(hist.numBins > 0);

  // if (is_equal(hist.min, hist.max)) return 0;

  histCheckValue(hist, value);
  if (value < hist.min || value > hist.max) return 1;

  if (hist.nsamp < hist.capacity) { FLT_PTR(hist.ptr)[hist.nsamp] = value; }
  else
  {
    if (hist.nsamp == hist.capacity) histBin(hist);

    if (hist.isUint32)
      histBinAddValue(hist, INT_PTR(hist.ptr), value);
    else
      histBinAddValue(hist, SHR_PTR(hist.ptr), value);
  }

  hist.nsamp++;

  return 0;
}

static void
histRemoveValue(HistogramEntry &hist, float value)
{
  auto fltptr = FLT_PTR(hist.ptr);

  int i = 0;
  for (i = 0; i < hist.nsamp; ++i)
  {
    if (is_equal(fltptr[i], value))
    {
      if (i != hist.nsamp - 1) fltptr[i] = fltptr[hist.nsamp - 1];
      break;
    }
  }
  if (i == hist.nsamp)
    cdo_warning("'%f' not found in histogram!", value);
  else
    hist.nsamp--;
}

static int
histSubValue(HistogramEntry &hist, float value)
{
  assert(hist.numBins > 0);

  if (is_equal(hist.min, hist.max)) return 0;

  histCheckValue(hist, value);
  if (value < hist.min || value > hist.max) return 1;

  if (hist.nsamp < hist.capacity) { histRemoveValue(hist, value); }
  else if (hist.nsamp > hist.capacity)
  {
    if (hist.isUint32)
      histBinSubValue(hist, INT_PTR(hist.ptr), value);
    else
      histBinSubValue(hist, SHR_PTR(hist.ptr), value);

    hist.nsamp--;
  }
  else
    return 1;

  return 0;
}

template <typename T>
static double
histGetBin(int nbins, double s, const T *ptr)
{
  int i = 0, count = 0;

  do count += ptr[i++];
  while (count < s);

  assert(i > 0);
  assert((i - 1) < nbins);
  assert(ptr[i - 1] > 0);

  double t = (count - s) / ptr[i - 1];

  return (i - t);
}

static double
histGetPercentile(const HistogramEntry &hist, double p)
{
  assert(hist.nsamp > 0);
  assert(hist.numBins > 0);
  assert(p >= 0.0);
  assert(p <= 100.0);

  if (hist.nsamp > hist.capacity)
  {
    static auto lprint = true;
    if (lprint && Options::cdoVerbose)
    {
      lprint = false;
      cdo_print("Using percentile method: histogram with %d bins", hist.numBins);
    }

    double s = hist.nsamp * (p / 100.0);

    auto bin = hist.isUint32 ? histGetBin(hist.numBins, s, INT_PTR(hist.ptr)) : histGetBin(hist.numBins, s, SHR_PTR(hist.ptr));

    // assert(hist.step > 0.0f);

    return hist.min + bin * hist.step;
  }
  else
  {
    return percentile(FLT_PTR(hist.ptr), hist.nsamp, p);
  }
}

void
HistogramSet::createVarLevels(int varID, int numLevels, size_t numHists)
{
  auto numBins = histGetEnvNumBins();

  assert(numBins > 0);
  assert(numLevels > 0);
  assert(numHists > 0);

  // if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);
  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  this->var_numLevels[varID] = numLevels;
  this->var_numHists[varID] = numHists;
  this->histograms[varID].resize(numLevels);

  for (int levelID = 0; levelID < numLevels; ++levelID)
  {
    this->histograms[varID][levelID].resize(numHists);
    auto &hists = this->histograms[varID][levelID];

    for (size_t histID = 0; histID < numHists; histID++)
    {
      hists[histID].min = 0.0f;
      hists[histID].max = 0.0f;
      hists[histID].step = 0.0f;
      hists[histID].numBins = numBins;
      hists[histID].nsamp = 0;
      auto isUint32 = (this->numSteps >= USHRT_MAX);
      auto vsize = isUint32 ? sizeof(unsigned int) : sizeof(unsigned short);
      hists[histID].isUint32 = isUint32;
      hists[histID].capacity = FLT_CAPACITY(numBins, vsize);
      hists[histID].ptr = std::malloc(numBins * vsize);
      if (hists[histID].ptr == nullptr) cdo_abort("Not enough memory (%s)", __func__);
    }
  }
}

template <typename T1, typename T2>
static void
def_bounds(size_t numHists, std::vector<HistogramEntry> &hists, Varray<T1> const &v1, Varray<T2> const &v2, float mv1, float mv2)
{
  assert(!v1.empty());
  assert(!v2.empty());

  for (size_t i = 0; i < numHists; ++i)
  {
    float a = v1[i];
    float b = v2[i];

    if (fp_is_equal(a, mv1) || fp_is_equal(b, mv2) /*|| fp_is_equal(a, b)*/)
      histDefBounds(hists[i], 0.0f, 0.0f);
    else
      histDefBounds(hists[i], a, b);
  }
}

void
HistogramSet::defVarLevelBounds(int varID, int levelID, Field const &field1, Field const &field2)
{
  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto numLevels = this->var_numLevels[varID];
  if (levelID < 0 || levelID >= numLevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto numHists = this->var_numHists[varID];
  if (numHists != gridInqSize(field1.grid) || numHists != gridInqSize(field2.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  float mv1 = (float) field1.missval;
  float mv2 = (float) field2.missval;

  auto func = [&](auto const &v1, auto const &v2) { def_bounds(numHists, hists, v1, v2, mv1, mv2); };
  field_operation2(func, field1, field2);
}

template <typename T>
static int
histAddVarLevelValues(size_t numHists, std::vector<HistogramEntry> &hists, Varray<T> const &v, size_t numMissVals, double missval)
{
  T mv = missval;
  assert(!v.empty());

  int nign = 0;
  if (numMissVals)
  {
    for (size_t i = 0; i < numHists; ++i)
      if (fp_is_not_equal(v[i], mv)) nign += histAddValue(hists[i], v[i]);
  }
  else
  {
    for (size_t i = 0; i < numHists; ++i) nign += histAddValue(hists[i], v[i]);
  }

  return nign;
}

template <typename T>
static int
histSubVarLevelValues(size_t numHists, std::vector<HistogramEntry> &hists, Varray<T> const &v, size_t numMissVals, double missval)
{
  T mv = missval;
  assert(!v.empty());

  int nign = 0;
  if (numMissVals)
  {
    for (size_t i = 0; i < numHists; ++i)
      if (fp_is_not_equal(v[i], mv)) nign += histSubValue(hists[i], v[i]);
  }
  else
  {
    for (size_t i = 0; i < numHists; ++i) nign += histSubValue(hists[i], v[i]);
  }

  return nign;
}

int
HistogramSet::addVarLevelValues(int varID, int levelID, Field const &field)
{
  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto numLevels = this->var_numLevels[varID];
  if (levelID < 0 || levelID >= numLevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto numHists = this->var_numHists[varID];
  if (numHists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  auto func = [&](auto &v, auto numMissVals, double mv) { return histAddVarLevelValues(numHists, hists, v, numMissVals, mv); };
  auto nign = field_operation(func, field, field.numMissVals, field.missval);
  if (nign)
  {
    cdo_warning("%d out of %d grid values are out of bounds and have been ignored (%s)", nign, numHists, __func__);
    return 1;
  }

  return 0;
}

int
HistogramSet::subVarLevelValues(int varID, int levelID, Field const &field)
{
  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto numLevels = this->var_numLevels[varID];
  if (levelID < 0 || levelID >= numLevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto numHists = this->var_numHists[varID];
  if (numHists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto &hists = this->histograms[varID][levelID];

  auto func = [&](auto &v, auto numMissVals, double mv) { return histSubVarLevelValues(numHists, hists, v, numMissVals, mv); };
  auto nign = field_operation(func, field, field.numMissVals, field.missval);
  if (nign)
  {
    cdo_warning("%d out of %d grid values are out of bounds and have been ignored (%s)", nign, numHists, __func__);
    return 1;
  }

  return 0;
}
/* unused
void
HistogramSet::reset(int varID, int levelID)
{
  assert(numVars > 0);

  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto numLevels = this->var_numLevels[varID];
  assert(numLevels > 0);

  if (levelID < 0 || levelID >= numLevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto numHists = this->var_numHists[varID];
  assert(numHists > 0);

  auto &hists = this->histograms[varID][levelID];
  for (size_t i = 0; i < numHists; ++i) histReset(hists[i]);
}
*/

template <typename T>
static size_t
calcPercentile(size_t numHists, std::vector<HistogramEntry> const &hists, double p, Varray<T> &v, double missval)
{
  T mv = missval;
  assert(!v.empty());

  size_t numMissVals = 0;

  for (size_t i = 0; i < numHists; ++i)
  {
    if (hists[i].nsamp) { v[i] = histGetPercentile(hists[i], p); }
    else
    {
      v[i] = mv;
      numMissVals++;
    }
  }

  return numMissVals;
}

void
HistogramSet::getVarLevelPercentiles(Field &field, int varID, int levelID, double p)
{
  if (varID < 0 || varID >= numVars) cdo_abort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  auto numLevels = this->var_numLevels[varID];
  if (levelID < 0 || levelID >= numLevels) cdo_abort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  auto numHists = this->var_numHists[varID];
  if (numHists != gridInqSize(field.grid)) cdo_abort("Grids are different (%s)", __func__);

  auto const &hists = this->histograms[varID][levelID];

  auto func = [&](auto &v, double mv) { return calcPercentile(numHists, hists, p, v, mv); };
  field.numMissVals = field_operation(func, field, field.missval);
}
