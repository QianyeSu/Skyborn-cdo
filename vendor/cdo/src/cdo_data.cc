#include <cstddef>
// #include <random>

#include "cdo_omp.h"
#include "cdo_data.h"
#include "constants.h"
#include "mpim_grid/grid_convert.h"

static constexpr double etopoScale = 3.0;
static constexpr double etopoOffset = 11000.0;
static constexpr unsigned short etopo[] = {
#include "etopo.dat"
};

static constexpr double tempScale = 500.0;
static constexpr double tempOffset = -220.0;
static constexpr unsigned short temp[] = {
#include "temp.dat"
};

static constexpr double maskScale = 1.0;
static constexpr double maskOffset = 0.0;
static constexpr unsigned short mask[] = {
#include "mask.dat"
};

struct PackedData
{
  double scale{};
  double offset{};
  const unsigned short *pdata = nullptr;
  size_t size{ 0 };
  PackedData(double _scale, double _offset, const unsigned short *_pdata, size_t _size)
      : scale(_scale), offset(_offset), pdata(_pdata), size(_size)
  {
  }
};

namespace cdo
{
const PackedData topoData(etopoScale, etopoOffset, etopo, sizeof(etopo) / sizeof(unsigned short));
const PackedData tempData(tempScale, tempOffset, temp, sizeof(temp) / sizeof(unsigned short));
const PackedData maskData(maskScale, maskOffset, mask, sizeof(mask) / sizeof(unsigned short));

Varray<float>
unpack_data(const PackedData &packedData)
{
  auto n = packedData.size;
  Varray<float> data(n);
  for (size_t i = 0; i < n; ++i) data[i] = packedData.pdata[i] / packedData.scale - packedData.offset;
  return data;
}

void
fill_random(Varray<float> &varray)
{
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> dis(0.0, 1.0);
  // for (auto &v : varray) v = dis(gen);
  for (auto &v : varray) v = ((double) std::rand()) / ((double) RAND_MAX);
}

void
fill_sincos(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals)
{
  auto n = varray.size();
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) { varray[i] = std::cos(1.0 * xvals[i]) * std::sin(2.0 * yvals[i]); }
}

void
fill_coshill(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals)
{
  auto n = varray.size();
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) { varray[i] = 2.0 - std::cos(std::acos(std::cos(xvals[i]) * std::cos(yvals[i])) / 1.2); }
}

void
fill_testfield(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals)
{
  auto n = varray.size();
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i)
  {
    double xyz[3];
    gcLLtoXYZ(xvals[i], yvals[i], xyz);
    auto x = xyz[0];
    auto y = xyz[1];
    auto z = xyz[2];
    varray[i] = 1.0 + std::pow(x, 8.0) + std::exp(2.0 * y * y * y) + std::exp(2.0 * x * x) + 10.0 * x * y * z;
  }
}

// Some Constants for creating temperatur and pressure for the standard atmosphere
constexpr double T_ZERO = 213.0;
constexpr double T_DELTA = 75.0;
constexpr double SCALEHEIGHT = 10000.0;  // [m]

double
std_atm_temperatur(double height)
{
  // Compute the temperatur for the given height (in meters) according to the solution of the hydrostatic atmosphere
  return (T_ZERO + T_DELTA * std::exp((-1) * (height / SCALEHEIGHT)));
}

double
std_atm_pressure(double height)
{
  constexpr double P_ZERO = 1013.25;  // surface pressure [hPa]
  constexpr double CC_R = 287.05;     // specific gas constant for air
  constexpr double TMP4PRESSURE = (C_EARTH_GRAV * SCALEHEIGHT) / (CC_R * T_ZERO);

  // Compute the pressure for the given height (in meters) according to the solution of the hydrostatic atmosphere
  return (P_ZERO
          * std::exp((-1) * TMP4PRESSURE * std::log((std::exp(height / SCALEHEIGHT) * T_ZERO + T_DELTA) / (T_ZERO + T_DELTA))));
}

}  // namespace cdo
