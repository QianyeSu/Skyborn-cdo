/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>

#include <vector>

#include <cdi.h>
#include "cdo_cdi_wrapper.h"
#include "cdo_omp.h"
#include "mpim_grid.h"
#include "grid_proj.h"
#include "grid_convert.h"
#include "grid_rot.h"
#include "grid_healpix.h"
#include "gridreference.h"
#include "grid_reduced.h"

#include "compare.h"
#include "cdo_output.h"

bool gridVerbose = false;

void
gridEnableVerbose(bool enable)
{
  gridVerbose = enable;
}

int
nfc_to_nlat(int nfc, int ntr)
{
  return (nfc / (ntr + 1)) / 2;
}

int
nlat_to_ntr(int nlat)
{
  return (nlat * 2 - 1) / 3;
}

int
nlat_to_ntr_linear(int nlat)
{
  return (nlat * 2 - 1) / 2;
}

int
nlat_to_ntr_cubic(int nlat)
{
  return (nlat * 2 - 1) / 4;
}

int
ntr_to_nlat(int ntr)
{
  auto nlat = (int) std::lround((ntr * 3.0 + 1.0) / 2.0);
  if ((nlat % 2) > 0) nlat++;

  return nlat;
}

int
ntr_to_nlat_linear(int ntr)
{
  auto nlat = (int) std::lround((ntr * 2.0 + 1.0) / 2.0);
  if ((nlat % 2) > 0) nlat++;

  return nlat;
}

int
ntr_to_nlat_cubic(int ntr)
{
  auto nlat = (int) std::lround((ntr * 4.0 + 1.0) / 2.0);
  if ((nlat % 2) > 0) nlat++;

  return nlat;
}

int
nlat_to_nlon(int nlat)
{
  return 2 * nlat;
}

int
nlat_to_nlon_cubic(int nlat)
{
  return 2 * nlat + 16;
}

static void
scale_vec(double scalefactor, size_t n, double *values)
{
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) values[i] *= scalefactor;
}

void
grid_copy_names(int gridID1, int gridID2)
{
  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_VDIMNAME, gridID2);
  cdiCopyKey(gridID1, CDI_XAXIS, CDI_KEY_DIMNAME, gridID2);
  cdiCopyKey(gridID1, CDI_YAXIS, CDI_KEY_DIMNAME, gridID2);
  cdiCopyKey(gridID1, CDI_XAXIS, CDI_KEY_NAME, gridID2);
  cdiCopyKey(gridID1, CDI_YAXIS, CDI_KEY_NAME, gridID2);
  cdiCopyKey(gridID1, CDI_XAXIS, CDI_KEY_LONGNAME, gridID2);
  cdiCopyKey(gridID1, CDI_YAXIS, CDI_KEY_LONGNAME, gridID2);
  cdiCopyKey(gridID1, CDI_XAXIS, CDI_KEY_UNITS, gridID2);
  cdiCopyKey(gridID1, CDI_YAXIS, CDI_KEY_UNITS, gridID2);
}

void
grid_copy_mapping(int gridID1, int gridID2)
{
  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, gridID2);
  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gridID2);
  cdiCopyAtts(gridID1, CDI_GLOBAL, gridID2, CDI_GLOBAL);
}

void
grid_to_radian(Varray<double> &values)
{
  scale_vec(DEG2RAD, values.size(), values.data());
}

static void
grid_to_degree(size_t nvals, double *values)
{
  scale_vec(RAD2DEG, nvals, values);
}

void
cdo_grid_to_radian(int gridID, int varID, Varray<double> &values, std::string const &description)
{
  auto llUnits = cdo_grid_get_units(gridID, varID, description);
  if (llUnits == LonLatUnits::Deg) grid_to_radian(values);
}

void
cdo_grid_to_degree(int gridID, int varID, size_t nvals, double *values, std::string const &description)
{
  auto llUnits = cdo_grid_get_units(gridID, varID, description);
  if (llUnits == LonLatUnits::Rad) grid_to_degree(nvals, values);
}

void
cdo_grid_to_degree(int gridID, int varID, Varray<double> &values, std::string const &description)
{
  cdo_grid_to_degree(gridID, varID, values.size(), values.data(), description);
}

int
gridToZonal(int gridID1)
{
  int gridID2 = CDI_UNDEFID;

  auto gridType = gridInqType(gridID1);
  if (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_GAUSSIAN_REDUCED || gridType == GRID_GENERIC)
  {
    if (gridType == GRID_GAUSSIAN_REDUCED) gridType = GRID_GAUSSIAN;

    auto numLats = gridInqYsize(gridID1);
    gridID2 = gridCreate(gridType, numLats);

    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, numLats);

    if (gridType == GRID_GAUSSIAN) gridDefNP(gridID2, gridInqNP(gridID1));

    double xval = 0.0;
    gridDefXvals(gridID2, &xval);

    if (gridInqYvals(gridID1, nullptr))
    {
      std::vector<double> yvals(numLats);
      gridInqYvals(gridID1, yvals.data());
      gridDefYvals(gridID2, yvals.data());
    }
  }
  else if (is_healpix_grid(gridID1))
  {
    auto hpParams = cdo::get_healpix_params(gridID1);
    auto nside = hpParams.nside();

    auto numLats = 4 * nside - 1;
    std::vector<double> latitudes(numLats);

    hp_generate_latitudes(nside, latitudes);
    for (auto &lat : latitudes) { lat *= RAD2DEG; }

    gridID2 = gridCreate(GRID_LONLAT, numLats);

    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, numLats);

    double xval = 0.0;
    gridDefXvals(gridID2, &xval);
    gridDefYvals(gridID2, latitudes.data());
  }
  else
  {
    cdo_abort("Gridtype %s unsupported!", gridNamePtr(gridType));
  }

  return gridID2;
}

int
gridToMeridional(int gridID1)
{
  auto gridType = gridInqType(gridID1);
  auto gridSize = gridInqXsize(gridID1);
  auto gridID2 = gridCreate(gridType, gridSize);

  if (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_GENERIC)
  {
    gridDefXsize(gridID2, gridSize);
    gridDefYsize(gridID2, 1);

    if (gridInqXvals(gridID1, nullptr))
    {
      std::vector<double> xvals(gridSize);
      gridInqXvals(gridID1, xvals.data());
      gridDefXvals(gridID2, xvals.data());
    }

    double yval = 0.0;
    gridDefYvals(gridID2, &yval);
  }
  else
  {
    cdo_abort("Gridtype %s unsupported!", gridNamePtr(gridType));
  }

  return gridID2;
}

void
grid_gen_corners(size_t n, Varray<double> const &vals, Varray<double> &corners)
{
  if (n == 1)
  {
    corners[0] = vals[0];
    corners[1] = vals[0];
  }
  else
  {
    for (size_t i = 0; i < n - 1; ++i) corners[i + 1] = 0.5 * (vals[i] + vals[i + 1]);

    corners[0] = 2 * vals[0] - corners[1];
    corners[n] = 2 * vals[n - 1] - corners[n - 1];
  }
}

void
grid_gen_bounds(size_t n, Varray<double> const &vals, Varray<double> &bounds)
{
  auto lrev = (vals[0] > vals[n - 1]);
  if (lrev)
  {
    for (size_t i = 0; i < n - 1; ++i)
    {
      bounds[2 * i] = 0.5 * (vals[i] + vals[i + 1]);
      bounds[2 * (i + 1) + 1] = 0.5 * (vals[i] + vals[i + 1]);
    }

    bounds[1] = 2 * vals[0] - bounds[0];
    bounds[2 * n - 2] = 2 * vals[n - 1] - bounds[2 * n - 1];
  }
  else
  {
    for (size_t i = 0; i < n - 1; ++i)
    {
      bounds[2 * i + 1] = 0.5 * (vals[i] + vals[i + 1]);
      bounds[2 * (i + 1)] = 0.5 * (vals[i] + vals[i + 1]);
    }

    bounds[0] = 2 * vals[0] - bounds[1];
    bounds[2 * n - 1] = 2 * vals[n - 1] - bounds[2 * (n - 1)];
  }
}

void
grid_check_lat_borders(int n, Varray<double> &ybounds)
{
  constexpr double YMAX = 90.0;
  constexpr double YLIM = 88.0;
  auto lrev = (ybounds[0] > ybounds[n - 1]);
  if (lrev)
  {
    if (ybounds[0] > ybounds[1])
    {
      if (ybounds[0] > YLIM) ybounds[0] = YMAX;
      if (ybounds[n - 1] < -YLIM) ybounds[n - 1] = -YMAX;
    }
    else
    {
      if (ybounds[1] > YLIM) ybounds[1] = YMAX;
      if (ybounds[n - 2] < -YLIM) ybounds[n - 2] = -YMAX;
    }
  }
  else
  {
    if (ybounds[0] < ybounds[1])
    {
      if (ybounds[0] < -YLIM) ybounds[0] = -YMAX;
      if (ybounds[n - 1] > YLIM) ybounds[n - 1] = YMAX;
    }
    else
    {
      if (ybounds[1] < -YLIM) ybounds[1] = -YMAX;
      if (ybounds[n - 2] > YLIM) ybounds[n - 2] = YMAX;
    }
  }
}

/*****************************************************************************/

static void
gridGenCenterRLL(int gridID, size_t nx, size_t ny, Varray<double> const &xvals, Varray<double> const &yvals,
                 Varray<double> &xvals2D, Varray<double> &yvals2D)
{
  double xpole = 0.0, ypole = 0.0, angle = 0.0;
  gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      xvals2D[j * nx + i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
      yvals2D[j * nx + i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
    }
}

static void
gridGenBoundsRLL(int gridID, size_t nx, size_t ny, Varray<double> const &xbounds, Varray<double> const &ybounds,
                 Varray<double> &xbounds2D, Varray<double> &ybounds2D)
{
  double xpole = 0.0, ypole = 0.0, angle = 0.0;
  gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  double minlon, maxlon;
  double minlat, maxlat;

  for (size_t j = 0; j < ny; ++j)
  {
    if (ybounds[0] > ybounds[1])
    {
      maxlat = ybounds[2 * j];
      minlat = ybounds[2 * j + 1];
    }
    else
    {
      maxlat = ybounds[2 * j + 1];
      minlat = ybounds[2 * j];
    }

    for (size_t i = 0; i < nx; ++i)
    {
      minlon = xbounds[2 * i];
      maxlon = xbounds[2 * i + 1];

      size_t index = j * 4 * nx + 4 * i;
      xbounds2D[index + 0] = lamrot_to_lam(minlat, minlon, ypole, xpole, angle);
      xbounds2D[index + 1] = lamrot_to_lam(minlat, maxlon, ypole, xpole, angle);
      xbounds2D[index + 2] = lamrot_to_lam(maxlat, maxlon, ypole, xpole, angle);
      xbounds2D[index + 3] = lamrot_to_lam(maxlat, minlon, ypole, xpole, angle);

      ybounds2D[index + 0] = phirot_to_phi(minlat, minlon, ypole, angle);
      ybounds2D[index + 1] = phirot_to_phi(minlat, maxlon, ypole, angle);
      ybounds2D[index + 2] = phirot_to_phi(maxlat, maxlon, ypole, angle);
      ybounds2D[index + 3] = phirot_to_phi(maxlat, minlon, ypole, angle);
    }
  }
}

void
grid_gen_xbounds2D(size_t nx, size_t ny, Varray<double> const &xbounds, Varray<double> &xbounds2D)
{
#ifdef _OPENMP
  auto n = nx * ny;
#pragma omp parallel for if (n > cdoMinLoopSize) default(none) shared(nx, ny, xbounds, xbounds2D) schedule(static)
#endif
  for (size_t i = 0; i < nx; ++i)
  {
    auto minlon = (xbounds[0] > xbounds[1]) ? xbounds[2 * i + 1] : xbounds[2 * i];
    auto maxlon = (xbounds[0] > xbounds[1]) ? xbounds[2 * i] : xbounds[2 * i + 1];

    for (size_t j = 0; j < ny; ++j)
    {
      auto index = 4 * (j * nx + i);
      xbounds2D[index] = minlon;
      xbounds2D[index + 1] = maxlon;
      xbounds2D[index + 2] = maxlon;
      xbounds2D[index + 3] = minlon;
    }
  }
}

void
grid_gen_ybounds2D(size_t nx, size_t ny, Varray<double> const &ybounds, Varray<double> &ybounds2D)
{
#ifdef _OPENMP
  auto n = nx * ny;
#pragma omp parallel for if (n > cdoMinLoopSize) default(none) shared(nx, ny, ybounds, ybounds2D) schedule(static)
#endif
  for (size_t j = 0; j < ny; ++j)
  {
    auto minlat = (ybounds[0] > ybounds[1]) ? ybounds[2 * j + 1] : ybounds[2 * j];
    auto maxlat = (ybounds[0] > ybounds[1]) ? ybounds[2 * j] : ybounds[2 * j + 1];

    for (size_t i = 0; i < nx; ++i)
    {
      auto index = 4 * (j * nx + i);
      ybounds2D[index] = minlat;
      ybounds2D[index + 1] = minlat;
      ybounds2D[index + 2] = maxlat;
      ybounds2D[index + 3] = maxlat;
    }
  }
}

static void
get_xfirst_and_xlast(int gridID, double &xfirst, double &xlast)
{
  double xfirstandlast[2] = { 0.0, 0.0 };
  gridInqXvals(gridID, xfirstandlast);
  if (is_not_equal(xfirstandlast[0], xfirstandlast[1]))
  {
    xfirst = xfirstandlast[0];
    xlast = xfirstandlast[1];
    if (xfirst > xlast && xfirst > 180.0) xfirst -= 360.0;
  }
}

static bool
reduced_grid_is_global(int np, int nxmax, double xfirst, double xlast)
{
  auto dx_global = (np > 0) ? (90.0 / np) : 999.0;
  auto dx_data = 360.0 - (xlast - xfirst);
  if ((dx_data > dx_global) && (dx_data * nxmax > 360.0)) dx_data = 360.0 / nxmax;
  return !(dx_data > dx_global);
}

void
field2regular(int gridID1, int gridID2, double missval, Varray<double> &array, size_t numMissVals, int lnearest)
{
  auto gridType = gridInqType(gridID1);
  if (gridType != GRID_GAUSSIAN_REDUCED) cdo_abort("Not a reduced Gaussian grid!");

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);
  auto np = gridInqNP(gridID1);

  std::vector<int> reducedPoints(ny);
  gridInqReducedPoints(gridID1, reducedPoints.data());

  double xfirst = 0.0, xlast = 359.9999;
  if (nx == 2) get_xfirst_and_xlast(gridID1, xfirst, xlast);

  int nxmax = 0;
  for (size_t i = 0; i < ny; ++i) nxmax = std::max(nxmax, reducedPoints[i]);

  int lmiss = (numMissVals > 0);
  int lperio = 1;

  int iret;
  if (reduced_grid_is_global(np, nxmax, xfirst, xlast))
  {
    nx = gridInqXsize(gridID2);
    qu2reg3_double(array.data(), reducedPoints.data(), ny, nx, missval, &iret, lmiss, lperio, lnearest);
  }
  else
  {
    nx = qu2reg_subarea(gridInqSize(gridID1), np, xfirst, xlast, array.data(), reducedPoints.data(), ny, missval, &iret, lmiss,
                        lperio, lnearest);
  }

  if (gridInqSize(gridID2) != nx * ny) cdo_abort("Gridsize differ!");
}

int
gridToRegular(int gridID1)
{
  auto gridType = gridInqType(gridID1);
  if (gridType != GRID_GAUSSIAN_REDUCED) cdo_abort("Not a reduced Gaussian grid!");

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);
  auto np = gridInqNP(gridID1);

  std::vector<double> xvals, yvals(ny);
  gridInqYvals(gridID1, yvals.data());

  std::vector<int> reducedPoints(ny);
  gridInqReducedPoints(gridID1, reducedPoints.data());

  double xfirst = 0.0, xlast = 359.9999;
  if (nx == 2) get_xfirst_and_xlast(gridID1, xfirst, xlast);

  int nxmax = 0;
  for (size_t i = 0; i < ny; ++i) nxmax = std::max(nxmax, reducedPoints[i]);

  if (reduced_grid_is_global(np, nxmax, xfirst, xlast))
  {
    nx = reducedPoints[ny / 2];
    if (nx < (2 * ny)) cdo_abort("Number of longitudes %zu is less than 2*ny=%zu!", nx, ny * 2);
    xvals.resize(nx);
    for (size_t i = 0; i < nx; ++i) xvals[i] = xfirst + i * 360.0 / nx;
  }
  else
  {
    if (np <= 0) cdo_abort("Number of values between pole and equator missing!");

    // sub area (longitudes)
    auto np4 = np * 4;
    long ilon_first, ilon_last, row_count;
    grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_first, &ilon_last);

    nx = row_count;
    xvals.resize(nx);
    for (size_t i = 0; i < nx; ++i)
    {
      xvals[i] = ((ilon_first + i) * 360.0) / np4;
      if (xfirst > xlast) xvals[i] -= 360.0;
    }
  }

  auto gridsize = nx * ny;
  auto gridID2 = gridCreate(GRID_GAUSSIAN, gridsize);

  gridDefXsize(gridID2, nx);
  gridDefYsize(gridID2, ny);

  gridDefXvals(gridID2, xvals.data());
  gridDefYvals(gridID2, yvals.data());
  gridDefNP(gridID2, np);

  return gridID2;
}

static void
grid_copy_mask(int gridID1, int gridID2, long gridsize)
{
  if (gridInqMask(gridID1, nullptr))
  {
    std::vector<int> mask(gridsize);
    gridInqMask(gridID1, mask.data());
    gridDefMask(gridID2, mask.data());
  }
}

static bool
check_range(long n, Varray<double> const &vals, double valid_min, double valid_max)
{
  bool status = false;

  for (long i = 0; i < n; ++i)
  {
    if (vals[i] < valid_min || vals[i] > valid_max)
    {
      status = true;
      break;
    }
  }

  return status;
}

bool
grid_has_proj_params(int gridID)
{
  bool hasProjParams = false;

  auto gridType = gridInqType(gridID);
  if (gridType == GRID_PROJECTION)
  {
    int atttype, attlen;
    char attname[CDI_MAX_NAME + 1];

    int natts;
    cdiInqNatts(gridID, CDI_GLOBAL, &natts);

    for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);

      if (atttype == CDI_DATATYPE_TXT)
      {
        if (cdo_cmpstr(attname, "proj_params") || cdo_cmpstr(attname, "proj4_params"))
        {
          hasProjParams = true;
          break;
        }
      }
    }
  }

  return hasProjParams;
}

std::string
grid_get_proj_params(int gridID)
{
  std::string projParams;

  auto gridType = gridInqType(gridID);
  if (gridType == GRID_PROJECTION)
  {
    char attname[CDI_MAX_NAME + 1];

    int natts;
    cdiInqNatts(gridID, CDI_GLOBAL, &natts);

    for (int iatt = 0; iatt < natts; ++iatt)
    {
      int atttype, attlen;
      cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);

      if (atttype == CDI_DATATYPE_TXT)
      {
        std::vector<char> atttxt(attlen + 1);
        cdiInqAttTxt(gridID, CDI_GLOBAL, attname, attlen, atttxt.data());
        atttxt[attlen] = 0;
        if (cdo_cmpstr(attname, "proj_params") || cdo_cmpstr(attname, "proj4_params"))
        {
          projParams = atttxt.data();
          break;
        }
      }
    }
  }

  return projParams;
}

static void
check_units(std::string const &name, std::string const &units)
{
  auto len = units.size();
  auto unitsIsValid = (units == "m" || units == "km" || units.rfind("meter", 0) == 0);

  if (!unitsIsValid)
    cdo_warning("Possibly wrong result! %s %s-coordinate units: %s%s%s (expected \"m\" or \"km\"; default \"m\")",
                len ? "Invalid" : "Missing", name, len ? "\"" : "", units, len ? "\"" : "");
}

static void
center_1D_to_2D(size_t nx, size_t ny, Varray<double> const &xvals, Varray<double> const &yvals, Varray<double> &xvals2D,
                Varray<double> &yvals2D, double xscale, double yscale)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      xvals2D[j * nx + i] = xscale * xvals[i];
      yvals2D[j * nx + i] = yscale * yvals[j];
    }
}

static void
bounds_1D_to_2D(size_t nx, size_t ny, Varray<double> const &xbounds, Varray<double> const &ybounds, Varray<double> &xbounds2D,
                Varray<double> &ybounds2D, double xscale, double yscale)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      auto index = 4 * (j * nx + i);
      xbounds2D[index + 0] = xscale * xbounds[2 * i];
      ybounds2D[index + 0] = yscale * ybounds[2 * j];
      xbounds2D[index + 1] = xscale * xbounds[2 * i];
      ybounds2D[index + 1] = yscale * ybounds[2 * j + 1];
      xbounds2D[index + 2] = xscale * xbounds[2 * i + 1];
      ybounds2D[index + 2] = yscale * ybounds[2 * j + 1];
      xbounds2D[index + 3] = xscale * xbounds[2 * i + 1];
      ybounds2D[index + 3] = yscale * ybounds[2 * j];
    }
}

enum struct Projection
{
  NONE,
  PARAMS,
  RLL,
  LAEA,
  LCC,
  SINU,
  STERE,
  HEALPIX
};

Projection
get_projection(int gridID1)
{
  auto projType = gridInqProjType(gridID1);
  if (projType == CDI_PROJ_RLL) return Projection::RLL;
  if (projType == CDI_PROJ_LAEA) return Projection::LAEA;
  if (projType == CDI_PROJ_LCC) return Projection::LCC;
  if (projType == CDI_PROJ_SINU) return Projection::SINU;
  if (projType == CDI_PROJ_STERE) return Projection::STERE;
  if (projType == CDI_PROJ_HEALPIX) return Projection::HEALPIX;

  auto gmapname = cdo::inq_key_string(gridID1, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME);
  if (gmapname.size())
    cdo_abort("Projection type >%s< unsupported!", gmapname);
  else
    cdo_abort("Projection parameter missing!");

  return Projection::NONE;
}

static void
apply_projection(Projection projection, std::string const &projParams, int gridID, size_t n, Varray<double> &x, Varray<double> &y)
{
  // clang-format off
  if      (projection == Projection::SINU)    cdo_sinu_to_lonlat(n, x, y);
  else if (projection == Projection::LAEA)    cdo_laea_to_lonlat(gridID, n, x, y);
  else if (projection == Projection::LCC)     cdo_lcc_to_lonlat(gridID, n, x, y);
  else if (projection == Projection::STERE)   cdo_stere_to_lonlat(gridID, n, x, y);
  else if (projection == Projection::PARAMS)  cdo_proj_to_lonlat(projParams, n, x, y);
  // clang-format on
}

static bool
is_valid_projection(Projection projection)
{
  return (projection == Projection::LAEA || projection == Projection::LCC || projection == Projection::SINU
          || projection == Projection::STERE);
}

int
gridToCurvilinear(int gridID1, NeedCorners needCorners)
{
  auto withBounds = (needCorners != NeedCorners::No);
  auto genBounds = (needCorners == NeedCorners::Yes);
  auto gridType = gridInqType(gridID1);
  if (!(gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_PROJECTION))
    cdo_abort("%s: Grid type >%s< unsupported!", __func__, gridNamePtr(gridType));

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);

  bool hasGridCoords = gridHasCoordinates(gridID1);
  if (!hasGridCoords) cdo_abort("Grid coordinates missing!");

  auto gridsize = gridInqSize(gridID1);
  auto gridID2 = gridCreate(GRID_CURVILINEAR, gridsize);
  cdiDefKeyInt(gridID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT32);

  Projection projection(Projection::NONE);
  std::string projParams;
  if (gridType == GRID_PROJECTION && gridsize == nx * ny)
  {
    projParams = grid_get_proj_params(gridID1);
    projection = !projParams.empty() ? Projection::PARAMS : get_projection(gridID1);
    // if (projection != Projection::NONE) gridtype = GRID_LONLAT;
  }

  auto isProjection = is_valid_projection(projection);

  auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);
  if (isProjection) check_units("x", xunits);
  if (isProjection) check_units("y", yunits);

  if (isProjection || projection == Projection::PARAMS)
  {
    auto xname = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_NAME);
    auto yname = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_NAME);
    if (xname.size() && yname.size())
    {
      cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_DIMNAME, xname.c_str());
      cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_DIMNAME, yname.c_str());
    }
  }

  double xscale = (xunits[0] == 'k' && xunits[1] == 'm') ? 1000.0 : 1.0;
  double yscale = (yunits[0] == 'k' && yunits[1] == 'm') ? 1000.0 : 1.0;

  gridDefXsize(gridID2, nx);
  gridDefYsize(gridID2, ny);

  Varray<double> xvals2D(gridsize), yvals2D(gridsize);

  if (nx == 0) nx = 1;
  if (ny == 0) ny = 1;

  Varray<double> xvals(nx, 0.0), yvals(ny, 0.0);
  if (gridInqXvals(gridID1, nullptr)) gridInqXvals(gridID1, xvals.data());
  if (gridInqYvals(gridID1, nullptr)) gridInqYvals(gridID1, yvals.data());

  if (projection == Projection::RLL)
  {
    gridDefProj(gridID2, gridID1);
    gridGenCenterRLL(gridID1, nx, ny, xvals, yvals, xvals2D, yvals2D);
  }
  else
  {
    center_1D_to_2D(nx, ny, xvals, yvals, xvals2D, yvals2D, xscale, yscale);
    if (projection != Projection::NONE)
    {
      gridDefProj(gridID2, gridID1);
      apply_projection(projection, projParams, gridID1, gridsize, xvals2D, yvals2D);
    }
  }

  gridDefXvals(gridID2, xvals2D.data());
  gridDefYvals(gridID2, yvals2D.data());

  if (withBounds)
  {
    auto nvertex = (size_t) gridInqNvertex(gridID1);
    Varray<double> xbounds, ybounds;

    if (nvertex == 2 && gridInqXbounds(gridID1, nullptr))
    {
      xbounds.resize(2 * nx);
      gridInqXbounds(gridID1, xbounds.data());
      if (check_range(2 * nx, xbounds, -720, 720))
      {
        cdo_warning("longitude bounds out of range, skipped!");
        xbounds.clear();
      }
    }
    else if (genBounds && nx > 1)
    {
      xbounds.resize(2 * nx);
      grid_gen_bounds(nx, xvals, xbounds);
    }

    if (nvertex == 2 && gridInqYbounds(gridID1, nullptr))
    {
      ybounds.resize(2 * ny);
      gridInqYbounds(gridID1, ybounds.data());
      if (check_range(2 * ny, ybounds, -180, 180))
      {
        cdo_warning("latitude bounds out of range, skipped!");
        ybounds.clear();
      }
    }
    else if (genBounds && ny > 1)
    {
      ybounds.resize(2 * ny);
      if (isProjection || projection == Projection::PARAMS) { grid_gen_bounds(ny, yvals, ybounds); }
      else
      {
        grid_gen_bounds(ny, yvals, ybounds);
        grid_check_lat_borders(2 * ny, ybounds);
      }
    }

    if (xbounds.size() && ybounds.size())
    {
      Varray<double> xbounds2D(4 * gridsize), ybounds2D(4 * gridsize);

      if (projection == Projection::RLL) { gridGenBoundsRLL(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D); }
      else if (isProjection || projection == Projection::PARAMS)
      {
        bounds_1D_to_2D(nx, ny, xbounds, ybounds, xbounds2D, ybounds2D, xscale, yscale);
        apply_projection(projection, projParams, gridID1, 4 * gridsize, xbounds2D, ybounds2D);
      }
      else
      {
        grid_gen_xbounds2D(nx, ny, xbounds, xbounds2D);
        grid_gen_ybounds2D(nx, ny, ybounds, ybounds2D);
      }

      gridDefXbounds(gridID2, xbounds2D.data());
      gridDefYbounds(gridID2, ybounds2D.data());
    }
  }

  grid_copy_mask(gridID1, gridID2, gridsize);

  return gridID2;
}

int
gridToUnstructuredSelecton(int gridID1, std::vector<size_t> const &selectionIndexList, int nocoords, int nobounds)
{
  auto selectionSize = selectionIndexList.size();

  NeedCorners needCorners = nobounds ? NeedCorners::No : NeedCorners::Yes;
  // transform input grid into a unstructured Version if necessary
  auto unstructuredGridID = (GRID_UNSTRUCTURED == gridInqType(gridID1)) ? gridID1 : gridToUnstructured(gridID1, needCorners);

  auto unstructuredGridSize = gridInqSize(unstructuredGridID);

  auto unstructuredSelectionGridID = gridCreate(GRID_UNSTRUCTURED, selectionSize);

  if (nocoords) return unstructuredSelectionGridID;

  // copy meta data of coordinates
  grid_copy_names(unstructuredGridID, unstructuredSelectionGridID);

  // TODO: select bounds

  // copy relevant coordinate
  std::vector<double> xvalsUnstructured(unstructuredGridSize), yvalsUnstructured(unstructuredGridSize);
  gridInqXvals(unstructuredGridID, xvalsUnstructured.data());
  gridInqYvals(unstructuredGridID, yvalsUnstructured.data());

  gridDefXsize(unstructuredSelectionGridID, selectionSize);
  gridDefYsize(unstructuredSelectionGridID, selectionSize);
  std::vector<double> xvals(selectionSize), yvals(selectionSize);

  for (size_t i = 0; i < selectionSize; ++i) xvals[i] = xvalsUnstructured[selectionIndexList[i]];
  for (size_t i = 0; i < selectionSize; ++i) yvals[i] = yvalsUnstructured[selectionIndexList[i]];

  gridDefXvals(unstructuredSelectionGridID, xvals.data());
  gridDefYvals(unstructuredSelectionGridID, yvals.data());

  // copy bounds if requested
  if (!nobounds)
  {
    size_t nvertex = gridInqNvertex(unstructuredGridID);
    std::vector<double> xbounds(nvertex * selectionSize), ybounds(nvertex * selectionSize);
    std::vector<double> xboundsUnstructured(nvertex * unstructuredGridSize), yboundsUnstructured(nvertex * unstructuredGridSize);
    gridInqXbounds(unstructuredGridID, xboundsUnstructured.data());
    gridInqYbounds(unstructuredGridID, yboundsUnstructured.data());
    for (size_t i = 0; i < selectionSize; ++i)
    {
      auto offset = selectionIndexList[i] * nvertex;
      for (size_t k = 0; k < nvertex; ++k) xbounds[i * nvertex + k] = xboundsUnstructured[offset + k];
      for (size_t k = 0; k < nvertex; ++k) ybounds[i * nvertex + k] = yboundsUnstructured[offset + k];
    }
    gridDefNvertex(unstructuredSelectionGridID, nvertex);
    gridDefXbounds(unstructuredSelectionGridID, xbounds.data());
    gridDefYbounds(unstructuredSelectionGridID, ybounds.data());
  }

  return unstructuredSelectionGridID;
}

static void
gridToUnstructuredRegular(int gridID1, int gridID2, size_t gridsize, int withBounds, bool isProjRLL)
{
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_NAME, "lon");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_NAME, "lat");
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_LONGNAME, "longitude");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_LONGNAME, "latitude");
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_UNITS, "degrees_east");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_UNITS, "degrees_north");

  gridDefNvertex(gridID2, 4);

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);

  constexpr double xscale = 1.0;
  constexpr double yscale = 1.0;

  gridDefXsize(gridID2, gridsize);
  gridDefYsize(gridID2, gridsize);

  Varray<double> xvals(nx, 0), yvals(ny, 0);
  if (gridInqXvals(gridID1, nullptr)) gridInqXvals(gridID1, xvals.data());
  if (gridInqXvals(gridID1, nullptr)) gridInqYvals(gridID1, yvals.data());

  {
    Varray<double> xvals2D(gridsize), yvals2D(gridsize);

    if (isProjRLL)
      gridGenCenterRLL(gridID1, nx, ny, xvals, yvals, xvals2D, yvals2D);
    else
      center_1D_to_2D(nx, ny, xvals, yvals, xvals2D, yvals2D, xscale, yscale);

    gridDefXvals(gridID2, xvals2D.data());
    gridDefYvals(gridID2, yvals2D.data());
  }

  if (withBounds)
  {
    size_t nvertex = (size_t) gridInqNvertex(gridID1);
    Varray<double> xbounds, ybounds;

    if (nvertex == 2 && gridInqXbounds(gridID1, nullptr))
    {
      xbounds.resize(2 * nx);
      gridInqXbounds(gridID1, xbounds.data());
    }
    else if (nx > 1)
    {
      xbounds.resize(2 * nx);
      grid_gen_bounds(nx, xvals, xbounds);
    }

    if (nvertex == 2 && gridInqYbounds(gridID1, nullptr))
    {
      ybounds.resize(2 * ny);
      gridInqYbounds(gridID1, ybounds.data());
    }
    else if (ny > 1)
    {
      ybounds.resize(2 * ny);
      grid_gen_bounds(ny, yvals, ybounds);
      grid_check_lat_borders(2 * ny, ybounds);
    }

    if (xbounds.size() && ybounds.size())
    {
      Varray<double> xbounds2D(4 * gridsize), ybounds2D(4 * gridsize);

      if (isProjRLL) { gridGenBoundsRLL(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D); }
      else
      {
        grid_gen_xbounds2D(nx, ny, xbounds, xbounds2D);
        grid_gen_ybounds2D(nx, ny, ybounds, ybounds2D);
      }

      gridDefXbounds(gridID2, xbounds2D.data());
      gridDefYbounds(gridID2, ybounds2D.data());
    }
  }

  grid_copy_mask(gridID1, gridID2, gridsize);
}

static void
gridToUnstructuredGaussianReduced(int gridID1, int gridID2, size_t gridsize, int withBounds)
{
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_NAME, "lon");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_NAME, "lat");
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_LONGNAME, "longitude");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_LONGNAME, "latitude");
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_UNITS, "degrees_east");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_UNITS, "degrees_north");

  gridDefNvertex(gridID2, 4);

  auto nlat = gridInqYsize(gridID1);
  std::vector<int> reducedPoints(nlat);
  gridInqReducedPoints(gridID1, reducedPoints.data());

  gridDefXsize(gridID2, gridsize);
  gridDefYsize(gridID2, gridsize);

  if (nlat == gridInqYvals(gridID1, nullptr))
  {
    std::vector<double> yvals2D(gridsize), yvals(nlat);
    gridInqYvals(gridID1, yvals.data());

    for (size_t j = 0, ij = 0; j < nlat; ++j)
    {
      size_t nlon = reducedPoints[j];
      for (size_t i = 0; i < nlon; ++i) yvals2D[ij++] = yvals[j];
    }

    gridDefYvals(gridID2, yvals2D.data());
  }
  else
  {
    cdo_abort("%s: latitude coordinates missing!", gridNamePtr(gridInqType(gridID1)));
  }

  std::vector<double> xvals(gridsize);

  if (gridsize == gridInqXvals(gridID1, nullptr)) { gridInqXvals(gridID1, xvals.data()); }
  else
  {
    for (size_t j = 0, ij = 0; j < nlat; ++j)
    {
      size_t nlon = reducedPoints[j];
      for (size_t i = 0; i < nlon; ++i) xvals[ij++] = i * 360.0 / nlon;
    }
  }

  gridDefXvals(gridID2, xvals.data());

  if (withBounds)
  {
    auto nvertex = (size_t) gridInqNvertex(gridID1);
    std::vector<double> ybounds;
    if (nvertex == 2 && gridInqYbounds(gridID1, nullptr))
    {
      ybounds.resize(2 * nlat);
      gridInqYbounds(gridID1, ybounds.data());
    }

    if (ybounds.size())
    {
      std::vector<double> xbounds2D(4 * gridsize), ybounds2D(4 * gridsize);

      for (size_t j = 0, ij = 0; j < nlat; ++j)
      {
        for (size_t i = 0, nlon = reducedPoints[j]; i < nlon; ++i)
        {
          xbounds2D[ij + 0] = (i + 0.5) * 360.0 / nlon;
          xbounds2D[ij + 1] = (i + 0.5) * 360.0 / nlon;
          xbounds2D[ij + 2] = (i - 0.5) * 360.0 / nlon;
          xbounds2D[ij + 3] = (i - 0.5) * 360.0 / nlon;
          ybounds2D[ij + 0] = ybounds[j * 2];
          ybounds2D[ij + 1] = ybounds[j * 2 + 1];
          ybounds2D[ij + 2] = ybounds[j * 2 + 1];
          ybounds2D[ij + 3] = ybounds[j * 2];
          ij += 4;
        }
      }

      gridDefXbounds(gridID2, xbounds2D.data());
      gridDefYbounds(gridID2, ybounds2D.data());
    }
  }

  grid_copy_mask(gridID1, gridID2, gridsize);
}

static void
gridToUnstructuredGME(int gridID1, int gridID2, size_t gridsize, int withBounds)
{
  constexpr size_t nv = 6;

  int nd, ni, ni2, ni3;
  gridInqParamGME(gridID1, &nd, &ni, &ni2, &ni3);

  std::vector<int> imask(gridsize);
  std::vector<double> xvals(gridsize), yvals(gridsize);
  std::vector<double> xbounds, ybounds;
  if (withBounds) xbounds.resize(nv * gridsize);
  if (withBounds) ybounds.resize(nv * gridsize);

  gme_grid(withBounds, gridsize, xvals.data(), yvals.data(), xbounds.data(), ybounds.data(), imask.data(), ni, nd, ni2, ni3);

  for (size_t i = 0; i < gridsize; ++i)
  {
    xvals[i] *= RAD2DEG;
    yvals[i] *= RAD2DEG;

    if (withBounds)
      for (size_t j = 0; j < nv; ++j)
      {
        xbounds[i * nv + j] *= RAD2DEG;
        ybounds[i * nv + j] *= RAD2DEG;
      }
    // printf("%d %g %g\n", i, xvals[i], yvals[i]);
  }

  gridDefXsize(gridID2, gridsize);
  gridDefYsize(gridID2, gridsize);

  gridDefXvals(gridID2, xvals.data());
  gridDefYvals(gridID2, yvals.data());

  gridDefMaskGME(gridID2, imask.data());

  gridDefNvertex(gridID2, nv);

  if (withBounds) gridDefXbounds(gridID2, xbounds.data());
  if (withBounds) gridDefYbounds(gridID2, ybounds.data());

  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_UNITS, "degrees_east");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_UNITS, "degrees_north");

  grid_copy_mask(gridID1, gridID2, gridsize);
}

static void
gridToUnstructuredHealpix(int gridID1, int gridID2, size_t gridsize, int withBounds)
{
  constexpr size_t nv = 4;

  Varray<double> xvals(gridsize), yvals(gridsize);
  Varray<double> xbounds, ybounds;
  if (withBounds) xbounds.resize(nv * gridsize);
  if (withBounds) ybounds.resize(nv * gridsize);

  cdo_healpix_to_lonlat(gridID1, gridsize, xvals, yvals, withBounds, xbounds, ybounds);

  gridDefXvals(gridID2, xvals.data());
  gridDefYvals(gridID2, yvals.data());

  if (withBounds)
  {
    gridDefNvertex(gridID2, nv);
    gridDefXbounds(gridID2, xbounds.data());
    gridDefYbounds(gridID2, ybounds.data());
  }

  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_UNITS, "radian");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_UNITS, "radian");
}

int
gridToUnstructured(int gridID1, NeedCorners needCorners)
{
  auto withBounds = (needCorners == NeedCorners::Yes);
  auto gridID2 = CDI_UNDEFID;
  auto gridType = gridInqType(gridID1);
  auto gridSize = gridInqSize(gridID1);

  auto isProjRLL = false;
  auto isProjHealpix = false;
  if (gridType == GRID_PROJECTION)
  {
    auto projType = gridInqProjType(gridID1);
    if (projType == CDI_PROJ_RLL)
    {
      gridType = GRID_LONLAT;
      isProjRLL = true;
    }
    else if (projType == CDI_PROJ_HEALPIX) { isProjHealpix = true; }
    else
    {
      cdo_abort("Projection unsupported!");
    }
  }

  switch (gridType)
  {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    {
      gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
      gridToUnstructuredRegular(gridID1, gridID2, gridSize, withBounds, isProjRLL);
      break;
    }
    case GRID_GAUSSIAN_REDUCED:
    {
      gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
      gridToUnstructuredGaussianReduced(gridID1, gridID2, gridSize, withBounds);
      break;
    }
    case GRID_GME:
    {
      gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
      gridToUnstructuredGME(gridID1, gridID2, gridSize, withBounds);
      break;
    }
    case GRID_CURVILINEAR:
    {
      gridID2 = gridDuplicate(gridID1);
      gridChangeType(gridID2, GRID_UNSTRUCTURED);
      gridDefXsize(gridID2, gridSize);
      gridDefYsize(gridID2, gridSize);
      break;
    }
    case GRID_HEALPIX:
    {
      gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
      gridToUnstructuredHealpix(gridID1, gridID2, gridSize, withBounds);
      break;
    }
    case GRID_PROJECTION:
    {
      if (isProjHealpix)
      {
        gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
        gridToUnstructuredHealpix(gridID1, gridID2, gridSize, withBounds);
      }
      break;
    }
    case GRID_GENERIC:
    {
      gridID2 = gridCreate(GRID_UNSTRUCTURED, gridSize);
      break;
    }
    default:
    {
      cdo_abort("Grid type >%s< unsupported!", gridNamePtr(gridType));
      break;
    }
  }

  cdiDefKeyInt(gridID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT32);

  return gridID2;
}

int
gridCurvilinearToRegular(int gridID1)
{
  int gridID2 = -1;
  auto lx = true, ly = true;

  auto gridType = gridInqType(gridID1);
  auto gridSize = gridInqSize(gridID1);

  if (gridType != GRID_CURVILINEAR) return gridID2;

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);

  std::vector<double> xvals2D(gridSize), yvals2D(gridSize);
  gridInqXvals(gridID1, xvals2D.data());
  gridInqYvals(gridID1, yvals2D.data());

  std::vector<double> xvals(nx), yvals(ny);

  for (size_t i = 0; i < nx; ++i) xvals[i] = xvals2D[i];
  for (size_t j = 0; j < ny; ++j) yvals[j] = yvals2D[j * nx];

  for (size_t j = 1; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      if (std::fabs(xvals[i] - xvals2D[j * nx + i]) > 1.e-6)
      {
        lx = false;
        j = ny;
        break;
      }
    }

  for (size_t i = 1; i < nx; ++i)
    for (size_t j = 0; j < ny; ++j)
    {
      if (std::fabs(yvals[j] - yvals2D[j * nx + i]) > 1.e-6)
      {
        ly = false;
        i = nx;
        break;
      }
    }

  if (lx && ly)
  {
    gridID2 = gridCreate(GRID_LONLAT, gridSize);
    gridDefXsize(gridID2, nx);
    gridDefYsize(gridID2, ny);

    // cdiDefKeyInt(gridID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT32);

    cdo_grid_to_degree(gridID2, CDI_XAXIS, nx, xvals.data(), "grid1 center lon");
    cdo_grid_to_degree(gridID2, CDI_YAXIS, ny, yvals.data(), "grid1 center lat");

    gridDefXvals(gridID2, xvals.data());
    gridDefYvals(gridID2, yvals.data());
  }

  return gridID2;
}

int
gridProjectionToRegular(int gridID1)
{
  int gridID2 = -1;

  auto gridType = gridInqType(gridID1);
  auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);

  if (gridType == GRID_PROJECTION && xunits.rfind("deg", 0) == 0 && yunits.rfind("deg", 0) == 0)
  {
    auto gridsize = gridInqSize(gridID1);
    auto nx = gridInqXsize(gridID1);
    auto ny = gridInqYsize(gridID1);

    std::vector<double> xvals(nx), yvals(ny);
    gridInqXvals(gridID1, xvals.data());
    gridInqYvals(gridID1, yvals.data());

    gridID2 = gridCreate(GRID_LONLAT, gridsize);
    gridDefXsize(gridID2, nx);
    gridDefYsize(gridID2, ny);

    gridDefXvals(gridID2, xvals.data());
    gridDefYvals(gridID2, yvals.data());
  }

  return gridID2;
}

static int
compute_gridcell_weights(int gridID, Varray<double> const &gridCellArea, Varray<double> &gridCellWeights)
{
  int status = 0;
  std::vector<int> gridMask;

  auto gridType = gridInqType(gridID);
  auto gridSize = gridInqSize(gridID);

  if (gridType == GRID_GME)
  {
    gridID = gridToUnstructured(gridID, NeedCorners::Yes);
    gridMask.resize(gridSize);
    gridInqMaskGME(gridID, gridMask.data());
  }

  double totalArea = 0.0;
  // int nvals = 0;
  for (size_t i = 0; i < gridSize; ++i)
  {
    if (gridMask.size())
      if (gridMask[i] == 0) continue;
    totalArea += gridCellArea[i];
    //  nvals++;
  }

  if (gridVerbose) cdo_print("Total area = %g", totalArea);

  for (size_t i = 0; i < gridSize; ++i)
  {
    if (gridMask.size())
      if (gridMask[i] == 0)
      {
        gridCellWeights[i] = 0.0;
        continue;
      }

    gridCellWeights[i] = gridCellArea[i] / totalArea;
  }

  return status;
}

int
gridcell_weights(int gridID, Varray<double> &gridCellWeights)
{
  auto weightStatus = 1;
  auto areaStatus = 0;

  auto gridSize = gridInqSize(gridID);
  Varray<double> gridCellArea(gridSize);

  if (gridHasArea(gridID))
  {
    if (gridVerbose) cdo_print("Using existing grid cell area!");
    gridInqArea(gridID, gridCellArea.data());
  }
  else
  {
    auto gridType = gridInqType(gridID);
    if (gridProjIsSupported(gridID) || gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_GME
        || gridType == GRID_CURVILINEAR || gridType == GRID_UNSTRUCTURED || gridType == GRID_HEALPIX)
    {
      areaStatus = gridGenArea(gridID, gridCellArea);
      if (areaStatus != 0 && (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN))
        areaStatus = gridGenAreaReg2Dweights(gridID, gridCellArea);
    }
    else
    {
      areaStatus = 1;
    }
  }

  if (areaStatus == 0) { weightStatus = compute_gridcell_weights(gridID, gridCellArea, gridCellWeights); }
  else
  {
    for (size_t i = 0; i < gridSize; ++i) gridCellWeights[i] = 1.0 / gridSize;
  }
  /*
  for (i = 0; i < gridsize; ++i)
    printf("weights: %d %d %d %g %g\n", areaStatus, weightStatus, i, gridCellArea[i], gridCellWeights[i]);
  */

  return weightStatus;
}

bool
grid_is_distance_generic(int gridID)
{
  auto status = false;

  if (gridInqType(gridID) == GRID_GENERIC)
  {
    auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
    auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
    if (xunits == "m" && yunits == "m" && gridHasCoordinates(gridID)) status = true;
  }

  return status;
}

bool
is_point_grid(int gridID)
{
  auto gridType = gridInqType(gridID);
  auto projType = gridInqProjType(gridID);

  auto isProjection
      = (gridType == GRID_PROJECTION && (projType == CDI_PROJ_RLL || projType == CDI_PROJ_LCC || projType == CDI_PROJ_STERE));

  return (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_CURVILINEAR || isProjection
          || gridType == GRID_UNSTRUCTURED || gridType == GRID_GME || is_healpix_grid(gridID));
}

static int
generate_full_grid(int gridID, NeedCorners needCorners)
{
  auto gridType = gridInqType(gridID);

  if (gridType == GRID_GME || gridType == GRID_GAUSSIAN_REDUCED || is_healpix_grid(gridID))
  {
    gridID = gridToUnstructured(gridID, needCorners);
    gridType = GRID_UNSTRUCTURED;
  }

  if (gridType != GRID_UNSTRUCTURED && gridType != GRID_CURVILINEAR)
  {
    gridID = gridToCurvilinear(gridID, needCorners);
    gridType = GRID_CURVILINEAR;
  }

  if (gridType == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID))
  {
    auto reference = dereferenceGrid(gridID);
    if (reference.isValid) gridID = reference.gridID;
    if (reference.notFound) cdo_abort("Reference to source grid not found!");
  }

  return gridID;
}

int
generate_full_grid(int gridID)
{
  constexpr auto needCorners = NeedCorners::IfAvail;
  return generate_full_grid(gridID, needCorners);
}

int
generate_full_point_grid(int gridID)
{
  constexpr auto needCorners = NeedCorners::No;
  return generate_full_grid(gridID, needCorners);
}

int
generate_full_cell_grid(int gridID)
{
  constexpr auto needCorners = NeedCorners::Yes;
  return generate_full_grid(gridID, needCorners);
}

static CoordRange
get_range(Varray<double> const &vals)
{
  double vmin = 1.e300;
  double vmax = -1.e300;
  auto hasNan = false;
  auto n = vals.size();

#ifdef HAVE_OPENMP4
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
  for (size_t i = 0; i < n; ++i)
  {
    auto val = vals[i];
    if (std::isnan(val)) { hasNan = true; }
    else
    {
      vmin = std::min(vmin, val);
      vmax = std::max(vmax, val);
    }
  }

  return CoordRange{ vmin, vmax, hasNan };
}

void
check_longitude_range(Varray<double> const &xvals, std::string const &txt, LonLatUnits units)
{
  auto range = get_range(xvals);
  if (range.hasNan) cdo_warning("Grid cell %s longitudes have NaNs!", txt);

  double minLon = (units == LonLatUnits::Deg) ? -360.0 : -2.0 * M_PI;
  double maxLon = (units == LonLatUnits::Deg) ? 720.0 : 2.0 * M_PI;
  if (range.minVal < minLon || range.maxVal > maxLon)
    cdo_warning("Grid cell %s longitudes out of range (min=%g/max=%g)!", txt, range.minVal, range.maxVal);
}

void
check_latitude_range(Varray<double> const &yvals, std::string const &txt, LonLatUnits units)
{
  auto range = get_range(yvals);
  if (range.hasNan) cdo_warning("Grid cell %s latitudes have NaNs!", txt);

  double minLat = (units == LonLatUnits::Deg) ? -90.0001 : -M_PI;
  double maxLat = (units == LonLatUnits::Deg) ? 90.0001 : M_PI;
  if (range.minVal < minLat || range.maxVal > maxLat)
    cdo_warning("Grid cell %s latitudes out of range (min=%g/max=%g)!", txt, range.minVal, range.maxVal);
}

LonLatUnits
string_to_LonLatUnits(std::string const &units, std::string const &description)
{
  if (units.rfind("rad", 0) == 0) { return LonLatUnits::Rad; }
  else if (units.rfind("deg", 0) == 0) { return LonLatUnits::Deg; }

  static bool warn = true;
  if (warn)
  {
    warn = false;
    cdo_warning("Unknown units [%s] supplied for %s; proceeding assuming degress!", units, description);
  }

  return LonLatUnits::Deg;
}

LonLatUnits
cdo_grid_get_units(int gridID, int varID, std::string const &description)
{
  return string_to_LonLatUnits(cdo::inq_key_string(gridID, varID, CDI_KEY_UNITS), description);
}
