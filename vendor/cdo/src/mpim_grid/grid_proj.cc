/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPROJ
#include "proj.h"
#endif

#include <cstdio>
#include <cstdarg> /* va_list */

#include <atomic>
#include <vector>

#include <cdi.h>
#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"
#include "grid_proj.h"
#include "grid_healpix.h"
#include "cdo_output.h"
#include "compare.h"

static void
set_xyvals(double val, size_t nvals, double *xvals, double *yvals)
{
  for (size_t i = 0; i < nvals; ++i)
  {
    xvals[i] = val;
    yvals[i] = val;
  }
}

#ifdef HAVE_LIBPROJ
static void
check_xyvals(size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
  for (size_t i = 0; i < nvals; ++i)
  {
    if (xvals[i] < -9000. || xvals[i] > 9000.) xvals[i] = -9999.;
    if (yvals[i] < -9000. || yvals[i] > 9000.) yvals[i] = -9999.;
  }
}
#endif

#ifdef HAVE_LIBPROJ
static std::string
gen_param(const char *fmt, ...)
{
  constexpr size_t len = 256;
  char str[len];

  va_list args;
  va_start(args, fmt);

  vsnprintf(str, len, fmt, args);

  va_end(args);

  std::string res(str);
  return res;
}

static void
proj_fwd_xyvals(PJ *proj, size_t nvals, double *xvals, double *yvals)
{
  std::atomic<size_t> atomicCountNans{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nvals, proj, atomicCountNans, xvals, yvals) schedule(static)
#endif
  for (size_t i = 0; i < nvals; ++i)
  {
    PJ_COORD p;
    p.uv.u = proj_torad(xvals[i]);
    p.uv.v = proj_torad(yvals[i]);
    p = proj_trans(proj, PJ_FWD, p);
    if (std::isnan(p.uv.u) || std::isnan(p.uv.v)) atomicCountNans++;
    xvals[i] = p.uv.u;
    yvals[i] = p.uv.v;
  }

  size_t countNans = atomicCountNans;
  if (countNans) cdo_warning("%s: %zu of %zu projection coodinates are NaN!", __func__, countNans, 2 * nvals);
}

static void
proj_inv_xyvals(PJ *proj, size_t nvals, double *xvals, double *yvals)
{
  std::atomic<size_t> atomicCountNans{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nvals, proj, atomicCountNans, xvals, yvals) schedule(static)
#endif
  for (size_t i = 0; i < nvals; ++i)
  {
    PJ_COORD p;
    p.uv.u = xvals[i];
    p.uv.v = yvals[i];
    p = proj_trans(proj, PJ_INV, p);
    if (std::isnan(p.uv.u) || std::isnan(p.uv.v)) atomicCountNans++;
    xvals[i] = proj_todeg(p.uv.u);
    yvals[i] = proj_todeg(p.uv.v);
  }

  size_t countNans = atomicCountNans;
  if (countNans) cdo_warning("%s: %zu of %zu projection coodinates are NaN!", __func__, countNans, 2 * nvals);
}

static int
do_proj_fwd(const char *params, size_t nvals, double *xvals, double *yvals)
{
  if (Options::cdoVerbose) cdo_print("Proj fwd: %s", params);

  auto proj = proj_create(PJ_DEFAULT_CTX, params);
  auto status = proj_errno(proj);
  if (status == 0)
  {
    proj_fwd_xyvals(proj, nvals, xvals, yvals);
    proj_destroy(proj);
  }

  return status;
}

static int
do_proj_inv(std::string const &params, size_t nvals, double *xvals, double *yvals)
{
  if (Options::cdoVerbose) cdo_print("Proj inv: %s", params);

  auto proj = proj_create(PJ_DEFAULT_CTX, params.c_str());
  auto status = proj_errno(proj);
  if (status == 0)
  {
    proj_inv_xyvals(proj, nvals, xvals, yvals);
    proj_destroy(proj);
  }

  return status;
}

#endif

static inline bool
checkValIsMiss(std::string const &projection, std::string const &name, double val, double missval)
{
  if (is_equal(val, missval))
  {
    cdo_warning("%s mapping parameter %s missing!", projection, name);
    return true;
  }

  return false;
}

static inline void
checkRangeWarning(const char *projection, const char *name)
{
  cdo_warning("%s mapping parameter %s out of bounds!", projection, name);
}

static inline void
checkLonRange(const char *projection, const char *name, double lon)
{
  if (lon < -360 || lon > 360) checkRangeWarning(projection, name);
}

static inline void
checkLatRange(const char *projection, const char *name, double lat)
{
  if (lat < -90 || lat > 90) checkRangeWarning(projection, name);
}

static inline void
checkRange(const char *projection, const char *name, double val, double missval, double rmin, double rmax)
{
  if (is_not_equal(val, missval) && (val < rmin || val > rmax)) checkRangeWarning(projection, name);
}

static inline void
checkUpperRange(const char *projection, const char *name, double val, double missval, double rmax)
{
  if (is_not_equal(val, missval) && (val > rmax)) checkRangeWarning(projection, name);
}

static void
verify_lcc_parameter(CDI_GridProjParams const &gpp)
{
  constexpr auto projection = "lambert_conformal_conic";

  checkUpperRange(projection, "earth_radius", gpp.a, gpp.mv, 1.e10);
  checkUpperRange(projection, "inverse_flattening", gpp.rf, gpp.mv, 400);

  checkLonRange(projection, "longitude_of_central_meridian", gpp.lon_0);

  checkLatRange(projection, "latitude_of_central_meridian", gpp.lat_0);
  checkLatRange(projection, "standard_parallel", gpp.lat_1);
  checkLatRange(projection, "standard_parallel", gpp.lat_2);

  checkRange(projection, "false_easting", gpp.x_0, gpp.mv, -1.e20, 1.e20);
  checkRange(projection, "false_northing", gpp.y_0, gpp.mv, -1.e20, 1.e20);
}

int
proj_lonlat_to_lcc(struct CDI_GridProjParams gpp, size_t nvals, double *xvals, double *yvals)
{
#ifdef HAVE_LIBPROJ
  std::string params = "+proj=lcc ";
  if (is_not_equal(gpp.a, gpp.mv) && gpp.a > 0) params += gen_param("+a=%.15g ", gpp.a);
  if (is_not_equal(gpp.b, gpp.mv) && gpp.b > 0) params += gen_param("+b=%.15g ", gpp.b);
  if (is_not_equal(gpp.rf, gpp.mv) && gpp.rf > 0) params += gen_param("+rf=%.15g ", gpp.rf);
  params += gen_param("+lon_0=%.15g ", gpp.lon_0);
  params += gen_param("+lat_0=%.15g ", gpp.lat_0);
  params += gen_param("+lat_1=%.15g ", gpp.lat_1);
  params += gen_param("+lat_2=%.15g ", gpp.lat_2);
  params += gen_param("+units=m ");
  //  params += gen_param("+no_defs ");

  int status = do_proj_fwd(params.c_str(), nvals, xvals, yvals);
#else
  int status = 1;
#endif

  if (status == 1) set_xyvals(gpp.mv, nvals, xvals, yvals);

  return status;
}

#ifdef HAVE_LIBPROJ
static void
lonlat_to_lcc(CDI_GridProjParams const &gpp, size_t nvals, double *xvals, double *yvals)
{
  auto status = proj_lonlat_to_lcc(gpp, nvals, xvals, yvals);
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));
}
#endif

int
proj_lcc_to_lonlat(struct CDI_GridProjParams gpp, double x_0, double y_0, size_t nvals, double *xvals, double *yvals)
{
#ifdef HAVE_LIBPROJ
  std::string params = "+proj=lcc ";
  if (is_not_equal(gpp.a, gpp.mv) && gpp.a > 0) params += gen_param("+a=%.15g ", gpp.a);
  if (is_not_equal(gpp.b, gpp.mv) && gpp.b > 0) params += gen_param("+b=%.15g ", gpp.b);
  if (is_not_equal(gpp.rf, gpp.mv) && gpp.rf > 0) params += gen_param("+rf=%.15g ", gpp.rf);
  params += gen_param("+lon_0=%.15g ", gpp.lon_0);
  params += gen_param("+lat_0=%.15g ", gpp.lat_0);
  params += gen_param("+lat_1=%.15g ", gpp.lat_1);
  params += gen_param("+lat_2=%.15g ", gpp.lat_2);
  if (is_not_equal(x_0, gpp.mv)) params += gen_param("+x_0=%.15g ", x_0);
  if (is_not_equal(y_0, gpp.mv)) params += gen_param("+y_0=%.15g ", y_0);

  int status = do_proj_inv(params, nvals, xvals, yvals);
#else
  int status = 1;
#endif

  if (status == 1) set_xyvals(gpp.mv, nvals, xvals, yvals);

  return status;
}

static void
lcc_to_lonlat(CDI_GridProjParams const &gpp, size_t nvals, double *xvals, double *yvals)
{
  auto status = proj_lcc_to_lonlat(gpp, gpp.x_0, gpp.y_0, nvals, xvals, yvals);
#ifdef HAVE_LIBPROJ
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));
#else
  if (status == 1) cdo_abort("proj library support not compiled in!");
#endif
}

int
cdo_lcc_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
  constexpr auto projection = "lambert_conformal_conic";

  CDI_GridProjParams gpp;
  gridInqParamsLCC(gridID, &gpp);

  auto paramIsMissing = checkValIsMiss(projection, "longitude_of_central_meridian", gpp.lon_0, gpp.mv)
                        || checkValIsMiss(projection, "latitude_of_projection_origin", gpp.lat_0, gpp.mv)
                        || checkValIsMiss(projection, "standard_parallel", gpp.lat_1, gpp.mv);

  if (!paramIsMissing && is_equal(gpp.x_0, gpp.mv) && is_equal(gpp.y_0, gpp.mv) && is_not_equal(gpp.xval_0, gpp.mv)
      && is_not_equal(gpp.yval_0, gpp.mv))
  {
#ifdef HAVE_LIBPROJ
    gpp.x_0 = gpp.xval_0;
    gpp.y_0 = gpp.yval_0;
    lonlat_to_lcc(gpp, 1, &gpp.x_0, &gpp.y_0);
    gpp.x_0 = -gpp.x_0;
    gpp.y_0 = -gpp.y_0;
#else
    paramIsMissing = true;
    cdo_warning("%s mapping parameter %s missing!", projection, "false_easting and false_northing");
    cdo_abort("proj library support not compiled in!");
#endif
  }

  if (paramIsMissing) cdo_abort("%s mapping parameter missing!", projection);

  verify_lcc_parameter(gpp);

  lcc_to_lonlat(gpp, nvals, xvals.data(), yvals.data());

  return 0;
}

static void
verify_stere_parameter(CDI_GridProjParams const &gpp)
{
  constexpr auto projection = "polar_stereographic";

  checkUpperRange(projection, "earth_radius", gpp.a, gpp.mv, 1.e10);

  checkLonRange(projection, "straight_vertical_longitude_from_pole", gpp.lon_0);

  checkLatRange(projection, "latitude_of_projection_origin", gpp.lat_0);
  checkLatRange(projection, "standard_parallel", gpp.lat_1);

  checkRange(projection, "false_easting", gpp.x_0, gpp.mv, -1.e20, 1.e20);
  checkRange(projection, "false_northing", gpp.y_0, gpp.mv, -1.e20, 1.e20);
}

int
proj_lonlat_to_stere(struct CDI_GridProjParams gpp, size_t nvals, double *xvals, double *yvals)
{
#ifdef HAVE_LIBPROJ
  std::string params = "+proj=stere ";
  if (is_not_equal(gpp.a, gpp.mv) && gpp.a > 0) params += gen_param("+a=%.15g ", gpp.a);
  if (is_not_equal(gpp.b, gpp.mv) && gpp.b > 0) params += gen_param("+b=%.15g ", gpp.b);
  if (is_not_equal(gpp.rf, gpp.mv) && gpp.rf > 0) params += gen_param("+rf=%.15g ", gpp.rf);
  params += gen_param("+lon_0=%.15g ", gpp.lon_0);
  params += gen_param("+lat_ts=%.15g ", gpp.lat_1);
  params += gen_param("+lat_0=%.15g ", gpp.lat_0);
  params += gen_param("+units=m ");
  //  params += gen_param("+no_defs ");

  int status = do_proj_fwd(params.c_str(), nvals, xvals, yvals);
#else
  int status = 1;
#endif

  if (status == 1) set_xyvals(gpp.mv, nvals, xvals, yvals);

  return status;
}

#ifdef HAVE_LIBPROJ
static void
lonlat_to_stere(CDI_GridProjParams const &gpp, size_t nvals, double *xvals, double *yvals)
{
  auto status = proj_lonlat_to_stere(gpp, nvals, xvals, yvals);
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));
}
#endif

int
proj_stere_to_lonlat(struct CDI_GridProjParams gpp, double x_0, double y_0, size_t nvals, double *xvals, double *yvals)
{
#ifdef HAVE_LIBPROJ
  std::string params = "+proj=stere ";
  if (is_not_equal(gpp.a, gpp.mv) && gpp.a > 0) params += gen_param("+a=%.15g ", gpp.a);
  if (is_not_equal(gpp.b, gpp.mv) && gpp.b > 0) params += gen_param("+b=%.15g ", gpp.b);
  if (is_not_equal(gpp.rf, gpp.mv) && gpp.rf > 0) params += gen_param("+rf=%.15g ", gpp.rf);
  params += gen_param("+lon_0=%.15g ", gpp.lon_0);
  params += gen_param("+lat_ts=%.15g ", gpp.lat_1);
  params += gen_param("+lat_0=%.15g ", gpp.lat_0);
  if (is_not_equal(x_0, gpp.mv)) params += gen_param("+x_0=%.15g ", x_0);
  if (is_not_equal(y_0, gpp.mv)) params += gen_param("+y_0=%.15g ", y_0);

  int status = do_proj_inv(params, nvals, xvals, yvals);
#else
  int status = 1;
#endif

  if (status == 1) set_xyvals(gpp.mv, nvals, xvals, yvals);

  return status;
}

static void
stere_to_lonlat(CDI_GridProjParams const &gpp, size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
  auto status = proj_stere_to_lonlat(gpp, gpp.x_0, gpp.y_0, nvals, xvals.data(), yvals.data());
#ifdef HAVE_LIBPROJ
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));
#else
  if (status == 1) cdo_abort("proj library support not compiled in!");
#endif
}

int
cdo_stere_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
  constexpr auto projection = "polar_stereographic";

  CDI_GridProjParams gpp;
  gridInqParamsSTERE(gridID, &gpp);

  auto paramIsMissing = checkValIsMiss(projection, "straight_vertical_longitude_from_pole", gpp.lon_0, gpp.mv)
                        || checkValIsMiss(projection, "latitude_of_projection_origin", gpp.lat_0, gpp.mv)
                        || checkValIsMiss(projection, "standard_parallel", gpp.lat_1, gpp.mv);

  if (!paramIsMissing && is_equal(gpp.x_0, gpp.mv) && is_equal(gpp.y_0, gpp.mv) && is_not_equal(gpp.xval_0, gpp.mv)
      && is_not_equal(gpp.yval_0, gpp.mv))
  {
#ifdef HAVE_LIBPROJ
    gpp.x_0 = gpp.xval_0;
    gpp.y_0 = gpp.yval_0;
    lonlat_to_stere(gpp, 1, &gpp.x_0, &gpp.y_0);
    gpp.x_0 = -gpp.x_0;
    gpp.y_0 = -gpp.y_0;
#else
    paramIsMissing = true;
    cdo_warning("%s mapping parameter %s missing!", projection, "false_easting and false_northing");
    cdo_abort("proj library support not compiled in!");
#endif
  }

  if (paramIsMissing) cdo_abort("%s mapping parameter missing!", projection);

  verify_stere_parameter(gpp);

  stere_to_lonlat(gpp, nvals, xvals, yvals);

  return 0;
}

int
cdo_healpix_to_lonlat(int gridID, size_t numVals, Varray<double> &xvals, Varray<double> &yvals, bool withBounds,
                      Varray<double> &xbounds, Varray<double> &ybounds)
{
  constexpr auto projection = "healpix";

  auto gridType = gridInqType(gridID);
  auto refinementLevelName = (gridType == GRID_HEALPIX) ? "refinement_level" : "healpix_nside";
  auto indexingSchemeName = (gridType == GRID_HEALPIX) ? "indexing_scheme" : "healpix_order";
  auto refinementLevel = cdo::inq_att_int(gridID, CDI_GLOBAL, refinementLevelName);
  size_t nside = (gridType == GRID_HEALPIX) ? std::lround(std::pow(2.0, refinementLevel)) : refinementLevel;
  auto order = cdo::inq_att_string(gridID, CDI_GLOBAL, indexingSchemeName);
  if (refinementLevel == -1 || order.empty())
  {
    if (order.empty()) cdo_warning("%s mapping parameter %s missing!", projection, indexingSchemeName);
    if (refinementLevel == -1) cdo_warning("%s mapping parameter %s missing!", projection, refinementLevelName);
    cdo_abort("%s mapping parameter missing!", projection);
  }

  auto hpOrder = hp_get_order(order);
  if (hpOrder == HpOrder::Undef) cdo_abort("%s mapping parameter %s=%s unsupported!", projection, indexingSchemeName, order);

  size_t gridSize = 12 * nside * nside;
  if (gridSize == numVals)
  {
    hp_generate_coords(hpOrder, nside, numVals, xvals.data(), yvals.data(), withBounds, xbounds.data(), ybounds.data());
  }
  else
  {
    auto numIndices = gridInqIndices(gridID, nullptr);
    if (numIndices != numVals)
      cdo_abort("Not a global healpix grid (refinementLevel=%d  numCells=%zu/%zu)!", nside, numVals, gridSize);

    std::vector<int64_t> globalIndices(numVals);
    gridInqIndices(gridID, globalIndices.data());
    hp_generate_coords_indices(hpOrder, nside, numVals, xvals.data(), yvals.data(), withBounds, xbounds.data(), ybounds.data(),
                               globalIndices.data());
  }

  return 0;
}

void
grid_def_params_sinu(int gridID)
{
  constexpr auto projection = "sinusoidal";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, projection);
  constexpr auto gmapvarname = "Sinusoidal";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, gmapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) std::strlen(projection), projection);
}

void
grid_def_params_laea(int gridID, double a, double lon_0, double lat_0)
{
  constexpr auto projection = "lambert_azimuthal_equal_area";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, projection);
  constexpr auto gmapvarname = "Lambert_AEA";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, gmapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) std::strlen(projection), projection);

  cdiDefAttFlt(gridID, CDI_GLOBAL, "earth_radius", CDI_DATATYPE_FLT64, 1, &a);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "longitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lat_0);
}

void
cdo_sinu_to_lonlat(size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
#ifdef HAVE_LIBPROJ
  std::string params = "+proj=sinu +ellps=WGS84 ";

  auto status = do_proj_inv(params, nvals, xvals.data(), yvals.data());
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));

  check_xyvals(nvals, xvals, yvals);
#else
  cdo_abort("proj library support not compiled in!");
#endif
}

#ifdef HAVE_LIBPROJ
static bool
cdiInqAttConvertedToFloat(int gridID, int atttype, const char *attname, int attlen, double *attflt)
{
  auto status = true;

  if (atttype == CDI_DATATYPE_INT32)
  {
    std::vector<int> attint(attlen);
    cdiInqAttInt(gridID, CDI_GLOBAL, attname, attlen, &attint[0]);
    for (int i = 0; i < attlen; ++i) attflt[i] = (double) attint[i];
  }
  else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
  {
    cdiInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, attflt);
  }
  else { status = false; }

  return status;
}

static CDI_GridProjParams
grid_inq_params_laea(int gridID)
{
  CDI_GridProjParams gpp;
  gridProjParamsInit(&gpp);

  auto gridtype = gridInqType(gridID);
  if (gridtype == GRID_PROJECTION)
  {
    constexpr auto projection = "lambert_azimuthal_equal_area";
    auto gmapname = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME);
    if (gmapname == projection)
    {
      char attname[CDI_MAX_NAME + 1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      for (int iatt = 0; iatt < natts; ++iatt)
      {
        int atttype, attlen;
        cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
        if (attlen != 1) continue;

        double attflt;
        if (cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, &attflt))
        {
          // clang-format off
                  if      (cdo_cmpstr(attname, "earth_radius"))                   gpp.a     = attflt;
                  else if (cdo_cmpstr(attname, "semi_major_axis"))                gpp.a     = attflt;
                  else if (cdo_cmpstr(attname, "semi_minor_axis"))                gpp.b     = attflt;
                  else if (cdo_cmpstr(attname, "inverse_flattening"))             gpp.rf    = attflt;
                  else if (cdo_cmpstr(attname, "longitude_of_projection_origin")) gpp.lon_0 = attflt;
                  else if (cdo_cmpstr(attname, "latitude_of_projection_origin"))  gpp.lat_0 = attflt;
                  else if (cdo_cmpstr(attname, "false_easting"))                  gpp.x_0   = attflt;
                  else if (cdo_cmpstr(attname, "false_northing"))                 gpp.y_0   = attflt;
          // clang-format on
        }
      }
    }
    else
      cdo_warning("%s mapping parameter missing!", projection);
  }

  return gpp;
}
#endif

void
cdo_laea_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
#ifdef HAVE_LIBPROJ
  auto gpp = grid_inq_params_laea(gridID);

  std::string params = "+proj=laea ";
  if (is_not_equal(gpp.a, gpp.mv) && gpp.a > 0) params += gen_param("+a=%.15g ", gpp.a);
  if (is_not_equal(gpp.b, gpp.mv) && gpp.b > 0) params += gen_param("+b=%.15g ", gpp.b);
  if (is_not_equal(gpp.rf, gpp.mv) && gpp.rf > 0) params += gen_param("+rf=%.15g ", gpp.rf);
  params += gen_param("+lon_0=%.15g ", gpp.lon_0);
  params += gen_param("+lat_0=%.15g ", gpp.lat_0);
  if (is_not_equal(gpp.x_0, gpp.mv)) params += gen_param("+x_0=%.15g ", gpp.x_0);
  if (is_not_equal(gpp.y_0, gpp.mv)) params += gen_param("+y_0=%.15g ", gpp.y_0);

  auto status = do_proj_inv(params, nvals, xvals.data(), yvals.data());
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));

  check_xyvals(nvals, xvals, yvals);
#else
  cdo_abort("proj library support not compiled in!");
#endif
}

void
cdo_proj_to_lonlat(std::string const &projParams, size_t nvals, Varray<double> &xvals, Varray<double> &yvals)
{
#ifdef HAVE_LIBPROJ
  auto status = do_proj_inv(projParams, nvals, xvals.data(), yvals.data());
  if (status) cdo_abort("proj library error: %s", proj_errno_string(status));

  check_xyvals(nvals, xvals, yvals);
#else
  cdo_abort("proj library support not compiled in!");
#endif
}

double
gridGetPlanetRadius(int gridID)
{
  double radius{ 0 };
  return (CDI_NOERR == cdiInqAttFlt(gridID, CDI_GLOBAL, "earth_radius", 1, &radius) && radius > 1) ? radius : 0.0;
}
