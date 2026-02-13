/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_output.h"
#include "mpim_grid.h"
#include "grid_healpix.h"
#include "gaussian_latitudes.h"
#include "griddes.h"
#include "util_string.h"
#include "dcw_reader.h"

size_t gen_icosphere_coords(int subdivisions, bool withBounds, Varray<double> &xvals, Varray<double> &yvals,
                            Varray<double> &xbounds, Varray<double> &ybounds);

static void
generate_grid_icosphere(GridDesciption &grid, std::string const &gridname)
{
  int gridtype = GRID_UNSTRUCTURED;

  char endChar = '?';
  int bVal = -1, addBounds = -1;
  auto numVals = std::sscanf(gridname.c_str(), "icor2b%d_%d%c", &bVal, &addBounds, &endChar);
  // printf("%s: %d R2B%d_%d<%c>\n", gridname.c_str(), numVals, bVal, addBounds, endChar);

  if (numVals == 0 || numVals > 2) return;
  if (bVal < 0 || bVal > 99) return;
  if (numVals == 2 && addBounds != 0) return;

  auto withBounds = (addBounds == 0) ? false : true;

  grid.type = gridtype;
  if (withBounds) grid.nvertex = 3;

  auto ncells = gen_icosphere_coords(bVal + 1, withBounds, grid.xvals, grid.yvals, grid.xbounds, grid.ybounds);
  grid.xsize = ncells;
  grid.ysize = ncells;
  grid.xname = "clon";
  grid.yname = "clat";
  grid.xunits = "radian";
  grid.yunits = "radian";
}

static void
generate_grid_zonal(GridDesciption &grid, std::string const &gridname, double inc, double lon1, double lon2, double lat1,
                    double lat2)
{
  int gridtype = GRID_LONLAT;

  char endChar = '?';
  double increment = 0.0;
  auto numVals = std::sscanf(gridname.c_str(), "zonal_%lf%c", &increment, &endChar);
  // printf("%s: %d zonal_%g<%c>\n", gridname.c_str(), numVals, increment, endChar);

  if (numVals == 0 || numVals > 1) return;
  if (numVals == 1) inc = increment;

  bool withBounds = true;

  if (inc < 1.0e-9) inc = 1.0;
  if (inc > 180.0) cdo_abort("Zonal latitude increment out if range (max=180)!");

  grid.type = gridtype;

  if (lon1 >= lon2 || lat1 >= lat2) cdo_abort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  auto nlon = 1;
  auto nlat = (size_t) ((lat2 - lat1) / inc + 0.5);

  grid.xvals.resize(nlon, 0.0);
  grid.yvals.resize(nlat);

  for (size_t i = 0; i < nlat; ++i) grid.yvals[i] = lat1 + inc / 2.0 + i * inc;

  grid.xsize = nlon;
  grid.ysize = nlat;

  if (withBounds)
  {
    grid.xbounds.resize(2);
    grid.xbounds[0] = lon1;
    grid.xbounds[1] = lon2;
    grid.ybounds.resize(2 * nlat);
    grid.ybounds[0] = lat1;
    grid.ybounds[1] = lat2;
    if (nlat > 1) grid_gen_bounds(nlat, grid.yvals, grid.ybounds);
  }
}

static void
generate_grid_lonlat(GridDesciption &grid, std::string const &gridname, double inc, double lon1, double lon2, double lat1,
                     double lat2)
{
  int gridtype = GRID_LONLAT;
  bool withBounds = true;

  if (gridname.size())
  {
    char endChar = '?', gridChar = '?';
    double offset = 0.0, increment = 0.0;
    int addBounds = -1;
    auto numVals = std::sscanf(gridname.c_str(), "%lf_%lf_%c%d%c", &offset, &increment, &gridChar, &addBounds, &endChar);
    // printf("%s: %d lonlat%g_%lf_%c%d<%c>\n", gridname.c_str(), numVals, offset, increment, gridChar, addBounds, endChar);

    if (numVals == 0)
    {
      numVals = std::sscanf(gridname.c_str(), "_%lf_%c%d%c", &increment, &gridChar, &addBounds, &endChar);
      // printf("%s: %d lonlat_%lf_%c%d<%c>\n", gridname.c_str(), numVals, increment, gridChar, addBounds, endChar);
      if (numVals < 1 || numVals > 3) return;
    }

    if (numVals < 1 || numVals > 5) return;
    if (gridChar != '?' && gridChar != 'c' && gridChar != 'u') return;
    if (addBounds != -1 && addBounds != 0) return;

    if (addBounds == 0) withBounds = false;

    if (is_not_equal(offset, 0.0))
    {
      lon1 -= offset;
      lon2 += offset;
      lat1 -= offset;
      lat2 += offset;
      lat1 = std::max(lat1, -90.0);
      lat2 = std::min(lat2, 90.0);
    }

    if (is_not_equal(increment, 0.0)) inc = increment;

    if (gridChar == 'c') gridtype = GRID_CURVILINEAR;
    if (gridChar == 'u') gridtype = GRID_UNSTRUCTURED;
  }

  if (inc < 1.0e-9) cdo_abort("LonLat increment >%g< to small (min=1.0e-9)!", inc);
  if (inc > 180.0) cdo_abort("LonLat increment >%g< out if range (max=180)!", inc);

  grid.type = gridtype;

  if (lon1 >= lon2 || lat1 >= lat2) cdo_abort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  auto nlon = (size_t) ((lon2 - lon1) / inc + 0.5);
  auto nlat = (size_t) ((lat2 - lat1) / inc + 0.5);

  grid.xvals.resize(nlon);
  grid.yvals.resize(nlat);

  for (size_t i = 0; i < nlon; ++i) grid.xvals[i] = lon1 + inc * 0.5 + i * inc;
  for (size_t i = 0; i < nlat; ++i) grid.yvals[i] = lat1 + inc * 0.5 + i * inc;

  if (gridtype == GRID_LONLAT)
  {
    grid.xsize = nlon;
    grid.ysize = nlat;
  }
  else
  {
    Varray<double> yvals(nlat);
    for (size_t j = 0; j < nlat; ++j) yvals[j] = grid.yvals[j];
    auto gridsize = nlon * nlat;
    grid.xvals.resize(gridsize);
    grid.yvals.resize(gridsize);
    for (size_t j = 0; j < nlat; ++j)
      for (size_t i = 0; i < nlon; ++i)
      {
        grid.xvals[j * nlon + i] = grid.xvals[i];
        grid.yvals[j * nlon + i] = yvals[j];
      }

    if (gridtype == GRID_CURVILINEAR)
    {
      grid.xsize = nlon;
      grid.ysize = nlat;
    }
    else
    {
      grid.xsize = gridsize;
      grid.ysize = gridsize;
      if (withBounds) grid.nvertex = 4;
    }

    if (withBounds && nlon > 1 && nlat > 1)
    {
      Varray<double> xbounds(2 * nlon), ybounds(2 * nlat);

      grid_gen_bounds(nlon, grid.xvals, xbounds);
      grid_gen_bounds(nlat, yvals, ybounds);
      grid_check_lat_borders(2 * nlat, ybounds);

      grid.xbounds.resize(4 * gridsize);
      grid.ybounds.resize(4 * gridsize);
      grid_gen_xbounds2D(nlon, nlat, xbounds, grid.xbounds);
      grid_gen_ybounds2D(nlon, nlat, ybounds, grid.ybounds);
    }
  }
}

static void
generate_grid_dcw(GridDesciption &grid, std::string const &gridname, double inc)
{
  auto param1 = gridname.c_str();
  auto param2 = std::strstr(param1, "_");
  auto param1len = param2 ? param2 - param1 : gridname.size();

  if (param2)
  {
    char endChar = '?';
    double increment = 0.0;
    auto numVals = std::sscanf(param2, "_%lf%c", &increment, &endChar);
    // printf("%s: %d _%g<%c>\n", param2, numVals, increment, endChar);
    if (numVals == 0 || numVals > 1) return;
    if (numVals == 1) inc = increment;
  }

  const std::string codeNames(string_to_upper({ param1, param1len }));

  DCW_Lists dcw_lists;
  if (dcw_load_lists(dcw_lists)) cdo_abort("dcw_load_lists failed!");

  auto codeList = split_string(codeNames, "\\+");

  dcw_sort_countries(dcw_lists);

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  Region region;
  if (dcw_get_region(dcw_lists, codeList, region)) cdo_abort("dcw_get_region() failed!");

  // printf("lon1, lon2, lat1, lat2 %g %g %g %g\n", region.west, region.east, region.south, region.north);
  auto lon1 = std::round(region.west / inc - 0.5) * inc;
  auto lon2 = std::round(region.east / inc + 0.5) * inc;
  auto lat1 = std::round(region.south / inc - 0.5) * inc;
  auto lat2 = std::round(region.north / inc + 0.5) * inc;
  // printf("lon1, lon2, lat1, lat2 %g %g %g %g\n", lon1, lon2, lat1, lat2);

  const char *param = param2 ? param2 : "";
  generate_grid_lonlat(grid, param, inc, lon1 - inc * 0.5, lon2 + inc * 0.5, lat1 - inc * 0.5, lat2 + inc * 0.5);
}

static void
generate_grid_gme(GridDesciption &grid, std::string const &gridname)
{
  if (gridname.size() == 0) return;

  char endChar = '?';
  int intVal = 0;
  auto numVals = std::sscanf(gridname.c_str(), "%d%c", &intVal, &endChar);
  if (numVals == 0 || numVals > 1) return;
  if (numVals == 1)
  {
    grid.type = GRID_GME;
    grid.ni = intVal;
    grid.nd = 10;
    gme_factorni(grid.ni, &grid.ni2, &grid.ni3);
    grid.size = (grid.ni + 1) * (grid.ni + 1) * 10;
  }
}

static int
scan_healpix_params(std::string const &gridname, bool isZoom, int &refinementLevel, HpOrder &hpOrder, bool &createIndices)
{
  char underChar = '?';
  std::vector<char> orderString(gridname.size(), 0);
  auto numVals = std::sscanf(gridname.c_str(), "%d%c%s", &refinementLevel, &underChar, orderString.data());
  if (refinementLevel < 0) return -1;
  if (isZoom)
  {
    if (numVals == 2 && underChar == 'i')
    {
      createIndices = true;
      return 0;
    }
    if (numVals == 1) return (numVals == 1) ? 0 : -1;
  }

  if (numVals == 0 || numVals == 2 || numVals > 3) return -1;
  if (numVals == 3 && underChar != '_') return -1;

  if (orderString[0])
  {
    hpOrder = hp_get_order(orderString.data());
    if (hpOrder == HpOrder::Undef || hpOrder == HpOrder::XY) return -1;
  }

  return 0;
}

static void
generate_proj_healpix(GridDesciption &grid, std::string const &gridname, bool isZoom = false)
{
  if (gridname.size() == 0) return;

  int refinementLevel = -1;
  HpOrder hpOrder(HpOrder::Nested);
  bool createIndices{ false };
  auto status = scan_healpix_params(gridname, isZoom, refinementLevel, hpOrder, createIndices);
  if (status == -1) return;

  size_t nside = isZoom ? std::lround(std::pow(2, refinementLevel)) : refinementLevel;

  size_t ncells = 12 * nside * nside;
  grid.type = GRID_PROJECTION;
  grid.projection = "healpix";
  grid.size = ncells;
  grid.healpixNside = nside;
  grid.healpixOrder = (hpOrder == HpOrder::Ring) ? "ring" : "nested";
}

static void
generate_grid_healpix(GridDesciption &grid, std::string const &gridname, bool isZoom = false)
{
  if (gridname.size() == 0) return;

  int refinementLevel = -1;
  HpOrder hpOrder(HpOrder::Nested);
  bool createIndices{ false };
  auto status = scan_healpix_params(gridname, isZoom, refinementLevel, hpOrder, createIndices);
  if (status == -1) return;

  size_t nside = isZoom ? std::lround(std::pow(2, refinementLevel)) : refinementLevel;

  size_t ncells = 12 * nside * nside;
  grid.type = GRID_HEALPIX;
  grid.size = ncells;
  grid.refinementLevel = refinementLevel;
  grid.healpixOrder = (hpOrder == HpOrder::Ring) ? "ring" : "nested";
  if (createIndices)
  {
    grid.indices.resize(ncells);
    for (size_t i = 0; i < ncells; ++i) { grid.indices[i] = i; }
  }
}

void
gaussian_latitudes_in_degrees(Varray<double> &lats, Varray<double> &lat_bounds)
{
  auto nlat = lats.size();
  // lats(nlat)
  // lat_bounds(nlat+1)
  std::vector<double> latw(nlat), latw_cumsum(nlat);

  gaussian_latitudes(nlat, lats.data(), latw.data());

  for (size_t j = 0; j < nlat; ++j) lats[j] = RAD2DEG * std::asin(lats[j]);

  latw_cumsum[0] = latw[0];
  for (size_t j = 1; j < nlat; ++j) latw_cumsum[j] = latw_cumsum[j - 1] + latw[j];

  lat_bounds[0] = 1.0;
  for (size_t j = 1; j < nlat; ++j) lat_bounds[j] = 1.0 - latw_cumsum[j - 1];
  lat_bounds[nlat] = -1.0;

  for (size_t j = 0; j < nlat + 1; ++j) lat_bounds[j] = RAD2DEG * std::asin(lat_bounds[j]);
}

static void
generate_grid_gea(GridDesciption &grid, std::string const &gridname)
{
  char endChar = '?';
  double dx = 0.0;
  auto numVals = std::sscanf(gridname.c_str(), "gea%lf%c", &dx, &endChar);
  // printf("%s: %d %lf<%c>\n", gridname.c_str(), numVals, dx, endChar);

  if (numVals != 1) return;
  if (dx < 1.0) return;

  auto dy = dx;
  constexpr auto re = 6378.137;
  constexpr auto f = 1.0 / 298.257223563;
  constexpr auto rp = re * (1.0 - f);
  constexpr auto polar_circumference = 2.0 * M_PI * rp;
  // constexpr auto equator_circumference = 2.0 * M_PI * re;

  size_t nlat = 0.5 * polar_circumference / dy;
  if (nlat % 2) nlat++;

  Varray<double> lats(nlat), lat_bounds(nlat + 1);
  gaussian_latitudes_in_degrees(lats, lat_bounds);

  std::vector<double> cell_height(nlat);
  for (size_t j = 0; j < nlat; ++j) cell_height[j] = 0.25 * polar_circumference * (lat_bounds[j] - lat_bounds[j + 1]) / 90.0;

  // size_t nlone = equator_circumference / dx;
  // if (nlone % 2) nlone++;

  std::vector<int> reducedPoints(nlat);
  size_t ncells = 0;
  for (size_t j = 0; j < nlat; ++j)
  {
    auto rlat = re * std::cos(DEG2RAD * lats[j]);
    auto circumference = 2.0 * M_PI * rlat;
    auto dx_to_use = dx * dy / cell_height[j];

    size_t nlon = std::max((int) std::lround(circumference / dx_to_use), 1);
    if (nlon % 2) nlon++;

    reducedPoints[j] = nlon;
    ncells += nlon;
  }

  // printf("%zu %zu %zu %zu %g\n", ncells, nlone, nlat, nlone*nlat, 100.0*ncells/(nlone*nlat));

  std::vector<double> lons(ncells);
  size_t ij = 0;
  for (size_t j = 0; j < nlat; ++j)
  {
    size_t nlon = reducedPoints[j];
    for (size_t i = 0; i < nlon; ++i) lons[ij++] = i * 360. / nlon;
  }

  grid.type = GRID_GAUSSIAN_REDUCED;
  grid.size = ncells;
  grid.xsize = ncells;
  grid.ysize = nlat;
  grid.numLPE = nlat / 2;
  grid.xvals.resize(ncells);
  grid.yvals.resize(nlat);
  grid.ybounds.resize(nlat * 2);
  grid.reducedPoints.resize(nlat);
  for (size_t i = 0; i < ncells; ++i) grid.xvals[i] = lons[i];
  for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
  for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
  for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
  for (size_t j = 0; j < nlat; ++j) grid.reducedPoints[j] = reducedPoints[j];
}

static void
generate_grid_zonal(GridDesciption &grid, std::string const &gridname)
{
  char endChar = '?';
  int intVal = 0;
  auto numVals = std::sscanf(gridname.c_str(), "z%d%c", &intVal, &endChar);
  // printf("%s: %d %d<%c>\n", gridname.c_str(), numVals, intVal, endChar);

  if (numVals != 1) return;
  if (intVal < 1) return;

  constexpr size_t nextra = 1;
  grid.type = GRID_UNSTRUCTURED;
  auto nlats = (size_t) intVal;
  auto nlons = nlats * 2;
  auto gridsize = nlats;
  grid.size = gridsize;
  grid.xsize = gridsize;
  grid.ysize = gridsize;
  grid.xvals.resize(gridsize);
  grid.yvals.resize(gridsize);
  for (size_t i = 0; i < nlats; ++i) grid.xvals[i] = 180.0;
  auto dlat = 180.0 / nlats;
  auto dlon = 360.0 / nlons;
  // printf("dlat %g dlon %g\n", dlat, dlon);
  for (size_t i = 0; i < nlats; ++i) grid.yvals[i] = -90.0 + i * dlat + dlat / 2.0;
  auto nv = (nlons + 1) * 2;
  grid.nvertex = nv;
  grid.xbounds.resize(nv * gridsize);
  grid.ybounds.resize(nv * gridsize);
  std::vector<double> xbounds(nlons + 1), ybounds(nlats + 1);
  for (size_t i = 0; i <= nlons; ++i) xbounds[i] = 0.0 + i * dlon;
  for (size_t i = 0; i <= nlats; ++i) ybounds[i] = -90.0 + i * dlat;
  // for (size_t i = 0; i <= nlons; ++i)  printf("lon %zu %g\n", i, xbounds[i]);
  // for (size_t i = 0; i <= nlats; ++i)  printf("lat %zu %g\n", i, ybounds[i]);
  size_t k = 0;
  for (size_t j = 0; j < nlats; ++j)
  {
    for (size_t i = nlons; i > 0; i--)
    {
      grid.xbounds[k] = xbounds[i];
      grid.ybounds[k] = ybounds[j + 1];
      k++;
    }
    for (size_t i = 0; i < nextra; ++i)
    {
      grid.xbounds[k] = xbounds[0];
      grid.ybounds[k] = ybounds[j + 1];
      k++;
    }
    for (size_t i = 0; i < nlons; ++i)
    {
      grid.xbounds[k] = xbounds[i];
      grid.ybounds[k] = ybounds[j];
      k++;
    }
    for (size_t i = 0; i < nextra; ++i)
    {
      grid.xbounds[k] = xbounds[nlons];
      grid.ybounds[k] = ybounds[j];
      k++;
    }
  }
}

static void
generate_grid_reg2d(GridDesciption &grid, std::string const &gridname)
{
  char endChar = '?';
  char sepChar = '?';
  int nlon = 0, nlat = 0;
  auto numVals = std::sscanf(gridname.c_str(), "r%d%c%d%c", &nlon, &sepChar, &nlat, &endChar);
  // printf("%s: %d %d%c%d<%c>\n", gridname.c_str(), numVals, nlon, sepChar, nlat, endChar);

  if (numVals != 3 || (sepChar != 'x' && sepChar != '/' && sepChar != '_')) return;

  grid.type = GRID_LONLAT;
  grid.xsize = nlon;
  grid.ysize = nlat;
  grid.xfirst = 0.0;
  grid.yfirst = 0.0;
}

static void
generate_grid_point(GridDesciption &grid, std::string const &gridname)
{
  char endChar = '?';
  char sepChar = '?';
  double lon = 0.0, lat = 0.0;
  auto numVals = std::sscanf(gridname.c_str(), "lon=%lf%clat=%lf%c", &lon, &sepChar, &lat, &endChar);
  // printf("%s: %d lon%lf%clat=%lf<%c>\n", gridname.c_str(), numVals, lon, sepChar, lat, endChar);

  if (numVals != 3 || (sepChar != 'x' && sepChar != '/' && sepChar != '_')) return;

  grid.type = GRID_LONLAT;
  grid.xsize = 1;
  grid.ysize = 1;
  grid.xvals.resize(1);
  grid.yvals.resize(1);
  grid.xvals[0] = lon;
  grid.yvals[0] = lat;
}

static void
generate_grid_generic(GridDesciption &grid, std::string const &gridname)
{
  char endChar = '?';
  char sepChar = '?';
  int nlon = 0, nlat = 0;
  auto numVals = std::sscanf(gridname.c_str(), "g%d%c%d%c", &nlon, &sepChar, &nlat, &endChar);
  // printf("%s: %d %d%c%d<%c>\n", gridname.c_str(), numVals, nlon, sepChar, nlat, endChar);

  if (numVals != 1 && numVals != 3) return;
  if (numVals == 3 && (sepChar != 'x' && sepChar != '/' && sepChar != '_')) return;

  grid.type = GRID_GENERIC;
  grid.xsize = nlon;
  if (numVals == 3) { grid.ysize = nlat; }
  else if (grid.xsize == 1)
  {
    grid.size = 1;
    grid.xsize = 0;
  }
}

static void
generate_grid_gaussian(GridDesciption &grid, std::string gridname)
{
  if (gridname.size())
  {
    int type = 'q';
    if (gridname[0] == 'l')
    {
      type = 'l';
      gridname.erase(0, 1);
    }
    else if (gridname[0] == 'c')
    {
      type = 'c';
      gridname.erase(0, 1);
    }

    if (gridname.size())
    {
      int intValue = 0;
      std::vector<char> typeString(gridname.size(), 0);
      auto numValues = std::sscanf(gridname.c_str(), "%d%s", &intValue, typeString.data());
      // printf("%s: %d %d%s\n", gridname.c_str(), numVals, intVal, typeString.data());

      if (numValues <= 0 || numValues > 2) return;
      if (intValue < 0) return;

      grid.ntr = intValue;
      // clang-format off
      if      (cdo_cmpstr(typeString.data(), "grid")) grid.type = GRID_GAUSSIAN;
      else if (cdo_cmpstr(typeString.data(), "zon"))  grid.type = GRID_GAUSSIAN;
      else if (cdo_cmpstr(typeString.data(), "spec")) grid.type = GRID_SPECTRAL;
      else if (cdo_cmpstr(typeString.data(), ""))     grid.type = GRID_SPECTRAL;
      else return;
      // clang-format on

      if (grid.type == GRID_GAUSSIAN)
      {
        // clang-format off
        if      (type == 'l') grid.ysize = ntr_to_nlat_linear(grid.ntr);
        else if (type == 'c') grid.ysize = ntr_to_nlat_cubic(grid.ntr);
        else                  grid.ysize = ntr_to_nlat(grid.ntr);
        // clang-format on

        grid.numLPE = grid.ysize / 2;
        grid.xsize = (cdo_cmpstr(typeString.data(), "zon"))
                         ? 1
                         : ((type == 'c') ? nlat_to_nlon_cubic(grid.ysize) : nlat_to_nlon(grid.ysize));

        grid.xfirst = 0.0;
        grid.yfirst = 0.0;
        grid.yvals.resize(grid.ysize);
        grid.ybounds.resize(grid.ysize * 2);

        auto nlat = grid.ysize;
        Varray<double> lats(nlat), lat_bounds(nlat + 1);
        gaussian_latitudes_in_degrees(lats, lat_bounds);

        for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
        for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
        for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
      }
    }
  }
}

static void
generate_grid_gaussian_full(GridDesciption &grid, std::string const &gridname)
{
  if (gridname.size() == 0) return;

  // F<N> - full (regular) Gaussian grid with N latitude lines between the pole and equator
  int intVal = 0;
  std::vector<char> typeString(gridname.size(), 0);
  auto numVals = std::sscanf(gridname.c_str(), "%d%s", &intVal, typeString.data());
  // printf("%s: %d %d%s\n", gridname.c_str(), numVals, intVal, typeString.data());

  if (numVals <= 0 || numVals > 2) return;
  if (intVal < 0) return;

  auto numLPE = intVal;

  int index = 0;
  if (typeString[index] == 'b')
  {
    grid.genBounds = true;
    index++;
  }

  if (typeString[index] == '_') index++;

  if (cdo_cmpstr(&typeString[index], "zon"))
  {
    grid.xsize = 1;
    index += 3;
  }

  if (typeString[index] == 0)
  {
    grid.type = GRID_GAUSSIAN;
    grid.numLPE = numLPE;
    grid.ysize = numLPE * 2;
    if (!grid.xsize) grid.xsize = nlat_to_nlon(grid.ysize);

    grid.xfirst = 0.0;
    grid.yfirst = 0.0;
    /* this will change the result of remapcon
    grid.yvals.resize(grid.ysize);
    grid.ybounds.resize(grid.ysize * 2);

    size_t nlat = grid.ysize;
    Varray<double> lats(nlat), lat_bounds(nlat + 1);
    gaussian_latitudes_in_degrees(lats, lat_bounds);

    for (size_t j = 0; j < nlat; ++j) grid.yvals[j] = lats[j];
    for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2 + 1] = lat_bounds[j];
    for (size_t j = 0; j < nlat; ++j) grid.ybounds[j * 2] = lat_bounds[j + 1];
    */
  }
}

static void
generate_grid_gaussian_o(GridDesciption &grid, std::string const &gridname)
{
  // O<N> - octahedral (regular) Gaussian grid with N latitude lines between the pole and equator
  int intVal = 0;
  std::vector<char> typeString(gridname.size(), 0);
  auto numVals = std::sscanf(gridname.c_str(), "o%d%s", &intVal, typeString.data());
  // printf("%s: %d %d%s\n", gridname.c_str(), numVals, intVal, typeString.data());

  if (numVals <= 0 || numVals > 2) return;
  if (intVal < 0) return;

  auto numLPE = intVal;

  int index = 0;
  if (typeString[index] == 'b')
  {
    grid.genBounds = true;
    index++;
  }

  if (typeString[index] == '_') index++;

  if (cdo_cmpstr(&typeString[index], "zon"))
  {
    grid.xsize = 1;
    index += 3;
  }

  if (typeString[index] == 0)
  {
    grid.type = GRID_GAUSSIAN;
    grid.numLPE = numLPE;
    grid.ysize = numLPE * 2;
    if (!grid.xsize) grid.xsize = nlat_to_nlon(grid.ysize) + 16;

    grid.xfirst = 0.0;
    grid.yfirst = 0.0;
  }
}

int
grid_from_name(std::string const &gridname)
{
  int gridID = CDI_UNDEFID;
  if (gridname.size() < 2) return gridID;

  auto gridNameLC = string_to_lower(gridname);

  GridDesciption grid;

  // t<RES>grid or t<RES>spec
  if (gridNameLC[0] == 't') { generate_grid_gaussian(grid, gridNameLC.substr(1)); }
  // r<LON>x<LAT>; regular 2D grid
  else if (gridNameLC[0] == 'r') { generate_grid_reg2d(grid, gridNameLC); }
  // lon=<LON>_lat=<LAT>; one gridpoint
  else if (gridNameLC.starts_with("lon=")) { generate_grid_point(grid, gridNameLC); }
  // gme<NI>
  else if (gridNameLC.starts_with("gme")) { generate_grid_gme(grid, gridNameLC.substr(3)); }
  // ni<NI>
  else if (gridNameLC.starts_with("ni")) { generate_grid_gme(grid, gridNameLC.substr(2)); }
  // healpix  hpr<refinementLevel>[_order] (order=[nested|ring|xy])
  else if (gridNameLC.starts_with("hpr")) { generate_grid_healpix(grid, gridNameLC.substr(3), true); }
  // healpix  hpz<zoom>[_order] (order=[nested|ring|xy])
  else if (gridNameLC.starts_with("hpz")) { generate_proj_healpix(grid, gridNameLC.substr(3), true); }
  // healpix  hp<nside>[_order] (order=[nested|ring|xy])
  else if (gridNameLC.starts_with("hp")) { generate_proj_healpix(grid, gridNameLC.substr(2)); }
  // gea<DX>: gaussian reduced equal area; DX in km
  else if (gridNameLC.starts_with("gea")) { generate_grid_gea(grid, gridNameLC); }
  // F<N> - full (regular) Gaussian grid with N latitude lines between the pole and equator
  else if ((gridNameLC[0] == 'f' || gridNameLC[0] == 'n') && std::isdigit((int) gridNameLC[1]))
  {
    generate_grid_gaussian_full(grid, gridNameLC.substr(1));
  }
  // O<N>
  else if (gridNameLC[0] == 'o' && std::isdigit((int) gridNameLC[1])) { generate_grid_gaussian_o(grid, gridNameLC); }
  // g<LON>x<LAT> or g<SIZE>
  else if (gridNameLC[0] == 'g' && std::isdigit(gridNameLC[1])) { generate_grid_generic(grid, gridNameLC); }
  // dcw:code_Xdeg
  else if (gridNameLC.starts_with("dcw:")) { generate_grid_dcw(grid, gridNameLC.substr(4), 0.1); }
  // germany_Xdeg
  else if (gridNameLC.starts_with("germany")) { generate_grid_lonlat(grid, gridNameLC.substr(7), 0.1, 5.6, 15.2, 47.1, 55.1); }
  // europe_Xdeg
  else if (gridNameLC.starts_with("europe")) { generate_grid_lonlat(grid, gridNameLC.substr(6), 1, -30, 60, 30, 80); }
  // africa_Xdeg
  else if (gridNameLC.starts_with("africa")) { generate_grid_lonlat(grid, gridNameLC.substr(6), 1, -20, 60, -40, 40); }
  // global_Xdeg
  else if (gridNameLC.starts_with("global")) { generate_grid_lonlat(grid, gridNameLC.substr(6), 1, -180, 180, -90, 90); }
  // zonal_Xdeg
  else if (gridNameLC.starts_with("zonal")) { generate_grid_zonal(grid, gridNameLC, 1, -180, 180, -90, 90); }
  // z<LAT>; zonal unstructured grid with <LAT> latitudes
  else if (gridNameLC[0] == 'z') { generate_grid_zonal(grid, gridNameLC); }
  // icoR02BXX
  else if (gridNameLC.starts_with("ico")) { generate_grid_icosphere(grid, gridNameLC); }

  if (grid.type != -1) gridID = grid_define(grid);

  return gridID;
}
