/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Getgridcell     gridcellindex    Get grid cell index
*/

// #include <algorithm>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include <mpim_grid.h>
#include "grid_healpix.h"
#include "grid_pointsearch.h"

namespace
{
struct Parameter
{
  double lon = 0.0;
  double lat = 0.0;
  double arc_radius = 0.0;
  double radius = 1.0;
};
}  // namespace

static size_t
lonlat_to_index(int gridID, Parameter const &params)
{
  auto gridID0 = gridID;
  auto gridsize = gridInqSize(gridID);

  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yvals, "grid center lat");

  GridPointsearch gps;
  gps.set_radius((params.arc_radius > 0.0) ? arc_to_chord_length(params.arc_radius) : params.radius);
  grid_pointsearch_create_unstruct(gps, xvals, yvals, true);

  constexpr size_t numNeighbors = 1;
  KnnData knnData(numNeighbors);
  grid_search_point_unstruct(gps, PointLonLat{ deg_to_rad(params.lon), deg_to_rad(params.lat) }, knnData);
  auto cellIdx = knnData.m_indices[0];

  if (gridID0 != gridID) gridDestroy(gridID);

  return cellIdx;
}

static Parameter
get_parameter()
{
  Parameter params{};

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "lon")        params.lon = parameter_to_double(value);
      else if (key == "lat")        params.lat = parameter_to_double(value);
      else if (key == "radius")     params.radius = radius_str_to_deg(value);
      else if (key == "arc_radius") params.arc_radius = radius_str_to_deg(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
print_parameter(Parameter const &params)
{
  std::stringstream outbuffer;
  outbuffer << "lon=" << params.lon << ", lat=" << params.lat;
  cdo_print("%s", outbuffer.str());
}

static void
check_radius_range(double radius, const char *name)
{
  if (radius < 0.0 || radius > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", name, radius);
}

class Getgridcell : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Getgridcell",
    .operators = { { "gridcellindex", 0, 0, "lon/lat coordinate of a single cell", GetgridcellHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Getgridcell> registration = RegisterEntry<Getgridcell>();

  CdoStreamID streamID1{};
  Parameter params{};
  int gridID1{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }

    params = get_parameter();

    check_radius_range(params.radius, "radius");
    check_radius_range(params.arc_radius, "arc_radius");

    if (Options::cdoVerbose) print_parameter(params);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids != 1) cdo_abort("Too many different grids!");

    gridID1 = vlistGrid(vlistID1, 0);
  }

  void
  run() override
  {
    int64_t cellIdx = -1;
    if (is_healpix_grid(gridID1))
    {
      cellIdx = hp_lonlat_to_index(cdo::get_healpix_params(gridID1), deg_to_rad(params.lon), deg_to_rad(params.lat));
    }
    else { cellIdx = lonlat_to_index(gridID1, params); }
    {
    }

    printf("%ld\n", (long) cellIdx + 1);
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
  }
};
