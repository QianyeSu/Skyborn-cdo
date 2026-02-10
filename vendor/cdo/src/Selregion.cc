/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include <cstdint>

#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "grid_pointsearch.h"
#include "mpim_grid.h"
#include "region.h"

int gengridcell(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices);
void window_cell(Field const &field1, Field &field2, std::vector<int64_t> const &cellIndices);
double radiusDegToKm(double radiusInDeg);

namespace
{
struct CirclePoint
{
  double radius{ 1.0 };
  double lon{ 0.0 };
  double lat{ 0.0 };
  size_t maxpoints{ SIZE_MAX };
};

struct RegionInfo
{
  std::vector<int64_t> cellIndices;
  long nvals{ 0 };
  int gridtype{ -1 };
  int gridID1{ -1 };
  int gridID2{ -1 };
};
}  // namespace

static inline bool
is_point_inside(double xval, double yval, double xi, double xj, double yi, double yj)
{
  return (((yval >= yi && yval < yj) || (yval > yj && yval <= yi)) && (xval < ((xj - xi) * (yval - yi) / (yj - yi) + xi)));
}

static bool
point_is_inside(double xval, double yval, size_t n, const double *xcoords, const double *ycoords)
{
  auto c = false;

  for (size_t i = 0, j = n - 1; i < n; j = i++)
  {
    if (is_point_inside(xval, yval, xcoords[i], xcoords[j], ycoords[i], ycoords[j])) c = !c;
  }

  return c;
}

static bool
point_is_inside(double xval, double yval, double xmin, double xmax, const double *xcoords, const double *ycoords, size_t nofcoords)
{
  auto c = false;

  // clang-format off
  if      (xval >= xmin && xval <= xmax)
    c = point_is_inside(xval,         yval, nofcoords, xcoords, ycoords);
  else if (xval > 180.0 && xval - 360.0 >= xmin && xval - 360.0 <= xmax)
    c = point_is_inside(xval - 360.0, yval, nofcoords, xcoords, ycoords);
  else if (xval <   0.0 && xval + 360.0 >= xmin && xval + 360.0 <= xmax)
    c = point_is_inside(xval + 360.0, yval, nofcoords, xcoords, ycoords);
  // clang-format on

  return c;
}

static void
sel_region_cell(Vmask &mask, size_t gridsize, Varray<double> const &xvals, Varray<double> const &yvals, const double *xcoords,
                const double *ycoords, size_t segmentSize, std::vector<int64_t> &cellIndices)
{
  auto xmm = varray_min_max(segmentSize, xcoords);
  auto ymm = varray_min_max(segmentSize, ycoords);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (mask[i]) continue;

    auto yval = yvals[i];
    if (yval > ymm.min && yval < ymm.max)
    {
      if (point_is_inside(xvals[i], yval, xmm.min, xmm.max, xcoords, ycoords, segmentSize)) mask[i] = true;
    }
  }

  for (size_t i = 0; i < gridsize; ++i)
  {
    if (mask[i]) cellIndices.push_back(i);
  }
}

static int
generate_region_grid(int gridID1, long &gridsize2, std::vector<int64_t> &cellIndices, int numFiles)
{
  auto gridID0 = gridID1;

  gridID1 = generate_full_grid(gridID1);
  if (!gridHasCoordinates(gridID1)) cdo_abort("Cell center coordinates missing!");

  auto gridsize = gridInqSize(gridID1);
  Varray<double> xvals(gridsize), yvals(gridsize);

  gridInqXvals(gridID1, xvals.data());
  gridInqYvals(gridID1, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID1, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID1, CDI_YAXIS, yvals, "grid center lat");

  Vmask mask(gridsize, false);
  for (int i = 0; i < numFiles; ++i)
  {
    Regions regions;
    auto param = cdo_operator_argv(i);
    if (param.starts_with("dcw:"))
      read_regions_from_dcw(param.c_str() + 4, regions);
    else
      read_regions_from_file(param, regions);

    for (size_t k = 0; k < regions.numSegments; ++k)
    {
      auto segmentSize = regions.segmentSize[k];
      if (segmentSize < 3) continue;
      auto offset = regions.segmentOffset[k];
      const auto xcoords = &regions.x[offset];
      const auto ycoords = &regions.y[offset];
      sel_region_cell(mask, gridsize, xvals, yvals, xcoords, ycoords, segmentSize, cellIndices);
    }
  }

  gridsize2 = cellIndices.size();
  if (gridsize2 == 0) cdo_abort("No grid points found!");

  auto gridID2 = gridsize2 ? gengridcell(gridID1, gridsize2, cellIndices) : CDI_UNDEFID;

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static int
generate_circle_grid(int gridID1, long &gridsize2, std::vector<int64_t> &cellIndices, const CirclePoint &cpoint)
{
  auto gridID0 = gridID1;

  gridID1 = generate_full_grid(gridID1);
  if (!gridHasCoordinates(gridID1)) cdo_abort("Cell center coordinates missing!");

  {
    auto gridsize1 = gridInqSize(gridID1);

    Varray<double> xvals(gridsize1), yvals(gridsize1);
    gridInqXvals(gridID1, xvals.data());
    gridInqYvals(gridID1, yvals.data());

    // Convert lat/lon units if required
    cdo_grid_to_radian(gridID1, CDI_XAXIS, xvals, "grid center lon");
    cdo_grid_to_radian(gridID1, CDI_YAXIS, yvals, "grid center lat");

    GridPointsearch gps;
    gps.set_radius(arc_to_chord_length(cpoint.radius));
    grid_pointsearch_create_unstruct(gps, xvals, yvals);

    auto numNeighbors = cpoint.maxpoints;
    if (numNeighbors > gridsize1) numNeighbors = gridsize1;

    KnnData knnData(numNeighbors);
    grid_search_point_smooth(gps, PointLonLat{ cpoint.lon, cpoint.lat }, knnData);

    auto nvals = knnData.m_numNeighbors;
    cellIndices.resize(nvals);

    for (size_t i = 0; i < nvals; ++i) cellIndices[i] = knnData.m_indices[i];

    if (nvals == 0) cdo_abort("No grid points found!");

    gridsize2 = nvals;
  }

  auto gridID2 = gengridcell(gridID1, gridsize2, cellIndices);

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static void
selcircle_get_parameter(CirclePoint &cpoint)
{
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
      if      (key == "maxpoints") cpoint.maxpoints = parameter_to_size_t(value);
      else if (key == "lon")       cpoint.lon = parameter_to_double(value);
      else if (key == "lat")       cpoint.lat = parameter_to_double(value);
      else if (key == "radius")    cpoint.radius = radius_str_to_deg(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Selregion : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selregion",
    .operators = { { "selregion", 0, 0, "DCW region or the path to region file", SelregionHelp },
                   { "selcircle", 0, 0, "Longitude, latitude of the center and radius of the circle", SelregionHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Selregion>(module);

  int SELREGION{}, SELCIRCLE{};
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numGrids{};

  VarList varList1;
  VarList varList2;
  std::vector<bool> varIDs;
  std::vector<RegionInfo> regions;

public:
  void
  init() override
  {
    SELREGION = module.get_id("selregion");
    SELCIRCLE = module.get_id("selcircle");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numVars = varList1.numVars();
    varIDs = std::vector<bool>(numVars, false);

    numGrids = varList1.numGrids();
    regions = std::vector<RegionInfo>(numGrids);

    int numFiles = 0;
    CirclePoint cpoint;

    if (operatorID == SELREGION)
    {
      numFiles = cdo_operator_argc();
      if (numFiles == 0) cdo_abort("Region parameter missing!");
    }
    else if (operatorID == SELCIRCLE)
    {
      selcircle_get_parameter(cpoint);
      if (cpoint.radius < 0.0 || cpoint.radius > 180.0) cdo_abort("radius=%g out of bounds (0-180 deg)!", cpoint.radius);

      if (varList1.gridsizeMax() < cpoint.maxpoints) cpoint.maxpoints = varList1.gridsizeMax();
      if (Options::cdoVerbose)
        cdo_print("lon = %g, lat = %g, radius = %gdeg(%gkm)", cpoint.lon, cpoint.lat, cpoint.radius, radiusDegToKm(cpoint.radius));

      cpoint.radius *= DEG2RAD;
      cpoint.lon *= DEG2RAD;
      cpoint.lat *= DEG2RAD;
    }

    for (int index = 0; index < numGrids; ++index)
    {
      auto &region = regions[index];
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      if (is_point_grid(gridID1))
      {
        auto gridsize = gridInqSize(gridID1);
        if (gridsize == 1) continue;

        region.cellIndices.reserve(gridsize);

        int gridID2 = CDI_UNDEFID;
        if (operatorID == SELREGION)
          gridID2 = generate_region_grid(gridID1, region.nvals, region.cellIndices, numFiles);
        else if (operatorID == SELCIRCLE)
          gridID2 = generate_circle_grid(gridID1, region.nvals, region.cellIndices, cpoint);

        region.cellIndices.shrink_to_fit();

        if (gridID2 != CDI_UNDEFID)
        {
          region.gridtype = gridtype;
          region.gridID1 = gridID1;
          region.gridID2 = gridID2;

          vlistChangeGridIndex(vlistID2, index, gridID2);

          for (auto const &var : varList1.vars)
            if (gridID1 == var.gridID) varIDs[var.ID] = true;
        }
      }
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

    varList2 = VarList(vlistID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field1, field2;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID1, field1);

        cdo_def_field(streamID2, varID, levelID);

        if (varIDs[varID])
        {
          auto gridID1 = var.gridID;
          int index;
          for (index = 0; index < numGrids; ++index)
            if (gridID1 == regions[index].gridID1) break;
          if (index == numGrids) cdo_abort("Internal problem, grid not found!");

          field2.init(varList2.vars[varID]);
          window_cell(field1, field2, regions[index].cellIndices);

          if (field1.numMissVals) field_num_mv(field2);

          cdo_write_field(streamID2, field2);
        }
        else { cdo_write_field(streamID2, field1); }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
