/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include <cstddef>
#include <cstdint>

#include "cdo_output.h"
#include "cdo_omp.h"
#include "cdo_timer.h"
#include "griddes.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "grid_pointsearch.h"
#include "mpim_grid.h"
#include "util_string.h"

int gengridcell(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices);
void window_cell(Field const &field1, Field &field2, std::vector<int64_t> const &cellIndices);
double radiusDegToKm(double radiusInDeg);

enum struct LatOrientation
{
  Undefined,
  Positive,
  Negative
};

namespace
{
struct Parameter
{
  LatOrientation data{ LatOrientation::Undefined };
  LatOrientation lat{ LatOrientation::Positive };
  std::string gridFile;
};

struct GridInfo
{
  std::vector<int64_t> cellIndices;
  size_t gridSize{ 0 };
  int gridType{ -1 };
  int gridID{ -1 };
};
}  // namespace

static void
mirror_upper_half(std::vector<int64_t> &cellIndices, size_t xSize, size_t ySize)
{
  auto ySizeHalf = ySize / 2;
  for (size_t j = 0, index = 0; j < ySizeHalf; ++j)
  {
    for (size_t i = 0; i < xSize; ++i) { cellIndices[index++] = (ySize - j - 1) * xSize + i; }
  }
}

static void
mirror_lower_half(std::vector<int64_t> &cellIndices, size_t xSize, size_t ySize)
{
  auto ySizeHalf = ySize / 2;
  for (size_t j = 0, index = ySizeHalf * xSize; j < ySizeHalf; ++j)
  {
    for (size_t i = 0; i < xSize; ++i) { cellIndices[index++] = (ySizeHalf - j - 1) * xSize + i; }
  }
}

static int
generate_indices_reg2d(GridInfo &gridInfo, Parameter const &params)
{
  auto &cellIndices = gridInfo.cellIndices;
  auto gridSize = gridInfo.gridSize;
  int gridID = gridInfo.gridID;

  for (size_t i = 0; i < gridSize; ++i) { cellIndices[i] = i; }

  auto xSize = gridInqXsize(gridID);
  auto ySize = gridInqYsize(gridID);
  if (ySize <= 2) return -1;
  // if (ySize % 2 != 0) return -1;

  auto isSouthNorth = (params.data != LatOrientation::Undefined) ? (params.data == LatOrientation::Negative) : true;

  {
    std::vector<double> yVals(ySize);
    gridInqYvals(gridID, yVals.data());
    if (yVals[0] > yVals[ySize - 1]) { isSouthNorth = false; }
  }

  if (params.lat == LatOrientation::Positive)
  {
    if (isSouthNorth) { mirror_upper_half(cellIndices, xSize, ySize); }
    else { mirror_lower_half(cellIndices, xSize, ySize); }
  }
  else
  {
    if (isSouthNorth) { mirror_lower_half(cellIndices, xSize, ySize); }
    else { mirror_upper_half(cellIndices, xSize, ySize); }
  }

  return 0;
}

static int
generate_indices_unstruct(GridInfo &gridInfo, Parameter const &params)
{
  auto &cellIndices = gridInfo.cellIndices;
  auto gridSize = gridInfo.gridSize;

  auto gridID = (params.gridFile.size() > 0) ? cdo_define_grid(params.gridFile) : gridInfo.gridID;

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  for (size_t i = 0; i < gridSize; ++i) { cellIndices[i] = i; }

  std::vector<double> xVals(gridSize);
  std::vector<double> yVals(gridSize);
  gridInqXvals(gridID, xVals.data());
  gridInqYvals(gridID, yVals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xVals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yVals, "grid center lat");

  size_t numNeighbors = 1;

  std::vector<KnnData> knnDataList;
  knnDataList.reserve(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(numNeighbors);

  cdo::timer timer;

  GridPointsearch gps;
  grid_pointsearch_create_unstruct(gps, xVals, yVals, true);

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds (%zu points)", timer.elapsed(), gridSize);

  timer.reset();

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
  for (size_t i = 0; i < gridSize; ++i)
  {
    auto ompthID = cdo_omp_get_thread_num();

    if ((params.lat == LatOrientation::Positive && yVals[i] < 0) || (params.lat == LatOrientation::Negative && yVals[i] > 0))
    {
      auto &knnData = knnDataList[ompthID];
      grid_search_point_unstruct(gps, PointLonLat{ xVals[i], -yVals[i] }, knnData);
      cellIndices[i] = knnData.m_indices[0];
    }
  }

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", timer.elapsed());

  return 0;
}

static int
generate_indices(GridInfo &gridInfo, Parameter const &params)
{
  auto gridType = gridInfo.gridType;
  if (gridType == GRID_GAUSSIAN || gridType == GRID_LONLAT) { return generate_indices_reg2d(gridInfo, params); }
  else if (gridType == GRID_UNSTRUCTURED || gridType == GRID_CURVILINEAR) { return generate_indices_unstruct(gridInfo, params); }
  else { cdo_abort("Unsupported grid!"); }

  return 0;
}

LatOrientation
get_lat_orientation(std::string const &value)
{
  auto type = string_to_lower(value);
  if (type.starts_with("positive")) return LatOrientation::Positive;
  if (type.starts_with("negativ")) return LatOrientation::Negative;

  cdo_abort("parameter type=%s: invalid value!", value);
  return LatOrientation::Undefined;
}

static Parameter
get_parameter()
{
  Parameter params;

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
      if      (key == "lat")   { params.lat = get_lat_orientation(value); }
      else if (key == "data")  { params.data = get_lat_orientation(value); }
      else if (key == "grid")  { params.gridFile = parameter_to_word(value); }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

class Symmetrize : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Symmetrize",
    .operators = { { "symmetrize", SymmetrizeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Symmetrize>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numGrids{};

  VarList varList1;
  std::vector<bool> varIDs;
  std::vector<GridInfo> gridInfoList;

public:
  void
  init() override
  {
    auto params = get_parameter();

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numVars = varList1.numVars();
    varIDs.resize(numVars, false);

    numGrids = varList1.numGrids();
    gridInfoList.resize(numGrids);

    for (int index = 0; index < numGrids; ++index)
    {
      auto &gridInfo = gridInfoList[index];
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridType = gridInqType(gridID1);
      if (is_point_grid(gridID1))
      {
        auto gridSize = gridInqSize(gridID1);
        if (gridSize == 1) continue;

        gridInfo.cellIndices.resize(gridSize);
        gridInfo.gridSize = gridSize;
        gridInfo.gridType = gridType;
        gridInfo.gridID = gridID1;

        auto status = generate_indices(gridInfo, params);
        if (status != CDI_UNDEFID)
        {
          for (auto const &var : varList1.vars)
            if (gridID1 == var.gridID) { varIDs[var.ID] = true; }
          {
          }
        }
      }
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridType)); }
    }

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
            if (gridID1 == gridInfoList[index].gridID) break;
          if (index == numGrids) cdo_abort("Internal problem, grid not found!");

          field2.init(varList1.vars[varID]);
          window_cell(field1, field2, gridInfoList[index].cellIndices);

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
