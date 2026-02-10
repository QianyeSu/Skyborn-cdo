/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "util_files.h"
#include "cdo_options.h"
#include "cdi_lockedIO.h"
#include "cdo_omp.h"
#include "pmlist.h"

static int globalGridType = CDI_UNDEFID;

namespace
{
struct GridInfo2
{
  size_t nx{ 0 };
  int globalIndicesID{ -1 };
  bool needed{ false };
  bool isReg2D{ false };
};

struct CollgridInfo
{
  std::vector<std::vector<long>> cellIndices;
  VarList varList;
  CdoStreamID streamID;
  size_t nx{ 0 }, ny{ 0 };
  size_t offset{ 0 };
};

struct xyinfoType
{
  double x{ 0.0 };
  double y{ 0.0 };
  int id{ -1 };
};
}  // namespace

static bool
cmpxy_lt(xyinfoType const &a, xyinfoType const &b)
{
  return (a.y < b.y || (std::fabs(a.y - b.y) <= 0 && a.x < b.x));
}

static bool
cmpxy_gt(xyinfoType const &a, xyinfoType const &b)
{
  return (a.y > b.y || (std::fabs(a.y - b.y) <= 0 && a.x < b.x));
}

static int
gen_coll_grid(int numGrids, int numFiles, GridInfo2 const &gridInfo2, std::vector<CollgridInfo> &collgridInfo, int gindex,
              long nxblocks)
{
  auto isReg2D = gridInfo2.isReg2D;
  auto isSouthNorth = true;
  auto isRegular = false;
  auto isCurvilinear = false;

  long nx = (nxblocks != -1) ? nxblocks : -1;

  auto gridID = vlistGrid(collgridInfo[0].varList.vlistID, gindex);
  auto gridtype0 = (globalGridType != CDI_UNDEFID) ? globalGridType : gridInqType(gridID);
  if (numGrids > 1 && gridtype0 == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0) return -1;

  auto isUnstructured = (gridtype0 == GRID_UNSTRUCTURED);
  auto nv = isUnstructured ? gridInqNvertex(gridID) : 0;
  auto withCenter = (globalGridType == CDI_UNDEFID && gridHasCoordinates(gridID));
  auto withBounds = (isUnstructured && globalGridType == CDI_UNDEFID && gridHasBounds(gridID));

  std::vector<xyinfoType> xyinfoList(numFiles);
  std::vector<long> xsize(numFiles), ysize(numFiles);
  Varray2D<double> xvals(numFiles), yvals(numFiles);
  Varray2D<double> xbounds(numFiles), ybounds(numFiles);

  for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
  {
    gridID = vlistGrid(collgridInfo[fileIdx].varList.vlistID, gindex);
    auto gridType = (globalGridType != CDI_UNDEFID) ? globalGridType : gridInqType(gridID);
    if (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_PROJECTION)
      isRegular = true;
    else if (gridType == GRID_CURVILINEAR)
      isCurvilinear = true;
    else if (gridType == GRID_UNSTRUCTURED)
      isUnstructured = true;
    else if (gridType == GRID_GENERIC /*&& gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0*/)
      isRegular = withCenter;
    else
      cdo_abort("Unsupported grid type: %s!", gridNamePtr(gridType));

    xsize[fileIdx] = isUnstructured ? gridInqSize(gridID) : gridInqXsize(gridID);
    ysize[fileIdx] = isUnstructured ? 1 : gridInqYsize(gridID);
    if (xsize[fileIdx] == 0) xsize[fileIdx] = 1;
    if (ysize[fileIdx] == 0) ysize[fileIdx] = 1;

    if (isRegular)
    {
      xvals[fileIdx].resize(xsize[fileIdx]);
      yvals[fileIdx].resize(ysize[fileIdx]);
    }
    else if (isCurvilinear || isUnstructured)
    {
      if (withCenter) xvals[fileIdx].resize(xsize[fileIdx] * ysize[fileIdx]);
      if (withCenter) yvals[fileIdx].resize(xsize[fileIdx] * ysize[fileIdx]);
      if (withBounds) xbounds[fileIdx].resize(nv * xsize[fileIdx] * ysize[fileIdx]);
      if (withBounds) ybounds[fileIdx].resize(nv * xsize[fileIdx] * ysize[fileIdx]);
    }

    if (isRegular || isCurvilinear || isUnstructured)
    {
      if (withCenter) gridInqXvals(gridID, xvals[fileIdx].data());
      if (withCenter) gridInqYvals(gridID, yvals[fileIdx].data());
      if (withBounds) gridInqXbounds(gridID, xbounds[fileIdx].data());
      if (withBounds) gridInqYbounds(gridID, ybounds[fileIdx].data());
    }
    // printf("fileIdx %d, gridID %d\n", fileIdx, gridID);

    xyinfoList[fileIdx].id = fileIdx;
    if (isRegular)
    {
      xyinfoList[fileIdx].x = xvals[fileIdx][0];
      xyinfoList[fileIdx].y = yvals[fileIdx][0];
      if (ysize[fileIdx] > 1 && yvals[fileIdx][0] > yvals[fileIdx][ysize[fileIdx] - 1]) isSouthNorth = false;
    }
  }

  if (isRegular)
  {
    if (Options::cdoVerbose)
      for (auto const &xyinfo : xyinfoList) printf("1 %d %g %g \n", xyinfo.id, xyinfo.x, xyinfo.y);

    std::ranges::sort(xyinfoList, {}, &xyinfoType::x);

    if (Options::cdoVerbose)
      for (auto const &xyinfo : xyinfoList) printf("2 %d %g %g \n", xyinfo.id, xyinfo.x, xyinfo.y);

    std::ranges::sort(xyinfoList, isSouthNorth ? cmpxy_lt : cmpxy_gt);

    if (Options::cdoVerbose)
      for (auto const &xyinfo : xyinfoList) printf("3 %d %g %g \n", xyinfo.id, xyinfo.x, xyinfo.y);

    if (nx <= 0)
    {
      nx = 1;
      for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      {
        if (fp_is_equal(xyinfoList[0].y, xyinfoList[fileIdx].y))
          nx++;
        else
          break;
      }
    }
  }
  else
  {
    if (nx <= 0) nx = numFiles;
  }

  long ny = numFiles / nx;
  if (nx * ny != numFiles) cdo_abort("Number of input files (%ld) and number of blocks (%ldx%ld) differ!", numFiles, nx, ny);

  long xsize2 = 0;
  for (long i = 0; i < nx; ++i) xsize2 += xsize[xyinfoList[i].id];
  long ysize2 = 0;
  for (long j = 0; j < ny; ++j) ysize2 += ysize[xyinfoList[j * nx].id];
  if (Options::cdoVerbose) cdo_print("xsize2 %ld  ysize2 %ld", xsize2, ysize2);

  {  // verify size of data
    auto xs = xsize[xyinfoList[0].id];
    for (long j = 1; j < ny; ++j)
      if (xsize[xyinfoList[j * nx].id] != xs)
        cdo_abort("xsize=%ld differ from first file (xsize=%ld)!", xsize[xyinfoList[j * nx].id], xs);
    auto ys = ysize[xyinfoList[0].id];
    for (long i = 1; i < nx; ++i)
      if (ysize[xyinfoList[i].id] != ys) cdo_abort("ysize=%ld differ from first file (ysize=%ld)!", ysize[xyinfoList[i].id], ys);
  }

  Varray<double> xvals2, yvals2;
  Varray<double> xbounds2, ybounds2;
  if (isRegular)
  {
    xvals2.resize(xsize2);
    yvals2.resize(ysize2);
  }
  else if (isCurvilinear || isUnstructured)
  {
    if (withCenter) xvals2.resize(xsize2 * ysize2);
    if (withCenter) yvals2.resize(xsize2 * ysize2);
    if (withBounds) xbounds2.resize(nv * xsize2 * ysize2);
    if (withBounds) ybounds2.resize(nv * xsize2 * ysize2);
  }

  std::vector<long> xoff(nx + 1), yoff(ny + 1);

  xoff[0] = 0;
  for (long i = 0; i < nx; ++i)
  {
    auto idx = xyinfoList[i].id;
    if (isRegular) array_copy(xsize[idx], xvals[idx].data(), &xvals2[xoff[i]]);
    xoff[i + 1] = xoff[i] + xsize[idx];
  }

  yoff[0] = 0;
  for (long j = 0; j < ny; ++j)
  {
    auto idx = xyinfoList[j * nx].id;
    if (isRegular) array_copy(ysize[idx], yvals[idx].data(), &yvals2[yoff[j]]);
    yoff[j + 1] = yoff[j] + ysize[idx];
  }

  for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
  {
    auto idx = xyinfoList[fileIdx].id;
    long iy = fileIdx / nx;
    long ix = fileIdx - iy * nx;

    long offset = yoff[iy] * xsize2 + xoff[ix];
    // printf("fileIdx %d %d, iy %d, ix %d, offset %d\n", fileIdx, xyinfo[fileIdx].id, iy, ix, offset);

    if (isReg2D)
    {
      collgridInfo[idx].nx = xsize[idx];
      collgridInfo[idx].ny = ysize[idx];
      collgridInfo[idx].offset = offset;
    }

    long ij = 0;
    for (long j = 0; j < ysize[idx]; ++j)
      for (long i = 0; i < xsize[idx]; ++i)
      {
        if (isCurvilinear || isUnstructured)
        {
          if (withCenter) { xvals2[offset + j * xsize2 + i] = xvals[idx][ij]; }
          if (withCenter) { yvals2[offset + j * xsize2 + i] = yvals[idx][ij]; }
          if (withBounds)
          {
            for (long k = 0; k < nv; ++k) xbounds2[(offset + j * xsize2 + i) * nv + k] = xbounds[idx][ij * nv + k];
          }
          if (withBounds)
          {
            for (long k = 0; k < nv; ++k) ybounds2[(offset + j * xsize2 + i) * nv + k] = ybounds[idx][ij * nv + k];
          }
        }
        if (!isReg2D) collgridInfo[idx].cellIndices[gindex][ij++] = offset + j * xsize2 + i;
      }
  }

  auto gridID2 = gridCreate(gridtype0, xsize2 * ysize2);
  if (!isUnstructured)
  {
    gridDefXsize(gridID2, xsize2);
    gridDefYsize(gridID2, ysize2);
  }
  else if (nv > 0) { gridDefNvertex(gridID2, nv); }

  if (isRegular || isCurvilinear || isUnstructured)
  {
    if (withCenter) gridDefXvals(gridID2, xvals2.data());
    if (withCenter) gridDefYvals(gridID2, yvals2.data());
    if (withBounds) gridDefXbounds(gridID2, xbounds2.data());
    if (withBounds) gridDefYbounds(gridID2, ybounds2.data());
  }

  gridID = vlistGrid(collgridInfo[0].varList.vlistID, gindex);

  grid_copy_names(gridID, gridID2);

  if (gridtype0 == GRID_PROJECTION) grid_copy_mapping(gridID, gridID2);

  return gridID2;
}

static void
collect_cells_reg2d(Field const &field1, Field &field2, CollgridInfo const &collgridInfo, size_t nlon)
{
  auto nx = collgridInfo.nx;
  auto ny = collgridInfo.ny;

  for (size_t j = 0; j < ny; ++j)
  {
    auto offset1 = j * nx;
    auto offset2 = collgridInfo.offset + j * nlon;

    auto func = [&](auto const &v1, auto &v2)
    {
      for (size_t i = 0; i < nx; ++i) { v2[offset2 + i] = v1[offset1 + i]; }
    };
    field_operation2(func, field1, field2);
  }
}

static void
collect_cells(Field const &field1, Field &field2, std::vector<long> const &cellIndex)
{
  auto func = [&](auto const &v1, auto &v2, auto n)
  {
    for (size_t i = 0; i < n; ++i) { v2[cellIndex[i]] = v1[i]; }
  };
  field_operation2(func, field1, field2, field1.size);
}

static std::vector<int>
get_var_gridindex(VarList const &varList)
{
  auto numVars = varList.numVars();
  auto numGrids = varList.numGrids();

  std::vector<int> varGridIndex(numVars, 0);
  for (auto const &var : varList.vars)
  {
    for (int index = 0; index < numGrids; ++index)
    {
      if (var.gridID == vlistGrid(varList.vlistID, index))
      {
        varGridIndex[var.ID] = index;
        break;
      }
    }
  }

  return varGridIndex;
}

static std::vector<GridInfo2>
get_gridinfo(VarList const &varList, std::vector<int> const &varGridIndex, std::vector<bool> &selectedVars)
{
  std::vector<GridInfo2> gridInfo(varList.numGrids());

  int globalCellIndicesID = -1;
  int globalVertIndicesID = -1;
  int globalEdgeIndicesID = -1;
  for (auto const &var : varList.vars)
  {
    // clang-format off
      if      (var.name == "global_cell_indices") globalCellIndicesID = var.ID;
      else if (var.name == "global_vert_indices") globalVertIndicesID = var.ID;
      else if (var.name == "global_edge_indices") globalEdgeIndicesID = var.ID;
    // clang-format on
  }
  if (globalCellIndicesID != -1) selectedVars[globalCellIndicesID] = false;
  if (globalVertIndicesID != -1) selectedVars[globalVertIndicesID] = false;
  if (globalEdgeIndicesID != -1) selectedVars[globalEdgeIndicesID] = false;

  if (globalCellIndicesID != -1) gridInfo[varGridIndex[globalCellIndicesID]].globalIndicesID = globalCellIndicesID;
  if (globalVertIndicesID != -1) gridInfo[varGridIndex[globalVertIndicesID]].globalIndicesID = globalVertIndicesID;
  if (globalEdgeIndicesID != -1) gridInfo[varGridIndex[globalEdgeIndicesID]].globalIndicesID = globalEdgeIndicesID;

  for (auto const &var : varList.vars)
    if (selectedVars[var.ID]) gridInfo[varGridIndex[var.ID]].needed = true;

  for (auto const &grid : gridInfo)
  {
    if (grid.needed && grid.globalIndicesID != CDI_UNDEFID)
    {
      if (grid.globalIndicesID == globalCellIndicesID) cdo_print("Using global_cell_indices array for indexing");
      if (grid.globalIndicesID == globalVertIndicesID) cdo_print("Using global_vert_indices array for indexing");
      if (grid.globalIndicesID == globalEdgeIndicesID) cdo_print("Using global_edge_indices array for indexing");
    }
  }

  for (int index = 0; index < varList.numGrids(); ++index)
  {
    auto gridID = vlistGrid(varList.vlistID, index);
    auto gridType = gridInqType(gridID);
    gridInfo[index].isReg2D = (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_PROJECTION);
  }

  return gridInfo;
}

static std::vector<bool>
get_selected_vars(VarList const &varList1, std::vector<std::string> const &nameList)
{
  auto numVars = varList1.numVars();
  std::vector<bool> selectedVars(numVars, false);

  int numNames = (int) nameList.size();
  if (numNames == 0)
  {
    for (int varID = 0; varID < numVars; ++varID) selectedVars[varID] = true;
  }
  else
  {
    if (Options::cdoVerbose)
      for (int i = 0; i < numNames; ++i) cdo_print("name %d = %s", i + 1, nameList[i]);

    std::vector<bool> selfound(numNames, false);

    for (int varID = 0; varID < numVars; ++varID)
    {
      for (int i = 0; i < numNames; ++i)
      {
        if (nameList[i] == varList1.vars[varID].name)
        {
          selfound[i] = true;
          selectedVars[varID] = true;
        }
      }
    }

    int notFound = 0;
    for (int i = 0; i < numNames; ++i)
    {
      if (selfound[i] == false)
      {
        notFound++;
        cdo_warning("Variable name %s not found!", nameList[i]);
      }
    }
    if (notFound) cdo_abort("Could not find all requested variables: (%d/%d)", numNames - notFound, numNames);
  }

  return selectedVars;
}

static std::vector<bool>
get_selected_levels(VarList const &varList1, std::vector<bool> const &selectedVars, std::vector<int> const &levelIndices)
{
  auto numVars = varList1.numVars();
  int maxLevels = 0;
  for (auto const &var : varList1.vars) { maxLevels = std::max(maxLevels, var.nlevels); }
  std::vector<bool> selectedLevels(maxLevels, false);

  int numLevelIndices = (int) levelIndices.size();
  if (numLevelIndices == 0)
  {
    for (int i = 0; i < maxLevels; ++i) selectedLevels[i] = true;
  }
  else
  {
    if (Options::cdoVerbose)
      for (int i = 0; i < numLevelIndices; ++i) cdo_print("levidx %d = %d", i + 1, levelIndices[i]);

    std::vector<bool> selfound(numLevelIndices, false);

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (selectedVars[varID])
      {
        for (int i = 0; i < numLevelIndices; ++i)
        {
          if (levelIndices[i] > 0 && levelIndices[i] < varList1.vars[varID].nlevels)
          {
            selfound[i] = true;
            selectedLevels[levelIndices[i] - 1] = true;
          }
        }
      }
    }

    int notFound = 0;
    for (int i = 0; i < numLevelIndices; ++i)
    {
      if (selfound[i] == false)
      {
        notFound++;
        cdo_warning("Level index %d not found!", levelIndices[i]);
      }
    }
    if (notFound) cdo_abort("Could not find all requested level indices: (%d/%d)", numLevelIndices - notFound, numLevelIndices);
  }

  return selectedLevels;
}

static void
select_vars(VarList const &varList1, std::vector<bool> const &selectedVars, std::vector<bool> const &selectedLevels)
{
  auto numVars = varList1.numVars();
  int numVarsFound = 0;

  for (int varID = 0; varID < numVars; ++varID)
  {
    if (selectedVars[varID])
    {
      numVarsFound++;
      auto numLevels = varList1.vars[varID].nlevels;
      for (int levelID = 0; levelID < numLevels; levelID++)
      {
        if (selectedLevels[levelID] || numLevels == 1) { vlistDefFlag(varList1.vlistID, varID, levelID, true); }
      }
    }
  }

  if (numVarsFound == 0) cdo_abort("No variables selected!");
}

namespace
{
struct Parameter
{
  std::vector<std::string> nameList;
  std::vector<int> levelIndices;
  int nx{ 0 };
  int gridType{ CDI_UNDEFID };
};
}  // namespace

static int
get_gridType(std::string_view gridTypeStr)
{
  if (gridTypeStr == "unstructured") return GRID_UNSTRUCTURED;
  cdo_abort("gridtype=%d unsupported!", gridTypeStr);
  return -1;
}

static Parameter
get_parameter(void)
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
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &values = kv.values;
      auto const &value = kv.values[0];
      int numValues = kv.nvalues;
      if (numValues == 1 && value.empty()) numValues = 0;

      // clang-format off
      if      (key == "nx")       params.nx = parameter_to_int(value);
      else if (key == "gridtype") params.gridType = get_gridType(parameter_to_word(value));
      else if (key == "name")
      {
        params.nameList.resize(numValues);
        for (int i = 0; i < numValues; ++i) { params.nameList[i] = parameter_to_word(values[i]); }
      }
      else if (key == "levidx")
      {
        params.levelIndices = cdo_argv_to_intarr(values);;
      }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

class Collgrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Collgrid",
    .operators = { { "collgrid", CollgridHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Collgrid> registration = RegisterEntry<Collgrid>(module);

  int nxblocks{ -1 };

  CdoStreamID streamID1;
  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };
  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int numFiles{};

  std::vector<GridInfo2> gridInfo;
  std::vector<CollgridInfo> collgridInfo;
  std::vector<bool> collectVars2;
  std::vector<int> gridID2s;
  std::vector<int> varGridIndex;

  VarList varList2;
  std::vector<size_t> targetGridsize;

public:
  void
  init() override
  {
    numFiles = cdo_stream_cnt() - 1;
    std::string ofilename = cdo_get_stream_name(numFiles);

    if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
      cdo_abort("Outputfile %s already exists!", ofilename);

    if (Options::cdoVerbose) cdo_print("Number of patches: %d", numFiles);
    collgridInfo.resize(numFiles);

    cdo::set_numfiles(numFiles + 8);

    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      auto streamID = cdo_open_read(fileIdx);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      collgridInfo[fileIdx].streamID = streamID;
      collgridInfo[fileIdx].varList = VarList(vlistID);
    }

    auto const &varList1 = collgridInfo[0].varList;
    vlistID1 = varList1.vlistID;
    vlistClearFlag(vlistID1);

    // check that the contents is always the same
    for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      varList_compare(varList1, collgridInfo[fileIdx].varList, CmpVarList::Name | CmpVarList::NumLevels);

    std::vector<std::string> nameList;
    std::vector<int> levelIndices;
    auto numArgs = cdo_operator_argc();
    if (numArgs > 0)
    {
      auto const &argList = cdo_get_oper_argv();
      if (numArgs == 1 && std::isdigit((int) argList[0][0])) { nxblocks = parameter_to_int(argList[0]); }
      else
      {
        auto params = get_parameter();
        nxblocks = params.nx;
        globalGridType = params.gridType;
        nameList = params.nameList;
        levelIndices = params.levelIndices;
      }
    }

    auto selectedVars = get_selected_vars(varList1, nameList);
    auto selectedLevels = get_selected_levels(varList1, selectedVars, levelIndices);
    varGridIndex = get_var_gridindex(varList1);
    gridInfo = get_gridinfo(varList1, varGridIndex, selectedVars);
    select_vars(varList1, selectedVars, selectedLevels);

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids1 = varList1.numGrids();
    auto numGrids2 = vlistNumGrids(vlistID2);
    targetGridsize.resize(numGrids1, 0);

    // allocate cellIndices
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      auto &collGrid = collgridInfo[fileIdx];
      collGrid.cellIndices.resize(numGrids1);
      for (int gindex = 0; gindex < numGrids1; ++gindex)
      {
        if (gridInfo[gindex].needed && gridInfo[gindex].isReg2D == false)
        {
          auto patchSize = gridInqSize(vlistGrid(collGrid.varList.vlistID, gindex));
          collGrid.cellIndices[gindex].resize(patchSize);
        }
      }
    }

    gridID2s.resize(numGrids2);

    for (int i2 = 0; i2 < numGrids2; ++i2)
    {
      int i1;
      for (i1 = 0; i1 < numGrids1; ++i1)
        if (vlistGrid(vlistID1, i1) == vlistGrid(vlistID2, i2)) break;

      gridID2s[i2] = gen_coll_grid(numGrids2, numFiles, gridInfo[i1], collgridInfo, i1, nxblocks);
      targetGridsize[i1] = gridInqSize(gridID2s[i2]);
      if (Options::cdoVerbose) cdo_print("Target grid%d size: %zu", i2 + 1, targetGridsize[i1]);
      if (gridInfo[i1].isReg2D) { gridInfo[i1].nx = gridInqXsize(gridID2s[i2]); }
    }

    for (int i = 0; i < numGrids2; ++i)
    {
      if (gridID2s[i] != -1) vlistChangeGridIndex(vlistID2, i, gridID2s[i]);
    }

    varList2 = VarList(vlistID2);
    auto numVars2 = varList2.numVars();

    collectVars2 = std::vector<bool>(numVars2, false);
    for (int varID = 0; varID < numVars2; ++varID)
    {
      auto gridID = varList2.vars[varID].gridID;
      for (int i = 0; i < numGrids2; ++i)
      {
        if (gridID2s[i] != -1 && gridID == vlistGrid(vlistID2, i))
        {
          collectVars2[varID] = true;
          break;
        }
      }
    }

    streamID2 = cdo_open_write(numFiles);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    std::vector<Field> field1vec(Threading::ompNumMaxThreads);
    Field field2;

    int numFields0 = 0;
    int tsID = 0;
    do {
      numFields0 = cdo_stream_inq_timestep(collgridInfo[0].streamID, tsID);
      for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      {
        auto numFields = cdo_stream_inq_timestep(collgridInfo[fileIdx].streamID, tsID);
        if (numFields != numFields0)
          cdo_abort("Number of fields at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(0),
                    cdo_get_stream_name(fileIdx));
      }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      if (numFields0 > 0) cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields0; ++fieldID)
      {
        int varID = 0, levelID = 0;
        for (int fileIdx = numFiles - 1; fileIdx >= 0; fileIdx--)
        {
          std::tie(varID, levelID) = cdo_inq_field(collgridInfo[fileIdx].streamID);
        }

        auto gindex = varGridIndex[varID];
        if (gridInfo[gindex].needed && gridInfo[gindex].isReg2D == false && gridInfo[gindex].globalIndicesID == varID)
        {
          std::vector<Varray<double>> cellIndicesMem(Threading::ompNumMaxThreads);
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
          {
            auto ompthID = cdo_omp_get_thread_num();
            auto &cellIndices = cellIndicesMem[ompthID];
            auto &collgrid = collgridInfo[fileIdx];
            auto patchSize = collgrid.varList.vars[varID].gridsize;
            if (cellIndices.size() < patchSize) cellIndices.resize(patchSize);
            size_t numMissVals;
            cdo_read_field(collgrid.streamID, cellIndices.data(), &numMissVals);
            for (size_t i = 0; i < patchSize; ++i)
            {
              auto index = std::lround(cellIndices[i]);
              if (index < 1 || index > (long) targetGridsize[gindex])
                cdo_abort("Global cell index out of range (%ld/%zu)", index, targetGridsize[gindex]);
              collgrid.cellIndices[gindex][i] = index - 1;
            }
          }
        }

        if (vlistInqFlag(vlistID1, varID, levelID) == true)
        {
          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
          // if (Options::cdoVerbose && tsID == 0) printf("varID %d %d levelID %d %d\n", varID, varID2, levelID, levelID2);

          field2.init(varList2.vars[varID2]);
          field_fill(field2, field2.missval);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
          {
            auto &collgrid = collgridInfo[fileIdx];
            auto ompthID = cdo_omp_get_thread_num();
            auto &field1 = field1vec[ompthID];
            field1.init(collgrid.varList.vars[varID]);
            cdo_read_field(collgrid.streamID, field1);

            if (collectVars2[varID2])
            {
              if (gridInfo[gindex].isReg2D)
                collect_cells_reg2d(field1, field2, collgrid, gridInfo[0].nx);
              else
                collect_cells(field1, field2, collgrid.cellIndices[gindex]);
            }
          }

          cdo_def_field(streamID2, varID2, levelID2);

          if (collectVars2[varID2])
          {
            field_num_mv(field2);
            cdo_write_field(streamID2, field2);
          }
          else { cdo_write_field(streamID2, field1vec[0]); }
        }
      }

      tsID++;
    } while (numFields0 > 0);
  }

  void
  close() override
  {
    for (auto const &collgrid : collgridInfo) cdo_stream_close(collgrid.streamID);

    cdo_stream_close(streamID2);
  }
};
