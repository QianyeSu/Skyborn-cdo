/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selgridcell     selgridcell    Select grid cells by indices
*/

#include <algorithm>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "util_files.h"
#include "param_conversion.h"
#include "mpim_grid/grid_healpix.h"

int gengridcell(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices);

static int
gen_grid_healpix_index(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices)
{
  auto hpParams = cdo::get_healpix_params(gridID1);
  auto refinementLevel = hpParams.level();
  auto healpixOrder = (hpParams.order() == HpOrder::Ring) ? "ring" : "nested";

  auto gridID2 = gridCreate(GRID_HEALPIX, gridsize2);
  cdiDefKeyString(gridID2, CDI_GLOBAL, CDI_KEY_DIMNAME, "cell");
  cdiDefKeyString(gridID2, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "crs");
  const std::string gridmapName = "healpix";
  cdiDefKeyString(gridID2, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gridmapName.c_str());
  cdiDefAttTxt(gridID2, CDI_GLOBAL, "grid_mapping_name", (int) gridmapName.size(), gridmapName.c_str());
  cdiDefAttInt(gridID2, CDI_GLOBAL, "refinement_level", CDI_DATATYPE_INT32, 1, &refinementLevel);
  cdiDefAttTxt(gridID2, CDI_GLOBAL, "indexing_scheme", (int) std::strlen(healpixOrder), healpixOrder);

  if (refinementLevel < 14) cdiDefKeyInt(gridID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_INT32);

  gridDefIndices(gridID2, cellIndices.data());

  return gridID2;
}

static int
gen_grid_unstr_from_healpix(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices)
{
  Varray<double> xvals(gridsize2), yvals(gridsize2);

  auto hpParams = cdo::get_healpix_params(gridID1);

  for (size_t i = 0; i < gridsize2; ++i) { hp_index_to_lonlat(hpParams, cellIndices[i], &xvals[i], &yvals[i]); }

  auto gridID2 = gridCreate(GRID_UNSTRUCTURED, gridsize2);
  cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_UNITS, "radians");
  cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_UNITS, "radians");
  gridDefXvals(gridID2, xvals.data());
  gridDefYvals(gridID2, yvals.data());

  return gridID2;
}

static int
genindexgrid(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices)
{
  auto gridID0 = gridID1;
  auto gridtype1 = gridInqType(gridID1);

  if (is_healpix_grid(gridID1)) {}
  else if (gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN || gridtype1 == GRID_PROJECTION)
  {
    gridID1 = gridToCurvilinear(gridID1);
    gridtype1 = GRID_CURVILINEAR;
  }
  else if (gridtype1 == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID1))
  {
    auto reference = dereferenceGrid(gridID1);
    if (reference.isValid) gridID1 = reference.gridID;
    if (reference.notFound) cdo_warning("Reference to source grid not found!");
  }

  int gridID2 = -1;
  if (gridtype1 == GRID_UNSTRUCTURED || gridtype1 == GRID_CURVILINEAR)
    gridID2 = gengridcell(gridID1, gridsize2, cellIndices);
  else if (gridtype1 == GRID_GENERIC && gridInqYsize(gridID1) == 0)
    gridID2 = gengridcell(gridID1, gridsize2, cellIndices);
  else if (gridtype1 == GRID_HEALPIX)
    gridID2 = gen_grid_healpix_index(gridID1, gridsize2, cellIndices);
  else if (is_healpix_grid(gridID1))
    gridID2 = gen_grid_unstr_from_healpix(gridID1, gridsize2, cellIndices);

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static void
select_index(Field const &field1, Field &field2, long nind, std::vector<int64_t> const &cellIndex)
{
  auto func = [&](auto const &v1, auto &v2)
  {
    for (long i = 0; i < nind; ++i) v2[i] = v1[cellIndex[i]];
  };
  field_operation2(func, field1, field2);
}

class Selgridcell : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selgridcell",
    .operators = { { "selgridcell", 0, 0, "gridcell indices(1-N)", SelgridcellHelp },
                   { "delgridcell", 0, 0, "gridcell indices(1-N)", SelgridcellHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Selgridcell> registration = RegisterEntry<Selgridcell>(module);

private:
  struct sindex_t
  {
    int gridID1, gridID2;
  };
  int DELGRIDCELL{};
  int gridID1 = -1, gridID2{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  long numCells = 0;

  VarList varList1{};
  VarList varList2{};

  std::vector<int64_t> cellIndices{};
  std::vector<sindex_t> sindex{};
  std::vector<bool> processVars{};

public:
  void
  init() override
  {
    int numIndices = 0;
    std::vector<int> indices;

    DELGRIDCELL = module.get_id("delgridcell");

    operator_input_arg(cdo_operator_enter(0));

    auto operatorID = cdo_operator_id();

    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

    if (cdo_operator_argc() == 1)
    {
      bool maskfile = true, indexfile = false;
      auto const &filename = cdo_operator_argv(0);
      const auto *filenameCstr = filename.c_str();
      if (filename.starts_with("index="))
      {
        filenameCstr += 6;
        indexfile = true;
      }
      else if (filename.starts_with("mask="))
      {
        filenameCstr += 5;
        maskfile = true;
      }

      if (FileUtils::file_exists(filenameCstr))
      {
        if (indexfile)
        {
          std::vector<int> cdo_read_index(const char *indexFileName);
          indices = cdo_read_index(filenameCstr);
          numIndices = indices.size();
          if (numIndices == 0) cdo_abort("Index file %s generates no input!", cdo_operator_argv(0));
        }
        else if (maskfile)
        {
          std::vector<bool> cdo_read_mask(const char *maskFileName);
          auto mask = cdo_read_mask(filenameCstr);
          numIndices = 0;
          auto maskSize = mask.size();
          for (size_t i = 0; i < maskSize; ++i)
            if (mask[i]) { numIndices++; }
          if (numIndices == 0) cdo_abort("Mask is empty!");

          indices.resize(numIndices);
          numIndices = 0;
          for (size_t i = 0; i < maskSize; ++i)
            if (mask[i]) { indices[numIndices++] = i; }
          if (numIndices == 0) cdo_abort("Mask file %s generates no input!", cdo_operator_argv(0));
        }
      }
    }

    if (numIndices == 0)
    {
      indices = cdo_argv_to_intarr(cdo_get_oper_argv());
      numIndices = indices.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < numIndices; ++i) cdo_print("int %d = %d", i + 1, indices[i]);

      for (int i = 0; i < numIndices; ++i) indices[i] -= 1;
    }

    if (numIndices == 0) cdo_abort("Argument %s generates no input!", cdo_operator_argv(0));

    auto [indmin, indmax] = std::ranges::minmax_element(indices);
    if (*indmin < 0) cdo_abort("Index < 1 not allowed!");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    auto numVars = varList1.numVars();
    processVars = std::vector<bool>(numVars, false);

    auto numGrids = vlistNumGrids(vlistID1);
    sindex = std::vector<sindex_t>(numGrids);

    numCells = numIndices;
    if (operatorID == DELGRIDCELL)
    {
      auto gridsizeMax = varList1.gridsizeMax();
      numCells = gridsizeMax - numIndices;
      cellIndices.resize(gridsizeMax, 1);
      for (long i = 0; i < numIndices; ++i) cellIndices[indices[i]] = 0;
      long j = 0;
      for (size_t i = 0; i < gridsizeMax; ++i)
        if (cellIndices[i] == 1) cellIndices[j++] = i;
      if (j != numCells) cdo_abort("Internal error; number of cells differ");
    }
    else
    {
      cellIndices.resize(numIndices);
      for (int i = 0; i < numIndices; ++i) cellIndices[i] = indices[i];
    }

    if (numCells == 0) cdo_abort("Mask is empty!");

    for (int index = 0; index < numGrids; ++index)
    {
      gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      auto gridsize = gridInqSize(gridID1);
      if (gridsize == 1) continue;
      if ((size_t) *indmax >= gridsize)
      {
        cdo_warning("Max grid index is greater than grid size, skipped grid %d!", index + 1);
        continue;
      }

      gridID2 = genindexgrid(gridID1, numCells, cellIndices);
      if (gridID2 == -1)
      {
        cdo_warning("Unsupported grid type >%s<, skipped grid %d!", gridNamePtr(gridtype), index + 1);
        continue;
      }

      sindex[index].gridID1 = gridID1;
      sindex[index].gridID2 = gridID2;

      vlistChangeGridIndex(vlistID2, index, gridID2);

      for (int varID = 0; varID < numVars; ++varID)
        if (gridID1 == varList1.vars[varID].gridID) processVars[varID] = true;
    }

    {
      int varID = 0;
      for (; varID < numVars; ++varID)
        if (processVars[varID]) break;

      if (varID >= numVars) cdo_abort("No variables selected!");
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
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        cdo_def_field(streamID2, varID, levelID);

        if (processVars[varID])
        {
          field2.init(varList1.vars[varID]);
          select_index(field1, field2, numCells, cellIndices);

          if (field1.numMissVals) field2.numMissVals = field_num_mv(field2);

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

    vlistDestroy(vlistID2);
  }
};
