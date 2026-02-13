/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "gridreference.h"
#include "util_files.h"

static void
change_grid(int vlistID1, int vlistID2, int gridID2)
{
  auto gridsize2 = gridInqSize(gridID2);
  int found = 0;
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    if (gridsize2 == gridInqSize(vlistGrid(vlistID1, index)))
    {
      vlistChangeGridIndex(vlistID2, index, gridID2);
      found++;
    }
  }
  if (!found) cdo_warning("No horizontal grid with %zu cells found!", gridsize2);
}

static void
grid_set_type(int vlistID1, int vlistID2, int gridType, std::string const &gridname, bool lregular, bool lregularnn,
              bool ldereference, bool withBounds, std::vector<int> &grid2_vgpm)
{
  auto needCorners = withBounds ? NeedCorners::Yes : NeedCorners::No;
  int gridID2;
  auto lrgrid = false;
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    auto gridType1 = gridInqType(gridID1);
    gridID2 = -1;

    if (gridType1 == GRID_GENERIC && gridInqSize(gridID1) == 1) continue;

    if (lregular || lregularnn)
    {
      if (gridType1 == GRID_GAUSSIAN_REDUCED) gridID2 = gridToRegular(gridID1);
    }
    else if (ldereference)
    {
      auto reference = dereferenceGrid(gridID1);
      if (reference.isValid) gridID2 = reference.gridID;
      if (reference.notFound) cdo_abort("Reference to horizontal grid not found!");
    }
    else
    {
      if (gridType == GRID_CURVILINEAR)
      {
        gridID2 = (gridType1 == GRID_CURVILINEAR) ? gridID1 : gridToCurvilinear(gridID1, needCorners);
      }
      else if (gridType == GRID_UNSTRUCTURED)
      {
        gridID2 = gridToUnstructured(gridID1, needCorners);
        if (gridType1 == GRID_GME)
        {
          auto grid2_nvgp = gridInqSize(gridID2);
          grid2_vgpm.resize(grid2_nvgp);
          gridInqMaskGME(gridID2, grid2_vgpm.data());
          gridCompress(gridID2);
        }
      }
      else if (gridType == GRID_LONLAT && gridType1 == GRID_CURVILINEAR)
      {
        gridID2 = gridCurvilinearToRegular(gridID1);
        if (gridID2 == -1) cdo_warning("Conversion of curvilinear grid to regular grid failed!");
      }
      else if (gridType == GRID_LONLAT && gridType1 == GRID_UNSTRUCTURED)
      {
        gridID2 = -1;
        cdo_warning("Conversion of unstructured grid to regular grid failed!");
      }
      else if (gridType == GRID_LONLAT && gridType1 == GRID_GENERIC)
      {
        gridID2 = -1;
        cdo_warning("Conversion of generic grid to regular grid failed!");
      }
      else if (gridType == GRID_LONLAT && gridType1 == GRID_PROJECTION)
      {
        gridID2 = gridProjectionToRegular(gridID1);
        if (gridID2 == -1) cdo_warning("Conversion of projection to regular grid failed!");
      }
      else if (gridType == GRID_LONLAT && gridType1 == GRID_LONLAT) { gridID2 = gridID1; }
      else if (gridType == GRID_PROJECTION)
      {
        if (gridInqType(gridID1) == GRID_PROJECTION)
          gridID2 = gridID1;
        else
        {
          auto projID = gridInqProj(gridID1);
          if (projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION) gridID2 = projID;
        }
        if (gridID2 == -1) cdo_warning("Projection not found!");
      }
      else
        cdo_abort("Unsupported grid name: %s", gridname);
    }

    if (gridID2 == -1)
    {
      if (!(lregular || lregularnn)) cdo_abort("Unsupported grid type!");
    }

    if (gridID2 != -1)
    {
      if (lregular || lregularnn) lrgrid = true;
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }
  }

  if ((lregular || lregularnn) && !lrgrid) cdo_warning("No reduced Gaussian grid found!");
}

static void
grid_set_cellarea(int vlistID1, int vlistID2, Varray<double> &gridcellArea)
{
  auto areasize = gridcellArea.size();
  int numGridsChanges = 0;
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    auto gridsize = gridInqSize(gridID1);
    if (gridsize == areasize)
    {
      auto gridID2 = gridDuplicate(gridID1);
      gridDefArea(gridID2, gridcellArea.data());
      vlistChangeGridIndex(vlistID2, index, gridID2);
      numGridsChanges++;
    }
  }
  if (!numGridsChanges) cdo_warning("No grid with %zu cells found!", areasize);
}

static void
grid_set_mask(int vlistID1, int vlistID2, Varray<double> const &gridmask)
{
  auto masksize = gridmask.size();
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    auto gridsize = gridInqSize(gridID1);
    if (gridsize == masksize)
    {
      std::vector<int> mask(masksize);
      for (size_t i = 0; i < masksize; ++i)
      {
        mask[i] = (gridmask[i] < 0 || gridmask[i] > 255) ? 0 : (int) std::lround(gridmask[i]);
      }
      auto gridID2 = gridDuplicate(gridID1);
      gridDefMask(gridID2, mask.data());
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }
  }
}

static void
grid_unset_mask(int vlistID1, int vlistID2)
{
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    auto gridID2 = gridDuplicate(gridID1);
    gridDefMask(gridID2, nullptr);
    vlistChangeGridIndex(vlistID2, index, gridID2);
  }
}

static void
grid_set_proj_params(int vlistID1, int vlistID2, std::string const &projparams)
{
  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    if (gridInqType(gridID1) == GRID_PROJECTION)
    {
      auto gridID2 = gridDuplicate(gridID1);
      cdiDefAttTxt(gridID2, CDI_GLOBAL, "proj_params", (int) projparams.size(), projparams.c_str());
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }
  }
}

static void
grid_read_cellarea(std::string const &areafile, Varray<double> &gridcellArea)
{
  auto searchName = false;
  std::string filename = areafile;
  std::string varname;

  if (!FileUtils::file_exists(areafile))
  {
    auto pos = filename.find_last_of(':');
    if (pos > 1 && pos < (filename.size() - 1))
    {
      varname = filename.substr(pos + 1);
      filename.resize(pos);
      searchName = true;
    }
  }

  auto streamID = stream_open_read_locked(filename.c_str());
  auto vlistID = streamInqVlist(streamID);

  VarList varList(vlistID);

  int svarID = 0;
  if (searchName)
  {
    auto numVars = varList.numVars();
    int varID = 0;
    for (; varID < numVars; ++varID)
    {
      if (varList.vars[varID].name == varname)
      {
        svarID = varID;
        break;
      }
    }
    if (varID == numVars) cdo_abort("Variable %s not found in %s\n", varname, filename);
  }

  auto numFields = streamInqTimestep(streamID, 0);
  for (int fieldID = 0; fieldID < numFields; ++fieldID)
  {
    int varID, levelID;
    streamInqField(streamID, &varID, &levelID);
    if (varID == svarID)
    {
      gridcellArea.resize(varList.vars[varID].gridsize);

      size_t numMissVals;
      streamReadField(streamID, gridcellArea.data(), &numMissVals);
      break;
    }
  }

  streamClose(streamID);
}

class Setgrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setgrid",
    // clang-format off
    .operators = { { "setgrid", 0, 0, "grid description file or name", SetgridHelp },
                   { "setgridtype", 0, 0, "gridtype", SetgridHelp },
                   { "setgridarea", 0, 0, "filename with areaweights", SetgridHelp },
                   { "setgridmask", 0, 0, "filename with gridmask", SetgridHelp },
                   { "unsetgridmask", SetgridHelp },
                   { "setgridnumber", 0, 0, "gridnumber and optionally grid position", SetgridHelp },
                   { "setgriduri", 0, 0, "reference URI of the horizontal grid", SetgridHelp },
                   { "usegridnumber", 0, 0, "use existing grid identified by gridnumber", SetgridHelp },
                   { "setprojparams", 0, 0, "proj library parameter (e.g.:+init=EPSG:3413)", SetgridHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Setgrid>();

private:
  int SETGRID{}, SETGRIDTYPE{}, SETGRIDAREA{}, SETGRIDMASK{}, UNSETGRIDMASK{}, SETGRIDNUMBER{}, SETGRIDURI{}, USEGRIDNUMBER{},
      SETPROJPARAMS{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool lregular = false;
  bool lregularnn = false;

  std::vector<int> grid2_vgpm;
  Varray<double> gridmask, gridcellArea;
  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    int number = 0, position = 0;
    std::string griduri;
    std::string projparams;
    std::string gridname;
    std::string gridfile;

    SETGRID = module.get_id("setgrid");
    SETGRIDTYPE = module.get_id("setgridtype");
    SETGRIDAREA = module.get_id("setgridarea");
    SETGRIDMASK = module.get_id("setgridmask");
    UNSETGRIDMASK = module.get_id("unsetgridmask");
    SETGRIDNUMBER = module.get_id("setgridnumber");
    SETGRIDURI = module.get_id("setgriduri");
    USEGRIDNUMBER = module.get_id("usegridnumber");
    SETPROJPARAMS = module.get_id("setprojparams");

    auto operatorID = cdo_operator_id();

    if (operatorID != UNSETGRIDMASK) operator_input_arg(cdo_operator_enter(operatorID));

    if (operatorID == SETGRID)
    {
      operator_check_argc(1);
      gridfile = cdo_operator_argv(0);
    }
    else if (operatorID == SETGRIDTYPE)
    {
      operator_check_argc(1);
      gridname = cdo_operator_argv(0);
    }
    else if (operatorID == SETGRIDAREA)
    {
      operator_check_argc(1);

      grid_read_cellarea(cdo_operator_argv(0), gridcellArea);

      if (Options::cdoVerbose)
      {
        auto areasize = gridcellArea.size();
        auto mmm = varray_min_max_mean(gridcellArea, areasize);
        cdo_print("gridcellAreas: %zu %#12.5g%#12.5g%#12.5g", areasize, mmm.min, mmm.mean, mmm.max);
      }
    }
    else if (operatorID == SETGRIDMASK)
    {
      operator_check_argc(1);
      auto maskfile = cdo_operator_argv(0).c_str();
      auto streamID = stream_open_read_locked(maskfile);
      auto vlistID = streamInqVlist(streamID);
      VarList varList(vlistID);

      (void) streamInqTimestep(streamID, 0);
      int varID, levelID;
      streamInqField(streamID, &varID, &levelID);
      auto missval = varList.vars[varID].missval;
      auto gridID = varList.vars[varID].gridID;
      auto masksize = gridInqSize(gridID);
      gridmask.resize(masksize);

      size_t numMissVals;
      streamReadField(streamID, gridmask.data(), &numMissVals);
      streamClose(streamID);

      for (size_t i = 0; i < masksize; ++i)
        if (fp_is_equal(gridmask[i], missval)) gridmask[i] = 0;
    }
    else if (operatorID == USEGRIDNUMBER)
    {
      operator_check_argc(1);
      number = parameter_to_int(cdo_operator_argv(0));
    }
    else if (operatorID == SETGRIDNUMBER)
    {
      if (cdo_operator_argc() >= 1 && cdo_operator_argc() <= 2)
      {
        number = parameter_to_int(cdo_operator_argv(0));
        if (cdo_operator_argc() == 2) position = parameter_to_int(cdo_operator_argv(1));
      }
      else { operator_check_argc(1); }
    }
    else if (operatorID == SETGRIDURI)
    {
      operator_check_argc(1);
      griduri = cdo_operator_argv(0);
    }
    else if (operatorID == SETPROJPARAMS)
    {
      operator_check_argc(1);
      projparams = cdo_operator_argv(0);
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if (operatorID == SETGRID)
    {
      auto gridID2 = cdo_define_grid(gridfile);
      change_grid(vlistID1, vlistID2, gridID2);
    }
    else if (operatorID == SETGRIDNUMBER || operatorID == SETGRIDURI || operatorID == USEGRIDNUMBER)
    {
      int gridID2x = -1;
      if (operatorID == SETGRIDNUMBER)
      {
        auto gridID1 = vlistGrid(vlistID1, 0);
        gridID2x = gridCreate(GRID_UNSTRUCTURED, gridInqSize(gridID1));
        cdiDefKeyInt(gridID2x, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, number);
        cdiDefKeyInt(gridID2x, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, position);
      }
      else if (operatorID == USEGRIDNUMBER)
      {
        if (number < 1 || number > vlistNumGrids(vlistID1))
          cdo_abort("Invalid grid number: %d (max = %d)!", number, vlistNumGrids(vlistID1));

        gridID2x = vlistGrid(vlistID1, number - 1);
      }
      else
      {
        auto gridID1 = vlistGrid(vlistID1, 0);
        gridID2x = gridDuplicate(gridID1);
        cdiDefKeyString(gridID2x, CDI_GLOBAL, CDI_KEY_REFERENCEURI, griduri.c_str());
      }

      change_grid(vlistID1, vlistID2, gridID2x);
    }
    else if (operatorID == SETGRIDTYPE)
    {
      int gridType = -1;
      bool ldereference = false;
      bool withBounds = true;

      // clang-format off
      if      (gridname == "curvilinear0")  { gridType = GRID_CURVILINEAR; withBounds = false; }
      else if (gridname == "curvilinear")   { gridType = GRID_CURVILINEAR; withBounds = true; }
      else if (gridname == "cell")          { gridType = GRID_UNSTRUCTURED; }
      else if (gridname == "unstructured0") { gridType = GRID_UNSTRUCTURED; withBounds = false; }
      else if (gridname == "unstructured")  { gridType = GRID_UNSTRUCTURED; withBounds = true; }
      else if (gridname == "generic")       { gridType = GRID_GENERIC; }
      else if (gridname == "dereference")   { ldereference = true; }
      else if (gridname == "lonlat")        { gridType = GRID_LONLAT; }
      else if (gridname == "gaussian")      { gridType = GRID_GAUSSIAN; }
      else if (gridname == "regularnn")     { gridType = GRID_GAUSSIAN; lregularnn = true; }
      else if (gridname == "regular")       { gridType = GRID_GAUSSIAN; lregular = true; }
      else if (gridname == "projection")    { gridType = GRID_PROJECTION; }
      else cdo_abort("Unsupported grid name: %s", gridname);
      // clang-format on

      grid_set_type(vlistID1, vlistID2, gridType, gridname, lregular, lregularnn, ldereference, withBounds, grid2_vgpm);
    }
    else if (operatorID == SETGRIDAREA) { grid_set_cellarea(vlistID1, vlistID2, gridcellArea); }
    else if (operatorID == SETGRIDMASK) { grid_set_mask(vlistID1, vlistID2, gridmask); }
    else if (operatorID == UNSETGRIDMASK) { grid_unset_mask(vlistID1, vlistID2); }
    else if (operatorID == SETPROJPARAMS) { grid_set_proj_params(vlistID1, vlistID2, projparams); }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
    // vlistPrint(vlistID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    if (lregular || lregularnn)
      for (auto &var : varList1.vars) var.memType = MemType::Double;
    if (lregular || lregularnn)
      for (auto &var : varList2.vars) var.memType = MemType::Double;
  }

  void
  run() override
  {
    Field field;

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
        cdo_def_field(streamID2, varID, levelID);

        auto const &var1 = varList1.vars[varID];
        auto const &var2 = varList2.vars[varID];
        field.init((lregular || lregularnn) ? var2 : var1);
        cdo_read_field(streamID1, field);
        auto numMissVals = field.numMissVals;

        auto gridID1 = var1.gridID;
        if (lregular || lregularnn)
        {
          if (gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED)
          {
            auto missval = var1.missval;
            int lnearst = lregularnn ? 1 : 0;
            field2regular(gridID1, var2.gridID, missval, field.vec_d, numMissVals, lnearst);
          }
        }
        else if (gridInqType(gridID1) == GRID_GME)
        {
          auto n = var1.gridsize;
          auto func = [&](auto &v)
          {
            for (size_t i = 0, j = 0; i < n; ++i)
              if (grid2_vgpm[i]) v[j++] = v[i];
          };
          field_operation(func, field);
        }

        cdo_write_field(streamID2, field);
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
