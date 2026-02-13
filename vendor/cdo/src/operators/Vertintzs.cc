/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    zs2zl           Model depth level to depth level interpolation
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "util_string.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"

static bool
is_depth_axis(int zaxisID)
{
  return (zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA);
}

static int
create_zaxis_depth(Varray<double> &depthLevels)
{
  int zaxisID = CDI_UNDEFID;
  auto &arg1 = cdo_operator_argv(0);
  if (cdo_operator_argc() == 1 && !std::isdigit(arg1[0]))
  {
    auto const &zfilename = arg1;
    auto fobj = c_fopen(zfilename, "r");
    if (fobj.get())
    {
      zaxisID = zaxis_from_file(fobj.get(), zfilename);
      if (zaxisID == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zfilename);
      auto numLevels = zaxisInqSize(zaxisID);
      depthLevels.resize(numLevels);
      zaxisInqLevels(zaxisID, depthLevels.data());
    }
    else if (arg1 == "default")
      depthLevels = { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 };
    else
      cdo_abort("Open failed on %s", zfilename);
  }
  else { depthLevels = cdo_argv_to_fltarr(cdo_get_oper_argv()); }

  if (zaxisID == CDI_UNDEFID)
  {
    zaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, depthLevels.size());
    zaxisDefLevels(zaxisID, depthLevels.data());
  }

  return zaxisID;
}

class Vertintzs : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertintzs",
    .operators = { { "zs2zl", 0, 0, "depth levels in meter" }, { "zs2zlx" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertintzs> registration = RegisterEntry<Vertintzs>();

private:
  int numVars{};
  std::vector<bool> processVars;
  std::vector<bool> interpVars;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;
  Field3DVector vardata1;
  Field3DVector vardata2;

  Varray2D<size_t> varnumMissVals;

  int depthID = -1;

  Field3D fullDepth;

  std::vector<int> vertIndexFull;
  Varray<double> depthLevels;

  int gridsize{};
  bool extrapolate = true;  // do not use missing values
                            //
  Field depthBottom;
  int numFullLevels{};

  Varray<size_t> pnumMissVals;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto zaxisID2 = create_zaxis_depth(depthLevels);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    gridsize = vlist_check_gridsize(vlistID1);

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);
    auto memType = varList1.vars[0].memType;

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      if (string_to_lower(varList1.vars[varID].name) == "depth_c") depthID = varID;
    }

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      if (-1 != depthID) { cdo_print("  %s -> %s", "zstar depth at cell center", varList1.vars[depthID].name); }
    }

    if (-1 == depthID) cdo_abort("depth_c not found!");

    auto zaxisIDfull = (-1 == depthID) ? -1 : varList1.vars[depthID].zaxisID;
    numFullLevels = (-1 == zaxisIDfull) ? 0 : zaxisInqSize(zaxisIDfull);

    auto numZaxes = vlistNumZaxis(vlistID1);
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisIDfull || (is_depth_axis(zaxisID) && nlevels == numFullLevels))
        vlistChangeZaxis(vlistID2, zaxisID, zaxisID2);
    }

    varList2 = VarList(vlistID2);
    varList_set_memtype(varList2, memType);

    if (!extrapolate) pnumMissVals.resize(depthLevels.size());

    vertIndexFull.resize(gridsize * depthLevels.size());

    processVars.resize(numVars);
    interpVars.resize(numVars);
    varnumMissVals.resize(numVars);
    vardata1.resize(numVars);
    vardata2.resize(numVars);

    auto maxLevels = std::max(numFullLevels, (int) depthLevels.size());

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];

      if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var);
      varnumMissVals[varID].resize(maxLevels, 0);

      interpVars[varID] = (var.zaxisID == zaxisIDfull || (is_depth_axis(var.zaxisID) && (var.nlevels == numFullLevels)));

      if (interpVars[varID]) { vardata2[varID].init(varList2.vars[varID]); }
      else if (is_depth_axis(var.zaxisID) && var.nlevels > 1)
      {

        cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1, var.name, var.nlevels);
      }
    }

    fullDepth.init(varList1.vars[depthID]);

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (interpVars[varID] && varList1.vars[varID].isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int varID = 0; varID < numVars; ++varID)
      {
        processVars[varID] = false;
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID) varnumMissVals[varID][levelID] = 0;
      }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, vardata1[varID], levelID, &varnumMissVals[varID][levelID]);
        processVars[varID] = true;
      }

      for (int varID = 0; varID < numVars; ++varID)
        if (interpVars[varID]) processVars[varID] = true;

      if (tsID == 0 || !varList1.vars[depthID].isConstant)
      {
        constexpr auto lreverse = false;
        field_copy(vardata1[depthID], fullDepth);
        gen_vert_index(vertIndexFull, depthLevels, fullDepth, gridsize, lreverse);
        if (!extrapolate)
        {
          depthBottom.init(varList1.vars[depthID]);
          field_copy(fullDepth, numFullLevels - 1, depthBottom);
          gen_vert_index_mv(vertIndexFull, depthLevels, gridsize, depthBottom, pnumMissVals, lreverse);
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        if (processVars[varID])
        {
          auto const &var = varList1.vars[varID];

          if (tsID > 0 && !interpVars[varID] && var.isConstant) continue;

          if (interpVars[varID])
          {
            if (var.nlevels != numFullLevels) cdo_abort("Number of depth level differ from full level (param=%s)!", var.name);

            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
            }

            vertical_interp_X(fullDepth, vardata1[varID], vardata2[varID], vertIndexFull, depthLevels, gridsize);

            if (!extrapolate) varray_copy(depthLevels.size(), pnumMissVals, varnumMissVals[varID]);
          }

          for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
          {
            cdo_def_field(streamID2, varID, levelID);
            auto varout = (interpVars[varID] ? vardata2[varID] : vardata1[varID]);
            cdo_write_field(streamID2, varout, levelID, varnumMissVals[varID][levelID]);
          }
        }
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
