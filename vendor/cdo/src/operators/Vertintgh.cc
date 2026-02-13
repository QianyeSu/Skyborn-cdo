/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     gh2hl           Model geometric height level to height level interpolation
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "util_string.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "vertint_util.h"

static bool
is_height_axis(int zaxisID)
{
  auto isHeight = false;
  if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
  {
    auto units = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
    auto stdname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_STDNAME);
    if (stdname == "height" && units.empty()) isHeight = true;
  }
  return isHeight;
}

static int
create_zaxis_height(Varray<double> &heightLevels)
{
  int zaxisID = CDI_UNDEFID;
  auto const &arg1 = cdo_operator_argv(0);
  if (cdo_operator_argc() == 1 && !std::isdigit(arg1[0]))
  {
    auto const &zfilename = arg1;
    auto fobj = c_fopen(zfilename, "r");
    if (fobj.get())
    {
      zaxisID = zaxis_from_file(fobj.get(), zfilename);
      if (zaxisID == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zfilename);
      auto numLevels = zaxisInqSize(zaxisID);
      heightLevels.resize(numLevels);
      zaxisInqLevels(zaxisID, heightLevels.data());
    }
    else if (arg1 == "default")
      heightLevels = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
    else
      cdo_abort("Open failed on %s", zfilename);
  }
  else { heightLevels = cdo_argv_to_fltarr(cdo_get_oper_argv()); }

  if (zaxisID == CDI_UNDEFID)
  {
    zaxisID = zaxisCreate(ZAXIS_HEIGHT, heightLevels.size());
    zaxisDefLevels(zaxisID, heightLevels.data());
  }

  return zaxisID;
}

class Vertintgh : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertintgh",
    // clang-format off
    .operators = { { "gh2hl", 0, 0, "height levels in meter", VertintghHelp },
                   { "gh2hlx", 0, 0, "height levels in meter", VertintghHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertintgh> registration = RegisterEntry<Vertintgh>();

  int GH2HLX{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  size_t gridsize{};
  int numVars{};
  int numFullLevels{};
  int numHalfLevels{};

  bool extrapolate{};

  std::vector<bool> processVars;
  std::vector<bool> interpVars;

  int heightID_FL = -1, heightID_HL = -1;

  VarList varList1{};
  VarList varList2{};
  Varray2D<size_t> varnumMissVals;
  Field3DVector vardata1, vardata2;

  Varray<double> heightLevels;

  Varray<size_t> numMiss_FL, numMiss_HL;
  std::vector<int> vertIndex_FL;
  std::vector<int> vertIndex_HL;

public:
  void
  init() override
  {
    GH2HLX = module.get_id("gh2hlx");

    auto operatorID = cdo_operator_id();

    extrapolate = (operatorID == GH2HLX);
    if (extrapolate == false) extrapolate = getenv_extrapolate();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto zaxisID2 = create_zaxis_height(heightLevels);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    gridsize = vlist_check_gridsize(vlistID1);

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);
    auto memtype = varList1.vars[0].memType;

    auto stdnameHeight_FL = var_stdname(geometric_height_at_full_level_center);
    auto stdnameHeight_HL = var_stdname(geometric_height_at_half_level_center);

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto stdname = string_to_lower(varList1.vars[varID].stdname);
      if (stdname == stdnameHeight_FL) heightID_FL = varID;
      if (stdname == stdnameHeight_HL) heightID_HL = varID;
    }

    if (-1 == heightID_FL && -1 == heightID_HL)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        auto stdname = string_to_lower(var.stdname);
        if (stdname == "height" && var.units == "m")
        {
          // clang-format off
          if      (var.longname == "geometric height at full level center") heightID_FL = varID;
          else if (var.longname == "geometric height at half level center") heightID_HL = varID;
          // clang-format on
        }
      }
    }

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      if (-1 != heightID_FL) cdo_print("  %s -> %s", stdnameHeight_FL, varList1.vars[heightID_FL].name);
      if (-1 != heightID_HL) cdo_print("  %s -> %s", stdnameHeight_HL, varList1.vars[heightID_HL].name);
    }

    if (-1 == heightID_FL && -1 == heightID_HL) cdo_abort("%s not found!", stdnameHeight_FL);

    auto zaxisID_FL = (-1 == heightID_FL) ? -1 : varList1.vars[heightID_FL].zaxisID;
    auto zaxisID_HL = (-1 == heightID_HL) ? -1 : varList1.vars[heightID_HL].zaxisID;
    numFullLevels = (-1 == zaxisID_FL) ? 0 : varList1.vars[heightID_FL].nlevels;
    numHalfLevels = (-1 == zaxisID_HL) ? 0 : varList1.vars[heightID_HL].nlevels;

    auto numZaxes = vlistNumZaxis(vlistID1);
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto numLevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisID_FL || zaxisID == zaxisID_HL
          || (is_height_axis(zaxisID) && (numLevels == numHalfLevels || numLevels == numFullLevels)))
        vlistChangeZaxis(vlistID2, zaxisID, zaxisID2);
    }

    varList2 = VarList(vlistID2);
    varList_set_memtype(varList2, memtype);

    if (!extrapolate) numMiss_FL.resize(heightLevels.size());
    if (!extrapolate) numMiss_HL.resize(heightLevels.size());

    if (-1 != heightID_FL) vertIndex_FL.resize(gridsize * heightLevels.size());
    if (-1 != heightID_HL) vertIndex_HL.resize(gridsize * heightLevels.size());

    processVars.resize(numVars);
    interpVars.resize(numVars);
    varnumMissVals.resize(numVars);
    vardata1.resize(numVars);
    vardata2.resize(numVars);

    auto maxLevels = std::max(std::max(numFullLevels, numHalfLevels), (int) heightLevels.size());

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto isHeightAxis = is_height_axis(var.zaxisID);

      if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var);

      interpVars[varID] = (var.zaxisID == zaxisID_FL || var.zaxisID == zaxisID_HL
                           || (isHeightAxis && (var.nlevels == numHalfLevels || var.nlevels == numFullLevels)));

      if (interpVars[varID])
      {
        varnumMissVals[varID].resize(maxLevels, 0);
        vardata2[varID].init(varList2.vars[varID]);
      }
      else
      {
        if (isHeightAxis && var.nlevels > 1)
        {
          if (-1 == heightID_FL && -1 != heightID_HL && var.nlevels == (numHalfLevels - 1))
            cdo_abort("%s not found (needed for %s)!", stdnameHeight_FL, var.name);
          else if (-1 != heightID_FL && -1 == heightID_HL && var.nlevels == (numFullLevels + 1))
            cdo_abort("%s not found (needed for %s)!", stdnameHeight_HL, var.name);
          else
            cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1, var.name, var.nlevels);
        }
        varnumMissVals[varID].resize(var.nlevels);
      }
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      if (interpVars[varID] && var.isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field heightBottom;

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

      auto lreverse = true;
      if (-1 != heightID_FL && (tsID == 0 || !varList1.vars[heightID_FL].isConstant))
      {
        gen_vert_index(vertIndex_FL, heightLevels, vardata1[heightID_FL], gridsize, lreverse);
        if (!extrapolate)
        {
          heightBottom.init(varList1.vars[heightID_FL]);
          field_copy(vardata1[heightID_FL], numFullLevels - 1, heightBottom);
          gen_vert_index_mv(vertIndex_FL, heightLevels, gridsize, heightBottom, numMiss_FL, lreverse);
        }
      }

      if (-1 != heightID_HL && (tsID == 0 || !varList1.vars[heightID_HL].isConstant))
      {
        gen_vert_index(vertIndex_HL, heightLevels, vardata1[heightID_HL], gridsize, lreverse);
        if (!extrapolate)
        {
          heightBottom.init(varList1.vars[heightID_HL]);
          field_copy(vardata1[heightID_HL], numHalfLevels - 1, heightBottom);
          gen_vert_index_mv(vertIndex_HL, heightLevels, gridsize, heightBottom, numMiss_HL, lreverse);
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
            if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
              cdo_abort("Number of generalized height level differ from full/half level (param=%s)!", var.name);

            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
            }

            auto const &height3D = (var.nlevels == numFullLevels) ? vardata1[heightID_FL] : vardata1[heightID_HL];
            auto const &vertIndex3D = (var.nlevels == numFullLevels) ? vertIndex_FL : vertIndex_HL;
            vertical_interp_X(height3D, vardata1[varID], vardata2[varID], vertIndex3D, heightLevels, gridsize);

            if (!extrapolate)
            {
              auto numMiss = (var.nlevels == numFullLevels) ? numMiss_FL : numMiss_HL;
              varray_copy(heightLevels.size(), numMiss, varnumMissVals[varID]);
            }
          }

          for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                            varnumMissVals[varID][levelID]);
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
