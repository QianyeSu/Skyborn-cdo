/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     setcodetab      Set parameter code table
     setcode         Set code number
     setparam        Set parameter identifier
     setname         Set variable name
     setstdname      Set standard name
     setunit         Set variable unit
     setlevel        Set level
     setltype        Set GRIB level type
     setmaxsteps     Set max timesteps
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"

static void
set_level(int vlistID2, double newlevel)
{
  auto numZaxes = vlistNumZaxis(vlistID2);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID1 = vlistZaxis(vlistID2, index);
    auto zaxisID2 = zaxisDuplicate(zaxisID1);
    auto nlevs = zaxisInqSize(zaxisID2);
    Varray<double> levels(nlevs);
    cdo_zaxis_inq_levels(zaxisID2, levels.data());
    levels[0] = newlevel;
    zaxisDefLevels(zaxisID2, levels.data());
    vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
  }
}

static void
set_ltype(int vlistID2, double newval)
{
  auto numZaxes = vlistNumZaxis(vlistID2);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID1 = vlistZaxis(vlistID2, index);
    auto zaxisID2 = zaxisDuplicate(zaxisID1);
    auto zaxistype = ZAXIS_GENERIC;
    zaxisChangeType(zaxisID2, zaxistype);
    cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, newval);
    vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
  }
}

class Set : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Set",
    .operators = { { "setcode", 0, 0, "code number", SetHelp },
                   { "setparam", 0, 0, "parameter identifier (format:code[.tabnum]ornum[.cat[.dis]])", SetHelp },
                   { "setname", 0, 0, "variable name", SetHelp },
                   { "setstdname", 0, 0, "standard name", SetHelp },
                   { "setunit", 0, 0, "variable unit", SetHelp },
                   { "setlevel", 0, 0, "level", SetHelp },
                   { "setltype", 0, 0, "GRIB level type", SetHelp },
                   { "settabnum", 0, 0, "GRIB table number", SetHelp },
                   { "setmaxsteps", 0, 0, "max. number of timesteps", SetHelp } },
    .aliases = { { "setvar", "setname" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Set> registration = RegisterEntry<Set>(module);

  int SETCODE{}, SETPARAM{}, SETNAME{}, SETSTDNAME{}, SETUNIT{}, SETLEVEL{}, SETLTYPE{}, SETTABNUM{}, SETMAXSTEPS{};
  int maxSteps = -1;
  int newval = -1, tabnum = 0;
  int newparam = 0;
  const char *newname = nullptr;
  const char *newstdname = nullptr;
  const char *newunit = nullptr;
  double newlevel = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

public:
  void
  init() override
  {
    SETCODE = module.get_id("setcode");
    SETPARAM = module.get_id("setparam");
    SETNAME = module.get_id("setname");
    SETSTDNAME = module.get_id("setstdname");
    SETUNIT = module.get_id("setunit");
    SETLEVEL = module.get_id("setlevel");
    SETLTYPE = module.get_id("setltype");
    SETTABNUM = module.get_id("settabnum");
    SETMAXSTEPS = module.get_id("setmaxsteps");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));
    if (operatorID == SETCODE || operatorID == SETLTYPE) { newval = parameter_to_int(cdo_operator_argv(0)); }
    else if (operatorID == SETPARAM) { newparam = string_to_param(cdo_operator_argv(0)); }
    else if (operatorID == SETNAME) { newname = cdo_operator_argv(0).c_str(); }
    else if (operatorID == SETSTDNAME) { newstdname = cdo_operator_argv(0).c_str(); }
    else if (operatorID == SETUNIT) { newunit = cdo_operator_argv(0).c_str(); }
    else if (operatorID == SETTABNUM) { tabnum = parameter_to_int(cdo_operator_argv(0)); }
    else if (operatorID == SETLEVEL) { newlevel = parameter_to_double(cdo_operator_argv(0)); }
    else if (operatorID == SETMAXSTEPS) { maxSteps = parameter_to_int(cdo_operator_argv(0)); }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);
    // vlistPrint(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    if (operatorID == SETCODE)
    {
      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID) vlistDefVarCode(vlistID2, varID, newval);
    }
    else if (operatorID == SETPARAM) { vlistDefVarParam(vlistID2, 0, newparam); }
    else if (operatorID == SETNAME) { cdiDefKeyString(vlistID2, 0, CDI_KEY_NAME, newname); }
    else if (operatorID == SETSTDNAME) { cdiDefKeyString(vlistID2, 0, CDI_KEY_STDNAME, newstdname); }
    else if (operatorID == SETUNIT) { cdiDefKeyString(vlistID2, 0, CDI_KEY_UNITS, newunit); }
    else if (operatorID == SETTABNUM)
    {
      auto tableID = tableDef(-1, tabnum, nullptr);
      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID) vlistDefVarTable(vlistID2, varID, tableID);
    }
    else if (operatorID == SETLEVEL) { set_level(vlistID2, newlevel); }
    else if (operatorID == SETLTYPE) { set_ltype(vlistID2, newval); }
    else if (operatorID == SETMAXSTEPS) { vlistDefNtsteps(vlistID2, maxSteps); }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
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

        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
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
