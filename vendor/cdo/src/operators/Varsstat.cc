/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_stepstat.h"
#include "process_int.h"
#include "field_functions.h"

static void
check_unique_zaxis(int vlistID)
{
  auto numZaxes = vlistNumZaxis(vlistID);
  auto zaxisID = vlistZaxis(vlistID, 0);
  auto nlevels = zaxisInqSize(zaxisID);
  for (int index = 1; index < numZaxes; ++index)
  {
    if (nlevels != zaxisInqSize(vlistZaxis(vlistID, index))) cdo_abort("Number of level differ!");
  }
}

static void
check_unique_gridsize(int vlistID)
{
  auto numGrids = vlistNumGrids(vlistID);
  auto gridID = vlistGrid(vlistID, 0);
  auto gridsize = gridInqSize(gridID);
  for (int index = 0; index < numGrids; ++index)
  {
    if (gridsize != gridInqSize(vlistGrid(vlistID, index))) cdo_abort("Horizontal gridsize differ!");
  }
}

static void
set_attributes(const CdoVars &cdoVars1, int vlistID2, int varID2, int operatorID)
{
  auto const &var0 = cdoVars1[0];
  auto paramIsEqual = true;
  auto name = var0.name;
  auto param = var0.param;
  int nvars = cdoVars1.size();
  for (int varID = 1; varID < nvars; ++varID)
  {
    if (param != cdoVars1[varID].param || name != cdoVars1[varID].name)
    {
      paramIsEqual = false;
      break;
    }
  }

  if (!paramIsEqual) name = cdo_operator_name(operatorID);
  cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, name.c_str());
  vlistDefVarMissval(vlistID2, varID2, var0.missval);

  if (paramIsEqual)
  {
    if (param >= 0) vlistDefVarParam(vlistID2, varID2, param);
    if (var0.longname.size()) cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, var0.longname.c_str());
    if (var0.units.size()) cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, var0.units.c_str());
  }
}

class Varsstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Varsstat",
    .operators = { { "varsrange", FieldFunc_Range, 0, VarsstatHelp },
                   { "varsmin", FieldFunc_Min, 0, VarsstatHelp },
                   { "varsmax", FieldFunc_Max, 0, VarsstatHelp },
                   { "varssum", FieldFunc_Sum, 0, VarsstatHelp },
                   { "varsmean", FieldFunc_Mean, 0, VarsstatHelp },
                   { "varsavg", FieldFunc_Avg, 0, VarsstatHelp },
                   { "varsstd", FieldFunc_Std, 0, VarsstatHelp },
                   { "varsstd1", FieldFunc_Std1, 0, VarsstatHelp },
                   { "varsvar", FieldFunc_Var, 0, VarsstatHelp },
                   { "varsvar1", FieldFunc_Var1, 0, VarsstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Varsstat> registration = RegisterEntry<Varsstat>();

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int numLevels{};

  VarList varList1;

  cdo::StepStat1Dlevels stepStat;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    stepStat.init(operfunc);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    check_unique_zaxis(vlistID1);
    auto zaxisID = vlistZaxis(vlistID1, 0);
    numLevels = zaxisInqSize(zaxisID);

    check_unique_gridsize(vlistID1);

    auto timetype = varList1.vars[0].timeType;
    for (auto const &var : varList1.vars)
    {
      if (timetype != var.timeType) cdo_abort("Number of timesteps differ!");
    }

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    auto gridID = vlistGrid(vlistID1, 0);
    auto varID2 = vlistDefVar(vlistID2, gridID, zaxisID, timetype);
    set_attributes(varList1.vars, vlistID2, varID2, operatorID);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    stepStat.alloc(varList1, VARS_MEMTYPE);
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

        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        auto numSets = stepStat.var1(levelID).nsamp;
        stepStat.add_field(field, levelID, numSets);

        if (varID == 0 && stepStat.lvarstd) stepStat.moq(levelID);
      }

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        auto numSets = stepStat.var1(levelID).nsamp;
        if (numSets)
        {
          stepStat.process(levelID, numSets);

          cdo_def_field(streamID2, 0, levelID);
          cdo_write_field(streamID2, stepStat.var1(levelID));
          stepStat.var1(levelID).nsamp = 0;
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

    vlistDestroy(vlistID2);
  }
};
