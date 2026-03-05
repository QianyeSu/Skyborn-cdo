/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include <cdi.h>

#include "cdo_stepstat.h"
#include "cdo_omp.h"
#include "process_int.h"
#include "field_functions.h"
#include "param_conversion.h"

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
                   { "varsvar1", FieldFunc_Var1, 0, VarsstatHelp },
                   { "varsskew", FieldFunc_Skew, 1, VarsstatHelp },
                   { "varskurt", FieldFunc_Kurt, 1, VarsstatHelp },
                   { "varsmedian", FieldFunc_Median, 1, VarsstatHelp },
                   { "varspctl", FieldFunc_Pctl, 1, VarsstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Varsstat>();

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operFunc{ 0 };
  int operMethod{ 0 };

  double pn{ 0 };

  int numLevels{};

  VarList varList1;

  cdo::StepStat1Dlevels stepStat;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operFunc = cdo_operator_f1(operatorID);
    operMethod = cdo_operator_f2(operatorID);

    auto lpctl = (operFunc == FieldFunc_Pctl);

    auto argc = cdo_operator_argc();

    if (lpctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
      argc--;
    }
    else { operator_check_argc(0); }

    if (operMethod == 0) stepStat.init(operFunc);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    check_unique_zaxis(vlistID1);
    auto zaxisID = vlistZaxis(vlistID1, 0);
    numLevels = zaxisInqSize(zaxisID);

    check_unique_gridsize(vlistID1);

    auto timeType = varList1.vars[0].timeType;
    for (auto const &var : varList1.vars)
    {
      if (timeType != var.timeType) cdo_abort("Number of timesteps differ!");
    }

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    auto gridID = vlistGrid(vlistID1, 0);
    auto varID2 = vlistDefVar(vlistID2, gridID, zaxisID, timeType);
    set_attributes(varList1.vars, vlistID2, varID2, operatorID);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    if (operMethod == 0) stepStat.alloc(varList1, VARS_MEMTYPE);
  }

  void
  run() override
  {
    if (operMethod == 0)
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
    else
    {
      Varray<double> array2(varList1.vars[0].gridsize);

      auto numVars = varList1.numVars();
      FieldVector2D fieldList2D(numVars);
      for (auto &f : fieldList2D) { f.resize(numLevels); }

      FieldVector workFields(Threading::ompNumMaxThreads);
      for (auto &work : workFields) work.resize(numVars);

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

          fieldList2D[varID][levelID].init(varList1.vars[varID]);
          cdo_read_field(streamID1, fieldList2D[varID][levelID]);
        }

        for (int levelID = 0; levelID < numLevels; ++levelID)
        {
          auto hasMissvals = false;
          for (int varID = 0; varID < numVars; ++varID)
            if (fieldList2D[varID][levelID].numMissVals > 0) hasMissvals = true;

          auto gridsize = varList1.vars[0].gridsize;
          auto missval = varList1.vars[0].missval;
          auto memType = varList1.vars[0].memType;

          std::atomic<size_t> atomicNumMiss{ 0 };
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < gridsize; ++i)
          {
            auto ompthID = cdo_omp_get_thread_num();

            auto &work = workFields[ompthID];
            work.missval = missval;
            work.numMissVals = 0;
            if (memType == MemType::Float)
              for (int k = 0; k < numVars; ++k) work.vec_d[k] = fieldList2D[k][levelID].vec_f[i];
            else
              for (int k = 0; k < numVars; ++k) work.vec_d[k] = fieldList2D[k][levelID].vec_d[i];

            if (hasMissvals)
              for (int k = 0; k < numVars; ++k)
              {
                if (fp_is_equal(work.vec_d[k], varList1.vars[k].missval))
                {
                  work.vec_d[k] = missval;
                  work.numMissVals++;
                }
              }

            auto lpctl = (operFunc == FieldFunc_Pctl);
            array2[i] = lpctl ? field_pctl(work, pn) : field_function(work, operFunc);

            if (fp_is_equal(array2[i], work.missval)) atomicNumMiss++;
          }

          size_t numMissVals = atomicNumMiss;

          cdo_def_field(streamID2, 0, levelID);
          cdo_write_field(streamID2, array2.data(), numMissVals);
        }

        tsID++;
      }
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
