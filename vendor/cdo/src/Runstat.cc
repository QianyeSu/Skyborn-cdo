/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Runstat    runrange        Running range
      Runstat    runmin          Running minimum
      Runstat    runmax          Running maximum
      Runstat    runsum          Running sum
      Runstat    runmean         Running mean
      Runstat    runavg          Running average
      Runstat    runvar          Running variance
      Runstat    runvar1         Running variance [Normalize by (n-1)]
      Runstat    runstd          Running standard deviation
      Runstat    runstd1         Running standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_stepstat.h"
#include "process_int.h"
#include "param_conversion.h"
#include "datetime.h"
#include "field_functions.h"

class Runstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Runstat",
    .operators = { { "runrange", FieldFunc_Range, 0, RunstatHelp },
                   { "runmin", FieldFunc_Min, 0, RunstatHelp },
                   { "runmax", FieldFunc_Max, 0, RunstatHelp },
                   { "runsum", FieldFunc_Sum, 0, RunstatHelp },
                   { "runmean", FieldFunc_Mean, 0, RunstatHelp },
                   { "runavg", FieldFunc_Avg, 0, RunstatHelp },
                   { "runstd", FieldFunc_Std, 0, RunstatHelp },
                   { "runstd1", FieldFunc_Std1, 0, RunstatHelp },
                   { "runvar", FieldFunc_Var, 0, RunstatHelp },
                   { "runvar1", FieldFunc_Var1, 0, RunstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Runstat> registration = RegisterEntry<Runstat>(module);

private:
  bool runstat_nomiss = false;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

  int ndates{};

  bool varsData2needed{};

  cdo::StepStat3D stepStat{};

public:
  void
  init() override
  {
    const auto envstr = getenv("RUNSTAT_NOMISS");
    if (envstr)
    {
      char *endptr;
      auto envval = (int) std::strtol(envstr, &endptr, 10);
      if (envval == 1) runstat_nomiss = true;
    }

    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);  // used in omp loop

    stepStat.init(operfunc);

    varsData2needed = (stepStat.lvarstd || stepStat.lrange);

    operator_input_arg("number of timesteps");
    operator_check_argc(1);
    ndates = parameter_to_int(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!stepStat.lminmax) vlist_unpack(vlistID2);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);
    // Number of timestep will be reduced compared to the input error handling in case of not enough timesteps is done per field
    auto numSteps = varList1.numSteps();
    if (numSteps != -1)
    {
      numSteps -= ndates - 1;
      if (numSteps > 0) vlistDefNtsteps(vlistID2, numSteps);
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    stepStat.set_dimlen0(ndates + 1);

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    for (int its = 0; its < ndates; its++)
    {
      field2D_init(stepStat.samp(its), varList1, !runstat_nomiss ? FIELD_VEC : 0);
      field2D_init(stepStat.var1(its), varList1, FIELD_VEC | VARS_MEMTYPE);
      field2D_init(stepStat.var2(its), varList1, varsData2needed ? FIELD_VEC : 0);
    }
  }

  void
  run() override
  {
    TimeStat timestatDate{ TimeStat::MEAN };
    DateTimeList dtlist;
    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    Vmask imask;

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
    FILL_FIRST_NDATES:
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0)
      {
        if (tsID < ndates)
          cdo_abort("File has less then %d timesteps!", ndates);
        else
          break;
      }

      auto numSteps = (tsID < ndates) ? tsID : ndates - 1;

      dtlist.taxis_inq_timestep(taxisID1, numSteps);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varIDx, levelIDx] = cdo_inq_field(streamID1);
        int varID = varIDx;  // needed for omp loop with intel icpx 2022.0.0
        int levelID = levelIDx;

        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

        auto &rsamp = stepStat.samp(numSteps, varID, levelID);
        auto &rvar1 = stepStat.var1(numSteps, varID, levelID);
        auto &rvar2 = stepStat.var2(numSteps, varID, levelID);

        auto fieldsize = rvar1.size;  // used in omp loop

        cdo_read_field(streamID1, rvar1);

        if (runstat_nomiss && rvar1.numMissVals) cdo_abort("Missing values supported was swichted off by env. RUNSTAT_NOMISS!");

        if (stepStat.lrange) field_copy(rvar1, rvar2);

        if (!runstat_nomiss)
        {
          imask.resize(fieldsize);

          auto func = [&](auto const &v, std::decay_t<decltype(v[0])> missval)
          {
            for (size_t i = 0; i < fieldsize; ++i) { imask[i] = fp_is_not_equal(v[i], missval); }
          };
          field_operation(func, rvar1, rvar1.missval);

          for (size_t i = 0; i < fieldsize; ++i) rsamp.vec_d[i] = (double) imask[i];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int inp = 0; inp < numSteps; ++inp)
          {
            auto &samp = stepStat.samp(inp, varID, levelID).vec_d;
            for (size_t i = 0; i < fieldsize; ++i)
              if (imask[i]) samp[i]++;
          }
        }

        if (stepStat.lvarstd)
        {
          field2_moq(stepStat.var2(numSteps, varID, levelID), stepStat.var1(numSteps, varID, levelID));
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int inp = 0; inp < numSteps; ++inp)
          {
            field2_sumsumq(stepStat.var1(inp, varID, levelID), stepStat.var2(inp, varID, levelID), rvar1);
          }
        }
        else if (stepStat.lrange)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int inp = 0; inp < numSteps; ++inp)
          {
            field2_maxmin(stepStat.var1(inp, varID, levelID), stepStat.var2(inp, varID, levelID), rvar1);
          }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (int inp = 0; inp < numSteps; ++inp)
          {
            field2_function(stepStat.var1(inp, varID, levelID), rvar1, stepStat.operfunc);
          }
        }
      }

      tsID++;  // don't move this line

      if (tsID < ndates) goto FILL_FIRST_NDATES;

      auto numSets = ndates;
      cdo::fields_process_3D(0, fieldInfoList, varList1, stepStat, numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, ndates);
      cdo_def_timestep(streamID2, otsID);

      for (int fieldID = 0; fieldID < maxFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        if (otsID && varList1.vars[varID].isConstant) continue;

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, stepStat.var1(0, varID, levelID));
      }

      otsID++;

      dtlist.shift();

      stepStat.var1(ndates) = stepStat.var1(0);
      if (!runstat_nomiss) stepStat.samp(ndates) = stepStat.samp(0);
      if (varsData2needed) stepStat.var2(ndates) = stepStat.var2(0);

      for (int inp = 0; inp < ndates; ++inp)
      {
        stepStat.var1(inp) = stepStat.var1(inp + 1);
        if (!runstat_nomiss) stepStat.samp(inp) = stepStat.samp(inp + 1);
        if (varsData2needed) stepStat.var2(inp) = stepStat.var2(inp + 1);
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
