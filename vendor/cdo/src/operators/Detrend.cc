/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Detrend    detrend         Detrend
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "workerthread.h"
#include "field_trend.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "field_functions.h"
#include "arithmetic.h"

static void
get_parameter(bool &tstepIsEqual)
{
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

      if (key == "equal") { tstepIsEqual = parameter_to_bool(value); }
      else { cdo_abort("Invalid parameter key >%s<!", key); }
    }
  }
}

class Detrend : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Detrend",
    .operators = { { "detrend", DetrendHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Detrend> registration = RegisterEntry<Detrend>();

  static const int numWork = 5;

  DateTimeList dtlist{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  VarList varList1{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };

  bool tstepIsEqual{ true };

public:
  void
  init() override
  {
    get_parameter(tstepIsEqual);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  vars_calc_trend_param(FieldVector3D &work)
  {
    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        auto gridsize = var.gridsize;
        auto missval1 = var.missval;
        auto missval2 = var.missval;

        auto &paramA = work[0][varID][levelID].vec_d;
        auto &paramB = work[1][varID][levelID].vec_d;
        auto &sumj = work[0][varID][levelID].vec_d;
        auto &sumjj = work[1][varID][levelID].vec_d;
        auto const &sumjx = work[2][varID][levelID].vec_d;
        auto const &sumx = work[3][varID][levelID].vec_d;
        auto const &zn = work[4][varID][levelID].vec_d;

        auto trend_kernel = [&](auto i, auto is_EQ)
        {
          auto temp1 = SUBM(sumjx[i], DIVMX(MULM(sumj[i], sumx[i]), zn[i]));
          auto temp2 = SUBM(sumjj[i], DIVMX(MULM(sumj[i], sumj[i]), zn[i]));
          auto temp3 = DIVM(temp1, temp2);

          paramA[i] = SUBM(DIVMX(sumx[i], zn[i]), MULM(DIVMX(sumj[i], zn[i]), temp3));
          paramB[i] = temp3;
        };

        if (std::isnan(var.missval))
          for (size_t i = 0; i < gridsize; ++i) trend_kernel(i, fp_is_equal);
        else
          for (size_t i = 0; i < gridsize; ++i) trend_kernel(i, is_equal);
      }
    }
  }

  static void
  vars_sub_trend(FieldVector3D &work, FieldVector2D &varsData, VarList const &varList, double zj)
  {
    auto numVars = varList.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList.vars[varID];
      if (var.isConstant) continue;
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        auto &field = varsData[varID][levelID];
        auto const &paramA = work[0][varID][levelID];
        auto const &paramB = work[1][varID][levelID];
        sub_trend(zj, field, paramA, paramB);
      }
    }
  }

  static void
  vars_trend_sum(FieldVector3D &work, FieldVector2D const &varsData, VarList const &varList, double zj)
  {
    auto numVars = varList.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList.vars[varID];
      if (var.isConstant) continue;
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        calc_trend_sum(work, varsData[varID][levelID], zj, varID, levelID);
      }
    }
  }

  void
  run() override
  {
    auto runAsync = (Options::CDO_Async_Read > 0);
    auto workerThread = runAsync ? std::make_unique<WorkerThread>() : nullptr;

    auto calendar = taxisInqCalendar(taxisID1);
    CheckTimeIncr checkTimeIncr;
    JulianDate julianDate0;
    double deltat1 = 0.0;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    FieldVector3D varsData{};
    if (numSteps > 0) varsData.resize(numSteps);

    FieldVector3D work(numWork);
    for (auto &w : work) field2D_init(w, varList1, FIELD_VEC, 0);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (numSteps > 1) progress.update((tsID + 1.0) / numSteps, 0.0, 0.5);

      dtlist.taxis_inq_timestep(taxisID1, tsID);
      auto vDateTime = dtlist.vDateTime(tsID);
      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);
      field2D_init(varsData[tsID], varList1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field = varsData[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }

      if (runAsync && tsID > 0) workerThread->wait();

      std::function<void()> vars_trend_sum_task
          = std::bind(vars_trend_sum, std::ref(work), std::cref(varsData[tsID]), std::cref(varList1), zj);

      runAsync ? workerThread->doAsync(vars_trend_sum_task) : vars_trend_sum_task();

      tsID++;
    }

    if (runAsync) workerThread->wait();

    numSteps = tsID;

    vars_calc_trend_param(work);

    if (runAsync)
    {
      auto step = 0;
      auto vDateTime = dtlist.vDateTime(step);
      auto zj = tstepIsEqual ? (double) step : delta_time_step_0(step, calendar, vDateTime, julianDate0, deltat1);
      std::function<void()> vars_sub_trend_func
          = std::bind(vars_sub_trend, std::ref(work), std::ref(varsData[step]), std::ref(varList1), zj);
      workerThread->doAsync(vars_sub_trend_func);
    }

    for (tsID = 0; tsID < numSteps; ++tsID)
    {
      progress.update((tsID + 1.0) / numSteps, 0.5, 0.5);

      if (runAsync) workerThread->wait();
      auto step = runAsync ? tsID + 1 : tsID;

      if (step < numSteps)
      {
        auto vDateTime = dtlist.vDateTime(step);
        auto zj = tstepIsEqual ? (double) step : delta_time_step_0(step, calendar, vDateTime, julianDate0, deltat1);

        std::function<void()> vars_sub_trend_func
            = std::bind(vars_sub_trend, std::ref(work), std::ref(varsData[step]), std::cref(varList1), zj);

        runAsync ? workerThread->doAsync(vars_sub_trend_func) : vars_sub_trend_func();
      }

      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (tsID && var.isConstant) continue;
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto &field = varsData[tsID][varID][levelID];
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, field);
        }
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
