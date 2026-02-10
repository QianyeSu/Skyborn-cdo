/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Trend      trend           Trend
*/

#include <cdi.h>

#include "field.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "workerthread.h"
#include "field_trend.h"
#include "cdo_omp.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "field_functions.h"

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

      // clang-format off
      if (key == "equal") tstepIsEqual = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Trend : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Trend",
    .operators = { { "trend", TrendHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 2, OnlyFirst },
  };
  inline static RegisterEntry<Trend> registration = RegisterEntry<Trend>(module);

  static const int numWork = 5;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int maxFields{};

  bool tstepIsEqual = true;

  VarList varList1{};
  VarList varList2{};
  std::vector<FieldInfo> fieldInfoList;

public:
  void
  init() override
  {
    get_parameter(tstepIsEqual);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    vlistDefNtsteps(vlistID2, 1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    for (auto &var : varList2.vars) var.memType = MemType::Double;

    maxFields = varList1.maxFields();
    fieldInfoList = std::vector<FieldInfo>(maxFields);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);

    streamID2 = cdo_open_write(1);
    streamID3 = cdo_open_write(2);

    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_vlist(streamID3, vlistID2);
  }

  void
  write_output(FieldVector3D const &work)
  {
    Field field2, field3;

    cdo_def_timestep(streamID2, 0);
    cdo_def_timestep(streamID3, 0);

    for (int fieldID = 0; fieldID < maxFields; ++fieldID)
    {
      auto [varID, levelID] = fieldInfoList[fieldID].get();

      auto const &var = varList2.vars[varID];
      field2.init(var);
      field3.init(var);

      calc_trend_param(work, field2, field3, varID, levelID);

      field_num_mv(field2);
      field_num_mv(field3);

      cdo_def_field(streamID2, varID, levelID);
      cdo_write_field(streamID2, field2);

      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, field3);
    }
  }

  void
  run_sync()
  {
    auto calendar = taxisInqCalendar(taxisID1);
    CheckTimeIncr checkTimeIncr;
    JulianDate julianDate0;
    CdiDateTime vDateTime{};
    double deltat1 = 0.0;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());
    Field field1;

    FieldVector3D work(numWork);
    for (auto &w : work) field2D_init(w, varList1, FIELD_VEC, 0);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
        if (numSteps > 0) progress.update(fstatus);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        fieldInfoList[fieldID].set(varID, levelID);
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        calc_trend_sum(work, field1, zj, varID, levelID);
      }

      tsID++;
    }

    taxisDefVdatetime(taxisID2, vDateTime);
    write_output(work);
  }

  static void
  fields_calc_trend_sum(FieldVector3D &work, FieldVector2D const &fields2D, std::vector<FieldInfo> const &fieldInfoList,
                        double zj) noexcept
  {
    for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      calc_trend_sum(work, fields2D[varID][levelID], zj, varID, levelID);
    }
  }

  void
  run_async()
  {
    auto calendar = taxisInqCalendar(taxisID1);
    CheckTimeIncr checkTimeIncr;
    JulianDate julianDate0;
    CdiDateTime vDateTime{};
    double deltat1 = 0.0;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    FieldVector3D work(numWork);
    for (auto &w : work) field2D_init(w, varList1, FIELD_VEC, 0);

    FieldVector3D fields3D(2);
    field2D_init(fields3D[0], varList1, FIELD_VEC | FIELD_NAT);
    field2D_init(fields3D[1], varList1, FIELD_VEC | FIELD_NAT);

    bool useTask = true;
    auto workerThread = useTask ? std::make_unique<WorkerThread>() : nullptr;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
        if (numSteps > 0) progress.update(fstatus);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        fieldInfoList[fieldID].set(varID, levelID);
        cdo_read_field(streamID1, fields3D[tsID % 2][varID][levelID]);
      }

      if (useTask && tsID > 0) workerThread->wait();

      std::function<void()> fields_calc_trend_sum_task
          = std::bind(fields_calc_trend_sum, std::ref(work), std::ref(fields3D[tsID % 2]), std::cref(fieldInfoList), zj);

      if (useTask) { workerThread->doAsync(fields_calc_trend_sum_task); }
      else { fields_calc_trend_sum_task(); }

      tsID++;
    }

    if (useTask) workerThread->wait();

    taxisDefVdatetime(taxisID2, vDateTime);
    write_output(work);
  }

  void
  run() override
  {
    auto runAsync = (Options::CDO_Async_Read > 0);
    runAsync ? run_async() : run_sync();
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
