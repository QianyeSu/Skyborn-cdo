/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydaystat   ydayrange       Multi-year daily range
      Ydaystat   ydaymin         Multi-year daily minimum
      Ydaystat   ydaymax         Multi-year daily maximum
      Ydaystat   ydaysum         Multi-year daily sum
      Ydaystat   ydaymean        Multi-year daily mean
      Ydaystat   ydayavg         Multi-year daily average
      Ydaystat   ydayvar         Multi-year daily variance
      Ydaystat   ydayvar1        Multi-year daily variance [Normalize by (n-1)]
      Ydaystat   ydaystd         Multi-year daily standard deviation
      Ydaystat   ydaystd1        Multi-year daily standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "printinfo.h"
#include "progress.h"
#include "field_functions.h"

struct YstatParam
{
  int year{};
  bool yearMode{};
};

static YstatParam
setParameter(void)
{
  YstatParam params;
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

      if (key == "yearMode")
        params.yearMode = parameter_to_bool(value);
      else if (key == "year")
        params.year = parameter_to_int(value);
      else
        cdo_abort("Invalid parameter key >%s<!", key);
    }
  }

  return params;
}

class Ydaystat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ydaystat",
    .operators = { { "ydayrange", FieldFunc_Range, 0, YdaystatHelp },
                   { "ydaymin", FieldFunc_Min, 0, YdaystatHelp },
                   { "ydaymax", FieldFunc_Max, 0, YdaystatHelp },
                   { "ydaysum", FieldFunc_Sum, 0, YdaystatHelp },
                   { "ydaymean", FieldFunc_Mean, 0, YdaystatHelp },
                   { "ydayavg", FieldFunc_Avg, 0, YdaystatHelp },
                   { "ydaystd", FieldFunc_Std, 0, YdaystatHelp },
                   { "ydaystd1", FieldFunc_Std1, 0, YdaystatHelp },
                   { "ydayvar", FieldFunc_Var, 0, YdaystatHelp },
                   { "ydayvar1", FieldFunc_Var1, 0, YdaystatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Ydaystat> registration = RegisterEntry<Ydaystat>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  int maxFields{};
  std::vector<FieldInfo> fieldInfoList{};

  cdo::StepStat3D stepStat{};
  YstatParam params{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    stepStat.init(operfunc);

    params = setParameter();

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    if (!stepStat.lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxFields = varList1.maxFields();
    fieldInfoList = std::vector<FieldInfo>(maxFields);
  }

  void
  run() override
  {
    constexpr auto timestatDate{ TimeStat::LAST };
    constexpr int MaxDays = 373;
    constexpr int MaxSteps = MaxDays;
    std::vector<DateTimeList> dtLists(MaxSteps);
    std::vector<int> rangeNumSets(MaxSteps, 0);
    Field field;

    stepStat.set_dimlen0(MaxSteps);
    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;

    auto calendar = taxisInqCalendar(taxisID1);
    for (int stepIndex = 0; stepIndex < MaxSteps; ++stepIndex)
    {
      dtLists[stepIndex].set_stat(timestatDate);
      dtLists[stepIndex].set_calendar(calendar);
    }

    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

      auto vDateTime = taxisInqVdatetime(taxisID1);
      if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

      auto stepIndex = decode_day_of_year(vDateTime.date);
      if (stepIndex < 0 || stepIndex >= MaxSteps)
        cdo_abort("Day of year %d out of range (%s)!", stepIndex, datetime_to_string(vDateTime));

      dtLists[stepIndex].taxis_set_next_timestep(taxisID1);

      if (!stepStat.var1(stepIndex).size()) { stepStat.alloc(stepIndex, varList1, VARS_MEMTYPE); }

      auto numSets = rangeNumSets[stepIndex];
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
        stepStat.add_field(field, stepIndex, varID, levelID, numSets);
      }

      rangeNumSets[stepIndex]++;
      tsID++;
    }

    // set the year to the minimum of years found on output timestep
    if (params.yearMode)
    {
      int outyear = 1e9;
      for (int stepIndex = 0; stepIndex < MaxSteps; stepIndex++)
      {
        if (rangeNumSets[stepIndex])
        {
          auto numEntries = dtLists[stepIndex].size();
          auto const &dtInfo = dtLists[stepIndex].info();
          outyear = std::min(outyear, dtInfo[numEntries - 1].v.date.year);
        }
      }
      params.year = outyear;
    }

    for (int stepIndex = 0; stepIndex < MaxSteps; stepIndex++)
    {
      auto numSets = rangeNumSets[stepIndex];
      if (numSets)
      {
        cdo::fields_process_3D(stepIndex, fieldInfoList, varList1, stepStat, numSets);

        if (params.year) dtLists[stepIndex].set_year(params.year);
        dtLists[stepIndex].stat_taxis_def_timestep(taxisID2);
        cdo_def_timestep(streamID2, otsID);

        for (int fieldID = 0; fieldID < maxFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();
          if (otsID && varList1.vars[varID].isConstant) continue;

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, stepStat.var1(stepIndex, varID, levelID));
        }

        otsID++;
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
