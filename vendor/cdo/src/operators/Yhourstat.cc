/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

   Yhourstat     yhourmin        Multi-year hourly minimum
   Yhourstat     yhourmax        Multi-year hourly maximum
   Yhourstat     yhourrange      Multi-year hourly range
   Yhourstat     yhoursum        Multi-year hourly sum
   Yhourstat     yhourmean       Multi-year hourly mean
   Yhourstat     yhouravg        Multi-year hourly average
   Yhourstat     yhourstd        Multi-year hourly standard deviation
   Yhourstat     yhourstd1       Multi-year hourly standard deviation (n-1)
   Yhourstat     yhourvar        Multi-year hourly variance
   Yhourstat     yhourvar1       Multi-year hourly variance (n-1)

   Dhourstat     dhourmin        Multi-day hourly minimum
   Dhourstat     dhourmax        Multi-day hourly maximum
   Dhourstat     dhourrange      Multi-day hourly range
   Dhourstat     dhoursum        Multi-day hourly sum
   Dhourstat     dhourmean       Multi-day hourly mean
   Dhourstat     dhouravg        Multi-day hourly average
   Dhourstat     dhourstd        Multi-day hourly standard deviation
   Dhourstat     dhourstd1       Multi-day hourly standard deviation (n-1)
   Dhourstat     dhourvar        Multi-day hourly variance
   Dhourstat     dhourvar1       Multi-day hourly variance (n-1)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "datetime.h"
#include "process_int.h"
#include "printinfo.h"
#include "progress.h"
#include "field_functions.h"

class Yhourstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Yhourstat",
    // clang-format off
    .operators = { { "yhourrange", FieldFunc_Range, 0, YhourstatHelp },
                   { "yhourmin", FieldFunc_Min, 0, YhourstatHelp },
                   { "yhourmax", FieldFunc_Max, 0, YhourstatHelp },
                   { "yhoursum", FieldFunc_Sum, 0, YhourstatHelp },
                   { "yhourmean", FieldFunc_Mean, 0, YhourstatHelp },
                   { "yhouravg", FieldFunc_Avg, 0, YhourstatHelp },
                   { "yhourstd", FieldFunc_Std, 0, YhourstatHelp },
                   { "yhourstd1", FieldFunc_Std1, 0, YhourstatHelp },
                   { "yhourvar", FieldFunc_Var, 0, YhourstatHelp },
                   { "yhourvar1", FieldFunc_Var1, 0, YhourstatHelp },
                   { "dhourrange", FieldFunc_Range, 1, DhourstatHelp },
                   { "dhourmin", FieldFunc_Min, 1, DhourstatHelp },
                   { "dhourmax", FieldFunc_Max, 1, DhourstatHelp },
                   { "dhoursum", FieldFunc_Sum, 1, DhourstatHelp },
                   { "dhourmean", FieldFunc_Mean, 1, DhourstatHelp },
                   { "dhouravg", FieldFunc_Avg, 1, DhourstatHelp },
                   { "dhourstd", FieldFunc_Std, 1, DhourstatHelp },
                   { "dhourstd1", FieldFunc_Std1, 1, DhourstatHelp },
                   { "dhourvar", FieldFunc_Var, 1, DhourstatHelp },
                   { "dhourvar1", FieldFunc_Var1, 1, DhourstatHelp },
                   { "dminuterange", FieldFunc_Range, 3, DminutestatHelp },
                   { "dminutemin", FieldFunc_Min, 3, DminutestatHelp },
                   { "dminutemax", FieldFunc_Max, 3, DminutestatHelp },
                   { "dminutesum", FieldFunc_Sum, 3, DminutestatHelp },
                   { "dminutemean", FieldFunc_Mean, 3, DminutestatHelp },
                   { "dminuteavg", FieldFunc_Avg, 3, DminutestatHelp },
                   { "dminutestd", FieldFunc_Std, 3, DminutestatHelp },
                   { "dminutestd1", FieldFunc_Std1, 3, DminutestatHelp },
                   { "dminutevar", FieldFunc_Var, 3, DminutestatHelp },
                   { "dminutevar1", FieldFunc_Var1, 3, DminutestatHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Yhourstat> registration = RegisterEntry<Yhourstat>();

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1;

  cdo::StepStat3D stepStat;

  bool ldaily{};
  bool lminute{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    auto f2 = cdo_operator_f2(operatorID);
    ldaily = (f2 == 1 || f2 == 3);
    lminute = (f2 == 3);

    stepStat.init(operfunc);

    operator_check_argc(0);

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
  }

  int
  get_stepIndex(CdiDateTime vDateTime, int MaxSteps)
  {
    return lminute ? decode_minute_of_day(vDateTime, MaxSteps)
                   : (ldaily ? decode_hour_of_day(vDateTime, MaxSteps) : decode_hour_of_year(vDateTime, MaxSteps));
  }

  void
  run() override
  {
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    constexpr auto timestatDate{ TimeStat::LAST };
    int MaxHours = ldaily ? 25 : 9301;  // year: 31*12*25 + 1
    int MaxMinutes = 1501;              // 25*60 + 1;
    int MaxSteps = lminute ? MaxMinutes : MaxHours;
    std::vector<DateTimeList> dtlist(MaxSteps);
    std::vector<int> rangeNumSets(MaxSteps, 0);
    Field field;

    stepStat.set_dimlen0(MaxSteps);
    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;

    auto calendar = taxisInqCalendar(taxisID1);
    for (int stepIndex = 0; stepIndex < MaxHours; ++stepIndex)
    {
      dtlist[stepIndex].set_stat(timestatDate);
      dtlist[stepIndex].set_calendar(calendar);
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

      auto stepIndex = get_stepIndex(vDateTime, MaxSteps);

      dtlist[stepIndex].taxis_inq_timestep(taxisID1, rangeNumSets[stepIndex]);

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

    for (int stepIndex = 0; stepIndex < MaxSteps; stepIndex++)
    {
      auto numSets = rangeNumSets[stepIndex];
      if (numSets)
      {
        cdo::fields_process_3D(stepIndex, fieldInfoList, varList1, stepStat, numSets);

        dtlist[stepIndex].stat_taxis_def_timestep(taxisID2, rangeNumSets[stepIndex]);
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
