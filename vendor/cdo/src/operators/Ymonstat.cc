/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymonstat   ymonrange       Multi-year monthly range
      Ymonstat   ymonmin         Multi-year monthly minimum
      Ymonstat   ymonmax         Multi-year monthly maximum
      Ymonstat   ymonsum         Multi-year monthly sum
      Ymonstat   ymonmean        Multi-year monthly mean
      Ymonstat   ymonavg         Multi-year monthly average
      Ymonstat   ymonvar         Multi-year monthly variance
      Ymonstat   ymonvar1        Multi-year monthly variance [Normalize by (n-1)]
      Ymonstat   ymonstd         Multi-year monthly standard deviation
      Ymonstat   ymonstd1        Multi-year monthly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "datetime.h"
#include "process_int.h"
#include "printinfo.h"
#include "progress.h"
#include "field_functions.h"

class Ymonstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ymonstat",
    .operators = { { "ymonrange", FieldFunc_Range, 0, YmonstatHelp },
                   { "ymonmin", FieldFunc_Min, 0, YmonstatHelp },
                   { "ymonmax", FieldFunc_Max, 0, YmonstatHelp },
                   { "ymonsum", FieldFunc_Sum, 0, YmonstatHelp },
                   { "ymonmean", FieldFunc_Mean, 0, YmonstatHelp },
                   { "ymonavg", FieldFunc_Avg, 0, YmonstatHelp },
                   { "ymonstd", FieldFunc_Std, 0, YmonstatHelp },
                   { "ymonstd1", FieldFunc_Std1, 0, YmonstatHelp },
                   { "ymonvar", FieldFunc_Var, 0, YmonstatHelp },
                   { "ymonvar1", FieldFunc_Var1, 0, YmonstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Ymonstat> registration = RegisterEntry<Ymonstat>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  cdo::StepStat3D stepStat{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

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

  void
  run() override
  {
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    constexpr auto timestatDate{ TimeStat::LAST };
    constexpr int MaxMonths = 17;
    constexpr int MaxSteps = MaxMonths;
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

      int stepIndex = vDateTime.date.month;
      if (stepIndex < 0 || stepIndex >= MaxSteps) cdo_abort("Month %d out of range!", stepIndex);

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

    for (int stepIndex = 0; stepIndex < MaxSteps; stepIndex++)
    {
      auto numSets = rangeNumSets[stepIndex];
      if (numSets)
      {
        cdo::fields_process_3D(stepIndex, fieldInfoList, varList1, stepStat, numSets);

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
