/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seasstat   seasrange       Seasonal range
      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasvar1        Seasonal variance [Normalize by (n-1)]
      Seasstat   seasstd         Seasonal standard deviation
      Seasstat   seasstd1        Seasonal standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_season.h"
#include "progress.h"
#include "field_functions.h"

class Seasstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seasstat",
    // clang-format off
    .operators = { { "seasrange", FieldFunc_Range, 0, SeasstatHelp },
                   { "seasmin",   FieldFunc_Min,   0, SeasstatHelp },
                   { "seasmax",   FieldFunc_Max,   0, SeasstatHelp },
                   { "seassum",   FieldFunc_Sum,   0, SeasstatHelp },
                   { "seasmean",  FieldFunc_Mean,  0, SeasstatHelp },
                   { "seasavg",   FieldFunc_Avg,   0, SeasstatHelp },
                   { "seasstd",   FieldFunc_Std,   0, SeasstatHelp },
                   { "seasstd1",  FieldFunc_Std1,  0, SeasstatHelp },
                   { "seasvar",   FieldFunc_Var,   0, SeasstatHelp },
                   { "seasvar1",  FieldFunc_Var1,  0, SeasstatHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Seasstat> registration = RegisterEntry<Seasstat>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  TimeStat timestatDate{ TimeStat::MEAN };

  int seas0 = 0;
  int oldmon = 0;
  int nseason = 0;

  cdo::StepStat2D stepStat{};

  VarList varList1{};

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

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    stepStat.alloc(varList1, VARS_MEMTYPE);
  }

  void
  run() override
  {
    std::vector<FieldInfo> fieldInfoList(varList1.maxFields());

    DateTimeList dtlist{};
    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    Field field;

    CdiDateTime vDateTime0{};
    CdiDateTime vDateTime1{};
    auto seasonStartIsDecember = (get_season_start() == SeasonStart::DEC);
    auto seasonNames = get_season_name();

    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numSets = 0;
      bool newseas = false;
      int numFields = 0;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime = dtlist.vDateTime(numSets);

        auto month = decode_month(vDateTime.date);
        auto newmon = month;
        if (seasonStartIsDecember && newmon == 12) newmon = 0;

        auto seas = month_to_season(month);

        if (numSets == 0)
        {
          nseason++;
          vDateTime0 = vDateTime;
          seas0 = seas;
          oldmon = newmon;
        }

        if (newmon < oldmon) newseas = true;

        if ((seas != seas0) || newseas)
        {
          cdo_add_steps(-1);
          break;
        }

        oldmon = newmon;

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          stepStat.add_field(field, varID, levelID, numSets);
        }

        vDateTime1 = vDateTime;
        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      cdo::fields_process(fieldInfoList, varList1, stepStat, numSets);

      if (Options::cdoVerbose)
        cdo_print("season: %3d %3s  start: %s  end: %s ntimesteps: %ld", nseason, seasonNames[seas0],
                  datetime_to_string(vDateTime0), datetime_to_string(vDateTime1), numSets);

      if (numSets < 3)
        cdo_warning("Season %3d (%s) has only %d input time step%s!", otsID + 1, date_to_string(vDateTime0.date), numSets,
                    (numSets == 1) ? "" : "s");

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo::write_out_stream(streamID2, fieldInfoList, varList1, stepStat, otsID);

      if (numFields == 0) break;
      otsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
