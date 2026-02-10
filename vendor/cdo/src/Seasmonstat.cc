/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seasmonstat   seasmonmean        Seasonal mean from monthly data
      Seasmonstat   seasmonavg         Seasonal average from monthly data
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "cdo_season.h"
#include "field_functions.h"

class Seasmonstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seasmonstat",
    .operators = { { "seasmonmean", FieldFunc_Mean, 0, nullptr }, { "seasmonavg", FieldFunc_Avg, 0, nullptr } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Seasmonstat> registration = RegisterEntry<Seasmonstat>(module);

  TimeStat timestatDate{ TimeStat::MEAN };
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTime1{};
  int numFields;
  int seas0 = 0;
  int oldmon = 0;
  int nseason = 0;
  int month0 = 0;
  int year{}, month{}, day{};

  DateTimeList dtlist{};
  int calendar{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  FieldVector2D samp1{}, varsData1{};

  int operfunc{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    calendar = taxisInqCalendar(taxisID1);
    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(calendar);

    varList1 = VarList(vlistID1);

    field2D_init(samp1, varList1);
    field2D_init(varsData1, varList1, FIELD_VEC);
  }

  void
  run() override
  {
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);
    Field field;

    auto seasonStart = get_season_start();
    auto seasonNames = get_season_name();

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      long numSets = 0;
      double dsets = 0.0;
      auto newseas = false;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields <= 0) break;

        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime = dtlist.vDateTime(numSets);
        cdiDate_decode(vDateTime.date, &year, &month, &day);

        auto newmon = month;
        if (seasonStart == SeasonStart::DEC && newmon == 12) newmon = 0;

        auto seas = month_to_season(month);

        if (numSets > 0 && month == month0)
        {
          cdo_warning("   last timestep: %s", datetime_to_string(vDateTime0));
          cdo_warning("current timestep: %s", datetime_to_string(vDateTime));
          cdo_abort("Month does not change!");
        }

        auto dpm = days_per_month(calendar, year, month);

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
          auto const &var = varList1.vars[varID];

          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

          auto &rsamp1 = samp1[varID][levelID];
          auto &rvars1 = varsData1[varID][levelID];

          auto fieldsize = rvars1.size;

          if (numSets == 0)
          {
            cdo_read_field(streamID1, rvars1);

            fieldc_mul(rvars1, dpm);

            if (rvars1.numMissVals || !rsamp1.empty())
            {
              if (rsamp1.empty()) rsamp1.resize(fieldsize);
              field2_vinit(rsamp1, rvars1, dpm);
            }
          }
          else
          {
            field.init(var);
            cdo_read_field(streamID1, field);

            fieldc_mul(field, dpm);

            if (field.numMissVals || !rsamp1.empty())
            {
              if (rsamp1.empty()) rsamp1.resize(fieldsize, dsets);
              field2_vincr(rsamp1, field, dpm);
            }

            field2_function(rvars1, field, operfunc);
          }
        }

        month0 = month;
        vDateTime1 = vDateTime;
        numSets++;
        dsets += dpm;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (var.isConstant) continue;
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto const &rsamp1 = samp1[varID][levelID];
          auto &rvars1 = varsData1[varID][levelID];
          if (!rsamp1.empty())
            field2_div(rvars1, rsamp1);
          else
            fieldc_div(rvars1, dsets);
        }
      }

      if (Options::cdoVerbose)
        cdo_print("season: %3d %3s  start: %s  end: %s ntimesteps: %ld", nseason, seasonNames[seas0],
                  datetime_to_string(vDateTime0), datetime_to_string(vDateTime1), numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

      if (numSets < 3)
        cdo_warning("Season %3d (%s) has only %d input time step%s!", otsID + 1, date_to_string(vDateTime0.date), numSets,
                    numSets == 1 ? "" : "s");

      for (int fieldID = 0; fieldID < maxFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        if (otsID && varList1.vars[varID].isConstant) continue;

        auto &rvars1 = varsData1[varID][levelID];
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, rvars1);
      }

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
