/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/

#include "cdi.h"
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "field_functions.h"

class Yearmonstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Yearmonstat",
    .operators = { { "yearmonmean", FieldFunc_Mean, 0, YearmonstatHelp }, { "yearmonavg", FieldFunc_Avg, 0, YearmonstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Yearmonstat> registration = RegisterEntry<Yearmonstat>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int year0{}, month0{};
  int year{}, month{}, day{};

  DateTimeList dtlist{};

  CdiDateTime vDateTime0{};

  int calendar{};
  int operfunc{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();

    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    calendar = taxisInqCalendar(taxisID1);
    constexpr auto timestatDate{ TimeStat::MEAN };
    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(calendar);
  }

  void
  run() override
  {
    Field field;
    VarList varList1(vlistID1);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    FieldVector2D samp1, varsData1;
    field2D_init(samp1, varList1);
    field2D_init(varsData1, varList1, FIELD_VEC);

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numFields = 0;
      int numSets = 0;
      double dsets = 0.0;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime = dtlist.vDateTime(numSets);
        cdiDate_decode(vDateTime.date, &year, &month, &day);

        if (numSets == 0) year0 = year;

        if (year != year0)
        {
          cdo_add_steps(-1);
          break;
        }

        if (numSets > 0 && month == month0)
        {
          cdo_warning("   last timestep: %s", datetime_to_string(vDateTime0));
          cdo_warning("current timestep: %s", datetime_to_string(vDateTime));
          cdo_abort("Month does not change!");
        }

        auto dpm = days_per_month(calendar, year, month);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          auto const &var = varList1.vars[varID];

          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

          auto &rsamp1 = samp1[varID][levelID];
          auto &rvars1 = varsData1[varID][levelID];

          if (numSets == 0)
          {
            cdo_read_field(streamID1, rvars1);

            fieldc_mul(rvars1, dpm);

            if (rvars1.numMissVals || !rsamp1.empty())
            {
              if (rsamp1.empty()) rsamp1.resize(rvars1.size);
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
              if (rsamp1.empty()) rsamp1.resize(rvars1.size, dsets);
              field2_vincr(rsamp1, field, dpm);
            }

            field2_function(rvars1, field, operfunc);
          }
        }

        month0 = month;
        vDateTime0 = vDateTime;
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

      if (Options::cdoVerbose) cdo_print("%s  numSets = %ld", datetime_to_string(vDateTime0), numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo_def_timestep(streamID2, otsID);

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
