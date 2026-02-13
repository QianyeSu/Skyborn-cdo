/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seaspctl   seaspctl        Seasonal percentiles
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "datetime.h"
#include "cdo_season.h"
#include "field_functions.h"

class Seaspctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seaspctl",
    .operators = { { "seaspctl", FieldFunc_Pctl, 0, SeaspctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Seaspctl> registration = RegisterEntry<Seaspctl>();

  TimeStat timestatDate{ TimeStat::MEAN };
  int seas0{ 0 };
  int oldmon{ 0 };

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};
  int taxisID4{};

  DateTimeList dtlist{};
  VarList varList1{};
  double pn{};

  HistogramSet hset{};

public:
  void
  init() override
  {
    operator_input_arg("percentile number");
    pn = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = cdo_stream_inq_vlist(streamID3);
    auto vlistID4 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID4);

    varList1 = VarList(vlistID1);
    VarList varList2(vlistID2);
    VarList varList3(vlistID3);

    varList_compare(varList1, varList2);
    varList_compare(varList1, varList3);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = vlistInqTaxis(vlistID3);
    // TODO - check that time axes 2 and 3 are equal

    taxisID4 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID4);
    vlistDefTaxis(vlistID4, taxisID4);

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);

    auto numVars = varList1.numVars();
    auto numSteps = varList1.numSteps();

    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    hset = HistogramSet(numVars, numSteps);

    for (auto const &var : varList1.vars) hset.createVarLevels(var.ID, var.nlevels, var.gridsize);
  }

  void
  run() override
  {
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);
    FieldVector constFields(maxFields);

    Field field1, field2;

    auto seasonStart = get_season_start();

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, otsID);
      if (numFields != cdo_stream_inq_timestep(streamID3, otsID))
        cdo_abort("Number of fields at time step %d of %s and %s differ!", otsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto vDateTime2 = taxisInqVdatetime(taxisID2);
      auto vDateTime3 = taxisInqVdatetime(taxisID3);
      if (cdiDateTime_isNE(vDateTime2, vDateTime3))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", otsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        auto const &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID2, field1);

        auto [varID2, levelID2] = cdo_inq_field(streamID3);
        field2.init(var);
        cdo_read_field(streamID3, field2);

        hset.defVarLevelBounds(varID2, levelID2, field1, field2);
      }

      int numSets = 0;
      auto newseas = false;
      while (numFields && (numFields = cdo_stream_inq_timestep(streamID1, tsID)))
      {
        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime1 = dtlist.vDateTime(numSets);

        auto month = decode_month(vDateTime1.date);
        auto newmon = month;
        if (seasonStart == SeasonStart::DEC && newmon == 12) newmon = 0;

        auto seas = month_to_season(month);

        if (numSets == 0)
        {
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

          if (tsID == 0 && var.isConstant)
          {
            constFields[fieldID].init(var);
            cdo_read_field(streamID1, constFields[fieldID]);
          }
          else
          {
            field1.init(var);
            cdo_read_field(streamID1, field1);

            hset.addVarLevelValues(varID, levelID, field1);
          }
        }

        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      dtlist.stat_taxis_def_timestep(taxisID4, numSets);
      cdo_def_timestep(streamID4, otsID);

      for (int fieldID = 0; fieldID < maxFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        auto const &var = varList1.vars[varID];
        if (otsID && var.isConstant) continue;

        cdo_def_field(streamID4, varID, levelID);

        if (var.isConstant) { cdo_write_field(streamID4, constFields[fieldID]); }
        else
        {
          field1.init(var);
          hset.getVarLevelPercentiles(field1, varID, levelID, pn);
          cdo_write_field(streamID4, field1);
        }
      }

      if (numFields == 0) break;
      otsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID4);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
