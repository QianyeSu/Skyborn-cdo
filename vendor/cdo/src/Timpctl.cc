/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timpctl    timpctl         Time percentiles
      Hourpctl   hourpctl        Hourly percentiles
      Daypctl    daypctl         Daily percentiles
      Monpctl    monpctl         Monthly percentiles
      Yearpctl   yearpctl        Yearly percentiles
*/

#include <cdi.h>

#include "util_date.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "datetime.h"
#include "field_functions.h"

class Timpctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timpctl",
    .operators = { { "timpctl", FieldFunc_Pctl, CMP_DATE, TimpctlHelp },
                   { "yearpctl", FieldFunc_Pctl, CMP_YEAR, YearpctlHelp },
                   { "monpctl", FieldFunc_Pctl, CMP_MONTH, MonpctlHelp },
                   { "daypctl", FieldFunc_Pctl, CMP_DAY, DaypctlHelp },
                   { "hourpctl", FieldFunc_Pctl, CMP_HOUR, HourpctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Timpctl> registration = RegisterEntry<Timpctl>(module);

  CdiDateTime vDateTime0{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};
  int taxisID4{};

  int compareDate{};
  double pn{};

  DateTimeList dtlist{};
  Field field1{}, field2{};

  VarList varList1{};
  HistogramSet hset{};

public:
  void
  init() override
  {
    operator_input_arg("percentile number");
    pn = parameter_to_double(cdo_operator_argv(0));

    auto operatorID = cdo_operator_id();
    compareDate = cdo_operator_f2(operatorID);

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

    if (cdo_operator_f2(operatorID) == 16) vlistDefNtsteps(vlistID4, 1);

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

    dtlist.set_stat(TimeStat::MEAN);
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
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID2, field1);

        (void) cdo_inq_field(streamID3);
        field2.init(varList1.vars[varID]);
        cdo_read_field(streamID3, field2);

        hset.defVarLevelBounds(varID, levelID, field1, field2);
      }

      int numSets = 0;
      while (numFields && (numFields = cdo_stream_inq_timestep(streamID1, tsID)))
      {
        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime1 = dtlist.vDateTime(numSets);
        if (numSets == 0) vDateTime0 = vDateTime1;

        if (date_is_neq(vDateTime1, vDateTime0, compareDate))
        {
          cdo_add_steps(-1);
          break;
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          auto const &var1 = varList1.vars[varID];

          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

          if (tsID == 0 && var1.isConstant)
          {
            constFields[fieldID].init(var1);
            cdo_read_field(streamID1, constFields[fieldID]);
          }
          else
          {
            field1.init(var1);
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
        auto const &var1 = varList1.vars[varID];
        if (otsID && var1.isConstant) continue;

        cdo_def_field(streamID4, varID, levelID);

        if (var1.isConstant) { cdo_write_field(streamID4, constFields[fieldID]); }
        else
        {
          field1.init(var1);
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
