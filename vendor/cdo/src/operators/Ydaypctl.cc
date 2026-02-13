/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydaypctl   ydaypctl        Multi-year daily percentiles
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "field_functions.h"

class Ydaypctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ydaypctl",
    .operators = { { "ydaypctl", FieldFunc_Pctl, 0, YdaypctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Ydaypctl> registration = RegisterEntry<Ydaypctl>();

  static const int MaxDays = 373;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};
  int taxisID4{};

  double pn{};

  VarList varList1{};

public:
  void
  init() override
  {
    operator_input_arg("percentile number");
    pn = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
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
    if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
    vlistDefTaxis(vlistID4, taxisID4);

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);
  }

  void
  run() override
  {
    Field field1, field2;
    std::vector<bool> vars1(MaxDays, false);
    CdiDateTime vDateTimes1[MaxDays]{};
    CdiDateTime vDateTimes2[MaxDays]{};
    HistogramSet hsets[MaxDays];
    long numSets[MaxDays] = { 0 };

    auto numVars = varList1.numVars();
    auto numSteps = varList1.numSteps();
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    FieldVector constFields(maxFields);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      if (numFields != cdo_stream_inq_timestep(streamID3, tsID))
        cdo_abort("Number of fields at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto vDateTime = taxisInqVdatetime(taxisID2);

      if (cdiDateTime_isNE(vDateTime, taxisInqVdatetime(taxisID3)))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

      auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes2[dayOfYear] = vDateTime;

      if (!vars1[dayOfYear])
      {
        vars1[dayOfYear] = true;
        hsets[dayOfYear].create(numVars, numSteps);

        for (auto const &var : varList1.vars) hsets[dayOfYear].createVarLevels(var.ID, var.nlevels, var.gridsize);
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        auto const &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID2, field1);

        (void) cdo_inq_field(streamID3);
        field2.init(var);
        cdo_read_field(streamID3, field2);

        hsets[dayOfYear].defVarLevelBounds(varID, levelID, field1, field2);
      }

      tsID++;
    }

    tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      const auto vDateTime = taxisInqVdatetime(taxisID1);

      // if (Options::cdoVerbose) cdo_print("process timestep: %d %s", tsID + 1, datetime_to_string(vDateTime));

      auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes1[dayOfYear] = vDateTime;

      if (!vars1[dayOfYear])
        cdo_abort("No data for day %d in %s and %s", dayOfYear, cdo_get_stream_name(1), cdo_get_stream_name(2));

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

          hsets[dayOfYear].addVarLevelValues(varID, levelID, field1);
        }
      }

      numSets[dayOfYear]++;
      tsID++;
    }

    int otsID = 0;
    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
      {
        if (decode_month_and_day(vDateTimes1[dayOfYear].date) != decode_month_and_day(vDateTimes2[dayOfYear].date))
          cdo_abort("Verification dates for the day %d of %s and %s are different!", dayOfYear, cdo_get_stream_name(0),
                    cdo_get_stream_name(1));

        taxisDefVdatetime(taxisID4, vDateTimes1[dayOfYear]);
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
            hsets[dayOfYear].getVarLevelPercentiles(field1, varID, levelID, pn);
            cdo_write_field(streamID4, field1);
          }
        }

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
