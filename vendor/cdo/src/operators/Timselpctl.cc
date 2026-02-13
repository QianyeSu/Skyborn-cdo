/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timselpctl    timselpctl         Time range percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "datetime.h"
#include "field_functions.h"

class Timselpctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timselpctl",
    .operators = { { "timselpctl", FieldFunc_Pctl, 0, TimselpctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Timselpctl> registration = RegisterEntry<Timselpctl>();

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };

  CdoStreamID streamID3;
  int taxisID3;

  CdoStreamID streamID4;
  int taxisID4;

  int noffset = 0, nskip = 0;
  int ndates;

  double pn;

  VarList varList1;

public:
  void
  init() override
  {
    operator_input_arg("percentile number, numSets <,noffset <,nskip>>");

    auto nargc = cdo_operator_argc();
    if (nargc < 2) cdo_abort("Too few arguments! Need %d found %d.", 2, nargc);

    pn = parameter_to_double(cdo_operator_argv(0));
    ndates = parameter_to_int(cdo_operator_argv(1));
    if (nargc > 2) noffset = parameter_to_int(cdo_operator_argv(2));
    if (nargc > 3) nskip = parameter_to_int(cdo_operator_argv(3));

    if (Options::cdoVerbose) cdo_print("numSets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip);

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
  }

  void
  run() override
  {
    DateTimeList dtlist;
    dtlist.set_stat(TimeStat::MEAN);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    HistogramSet hset(varList1.numVars(), varList1.numSteps());
    for (auto const &var : varList1.vars) hset.createVarLevels(var.ID, var.nlevels, var.gridsize);

    Field field1, field2;

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    FieldVector constFields(maxFields);

    int tsID;
    for (tsID = 0; tsID < noffset; ++tsID)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
      }
    }
    int otsID = 0;
    if (tsID < noffset)
    {
      cdo_warning("noffset is larger than number of timesteps!");
      return;
    }

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
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);
          field1.init(varList1.vars[varID]);
          cdo_read_field(streamID2, field1);
        }

        auto [varID, levelID] = cdo_inq_field(streamID3);
        field2.init(varList1.vars[varID]);
        cdo_read_field(streamID3, field2);

        hset.defVarLevelBounds(varID, levelID, field1, field2);
      }

      int numSets = 0;
      if (numFields)
        for (numSets = 0; numSets < ndates; numSets++)
        {
          numFields = cdo_stream_inq_timestep(streamID1, tsID);
          if (numFields == 0) break;

          dtlist.taxis_inq_timestep(taxisID1, numSets);

          for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = cdo_inq_field(streamID1);

            if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

            auto const &var1 = varList1.vars[varID];
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

      for (int i = 0; i < nskip; ++i)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;
        tsID++;
      }

      if (numFields == 0) break;
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
