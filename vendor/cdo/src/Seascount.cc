/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2007 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seascount   seascount         Seasonal counts
*/

#include <cdi.h>

#include "cdo_season.h"
#include "datetime.h"
#include "process_int.h"
#include "field_functions.h"

class Seascount : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seascount",
    .operators = { { "seascount" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Seascount> registration = RegisterEntry<Seascount>(module);

  CdiDateTime vDateTime0{};
  int seas0 = 0;
  int oldmon = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector2D varsData1;
    field2D_init(varsData1, varList1, FIELD_VEC);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);
    Field field;

    auto seasonStart = get_season_start();

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numFields = 0;
      int numSets = 0;
      auto newseas = false;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        auto month = decode_month(vDateTime.date);
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

          auto fieldsize = varsData1[varID][levelID].size;

          if (numSets == 0)
          {
            for (size_t i = 0; i < fieldsize; ++i) varsData1[varID][levelID].vec_d[i] = varsData1[varID][levelID].missval;
            varsData1[varID][levelID].numMissVals = fieldsize;
          }

          field.init(var);
          cdo_read_field(streamID1, field);

          field2_count(varsData1[varID][levelID], field);
        }

        vDateTime0 = vDateTime;
        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      taxisDefVdatetime(taxisID2, vDateTime0);
      cdo_def_timestep(streamID2, otsID);

      for (int fieldID = 0; fieldID < maxFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        if (otsID && varList1.vars[varID].isConstant) continue;

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, varsData1[varID][levelID]);
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
