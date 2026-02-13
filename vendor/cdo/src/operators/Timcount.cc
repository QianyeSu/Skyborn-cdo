/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timcount    timcount          Time counts
      Hourcount   hourcount         Hourly counts
      Daycount    daycount          Daily counts
      Moncount    moncount          Monthly counts
      Yearcount   yearcount         Yearly counts
*/

#include <cdi.h>

#include "process_int.h"
#include "util_date.h"
#include "field_functions.h"

class Timcount : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timcount",
    .operators = { { "timcount", 0, CMP_DATE, nullptr },
                   { "yearcount", 0, CMP_YEAR, nullptr },
                   { "moncount", 0, CMP_MONTH, nullptr },
                   { "daycount", 0, CMP_DAY, nullptr },
                   { "hourcount", 0, CMP_HOUR, nullptr } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timcount> registration = RegisterEntry<Timcount>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int compareDate{};

  CdiDateTime vDateTime0{};
  CdiDateTime vDateTimeN{};

  VarList varList1;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    compareDate = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "No.");

    vlistDefNtsteps(vlistID2, (compareDate == CMP_DATE) ? 1 : -1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    FieldVector2D varsData1;
    field2D_init(varsData1, varList1, FIELD_VEC);

    Field field;
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numFields = 0;
      int numSets = 0;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID1);

        if (numSets == 0) vDateTime0 = vDateTime;

        if (date_is_neq(vDateTime, vDateTime0, compareDate))
        {
          cdo_add_steps(-1);
          break;
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);

          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

          auto fieldsize = varsData1[varID][levelID].size;

          if (numSets == 0)
          {
            for (size_t i = 0; i < fieldsize; ++i) varsData1[varID][levelID].vec_d[i] = varsData1[varID][levelID].missval;
            varsData1[varID][levelID].numMissVals = fieldsize;
          }

          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);

          field2_count(varsData1[varID][levelID], field);
        }

        vDateTimeN = vDateTime;
        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      taxisDefVdatetime(taxisID2, vDateTimeN);
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
