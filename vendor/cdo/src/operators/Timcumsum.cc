/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timcumsum    timcumsum         Cumulative sum over time
*/

#include <cdi.h>

#include "process_int.h"
#include "field_functions.h"

void set_miss_to_const(Field &field, double rconst);

class Timcumsum : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timcumsum",
    .operators = { { "timcumsum", TimcumsumHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timcumsum> registration = RegisterEntry<Timcumsum>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;

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
    Field field;
    FieldVector2D varsData1;
    field2D_init(varsData1, varList1, FIELD_VEC);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &rvars1 = varsData1[varID][levelID];

        if (tsID == 0)
        {
          cdo_read_field(streamID1, rvars1);
          if (rvars1.numMissVals) { set_miss_to_const(rvars1, 0.0); }
        }
        else
        {
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          if (field.numMissVals) { set_miss_to_const(field, 0.0); }
          field2_sum(rvars1, field);
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, rvars1);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
