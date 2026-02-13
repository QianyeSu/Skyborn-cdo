/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

class Tee : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Tee",
    .operators = { { "tee", TeeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Tee> registration = RegisterEntry<Tee>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdiStreamID streamID3;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;
  VarList varList1;

public:
  void
  init() override
  {
    operator_check_argc(1);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    streamID2 = cdo_open_write(1);

    streamID3 = streamOpenWrite(cdo_operator_argv(0).c_str(), cdo_filetype());

    auto vlistID2 = vlistDuplicate(vlistID1);
    auto vlistID3 = vlistDuplicate(vlistID1);

    taxisID2 = taxisDuplicate(taxisID1);
    taxisID3 = taxisDuplicate(taxisID1);

    vlistDefTaxis(vlistID2, taxisID2);
    vlistDefTaxis(vlistID3, taxisID3);

    cdo_def_vlist(streamID2, vlistID2);
    streamDefVlist(streamID3, vlistID3);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID2, tsID);
      streamDefTimestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);

        streamDefField(streamID3, varID, levelID);
        if (field.memType == MemType::Float)
          streamWriteFieldF(streamID3, field.vec_f.data(), field.numMissVals);
        else
          streamWriteField(streamID3, field.vec_d.data(), field.numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
    streamClose(streamID3);
  }
};
