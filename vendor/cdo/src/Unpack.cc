/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

class Unpack : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Unpack",
    .operators = { { "unpack", UnpackHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Unpack> registration = RegisterEntry<Unpack>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

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

    varList1 = VarList(vlistID1);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      double addoffset = 0.0, scalefactor = 1.0;
      auto haveAddoffset = (cdiInqKeyFloat(vlistID1, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
      auto haveScalefactor = (cdiInqKeyFloat(vlistID1, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
      if (haveAddoffset || haveScalefactor)
      {
        auto datatype = vlistInqVarDatatype(vlistID1, varID);
        if (datatype != CDI_DATATYPE_FLT64) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
        if (haveAddoffset) cdiDeleteKey(vlistID2, varID, CDI_KEY_ADDOFFSET);
        if (haveScalefactor) cdiDeleteKey(vlistID2, varID, CDI_KEY_SCALEFACTOR);
      }
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
    vlistDestroy(vlistID2);
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
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
