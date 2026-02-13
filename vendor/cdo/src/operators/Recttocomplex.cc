/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

class Recttocomplex : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Recttocomplex",
    .operators = { { "recttocomplex" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Recttocomplex> registration = RegisterEntry<Recttocomplex>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};

  int vlistID3{};

  VarList varList1{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    streamID2 = cdo_open_read(1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    vlistID3 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    auto nvars = vlistNvars(vlistID3);
    for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array1(gridsizeMax);
    Varray<double> array2(gridsizeMax);
    Varray<double> array3(2 * gridsizeMax);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;
      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID3, varID, levelID);

        (void) cdo_inq_field(streamID2);

        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);
        cdo_read_field(streamID2, array2.data(), &numMissVals);

        auto gridsize = varList1.vars[varID].gridsize;
        for (size_t i = 0; i < gridsize; ++i)
        {
          array3[2 * i] = array1[i];
          array3[2 * i + 1] = array2[i];
        }

        cdo_write_field(streamID3, array3.data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID3);
  }
};
