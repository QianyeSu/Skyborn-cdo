/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Writerandom writerandom
*/

#include <cstdlib>
#include <cdi.h>

#include "process_int.h"

class Writerandom : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Writerandom",
    .operators = { { "writerandom" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Writerandom> registration = RegisterEntry<Writerandom>(module);

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
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      Varray2D<double> recdata(numFields);
      std::vector<int> recvarID(numFields), reclevelID(numFields), recindex(numFields);
      std::vector<size_t> recnumMissVals(numFields);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        recvarID[fieldID] = varID;
        reclevelID[fieldID] = levelID;
        recdata[fieldID].resize(varList1.vars[varID].gridsize);
        cdo_read_field(streamID1, recdata[fieldID].data(), &recnumMissVals[fieldID]);
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID) recindex[fieldID] = -1;

      for (int rindex = numFields - 1; rindex >= 0; rindex--)
      {
        auto index = (int) (rindex * ((double) std::rand()) / ((double) RAND_MAX));
        //	printf("rindex %d %d\n", rindex, index);
        int ipos = -1;
        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          if (recindex[fieldID] == -1) ipos++;
          if (recindex[fieldID] == -1 && ipos == index)
          {
            recindex[fieldID] = rindex;
            break;
          }
        }
      }

      // for ( int fieldID = 0; fieldID < numFields; fieldID++ ) printf("fieldID %d %d\n", fieldID, recindex[fieldID]);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
        if (recindex[fieldID] == -1) cdo_abort("Internal problem! Random initialize.");

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto rindex = recindex[fieldID];
        auto varID = recvarID[rindex];
        auto levelID = reclevelID[rindex];
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, recdata[rindex].data(), recnumMissVals[rindex]);
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
