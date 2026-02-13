/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selrec     selrec          Select records
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

class Selrec : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selrec",
    .operators = { { "selrec", SelvarHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, OnlyFirst },
  };
  inline static RegisterEntry<Selrec> registration = RegisterEntry<Selrec>();
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int nsel{};
  std::vector<int> intarr{};

public:
  void
  init() override
  {

    operator_input_arg("records");

    intarr = cdo_argv_to_intarr(cdo_get_oper_argv());
    nsel = intarr.size();

    if (Options::cdoVerbose)
    {
      for (int i = 0; i < nsel; ++i) cdo_print("intarr entry: %d %d", i, intarr[i]);
    }

    streamID1 = cdo_open_read(0);

    auto filetype = cdo_inq_filetype(streamID1);
    if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C
        || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR)
      cdo_abort("This operator does not work on NetCDF data!");

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    int recordID = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        recordID++;
        auto [varID, levelID] = cdo_inq_field(streamID1);

        for (int i = 0; i < nsel; ++i)
        {
          if (recordID == intarr[i])
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_copy_field(streamID1, streamID2);

            break;
          }
        }
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
