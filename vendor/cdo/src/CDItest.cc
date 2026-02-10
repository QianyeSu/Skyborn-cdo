/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:
*/

#include <cdi.h>

#include "param_conversion.h"
#include "cdo_timer.h"
#include "process_int.h"

class CDItest : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CDItest",
    .operators = { { "ncopy" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<CDItest> registration = RegisterEntry<CDItest>(module);

  int NCOPY{};
  bool dataIsUnchanged{};
  int max_copy{};

public:
  void
  init() override
  {
    dataIsUnchanged = false;
    // auto dataIsUnchanged = data_is_unchanged();

    NCOPY = module.get_id("ncopy");
    (void) (NCOPY);  // unused

    auto operatorID = cdo_operator_id();
    (void) (operatorID);  // unused

    //  operator_input_arg("Number of copies");
    max_copy = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 3;
  }

  void
  run() override
  {
    int n = 0;
    while (true)
    {
      cdo::timer timer;

      auto streamID1 = cdo_open_read(0);

      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      auto vlistID2 = vlistDuplicate(vlistID1);
      auto taxisID2 = taxisDuplicate(taxisID1);
      vlistDefTaxis(vlistID2, taxisID2);

      auto streamID2 = cdo_open_write(1);
      cdo_def_vlist(streamID2, vlistID2);

      Field field;

      VarList varList1(vlistID1);

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
          cdo_def_field(streamID2, varID, levelID);

          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            field.init(varList1.vars[varID]);
            cdo_read_field(streamID1, field);
            cdo_write_field(streamID2, field);
          }
        }

        tsID++;
      }

      cdo_stream_close(streamID1);
      cdo_stream_close(streamID2);

      vlistDestroy(vlistID2);
      taxisDestroy(taxisID2);

      n++;

      cdo_print("Copy number %d: %.2fs", n, timer.elapsed());

      if (n == max_copy) break;
    }
  }

  void
  close() override
  {
  }
};
