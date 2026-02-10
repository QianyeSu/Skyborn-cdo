/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Split      splitrec        Split records
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include "util_files.h"
#include "util_string.h"

class Splitrec : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Splitrec",
    .operators = { { "splitrec", SplitHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Splitrec> registration = RegisterEntry<Splitrec>(module);

  CdoStreamID streamID1{};
  int vlistID1{ CDI_UNDEFID };
  std::string fileSuffix{};
  Field field{};
  VarList varList1{};
  bool dataIsUnchanged{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    dataIsUnchanged = data_is_unchanged();

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    int index = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        vlistClearFlag(vlistID1);
        vlistDefFlag(vlistID1, varID, levelID, true);

        auto vlistID2 = vlistCreate();
        cdo_vlist_copy_flag(vlistID2, vlistID1);

        index++;
        auto fileName = cdo_get_obase() + string_format("%06d", index);
        if (fileSuffix.size() > 0) fileName += fileSuffix;

        if (Options::cdoVerbose) cdo_print("create file %s", fileName);

        auto streamID2 = open_write(fileName);
        cdo_def_vlist(streamID2, vlistID2);

        auto varID2 = vlistFindVar(vlistID2, varID);
        auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);

        cdo_def_timestep(streamID2, 0);
        cdo_def_field(streamID2, varID2, levelID2);
        if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
        else
        {
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          cdo_write_field(streamID2, field);
        }

        cdo_stream_close(streamID2);
        vlistDestroy(vlistID2);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
  }
};
