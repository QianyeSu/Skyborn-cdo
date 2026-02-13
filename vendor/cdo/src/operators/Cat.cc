/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Copy       cat             Concatenate datasets
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_timer.h"
#include "util_files.h"
#include "progress.h"
#include "cdo_options.h"

class Cat : public Process
{
  enum struct StreamMode
  {
    APPEND,
    CREATE
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Cat",
    .operators = { { "cat", CopyHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Cat> registration = RegisterEntry<Cat>();

  StreamMode streamMode{ StreamMode::APPEND };
  bool hasConstVars{ true };
  bool dataIsUnchanged{ false };
  int tsID2{ 0 };
  CdoStreamID streamID2;
  int vlistID2{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int numFiles{ 0 };

  void
  open_append(VarList const &varList1, int numSteps)
  {
    streamID2 = cdo_open_append(numFiles);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    taxisID2 = vlistInqTaxis(vlistID2);

    VarList varList2(vlistID2);
    varList_compare(varList1, varList2);

    tsID2 = varList2.numSteps();
    if (tsID2 == 0) tsID2 = 1;  // bug fix for time constant data only

    if (numSteps == 0) { hasConstVars = false; }
  }

  void
  open_write(int vlistID1, int taxisID1, VarList const &varList1, int numSteps, std::string const &ofilename)
  {
    if (Options::cdoVerbose) cdo_print("Output file doesn't exist, creating: %s", ofilename);

    streamMode = StreamMode::CREATE;
    streamID2 = cdo_open_write(numFiles);
    vlistID2 = vlistDuplicate(vlistID1);
    vlist_unpack(vlistID2);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if ((numSteps == 0 || numSteps == 1) && numFiles > 1)
    {
      if (numSteps == 0)
      {
        hasConstVars = false;
        auto numVars = varList1.numVars();
        for (int varID = 0; varID < numVars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
      }
      vlistDefNtsteps(vlistID2, -1);
    }

    cdo_def_vlist(streamID2, vlistID2);
  }

public:
  void
  init() override
  {
    operator_check_argc(0);

    dataIsUnchanged = data_is_unchanged();

    auto streamCnt = cdo_stream_cnt();
    numFiles = streamCnt - 1;
  }

  void
  run() override
  {
    Field field;
    cdo::Progress progress(get_id());

    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      cdo::timer timer;

      auto streamID1 = cdo_open_read(fileIdx);
      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1(vlistID1);
      auto numSteps = varList1.numSteps();

      if (fileIdx == 0)
      {
        if (numSteps == 1 && varList1.numVaryingVars() == 0) numSteps = 0;

        std::string ofilename = cdo_get_stream_name(numFiles);
        if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename)) { open_append(varList1, numSteps); }
        else { open_write(vlistID1, taxisID1, varList1, numSteps, ofilename); }
      }
      else
      {
        VarList varList2(vlistID2);
        varList_compare(varList1, varList2);
      }

      int tsID1 = 0;
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID1, tsID1);
        if (numFields == 0) break;

        auto fstatus = (numSteps > 1) ? fileIdx + (tsID1 + 1.0) / numSteps : fileIdx + 1.0;
        progress.update(fstatus / numFiles);

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID2);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);

          auto const &var1 = varList1.vars[varID];
          if (hasConstVars && tsID2 > 0 && tsID1 == 0 && var1.isConstant) continue;

          cdo_def_field(streamID2, varID, levelID);

          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            field.init(var1);
            cdo_read_field(streamID1, field);
            cdo_write_field(streamID2, field);
          }
        }

        tsID1++;
        tsID2++;
      }

      cdo_stream_close(streamID1);

      if (Options::cdoVerbose) cdo_print("Processed file: %s   %.2f seconds", cdo_get_stream_name(fileIdx), timer.elapsed());
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);

    if (streamMode == StreamMode::CREATE) { vlistDestroy(vlistID2); }
  }
};
