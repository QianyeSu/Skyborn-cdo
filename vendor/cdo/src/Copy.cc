/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Copy       copy            Copy datasets
*/

#include <cdi.h>

#include "process_int.h"
#include "progress.h"
#include "cdo_options.h"

bool
is_fdb_stream(std::string const &filename)
{
  return (filename.size() >= 4 && filename.starts_with("fdb:"));
}

bool
is_fdb_copy(bool dataIsUnchanged, int numFiles)
{
  auto isFdbCopy = false;

  if (dataIsUnchanged)
  {
    isFdbCopy = is_fdb_stream(cdo_get_stream_name(numFiles));
    if (numFiles == 1 && !isFdbCopy) isFdbCopy = is_fdb_stream(cdo_get_stream_name(0));
  }

  return isFdbCopy;
}

class Copy : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Copy",
    // clang-format off
    .operators = { { "copy", CopyHelp },
                   { "clone", CopyHelp },
                   { "szip", CopyHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Copy>(module);

  int CLONE{}, SZIP{};

  bool hasConstantFields{ true };
  CdoStreamID streamID2{ CDO_STREAM_UNDEF };
  int vlistID2{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool dataIsUnchanged{ false };

  int numFiles{ 0 };
  bool isFdbCopy{ false };

  int operatorID{ 0 };

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    CLONE = module.get_id("clone");
    SZIP = module.get_id("szip");

    operatorID = cdo_operator_id();
    if (operatorID == SZIP)
    {
      Options::cdoCompType = CDI_COMPRESS_SZIP;
      Options::cdoCompLevel = 0;
    }

    operator_check_argc(0);

    auto streamCnt = cdo_stream_cnt();
    numFiles = streamCnt - 1;

    isFdbCopy = is_fdb_copy(dataIsUnchanged, numFiles);
  }

  void
  run() override
  {
    Field field;

    cdo::Progress progress(get_id());

    int tsID2 = 0;
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      if (Options::cdoVerbose) cdo_print("Process file: %s", cdo_get_stream_name(fileIdx));

      auto streamID1 = cdo_open_read(fileIdx);
      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1(vlistID1);

      if (fileIdx == 0)
      {
        vlistID2 = vlistDuplicate(vlistID1);
        taxisID2 = taxisDuplicate(taxisID1);
        vlistDefTaxis(vlistID2, taxisID2);

        auto numSteps = varList1.numSteps();
        if (numSteps == 1 && varList1.numVaryingVars() == 0) numSteps = 0;

        if ((numSteps == 0 || numSteps == 1) && numFiles > 1)
        {
          if (numSteps == 0)
          {
            hasConstantFields = false;
            auto numVars = varList1.numVars();
            for (int varID = 0; varID < numVars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
          }
          vlistDefNtsteps(vlistID2, -1);
        }
      }
      else
      {
        VarList varList2(vlistID2);
        varList_compare(varList1, varList2);
      }

      if (streamID2 == CDO_STREAM_UNDEF)
      {
        streamID2 = cdo_open_write(numFiles);
        cdo_def_vlist(streamID2, vlistID2);
      }

      auto numSteps = varList1.numSteps();
      int tsID1 = 0;
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID1, tsID1);
        if (numFields == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID2);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          double fstatus = fileIdx + ((numSteps >= 0) ? (tsID1 + (fieldID + 1.0) / numFields) / numSteps : 1.0);
          progress.update(fstatus / numFiles);

          auto [varID, levelID] = cdo_inq_field(streamID1);

          auto const &var1 = varList1.vars[varID];
          if (hasConstantFields && tsID2 > 0 && tsID1 == 0 && var1.isConstant) continue;

          cdo_def_field(streamID2, varID, levelID);

          if (dataIsUnchanged && (isFdbCopy || operatorID == CLONE || operatorID == SZIP)) { cdo_copy_field(streamID1, streamID2); }
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
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);

    if (vlistID2 != CDI_UNDEFID) vlistDestroy(vlistID2);
  }
};
