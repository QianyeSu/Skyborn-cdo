/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Merge      mergetime       Merge datasets sorted by date and time
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_rlimit.h"
#include "process_int.h"
#include "util_string.h"
#include "util_files.h"
#include "printinfo.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "progress.h"

namespace
{
struct StreamInfo
{
  CdoStreamID streamID;
  CdiDateTime vDateTime{};
  int vlistID{ -1 };
  int taxisID{ -1 };
  int tsID{ -1 };
  int numFields{ 0 };
  VarList varList;
  std::map<int, int> mapOfVarIDs;
};
}  // namespace

bool
getenv_skip_same_time()
{
  auto envString = getenv_string("SKIP_SAME_TIME");
  if (envString.size())
  {
    auto ival = std::stoi(envString);
    if (ival == 1)
    {
      if (Options::cdoVerbose) cdo_print("Set SKIP_SAME_TIME to %d", ival);
      return true;
    }
  }

  return false;
}

static int
open_all_files(int numFiles, std::vector<StreamInfo> &streamInfoList)
{
  int numSteps = 0;
  for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
  {
    if (Options::cdoVerbose) cdo_print("process: %s", cdo_get_stream_name(fileIdx));

    auto &si = streamInfoList[fileIdx];
    si.streamID = cdo_open_read(fileIdx);
    si.vlistID = cdo_stream_inq_vlist(si.streamID);
    si.taxisID = vlistInqTaxis(si.vlistID);
    si.varList = VarList(si.vlistID);
    if (numSteps >= 0)
    {
      if (si.varList.numSteps() > 0) { numSteps += si.varList.numSteps(); }
      else { numSteps = -1; }
    }
  }

  return numSteps;
}

static void
read_first_timestep(int numFiles, std::vector<StreamInfo> &streamInfoList)
{
  for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
  {
    auto &si = streamInfoList[fileIdx];
    si.tsID = 0;
    si.numFields = cdo_stream_inq_timestep(si.streamID, si.tsID);
    if (si.numFields == 0)
    {
      cdo_stream_close(si.streamID);
      si.streamID = CDO_STREAM_UNDEF;
    }
    else { si.vDateTime = taxisInqVdatetime(si.taxisID); }
  }
}

static void
get_parameter(bool &skipSameTime, MapFlag &mapFlag)
{
  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      if (key == "skip_same_time") { skipSameTime = parameter_to_bool(value); }
      else if (key == "names")
      {
        if (value == "union") { mapFlag = MapFlag::Left; }
        else if (value == "intersect") { mapFlag = MapFlag::Intersect; }
        else { cdo_abort("Invalid value for key >%s< (names=<union/intersect>)", key, value); }
      }
      else { cdo_abort("Invalid parameter key >%s<!", key); }
    }
  }
}

class Mergetime : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Mergetime",
    .operators = { { "mergetime", MergeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Mergetime> registration = RegisterEntry<Mergetime>(module);

  int tsID2{ 0 };
  int vlistID2{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  CdiDateTime lastDateTime{};
  Field field;

  CdoStreamID streamID2{};

  int numFiles{ 0 };
  bool dataIsUnchanged{ false };
  bool skipSameTime{ false };
  std::vector<StreamInfo> streamInfoList;
  int vlistFileIDmin{ -1 };
  int vlistFileIDmax{ -1 };

  int numSteps{ 0 };

  MapFlag mapFlag{ MapFlag::Undefined };

public:
  void
  init() override
  {
    skipSameTime = getenv_skip_same_time();

    get_parameter(skipSameTime, mapFlag);

    dataIsUnchanged = data_is_unchanged();

    numFiles = cdo_stream_cnt() - 1;
    streamInfoList.resize(numFiles);

    cdo::set_numfiles(numFiles + 8);

    numSteps = open_all_files(numFiles, streamInfoList);

    // check that the contents is always the same
    if (mapFlag == MapFlag::Undefined)
    {
      for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      {
        auto &si = streamInfoList[fileIdx];
        varList_compare(streamInfoList[0].varList, si.varList);
        for (auto const &var : si.varList.vars) si.mapOfVarIDs[var.ID] = var.ID;
      }
    }
    else
    {
      vlistFileIDmin = 0;
      vlistFileIDmax = 0;
      for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      {
        auto numVars = streamInfoList[fileIdx].varList.numVars();
        if (numVars < streamInfoList[vlistFileIDmin].varList.numVars()) vlistFileIDmin = fileIdx;
        if (numVars > streamInfoList[vlistFileIDmax].varList.numVars()) vlistFileIDmax = fileIdx;
      }

      auto const &varList2 = streamInfoList[(mapFlag == MapFlag::Intersect) ? vlistFileIDmin : vlistFileIDmax].varList;
      for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
      {
        auto &si = streamInfoList[fileIdx];
        varList_map(si.varList, varList2, mapFlag, si.mapOfVarIDs);
      }
    }

    // read the first time step
    read_first_timestep(numFiles, streamInfoList);

    std::string ofilename = cdo_get_stream_name(numFiles);
    if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
      cdo_abort("Outputfile %s already exists!", ofilename);

    streamID2 = cdo_open_write(numFiles);
  }

  void
  fill_missing_fields(StreamInfo const &si1, StreamInfo &si2)
  {
    auto maxVars = si1.varList.numVars();
    auto numVars = si2.varList.numVars();
    if (numVars < maxVars)
    {
      std::vector<short> missingIDs(maxVars, 0);
      for (int varID = 0; varID < numVars; ++varID) { missingIDs[si2.mapOfVarIDs[varID]] = 1; }
      for (int varID = 0; varID < maxVars; ++varID)
      {
        if (missingIDs[varID] != 1)
        {
          auto const &var1 = si1.varList.vars[varID];
          field.init(var1);
          field_fill(field, var1.missval);
          field.numMissVals = field.gridsize;
          for (int levelID = 0; levelID < var1.nlevels; levelID++)
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field);
          }
        }
      }
    }
  }

  void
  run() override
  {
    cdo::Progress progress(get_id());

    while (true)
    {
      auto processTimestep = true;

      int nextFileID = -1;
      CdiDateTime vDateTime{};
      for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
      {
        if (streamInfoList[fileIdx].streamID != CDO_STREAM_UNDEF)
        {
          auto vdate = cdiDate_get(streamInfoList[fileIdx].vDateTime.date);
          auto vtime = cdiTime_get(streamInfoList[fileIdx].vDateTime.time);
          if (nextFileID == -1 || vdate < cdiDate_get(vDateTime.date)
              || (vdate == cdiDate_get(vDateTime.date) && vtime < cdiTime_get(vDateTime.time)))
          {
            nextFileID = fileIdx;
            vDateTime = streamInfoList[fileIdx].vDateTime;
          }
        }
      }

      auto fileIdx = nextFileID;
      if (Options::cdoVerbose) cdo_print("nextstep = %d  vDateTime = %s", fileIdx, datetime_to_string(vDateTime));
      if (fileIdx == -1) break;

      auto &si = streamInfoList[fileIdx];

      if (skipSameTime && cdiDateTime_isEQ(vDateTime, lastDateTime))
      {
        cdo_print("Timestep %4d in stream %d (%s) already exists, skipped!", si.tsID + 1, si.streamID->get_id(),
                  datetime_to_string(vDateTime));
        processTimestep = false;
      }

      if (processTimestep)
      {
        if (numSteps > 1) progress.update((tsID2 + 1.0) / numSteps);

        if (tsID2 == 0)
        {
          auto vlistID1 = (mapFlag == MapFlag::Undefined)
                              ? si.vlistID
                              : ((mapFlag == MapFlag::Intersect) ? streamInfoList[vlistFileIDmin].vlistID
                                                                 : streamInfoList[vlistFileIDmax].vlistID);
          vlistID2 = vlistDuplicate(vlistID1);
          vlist_unpack(vlistID2);
          auto taxisID1 = vlistInqTaxis(si.vlistID);
          taxisID2 = taxisDuplicate(taxisID1);
          vlistDefTaxis(vlistID2, taxisID2);
          vlistDefNtsteps(vlistID2, numSteps);
          cdo_def_vlist(streamID2, vlistID2);
        }

        lastDateTime = vDateTime;

        cdo_taxis_copy_timestep(taxisID2, si.taxisID);
        cdo_def_timestep(streamID2, tsID2);

        for (int fieldID = 0; fieldID < si.numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(si.streamID);

          auto const &var = si.varList.vars[varID];
          if (tsID2 > 0 && si.tsID == 0 && var.isConstant) continue;

          auto varID2 = varID;
          if (mapFlag != MapFlag::Undefined)
          {
            auto it = si.mapOfVarIDs.find(varID);
            if (it == si.mapOfVarIDs.end()) continue;
            varID2 = it->second;
          }

          cdo_def_field(streamID2, varID2, levelID);

          if (dataIsUnchanged) { cdo_copy_field(si.streamID, streamID2); }
          else
          {
            field.init(var);
            cdo_read_field(si.streamID, field);
            cdo_write_field(streamID2, field);
          }
        }

        if (mapFlag == MapFlag::Left) fill_missing_fields(streamInfoList[vlistFileIDmax], si);

        tsID2++;
      }

      si.numFields = cdo_stream_inq_timestep(si.streamID, ++si.tsID);
      if (si.numFields == 0)
      {
        cdo_stream_close(si.streamID);
        si.streamID = CDO_STREAM_UNDEF;
      }
      else { si.vDateTime = taxisInqVdatetime(si.taxisID); }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    if (vlistID2 != CDI_UNDEFID) vlistDestroy(vlistID2);
    if (taxisID2 != CDI_UNDEFID) taxisDestroy(taxisID2);
  }
};
