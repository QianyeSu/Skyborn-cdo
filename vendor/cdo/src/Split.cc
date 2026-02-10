/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Split      splitcode       Split codes
      Split      splitparam      Split parameters
      Split      splitname       Split variables
      Split      splitlevel      Split levels
      Split      splitgrid       Split grids
      Split      splitzaxis      Split zaxis
      Split      splittabnum     Split table numbers
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_history.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"
#include "util_files.h"
#include "util_string.h"

#include <cassert>

static void
gen_filename(std::string &fileName, bool swapObase, std::string const &obase, std::string const &suffix)
{
  if (swapObase) fileName += obase;
  if (suffix.size() > 0) fileName += suffix;
}

class Split : public Process
{
  int
  split_code(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto numVars = varList1.numVars();
    std::vector<int> codes(numVars);

    int nsplit = 0;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      int index;
      for (index = 0; index < varID; ++index)
        if (var.code == varList1.vars[index].code) break;

      if (index == varID) codes[nsplit++] = var.code;
    }

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (codes[index] == var.code)
        {
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto format = (codes[index] > 9999) ? "%05d" : ((codes[index] > 999) ? "%04d" : "%03d");
      auto formatted = fileName + string_format(format, codes[index]);
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_param(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    char paramstr[32];
    auto numVars = varList1.numVars();
    std::vector<int> params(numVars);

    int nsplit = 0;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      int index;
      for (index = 0; index < varID; ++index)
        if (var.param == varList1.vars[index].param) break;

      if (index == varID) params[nsplit++] = var.param;
    }

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (params[index] == var.param)
        {
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      cdiParamToString(params[index], paramstr, sizeof(paramstr));

      auto formatted = fileName + paramstr;
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_name(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto numVars = varList1.numVars();
    auto nsplit = numVars;

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      int varID = index;
      auto const &var = varList1.vars[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        vlistDefIndex(vlistID1, varID, levelID, index);
        vlistDefFlag(vlistID1, varID, levelID, true);
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto formatted = fileName + var.name;
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_level(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto numVars = varList1.numVars();
    auto numZaxes = vlistNumZaxis(vlistID1);
    Varray<double> ftmp(999, 0.0);

    int nsplit = 0;
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        auto level = cdo_zaxis_inq_level(zaxisID, levelID);
        int i;
        for (i = 0; i < nsplit; ++i)
          if (is_equal(level, ftmp[i])) break;
        if (i == nsplit) ftmp[nsplit++] = level;
      }
    }

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);
    Varray<double> levels(nsplit);
    for (int index = 0; index < nsplit; ++index) levels[index] = ftmp[index];

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
          if (is_equal(levels[index], level))
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto formatted = fileName + string_format("%06g", levels[index]);
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_grid(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto nsplit = vlistNumGrids(vlistID1);
    auto numVars = varList1.numVars();

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);
    std::vector<int> gridIDs(nsplit);
    for (int index = 0; index < nsplit; ++index) gridIDs[index] = vlistGrid(vlistID1, index);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (gridIDs[index] == var.gridID)
        {
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto formatted = fileName + string_format("%02d", vlistGridIndex(vlistID1, gridIDs[index]) + 1);
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_zaxis(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto nsplit = vlistNumZaxis(vlistID1);
    auto numVars = varList1.numVars();

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);
    std::vector<int> zaxisIDs(nsplit);
    for (int index = 0; index < nsplit; ++index) zaxisIDs[index] = vlistZaxis(vlistID1, index);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (zaxisIDs[index] == var.zaxisID)
        {
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto formatted = fileName + string_format("%02d", vlistZaxisIndex(vlistID1, zaxisIDs[index]) + 1);
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

  int
  split_tabnum(bool swapObase, std::string const &fileSuffix, std::string const &fileName)
  {
    auto numVars = varList1.numVars();
    std::vector<int> tabnums(numVars);

    int nsplit = 0;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      int index;
      for (index = 0; index < varID; ++index)
        if (tabnum == tableInqNum(vlistInqVarTable(vlistID1, index))) break;

      if (index == varID) tabnums[nsplit++] = tabnum;
    }

    vlistIDs.resize(nsplit);
    streamIDs.resize(nsplit);

    for (int index = 0; index < nsplit; ++index)
    {
      vlistClearFlag(vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
        if (tabnums[index] == tabnum)
        {
          auto const &var = varList1.vars[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            vlistDefIndex(vlistID1, varID, levelID, index);
            vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }
      }

      auto vlistID2 = vlistCreate();
      cdo_vlist_copy_flag(vlistID2, vlistID1);
      vlistIDs[index] = vlistID2;

      auto formatted = fileName + string_format("%03d", tabnums[index]);
      gen_filename(formatted, swapObase, cdo_get_obase(), fileSuffix);

      streamIDs[index] = open_write(formatted);
    }

    return nsplit;
  }

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Split",
    .operators = { { "splitcode", SplitHelp },
                   { "splitparam", SplitHelp },
                   { "splitname", SplitHelp },
                   { "splitlevel", SplitHelp },
                   { "splitgrid", SplitHelp },
                   { "splitzaxis", SplitHelp },
                   { "splittabnum", SplitHelp } },
    .aliases = { { "splitvar", "splitname" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Split> registration = RegisterEntry<Split>(module);

  int SPLITCODE{}, SPLITPARAM{}, SPLITNAME{}, SPLITLEVEL{}, SPLITGRID{}, SPLITZAXIS{}, SPLITTABNUM{};

  std::vector<int> vlistIDs{};
  std::vector<CdoStreamID> streamIDs{};
  VarList varList1{};

  CdoStreamID streamID1{};
  int vlistID1{ CDI_UNDEFID };
  bool dataIsUnchanged{};

  int numSplit = 0;

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    SPLITCODE = module.get_id("splitcode");
    SPLITPARAM = module.get_id("splitparam");
    SPLITNAME = module.get_id("splitname");
    SPLITLEVEL = module.get_id("splitlevel");
    SPLITGRID = module.get_id("splitgrid");
    SPLITZAXIS = module.get_id("splitzaxis");
    SPLITTABNUM = module.get_id("splittabnum");

    auto operatorID = cdo_operator_id();

    auto swapObase = false;
    const char *uuidAttribute = nullptr;
    for (int i = 0; i < cdo_operator_argc(); ++i)
    {
      if (cdo_operator_argv(i) == "swap")
        swapObase = true;
      else if (cdo_operator_argv(i).find("uuid=") == 0)
        uuidAttribute = &cdo_operator_argv(i)[0] + 5;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(0));
    }

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    std::string fileName;
    if (!swapObase) fileName = cdo_get_obase();

    auto fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    if (operatorID == SPLITCODE) { numSplit = split_code(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITPARAM) { numSplit = split_param(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITTABNUM) { numSplit = split_tabnum(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITNAME) { numSplit = split_name(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITLEVEL) { numSplit = split_level(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITGRID) { numSplit = split_grid(swapObase, fileSuffix, fileName); }
    else if (operatorID == SPLITZAXIS) { numSplit = split_zaxis(swapObase, fileSuffix, fileName); }
    else { cdo_abort("not implemented!"); }

    assert(numSplit > 0);

    for (int index = 0; index < numSplit; ++index)
    {
      if (uuidAttribute) cdo_def_tracking_id(vlistIDs[index], uuidAttribute);

      cdo_def_vlist(streamIDs[index], vlistIDs[index]);
    }
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

      for (int index = 0; index < numSplit; ++index) cdo_def_timestep(streamIDs[index], tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        auto index = vlistInqIndex(vlistID1, varID, levelID);
        auto vlistID2 = vlistIDs[index];
        auto varID2 = vlistFindVar(vlistID2, varID);
        auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
        // printf("%d %d %d %d %d %d\n", index, vlistID2, varID, levelID, varID2, levelID2);

        cdo_def_field(streamIDs[index], varID2, levelID2);
        if (dataIsUnchanged) { cdo_copy_field(streamID1, streamIDs[index]); }
        else
        {
          auto const &var = varList1.vars[varID];
          field.init(var);
          cdo_read_field(streamID1, field);
          cdo_write_field(streamIDs[index], field);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);

    for (auto const &streamID : streamIDs) cdo_stream_close(streamID);
    for (auto const &vlistID : vlistIDs) vlistDestroy(vlistID);
  }
};
