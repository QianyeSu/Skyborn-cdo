/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Sort sortcode  Sort by code number
*/

#include <ranges>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"
#include "field_functions.h"

namespace
{
struct LevInfo
{
  int levelID{};
  double level{};
};

struct VarInfo
{
  int varID{};
  int numLevels{};
  int code{};
  std::string param{};
  std::string name{};
  std::vector<LevInfo> levInfo{};
};
}  // namespace

static void
varInfo_init(VarInfo &varInfo, int &varID, int levelID, CdoVar const &var)
{
  varInfo.varID = varID;
  varInfo.code = var.code;
  varInfo.param = param_to_string_zerofilled(var.param);
  varInfo.name = var.name;
  varInfo.levInfo[levelID].levelID = levelID;
  varInfo.levInfo[levelID].level = cdo_zaxis_inq_level(var.zaxisID, levelID);
}

static void
varsInfo_print(std::vector<VarInfo> const &varsInfo, std::string const &txt)
{
  int numVars = (int) varsInfo.size();
  for (int varID = 0; varID < numVars; varID++)
  {
    auto const &varInfo = varsInfo[varID];
    for (int levelID = 0; levelID < varInfo.numLevels; ++levelID)
      printf("%s: %d %s %d %d %d %g\n", txt.c_str(), varID, varInfo.name.c_str(), varInfo.code, varInfo.numLevels,
             varInfo.levInfo[levelID].levelID, varInfo.levInfo[levelID].level);
  }
}

class Sort : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Sort",
    .operators = { { "sortcode" }, { "sortparam" }, { "sortname" }, { "sortlevel" } },
    .aliases = { { "sortvar", "sortname" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Sort>();

  int SORTCODE{}, SORTPARAM{}, SORTNAME{}, SORTLEVEL{};
  bool compareLess = true;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  int operatorID{};

public:
  void
  init() override
  {
    SORTCODE = module.get_id("sortcode");
    SORTPARAM = module.get_id("sortparam");
    SORTNAME = module.get_id("sortname");
    SORTLEVEL = module.get_id("sortlevel");

    operatorID = cdo_operator_id();

    if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

    if (operatorID == SORTLEVEL && cdo_operator_argc() == 1)
    {
      auto iarg = parameter_to_int(cdo_operator_argv(0));
      if (iarg < 0) compareLess = false;
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);
    /*
    if ( operatorID == SORTCODE ) { vlistSortCode(vlistID2); }
    else if ( operatorID == SORTNAME ) { ; }
    else if ( operatorID == SORTLEVEL ) { ; }
    */
    if (operatorID == SORTLEVEL)
    {
      auto numZaxes = vlistNumZaxis(vlistID2);
      for (int index = 0; index < numZaxes; ++index)
      {
        auto zaxisID1 = vlistZaxis(vlistID2, index);
        auto zaxisID2 = zaxisDuplicate(zaxisID1);
        auto numLevels = zaxisInqSize(zaxisID2);
        Varray<double> levels(numLevels);
        cdo_zaxis_inq_levels(zaxisID2, levels.data());

        compareLess ? std::ranges::sort(levels, std::ranges::less()) : std::ranges::sort(levels, std::ranges::greater());

        zaxisDefLevels(zaxisID2, levels.data());
        vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    auto numVars = varList1.numVars();
    std::vector<VarInfo> varsInfo(numVars);
    for (auto const &var : varList1.vars)
    {
      varsInfo[var.ID].numLevels = var.nlevels;
      varsInfo[var.ID].levInfo.resize(var.nlevels);
    }

    FieldVector2D varDataList;
    field2D_init(varDataList, varList1, FIELD_VEC | FIELD_NAT);

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
        auto const &var = varList1.vars[varID];

        if (tsID == 0) { varInfo_init(varsInfo[varID], varID, levelID, var); }

        cdo_read_field(streamID1, varDataList[varID][levelID]);
      }

      if (tsID == 0)
      {
        if (Options::cdoVerbose) { varsInfo_print(varsInfo, "sort in"); }

        if (operatorID == SORTCODE) { std::ranges::sort(varsInfo, {}, &VarInfo::code); }
        else if (operatorID == SORTPARAM) { std::ranges::sort(varsInfo, {}, &VarInfo::param); }
        else if (operatorID == SORTNAME) { std::ranges::sort(varsInfo, {}, &VarInfo::name); }
        else if (operatorID == SORTLEVEL)
        {
          for (int varID = 0; varID < numVars; varID++)
          {
            if (compareLess) { std::ranges::sort(varsInfo[varID].levInfo, std::ranges::less(), &LevInfo::level); }
            else { std::ranges::sort(varsInfo[varID].levInfo, std::ranges::greater(), &LevInfo::level); }
          }
        }

        if (Options::cdoVerbose) { varsInfo_print(varsInfo, "sort out"); }
      }

      for (int varID = 0; varID < numVars; varID++)
      {
        auto const &varInfo = varsInfo[varID];
        auto const &var = varList1.vars[varInfo.varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto levelID2 = varInfo.levInfo[levelID].levelID;
          if (tsID == 0 || !var.isConstant)
          {
            cdo_def_field(streamID2, varInfo.varID, levelID);
            cdo_write_field(streamID2, varDataList[varInfo.varID][levelID2]);
          }
        }
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
