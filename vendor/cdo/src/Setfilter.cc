/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "util_string.h"
#include "cdo_options.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"

#include <utility>

static void
string_replace_name(std::string &str, std::string const &name, std::string const &replace)
{
  auto pos = str.find(name);
  if (pos != std::string::npos) { str.replace(pos, name.size(), replace); }
}

void
expand_filter_names(std::string &filterSpec)
{
  // clang-format off
  const std::vector<std::pair<std::string, std::string>> filterList = {
      { "zip",        "1" },
      { "deflate",    "1" },
      { "shuffle",    "2" },
      { "fletcher32", "3" },
      { "bzip2",      "307" },
      { "blosc",      "32001" },
      { "lz4",        "32004" },
      { "bshuf",      "32008" },
      { "zfp",        "32013" },
      { "zstd",       "32015" },
      { "sz",         "32017" },
      { "sz3",        "32024" },
      { "blosc2",     "32026" }
  };
  // clang-format on

  for (auto const &[key, value] : filterList) string_replace_name(filterSpec, key, value);
}

namespace
{
struct Parameter
{
  std::string filename;
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params;

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

      // clang-format off
      if (key == "filename")  params.filename = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static std::vector<std::string>
get_vars_filter(VarList const &varList, std::string const &filename)
{
  auto numVars = varList.numVars();
  std::vector<std::string> varsFilter(numVars);

  if (filename.size())
  {
    auto fp = std::fopen(filename.c_str(), "r");
    if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
    PMList pmlist;
    pmlist.read_namelist(fp, filename.c_str());
    auto &kvlist = pmlist.front();
    std::fclose(fp);
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

      for (auto const &var : varList.vars)
      {
        if (key == var.name)
        {
          auto const &value = kv.values[0];
          auto filter = parameter_to_word(value);
          varsFilter[var.ID] = filter;
        }
      }
    }
  }

  return varsFilter;
}

static void
set_key_filterspec(int vlistID, int varID, std::string const &filterSpec)
{
  cdiDefKeyString(vlistID, varID, CDI_KEY_FILTERSPEC, filterSpec.c_str());
}

class Setfilter : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setfilter",
    .operators = { { "setfilter", SetfilterHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setfilter> registration = RegisterEntry<Setfilter>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1;

public:
  void
  init() override
  {
    auto params = get_parameter();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    varList1 = VarList(vlistID1);

    auto varsFilter = get_vars_filter(varList1, params.filename);

    vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      if (varsFilter[varID].size() > 0)
      {
        auto filterSpec = string_to_lower(varsFilter[varID]);
        expand_filter_names(filterSpec);
        set_key_filterspec(vlistID2, varID, filterSpec);
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        if (numSteps > 0)
        {
          auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
          progress.update(fstatus);
        }

        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
