/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdo_options.h"
#include "process_int.h"
#include "pmlist.h"

static void
dump_cmor_table(const PMList &pmlist)
{
  printf("# Number of lists: %zu\n", pmlist.size());
  int i = 0;
  for (auto const &kvlist : pmlist)
  {
    printf("# list ID: %d;   Number of elements: %zu\n", i, kvlist.size());
    printf("&%s\n", kvlist.name.c_str());
    for (auto const &kv : kvlist) { printf("  %s = %s\n", kv.key.c_str(), kv.values[0].c_str()); }
    printf("/\n");
    ++i;
  }
}

static void
conv_cmor_table(const PMList &pmlist)
{
  const char *hname = "Header";
  const char *vname = "variable";
  // const char *aname = "axis";

  bool hasmissval = false;
  double missval = 0;

  for (auto const &kvlist : pmlist)
  {
    auto const &listname = kvlist.name;

    if (listname.starts_with(hname))
    {
      for (auto const &kv : kvlist)
      {
        auto const &key = kv.key;
        auto const &value = kv.values[0];
        if (key == "missing_value")
        {
          missval = std::stof(value);
          hasmissval = true;
        }
      }
    }
    else if (listname.starts_with(vname))
    {
      printf("&%s\n", "parameter");
      for (auto const &kv : kvlist)
      {
        auto const &key = kv.key;
        auto const &value = kv.values[0];
        int vlen = value.size();

        int start = 0;
        if (vlen > 1 && value[0] == '"' && value[vlen - 1] == '"')
        {
          vlen -= 2;
          start++;
        }

        char *ovalue = strdup(value.c_str() + start);
        for (int i = 1; i < vlen; ++i)
        {
          if (ovalue[i - 1] == '"' && ovalue[i] == '"')
          {
            ovalue[i - 1] = '\'';
            for (int j = i + 1; j < vlen; ++j) ovalue[j - 1] = ovalue[j];
            vlen -= 1;
          }
        }

        if (vlen)
        {
          if (key == "name" || key == "standard_name" || key == "out_name" || key == "type" || key == "valid_min"
              || key == "valid_max" || key == "ok_min_mean_abs" || key == "ok_max_mean_abs")
            printf("  %-15s = %s\n", key.c_str(), ovalue);
          else if (key == "long_name" || key == "units" || key == "cell_methods" || key == "cell_measures" || key == "comment")
            printf("  %-15s = \"%.*s\"\n", key.c_str(), vlen, ovalue);
        }

        std::free(ovalue);
      }
      if (hasmissval) printf("  %-15s = %g\n", "missing_value", missval);
      printf("/\n");
    }
  }
}

class CMOR_table : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CMOR_table",
    .operators = { { "dump_cmor_table" }, { "conv_cmor_table" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<CMOR_table> registration = RegisterEntry<CMOR_table>(module);

public:
  void
  init() override
  {
    auto DUMP_CMOR_TABLE = module.get_id("dump_cmor_table");
    auto CONV_CMOR_TABLE = module.get_id("conv_cmor_table");

    auto operatorID = cdo_operator_id();

    if (cdo_operator_argc() != 1) cdo_abort("Too few arguments!");
    auto const &filename = cdo_operator_argv(0);

    if (Options::cdoVerbose) cdo_print("Parse file: %s", filename);

    auto fp = std::fopen(filename.c_str(), "r");
    if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);

    PMList pmlist;
    pmlist.read_cmor_table(fp, filename.c_str());
    std::fclose(fp);

    if (operatorID == DUMP_CMOR_TABLE)
      dump_cmor_table(pmlist);
    else if (operatorID == CONV_CMOR_TABLE)
      conv_cmor_table(pmlist);
  }

  void
  run() override
  {
  }

  void
  close() override
  {
  }
};
