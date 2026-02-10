/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "process_int.h"
#include "util_string.h"
#include "dcw_reader.h"

static void
print_polygons(DCW_Lists &dcw_lists, std::string const &codeNames)
{
  auto codeList = split_string(codeNames, "\\+");

  dcw_sort_countries(dcw_lists);

  printf("# Digital Chart of the World\n");
  printf("# Region for country:");
  for (auto const &code : codeList) printf(" %s", code.c_str());
  printf("\n");

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  Region region;
  if (dcw_get_region(dcw_lists, codeList, region)) cdo_abort("dcw_get_region() failed!");

  printf("#   West=%g  East=%g  South=%g  North=%g\n", region.west, region.east, region.south, region.north);
  printf("#\n");

  dcw_print_polygons(dcw_lists, codeList);
}

class DCW_util : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "DCW_util",
    .operators = { { "dcw" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<DCW_util> registration = RegisterEntry<DCW_util>(module);

private:
  int argc = 0;
  DCW_Lists dcw_lists;

public:
  void
  init() override
  {
    if (dcw_load_lists(dcw_lists)) cdo_abort("dcw_load_lists() failed!");

    argc = cdo_operator_argc();
  }

  void
  run() override
  {
    if (argc == 0) { cdo_abort("Parameter missing (available keywords are path/countries)"); }
    else if (argc > 1) { cdo_abort("Too many parameter, max=1!"); }
    else
    {
      const std::string paramCountry = "country=";
      const std::string paramDCW = "dcw:";
      auto &param = cdo_operator_argv(0);
      // clang-format off
      if (param == "path" || param == "dir")    dcw_print_path();
      else if (param == "countries")            dcw_print_countries(dcw_lists);
      else if (param == "states")               dcw_print_states(dcw_lists);
      else if (param.starts_with(paramCountry)) print_polygons(dcw_lists, param.substr(paramCountry.size()));
      else if (param.starts_with(paramDCW))     print_polygons(dcw_lists, param.substr(paramDCW.size()));
      else
        cdo_abort("Unsupported parameter: %s", param);
      // clang-format on
    }
  }

  void
  close() override
  {
  }
};
