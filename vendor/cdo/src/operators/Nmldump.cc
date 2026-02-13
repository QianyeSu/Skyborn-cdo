/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdi.h"

#include "cdo_options.h"
#include "process_int.h"
#include "parse_literals.h"
#include "pmlist.h"
#include "util_string.h"

static void
print_values(int numValues, std::vector<std::string> const &values)
{
  char fltstr[128];
  if (numValues)
  {
    auto dataType = literals_find_datatype(numValues, values);
    for (int i = 0; i < numValues; ++i)
    {
      auto const &value = values[i];
      if (i) printf(", ");
      switch (dataType)
      {
        case CDI_DATATYPE_INT8: printf("%db", literal_to_int(value)); break;
        case CDI_DATATYPE_INT16: printf("%ds", literal_to_int(value)); break;
        case CDI_DATATYPE_INT32: printf("%d", literal_to_int(value)); break;
        case CDI_DATATYPE_FLT32:
          printf("%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), literal_to_double(value)));
          break;
        case CDI_DATATYPE_FLT64:
          printf("%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), literal_to_double(value)));
          break;
        default: printf("\"%s\"", value.c_str());
      }
    }
  }
}

void
kvldump(const PMList &pmlist)
{
  for (auto const &kvlist : pmlist)
  {
    auto const &listname = kvlist.name;
    if (!listname.empty()) printf("&%s\n", listname.c_str());
    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (!listname.empty()) printf("  ");
      printf("%s = ", key.c_str());
      print_values(kv.nvalues, kv.values);
      printf("\n");
    }
    if (!listname.empty()) printf("/\n");
  }
}

class Nmldump : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Nmldump",
    .operators = { { "nmldump" }, { "kvldump" } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<Nmldump> registration = RegisterEntry<Nmldump>();

  int NMLDUMP, KVLDUMP;
  int operatorID;
  PMList pmlist;

public:
  void
  init() override
  {
    NMLDUMP = module.get_id("nmldump");
    KVLDUMP = module.get_id("kvldump");

    operatorID = cdo_operator_id();

    operator_check_argc(0);
  }

  void
  run() override
  {
    pmlist.read_namelist(stdin, "STDIN");

    if (operatorID == NMLDUMP)
      pmlist.print();
    else if (operatorID == KVLDUMP)
      kvldump(pmlist);
  }

  void
  close() override
  {
  }
};
