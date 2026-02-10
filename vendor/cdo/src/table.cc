/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <stdlib.h>

#include <cdi.h>
#include "table.h"
#include "cdo_output.h"
#include "util_files.h"
#include "util_string.h"

namespace cdo
{

int
define_table(std::string const &tablearg)
{
  auto tablename = tablearg.c_str();

  auto tableID = FileUtils::file_exists(tablename) ? tableRead(tablename) : CDI_UNDEFID;

  if (tableID == CDI_UNDEFID)
  {
    auto const &tablepath = getenv_string("CD_TABLEPATH");
    if (tablepath.size())
    {
      std::string tablefile = tablepath + "/" + tablename;
      if (FileUtils::file_exists(tablefile)) tableID = tableRead(tablefile.c_str());
    }
  }

  if (tableID == CDI_UNDEFID) tableID = tableInq(-1, 0, tablename);

  if (tableID == CDI_UNDEFID) cdo_abort("table <%s> not found", tablename);

  return tableID;
}

std::string
predefined_tables(int p_padding)
{
  const char *name;
  constexpr int id_padding = 4;
  int padding = p_padding + id_padding;
  int numTables = tableInqNumber();
  std::string tables{ "Predefined tables: " };
  for (int id = 0; id < numTables; id++)
  {
    if (id % 7 == 6) tables += std::string("\n") + std::string(padding, ' ');
    if ((name = tableInqNamePtr(id))) tables += std::string(name);
    if (id < numTables - 1) tables += ",";
  }
  return tables;
}

}  // namespace cdo
