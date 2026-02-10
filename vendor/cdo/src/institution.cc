/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <fstream>

#include <cdi.h>

#include "cdo_default_values.h"
#include "cdo_output.h"

static int
read_institution(std::string const &filename)
{
  int lnr = 0;
  int nvar = 0, maxvar = 4;
  int center = CDI_UNDEFID, subcenter = CDI_UNDEFID;
  std::string name, longname;

  std::ifstream file(filename);
  if (!file.is_open()) cdo_abort("Open failed on: %s\n", filename);

  std::string line;
  while (std::getline(file, line))
    {
      lnr++;
      if (line[0] == '#') continue;
      if (nvar == maxvar) break;
      nvar++;

      while (std::isspace((int) line[0])) line.erase(0, 1);

      if (nvar == 1) maxvar = std::isdigit((int) line[0]) ? 4 : 2;

      if (nvar == 1 && maxvar == 4) center = std::stoi(line);

      if (nvar == 2 && maxvar == 4)
        {
          if (!std::isdigit((int) line[0])) cdo_abort("wrong format in line %d. Missing subcenter!", lnr);

          subcenter = std::stoi(line);
        }

      if ((nvar == 3 && maxvar == 4) || (nvar == 1 && maxvar == 2)) name = line;
      if ((nvar == 4 && maxvar == 4) || (nvar == 2 && maxvar == 2)) longname = line;
    }

  file.close();

  auto instID = institutInq(center, subcenter, name.c_str(), longname.c_str());
  if (instID == CDI_UNDEFID) instID = institutDef(center, subcenter, name.c_str(), longname.c_str());

  return instID;
}

void
define_institution(std::string const &instString)
{
  int instID = read_institution(instString);

  if (instID == CDI_UNDEFID) instID = institutInq(0, 0, instString.c_str(), nullptr);
  if (instID == CDI_UNDEFID) cdo_abort("institution <%s> not found", instString);

  CdoDefault::InstID = instID;
}
