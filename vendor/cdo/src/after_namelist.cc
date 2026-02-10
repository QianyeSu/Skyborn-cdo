/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "afterburner.h"

static char *
amatch(char *msr, const char *sub)
{
  int nm = static_cast<int>(std::strlen(msr));
  int ns = static_cast<int>(std::strlen(sub));

  for (int i = 0; i < nm - ns; ++i)
    if (std::strncmp(msr + i, sub, ns) == 0) return (msr + i + ns);

  return nullptr;
}

int
scan_par_obsolete(char *namelist, const char *name, int def)
{
  int value;

  char *cp = amatch(namelist, name);

  if (cp == nullptr)
    value = def;
  else
    value = atoi(cp);
  /*
  fprintf(stdout, " %16.16s = %6d ", name, value);
  if ( value == def ) fprintf(stdout, " (default)\n");
  else                fprintf(stdout, "          \n");
  */
  return value;
}

int
scan_par(int verbose, char *namelist, const char *name, int def)
{
  int value;

  char *cp = amatch(namelist, name);

  if (cp == nullptr)
    value = def;
  else
    value = atoi(cp);

  if (verbose)
    {
      fprintf(stdout, " %16.16s = %6d ", name, value);
      if (value == def)
        fprintf(stdout, " (default)\n");
      else
        fprintf(stdout, "          \n");
    }

  return value;
}

int
scan_time(int verbose, char *namelist, int *hours, int max_hours)
{
  char *icp;
  int nrqh = 0;

  char *cp = amatch(namelist, "timesel");
  if (cp == nullptr)
    {
      hours[nrqh++] = -1;
      if (verbose) fprintf(stdout, " %16.16s = all\n", "timesel");
      return (nrqh);
    }

  int time = (int) std::strtol(cp, &icp, 10);

  while (icp != cp && nrqh < max_hours)
    {
      hours[nrqh++] = time;
      cp = icp;
      time = (int) std::strtol(cp, &icp, 10);
    }

  if (verbose)
    {
      fprintf(stdout, " %16.16s = ", "timesel");
      for (time = 0; time < nrqh; ++time) fprintf(stdout, " %02d", hours[time]);
      fprintf(stdout, "\n");
    }

  return nrqh;
}

void
scan_code(char *namelist, struct Variable *vars, int maxCodes, int *numCodes)
{
  char *icp;
  int ncodes = 0;

  char *cp = amatch(namelist, "code");
  if (cp != nullptr)
    {
      int code = (int) std::strtol(cp, &icp, 10);
      while (code > 0 && code < maxCodes)
        {
          ncodes++;
          vars[code].selected = 1;
          cp = icp;
          code = (int) std::strtol(cp, &icp, 10);
        }
    }

  *numCodes = ncodes;
}

void
scan_darray(char *namelist, const char *name, double *values, int maxValues, int *numValues)
{
  char *icp;
  double val;
  int nval = 0;

  char *cp = amatch(namelist, name);

  if (cp != nullptr)
    {
      val = std::strtod(cp, &icp);
      values[nval++] = val;
      cp = icp;
      val = std::strtod(cp, &icp);
      while (val > 0 && nval < maxValues)
        {
          values[nval++] = val;
          cp = icp;
          val = std::strtod(cp, &icp);
        }
    }

  *numValues = nval;
}
