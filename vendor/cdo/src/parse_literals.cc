/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <climits>
#include <cerrno>
#include <cstdlib>
#include <cfloat>

#include "cdi.h"
#include "parse_literals.h"

int
literal_get_datatype(std::string const &literal)
{
  if (!literal.empty())
  {
    char *endptr;
    errno = 0;
    long lval = std::strtol(literal.c_str(), &endptr, 10);
    if (errno == 0 && *endptr == 0)
      return CDI_DATATYPE_INT32;
    else if (errno == 0 && *(endptr + 1) == 0 && (*endptr == 's' || *endptr == 'b'))
    {
      if (*endptr == 's' && lval >= SHRT_MIN && lval <= SHRT_MAX)
        return CDI_DATATYPE_INT16;
      else if (*endptr == 'b' && lval >= SCHAR_MIN && lval <= SCHAR_MAX)
        return CDI_DATATYPE_INT8;
    }
    else
    {
      errno = 0;
      double dval = std::strtod(literal.c_str(), &endptr);
      if (errno == 0 && *endptr == 0)
        return CDI_DATATYPE_FLT64;
      else if (errno == 0 && *(endptr + 1) == 0)
      {
        if (*endptr == 'f' && dval >= -FLT_MAX && dval <= FLT_MAX) return CDI_DATATYPE_FLT32;
      }
    }
  }

  return -1;
}

int
literals_find_datatype(int n, std::vector<std::string> const &literals)
{
  int dtype = -1;

  if (n)
  {
    dtype = literal_get_datatype(literals[0]);
    if (dtype != -1)
      for (int i = 1; i < n; ++i)
      {
        auto xtype = literal_get_datatype(literals[i]);
        if (dtype != xtype)
        {
          if (xtype == CDI_DATATYPE_FLT32 || xtype == CDI_DATATYPE_FLT64)
          {
            dtype = (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64) ? std::max(dtype, xtype) : xtype;
          }
          else { dtype = (!(dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)) ? std::max(dtype, xtype) : xtype; }
        }
      }
  }

  return dtype;
}

int
literal_to_int(std::string const &literal)
{
  int ival = INT_MAX;

  if (!literal.empty())
  {
    char *endptr;
    ival = std::strtol(literal.c_str(), &endptr, 10);
  }

  return ival;
}

double
literal_to_double(std::string const &literal)
{
  double dval = DBL_MAX;

  if (!literal.empty())
  {
    char *endptr;
    dval = std::strtod(literal.c_str(), &endptr);
  }

  return dval;
}

#ifdef TEST_LITERAL
int
main(void)
{
  std::vector<std::string> literals
      = { "127b", "-32768s", "-2147483647", "-1.e+36f", "1.e+308", "temperature", "surface pressure", "1000." };
  int nliterals = literals.size();

  for (int i = 0; i < nliterals; ++i)
  {
    int dtype = literal_get_datatype(literals[i]);
    printf("%d %s type = %d", i + 1, literals[i].c_str(), dtype);
    if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
    {
      printf("  ival = %d", literal_to_int(literals[i]));
    }
    else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64) { printf("  dval = %g", literal_to_double(literals[i])); }
    else { printf("  sval = '%s'", literals[i].c_str()); }

    printf("\n");
  }

  return 0;
}
#endif
