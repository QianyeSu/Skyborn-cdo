/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "param_conversion.h"

#include <inttypes.h>
#include <climits>
#include <cassert>
#include <cctype>
#include <vector>
#include <cdi.h>

#include "cdo_output.h"
#include "util_string.h"
#include "constants.h"
#include "const.h"
#include "fill_1d.h"

const char *
parameter_to_word(const char *cstring)
{
  auto len = std::strlen(cstring);

  for (size_t i = 0; i < len; ++i)
  {
    int c = cstring[i];
    if (iscntrl(c) || isblank(c)) cdo_abort("Word parameter >%s< contains invalid character at position %zu!", cstring, i + 1);
  }

  if (len == 0) cdo_abort("Word parameter >%s< is empty!", cstring);

  return cstring;
}

static inline void
parameter_error(const char *name, const char *cstring, const char *endptr)
{
  cdo_abort("%s parameter >%s< contains invalid character at position %d!", name, cstring, (int) (endptr - cstring + 1));
}

long
parameter_to_bytes(std::string const &string)
{
  char *endptr = nullptr;
  auto numBytes = strtoimax(string.c_str(), &endptr, 10);
  if (*endptr)
  {
    switch (tolower((int) *endptr))
    {
        // clang-format off
      case 'k': numBytes *= 1024;       endptr++; break;
      case 'm': numBytes *= 1048576;    endptr++; break;
      case 'g': numBytes *= 1073741824; endptr++; break;
        // clang-format on
    }
  }
  if (*endptr) parameter_error("Bytes", string.c_str(), endptr);
  int len = string.size();
  if ((string[0] != '-' && len > 19) || len > 20)
    cdo_abort("Integer parameter %s too long (min=%ld max=%ld)!", string, LONG_MIN, LONG_MAX);
  return numBytes;
}

double
parameter_to_double(const char *const cstring)
{
  char *endptr = nullptr;
  auto fval = std::strtod(cstring, &endptr);
  if (*endptr == 'f') endptr++;
  if (*endptr) parameter_error("Float", cstring, endptr);
  return fval;
}

static intmax_t
parameter_to_imax(const char *const cstring)
{
  char *endptr = nullptr;
  auto ival = strtoimax(cstring, &endptr, 10);
  if (*endptr) parameter_error("Integer", cstring, endptr);
  int len = std::strlen(cstring);
  if ((cstring[0] != '-' && len > 19) || len > 20)
    cdo_abort("Integer parameter %s too long (min=%ld max=%ld)!", cstring, LONG_MIN, LONG_MAX);
  return ival;
}

int
parameter_to_int(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < INT_MIN || ival > INT_MAX)
    cdo_abort("Integer parameter >%s< out of range (min=%d max=%d)!", cstring, INT_MIN, INT_MAX);
  return static_cast<int>(ival);
}

long
parameter_to_long(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < LONG_MIN || ival > LONG_MAX)
    cdo_abort("Integer parameter >%s< out of range (min=%ld max=%ld)!", cstring, LONG_MIN, LONG_MAX);
  return static_cast<long>(ival);
}

size_t
parameter_to_size_t(const char *const cstring)
{
  auto ival = parameter_to_imax(cstring);
  if (ival < 0 || ival > LONG_MAX)
    cdo_abort("Unsigned integer parameter >%s< out of range (min=%ld max=%ld)!", cstring, 0L, LONG_MAX);
  return static_cast<size_t>(ival);
}

int
parameter_to_intlist(const char *const cstring)
{
  char *endptr = nullptr;
  auto ival = (int) std::strtol(cstring, &endptr, 10);
  if (*endptr && *endptr != '/' && (endptr - cstring) == 0) parameter_error("Integer", cstring, endptr);
  return ival;
}

std::string const &
parameter_to_word(std::string const &string)
{
  auto len = string.size();

  for (size_t i = 0; i < len; ++i)
  {
    int c = string[i];
    if (iscntrl(c) || isblank(c)) cdo_abort("Word parameter >%s< contains invalid character at position %zu!", string, i + 1);
  }

  if (len == 0) cdo_abort("Word parameter >%s< is empty!", string);

  return string;
}

bool
parameter_to_bool(std::string const &string)
{
  auto lstring = string_to_lower(string);

  if (lstring == "1" || lstring == "t" || lstring == "true") return true;
  if (lstring == "0" || lstring == "f" || lstring == "false") return false;

  cdo_abort("Boolean parameter >%s< contains invalid characters!", string);

  return false;
}

double
parameter_to_double(std::string const &string)
{
  return parameter_to_double(string.c_str());
}

int
parameter_to_int(std::string const &string)
{
  return parameter_to_int(string.c_str());
}

long
parameter_to_long(std::string const &string)
{
  return parameter_to_long(string.c_str());
}

size_t
parameter_to_size_t(std::string const &string)
{
  return parameter_to_size_t(string.c_str());
}

int
parameter_to_intlist(std::string const &string)
{
  return parameter_to_intlist(string.c_str());
}

double
radius_str_to_meter(std::string const &string)
{
  char *endptr = nullptr;
  auto radius = std::strtod(string.c_str(), &endptr);

  if (*endptr != 0)
  {
    if (std::strncmp(endptr, "km", 2) == 0)
      radius *= 1000;
    else if (std::strncmp(endptr, "m", 1) == 0)
      ;
    else
      cdo_abort("Float parameter >%s< contains invalid character at position %d!", string, (int) (endptr - string.c_str() + 1));
  }

  return radius;
}

double
radius_str_to_deg(std::string const &string)
{
  char *endptr = nullptr;
  auto radius = std::strtod(string.c_str(), &endptr);

  if (*endptr != 0)
  {
    if (std::strncmp(endptr, "km", 2) == 0) { radius = 360 * ((radius * 1000) / (2.0 * PlanetRadiusDefault * M_PI)); }
    else if (std::strncmp(endptr, "m", 1) == 0) { radius = 360 * (radius / (2.0 * PlanetRadiusDefault * M_PI)); }
    else if (std::strncmp(endptr, "deg", 3) == 0) { ; }
    else if (std::strncmp(endptr, "rad", 3) == 0) { radius *= RAD2DEG; }
    else
    {
      cdo_abort("Float parameter >%s< contains invalid character at position %d!", string, (int) (endptr - string.c_str() + 1));
    }
  }

  if (radius > 180.0) radius = 180.0;

  return radius;
}

static int
string_to_param(const char *const paramstr)
{
  int pnum = -1, pcat = 255, pdis = 255;
  std::sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);

  return cdiEncodeParam(pnum, pcat, pdis);
}

int
string_to_param(std::string const &paramstr)
{
  return string_to_param(paramstr.c_str());
}

std::string
param_to_string(int param)
{
  char paramstr[CDI_MAX_NAME];
  int dis, cat, num;
  cdiDecodeParam(param, &num, &cat, &dis);

  int maxlen = sizeof(paramstr);
  int len;
  if (dis == 255 && (cat == 255 || cat == 0))
    len = std::snprintf(paramstr, maxlen, "%03d", num);
  else if (dis == 255)
    len = std::snprintf(paramstr, maxlen, "%03d.%03d", num, cat);
  else
    len = std::snprintf(paramstr, maxlen, "%03d.%03d.%03d", num, cat, dis);

  if (len >= maxlen || len < 0) cdo_abort("Internal problem (%s): size of input string is too small!", __func__);

  return std::string{ paramstr };
}

/* time/date/season converisons */
/* =================================================================================== */
void
season_to_months(std::string const &season, int (&imonths)[13])
{
  std::string smons{ "JFMAMJJASONDJFMAMJJASOND" };
  const int imons[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
  assert(smons.size() == (sizeof(imons) / sizeof(int)));

  if (season == "ANN")
  {
    for (int k = 1; k < 13; ++k) imonths[k] = 1;
  }
  else
  {
    auto len = (int) season.size();
    if (len > 12) cdo_abort("Too many months %d (limit=12)!", len);
    auto pos = smons.find(season);
    if (pos == std::string::npos) { cdo_abort("Season %s not available!", season); }
    auto ks = (int) pos;
    auto ke = ks + len;
    for (int k = ks; k < ke; ++k) imonths[imons[k]]++;
  }
}

double
datestr_to_double(std::string const &datestr, int opt)
{
  for (size_t i = 0, n = datestr.size(); i < n; ++i)
  {
    int c = datestr[i];
    if (!(std::isdigit(c) || c == '-' || c == ':' || c == '.' || c == 'T'))
      cdo_abort("Date string >%s< contains invalid character at position %zu!", datestr, i + 1);
  }

  if (not string_contains(datestr, '-')) { return parameter_to_double(datestr); }
  if (datestr[0] == '-' && not string_contains(datestr.c_str() + 1, '-')) { return parameter_to_double(datestr); }

  int year = 1;
  unsigned month = 1, day = 1, hour = 0, minute = 0, second = 0;
  if (opt)
  {
    hour = 23;
    minute = 59;
    second = 59;
  }

  if (string_contains(datestr, 'T'))
  {
    auto status = std::sscanf(datestr.c_str(), "%d-%u-%uT%u:%u:%u", &year, &month, &day, &hour, &minute, &second);
    if (status != 6) cdo_abort("Invalid date string >%s<!", datestr);
  }
  else
  {
    auto status = std::sscanf(datestr.c_str(), "%d-%u-%u", &year, &month, &day);
    if (status != 3) cdo_abort("Invalid date string >%s<!", datestr);
  }

  if (month > 12 || day > 31) cdo_abort("Date components out of bounds >%s<!", datestr);
  if (hour > 24 || minute > 60 || second > 60) cdo_abort("Time components out of bounds >%s<!", datestr);

  double fval = cdiEncodeTime(hour, minute, second);
  if (std::fabs(fval) > 0.0) fval /= 1000000;
  fval += std::abs(cdiEncodeDate(year, month, day));
  if (year < 0) fval *= -1;

  return fval;
}

// argv conversions ==============================================
std::vector<int>
cdo_argv_to_intarr(std::vector<std::string> const &argv)
{
  std::vector<int> v;

  for (auto const &argument : argv)
  {
    int first, last, inc;
    split_intstring(argument, first, last, inc);

    if (inc >= 0)
    {
      for (auto ival = first; ival <= last; ival += inc) v.push_back(ival);
    }
    else
    {
      for (auto ival = first; ival >= last; ival += inc) v.push_back(ival);
    }
  }

  return v;
}

std::vector<double>
cdo_argv_to_fltarr(std::vector<std::string> const &argv)
{
  std::vector<double> v;

  for (auto const &argument : argv)
  {
    auto len = (int) argument.size();
    int i;
    for (i = 0; i < len; ++i)
      if (argument[i] != '/' && argument[i] != '-' && argument[i] != '.' && !std::isdigit(argument[i])) break;

    if (i != len) { v.push_back(parameter_to_double(argument)); }
    else
    {
      double first, last, inc;
      split_fltstring(argument, first, last, inc);

      if (inc >= 0)
      {
        for (double fval = first; fval <= last; fval += inc) v.push_back(fval);
      }
      else
      {
        for (double fval = first; fval >= last; fval += inc) v.push_back(fval);
      }
    }
  }

  return v;
}

static const char *
skip_word(const char *inStr, const char *word)
{
  auto len = std::strlen(word);
  if (len <= std::strlen(inStr) && std::strncmp(word, inStr, len) == 0) { inStr += len; }
  return inStr;
}

static void
split_intstring(const char *const intstr, int &first, int &last, int &inc)
{
  auto startptr = intstr;
  char *endptr = nullptr;
  auto ival = (int) std::strtol(startptr, &endptr, 10);
  if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
  first = ival;
  last = ival;
  inc = 1;

  if (*endptr == '/')
  {
    startptr = skip_word(endptr + 1, "to/");
    endptr = nullptr;
    ival = (int) std::strtol(startptr, &endptr, 10);
    if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
    last = ival;
    if (first > last) inc = -1;

    if (*endptr == '/')
    {
      startptr = skip_word(endptr + 1, "by/");
      endptr = nullptr;
      ival = (int) std::strtol(startptr, &endptr, 10);
      if (*endptr != 0 && (endptr - startptr) == 0) parameter_error("Integer", startptr, endptr);
      inc = ival;
    }
  }
}

void
split_intstring(std::string const &intstr, int &first, int &last, int &inc)
{
  split_intstring(intstr.c_str(), first, last, inc);
}

static void
split_fltstring(const char *const intstr, double &first, double &last, double &inc)
{
  auto startptr = intstr;
  char *endptr = nullptr;
  auto fval = std::strtod(startptr, &endptr);
  if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Float", startptr, endptr);
  first = fval;
  last = fval;
  inc = 1;

  if (*endptr == '/')
  {
    startptr = skip_word(endptr + 1, "to/");
    endptr = nullptr;
    fval = std::strtod(startptr, &endptr);
    if (*endptr != 0 && *endptr != '/' && (endptr - startptr) == 0) parameter_error("Float", startptr, endptr);
    last = fval;
    if (first > last) inc = -1;

    if (*endptr == '/')
    {
      startptr = skip_word(endptr + 1, "by/");
      endptr = nullptr;
      fval = std::strtod(startptr, &endptr);
      if (*endptr != 0 && (endptr - startptr) == 0) parameter_error("Float", startptr, endptr);
      inc = fval;
    }
  }
}

void
split_fltstring(std::string const &fltstr, double &first, double &last, double &inc)
{
  split_fltstring(fltstr.c_str(), first, last, inc);
}

template <>
int
convert(std::string const &str_value)
{
  return parameter_to_int(str_value.c_str());
}

template <>
double
convert(std::string const &str_value)
{
  return parameter_to_double(str_value.c_str());
}

template <>
bool
convert(std::string const &str_value)
{
  return parameter_to_bool(str_value);
}

template <>
long
convert(std::string const &str_value)
{
  return parameter_to_long(str_value.c_str());
}

template <>
size_t
convert(std::string const &str_value)
{
  return parameter_to_size_t(str_value.c_str());
}

template <>
std::string
convert(std::string const &str_value)
{
  return parameter_to_word(str_value.c_str());
}

template <>
FillMethod
convert(std::string const &str_value)
{
  auto fillMethod = string_to_fillmethod(parameter_to_word(str_value));

  if (fillMethod == FillMethod::Undefined) cdo_abort("method=%s unsupported!", fillMethod);

  return fillMethod;
}
