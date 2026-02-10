/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "convert_units.h"
#include "cdo_output.h"
#include "cdo_options.h"

#ifdef HAVE_UDUNITS2

#include <mutex>
static std::mutex udunitsMutex;
#define UDUNITS_LOCK() std::scoped_lock lock(udunitsMutex)

static ut_system *ut_read = nullptr;

static void *
get_converter(char *src_unit_str, char *tgt_unit_str, int *rstatus)
{
  cv_converter *ut_units_converter = nullptr;
  int status;

  *rstatus = -1;

  if (ut_read == nullptr)
    {
      ut_set_error_message_handler(ut_ignore);

      errno = 0;
      ut_read = ut_read_xml(nullptr);
      status = ut_get_status();
      if (status == UT_PARSE)
        {
          if (Options::cdoVerbose) cdo_warning("Udunits: Couldn't parse unit database!");
        }
      if (status == UT_OPEN_ENV || status == UT_OPEN_DEFAULT || status == UT_OS)
        {
          if (Options::cdoVerbose) cdo_warning("Udunits: %s", strerror(errno));
        }
      errno = 0;
      if (status != UT_SUCCESS)
        {
          if (Options::cdoVerbose) cdo_warning("Udunits: Error reading units system!");
          return nullptr;
        }
    }

  ut_trim(src_unit_str, UT_ASCII);
  ut_unit *src_unit = ut_parse(ut_read, src_unit_str, UT_ASCII);
  if (ut_get_status() != UT_SUCCESS)
    {
      if (Options::cdoVerbose) cdo_warning("Udunits: Error parsing units: [%s]", src_unit_str);
      return nullptr;
    }

  ut_trim(tgt_unit_str, UT_ASCII);
  ut_unit *tgt_unit = ut_parse(ut_read, tgt_unit_str, UT_ASCII);
  if (ut_get_status() != UT_SUCCESS)
    {
      if (Options::cdoVerbose) cdo_warning("Udunits: Error parsing units: [%s]", tgt_unit_str);
      return nullptr;
    }

  status = ut_compare(src_unit, tgt_unit);
  if (status == 0) *rstatus = -2;

  if (*rstatus == -1)
    {
      status = ut_are_convertible(src_unit, tgt_unit);
      if (status == 0) *rstatus = -3;
    }

  if (*rstatus == -1)
    {
      ut_units_converter = ut_get_converter(src_unit, tgt_unit);
      if (ut_units_converter == nullptr || ut_get_status() != UT_SUCCESS)
        {
          if (Options::cdoVerbose) cdo_warning("Udunits: Error getting converter from [%s] to [%s]", src_unit_str, tgt_unit_str);
        }
      else
        *rstatus = 0;
    }

  ut_free(src_unit);
  if (ut_get_status() != UT_SUCCESS)
    {
      if (Options::cdoVerbose) cdo_warning("Udunits: Error freeing units [%s]", src_unit_str);
      return nullptr;
    }

  ut_free(tgt_unit);
  if (ut_get_status() != UT_SUCCESS)
    {
      if (Options::cdoVerbose) cdo_warning("Udunits: Error freeing units [%s]", tgt_unit_str);
      return nullptr;
    }

  return (void *) ut_units_converter;
}

namespace cdo
{
void
convert_free(void *ut_converter)
{
  UDUNITS_LOCK();
  cv_free((cv_converter *) ut_converter);
}

void
convert_destroy()
{
  if (ut_read)
    {
      UDUNITS_LOCK();
      ut_free_system(ut_read);
      ut_read = nullptr;
    }
}
}  // namespace cdo
#endif

namespace cdo
{
void
convert_units(void **ut_converter, bool *changeunits, char *units, char *unitsOld, std::string const &name)
{
  (void) ut_converter;  // removes wrong warning, caused by ifdef

  if (*changeunits)
    {
#ifdef HAVE_UDUNITS2
      int status;
      {
        UDUNITS_LOCK();
        *ut_converter = get_converter(unitsOld, units, &status);
      }
      if (*ut_converter == nullptr)
        {
          if (status == -2)
            {
              if (Options::cdoVerbose) cdo_print("%s - not converted from  [%s] to [%s], units are equal!", name, unitsOld, units);
            }
          else if (status == -3)
            {
              cdo_warning("%s - converting units from [%s] to [%s] failed, not convertible!", name, unitsOld, units);
            }
          else { cdo_warning("%s - converting units from [%s] to [%s] failed!", name, unitsOld, units); }
          *changeunits = false;
        }
      else
        {
          // if ( Options::cdoVerbose )
          {
            char buf[64];
            cv_get_expression((const cv_converter *) *ut_converter, buf, sizeof(buf), name.c_str());
            cdo_print("%s - convert units from [%s] to [%s] (expression: %s).", name, unitsOld, units, buf);
          }
        }
#else
      static auto printWarning = true;
      if (printWarning)
        {
          cdo_warning("%s - converting units from [%s] to [%s] failed, UDUNITS2 support not compiled in!", name, unitsOld, units);
          *changeunits = false;
          printWarning = false;
        }
#endif
    }
}
}  // namespace cdo
