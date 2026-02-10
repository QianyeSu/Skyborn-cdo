#ifndef CONVERT_UNITS_H
#define CONVERT_UNITS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>

#if defined(HAVE_LIBUDUNITS2) && (defined(HAVE_UDUNITS2_H) || defined(HAVE_UDUNITS2_UDUNITS2_H))
#define HAVE_UDUNITS2
#endif

#ifdef HAVE_UDUNITS2
#ifdef HAVE_UDUNITS2_UDUNITS2_H
#include <udunits2/udunits2.h>
#else
#include <udunits2.h>
#endif

namespace cdo
{
void convert_free(void *ut_converter);
void convert_destroy();
}  // namespace cdo
#endif

namespace cdo
{
void convert_units(void **ut_converter, bool *changeunits, char *units, char *units_old, std::string const &name);
}
#endif  // CONVERT_UNITS_H
