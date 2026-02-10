/* DO NOT REMOVE the config.h include file under any circumstances,
 * it's very much needed on some platforms */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
/* DO NOT REMOVE the above config.h include file under any
 * circumstances as long as it's the autoconf configuration header
 * used to build this package. When it's missing on some platforms,
 * some poor person has to do long, tedious debugging sessions, where
 * struct offsets almost imperceptibly change from one file to the
 * next to find out what happened */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "cdi.h"
#include "cdf_int.h"

const char *
cdiStringError(int cdiErrno)
{
  // clang-format off
  static const char UnknownError[] = "Unknown Error";
  static const char _ETMOF[]       = "Too many open files";
  static const char _EINVAL[]      = "Invalid argument";
  static const char _EISDIR[]      = "Is a directory";
  static const char _EISEMPTY[]    = "File is empty";
  static const char _EUFTYPE[]     = "Unsupported file type";
  static const char _ELIBNAVAIL[]  = "Unsupported file type (library support not compiled in)";
  static const char _EUFSTRUCT[]   = "Unsupported file structure";
  static const char _EUNC4[]       = "Unsupported NetCDF4 structure";
  static const char _EDIMSIZE[]    = "Invalid dimension size";
  static const char _EQENF[]       = "Query entries not found";
  static const char _EQNAVAIL[]    = "Query not available for file type";
  static const char _ELIMIT[]      = "Internal limits exceeded";

  if (cdiErrno < -1000)
    {
      const char *cdfErrorString = cdf_strerror(cdiErrno + 1000);
      if (cdfErrorString) return cdfErrorString;
    }

  switch (cdiErrno) {
  case CDI_ESYSTEM:
    {
      const char *cp = strerror(errno);
      if (cp == NULL) break;
      return cp;
    }
  case CDI_ETMOF:      return _ETMOF;
  case CDI_EINVAL:     return _EINVAL;
  case CDI_EISDIR:     return _EISDIR;
  case CDI_EISEMPTY:   return _EISEMPTY;
  case CDI_EUFTYPE:    return _EUFTYPE;
  case CDI_ELIBNAVAIL: return _ELIBNAVAIL;
  case CDI_EUFSTRUCT:  return _EUFSTRUCT;
  case CDI_EUNC4:      return _EUNC4;
  case CDI_EDIMSIZE:   return _EDIMSIZE;
  case CDI_EQENF:      return _EQENF;
  case CDI_EQNAVAIL:   return _EQNAVAIL;
  case CDI_ELIMIT:     return _ELIMIT;
  }
  // clang-format on
  return UnknownError;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
