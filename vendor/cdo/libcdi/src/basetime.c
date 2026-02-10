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

#include <stddef.h>  // for NULL
#include "cdi.h"
#include "basetime.h"
#include "error.h"

void
basetimeInit(basetime_t *basetime)
{
  if (basetime == NULL) Error("Internal problem! Basetime not allocated.");

  if (basetime)
  {
    basetime->ncvarid = CDI_UNDEFID;
    basetime->ncdimid = CDI_UNDEFID;
    basetime->ncvarboundsid = CDI_UNDEFID;
    basetime->leadtimeid = CDI_UNDEFID;
    basetime->hasUnits = false;
    basetime->isWRF = false;
  }
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
