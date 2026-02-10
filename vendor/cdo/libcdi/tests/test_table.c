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

#include <assert.h>
#include <stdlib.h>

#include "cdi.h"

int
main(void)
{
  tableInqEntry(-1, -1, -1, NULL, NULL, NULL);
  // assert(tableInqEntry(-1, -1, -1, NULL, NULL, NULL) != 0);
  return EXIT_SUCCESS;
}
