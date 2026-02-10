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

#include <inttypes.h>
#include <sys/types.h>
#include <stdlib.h>

#include "cdi_cksum.h"
#include "cksum.h"
#include "error.h"
#include "serialize.h"

uint32_t
cdiCheckSum(int type, int count, const void *buffer)
{
  uint32_t s = 0U;
  xassert(count >= 0);
  size_t elemSize = (size_t) serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&s, (const unsigned char *) buffer, (size_t) count, elemSize);
  s = memcrc_finish(&s, (off_t) (elemSize * (size_t) count));
  return s;
}

void
cdiCheckSumRStart(struct cdiCheckSumState *state)
{
  state->sum = 0U;
  state->len = 0;
}

void
cdiCheckSumRAdd(struct cdiCheckSumState *state, int type, int count, const void *buffer)
{
  size_t elemSize = (size_t) serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&state->sum, (const unsigned char *) buffer, (size_t) count, elemSize);
  state->len += (off_t) (elemSize * (size_t) count);
}

uint32_t
cdiCheckSumRValue(struct cdiCheckSumState state)
{
  return memcrc_finish(&state.sum, state.len);
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
