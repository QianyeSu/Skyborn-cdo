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

#include "referenceCounting.h"

#include "dmemory.h"
#include "error.h"

void
cdiRefObject_construct(CdiReferencedObject *me)
{
  me->destructor = cdiRefObject_destruct;
  me->refCount = 1;
}

void
cdiRefObject_retain(CdiReferencedObject *me)
{
  size_t oldCount = me->refCount++;
  xassert(oldCount && "A reference counted object was used after it was destructed.");
}

void
cdiRefObject_release(CdiReferencedObject *me)
{
  size_t oldCount = me->refCount--;
  xassert(oldCount && "A reference counted object was released too often.");
  if (oldCount == 1)
  {
    me->destructor(me);
    Free(me);
  }
}

void
cdiRefObject_destruct(CdiReferencedObject *me)
{
  (void) me;
  /* Empty for now, but that's no reason not to call it! */
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
