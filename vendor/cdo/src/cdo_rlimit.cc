/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GETRLIMIT
#if __has_include(<sys/resource.h>)
#include <sys/time.h>      // getrlimit
#include <sys/resource.h>  // getrlimit
#endif
#endif

#include <cstdio>
#include <algorithm>

#include "cdo_rlimit.h"
#include "cdo_options.h"
#include "mpmo.h"

#ifdef SX
#define RLIM_T long long
#else
#define RLIM_T rlim_t
#endif

#define PRINT_RLIMIT(resource)                                                           \
  {                                                                                      \
    struct rlimit rlim;                                                                  \
    auto status = getrlimit(resource, &rlim);                                            \
    if (status == 0)                                                                     \
      {                                                                                  \
        if (sizeof(RLIM_T) > sizeof(long))                                               \
          {                                                                              \
            std::fprintf(stderr, "CUR %-15s = %llu\n", #resource, (long long) rlim.rlim_cur); \
            std::fprintf(stderr, "MAX %-15s = %llu\n", #resource, (long long) rlim.rlim_max); \
          }                                                                              \
        else                                                                             \
          {                                                                              \
            std::fprintf(stderr, "CUR %-15s = %lu\n", #resource, (long) rlim.rlim_cur);       \
            std::fprintf(stderr, "MAX %-15s = %lu\n", #resource, (long) rlim.rlim_max);       \
          }                                                                              \
      }                                                                                  \
  }

namespace cdo
{

void
print_rlimits(void)
{
#ifdef HAVE_GETRLIMIT
#ifdef RLIMIT_FSIZE
  PRINT_RLIMIT(RLIMIT_FSIZE);
#endif
#ifdef RLIMIT_NOFILE
  PRINT_RLIMIT(RLIMIT_NOFILE);
#endif
#ifdef RLIMIT_STACK
  PRINT_RLIMIT(RLIMIT_STACK);
#endif
#ifdef RLIMIT_CORE
  PRINT_RLIMIT(RLIMIT_CORE);
#endif
#ifdef RLIMIT_RSS
  PRINT_RLIMIT(RLIMIT_RSS);
#endif
#endif
}

#ifdef HAVE_GETRLIMIT
static void
set_rlimit(struct rlimit lim, int resource, const char *rname)
{
  auto stat = setrlimit(resource, &lim);
  if (MpMO::DebugLevel > 0)
    {
      if (stat == 0)
        {
          std::fprintf(stderr, "Set %s to %ld\n", rname, (long) lim.rlim_cur);
          PRINT_RLIMIT(resource);
        }
      else { std::fprintf(stderr, "Set %s to %ld failed!\n", rname, (long) lim.rlim_cur); }
      std::fprintf(stderr, "\n");
    }
}
#endif

static void
set_rlimit_min(long rsize, int resource, const char *rname)
{
#ifdef HAVE_GETRLIMIT
  struct rlimit lim;
  auto stat = getrlimit(resource, &lim);
  if (stat == 0)
    {
      RLIM_T minSize = rsize;
      minSize = std::min(minSize, lim.rlim_max);
      if (lim.rlim_cur < minSize)
        {
          lim.rlim_cur = minSize;
          set_rlimit(lim, resource, rname);
        }
    }
#endif
}

static void
set_rlimit_max(long rsize, int resource, const char *rname)
{
#ifdef HAVE_GETRLIMIT
  struct rlimit lim;
  auto stat = getrlimit(resource, &lim);
  if (stat == 0)
    {
      RLIM_T maxSize = rsize;
      if (maxSize < lim.rlim_cur)
        {
          lim.rlim_cur = maxSize;
          set_rlimit(lim, resource, rname);
        }
    }
#endif
}

void
set_numfiles(long numfiles)
{
#if defined(HAVE_GETRLIMIT) && defined(RLIMIT_NOFILE)
  set_rlimit_min(numfiles, RLIMIT_NOFILE, "numfiles");
#endif
}

void
set_stacksize(long stacksize)
{
#if defined(HAVE_GETRLIMIT) && defined(RLIMIT_STACK)
  set_rlimit_min(stacksize, RLIMIT_STACK, "stacksize");
#endif
}

void
set_coresize(long coresize)
{
#if defined(HAVE_GETRLIMIT) && defined(RLIMIT_CORE)
  set_rlimit_max(coresize, RLIMIT_CORE, "coresize");
#endif
}

long
get_rss_cur()
{
  long rss = 0;
#ifdef HAVE_GETRLIMIT
#ifdef RLIMIT_RSS
  struct rlimit lim;
  if (getrlimit(RLIMIT_RSS, &lim) == 0) rss = (long) lim.rlim_cur;
#endif
#endif
  return rss;
}

}  // namespace cdo
