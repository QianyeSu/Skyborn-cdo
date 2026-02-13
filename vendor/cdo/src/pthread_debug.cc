/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>

#ifdef HAVE_LIBPTHREAD

#include <pthread.h>
#include <errno.h>
#include "cdo_output.h"

#define POUT2(caller, x, a, b) pout2(caller, #x, x, #a, a, #b, b)
#define POUT3(caller, x, a, b, c) pout3(caller, #x, x, #a, a, #b, b, #c, c)

static void
pout2(const char *caller, const char *sval, int ival, const char *sval1, int oval1, const char *sval2, int oval2)
{
  if (ival == oval1)
    std::fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if (ival == oval2)
    std::fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else
    std::fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}

static void
pout3(const char *caller, const char *sval, int ival, const char *sval1, int oval1, const char *sval2, int oval2, const char *sval3,
      int oval3)
{
  if (ival == oval1)
    std::fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if (ival == oval2)
    std::fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else if (ival == oval3)
    std::fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval3);
  else
    std::fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}

void
print_pthread_attr(const char *caller, pthread_attr_t *attr)
{
  struct sched_param param;
  int detachstate, policy, inherit, scope, priority;
  size_t stacksize;

  pthread_attr_getdetachstate(attr, &detachstate);
  POUT2(caller, detachstate, PTHREAD_CREATE_JOINABLE, PTHREAD_CREATE_DETACHED);

#ifdef SCHED_FIFO
  pthread_attr_getschedpolicy(attr, &policy);
  POUT3(caller, policy, SCHED_FIFO, SCHED_RR, SCHED_OTHER);
  pthread_attr_getschedparam(attr, &param);
  priority = param.sched_priority;
  std::fprintf(stderr, "%-18s :  %-14s = %d\n", caller, "priority", priority);
#endif

#ifdef PTHREAD_INHERIT_SCHED
  pthread_attr_getinheritsched(attr, &inherit);
  POUT2(caller, inherit, PTHREAD_INHERIT_SCHED, PTHREAD_EXPLICIT_SCHED);
#endif

  pthread_attr_getscope(attr, &scope);
  POUT2(caller, scope, PTHREAD_SCOPE_SYSTEM, PTHREAD_SCOPE_PROCESS);

  pthread_attr_getstacksize(attr, &stacksize);
  std::fprintf(stderr, "%-18s :  %-14s = %ld\n", caller, "stacksize", (long) stacksize);
}

int
Pthread_create(const char *caller, pthread_t *th, pthread_attr_t *attr, void *(*start_routine)(void *), void *arg)
{
  int status;

  Debug(PTHREAD, "%s", caller);

  if (PTHREAD)
  {
    cdo_print("%s attributes:", caller);
    if (attr)
      print_pthread_attr(__func__, attr);
    else
      cdo_print("  default attributes");
  }

  status = pthread_create(th, attr, start_routine, arg);

  // Debug(PTHREAD,"-%s (thID = %ld, status = %d)", caller, (long) *th, status);

  return status;
}

int
Pthread_join(const char *caller, pthread_t th, void **thread_return)
{
  (void) caller;

  // Debug(PTHREAD,"+%s (thID = %ld)", caller, (void *) th);

  int status = pthread_join(th, thread_return);

  // Debug(PTHREAD,"-%s (thID = %ld, status = %d)", caller, (void *) th, status);

  return status;
}

#endif
