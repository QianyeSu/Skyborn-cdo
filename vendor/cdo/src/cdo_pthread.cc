/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include "pthread_debug.h"
#endif

#include <cstdio>

void
print_pthread_info()
{
#ifdef HAVE_LIBPTHREAD
  std::fprintf(stderr, "\n");

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  print_pthread_attr("Default pthread attr", &attr);
  pthread_attr_destroy(&attr);

  std::fprintf(stderr, "\n");
#endif
}
