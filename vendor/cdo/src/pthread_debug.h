/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PTHREAD_DEBUG_H
#define PTHREAD_DEBUG_H

#include <pthread.h>

void print_pthread_attr(const char *caller, pthread_attr_t *attr);

int Pthread_create(const char *caller, pthread_t *th, pthread_attr_t *attr, void *(*start_routine)(void *), void *arg);
int Pthread_join(const char *caller, pthread_t th, void **thread_return);

#define pthread_create(a, b, c, d) Pthread_create(__func__, a, b, c, d)
#define pthread_join(a, b) Pthread_join(__func__, a, b)

#endif /* PTHREAD_DEBUG_H */
