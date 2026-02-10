#ifndef CF_INTERFACE_H
#define CF_INTERFACE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_CF_INTERFACE

#ifndef __CFORTRAN_LOADED
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#endif
#include "cfortran.h"
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
#endif

#endif

#endif
