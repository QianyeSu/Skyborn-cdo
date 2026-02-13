#ifndef CDF_CONFIG_H
#define CDF_CONFIG_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include <netcdf.h>
#include <netcdf_meta.h>

#ifdef NC_FORMAT_64BIT_DATA
#define HAVE_NETCDF5 1
#endif

#if defined NC_HAS_NC2 && NC_HAS_NC2
#define HAVE_NETCDF2 1
#endif

#if defined NC_HAS_NCZARR && NC_HAS_NCZARR
#define HAVE_NCZARR 1
#endif

#if ((NC_VERSION_MAJOR == 4 && NC_VERSION_MINOR >= 9) || NC_VERSION_MAJOR >= 5)
#define HAVE_NC4FILTER 1
#endif

#endif  // HAVE_LIBNETCDF

#endif  // CDF_CONFIG_H
