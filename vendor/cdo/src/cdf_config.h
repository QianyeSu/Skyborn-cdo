#ifndef CDF_CONFIG_H
#define CDF_CONFIG_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include <netcdf_meta.h>

#if defined NC_HAS_HDF5 && NC_HAS_HDF5
#define HAVE_NC4HDF5 1
#endif

#if defined NC_HAS_S3_INTERNAL && NC_HAS_S3_INTERNAL
#define HAVE_NC4S3 1
#endif

#endif  // HAVE_LIBNETCDF

#endif  // CDF_CONFIG_H
