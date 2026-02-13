#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdf_config.h"

#ifdef HAVE_LIBNETCDF
#include <netcdf_meta.h>
#endif

int
cdiGetConfig(int confType)
{
#if defined HAVE_LIBCGRIBEX
  static const int cdi_has_cgribex = 1;
#else
  static const int cdi_has_cgribex = 0;
#endif
#if defined HAVE_NC4FILTER && HAVE_NC4FILTER
  static const int cdi_nc_has_filter = 1;
#else
  static const int cdi_nc_has_filter = 0;
#endif
#if defined NC_HAS_DAP2 && NC_HAS_DAP2
  static const int cdi_nc_has_dap = 1;
#else
  static const int cdi_nc_has_dap = 0;
#endif
#if defined NC_HAS_S3_INTERNAL && NC_HAS_S3_INTERNAL
  static const int cdi_nc_has_s3 = 1;
#else
  static const int cdi_nc_has_s3 = 0;
#endif
#if defined NC_HAS_HDF5 && NC_HAS_HDF5
  static const int cdi_nc_has_hdf5 = 1;
#else
  static const int cdi_nc_has_hdf5 = 0;
#endif

  if (confType == CDI_HAS_CGRIBEX) return cdi_has_cgribex;
  if (confType == CDI_NC_HAS_FILTER) return cdi_nc_has_filter;
  if (confType == CDI_NC_HAS_DAP) return cdi_nc_has_dap;
  if (confType == CDI_NC_HAS_S3) return cdi_nc_has_s3;
  if (confType == CDI_NC_HAS_HDF5) return cdi_nc_has_hdf5;

  return 0;
}
