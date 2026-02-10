#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdf_config.h"
#include "cdf_filter.h"
#include "error.h"

#ifdef HAVE_NC4FILTER
#include <netcdf.h>
#include <netcdf_filter.h>
#include <netcdf_aux.h>

#include <string.h>

static void
check_length(size_t maxLength, size_t len)
{
  if (len >= maxLength) Error("Internal error: size of filterSpec to small!");
}
#endif

bool
cdf_get_var_filter(int ncid, int varid, char *filterSpec, size_t maxLength)
{
  bool hasFilter = false;
#ifdef HAVE_NC4FILTER
  size_t numFilters = 0;
  nc_inq_var_filter_ids(ncid, varid, &numFilters, NULL);
  if (numFilters > 0 && numFilters < 16)
  {
    unsigned int filterids[16];
    nc_inq_var_filter_ids(ncid, varid, &numFilters, filterids);
    for (size_t k = 0; k < numFilters; ++k)
    {
      if (k > 0)
      {
        size_t len = strlen(filterSpec);
        check_length(maxLength, len);
        strncat(filterSpec, "|", maxLength - len - 1);
      }
      unsigned filterId = filterids[k];
      {
        size_t len = strlen(filterSpec);
        check_length(maxLength, len);
        snprintf(filterSpec + len, maxLength - len - 1, "%u", filterId);
      }
      size_t numParams;
      nc_inq_var_filter_info(ncid, varid, filterId, &numParams, NULL);
      if (numParams <= 16)
      {
        unsigned int params[16];
        nc_inq_var_filter_info(ncid, varid, filterId, &numParams, params);
        for (size_t i = 0; i < numParams; ++i)
        {
          size_t len = strlen(filterSpec);
          check_length(maxLength, len);
          snprintf(filterSpec + len, maxLength - len - 1, ",%u", params[i]);
        }
      }
    }
    if (filterSpec[0]) hasFilter = true;
  }
#else
  (void) ncid;
  (void) varid;
  (void) filterSpec;
  (void) maxLength;
#endif
  return hasFilter;
}

void
cdf_def_var_filter(int ncid, int ncvarID, const char *filterSpec)
{
  if (filterSpec)
  {
#ifdef HAVE_NC4FILTER
    size_t nfilters = 0;
    NC_H5_Filterspec **filters = NULL;
    int status = ncaux_h5filterspec_parselist(filterSpec, NULL, &nfilters, &filters);
    if (status != NC_NOERR)
    {
      Message("filterSpec=%s", filterSpec);
      Error("ncaux_h5filterspec_parselist failed: %s", nc_strerror(status));
    }

    if (filters != NULL)
    {
      for (size_t i = 0; i < nfilters; i++)
      {
        unsigned int filterid = filters[i]->filterid;
        // printf("filter %zu id:%d nparams:%zu param1 %d\n", i + 1, filterid, filters[i]->nparams, filters[i]->params[0]);
        status = nc_def_var_filter(ncid, ncvarID, filterid, filters[i]->nparams, filters[i]->params);
        if (status != NC_NOERR)
        {
          Message("filterid=%u  numParams=%zu", filterid, filters[i]->nparams);
          Error("nc_def_var_filter failed: %s", nc_strerror(status));
        }
      }

      for (size_t i = 0; i < nfilters; i++) ncaux_h5filterspec_free(filters[i]);
      free(filters);
    }
#else
    (void) ncid;
    (void) ncvarID;
    (void) filterSpec;
    Error("Filter failed, NetCDF4 function ncaux_h5filterspec_parselist() not available!");
#endif
  }
}
