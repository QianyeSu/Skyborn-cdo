#ifndef CDF_FILTER_H
#define CDF_FILTER_H

#include <stdio.h>

int cdf_get_var_filter(int ncid, int varid, char *filterSpec, size_t maxLength);

#endif
