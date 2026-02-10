#ifndef CDO_FILL_H
#define CDO_FILL_H

#include <cstdio>  // size_t
#include "varray.h"

void cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData);
void cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData, Varray2D<size_t> &p_varNmiss);

#endif
