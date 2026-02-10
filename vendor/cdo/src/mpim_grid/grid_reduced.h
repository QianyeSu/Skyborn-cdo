#ifndef GRID_REDUCED_H
#define GRID_REDUCED_H

#include "stdlib.h"

extern "C" int qu2reg3_double(double *pfield, int *kpoint, int klat, int klon, double msval, int *kret, int omisng, int operio,
                              int oveggy);

void grib_get_reduced_row(long pl, double lon_first, double lon_last, long *npoints, long *ilon_first, long *ilon_last);

int qu2reg_subarea(size_t gridsize, int np, double xfirst, double xlast, double *array, int *reducedPoints, int ny, double missval,
                   int *iret, int lmiss, int lperio, int lveggy);

#endif
