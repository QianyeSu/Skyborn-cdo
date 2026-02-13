#ifndef GRID_READ_PINGO_H
#define GRID_READ_PINGO_H

#include <cstdio>
#include "varray.h"

size_t input_darray(std::FILE *gfp, size_t n_values, Varray<double> &array);
int grid_read_pingo(std::FILE *gfp);

#endif
