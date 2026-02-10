#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H

#include <stddef.h>
#include "yac_types.h"

void yac_compute_weights_rbf(double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n, double *weights,
                             double const rbf_scale);

#endif
