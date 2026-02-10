#ifndef GRID_DEFINE_H
#define GRID_DEFINE_H

// Define a de-staggered grid for U and V
int cdo_define_destagered_grid(int gridID_u_stag, int gridID_v_stag, double *destagGridOffsets);

// Define a sampled grid of another grid
int cdo_define_sample_grid(int gridID, int sampleFactor);

// Define a sub-grid of another grid
int cdo_define_subgrid_grid(int gridSrcID, int subI0, int subI1, int subJ0, int subJ1);

#endif
