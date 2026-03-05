// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef AREA_H
#define AREA_H

#include "yac_types.h"
#include "clipping.h"

/** \file area.h
  * \brief Structs and interfaces for area calculations
  **/

/** an area of 20m x 20m on the Earth Surface is equivalent to an area on the
  * unit sphere:
 **/
#define YAC_AREA_TOL ((0.02 * 0.02) / (6371.2290 * 6371.2290))

/**
  * \brief Area calculation of a spherical cell
  *
  * The cell in split up into triangles that all have one corner in common,
  * then the area for each of the triangles is computed and summed up to build
  * the area of the cell.
  */
double yac_grid_cell_area(struct yac_grid_cell cell);
double yac_grid_cell_area_info(
  struct yac_grid_cell cell, double * barycenter, double sign);

#endif // AREA_H

