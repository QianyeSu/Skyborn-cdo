// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "geometry.h"
#include "utils_core.h"
#include "yac_types.h"
#include "ensure_array_size.h"

/** \file bnd_circle.c
 *  \brief Set of functions to calculate a bounding circle around a certain set of points.
 *
 **/

static inline double get_sin_vector_angle(
  double a[3], double b[3]) {

  double cross_ab[3];
  crossproduct_kahan(a, b, cross_ab);

  return sqrt(cross_ab[0]*cross_ab[0] +
              cross_ab[1]*cross_ab[1] +
              cross_ab[2]*cross_ab[2]);
}

// computes the bounding circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_bounding_circle_unstruct_triangle(
  double a[3], double b[3], double c[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0];
  middle_point[1] = a[1] + b[1] + c[1];
  middle_point[2] = a[2] + b[2] + c[2];

  normalise_vector(middle_point);

  double cos_angles[3] = {middle_point[0] * a[0] +
                          middle_point[1] * a[1] +
                          middle_point[2] * a[2],
                          middle_point[0] * b[0] +
                          middle_point[1] * b[1] +
                          middle_point[2] * b[2],
                          middle_point[0] * c[0] +
                          middle_point[1] * c[1] +
                          middle_point[2] * c[2]};

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  if (cos_angles[0] < cos_angles[1]) {
    if (cos_angles[0] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, a), cos_angles[0]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  } else {
    if (cos_angles[1] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, b), cos_angles[1]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes the bounding circle of a quad on the sphere
// it is assumed that all edges are great circles
static void yac_get_cell_bounding_circle_unstruct_quad(
  double a[3], double b[3], double c[3], double d[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0] + d[0];
  middle_point[1] = a[1] + b[1] + c[1] + d[1];
  middle_point[2] = a[2] + b[2] + c[2] + d[2];

  normalise_vector(middle_point);

  double cos_angle_a = middle_point[0] * a[0] +
                       middle_point[1] * a[1] +
                       middle_point[2] * a[2];
  double cos_angle_b = middle_point[0] * b[0] +
                       middle_point[1] * b[1] +
                       middle_point[2] * b[2];
  double cos_angle_c = middle_point[0] * c[0] +
                       middle_point[1] * c[1] +
                       middle_point[2] * c[2];
  double cos_angle_d = middle_point[0] * d[0] +
                       middle_point[1] * d[1] +
                       middle_point[2] * d[2];

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  double * temp_x, * temp_y;
  double temp_cos_angle_x, temp_cos_angle_y;

  if (cos_angle_a < cos_angle_b) {
    temp_x = a;
    temp_cos_angle_x = cos_angle_a;
  } else {
    temp_x = b;
    temp_cos_angle_x = cos_angle_b;
  }

  if (cos_angle_c < cos_angle_d) {
    temp_y = c;
    temp_cos_angle_y = cos_angle_c;
  } else {
    temp_y = d;
    temp_cos_angle_y = cos_angle_d;
  }

  if (temp_cos_angle_x < temp_cos_angle_y) {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_x), temp_cos_angle_x);
  } else {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_y), temp_cos_angle_y);
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes the bounding circle of a pentagon on the sphere
// it is assumed that all edges are great circles
static void yac_get_cell_bounding_circle_unstruct_penta(
  double a[3], double b[3], double c[3], double d[3], double e[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0] + d[0] + e[0];
  middle_point[1] = a[1] + b[1] + c[1] + d[1] + e[1];
  middle_point[2] = a[2] + b[2] + c[2] + d[2] + e[2];

  normalise_vector(middle_point);

  double cos_angle_a = middle_point[0] * a[0] +
                       middle_point[1] * a[1] +
                       middle_point[2] * a[2];
  double cos_angle_b = middle_point[0] * b[0] +
                       middle_point[1] * b[1] +
                       middle_point[2] * b[2];
  double cos_angle_c = middle_point[0] * c[0] +
                       middle_point[1] * c[1] +
                       middle_point[2] * c[2];
  double cos_angle_d = middle_point[0] * d[0] +
                       middle_point[1] * d[1] +
                       middle_point[2] * d[2];
  double cos_angle_e = middle_point[0] * e[0] +
                       middle_point[1] * e[1] +
                       middle_point[2] * e[2];

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  double * temp_x, * temp_y;
  double temp_cos_angle_x, temp_cos_angle_y;

  if (cos_angle_a < cos_angle_b) {
    temp_x = a;
    temp_cos_angle_x = cos_angle_a;
  } else {
    temp_x = b;
    temp_cos_angle_x = cos_angle_b;
  }

  if (cos_angle_c < cos_angle_d) {
    temp_y = c;
    temp_cos_angle_y = cos_angle_c;
  } else {
    temp_y = d;
    temp_cos_angle_y = cos_angle_d;
  }

  if (cos_angle_e < temp_cos_angle_y) {
    temp_y = e;
    temp_cos_angle_y = cos_angle_e;
  }

  if (temp_cos_angle_x < temp_cos_angle_y) {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_x), temp_cos_angle_x);
  } else {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_y), temp_cos_angle_y);
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes the bounding circle of a hexagon on the sphere
// it is assumed that all edges are great circles
static void yac_get_cell_bounding_circle_unstruct_hexa(
  double a[3], double b[3], double c[3], double d[3], double e[3], double f[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0] + d[0] + e[0] + f[0];
  middle_point[1] = a[1] + b[1] + c[1] + d[1] + e[1] + f[1];
  middle_point[2] = a[2] + b[2] + c[2] + d[2] + e[2] + f[2];

  normalise_vector(middle_point);

  double cos_angle_a = middle_point[0] * a[0] +
                       middle_point[1] * a[1] +
                       middle_point[2] * a[2];
  double cos_angle_b = middle_point[0] * b[0] +
                       middle_point[1] * b[1] +
                       middle_point[2] * b[2];
  double cos_angle_c = middle_point[0] * c[0] +
                       middle_point[1] * c[1] +
                       middle_point[2] * c[2];
  double cos_angle_d = middle_point[0] * d[0] +
                       middle_point[1] * d[1] +
                       middle_point[2] * d[2];
  double cos_angle_e = middle_point[0] * e[0] +
                       middle_point[1] * e[1] +
                       middle_point[2] * e[2];
  double cos_angle_f = middle_point[0] * f[0] +
                       middle_point[1] * f[1] +
                       middle_point[2] * f[2];

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  double * temp_x, * temp_y;
  double temp_cos_angle_x, temp_cos_angle_y;

  if (cos_angle_a < cos_angle_b) {
    temp_x = a;
    temp_cos_angle_x = cos_angle_a;
  } else {
    temp_x = b;
    temp_cos_angle_x = cos_angle_b;
  }

  if (cos_angle_c < cos_angle_d) {
    temp_y = c;
    temp_cos_angle_y = cos_angle_c;
  } else {
    temp_y = d;
    temp_cos_angle_y = cos_angle_d;
  }

  if (cos_angle_e < cos_angle_f) {
    if (cos_angle_e < temp_cos_angle_y) {
      temp_y = e;
      temp_cos_angle_y = cos_angle_e;
    }
  } else {
    if (cos_angle_f < temp_cos_angle_y) {
      temp_y = f;
      temp_cos_angle_y = cos_angle_f;
    }
  }

  if (temp_cos_angle_x < temp_cos_angle_y) {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_x), temp_cos_angle_x);
  } else {
    inc_angle =
      sin_cos_angle_new(
        get_sin_vector_angle(middle_point, temp_y), temp_cos_angle_y);
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes the bounding circle of a quad on the sphere
// it is assumed that the edges are circles of longitude and latitude
void yac_get_cell_bounding_circle_reg_quad(
   struct yac_grid_cell quad, struct bounding_circle * bnd_circle) {

  // compute edge middle point of the outwards pointing edge
  // (the one closest to the equator)
  int first_lat_edge_index = quad.edge_type[0] != YAC_LAT_CIRCLE_EDGE;
  int second_lat_edge_index = first_lat_edge_index + 2;
  int equator_lat_edge_start_index =
    (fabs(quad.coordinates_xyz[first_lat_edge_index][2]) <
     fabs(quad.coordinates_xyz[second_lat_edge_index][2]))?
       first_lat_edge_index:second_lat_edge_index;

  int equator_lat_edge_end_index =
    (equator_lat_edge_start_index == 3)?0:(equator_lat_edge_start_index+1);

  double edge_middle_point[3] =
    {quad.coordinates_xyz[equator_lat_edge_start_index][0] +
     quad.coordinates_xyz[equator_lat_edge_end_index][0],
     quad.coordinates_xyz[equator_lat_edge_start_index][1] +
     quad.coordinates_xyz[equator_lat_edge_end_index][1],
     quad.coordinates_xyz[equator_lat_edge_start_index][2]};
  double norm =
    sqrt((1.0 - edge_middle_point[2] * edge_middle_point[2]) /
         (edge_middle_point[0] * edge_middle_point[0] +
          edge_middle_point[1] * edge_middle_point[1]));
  edge_middle_point[0] *= norm;
  edge_middle_point[1] *= norm;

  yac_get_cell_bounding_circle_unstruct_penta(
    quad.coordinates_xyz[0], quad.coordinates_xyz[1],
    quad.coordinates_xyz[2], quad.coordinates_xyz[3],
    edge_middle_point, bnd_circle);
}

// computes the bounding circle of a polygon on the sphere
// it is assumed that all edges are great circles
static void yac_get_cell_bounding_circle_unstruct(
  size_t num_corners, double (* restrict coordinates_xyz)[3],
  struct bounding_circle * bnd_circle) {

  if (num_corners == 0) {

    bnd_circle->base_vector[0] = 1.0;
    bnd_circle->base_vector[1] = 0.0;
    bnd_circle->base_vector[2] = 0.0;
    bnd_circle->inc_angle = SIN_COS_ZERO;
    bnd_circle->sq_crd = DBL_MAX;
    return;
  }

  double middle_point[3] = {0.0, 0.0, 0.0};

  for (size_t i = 0; i < num_corners; ++i) {
    middle_point[0] += coordinates_xyz[i][0];
    middle_point[1] += coordinates_xyz[i][1];
    middle_point[2] += coordinates_xyz[i][2];
  }

  normalise_vector(middle_point);

  double min_cos_angle = DBL_MAX;
  double * coordinate_xyz = NULL;

  // find the biggest angle

  for (size_t i = 0; i < num_corners; ++i) {

    double * restrict temp_coordinate_xyz = coordinates_xyz[i];
    double temp_cos_angle = middle_point[0] * temp_coordinate_xyz[0] +
                            middle_point[1] * temp_coordinate_xyz[1] +
                            middle_point[2] * temp_coordinate_xyz[2];
    if (temp_cos_angle < min_cos_angle) {
      min_cos_angle = temp_cos_angle;
      coordinate_xyz = temp_coordinate_xyz;
    }
  }

  struct sin_cos_angle inc_angle =
    sin_cos_angle_new(
      get_sin_vector_angle(middle_point, coordinate_xyz), min_cos_angle);

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes a) the angle between the middle point of the edge with the middle
//             point of the polygon
//          b) half of the angle between between the two vertices of the edge
// returns the sum of both angles
static inline struct sin_cos_angle compute_edge_inc_angle(
  double * restrict a, double * restrict b, double * restrict middle_point) {

  double edge_middle_point[3] = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  normalise_vector(edge_middle_point);

  struct sin_cos_angle t1 = get_vector_angle_2(edge_middle_point, a);
  struct sin_cos_angle t2 = get_vector_angle_2(edge_middle_point, middle_point);

  return sum_angles_no_check(t1, t2);
}

void yac_get_cell_bounding_circle(struct yac_grid_cell cell,
                                  struct bounding_circle * bnd_circle) {

  int gc_edge_flag = 1;
  for (size_t i = 0; i < cell.num_corners; ++i)
    gc_edge_flag &= (cell.edge_type[i] == YAC_GREAT_CIRCLE_EDGE);

  // check for cells that only have gc edges
  if (gc_edge_flag) {

    // check for triangle
    if (cell.num_corners == 3)

      yac_get_cell_bounding_circle_unstruct_triangle(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[1],
        cell.coordinates_xyz[2], bnd_circle);

    // check for quad
    else if (cell.num_corners == 4)

      yac_get_cell_bounding_circle_unstruct_quad(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[1],
        cell.coordinates_xyz[2],
        cell.coordinates_xyz[3], bnd_circle);

    // check for pentagon
    else if (cell.num_corners == 5)

      yac_get_cell_bounding_circle_unstruct_penta(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[1],
        cell.coordinates_xyz[2],
        cell.coordinates_xyz[3],
        cell.coordinates_xyz[4], bnd_circle);

    // check for hexagon
    else if (cell.num_corners == 6)

      yac_get_cell_bounding_circle_unstruct_hexa(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[1],
        cell.coordinates_xyz[2],
        cell.coordinates_xyz[3],
        cell.coordinates_xyz[4],
        cell.coordinates_xyz[5], bnd_circle);

    else

      yac_get_cell_bounding_circle_unstruct(
        cell.num_corners,  cell.coordinates_xyz, bnd_circle);

  // check for regular quad
  } else if ((cell.num_corners == 4) &&
      (((cell.edge_type[0] == YAC_LAT_CIRCLE_EDGE) &&
        (cell.edge_type[1] == YAC_LON_CIRCLE_EDGE) &&
        (cell.edge_type[2] == YAC_LAT_CIRCLE_EDGE) &&
        (cell.edge_type[3] == YAC_LON_CIRCLE_EDGE)) ||
        ((cell.edge_type[0] == YAC_LON_CIRCLE_EDGE) &&
        (cell.edge_type[1] == YAC_LAT_CIRCLE_EDGE) &&
        (cell.edge_type[2] == YAC_LON_CIRCLE_EDGE) &&
        (cell.edge_type[3] == YAC_LAT_CIRCLE_EDGE)))) {

    yac_get_cell_bounding_circle_reg_quad(cell, bnd_circle);

  // general case application to all cell types
  } else {

    double middle_point[3];

    middle_point[0] = 0.0;
    middle_point[1] = 0.0;
    middle_point[2] = 0.0;

    size_t num_corners = cell.num_corners;

    // compute the coordinates in rad and 3d
    for (size_t i = 0; i < num_corners; ++i) {

      middle_point[0] += cell.coordinates_xyz[i][0];
      middle_point[1] += cell.coordinates_xyz[i][1];
      middle_point[2] += cell.coordinates_xyz[i][2];
    }

    normalise_vector(middle_point);

    // compute the angle required for the bounding circle
    struct sin_cos_angle max_angle =
      compute_edge_inc_angle(
        cell.coordinates_xyz[num_corners-1],
        cell.coordinates_xyz[0], middle_point);

    for (size_t i = 0; i < num_corners-1; ++i) {

      struct sin_cos_angle curr_angle =
        compute_edge_inc_angle(
          cell.coordinates_xyz[i], cell.coordinates_xyz[i+1], middle_point);

      if (compare_angles(max_angle, curr_angle) < 0) max_angle = curr_angle;
    }

    bnd_circle->base_vector[0] = middle_point[0];
    bnd_circle->base_vector[1] = middle_point[1];
    bnd_circle->base_vector[2] = middle_point[2];

    bnd_circle->inc_angle = sum_angles_no_check(max_angle, SIN_COS_TOL);
    bnd_circle->sq_crd = DBL_MAX;
  }
}

int yac_extents_overlap(struct bounding_circle * extent_a,
                        struct bounding_circle * extent_b) {

  struct sin_cos_angle base_vector_angle =
    get_vector_angle_2(extent_a->base_vector, extent_b->base_vector);

  struct sin_cos_angle tmp_angle, inc_angle_sum;
  int big_sum =
    sum_angles(extent_a->inc_angle, extent_b->inc_angle, &tmp_angle);

  if (big_sum ||
      sum_angles(tmp_angle, SIN_COS_TOL, &inc_angle_sum))
    return 1;

  return compare_angles(base_vector_angle, inc_angle_sum) <= 0;
}
