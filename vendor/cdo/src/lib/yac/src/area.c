// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "area.h"
#include "geometry.h"
#include "utils_core.h"
#include "ensure_array_size.h"

static inline double scalar_product(double a[], double b[]);

/** Compute the area of a spherical triangle using Eriksson's formula
  *
  * (Eriksson, F. (1990). On the Measure of Solid Angles. Mathematics Magazine,
  *  63(3), 184–187. https://doi.org/10.1080/0025570X.1990.11977515)
  *
  * This function computes the signed spherical excess (area) of a triangle
  * on the unit sphere using a numerically stable variant of Eriksson's formula:
  *
  *   area = 2 * atan2(a·(b×c), 1 + a·b + b·c + c·a)
  *
  * The implementation uses compensated arithmetic (error-free transformations
  * via FMA instructions) to maintain high accuracy for small triangles and
  * near-degenerate configurations.
  *
  * This formula was proposed for YAC by Hongyu Chen, who evaluated several
  * methods for computing the area of a spherical triangle for her paper:
  * "Accurate and Robust Algorithms for Conservative Regridding on the Sphere"
  * (full reference will follow). Erikksons formula significantly exceeds
  * the accuracy of L'Huilier's theorem based implementations.
  *
  * \param[in] a First vertex in Cartesian coordinates (x, y, z) on unit sphere
  * \param[in] b Second vertex in Cartesian coordinates (x, y, z) on unit sphere
  * \param[in] c Third vertex in Cartesian coordinates (x, y, z) on unit sphere
  *
  * \return Signed spherical excess (area) of the triangle. Positive if vertices
  *         are ordered counter-clockwise when viewed from outside the sphere,
  *         negative if clockwise.
  *
  * \remark All three vertices are assumed to lie on the unit sphere
  * \remark The triangle edges are great circle arcs
  */
static double tri_area(double const * restrict a,
                       double const * restrict b,
                       double const * restrict c) {

  /* --------------------------------------------------
    * Compensated cross product b x c
    * -------------------------------------------------- */

  double bc[3];   // b x c
  double bc_e[3]; // error from computation of bc
  det2_error(b[1], b[2], c[1], c[2], &bc[0], &bc_e[0]);
  det2_error(b[2], b[0], c[2], c[0], &bc[1], &bc_e[1]);
  det2_error(b[0], b[1], c[0], c[1], &bc[2], &bc_e[2]);

  /* --------------------------------------------------
    * Compensated scalar triple product
    * T = a * (b x c)
    * -------------------------------------------------- */

  double abc[3] =   // a * bc
    {a[0]*bc[0], a[1]*bc[1], a[2]*bc[2]};
  double abc_e[3] = // error from computation of abc
    {fma(a[0], bc[0], -abc[0]) + a[0] * bc_e[0],
     fma(a[1], bc[1], -abc[1]) + a[1] * bc_e[1],
     fma(a[2], bc[2], -abc[2]) + a[2] * bc_e[2]};

  double triple = 0.0, triple_e = 0.0;

  accu_eft(abc_e[0], &triple, &triple_e);
  accu_eft(abc_e[1], &triple, &triple_e);
  accu_eft(abc_e[2], &triple, &triple_e);
  accu_eft(abc[0], &triple, &triple_e);
  accu_eft(abc[1], &triple, &triple_e);
  accu_eft(abc[2], &triple, &triple_e);

  /* --------------------------------------------------
    * Compensated denominator
    * D = 1 + a·b + b·c + c·a
    * -------------------------------------------------- */

  double D = 1.0;   // denominator
  double D_e = 0.0; // error from computing D

  accu_product_eft(a[0], b[0], &D, &D_e);
  accu_product_eft(a[1], b[1], &D, &D_e);
  accu_product_eft(a[2], b[2], &D, &D_e);

  accu_product_eft(b[0], c[0], &D, &D_e);
  accu_product_eft(b[1], c[1], &D, &D_e);
  accu_product_eft(b[2], c[2], &D, &D_e);

  accu_product_eft(c[0], a[0], &D, &D_e);
  accu_product_eft(c[1], a[1], &D, &D_e);
  accu_product_eft(c[2], a[2], &D, &D_e);

  /* --------------------------------------------------
    * Signed spherical excess
    * -------------------------------------------------- */

  return 2.0 * atan2(triple + triple_e, D + D_e);
}

/* ----------------------------------- */

static inline int compute_norm_vector(double a[], double b[], double norm[]) {

  crossproduct_kahan(a, b, norm);

  double scale = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  if (scale <= yac_angle_tol) return 1;

  scale = 1.0 / scale;

  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;

  return 0;
}

static inline double XYZtoLon(double a[3]) {
  return atan2(a[1] , a[0]);
}

static double lat_edge_correction_(
  double base_vec[3], double a[3], double b[3], double * middle_lat) {

  //-------------------------------------------------------------------
  // compute the area that is enclosed between a great circle edges and
  // a lat circle edge using that go through the points a and b
  //-------------------------------------------------------------------

  double const tol = 1e-8;

  YAC_ASSERT(
    fabs(a[2]-b[2]) <= tol,
    "ERROR(lat_edge_correction_): "
    "latitude of both corners is not identical")

  double h = fabs(a[2]);

  // if we are at the equator or at a pole -> area is zero
  if (h < tol || fabs(1.0 - h) < tol) return 0.0;

  // compute the norm vector of the plane going through a and b
  // (if the angle between a and b is to small to compute a norm vector, then
  //  the area is negligible)
  double norm_ab[3];
  if (compute_norm_vector(a, b, norm_ab)) return 0.0;

  // compute the area of a triangle consisting of the lat edge and
  // the reference pole
  double lat_area = fabs((1.0 - h) * get_angle(XYZtoLon(a), XYZtoLon(b)));

  // compute the same area, but use a gc edge instead of a lat edge
  double pole[3] = {0, 0, (a[2] >= 0.0)?1.0:-1.0};
  double gc_area = fabs(tri_area(a, b, pole));

  // the difference between the two triangle areas is the area enclosed
  // between the points a and b using a lat and a gc edge
  double correction = MAX(lat_area - gc_area, 0.0);

  //-------------------------------------------------------------------
  // now we have to decide whether this correction needs to be added or
  // subtracted
  //-------------------------------------------------------------------

  // compute the middle point of the lat edge
  middle_lat[0] = a[0]+b[0];
  middle_lat[1] = a[1]+b[1];
  middle_lat[2] = a[2];
  double scale = sqrt(middle_lat[0]*middle_lat[0]+middle_lat[1]*middle_lat[1]);
  YAC_ASSERT(fabs(scale) >= 1e-18, "internal error")
  scale = sqrt(1.0 - a[2]*a[2]) / scale;
  middle_lat[0] *= scale;
  middle_lat[1] *= scale;

  // compute the cosine between norm vector and the base vector
  // --> their sign indicates the ordering of vertices
  double scalar_base = scalar_product(norm_ab, base_vec);

  int negative_correction_flag;

  // if the base point is on the same plane as a and b
  if (fabs(scalar_base) < 1e-11) {

    // compute the norm vector of the lat edge middle point and the associated
    // pole
    // (if the angle between lat edge middle point and the pols is to small to
    //  compute a norm vector, then the area is negligible)
    double norm_middle[3];
    if (compute_norm_vector(middle_lat, pole, norm_middle)) return 0.0;

    double scalar_a = scalar_product(norm_middle, a);
    negative_correction_flag = scalar_a <= 0;

  } else {

    // compute the cosine between the edge norm vector and the lat edge middle
    // point
    double scalar_middle_lat = scalar_product(norm_ab, middle_lat);
    negative_correction_flag = scalar_middle_lat >= 0;
  }

  return negative_correction_flag?-correction:correction;
}

static double lat_edge_correction(
  double base_vec[3], double a[3], double b[3]) {

  double middle_lat[3];
  return lat_edge_correction_(base_vec, a, b, middle_lat);
}

static double lat_edge_correction_info(
  double base_vec[3], double a[3], double b[3],
  double * barycenter, double sign) {

  double middle_lat[3];
  double correction = lat_edge_correction_(base_vec, a, b, middle_lat);

  if (correction == 0.0) return 0.0;

  double middle_gc[3] = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  normalise_vector(middle_gc);

  double temp_barycenter[3] =
    {middle_lat[0] + middle_gc[0],
     middle_lat[1] + middle_gc[1],
     middle_lat[2] + middle_gc[2]};
  normalise_vector(temp_barycenter);

  barycenter[0] += temp_barycenter[0] * correction * sign;
  barycenter[1] += temp_barycenter[1] * correction * sign;
  barycenter[2] += temp_barycenter[2] * correction * sign;

  return correction * sign;
}

/** Compute the area of a grid cell on the unit sphere
  *
  * This function computes the area of a spherical polygon that may contain
  * both great circle edges and latitude circle edges. The polygon is
  * triangulated using a fan triangulation from the first vertex, and the
  * signed areas of the triangles are summed to correctly handle concave cells.
  * For latitude circle edges, a correction term is added to account for the
  * difference between the great circle path and the latitude circle path.
  *
  * \param[in] cell Grid cell structure
  *
  * \return The unsigned area of the grid cell on the unit sphere
  *
  * \remark For triangles with only great circle edges, uses tri_area directly
  * \remark For polygons or cells with latitude edges, uses triangulation
  *         with corrections
  * \remark Returns 0.0 for degenerate cells with fewer than 2 corners
  * \remark The result is always non-negative (absolute value of signed area)
  */
double yac_grid_cell_area (struct yac_grid_cell cell) {

  size_t M = cell.num_corners;

  if (M < 2) return 0.0;

  // determine whether the cell contains at least one lat circle edge
  int lat_flag = 0;
  for (size_t i = 0; i < M; i++)
    lat_flag |= cell.edge_type[i] == YAC_LAT_CIRCLE_EDGE;

  // if this is a basic triangle only conisting of great circle edges
  // (e.g. ICON cell)
  if (M == 3 && !lat_flag) {

    // in case of triangle only consisting of great circle edges, return
    // absolut value of the signed area returned by tri_area
    return
      fabs(
        tri_area(cell.coordinates_xyz[0],
                 cell.coordinates_xyz[1],
                 cell.coordinates_xyz[2]));
  }

  // loop over all subtriangles generated by triangulation using
  // a fan approach
  double area = 0.0;
  for (size_t i = 1; i < M - 1; ++i) {

    // signed area summation, this ensures that this routine
    // also works for concave cells
    area += tri_area(cell.coordinates_xyz[0],
                     cell.coordinates_xyz[i],
                     cell.coordinates_xyz[i+1]);
  }

  // if there is at least one latitude circle edge
  if (lat_flag) {

    // loop over all edges of the polygon
    for (size_t i = 0, j = M - 1; i < M; j = i++) {

      if (cell.edge_type[j] == YAC_LAT_CIRCLE_EDGE) {

        // add correction for the area difference between the great circle
        // edge and the latitude circle edge from vertex j to vertex i
        // (uses vertex 0 as a reference to determine area sign)
        area += lat_edge_correction(cell.coordinates_xyz[0],
                                    cell.coordinates_xyz[j],
                                    cell.coordinates_xyz[i]);
      }
    }
  }

  return fabs(area);
}

/** Compute the signed area of a spherical triangle and update barycenter
  *
  * This function computes the signed area of a spherical triangle using
  * tri_area() and simultaneously updates a barycenter accumulator. The
  * barycenter contribution is computed using the edge normal vectors weighted
  * by half of their associated edge lengths. This approach is exact for
  * spherical triangles.
  *
  * \param[in] ref Reference (first) vertex in Cartesian coordinates on unit sphere
  * \param[in] a Second vertex in Cartesian coordinates on unit sphere
  * \param[in] b Third vertex in Cartesian coordinates on unit sphere
  * \param[in,out] barycenter Barycenter accumulator to be updated (3D vector)
  * \param[in] sign Sign multiplier for both area and barycenter contribution
  *
  * \return Signed area of the triangle multiplied by sign parameter
  *
  * \remark The barycenter is accumulated but not normalized by this function
  * \remark Short edges (sin_angle < yac_angle_tol) are skipped
  * \remark All vertices are assumed to lie on the unit sphere
  */
static double tri_area_info(
  double ref[3], double a[3], double b[3],
  double * barycenter, double sign) {

  // the barycenter of the triangle is given by the sum edge norm vector
  // scaled by half of the associated edge length
  // This approach is exact for spherical triangles. The weighted normalized
  // sum of the triangle vertices works for the planar triangle.
  double * corners[3] = {b, ref, a};
  for (int i = 0, j = 2; i < 3; j = i++) {

    double cross[3];
    crossproduct_kahan(corners[j], corners[i], cross);
    double sin_angle =
      sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

    // if the edge is too short
    if (sin_angle < yac_angle_tol) continue;

    double angle = asin(sin_angle);

    double scale = 0.5 * angle / sin_angle * sign;
    for (int k = 0; k < 3; ++k) barycenter[k] += cross[k] * scale;
  }

  // return signed area
  return tri_area(ref, a, b) * sign;
}

/** Compute the signed area of a grid cell and update barycenter
  *
  * This function computes the signed area of a spherical polygon that may
  * contain both great circle edges and latitude circle edges. It simultaneously
  * updates a barycenter accumulator that can be used to compute the cell's
  * center of mass. The polygon is triangulated using a fan triangulation from
  * the first vertex, and the signed areas of the triangles are summed to
  * correctly handle concave cells. For latitude circle edges, correction terms
  * are added to account for the difference between the great circle path and
  * the latitude circle path.
  *
  * \param[in] cell Grid cell structure
  * \param[in,out] barycenter Barycenter accumulator to be updated (3D vector)
  * \param[in] sign Sign multiplier for both area and barycenter contributions
  *
  * \return The signed area of the grid cell on the unit sphere multiplied by
  *         sign parameter
  *
  * \remark For triangles with only great circle edges, uses tri_area_info directly
  * \remark For polygons or cells with latitude edges, uses triangulation
  *         with corrections
  * \remark Returns 0.0 for degenerate cells with fewer than 2 corners
  * \remark The barycenter is accumulated but not normalized by this function
  * \remark The sign parameter allows for proper handling of overlapping regions
  *         in conservative remapping algorithms
  */
double yac_grid_cell_area_info(
  struct yac_grid_cell cell, double * barycenter, double sign) {

  size_t M = cell.num_corners;

  if (M < 2) return 0.0;

  // determine whether the cell contains at least one lat circle edge
  int lat_flag = 0;
  for (size_t i = 0; i < M; i++)
    lat_flag |= cell.edge_type[i] == YAC_LAT_CIRCLE_EDGE;

  // if this is a basic triangle only conisting of great circle edges
  // (e.g. ICON cell)
  if (M == 3 && !lat_flag) {

    // in case of triangle only consisting of great circle edges, return
    // value of the signed area returned by tri_area_info and update
    // barycenter
    return tri_area_info(cell.coordinates_xyz[0],
                         cell.coordinates_xyz[1],
                         cell.coordinates_xyz[2],
                         barycenter, sign);
  }

  // loop over all subtriangles generated by triangulation using
  // a fan approach
  double area = 0.0;
  for (size_t i = 1; i < M - 1; ++i) {

    // signed area summation, this ensures that this routine
    // also works for concave cells
    area +=
      tri_area_info(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[i],
        cell.coordinates_xyz[i+1],
        barycenter, sign);
  }

  // if there is at least one latitude circle edge
  if (lat_flag) {

    // loop over all edges of the polygon
    for (size_t i = 0, j = M - 1; i < M; j = i++) {

      if (cell.edge_type[j] == YAC_LAT_CIRCLE_EDGE) {

        // Add correction for the area difference between the great circle
        // edge and the latitude circle edge from vertex j to vertex i
        // (uses vertex 0 as a reference to determine area sign).
        // Also updates the barycenter accumulator with the contribution
        // from the latitude circle edge correction.
        area += lat_edge_correction_info(
          cell.coordinates_xyz[0],
          cell.coordinates_xyz[j],
          cell.coordinates_xyz[i],
          barycenter, sign);
      }
    }
  }

  return area;
}

/* ----------------------------------- */

static inline double scalar_product(double a[], double b[]) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* ----------------------------------- */

