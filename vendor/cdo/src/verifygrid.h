#ifndef VERIFYGRID_H
#define VERIFYGRID_H

#include "varray.h"

struct Point
{
  double x = 0.0, y = 0.0;
};

struct Point3D
{
  double X = 0.0, Y = 0.0, Z = 0.0;
};

int get_actual_number_of_corners(int ncorner, Varray<Point3D> const &cell_corners_xyz_open_cell);
int get_no_duplicates(int actual_number_of_corners, Varray<Point3D> const &cell_corners_xyz_open_cell,
                      std::vector<bool> &marked_duplicate_indices);
void copy_unique_corners(int actual_number_of_corners, Varray<Point3D> const &cell_corners_xyz_open_cell,
                         std::vector<bool> const &marked_duplicate_indices, Varray<Point3D> &cell_corners_xyz_without_duplicates);

void set_cell_corners_3D(int ncorner, const double *cellCornersLon, const double *cellCornersLat, Varray<Point3D> &cellCorners3D);
Point set_center_point_plane_projection(int coordinateToIgnore, const Point3D &centerPoint3D);
void set_cell_corners_plane_projection(int coordinateToIgnore, int ncorner, Varray<Point3D> const &cellCorners3D,
                                       Varray<Point> &cellCorners2D);
int find_coordinate_to_ignore(Varray<Point3D> const &cell_corners_xyz);
double polygon_area(Varray<Point> const &cellCorners, int numCorners);
bool are_polygon_vertices_arranged_in_clockwise_order(Varray<Point> const &cellCorners, int numCorners);
int winding_numbers_algorithm(Varray<Point> const &cell_corners, int number_corners, const Point &point);

#endif /* VERIFYGRID_H */
