/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Asela Rajapakse
          Uwe Schulzweida

*/

/*
   This module contains the following operators:
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "verifygrid.h"

constexpr double eps = 0.000000001;

static void
print_index_2D(size_t cell_no, size_t nx)
{
  auto iy = cell_no / nx;
  auto ix = cell_no - iy * nx;
  std::fprintf(stdout, " [i=%zu j=%zu]", ix + 1, iy + 1);
}

static void
quick_sort(double *array, size_t array_length)
{
  if (array_length < 2) return;
  auto p = array[array_length / 2];

  size_t i, j;
  for (i = 0, j = array_length - 1;; i++, j--)
  {
    while (array[i] < p) i++;
    while (p < array[j]) j--;
    if (i >= j) break;
    std::swap(array[i], array[j]);
  }
  quick_sort(array, i);
  quick_sort(array + i, array_length - i);
}

/* Quicksort is called with a pointer to the array of center points to be sorted and an integer indicating its length.
 * It sorts the array by its longitude coordinates */
static void
quick_sort_by_lon(double *array, size_t array_length)
{
  if (array_length < 4) return;

  auto p = ((array_length / 2) % 2) ? array[(array_length / 2) + 1] : array[array_length / 2];

  size_t i, j;
  for (i = 0, j = array_length - 2;; i += 2, j -= 2)
  {
    while (array[i] < p) i += 2;

    while (p < array[j]) j -= 2;

    if (i >= j) break;

    std::swap(array[i], array[j]);
    std::swap(array[i + 1], array[j + 1]);
  }

  quick_sort_by_lon(array, i);
  quick_sort_by_lon(array + i, array_length - i);
}

// This uses quicksort to sort the latitude coordinates in a subarray of all coordinates.
static void
quick_sort_of_subarray_by_lat(double *array, size_t subarray_start, size_t subarray_end)
{
  auto subarray_length = (subarray_end - subarray_start) / 2 + 1;
  Varray<double> subarray(subarray_length);
  size_t subarray_index = 0;

  for (size_t index = subarray_start + 1; index <= subarray_end + 1; index += 2)
  {
    subarray[subarray_index] = array[index];
    subarray_index++;
  }

  quick_sort(subarray.data(), subarray_length);

  subarray_index = 0;

  for (size_t index = subarray_start + 1; index <= subarray_end + 1; index += 2)
  {
    array[index] = subarray[subarray_index];
    subarray_index++;
  }
}

static double
determinant(const double (&matrix)[3][3])
{
  // Calculates the determinant for a 3 x 3 matrix.

  return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0]
         + matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0]
         - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1];
}

static Point3D
find_unit_normal(const Point3D &a, const Point3D &b, const Point3D &c)
{
  // Calculates the unit normal for a plane defined on three points a, b, c in Euclidean space.

  const double matrixX[3][3] = { { 1, a.Y, a.Z }, { 1, b.Y, b.Z }, { 1, c.Y, c.Z } };
  const double matrixY[3][3] = { { a.X, 1, a.Z }, { b.X, 1, b.Z }, { c.X, 1, c.Z } };
  const double matrixZ[3][3] = { { a.X, a.Y, 1 }, { b.X, b.Y, 1 }, { c.X, c.Y, 1 } };

  auto x = determinant(matrixX);
  auto y = determinant(matrixY);
  auto z = determinant(matrixZ);

  auto magnitude = std::sqrt(x * x + y * y + z * z);

  Point3D unitNormal;
  unitNormal.X = x / magnitude;
  unitNormal.Y = y / magnitude;
  unitNormal.Z = z / magnitude;

  return unitNormal;
}

int
find_coordinate_to_ignore(Varray<Point3D> const &cellCorners3D)
{
  // Takes the first three corners/vertices of the cell and calculates the unit normal via determinants.

  auto surfaceNormalOfTheCell = find_unit_normal(cellCorners3D[0], cellCorners3D[1], cellCorners3D[2]);

  // The surface normal is used to choose the coordinate to ignore.

  auto absX = std::fabs(surfaceNormalOfTheCell.X);
  auto absY = std::fabs(surfaceNormalOfTheCell.Y);
  auto absZ = std::fabs(surfaceNormalOfTheCell.Z);

  int coordinateToIgnore = 3;

  if (absX > absY)
  {
    if (absX > absZ) coordinateToIgnore = 1;
  }
  else
  {
    if (absY > absZ) coordinateToIgnore = 2;
  }

  return coordinateToIgnore;
}

static inline double
is_point_left_of_edge(const Point &p1, const Point &p2, const Point &p)
{
  /*
     Computes whether a point is left of the line through point1 and point2.
     This is part of the solution to the point in polygon problem.

     Returns: = 0 if the point is on the line through point1 and point2
              > 0 if the point is left of the line
              < 0 if the point is right of the line

     This algorithm is by Dan Sunday (geomalgorithms.com) and is completely free for use and modification.
  */

  return ((p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y));
}

int
winding_numbers_algorithm(Varray<Point> const &cellCorners, int numCorners, const Point &point)
{
  /*
     Computes whether a point is inside the bounds of a cell. This is the solution to the point in polygon problem.
     Returns 0 if the point is outside, returns 1 if the point is inside the cell. Based on an algorithm by Dan Sunday
     (geomalgorithms.com). His algorithm is completely free for use and modification.
  */

  int windingNumber = 0;

  for (int k = 0; k < numCorners - 1; ++k)
  {
    if (cellCorners[k].y <= point.y)
    {
      if (cellCorners[k + 1].y > point.y)
      {
        if (is_point_left_of_edge(cellCorners[k], cellCorners[k + 1], point) > 0.0) windingNumber++;
      }
    }
    else
    {
      if (cellCorners[k + 1].y <= point.y)
      {
        if (is_point_left_of_edge(cellCorners[k], cellCorners[k + 1], point) < 0.0) windingNumber--;
      }
    }
  }

  return windingNumber;
}

static bool
is_center_point_on_corner(Varray<Point> const &cellCorners, int numCorners, const Point &point)
{
  for (int k = 0; k < numCorners; ++k)
  {
    if (is_equal(cellCorners[k].x, point.x) && is_equal(cellCorners[k].y, point.y)) return true;
  }

  return false;
}

static int
sign(double x)
{
  // Is +1 if x is positive, -1 if x is negative and 0 if x is zero.

  return (x > 0.0) - (x < 0.0);
}

static bool
is_simple_polygon_convex(Varray<Point> const &cellCorners, int numCorners)
{
  // Tests in which direction the polygon winds when walking along its edges. Does so for all edges of the polygon.

  double direction = 0.0;

  for (int i = 0; i < numCorners - 2; ++i)
  {
    auto turnsTo = (cellCorners[i].x - cellCorners[i + 1].x) * (cellCorners[i + 1].y - cellCorners[i + 2].y)
                   - (cellCorners[i].y - cellCorners[i + 1].y) * (cellCorners[i + 1].x - cellCorners[i + 2].x);

    // In the first iteration the direction of winding of the entire polygon is set. Better not be 0.

    if (i == 1) direction = turnsTo;

    if (is_not_equal(sign(direction), sign(turnsTo)))
    {
      if (is_not_equal(direction, 0.0)) return false;
    }
    else { direction = turnsTo; }
  }

  return true;
}

double
polygon_area(Varray<Point> const &cellCorners, int numCorners)
{
  /* This algorithm is based on the calculation from Wolfram Mathworld Polygon Area. It results in the area of planar
   * non-self-intersecting polygon. */

  double twiceThePolygonArea = 0.0;

  for (int i = 0; i < numCorners - 1; ++i)
    twiceThePolygonArea += (cellCorners[i].x * cellCorners[i + 1].y) - (cellCorners[i + 1].x * cellCorners[i].y);

  return twiceThePolygonArea * 0.5;
}

bool
are_polygon_vertices_arranged_in_clockwise_order(Varray<Point> const &cellCorners, int numCorners)
{
  auto cellArea = polygon_area(cellCorners, numCorners);

  /* A negative area indicates a clockwise arrangement of vertices, a positive area a counterclockwise arrangement.
   * There should be an area to begin with. */

  return (cellArea < 0.0);
}

static void
create_sorted_point_array(size_t gridsize, Varray<double> const &grid_center_lon, Varray<double> const &grid_center_lat,
                          Varray<double> &center_point_array)
{
  for (size_t cell_no = 0; cell_no < gridsize; cell_no++)
  {
    center_point_array[cell_no * 2 + 0] = grid_center_lon[cell_no];
    center_point_array[cell_no * 2 + 1] = grid_center_lat[cell_no];
  }

  // The cell center points are sorted by their first coordinate (lon) with quicksort.

  quick_sort_by_lon(center_point_array.data(), gridsize * 2);

  // Now the lat coordinates in subarrays that reflect equal lon coordinates are sorted with quicksort.

  int subarray_start = 0;

  for (size_t cell_no = 0; cell_no < gridsize - 1; cell_no++)
  {
    if (cell_no == gridsize - 2)
    {
      int subarray_end = gridsize * 2 - 2;
      quick_sort_of_subarray_by_lat(center_point_array.data(), subarray_start, subarray_end);
    }

    if (std::fabs(center_point_array[cell_no * 2 + 0] - center_point_array[(cell_no + 1) * 2 + 0]) > eps)
    {
      int subarray_end = cell_no * 2;
      if ((subarray_end - subarray_start) > 1)
        quick_sort_of_subarray_by_lat(center_point_array.data(), subarray_start, subarray_end);

      subarray_start = subarray_end + 2;
    }
  }
}

static void
print_header(int gridtype, size_t gridsize, size_t nx, int gridno, int numGrids)
{
  if (numGrids == 1)
  {
    if (nx)
      cdo_print(Blue("Grid consists of %zu (%zux%zu) cells (type: %s), of which"), gridsize, nx, gridsize / nx,
                gridNamePtr(gridtype));
    else
      cdo_print(Blue("Grid consists of %zu cells (type: %s), of which"), gridsize, gridNamePtr(gridtype));
  }
  else
  {
    if (nx)
      cdo_print(Blue("Grid no %u (of %u) consists of %zu (%zux%zu) cells (type: %s), of which"), gridno + 1, numGrids, gridsize, nx,
                gridsize / nx, gridNamePtr(gridtype));
    else
      cdo_print(Blue("Grid no %u (of %u) consists of %zu cells (type: %s), of which"), gridno + 1, numGrids, gridsize,
                gridNamePtr(gridtype));
  }
}

static size_t
get_no_unique_center_points(size_t gridsize, size_t nx, Varray<double> const &grid_center_lon,
                            Varray<double> const &grid_center_lat)
{
  (void) nx;
  // For performing the first test, an array of all center point coordinates is built.
  Varray<double> center_points(gridsize * 2);
  create_sorted_point_array(gridsize, grid_center_lon, grid_center_lat, center_points);

  size_t no_unique_center_points = 1;
  for (size_t cell_no = 0; cell_no < gridsize - 1; cell_no++)
  {
    if (std::fabs(center_points[cell_no * 2 + 0] - center_points[(cell_no + 1) * 2 + 0]) < eps
        && std::fabs(center_points[cell_no * 2 + 1] - center_points[(cell_no + 1) * 2 + 1]) < eps)
    {
      if (Options::cdoVerbose)
      {
        std::fprintf(stdout, "Duplicate point [lon=%.5g lat=%.5g] was found", center_points[cell_no * 2 + 0],
                center_points[cell_no * 2 + 1]);
        /*
        std::fprintf(stdout, "Duplicate point [lon=%.5g lat=%.5g] was found in cell no %zu",
                center_points[cell_no * 2 + 0], center_points[cell_no * 2 + 1], cell_no + 1);
        if (nx) print_index_2D(cell_no, nx);
        */
        std::fprintf(stdout, "\n");
      }
    }
    else { no_unique_center_points++; }
  }

  return no_unique_center_points;
}

void
set_cell_corners_3D(int ncorner, const double *cellCornersLon, const double *cellCornersLat, Varray<Point3D> &cellCorners3D)
{
  double cornerCoordinates[3];
  for (int k = 0; k < ncorner; ++k)
  {
    // Conversion of corner spherical coordinates to Cartesian coordinates.
    gcLLtoXYZ(cellCornersLon[k], cellCornersLat[k], cornerCoordinates);

    // The components of the result vector are appended to the list of cell corner coordinates.
    cellCorners3D[k].X = cornerCoordinates[0];
    cellCorners3D[k].Y = cornerCoordinates[1];
    cellCorners3D[k].Z = cornerCoordinates[2];
  }
}

Point
set_center_point_plane_projection(int coordinateToIgnore, const Point3D &centerPoint3D)
{
  // clang-format off
  switch (coordinateToIgnore)
  {
    case 1: return Point {centerPoint3D.Y, centerPoint3D.Z};
    case 2: return Point {centerPoint3D.Z, centerPoint3D.X};
    case 3: return Point {centerPoint3D.X, centerPoint3D.Y};
  }
  // clang-format on

  return Point{ centerPoint3D.X, centerPoint3D.Y };
}

void
set_cell_corners_plane_projection(int coordinateToIgnore, int ncorner, Varray<Point3D> const &cellCorners3D,
                                  Varray<Point> &cellCorners2D)
{
  switch (coordinateToIgnore)
  {
    case 1:
      for (int k = 0; k <= ncorner; ++k) cellCorners2D[k] = Point{ cellCorners3D[k].Y, cellCorners3D[k].Z };
      break;
    case 2:
      for (int k = 0; k <= ncorner; ++k) cellCorners2D[k] = Point{ cellCorners3D[k].Z, cellCorners3D[k].X };
      break;
    case 3:
      for (int k = 0; k <= ncorner; ++k) cellCorners2D[k] = Point{ cellCorners3D[k].X, cellCorners3D[k].Y };
      break;
  }
}

int
get_actual_number_of_corners(int ncorner, Varray<Point3D> const &cellCorners3D_openCell)
{
  auto actualNumberOfCorners = ncorner;

  for (int k = ncorner - 1; k > 0; k--)
  {
    if (is_equal(cellCorners3D_openCell[k].X, cellCorners3D_openCell[k - 1].X)
        && is_equal(cellCorners3D_openCell[k].Y, cellCorners3D_openCell[k - 1].Y)
        && is_equal(cellCorners3D_openCell[k].Z, cellCorners3D_openCell[k - 1].Z))
      actualNumberOfCorners--;
    else
      break;
  }

  return actualNumberOfCorners;
}

int
get_no_duplicates(int actualNumberOfCorners, Varray<Point3D> const &cellCorners3D_openCell,
                  std::vector<bool> &markedDuplicateIndices)
{
  for (int i = 0; i < actualNumberOfCorners; ++i) markedDuplicateIndices[i] = false;

  int noDuplicates = 0;

  for (int i = 0; i < actualNumberOfCorners; ++i)
    for (int j = i + 1; j < actualNumberOfCorners; ++j)
      if (std::fabs(cellCorners3D_openCell[i].X - cellCorners3D_openCell[j].X) < eps
          && std::fabs(cellCorners3D_openCell[i].Y - cellCorners3D_openCell[j].Y) < eps
          && std::fabs(cellCorners3D_openCell[i].Z - cellCorners3D_openCell[j].Z) < eps)
      {
        noDuplicates++;
        markedDuplicateIndices[j] = true;
      }

  return noDuplicates;
}

void
copy_unique_corners(int actualNumberOfCorners, Varray<Point3D> const &cellCorners3D_openCell,
                    std::vector<bool> const &markedDuplicateIndices, Varray<Point3D> &cellCorners3D_without_duplicates)
{
  int uniqueCornerNumber = 0;

  for (int k = 0; k < actualNumberOfCorners; ++k)
  {
    if (markedDuplicateIndices[k] == false)
    {
      cellCorners3D_without_duplicates[uniqueCornerNumber] = cellCorners3D_openCell[k];
      uniqueCornerNumber++;
    }
  }
}

static inline void
normalise_vector(double v[])
{
  double norm = 1.0 / std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

  v[0] *= norm;
  v[1] *= norm;
  v[2] *= norm;
}

static void
compute_cell_barycenter(size_t numCorners, Varray<Point3D> const &cellCorners3D, double barycenter[3])
{
  barycenter[0] = 0.0;
  barycenter[1] = 0.0;
  barycenter[2] = 0.0;

  for (size_t i = 0; i < numCorners; ++i)
  {
    barycenter[0] += cellCorners3D[i].X;
    barycenter[1] += cellCorners3D[i].Y;
    barycenter[2] += cellCorners3D[i].Z;
  }

  normalise_vector(barycenter);
}

static void
verify_grid(size_t gridsize, size_t nx, int ncorner, Varray<double> const &grid_center_lon, Varray<double> const &grid_center_lat,
            Varray<double> const &grid_corner_lon, Varray<double> const &grid_corner_lat)
{
  /*
     First, this function performs the following test:

     1) it tests whether there are duplicate cells in the given grid by comparing their center point

     Additionally, on each cell of a given grid:

     2) it tests whether all cells are convex and all cell bounds have the same orientation,
        i.e. the corners of the cell are in clockwise or counterclockwise order

     3) it tests whether the center point is within the bounds of the cell

     The results of the tests are printed on stdout.
  */
  Varray<Point3D> cellCorners3D_openCell(ncorner);

  size_t no_of_cells_with_duplicates = 0;
  size_t no_usable_cells = 0;
  size_t no_convex_cells = 0;
  size_t no_clockwise_cells = 0;
  size_t no_counterclockwise_cells = 0;
  size_t no_of_cells_with_center_points_out_of_bounds = 0;
  size_t no_of_cells_with_center_points_on_corner = 0;

  std::vector<int> no_cells_with_a_specific_no_of_corners(ncorner, 0);

  // Checking for the number of unique center point coordinates.
  auto no_unique_center_points = get_no_unique_center_points(gridsize, nx, grid_center_lon, grid_center_lat);

  // used only actualNumberOfCorners
  std::vector<bool> markedDuplicateIndices(ncorner);
  Varray<Point3D> cellCorners3D(ncorner + 1);
  Varray<Point> cellCornersPlaneProjection(ncorner + 1);
  Varray<size_t> cells_with_duplicates;
  cells_with_duplicates.reserve(gridsize);

  /*
     Latitude and longitude are spherical coordinates on a unit circle. Each such coordinate tuple is transformed into a
     triple of Cartesian coordinates in Euclidean space. This is first done for the presumed center point of the cell
     and then for all the corners of the cell.
  */

  for (size_t cell_no = 0; cell_no < gridsize; cell_no++)
  {
    // Conversion of center point spherical coordinates to Cartesian coordinates.
    double centerCoordinates[3];
    gcLLtoXYZ(deg_to_rad(grid_center_lon[cell_no]), deg_to_rad(grid_center_lat[cell_no]), centerCoordinates);
    Point3D centerPoint3D;
    centerPoint3D.X = centerCoordinates[0];
    centerPoint3D.Y = centerCoordinates[1];
    centerPoint3D.Z = centerCoordinates[2];

    double cornerCoordinates[3];
    for (int k = 0; k < ncorner; ++k)
    {
      // Conversion of corner spherical coordinates to Cartesian coordinates.
      auto index = cell_no * ncorner + k;
      gcLLtoXYZ(deg_to_rad(grid_corner_lon[index]), deg_to_rad(grid_corner_lat[index]), cornerCoordinates);

      // The components of the result vector are appended to the list of cell corner coordinates.
      cellCorners3D_openCell[k].X = cornerCoordinates[0];
      cellCorners3D_openCell[k].Y = cornerCoordinates[1];
      cellCorners3D_openCell[k].Z = cornerCoordinates[2];
    }

    /*
       Not all cells have the same number of corners. The array, however, has ncorner * 3  values for each cell, where
       ncorner is the maximum number of corners. Unused values have been filled with the values of the final cell. The
       following identifies the surplus corners and gives the correct length of the cell.
    */

    auto actualNumberOfCorners = get_actual_number_of_corners(ncorner, cellCorners3D_openCell);

    no_cells_with_a_specific_no_of_corners[actualNumberOfCorners - 1]++;

    // If there are less than three corners in the cell, it is unusable and considered degenerate. No area can be computed.

    if (actualNumberOfCorners < 3)
    {
      if (Options::cdoVerbose)
      {
        std::fprintf(stdout, "Less than three vertices found in cell no %zu", cell_no + 1);
        if (nx) print_index_2D(cell_no, nx);
        std::fprintf(stdout, ", omitted!");
        std::fprintf(stdout, "\n");
      }
      continue;
    }

    no_usable_cells++;

    // Checks if there are any duplicate vertices in the list of corners. Note that the last (additional) corner has not been set
    // yet.

    auto noDuplicates = get_no_duplicates(actualNumberOfCorners, cellCorners3D_openCell, markedDuplicateIndices);

    if (Options::cdoVerbose)
    {
      for (int i = 0; i < actualNumberOfCorners; ++i)
      {
        if (markedDuplicateIndices[i])
        {
          std::fprintf(stdout, "Duplicate vertex [lon=%.5g lat=%.5g] was found in cell no %zu", grid_corner_lon[cell_no * ncorner + i],
                  grid_corner_lat[cell_no * ncorner + i], cell_no + 1);
          if (nx) print_index_2D(cell_no, nx);
          std::fprintf(stdout, "\n");
        }
      }
    }

    // Writes the unique corner vertices in a new array.

    copy_unique_corners(actualNumberOfCorners, cellCorners3D_openCell, markedDuplicateIndices, cellCorners3D);

    actualNumberOfCorners -= noDuplicates;

    // We are creating a closed polygon/cell by setting the additional last corner to be the same as the first one.

    cellCorners3D[actualNumberOfCorners] = cellCorners3D[0];

    if (noDuplicates != 0)
    {
      no_of_cells_with_duplicates++;
      cells_with_duplicates.push_back(cell_no);
    }

    /* If there are less than three corners in the cell left after removing duplicates, it is unusable and considered
     * degenerate. No area can be computed. */

    if (actualNumberOfCorners < 3)
    {
      if (Options::cdoVerbose)
        std::fprintf(stdout,
                "Less than three vertices found in cell no %zu. This cell is considered degenerate and "
                "will be omitted from further computation!\n",
                cell_no + 1);

      continue;
    }

    double barycenter[3];
    compute_cell_barycenter(actualNumberOfCorners, cellCorners3D, barycenter);

    auto coordinateToIgnore = find_coordinate_to_ignore(cellCorners3D);

    /* The remaining two-dimensional coordinates are extracted into one array for all the cell's corners and into one
     * array for the center point. */

    /* The following projection on the plane that two coordinate axes lie on changes the arrangement of the polygon
       vertices if the coordinate to be ignored along the third axis is smaller than 0. In this case, the result of
       the computation of the orientation of vertices needs to be inverted. Clockwise becomes counterclockwise and
       vice versa. */

    auto cval
        = (coordinateToIgnore == 1) ? cellCorners3D[0].X : ((coordinateToIgnore == 2) ? cellCorners3D[0].Y : cellCorners3D[0].Z);
    auto invertResult = (cval < 0.0);

    auto centerPoint2D = set_center_point_plane_projection(coordinateToIgnore, centerPoint3D);
    set_cell_corners_plane_projection(coordinateToIgnore, actualNumberOfCorners, cellCorners3D, cellCornersPlaneProjection);

    // Checking for convexity of the cell.

    if (is_simple_polygon_convex(cellCornersPlaneProjection, actualNumberOfCorners + 1)) { no_convex_cells++; }
    else if (Options::cdoVerbose)
    {
      std::fprintf(stdout, "Vertices are not convex in cell no %zu", cell_no + 1);
      if (nx) print_index_2D(cell_no, nx);
      std::fprintf(stdout, "\n");
    }

    // Checking the arrangement or direction of cell vertices.

    auto isClockwise = are_polygon_vertices_arranged_in_clockwise_order(cellCornersPlaneProjection, actualNumberOfCorners + 1);

    /* If the direction of the vertices was flipped during the projection onto the two-dimensional plane, the previous
     * result needs to be inverted now. */

    if (invertResult) isClockwise = !isClockwise;

    // The overall counter of (counter)clockwise cells is increased by one.
    if (isClockwise)
    {
      no_clockwise_cells++;
      if (Options::cdoVerbose)
      {
        std::fprintf(stdout, "Vertices arranged in a clockwise order in cell no %zu", cell_no + 1);
        if (nx) print_index_2D(cell_no, nx);
        std::fprintf(stdout, "\n");
      }
    }
    else { no_counterclockwise_cells++; }

    // The winding numbers algorithm is used to test whether the presumed center point is within the bounds of the cell.
    auto windingNumber = winding_numbers_algorithm(cellCornersPlaneProjection, actualNumberOfCorners + 1, centerPoint2D);

    if (windingNumber == 0)
    {
      if (is_center_point_on_corner(cellCornersPlaneProjection, actualNumberOfCorners, centerPoint2D))
      {
        no_of_cells_with_center_points_on_corner++;
      }
      else
      {
        no_of_cells_with_center_points_out_of_bounds++;

        if (Options::cdoVerbose)
        {
          std::fprintf(stdout, "Center point [lon=%.5g lat=%.5g] outside bounds [", grid_center_lon[cell_no], grid_center_lat[cell_no]);
          for (int k = 0; k < actualNumberOfCorners; ++k)
            std::fprintf(stdout, " %.5g/%.5g", grid_corner_lon[cell_no * ncorner + k], grid_corner_lat[cell_no * ncorner + k]);
          std::fprintf(stdout, "] in cell no %zu", cell_no + 1);
          if (nx) print_index_2D(cell_no, nx);
          std::fprintf(stdout, "\n");
        }
      }
    }
  }

  auto no_nonunique_cells = gridsize - no_unique_center_points;
  auto no_nonconvex_cells = gridsize - no_convex_cells;
  auto no_nonusable_cells = gridsize - no_usable_cells;

  for (int i = 2; i < ncorner; ++i)
    if (no_cells_with_a_specific_no_of_corners[i])
      cdo_print(Blue("%9d cells have %d vertices"), no_cells_with_a_specific_no_of_corners[i], i + 1);

  if (no_of_cells_with_duplicates) cdo_print(Blue("%9zu cells have duplicate vertices"), no_of_cells_with_duplicates);

  if (no_nonusable_cells) cdo_print(Red("%9zu cells have unusable vertices"), no_nonusable_cells);
  if (no_nonunique_cells) cdo_print(Red("%9zu cells are not unique"), no_nonunique_cells);
  if (no_nonconvex_cells) cdo_print(Red("%9zu cells are non-convex"), no_nonconvex_cells);
  if (no_clockwise_cells) cdo_print(Red("%9zu cells have their vertices arranged in a clockwise order"), no_clockwise_cells);

  if (no_of_cells_with_center_points_out_of_bounds)
    cdo_print(Red("%9zu cells have their center point located outside their boundaries"),
              no_of_cells_with_center_points_out_of_bounds);

  if (no_of_cells_with_center_points_on_corner)
    cdo_print(Red("%9zu cells have their center point located on their corners"), no_of_cells_with_center_points_on_corner);
}

static void
print_lonlat(int dig, int gridID, size_t gridsize, Varray<double> const &lon, Varray<double> const &lat)
{
  size_t num_lons_nan = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (std::isnan(lon[i])) num_lons_nan++;
  if (num_lons_nan) cdo_print(Red("%9zu longitudes are NaN"), num_lons_nan);

  size_t num_lats_nan = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (std::isnan(lat[i])) num_lats_nan++;
  if (num_lats_nan) cdo_print(Red("%9zu latitudes are NaN"), num_lats_nan);

  auto xname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
  auto xmm = varray_min_max(gridsize, lon);
  cdo_print("%10s : %.*g to %.*g degrees", xname, dig, xmm.min, dig, xmm.max);
  if (xmm.min < -720.0 || xmm.max > 720.0) cdo_warning("Grid cell center longitudes out of range!");

  auto yname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
  auto ymm = varray_min_max(gridsize, lat);
  cdo_print("%10s : %.*g to %.*g degrees", yname, dig, ymm.min, dig, ymm.max);
  if (ymm.min < -90.00001 || ymm.max > 90.00001) cdo_warning("Grid cell center latitudes out of range!");
}

class Verifygrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Verifygrid",
    .operators = { { "verifygrid", VerifygridHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Verifygrid> registration = RegisterEntry<Verifygrid>();

  CdoStreamID streamID{};
  int vlistID{};

  int numGrids{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID = cdo_open_read(0);
    vlistID = cdo_stream_inq_vlist(streamID);

    numGrids = vlistNumGrids(vlistID);
  }

  void
  run() override
  {
    for (int gridno = 0; gridno < numGrids; ++gridno)
    {
      auto lgeo = true;

      auto gridID = vlistGrid(vlistID, gridno);
      auto gridtype = gridInqType(gridID);
      auto isHealpixGrid = is_healpix_grid(gridID);

      if (gridtype == GRID_GENERIC || gridtype == GRID_SPECTRAL) lgeo = false;
      if (isHealpixGrid) lgeo = false;

      if (gridtype == GRID_GME) gridID = gridToUnstructured(gridID, NeedCorners::Yes);

      if (lgeo && gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
        gridID = gridToCurvilinear(gridID, NeedCorners::Yes);

      if (gridtype == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID))
      {
        auto reference = dereferenceGrid(gridID);
        if (reference.isValid) gridID = reference.gridID;
        if (reference.notFound)
        {
          cdo_print("Reference to source grid not found!");
          lgeo = false;
        }
      }

      auto gridsize = gridInqSize(gridID);
      auto nx = gridInqXsize(gridID);
      auto ny = gridInqYsize(gridID);
      nx = (nx * ny == gridsize) ? nx : 0;

      if (lgeo)
      {
        auto ncorner = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

        if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

        Varray<double> grid_center_lat(gridsize), grid_center_lon(gridsize);

        gridInqYvals(gridID, grid_center_lat.data());
        gridInqXvals(gridID, grid_center_lon.data());

        // Convert lat/lon units if required
        cdo_grid_to_degree(gridID, CDI_XAXIS, grid_center_lon, "grid center lon");
        cdo_grid_to_degree(gridID, CDI_YAXIS, grid_center_lat, "grid center lat");

        print_header(gridtype, gridsize, nx, gridno, numGrids);

        if (gridHasBounds(gridID))
        {
          if (ncorner == 0) cdo_abort("Number of cell corners undefined!");
          auto nalloc = ncorner * gridsize;
          Varray<double> grid_corner_lat(nalloc), grid_corner_lon(nalloc);

          gridInqYbounds(gridID, grid_corner_lat.data());
          gridInqXbounds(gridID, grid_corner_lon.data());

          cdo_grid_to_degree(gridID, CDI_XAXIS, grid_corner_lon, "grid corner lon");
          cdo_grid_to_degree(gridID, CDI_YAXIS, grid_corner_lat, "grid corner lat");

          verify_grid(gridsize, nx, ncorner, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);
        }
        else { cdo_warning("Grid cell corner coordinates missing!"); }

        int dig = Options::CDO_flt_digits;
        print_lonlat(dig, gridID, gridsize, grid_center_lon, grid_center_lat);
      }
      else
      {
        if (numGrids == 1)
          cdo_print(Blue("Grid consists of %zu points (type: %s)"), gridsize, isHealpixGrid ? "healpix" : gridNamePtr(gridtype));
        else
          cdo_print(Blue("Grid no %u (of %u) consists of %zu points (type: %s)"), gridno + 1, numGrids, gridsize,
                    gridNamePtr(gridtype));
        // cdo_print("");
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
