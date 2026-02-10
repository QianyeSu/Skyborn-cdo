/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>
#include <vector>

#include <cdi.h>
#include "cdo_cdi_wrapper.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "cdo_output.h"
#include "compare.h"
#include "cdo_omp.h"
#include "progress.h"
#include "cdo_omp.h"

static void
lonlat_to_xyz(double lon, double lat, double *xyz)
{
  auto coslat = std::cos(lat);
  xyz[0] = coslat * std::cos(lon);
  xyz[1] = coslat * std::sin(lon);
  xyz[2] = std::sin(lat);
}

static void
cross_product(const double *a, const double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/*
static
double scalar_product(const double *a, const double *b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}
*/

static double
norm(const double *a)
{
  return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
/*
static
double mod_cell_area(int num_corners, double *cell_corner_lon, double
*cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  // generalised version based on the ICON code, mo_base_geometry.f90
  // provided by Luis Kornblueh, MPI-M.

  int M = num_corners; // number of vertices
  int m; // loop index over number of vertices M
  int i; // loop index over the three dimensions

  double area = 0.0;
  double s[M];
  double ca[M];
  double a[M];

  double p[M][3];
  double u[M][3];

  // Convert into cartesian coordinates
  for ( m = 0; m < M; m++ )
    lonlat_to_xyz(cell_corner_lon[m], cell_corner_lat[m], p[m]);

  // First, compute cross products Uij = Vi x Vj.
  for ( m = 0; m < M; m++ )
    cross_product(p[m], p[(m+1)%M], u[m]);

  //  Normalize Uij to unit vectors.
  for ( m = 0; m < M; m++ )
    {
      s[m] = norm(u[m]);
      area += s[m];
    }

  // Test for a degenerated cells associated with collinear vertices.

  if ( std::fabs(area) > 0.0 )
    {
      for ( m = 0; m < M; m++ )
        s[m] = std::sqrt(s[m]);

      for ( m = 0; m < M; m++ )
        for ( i = 0; i < 3; i++ )
          u[m][i] = u[m][i]/s[m];

      //  Compute interior angles Ai as the dihedral angles between planes
      //  by using the definition of the scalar product
      //
      //	    ab = |a| |b| cos (phi)
      //
      //  As a and b are already normalised this reduces to
      //
      //            ab = cos (phi)

      //  There is no explanation so far for the - in the loop below.
      //  But otherwise we don't get the correct results for triangles
      //  and cells. Must have something to do with the theorem.

      for ( m = 0; m < M; m++ )
        {
          ca[m] = - scalar_product(u[m], u[(m+1)%M]);
          if ( ca[m] < -1.0 ) ca[m] = -1.0;
          if ( ca[m] >  1.0 ) ca[m] =  1.0;
          a[m] = std::acos(ca[m]);
        }

      //  Compute areas = a1 + a2 + a3 - (M-2) * pi. here for a unit sphere:
      area = - (double) (M-2) * M_PI;

      for ( m = 0; m < M; m++ )
        area += a[m];

      // area *= EarthRadius * EarthRadius;
      if ( area < 0.0 ) area = 0.0;
    }

  return area;
}
*/

/** area of a spherical triangle based on L'Huilier's Theorem
 *
 * source code is taken from code by Robert Oehmke of Earth System Modeling
 * Framework (www.earthsystemmodeling.org)
 *
 * the license statement for this routine is as follows:
 * Earth System Modeling Framework
 * Copyright 2002-2013, University Corporation for Atmospheric Research,
 * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
 * Laboratory, University of Michigan, National Centers for Environmental
 * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
 * NASA Goddard Space Flight Center.
 * Licensed under the University of Illinois-NCSA License.
 */
static double
mod_tri_area(const double *u, const double *v, const double *w)
{
  double tmp_vec[3];

  cross_product(u, v, tmp_vec);
  double sina = std::sqrt(norm(tmp_vec));
  double a = std::asin(sina);

  cross_product(u, w, tmp_vec);
  double sinb = std::sqrt(norm(tmp_vec));
  double b = std::asin(sinb);

  cross_product(w, v, tmp_vec);
  double sinc = std::sqrt(norm(tmp_vec));
  double c = std::asin(sinc);

  double s = 0.5 * (a + b + c);

  double t = std::tan(s * 0.5) * std::tan((s - a) * 0.5) * std::tan((s - b) * 0.5) * std::tan((s - c) * 0.5);

  double area = std::fabs(4.0 * std::atan(std::sqrt(std::fabs(t))));

  return area;
}

/*
 * source code is taken from code by Robert Oehmke of Earth System Modeling
 * Framework (www.earthsystemmodeling.org) and adjusted to CDO data structures
 *
 * the license statement for this routine is as follows:
 * Earth System Modeling Framework
 * Copyright 2002-2013, University Corporation for Atmospheric Research,
 * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
 * Laboratory, University of Michigan, National Centers for Environmental
 * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
 * NASA Goddard Space Flight Center.
 * Licensed under the University of Illinois-NCSA License.
 */
static double
mod_huiliers_area(int num_corners, double *cell_corner_lon, double *cell_corner_lat)
{
  if (num_corners < 3) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt1);
  lonlat_to_xyz(cell_corner_lon[1], cell_corner_lat[1], pnt2);

  for (int i = 2; i < num_corners; ++i)
  {
    // points that make up a side of cell
    lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);

    // compute angle for pnt2
    sum += mod_tri_area(pnt1, pnt2, pnt3);

    if (i < (num_corners - 1))
    {
      pnt2[0] = pnt3[0];
      pnt2[1] = pnt3[1];
      pnt2[2] = pnt3[2];
    }
  }

  return sum;
}

static double
mod_huiliers_area2(int num_corners, double *cell_corner_lon, double *cell_corner_lat, double cell_center_lon,
                   double cell_center_lat)
{
  if (num_corners < 3) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_center_lon, cell_center_lat, pnt1);
  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt2);

  for (int i = 1; i < num_corners; ++i)
  {
    if (is_equal(cell_corner_lon[i], cell_corner_lon[i - 1]) && is_equal(cell_corner_lat[i], cell_corner_lat[i - 1])) continue;

    // points that make up a side of cell
    lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);

    // compute angle for pnt2
    sum += mod_tri_area(pnt1, pnt2, pnt3);

    pnt2[0] = pnt3[0];
    pnt2[1] = pnt3[1];
    pnt2[2] = pnt3[2];
  }

  if (!(is_equal(cell_corner_lon[0], cell_corner_lon[num_corners - 1])
        && is_equal(cell_corner_lat[0], cell_corner_lat[num_corners - 1])))
  {
    lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt3);
    sum += mod_tri_area(pnt1, pnt2, pnt3);
  }

  return sum;
}

static void
getLonLatCorner(size_t nx, size_t idx, Varray<double> const &grid_corner_lon, Varray<double> const &grid_corner_lat,
                double (&lons)[4], double (&lats)[4])
{
  auto j = idx / nx;
  auto i = idx - j * nx;

  lons[0] = grid_corner_lon[2 * i];
  lons[1] = grid_corner_lon[2 * i + 1];
  lons[2] = grid_corner_lon[2 * i + 1];
  lons[3] = grid_corner_lon[2 * i];

  if (grid_corner_lat[2 * j + 1] > grid_corner_lat[2 * j])
  {
    lats[0] = grid_corner_lat[2 * j];
    lats[1] = grid_corner_lat[2 * j];
    lats[2] = grid_corner_lat[2 * j + 1];
    lats[3] = grid_corner_lat[2 * j + 1];
  }
  else
  {
    lats[0] = grid_corner_lat[2 * j + 1];
    lats[1] = grid_corner_lat[2 * j + 1];
    lats[2] = grid_corner_lat[2 * j];
    lats[3] = grid_corner_lat[2 * j];
  }
}

static int
gen_gridcellarea_reg2d(int gridID, Varray<double> &area, bool lweights)
{
  // lweights = true : area used only to calculate weights

  int status = 0;
  std::string unitstr;
  auto gridsize = gridInqSize(gridID);
  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);
  if (!nlon) nlon = 1;
  if (!nlat) nlat = 1;

  if (!(gridHasCoordinates(gridID) || gridHasBounds(gridID) || (nlon == 1 && gridInqYvals(gridID, nullptr))))
  {
    cdo_warning("Computation of grid cell area weights failed, grid cell center and bounds coordinates missing!");
    return 1;
  }

  Varray<double> grid_corner_lon(nlon * 2), grid_corner_lat(nlat * 2);

  if (gridInqXbounds(gridID, nullptr))
  {
    gridInqXbounds(gridID, grid_corner_lon.data());
    unitstr = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
  }
  else if (nlon == 1)
  {
    if (lweights)
    {
      grid_corner_lon[0] = 0;
      grid_corner_lon[1] = .01;
      unitstr = "radian";
    }
    else
      return 1;
  }
  else
  {
    Varray<double> grid_center_lon(nlon);
    gridInqXvals(gridID, grid_center_lon.data());
    grid_gen_bounds(nlon, grid_center_lon, grid_corner_lon);
    unitstr = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
  }

  if (gridInqYbounds(gridID, nullptr))
  {
    gridInqYbounds(gridID, grid_corner_lat.data());
    unitstr = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
  }
  else if (nlat == 1)
  {
    if (lweights)
    {
      grid_corner_lat[0] = 0;
      grid_corner_lat[1] = .01;
      unitstr = "radian";
    }
    else
      return 1;
  }
  else
  {
    Varray<double> grid_center_lat(nlat);
    gridInqYvals(gridID, grid_center_lat.data());
    grid_gen_bounds(nlat, grid_center_lat, grid_corner_lat);
    grid_check_lat_borders(nlat * 2, grid_corner_lat);
    unitstr = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
  }

  if (string_to_LonLatUnits(unitstr, "grid cell corner longitudes") == LonLatUnits::Deg) grid_to_radian(grid_corner_lon);
  if (string_to_LonLatUnits(unitstr, "grid cell corner latitudes") == LonLatUnits::Deg) grid_to_radian(grid_corner_lat);

  constexpr int nv = 4;
  std::atomic<size_t> atomicCount{ 0 };

  cdo::Progress progress;

#ifdef _OPENMP
#pragma omp parallel for if (gridsize > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && gridsize > progressMinSize) progress.update((double) atomicCount / gridsize);

    double lons[4], lats[4];
    getLonLatCorner(nlon, i, grid_corner_lon, grid_corner_lat, lons, lats);
    area[i] = mod_huiliers_area(nv, lons, lats);
  }

  return status;
}

int
gridGenAreaReg2Dweights(int gridID, Varray<double> &area)
{
  bool lweights = true;
  return gen_gridcellarea_reg2d(gridID, area, lweights);
}

static int
gen_gridcellarea_unstruct(int gridID, Varray<double> &area)
{
  int status = 0;
  auto lgriddestroy = false;

  auto gridsize = gridInqSize(gridID);
  auto gridtype = gridInqType(gridID);

  if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
  {
    if (gridtype == GRID_GME || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      lgriddestroy = true;
      gridID = gridToUnstructured(gridID, NeedCorners::Yes);
    }
    else
    {
      lgriddestroy = true;
      gridID = gridToCurvilinear(gridID, NeedCorners::Yes);
    }
  }

  if (gridtype == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID))
  {
    auto reference = dereferenceGrid(gridID);
    if (reference.exists)
    {
      if (!reference.isValid)
      {
        cdo_warning("Reference to source grid not found!");
        return 1;
      }
      gridID = reference.gridID;
      lgriddestroy = true;
    }
  }

  gridtype = gridInqType(gridID);

  size_t nv = (gridtype == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  if (!gridHasCoordinates(gridID))
  {
    cdo_warning("Computation of grid cell area weights failed, grid cell center coordinates missing!");
    return 1;
  }

  if (nv == 0)
  {
    cdo_warning("Computation of grid cell area weights failed, grid cell corner coordinates missing!");
    return 1;
  }

  Varray<double> grid_center_lon(gridsize), grid_center_lat(gridsize);
  gridInqXvals(gridID, grid_center_lon.data());
  gridInqYvals(gridID, grid_center_lat.data());

  Varray<double> grid_corner_lon(nv * gridsize), grid_corner_lat(nv * gridsize);

  if (!gridHasBounds(gridID)) return 1;

  gridInqXbounds(gridID, grid_corner_lon.data());
  gridInqYbounds(gridID, grid_corner_lat.data());

  cdo_grid_to_radian(gridID, CDI_XAXIS, grid_center_lon, "grid1 center longitudes");
  cdo_grid_to_radian(gridID, CDI_XAXIS, grid_corner_lon, "grid1 corner longitudes");

  cdo_grid_to_radian(gridID, CDI_YAXIS, grid_center_lat, "grid1 center latitudes");
  cdo_grid_to_radian(gridID, CDI_YAXIS, grid_corner_lat, "grid1 corner latitudes");

  if (lgriddestroy) gridDestroy(gridID);

  std::atomic<size_t> atomicCount{ 0 };

  cdo::Progress progress;

#ifdef _OPENMP
#pragma omp parallel for if (gridsize > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && gridsize > progressMinSize) progress.update((double) atomicCount / gridsize);

    if (nv <= 4)
      area[i] = mod_huiliers_area(nv, &grid_corner_lon[i * nv], &grid_corner_lat[i * nv]);
    else
      area[i] = mod_huiliers_area2(nv, &grid_corner_lon[i * nv], &grid_corner_lat[i * nv], grid_center_lon[i], grid_center_lat[i]);
  }

  return status;
}

int
gen_gridcellarea_healpix(int gridID, Varray<double> &area)
{
  auto gridsize = gridInqSize(gridID);
  auto cellarea = 4.0 * M_PI / gridsize;
  for (size_t i = 0; i < gridsize; ++i) { area[i] = cellarea; }
  return 0;
}

int
gridGenArea(int gridID, Varray<double> &area)
{
  int status = 0;

  auto gridsize = gridInqSize(gridID);
  auto gridtype = gridInqType(gridID);

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)
  {
    bool lweights = false;
    status = gen_gridcellarea_reg2d(gridID, area, lweights);
  }
  else if (is_healpix_grid(gridID)) { status = gen_gridcellarea_healpix(gridID, area); }
  else if (gridProjIsSupported(gridID) || gridtype == GRID_GME || gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED
           || gridtype == GRID_GAUSSIAN_REDUCED)
  {
    status = gen_gridcellarea_unstruct(gridID, area);
  }
  else { cdo_abort("Internal error! Unsupported gridtype: %s", gridNamePtr(gridtype)); }

  if (gridVerbose) cdo_print("Total area = %g steradians", varray_sum(gridsize, area));

  if (gridsize < 20 && is_equal(varray_sum(gridsize, area), 0.0)) status = 2;

  return status;
}
