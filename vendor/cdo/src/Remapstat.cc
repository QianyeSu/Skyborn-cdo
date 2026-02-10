/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "cdo_math.h"
#include "cdo_timer.h"
#include "process_int.h"
#include "griddes.h"
#include "grid_pointsearch.h"
#include <mpim_grid.h>
#include "cdo_omp.h"
#include "verifygrid.h"
#include "field_functions.h"

#define TESTIMPL 1
// #define USE_YAC 1
#ifdef TESTIMPL
#ifdef USE_YAC
extern "C"
{
#include "lib/yac/src/geometry.h"
}
#endif
#endif

constexpr double PI2 = M_PI * 2.0;

double radiusDegToKm(double radiusInDeg);

namespace
{
struct StatInfo
{
  size_t n = 0;
  size_t min_svals = SIZE_MAX, max_svals = 0, sum_svals = 0;
  size_t min_nvals = SIZE_MAX, max_nvals = 0, sum_nvals = 0;
  double min_radius = 1.e33, max_radius = 0.0, sum_radius = 0.0;

  void
  add(size_t svals, size_t nvals, double radius)
  {
    radius *= RAD2DEG;
    n++;
    sum_svals += svals;
    sum_nvals += nvals;
    sum_radius += radius;
    min_svals = std::min(min_svals, svals);
    min_nvals = std::min(min_nvals, nvals);
    min_radius = std::min(min_radius, radius);
    max_svals = std::max(max_svals, svals);
    max_nvals = std::max(max_nvals, nvals);
    max_radius = std::max(max_radius, radius);
  };

  void
  print()
  {
    cdo_print("N=%zu", n);
    cdo_print("Min:   svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", min_svals, min_nvals, min_radius, radiusDegToKm(min_radius));
    cdo_print("Mean:  svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", sum_svals / n, sum_nvals / n, sum_radius / n,
              radiusDegToKm(sum_radius / n));
    cdo_print("Max:   svals=%3zu  nvals=%2zu  radius=%.3gdeg(%.3gkm)", max_svals, max_nvals, max_radius, radiusDegToKm(max_radius));
  };
};
}  // namespace

static size_t
read_target_cell_bounds(int gridID, Varray<double> &xbounds, Varray<double> &ybounds)
{
  if (!gridHasBounds(gridID)) cdo_abort("Target cell corner coordinates missing!");

  auto nv = gridInqNvertex(gridID);
  auto gridsize = gridInqSize(gridID);
  xbounds.resize(nv * gridsize);
  ybounds.resize(nv * gridsize);

  gridInqXbounds(gridID, xbounds.data());
  gridInqYbounds(gridID, ybounds.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xbounds, "grid corner lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, ybounds, "grid corner lat");

  return nv;
}

static void
read_coordinates(int gridID, Varray<double> &xvals, Varray<double> &yvals)
{
  auto gridID0 = gridID;

  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  auto gridsize = gridInqSize(gridID);
  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yvals, "grid center lat");

  for (size_t i = 0; i < gridsize; ++i)
  {
    if (xvals[i] > PI2) xvals[i] -= PI2;
    if (xvals[i] < 0.0) xvals[i] += PI2;
  }

  if (gridID0 != gridID) gridDestroy(gridID);
}

static void
read_xbounds_reg2d(int gridID, Varray<double> &xbounds2d)
{
  auto nlon = gridInqXsize(gridID);
  xbounds2d.resize(nlon * 2);

  if (gridInqXbounds(gridID, nullptr)) { gridInqXbounds(gridID, xbounds2d.data()); }
  else
  {
    Varray<double> xvals(nlon);
    gridInqXvals(gridID, xvals.data());
    grid_gen_bounds(nlon, xvals, xbounds2d);
  }

  for (size_t i = 0; i < 2 * nlon; ++i) xbounds2d[i] *= DEG2RAD;

  for (size_t i = 0; i < nlon; ++i)
  {
    if (xbounds2d[2 * i + 1] > PI2)
    {
      xbounds2d[2 * i] -= PI2;
      xbounds2d[2 * i + 1] -= PI2;
    }
    if (xbounds2d[2 * i + 1] < 0.0)
    {
      xbounds2d[2 * i] += PI2;
      xbounds2d[2 * i + 1] += PI2;
    }
  }
}

static void
read_ybounds_reg2d(int gridID, Varray<double> &ybounds2d)
{
  auto nlat = gridInqYsize(gridID);
  ybounds2d.resize(nlat * 2);

  if (gridInqYbounds(gridID, nullptr)) { gridInqYbounds(gridID, ybounds2d.data()); }
  else
  {
    Varray<double> yvals(nlat);
    gridInqYvals(gridID, yvals.data());
    grid_gen_bounds(nlat, yvals, ybounds2d);
    grid_check_lat_borders(2 * nlat, ybounds2d);
  }

  for (size_t i = 0; i < 2 * nlat; ++i) ybounds2d[i] *= DEG2RAD;
}

static bool
is_zonal_grid(int gridID, Varray<double> const &xbounds2d, Varray<double> const &ybounds2d)
{
  bool isZonalGrid{ false };

  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);

  if (nlon == 1 && xbounds2d[1] > xbounds2d[0])
  {
    isZonalGrid = true;
    for (size_t i = 1; i < nlat; ++i)
    {
      if (is_not_equal(ybounds2d[2 * i - 1], ybounds2d[2 * i]))
      {
        isZonalGrid = false;
        break;
      }
    }
  }

  return isZonalGrid;
}

static double
calc_maxdist(size_t i, size_t nv, PointLonLat const &pointLL, Varray<double> const &xbounds, Varray<double> const &ybounds)
{
  double p1[3], p2[3];
  gcLLtoXYZ(pointLL.lon(), pointLL.lat(), p1);

  double maxdist = 0.0;
  for (size_t k = 0; k < nv; ++k)
  {
    auto lons = &xbounds[nv * i];
    auto lats = &ybounds[nv * i];
    gcLLtoXYZ(lons[k], lats[k], p2);
    auto sdist = cdo::sqr_distance(p1, p2);
    maxdist = std::max(maxdist, sdist);
  }

  return maxdist;
}

static double
calc_maxdist_rec2d(size_t i, size_t nlon, PointLonLat const &pointLL, Varray<double> const &xbounds, Varray<double> const &ybounds)
{
  double p1[3], p2[3];
  gcLLtoXYZ(pointLL.lon(), pointLL.lat(), p1);

  constexpr int nv = 4;
  double lons[nv], lats[nv];
  auto iy = i / nlon;
  auto ix = i - iy * nlon;
  lons[0] = xbounds[2 * ix];
  lons[1] = xbounds[2 * ix];
  lons[2] = xbounds[2 * ix + 1];
  lons[3] = xbounds[2 * ix + 1];
  lats[0] = ybounds[2 * iy + 1];
  lats[1] = ybounds[2 * iy];
  lats[2] = ybounds[2 * iy];
  lats[3] = ybounds[2 * iy + 1];

  double maxdist = 0.0;
  for (size_t k = 0; k < nv; ++k)
  {
    gcLLtoXYZ(lons[k], lats[k], p2);
    auto sdist = cdo::sqr_distance(p1, p2);
    maxdist = std::max(maxdist, sdist);
  }

  return maxdist;
}

#ifdef USE_YAC
static size_t
find_points_yac(Vmask &vmask, size_t cell_no, size_t ncorner, size_t numIndices, Varray<size_t> &indices,
                Varray<double> const &xvals1, Varray<double> const &yvals1, Varray<double> const &xbounds,
                Varray<double> const &ybounds)
{
#ifndef TESTIMPL
  cdo_abort("Internal error: find_points() not implemented!");
#endif

  yac_grid_cell cell;
  cell.array_size = ncorner;
  cell.num_corners = ncorner;
  cell.edge_type = new enum yac_edge_type[ncorner];
  cell.coordinates_xyz = new double[ncorner][3];

  // enum yac_edge_type lonlat_circle_type[] = { YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
  // YAC_LON_CIRCLE_EDGE
  // };

  for (size_t k = 0; k < ncorner; ++k)
  {
    auto lon = xbounds[cell_no * ncorner + k];
    auto lat = ybounds[cell_no * ncorner + k];
    cell.edge_type[k] = YAC_GREAT_CIRCLE_EDGE;
    // cell.edge_type[k] = lonlat_circle_type[k + 1];
    gcLLtoXYZ(lon, lat, cell.coordinates_xyz[k]);
  }

  double centerCoordinates[3];
  size_t nvalues = 0;
  for (size_t k = 0; k < numIndices; ++k)
  {
    auto index1 = indices[k];
    auto lon = xvals1[index1];
    auto lat = yvals1[index1];

    gcLLtoXYZ(lon, lat, centerCoordinates);

    if (yac_point_in_cell(centerCoordinates, cell))
    {
      if (vmask[index1] == 0)
      {
        indices[nvalues] = indices[k];
        nvalues++;
      }
      if (vmask[index1] < 127) vmask[index1]++;
    }
  }

  delete[] cell.edge_type;
  delete[] cell.coordinates_xyz;

  return nvalues;
}
#endif

static size_t
find_points(Vmask &vmask, size_t cell_no, size_t ncorner, size_t numIndices, Varray<size_t> &indices, Varray<double> const &xvals1,
            Varray<double> const &yvals1, Varray<double> const &xbounds, Varray<double> const &ybounds)
{
#ifndef TESTIMPL
  cdo_abort("Internal error: find_points() not implemented!");
#endif

  Point3D centerPoint3D;
  Varray<Point> cellCornersPlaneProjection(ncorner + 1);
  Varray<Point3D> cellCorners3D(ncorner + 1);
  Varray<Point3D> cell_corners_xyz_open_cell(ncorner);
  std::vector<bool> marked_duplicate_indices(ncorner);

  set_cell_corners_3D(ncorner, &xbounds[cell_no * ncorner], &ybounds[cell_no * ncorner], cell_corners_xyz_open_cell);

  auto actualNumberOfCorners = get_actual_number_of_corners(ncorner, cell_corners_xyz_open_cell);

  auto noDuplicates = get_no_duplicates(actualNumberOfCorners, cell_corners_xyz_open_cell, marked_duplicate_indices);

  copy_unique_corners(actualNumberOfCorners, cell_corners_xyz_open_cell, marked_duplicate_indices, cellCorners3D);

  actualNumberOfCorners -= noDuplicates;

  cellCorners3D[actualNumberOfCorners] = cellCorners3D[0];

  if (actualNumberOfCorners < 3) return 0;

  auto coordinateToIgnore = find_coordinate_to_ignore(cellCorners3D);

  auto cval
      = (coordinateToIgnore == 1) ? cellCorners3D[0].X : ((coordinateToIgnore == 2) ? cellCorners3D[0].Y : cellCorners3D[0].Z);
  auto invertResult = (cval < 0.0);

  set_cell_corners_plane_projection(coordinateToIgnore, actualNumberOfCorners, cellCorners3D, cellCornersPlaneProjection);

  auto isClockwise = are_polygon_vertices_arranged_in_clockwise_order(cellCornersPlaneProjection, actualNumberOfCorners + 1);

  if (invertResult) isClockwise = !isClockwise;
  if (isClockwise) return 0;

  double centerCoordinates[3];
  size_t nvalues = 0;
  for (size_t k = 0; k < numIndices; ++k)
  {
    auto index1 = indices[k];
    auto lon = xvals1[index1];
    auto lat = yvals1[index1];

    gcLLtoXYZ(lon, lat, centerCoordinates);
    centerPoint3D.X = centerCoordinates[0];
    centerPoint3D.Y = centerCoordinates[1];
    centerPoint3D.Z = centerCoordinates[2];

    auto centerPoint2D = set_center_point_plane_projection(coordinateToIgnore, centerPoint3D);

    auto windingNumber = winding_numbers_algorithm(cellCornersPlaneProjection, actualNumberOfCorners + 1, centerPoint2D);

    if (windingNumber != 0)
    {
      if (vmask[index1] == 0)
      {
        indices[nvalues] = indices[k];
        nvalues++;
      }
      if (vmask[index1] < 127) vmask[index1]++;
    }
  }

  return nvalues;
}

static size_t
find_points_rec2d(Vmask &vmask, size_t i, size_t nlon2, size_t numIndices, Varray<size_t> &indices, Varray<double> const &xvals1,
                  Varray<double> const &yvals1, Varray<double> const &xbounds2d, Varray<double> const &ybounds2d)
{
  auto iy = i / nlon2;
  auto ix = i - iy * nlon2;

  size_t nvalues = 0;
  for (size_t k = 0; k < numIndices; ++k)
  {
    auto index1 = indices[k];
    auto y = yvals1[index1];
    if (y >= ybounds2d[2 * iy] && y < ybounds2d[2 * iy + 1])
    {
      auto x = xvals1[index1];
      if ((x >= xbounds2d[2 * ix] && x < xbounds2d[2 * ix + 1])
          || ((x - PI2) >= xbounds2d[2 * ix] && (x - PI2) < xbounds2d[2 * ix + 1]))
      {
        if (vmask[index1] == 0)
        {
          indices[nvalues] = index1;
          nvalues++;
        }
        if (vmask[index1] < 127) vmask[index1]++;
      }
    }
  }

  return nvalues;
}

static long
get_bin_num(double lon1, double lat1, size_t numBins, Varray<double> const &xbounds2d, Varray<double> const &ybounds2)
{
  if ((lon1 >= xbounds2d[0] && lon1 <= xbounds2d[1]) || ((lon1 - PI2) >= xbounds2d[0] && (lon1 - PI2) <= xbounds2d[1]))
  {
    auto latIsUp = (ybounds2[1] > ybounds2[0]);
    if (latIsUp)
    {
      for (size_t i = 0; i < numBins; ++i)
      {
        if (lat1 > ybounds2[i * 2] && lat1 <= ybounds2[i * 2 + 1]) return (long) i;
      }
    }
    else
    {
      for (size_t i = 0; i < numBins; ++i)
      {
        if (lat1 <= ybounds2[i * 2] && lat1 > ybounds2[i * 2 + 1]) return (long) i;
      }
    }
  }

  return -1;
}

static void
check_vmask(Vmask const &vmask)
{
  constexpr int maxVals = 128;
  size_t vm[maxVals] = { 0 };
  auto size = vmask.size();
  for (size_t i = 0; i < size; ++i)
    if (vmask[i] >= 0) vm[vmask[i]]++;
  for (int i = 0; i < maxVals; ++i)
    if (vm[i]) cdo_print("Number of source values: %d --> %zu/%zu", i, vm[i], size);
  size_t sum = 0;
  for (int i = 1; i < maxVals; ++i) sum += vm[i];
  cdo_print("Sum of used source values:     %zu/%zu", sum, size);
}

static Varray2D<size_t>
gen_mapdata(int gridID1, int gridID2)
{
  auto gridsize1 = gridInqSize(gridID1);
  auto gridsize2 = gridInqSize(gridID2);

  Varray2D<size_t> mapdata(gridsize2);

  Varray<double> xvals1(gridsize1), yvals1(gridsize1);
  read_coordinates(gridID1, xvals1, yvals1);

  Varray<double> xvals2(gridsize2), yvals2(gridsize2);
  read_coordinates(gridID2, xvals2, yvals2);

  auto gridtype2 = gridInqType(gridID2);
  auto grid2IsReg2d = (gridtype2 == GRID_GAUSSIAN || gridtype2 == GRID_LONLAT);

  auto nlon2 = grid2IsReg2d ? gridInqXsize(gridID2) : 0;

  int gridIDdestroy = -1;

  Varray<double> xbounds, ybounds;
  Varray<double> xbounds2d, ybounds2d;
  size_t nv = 4;
  if (grid2IsReg2d)
  {
    read_xbounds_reg2d(gridID2, xbounds2d);
    read_ybounds_reg2d(gridID2, ybounds2d);
  }
  else
  {
    auto gridID = generate_full_cell_grid(gridID2);
    if (gridID != gridID2) gridIDdestroy = gridID2 = gridID;
    nv = read_target_cell_bounds(gridID2, xbounds, ybounds);
  }

  auto isZonalGrid = (grid2IsReg2d ? is_zonal_grid(gridID2, xbounds2d, ybounds2d) : false);

  if (isZonalGrid)
  {
    Varray<size_t> zonalBinsNum(gridsize2, 0);
    Varray<long> indexMap(gridsize1, -1);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for (size_t i = 0; i < gridsize1; ++i)
    {
      auto binNum = get_bin_num(xvals1[i], yvals1[i], gridsize2, xbounds2d, ybounds2d);
      if (binNum >= 0)
      {
        indexMap[i] = binNum;
#ifdef _OPENMP
#pragma omp atomic
#endif
        zonalBinsNum[binNum]++;
      }
    }

    for (size_t i = 0; i < gridsize2; ++i) mapdata[i].reserve(zonalBinsNum[i]);

    for (size_t i = 0; i < gridsize1; ++i)
    {
      if (indexMap[i] >= 0) mapdata[indexMap[i]].push_back(i);
    }
  }
  else
  {
    cdo::timer timer;

    GridPointsearch gps;
    grid_pointsearch_create_unstruct(gps, xvals1, yvals1, true);

    if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds", timer.elapsed());

    timer.reset();

    auto maxIndex = gridsize1;
    if (gridsize1 > 1000000) maxIndex /= 4;
    Vmask vmask(gridsize1, 0);
    Varray2D<size_t> indices_2D(Threading::ompNumMaxThreads, Varray<size_t>(maxIndex));
    Varray2D<double> dist_2D(Threading::ompNumMaxThreads, Varray<double>(maxIndex));

    StatInfo statInfo;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for (size_t i = 0; i < gridsize2; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &indices = indices_2D[ompthID];
      auto &dist = dist_2D[ompthID];

      PointLonLat pointLL2{ xvals2[i], yvals2[i] };

      auto searchRadius = grid2IsReg2d ? calc_maxdist_rec2d(i, nlon2, pointLL2, xbounds2d, ybounds2d)
                                       : calc_maxdist(i, nv, pointLL2, xbounds, ybounds);

      if (searchRadius < 2.0) { searchRadius = 1.01 * std::sqrt(searchRadius); }
      else { cdo_abort("Search radius of target grid cell[%zu] > 90 degrees!", i + 1); }

      gps.set_radius(searchRadius);
      auto numIndices = grid_pointsearch_qnearest(gps, pointLL2, maxIndex, indices.data(), dist.data());
      // printf("%zu numIndices %zu\n", i+1, numIndices);

      auto nvalues = grid2IsReg2d ? find_points_rec2d(vmask, i, nlon2, numIndices, indices, xvals1, yvals1, xbounds2d, ybounds2d)
#ifdef USE_YAC
                                  : find_points_yac(vmask, i, nv, numIndices, indices, xvals1, yvals1, xbounds, ybounds);
#else
                                  : find_points(vmask, i, nv, numIndices, indices, xvals1, yvals1, xbounds, ybounds);
#endif

      if (nvalues)
      {
        mapdata[i].resize(nvalues);
        for (size_t k = 0; k < nvalues; ++k) mapdata[i][k] = indices[k];
      }

      // if (Options::cdoVerbose) printf("%zu numIndices %zu nvalues %zu  maxdist %g\n", i+1, numIndices, nvalues, maxdist);
      if (Options::cdoVerbose && Threading::ompNumMaxThreads == 1) statInfo.add(numIndices, nvalues, searchRadius);
    }

    if (Options::cdoVerbose && Threading::ompNumMaxThreads == 1) statInfo.print();
    if (Options::cdoVerbose) check_vmask(vmask);
    if (Options::cdoVerbose) cdo_print("Point search qnearest: %.2f seconds", timer.elapsed());

    if (gridIDdestroy != -1) gridDestroy(gridIDdestroy);
  }

  return mapdata;
}

template <typename T>
static T
remap_kernel(int operfunc, Varray<size_t> const &indices, size_t &numMissVals2, Field &field, Varray<T> &fieldvec,
             Varray<T> const &vec1, double mv)
{
  T missval = mv;
  T value;
  auto nvalues = indices.size();
  if (nvalues)
  {
    for (size_t k = 0; k < nvalues; ++k) { fieldvec[k] = vec1[indices[k]]; }

    field.numMissVals = 0;
    for (size_t k = 0; k < nvalues; ++k)
    {
      if (fp_is_equal(fieldvec[k], missval)) field.numMissVals++;
    }

    field.size = nvalues;
    field.missval = missval;
    value = field_function(field, operfunc);
    if (fp_is_equal(value, missval)) numMissVals2++;
  }
  else
  {
    value = missval;
    numMissVals2++;
  }

  return value;
}

static void
remap_field(const Varray2D<size_t> &mapdata, Field const &field1, Field &field2, int operfunc)
{
  std::vector<Field> fields(Threading::ompNumMaxThreads);

  auto gridsize2 = gridInqSize(field2.grid);
  cdo::timer timer;

  size_t numMissVals2 = 0;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(+ : numMissVals2)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
  {
    auto const &indices = mapdata[i];
    auto nvalues = indices.size();
    auto missval = field1.missval;

    auto ompthID = cdo_omp_get_thread_num();
    auto &field = fields[ompthID];
    field.memType = field1.memType;
    if (field1.memType == MemType::Float)
      field.vec_f.resize(nvalues);
    else
      field.vec_d.resize(nvalues);

    double rvalue = 0.0;
    if (field1.memType == MemType::Float)
      rvalue = remap_kernel(operfunc, indices, numMissVals2, field, field.vec_f, field1.vec_f, missval);
    else
      rvalue = remap_kernel(operfunc, indices, numMissVals2, field, field.vec_d, field1.vec_d, missval);

    auto func = [&](auto &v) { v[i] = rvalue; };
    field_operation(func, field2);
  }

  field2.numMissVals = numMissVals2;

  if (Options::cdoVerbose) cdo_print("Remap: %.3f seconds", timer.elapsed());
}

class Remapstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Remapstat",
    .operators = { { "remaprange", FieldFunc_Range, 0, RemapstatHelp },
                   { "remapmin", FieldFunc_Min, 0, RemapstatHelp },
                   { "remapmax", FieldFunc_Max, 0, RemapstatHelp },
                   { "remapsum", FieldFunc_Sum, 0, RemapstatHelp },
                   { "remapmean", FieldFunc_Mean, 0, RemapstatHelp },
                   { "remapavg", FieldFunc_Avg, 0, RemapstatHelp },
                   { "remapstd", FieldFunc_Std, 0, RemapstatHelp },
                   { "remapstd1", FieldFunc_Std1, 0, RemapstatHelp },
                   { "remapvar", FieldFunc_Var, 0, RemapstatHelp },
                   { "remapvar1", FieldFunc_Var1, 0, RemapstatHelp },
                   { "remapskew", FieldFunc_Skew, 0, RemapstatHelp },
                   { "remapkurt", FieldFunc_Kurt, 0, RemapstatHelp },
                   { "remapmedian", FieldFunc_Median, 0, RemapstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Remapstat> registration = RegisterEntry<Remapstat>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  VarList varList2{};

  Varray2D<size_t> mapdata{};
  int operfunc{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);

    operator_input_arg("grid description file or name");
    auto gridID2 = cdo_define_grid(cdo_operator_argv(0));
    auto gridtype2 = gridInqType(gridID2);

    {
#ifdef TESTIMPL
      auto hasProjParams = ((gridtype2 == GRID_PROJECTION) && grid_has_proj_params(gridID2));
      if (!gridProjIsSupported(gridID2) && !hasProjParams && gridtype2 != GRID_LONLAT && gridtype2 != GRID_GAUSSIAN
          && gridtype2 != GRID_CURVILINEAR && gridtype2 != GRID_UNSTRUCTURED && gridtype2 != GRID_HEALPIX)
#else
      if (gridtype2 != GRID_LONLAT && gridtype2 != GRID_GAUSSIAN)
#endif
        cdo_abort("Remapping to %s data unsupported!", gridNamePtr(gridtype2));
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids = vlistNumGrids(vlistID1);
    auto gridID1 = vlistGrid(vlistID1, 0);
    for (int index = 0; index < numGrids; ++index)
    {
      if (index > 0) cdo_abort("Too many different grids!");

      // auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);
      auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID1));
      if (!gridProjIsSupported(gridID1) && !hasProjParams && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
          && gridtype != GRID_GME && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED && gridtype != GRID_HEALPIX)
        cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    mapdata = gen_mapdata(gridID1, gridID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    Field field1, field2;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        field2.init(varList2.vars[varID]);

        remap_field(mapdata, field1, field2, operfunc);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field2);
      }
      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
