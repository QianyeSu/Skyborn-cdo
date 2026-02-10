/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include <cdi.h>

#include <mpim_grid.h>
#include "cdo_output.h"
#include "interpol.h"
#include "field.h"
#include "progress.h"
#include "cdo_omp.h"

/**
* Find the interval i-1 .. i in which an element x fits and return i, the
* bigger one of the interval borders or x itself if it is an interval border.
*
* If no interval can be found return the length of the array.

* @param *array ascending or descending sorted list
* @param nelem  length of the sorted list
* @param x      the element to find a position for
*/
static long
find_element(double x, long nelem, Varray<double> const &v)
{
  long ii;
  long mid = 0;
  long first = 1;
  long last = nelem;

  if (v[0] < v[nelem - 1])  // ascending order
  {
    // return the length of the array if x is out of bounds
    if (x < v[0] || x > v[nelem - 1]) return nelem;

    // search for the interval in which x fits
    // implementation: binary search algorithm
    for (ii = 1; ii < nelem; ++ii)
    {
      // binary search: divide search room in the middle
      mid = (first + last) >> 1;

      // return the bigger interval border of the interval in which x fits
      if (!(x < v[mid - 1] || x > v[mid])) break;

      // binary search: ignore half of the search room
      if (x > v[mid])
        first = mid;
      else
        last = mid;
    }
  }
  else
  {
    // return the length of the array if x is out of bounds
    if (x < v[nelem - 1] || x > v[0]) return nelem;

    // search for the interval in which x fits
    // implementation: binary search algorithm
    for (ii = 1; ii < nelem; ++ii)
    {
      // binary search: divide search room in the middle
      mid = (first + last) >> 1;

      // return the bigger interval border of the interval in which x fits
      if (!(x < v[mid] || x > v[mid - 1])) break;

      // binary search: ignore half of the search room
      if (x < v[mid])
        first = mid;
      else
        last = mid;
    }
  }

  if (mid > 1 && is_equal(x, v[mid - 1])) mid--;

  return mid;
}

bool
rect_grid_search(size_t &ii, size_t &jj, double x, double y, size_t nxm, size_t nym, Varray<double> const &xm,
                 Varray<double> const &ym)
{
  constexpr double rtol = 1.e-12;
  auto pointFound{ false };

  jj = find_element(y, nym, ym);
  if (jj >= nym && std::fabs(ym[0] - y) < rtol) jj = 1;  // fix rounding errors

  if (jj < nym)
  {
    ii = find_element(x, nxm, xm);
    if (ii >= nxm && std::fabs(xm[0] - x) < rtol) ii = 1;  // fix rounding errors

    if (ii < nxm) pointFound = true;
  }

  return pointFound;
}

bool
rect_grid_search2(long &imin, long &imax, double xmin, double xmax, long nxm, Varray<double> const &xm)
{
  auto pointFound{ false };
  imin = nxm;
  imax = -1;

  auto isAscend = (xm[0] < xm[nxm - 1]);

  auto i1 = find_element(xmin, nxm, xm);
  auto i2 = find_element(xmax, nxm, xm);

  if (i1 > 0 && i1 < nxm)
  {
    pointFound = true;

    if (isAscend)
    {
      if (i1 > 1 && xmin <= xm[i1 - 1]) i1--;
      imin = i1 - 1;
      imax = i1 - 1;
    }
    else
    {
      if (i1 < nxm - 1 && xmin <= xm[i1]) i1++;
      imin = i1 - 1;
      imax = i1 - 1;
    }
  }

  if (i2 > 0 && i2 < nxm)
  {
    pointFound = true;

    if (isAscend)
    {
      if (i2 < nxm - 1 && xmax >= xm[i2]) i2++;
      imax = i2 - 1;
      if (imin == nxm) imin = imax;
    }
    else
    {
      if (i2 > 1 && xmax >= xm[i2 - 1]) i2--;
      imin = i2 - 1;
      if (imax == -1) imax = imin;
    }
  }

  return pointFound;
}
/*
double
intlinarr2p(long nxm, long nym, double **fieldm, Varray<double> const &xm, Varray<double> const &ym, double x, double y)
{
  long ii, jj;
  double value = 0;

  for (jj = 1; jj < nym; ++jj)
    if (y >= std::min(ym[jj - 1], ym[jj]) && y <= std::max(ym[jj - 1], ym[jj])) break;

  for (ii = 1; ii < nxm; ++ii)
    if (x >= xm[ii - 1] && x <= xm[ii]) break;

  if (jj < nym && ii < nxm)
    {
      value = fieldm[jj - 1][ii - 1] * (x - xm[ii]) * (y - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj - 1][ii] * (x - xm[ii - 1]) * (y - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj][ii - 1] * (x - xm[ii]) * (y - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]))
              + fieldm[jj][ii] * (x - xm[ii - 1]) * (y - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
    }

  return value;
}
*/

template <typename T>
static inline T
bilinear_remap(Varray<T> const &srcArray, const double (&wgt)[4], const size_t (&ind)[4])
{
  return srcArray[ind[0]] * wgt[0] + srcArray[ind[1]] * wgt[1] + srcArray[ind[2]] * wgt[2] + srcArray[ind[3]] * wgt[3];
}

template <typename T1, typename T2>
static void
intlinarr2(double mv, int lonIsCircular, size_t nxm, size_t nym, Varray<T1> const &varray1, Varray<double> const &xm,
           Varray<double> const &ym, size_t gridsize2, Varray<T2> &varray2, Varray<double> const &x, Varray<double> const &y)
{
  T1 missval = mv;
  auto nlon1 = nxm;
  std::atomic<size_t> atomicCount{ 0 };

  if (lonIsCircular) nlon1--;
  size_t gridsize1 = nlon1 * nym;
  Vmask gridMask1(gridsize1);
  for (size_t i = 0; i < gridsize1; ++i) gridMask1[i] = !fp_is_equal(varray1[i], missval);

  cdo::Progress progress;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
  {
    size_t srcIndices[4];  // indices for the four source points

    varray2[i] = missval;

    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && gridsize2 > progressMinSize) progress.update((double) atomicCount / gridsize2);

    auto plon = x[i];
    auto plat = y[i];
    size_t ii, jj;
    auto lfound = rect_grid_search(ii, jj, plon, y[i], nxm, nym, xm, ym);
    if (lfound)
    {
      size_t iix = (lonIsCircular && ii == (nxm - 1)) ? 0 : ii;
      srcIndices[0] = (jj - 1) * nlon1 + (ii - 1);
      srcIndices[1] = (jj - 1) * nlon1 + (iix);
      srcIndices[2] = (jj) *nlon1 + (ii - 1);
      srcIndices[3] = (jj) *nlon1 + (iix);

      // Check to see if points are missing values
      for (int n = 0; n < 4; ++n)
        if (!gridMask1[srcIndices[n]]) lfound = 0;
    }

    if (lfound)
    {
      double weights[4];
      weights[0] = (plon - xm[ii]) * (plat - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]));
      weights[1] = (plon - xm[ii - 1]) * (plat - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]));
      weights[3] = (plon - xm[ii - 1]) * (plat - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
      weights[2] = (plon - xm[ii]) * (plat - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]));

      varray2[i] = bilinear_remap(varray1, weights, srcIndices);
    }
  }
}

void
intlinarr(long nxm, double *ym, double *xm, int nx, double *y, double *x)
{
  /*
    intlinarr - lineare interpolation over 1D array

    Uwe Schulzweida  04/05/1995
  */
  for (long jj = 1; jj < nxm; ++jj)
    for (long j = 0; j < nx; ++j)
      if (x[j] >= xm[jj - 1] && x[j] <= xm[jj]) y[j] = intlin(x[j], ym[jj - 1], xm[jj - 1], ym[jj], xm[jj]);
}

void
intgrid_bil(Field const &field1, Field &field2)
{
  auto gridID1 = field1.grid;
  auto gridID2 = field2.grid;
  if (gridID1 == -1) cdo_abort("Source grid undefined!");
  if (gridID2 == -1) cdo_abort("Target grid undefined!");

  auto missval = field1.missval;

  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  int lonIsCircular = 0;

  bool isGeorefGrid = true;
  if (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2)) isGeorefGrid = false;

  if (isGeorefGrid)
  {
    if (!gridHasCoordinates(gridID1)) cdo_abort("Source grid has no coordinate values!");

    lonIsCircular = gridIsCircular(gridID1);
    if (lonIsCircular) nlon1 += 1;
  }

  Varray<double> lons1(nlon1), lats1(nlat1);
  gridInqXvals(gridID1, lons1.data());
  gridInqYvals(gridID1, lats1.data());

  if (isGeorefGrid)
  {
    if (lonIsCircular) lons1[nlon1 - 1] = 0;

    static bool doCheck = true;
    if (doCheck)
    {
      doCheck = false;
      check_longitude_range(lons1, "center", cdo_grid_get_units(gridID1, CDI_XAXIS, "grid1 center lon"));
      check_latitude_range(lats1, "center", cdo_grid_get_units(gridID1, CDI_YAXIS, "grid1 center lat"));
    }

    cdo_grid_to_radian(gridID1, CDI_XAXIS, lons1, "grid1 center lon");
    cdo_grid_to_radian(gridID1, CDI_YAXIS, lats1, "grid1 center lat");

    if (lonIsCircular) lons1[nlon1 - 1] = lons1[0] + 2 * M_PI;
  }

  auto xsize2 = gridInqXsize(gridID2);
  auto ysize2 = gridInqYsize(gridID2);

  if (isGeorefGrid)
  {
    gridID2 = generate_full_point_grid(gridID2);
    if (!gridHasCoordinates(gridID2)) cdo_abort("Target cell center coordinates missing!");
  }

  auto gridsize2 = gridInqSize(gridID2);
  Varray<double> xvals2(gridsize2), yvals2(gridsize2);

  if (isGeorefGrid)
  {
    gridInqXvals(gridID2, xvals2.data());
    gridInqYvals(gridID2, yvals2.data());

    cdo_grid_to_radian(gridID2, CDI_XAXIS, xvals2, "grid2 center lon");
    cdo_grid_to_radian(gridID2, CDI_YAXIS, yvals2, "grid2 center lat");

    for (size_t i = 0; i < gridsize2; ++i)
    {
      if (xvals2[i] < lons1[0]) xvals2[i] += 2 * M_PI;
      if (xvals2[i] > lons1[nlon1 - 1]) xvals2[i] -= 2 * M_PI;
    }
  }
  else
  {
    Varray<double> xcoord(xsize2), ycoord(ysize2);
    gridInqXvals(gridID2, xcoord.data());
    gridInqYvals(gridID2, ycoord.data());

    for (size_t j = 0; j < ysize2; ++j)
      for (size_t i = 0; i < xsize2; ++i)
      {
        xvals2[j * xsize2 + i] = xcoord[i];
        yvals2[j * xsize2 + i] = ycoord[j];
      }
  }

  if (field2.grid != gridID2) gridDestroy(gridID2);

  auto func = [&](auto &v1, auto &v2)
  { intlinarr2(missval, lonIsCircular, nlon1, nlat1, v1, lons1, lats1, gridsize2, v2, xvals2, yvals2); };
  field_operation2(func, field1, field2);

  field_num_mv(field2);
}
