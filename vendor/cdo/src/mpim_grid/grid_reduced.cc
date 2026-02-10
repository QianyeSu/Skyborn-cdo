#include "grid_reduced.h"
#include "cdo_output.h"

/*
 * grib_get_reduced_row: code from GRIB_API 1.10.4
 *
 * Description:
 *   computes the number of points within the range lon_first->lon_last and the zero based indexes ilon_first,ilon_last
 *   of the first and last point given the number of points along a parallel (pl)
 *
 */
void
grib_get_reduced_row(long pl, double lon_first, double lon_last, long *npoints, long *ilon_first, long *ilon_last)
{
  auto range = lon_last - lon_first;
  if (range < 0.0)
  {
    range += 360.0;
    lon_first -= 360.0;
  }

  // computing integer number of points and coordinates without using floating point resolution
  *npoints = (range * pl) / 360.0 + 1;
  *ilon_first = (lon_first * pl) / 360.0;
  *ilon_last = (lon_last * pl) / 360.0;

  auto irange = *ilon_last - *ilon_first + 1;

  if (irange != *npoints)
  {
    if (irange > *npoints)
    {
      // checking if the first point is out of range
      auto dlon_first = ((*ilon_first) * 360.0) / pl;
      if (dlon_first < lon_first)
      {
        (*ilon_first)++;
        irange--;
      }

      // checking if the last point is out of range
      auto dlon_last = ((*ilon_last) * 360.0) / pl;
      if (dlon_last > lon_last)
      {
        (*ilon_last)--;
        irange--;
      }
    }
    else
    {
      int ok = 0;
      // checking if the point before the first is in the range
      auto dlon_first = ((*ilon_first - 1) * 360.0) / pl;
      if (dlon_first > lon_first)
      {
        (*ilon_first)--;
        irange++;
        ok = 1;
      }

      // checking if the point after the last is in the range
      auto dlon_last = ((*ilon_last + 1) * 360.0) / pl;
      if (dlon_last < lon_last)
      {
        (*ilon_last)++;
        irange++;
        ok = 1;
      }

      // if neither of the two are triggered then npoints is too large
      if (!ok) (*npoints)--;
    }

    //   assert(*npoints==irange);
  }
  else
  {
    // checking if the first point is out of range
    auto dlon_first = ((*ilon_first) * 360.0) / pl;
    if (dlon_first < lon_first)
    {
      (*ilon_first)++;
      (*ilon_last)++;
    }
  }

  if (*ilon_first < 0) *ilon_first += pl;
}

int
qu2reg_subarea(size_t gridsize, int np, double xfirst, double xlast, double *array, int *reducedPoints, int ny, double missval,
               int *iret, int lmiss, int lperio, int lveggy)
{
  // sub area (longitudes)
  long ilon_firstx;
  long ilon_first, ilon_last;
  int i, j;
  long row_count;
  int rlon;
  int np4 = np * 4;
  size_t size = 0;
  int wlen;
  int ii;

  if (np <= 0) cdo_abort("Number of values between pole and equator missing!");

  grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_firstx, &ilon_last);
  int nx = row_count;
  // printf("nx %d  %ld %ld lon1 %g lon2 %g\n", nx, ilon_firstx, ilon_last,
  // (ilon_firstx*360.)/np4, (ilon_last*360.)/np4);

  // int nwork = 0;
  // for (j = 0; j < ny; ++j) nwork += reducedPoints[j];

  double **pwork = (double **) malloc(ny * sizeof(double *));
  double *work = (double *) malloc(ny * np4 * sizeof(double));
  wlen = 0;
  pwork[0] = work;
  for (j = 1; j < ny; ++j)
  {
    wlen += reducedPoints[j - 1];
    pwork[j] = work + wlen;
  }
  // printf(" ny, np4, nwork %d %d %d wlen %d\n", ny, np4, nwork, wlen);

  for (j = 0; j < ny; ++j)
  {
    rlon = reducedPoints[j];
    for (i = 0; i < rlon; ++i) pwork[j][i] = missval;
  }

  double *parray = array;
  for (j = 0; j < ny; ++j)
  {
    rlon = reducedPoints[j];
    row_count = 0;
    grib_get_reduced_row(rlon, xfirst, xlast, &row_count, &ilon_first, &ilon_last);
    // printf("j %d xfirst %g xlast %g reducedPoints %d %ld %ld %ld %g %g\n", j,
    // xfirst, xlast, rlon, row_count, ilon_first, ilon_last,
    // (ilon_first*360.)/rlon, (ilon_last*360.)/rlon);

    for (i = ilon_first; i < (ilon_first + row_count); ++i)
    {
      ii = i;
      if (ii >= rlon) ii -= rlon;
      pwork[j][ii] = *parray;
      parray++;
    }
    size += row_count;
  }

  if (gridsize != size) cdo_abort("gridsize1 inconsistent! (gridsize=%zu found=%zu)", gridsize, size);

  qu2reg3_double(work, reducedPoints, ny, np4, missval, iret, lmiss, lperio, lveggy);

  wlen = 0;
  pwork[0] = work;
  for (j = 1; j < ny; ++j)
  {
    wlen += np4;
    pwork[j] = work + wlen;
  }

  // printf("nx, ilon_firstx %d %ld\n", nx, ilon_firstx);
  parray = array;
  for (j = 0; j < ny; ++j)
  {
    for (i = ilon_firstx; i < (ilon_firstx + nx); ++i)
    {
      ii = i;
      if (ii >= np4) ii -= np4;
      *parray = pwork[j][ii];
      parray++;
    }
  }

  free(work);
  free(pwork);

  return nx;
}
