/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "grid_read_pingo.h"

#include <cdi.h>

#include "griddes.h"
#include <mpim_grid.h>
#include "cdo_output.h"
#include "compare.h"
#include "gaussian_latitudes.h"

static void
skip_nondigit_lines(std::FILE *gfp)
{
  int c;

  if (feof(gfp)) return;

  while (1)
  {
    do c = fgetc(gfp);
    while ((std::isspace(c) || c == ',') && c != EOF);

    if (c == EOF || std::isdigit(c) || c == '.' || c == '+' || c == '-' || c == 'N')
      break;  // N: for NaN
    else
      while (c != '\n' && c != '\r' && c != EOF) c = fgetc(gfp);
  }

  ungetc(c, gfp);
}

static int
input_ival(std::FILE *gfp, int &ival)
{
  skip_nondigit_lines(gfp);

  if (feof(gfp)) return 0;

  ival = 0;

  return std::fscanf(gfp, "%d", &ival);
}

size_t
input_darray(std::FILE *gfp, size_t n_values, Varray<double> &array)
{
  if (n_values <= 0) return 0;

  size_t read_items = 0;
  for (size_t i = 0; i < n_values; ++i)
  {
    skip_nondigit_lines(gfp);

    if (feof(gfp)) break;

    read_items += std::fscanf(gfp, "%lg", &array[i]);

    if (feof(gfp)) break;
  }

  return read_items;
}

int
grid_read_pingo(std::FILE *gfp)
{
  int gridID = -1;

  int nlon, nlat;
  if (!input_ival(gfp, nlon)) return gridID;
  if (!input_ival(gfp, nlat)) return gridID;

  GridDesciption grid;

  if (nlon > 0 && nlon < 99999 && nlat > 0 && nlat < 99999)
  {
    int i;

    grid.xsize = nlon;
    grid.ysize = nlat;

    grid.xvals.resize(grid.xsize);
    grid.yvals.resize(grid.ysize);

    if (!input_ival(gfp, nlon)) return gridID;
    if (nlon == 2)
    {
      if (input_darray(gfp, 2, grid.xvals) != 2) return gridID;
      grid.xvals[1] -= 360 * std::floor((grid.xvals[1] - grid.xvals[0]) / 360);

      if (grid.xsize > 1)
        if (is_equal(grid.xvals[0], grid.xvals[1])) grid.xvals[1] += 360;

      for (i = 0; i < (int) grid.xsize; ++i) grid.xvals[i] = grid.xvals[0] + i * (grid.xvals[1] - grid.xvals[0]);
    }
    else if (nlon == (int) grid.xsize)
    {
      if (input_darray(gfp, nlon, grid.xvals) != (size_t) nlon) return gridID;
      for (i = 0; i < nlon - 1; ++i)
        if (grid.xvals[i + 1] <= grid.xvals[i]) break;

      for (i++; i < nlon; ++i)
      {
        grid.xvals[i] += 360;
        if (i < nlon - 1 && grid.xvals[i + 1] + 360 <= grid.xvals[i])
        {
          cdo_print("Longitudes are not in ascending order!");
          return gridID;
        }
      }
    }
    else
      return gridID;

    if (!input_ival(gfp, nlat)) return gridID;
    if (nlat == 2)
    {
      if (input_darray(gfp, 2, grid.yvals) != 2) return gridID;
      for (i = 0; i < (int) grid.ysize; ++i) grid.yvals[i] = grid.yvals[0] + i * (grid.yvals[1] - grid.yvals[0]);
    }
    else if (nlat == (int) grid.ysize)
    {
      if (input_darray(gfp, nlat, grid.yvals) != (size_t) nlat) return gridID;
    }
    else
      return gridID;

    if (grid.yvals[0] > 90.001 || grid.yvals[nlat - 1] > 90.001 || grid.yvals[0] < -90.001 || grid.yvals[nlat - 1] < -90.001)
    {
      cdo_print("Latitudes must be between 90 and -90!");
      return gridID;
    }

    for (i = 0; i < nlat - 1; ++i)
      if (is_equal(grid.yvals[i + 1], grid.yvals[i])
          || (i < nlat - 2 && ((grid.yvals[i + 1] > grid.yvals[i]) != (grid.yvals[i + 2] > grid.yvals[i + 1]))))
      {
        cdo_print("Latitudes must be in descending or ascending order!");
        return gridID;
      }

    grid.type = is_gaussian_latitudes(nlat, grid.yvals.data()) ? GRID_GAUSSIAN : GRID_LONLAT;
  }

  if (grid.type != CDI_UNDEFID) gridID = grid_define(grid);

  return gridID;
}
