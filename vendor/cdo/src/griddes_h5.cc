/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define H5_USE_16_API

#ifdef HAVE_LIBHDF5
#include "hdf5.h"
#endif

#include <cdi.h>
#include "cdo_options.h"

#include "cdo_output.h"
#include "griddes.h"

#ifdef HAVE_LIBHDF5
static herr_t
obj_info(hid_t loc_id, const char *name, void *objname)
{
  herr_t lexist = 0;

  H5G_stat_t statbuf;
  H5Gget_objinfo(loc_id, name, false, &statbuf);

  if (std::strcmp(name, (char *) objname) == 0)
  {
    lexist = 1;

    H5G_obj_t obj_type = statbuf.type;

    switch (obj_type)
    {
      case H5G_GROUP:
        if (Options::cdoVerbose) cdo_print("HDF5 object '%s' is a group", name);
        break;
      case H5G_DATASET:
        if (Options::cdoVerbose) cdo_print("HDF5 object '%s' is a dataset", name);
        break;
      case H5G_TYPE:
        if (Options::cdoVerbose) cdo_print("HDF5 object '%s' is a named datatype", name);
        break;
      default:
        /*cdo_abort(" Unable to identify an object %s", name);*/
        break;
    }
  }

  return lexist;
}
#endif

#ifdef HAVE_LIBHDF5
static int
h5find_object(hid_t file_id, const char *name)
{
  return (int) H5Giterate(file_id, "/", nullptr, obj_info, (void *) name);
}
#endif

#ifdef HAVE_LIBHDF5
static void
fill_gridvals(size_t xsize, size_t ysize, double *xvals, double *yvals)
{
  size_t i, j, ii, jj;
  size_t index, index2;

  double xmin = -180;
  double xmax = 180;
  double ymin = -90;
  double ymax = 90;

  for (ii = 0; ii < xsize / 2; ++ii)
  {
    index2 = ysize / 2 * xsize + ii;
    if (xvals[index2] > -180 && xvals[index2] < 360)
    {
      xmin = xvals[index2];
      break;
    }
  }

  for (ii = xsize - 1; ii > xsize / 2; --ii)
  {
    index2 = ysize / 2 * xsize + ii;
    if (xvals[index2] > -180 && xvals[index2] < 360)
    {
      xmax = xvals[index2];
      break;
    }
  }
  /*
  for ( jj = 0; jj < ysize; ++jj )
    {
      index2 = jj*xsize + xsize/2;
      if ( xvals[index2] < -180 || xvals[index2] > 360 ) xvals[index2] = 0;
      index2 = jj*xsize + xsize/2-1;
      if ( xvals[index2] < -180 || xvals[index2] > 360 ) xvals[index2] = 0;
    }
  */
  for (jj = 0; jj < ysize / 2; ++jj)
  {
    index2 = jj * xsize + xsize / 2;
    if (yvals[index2] > -90 && yvals[index2] < 90)
    {
      ymax = yvals[index2];
      break;
    }
  }

  for (jj = ysize - 1; jj > ysize / 2; --jj)
  {
    index2 = jj * xsize + xsize / 2;
    if (yvals[index2] > -90 && yvals[index2] < 90)
    {
      ymin = yvals[index2];
      break;
    }
  }

  /* printf("xmin %g, xmax %g, ymin %g, ymax %g\n", xmin, xmax, ymin, ymax); */

  for (i = 0; i < xsize * ysize; ++i)
  {
    if (xvals[i] > -180 && xvals[i] < 360)
    {
      if (xvals[i] < xmin) xmin = xvals[i];
      if (xvals[i] > xmax) xmax = xvals[i];
    }

    if (yvals[i] > -90 && yvals[i] < 90)
    {
      if (yvals[i] < ymin) ymin = yvals[i];
      if (yvals[i] > ymax) ymax = yvals[i];
    }
  }

  for (j = 0; j < ysize; ++j)
    for (i = 0; i < xsize; ++i)
    {
      index = j * xsize + i;

      if (xvals[index] < -180 || xvals[index] > 360)
      {
        if (i < xsize / 2)
          xvals[index] = xmin;
        else
          xvals[index] = xmax;
        /*
        if ( j < ysize/2 )
          for ( jj = j+1; jj < ysize/2; ++jj )
            {
              index2 = jj*xsize + i;
              if ( xvals[index2] > -180 && xvals[index2] < 360 )
                {
                  xvals[index] = xvals[index2];
                  break;
                }
            }
        else
          for ( jj = j-1; jj > ysize/2; --jj )
            {
              index2 = jj*xsize + i;
              if ( xvals[index2] > -180 && xvals[index2] < 360 )
                {
                  xvals[index] = xvals[index2];
                  break;
                }
            }
        */
        /*
        if ( i < xsize/2 )
          {
            xvals[index] = xmin;
            for ( ii = i+1; ii < xsize/2; ++ii )
              {
                index2 = j*xsize + ii;
                if ( xvals[index2] > -180 && xvals[index2] < 360 )
                  {
                    xvals[index] = (xmin*(ii-i) + xvals[index2]*(i))/ii;
                    break;
                  }
              }
          }
        else
          {
            for ( ii = i-1; ii >= xsize/2; --ii )
              {
                index2 = j*xsize + ii;
                if ( xvals[index2] > -180 && xvals[index2] < 360 )
                  {
                    xvals[index] = (xmax*(i-ii) + xvals[index2]*((xsize-1)-i))/(xsize-1-ii); break;
                  }
              }
          }
        */
      }

      if (yvals[index] < -90 || yvals[index] > 90)
      {
        if (j < ysize / 2)
          yvals[index] = ymax;
        else
          yvals[index] = ymin;

        if (i < xsize / 2)
          for (ii = i + 1; ii < xsize / 2; ++ii)
          {
            index2 = j * xsize + ii;
            if (yvals[index2] > -90 && yvals[index2] < 90)
            {
              yvals[index] = yvals[index2];
              break;
            }
          }
        else
          for (ii = i - 1; ii > xsize / 2; --ii)
          {
            index2 = j * xsize + ii;
            if (yvals[index2] > -90 && yvals[index2] < 90)
            {
              yvals[index] = yvals[index2];
              break;
            }
          }
      }
    }
}

int
grid_from_h5_file(const char *gridfile)
{
  int gridID = -1;
  hid_t lon_id = -1;
  hid_t lat_id = -1;
  hid_t att_id;
  hid_t dataspace;
  hsize_t dims_out[9];  // dataset dimensions
  [[maybe_unused]] herr_t status;
  int rank;
  GridDesciption grid;

  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);

  // Open an existing file
  hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, fapl_id);

  H5Pclose(fapl_id);

  if (file_id < 0) return gridID;

  if (h5find_object(file_id, "lon") > 0 && h5find_object(file_id, "lat") > 0)
  {
    lon_id = H5Dopen(file_id, "/lon");
    lat_id = H5Dopen(file_id, "/lat");
  }
  else if (h5find_object(file_id, "Longitude") > 0 && h5find_object(file_id, "Latitude") > 0)
  {
    lon_id = H5Dopen(file_id, "/Longitude");
    lat_id = H5Dopen(file_id, "/Latitude");
  }

  if (lon_id >= 0 && lat_id >= 0)
  {
    dataspace = H5Dget_space(lon_id);  // dataspace handle
    rank = H5Sget_simple_extent_ndims(dataspace);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);

    if (rank != 2)
    {
      // if ( Options::cdoVerbose ) cdo_warning("Unexpected rank = %d!", rank);
      goto RETURN;
    }

    // check for netcdf4 attribute
    if (H5Aexists(lon_id, "DIMENSION_LIST")) goto RETURN;
    if (H5Aexists(lat_id, "DIMENSION_LIST")) goto RETURN;

    if (H5Aexists(lon_id, "bounds")) goto RETURN;
    if (H5Aexists(lat_id, "bounds")) goto RETURN;

    /*
    printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
           (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
    */

    hid_t type_id = H5Dget_type(lon_id);  // get datatype

    hid_t native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
    int ftype = 0;
    if (H5Tequal(native_type, H5T_NATIVE_SCHAR) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_UCHAR) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_SHORT) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_USHORT) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_INT) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_UINT) > 0) { ftype = 0; }
    else if (H5Tequal(native_type, H5T_NATIVE_FLOAT) > 0) { ftype = 1; }
    else if (H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0) { ftype = 1; }
    else
    {
      cdo_warning("Grid has unsupported native datatype!");
      goto RETURN;
    }
    H5Tclose(native_type);

    grid.xsize = dims_out[1];
    grid.ysize = dims_out[0];
    grid.size = grid.xsize * grid.ysize;

    grid.xvals.resize(grid.size);
    grid.yvals.resize(grid.size);

    if (ftype)
    {
      status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals.data());
      status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals.data());
    }
    else
    {
      std::vector<int> iarray(grid.size);
      status = H5Dread(lon_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray.data());
      for (size_t i = 0; i < grid.size; ++i) grid.xvals[i] = iarray[i];
      status = H5Dread(lat_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray.data());
      for (size_t i = 0; i < grid.size; ++i) grid.yvals[i] = iarray[i];
    }

    status = H5Sclose(dataspace);

    // Close the dataset
    status = H5Dclose(lon_id);
    status = H5Dclose(lat_id);

    fill_gridvals(grid.xsize, grid.ysize, grid.xvals.data(), grid.yvals.data());

    grid.type = GRID_CURVILINEAR;
    grid.datatype = CDI_DATATYPE_FLT32;

    gridID = grid_define(grid);
  }
  else if (h5find_object(file_id, "where") > 0)
  {
    double xscale = 1, yscale = 1;
    double xoffset = 0, yoffset = 0;
    hid_t grp_id;

    grp_id = H5Gopen(file_id, "/where/lon/what");
    if (grp_id >= 0)
    {
      att_id = H5Aopen_name(grp_id, "gain");
      if (att_id >= 0)
      {
        status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xscale);
        H5Aclose(att_id);
      }

      att_id = H5Aopen_name(grp_id, "offset");
      if (att_id >= 0)
      {
        status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xoffset);
        H5Aclose(att_id);
      }

      H5Gclose(grp_id);
    }

    grp_id = H5Gopen(file_id, "/where/lat/what");
    if (grp_id >= 0)
    {
      att_id = H5Aopen_name(grp_id, "gain");
      if (att_id >= 0)
      {
        status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yscale);
        H5Aclose(att_id);
      }

      att_id = H5Aopen_name(grp_id, "offset");
      if (att_id >= 0)
      {
        status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yoffset);
        H5Aclose(att_id);
      }

      H5Gclose(grp_id);
    }

    // Open an existing dataset
    lon_id = H5Dopen(file_id, "/where/lon/data");
    if (lon_id >= 0) lat_id = H5Dopen(file_id, "/where/lat/data");

    if (lon_id >= 0 && lat_id >= 0)
    {
      dataspace = H5Dget_space(lon_id);  // dataspace handle
      rank = H5Sget_simple_extent_ndims(dataspace);
      status = H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);

      if (rank != 2)
      {
        // if ( Options::cdoVerbose ) cdo_warning("Unexpected rank = %d!", rank);
        goto RETURN;
      }
      /*
      printf("\nRank: %d\nDimensions: %lu x %lu \n", rank, (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
      */

      hid_t type_id = H5Dget_type(lon_id);  // get datatype

      hid_t native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
      int ftype = 0;
      if (H5Tequal(native_type, H5T_NATIVE_SCHAR) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_UCHAR) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_SHORT) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_USHORT) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_INT) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_UINT) > 0) { ftype = 0; }
      else if (H5Tequal(native_type, H5T_NATIVE_FLOAT) > 0) { ftype = 1; }
      else if (H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0) { ftype = 1; }
      else
      {
        cdo_warning("Grid has unsupported native datatype!");
        goto RETURN;
      }
      H5Tclose(native_type);

      grid.xsize = dims_out[1];
      grid.ysize = dims_out[0];
      grid.size = grid.xsize * grid.ysize;

      grid.xvals.resize(grid.size);
      grid.yvals.resize(grid.size);

      if (ftype)
      {
        status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals.data());
        status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals.data());
      }
      else
      {
        std::vector<int> iarray(grid.size);
        status = H5Dread(lon_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray.data());
        for (size_t i = 0; i < grid.size; ++i) grid.xvals[i] = iarray[i];
        status = H5Dread(lat_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray.data());
        for (size_t i = 0; i < grid.size; ++i) grid.yvals[i] = iarray[i];
      }

      status = H5Sclose(dataspace);

      // Close the dataset.
      status = H5Dclose(lon_id);
      status = H5Dclose(lat_id);

      for (size_t i = 0; i < grid.size; ++i) grid.xvals[i] = grid.xvals[i] * xscale + xoffset;
      for (size_t i = 0; i < grid.size; ++i) grid.yvals[i] = grid.yvals[i] * yscale + yoffset;

      grid.type = GRID_CURVILINEAR;
      grid.datatype = CDI_DATATYPE_FLT32;

      gridID = grid_define(grid);
    }
  }

RETURN:

  // Close file
  if (file_id >= 0) status = H5Fclose(file_id);
  (void) status;

  if (gridID != -1 && Options::cdoVerbose) cdo_print("%s: grid created.", __func__);

  return gridID;
}
#else
int
grid_from_h5_file(const char *gridfile)
{
  (void) gridfile;
  cdo_warning("HDF5 support not compiled in!");
  return -1;
}
#endif
