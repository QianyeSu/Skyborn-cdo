/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef GRIDDES_H
#define GRIDDES_H

#include <string>
#include <cdi.h>

#include "varray.h"

constexpr double undefGridValue = 9.e20;

// clang-format off
struct  // GridDesciption
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
GridDesciption
// clang-format on
{
  std::vector<int64_t> indices;
  std::vector<int> mask;
  Varray<double> xvals;
  Varray<double> yvals;
  Varray<double> xbounds;
  Varray<double> ybounds;
  std::vector<double> area;
  std::vector<int> reducedPoints;
  char **xcvals = nullptr;
  char **ycvals = nullptr;
  double xfirst = undefGridValue, yfirst = undefGridValue;
  double xlast = undefGridValue, ylast = undefGridValue;
  double xinc = undefGridValue, yinc = undefGridValue;
  double xpole = 0.0, ypole = 0.0, angle = 0.0;  // rotated north pole
  int scanningMode = 64;
  /*
    scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
              64  = 128 * 0                + 64 *        1         + 32 * 0
              00  = 128 * 0                + 64 *        0         + 32 * 0
              96  = 128 * 0                + 64 *        1         + 32 * 1
    Default  implicit scanning mode is 64: i and j scan positively, i points are consecutive (row-major)
  */
  double a = 0.0;
  int isRotated = 0;  // true for rotated grids
  int datatype = CDI_UNDEFID;
  int type = CDI_UNDEFID;
  int ntr = 0;
  int nvertex = 0;
  size_t size = 0;
  size_t xsize = 0;
  size_t ysize = 0;
  int numLPE = 0;
  int lcomplex = 1;
  bool genBounds = false;
  int nd = 0, ni = 0, ni2 = 0, ni3 = 0;
  int number = 0, position = 0;
  unsigned char uuid[CDI_UUID_SIZE] = { 0 };
  std::string path;
  std::string xname;
  std::string xlongname;
  std::string xunits;
  std::string xdimname;
  std::string yname;
  std::string ylongname;
  std::string yunits;
  std::string ydimname;
  std::string vdimname;
  std::string projection;
  std::string healpixOrder;
  int healpixNside = 0;
  int refinementLevel = 0;
};

int grid_define(GridDesciption &grid);

int grid_from_nc_file(const char *gridfile);
int grid_from_h5_file(const char *gridfile);
int grid_from_name(std::string const &gridname);

void write_nc_grid(const char *gridfile, int gridID, int *imask);

int cdo_define_grid(std::string const &gridfile);

int grid_read(std::FILE *gfp, const char *dname);  // TODO: Find better place for this

int cdo_cdf_openread(const char *filename);
void cdo_cdf_close(int nc_file_id);
void cdo_set_grids(std::string const &gridarg);

void gaussian_latitudes_in_degrees(Varray<double> &lats, Varray<double> &lat_bounds);

#endif /* GRIDDES_H */
