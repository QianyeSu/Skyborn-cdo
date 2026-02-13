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
#include <cstdlib>  // realloc()

#include "cdo_options.h"
#include "varray.h"
#include "process_int.h"
#include "cdo_default_values.h"
#include <mpim_grid.h>

#define MAX_DSETS 1024

namespace
{
struct dset_obj_t
{
  char *name;
  char *description;
  char *units;
  char *title;
  char *time;
  int dtype;
  int nx;
  int ny;
  int nz;
  int nt;
  size_t gridsize;
  bool haveScalefactor;
  bool haveAddoffset;
  bool haveMissvals;
  double scalefactor;
  double addoffset;
  double missval;
  double *array;
};

struct datasets_t
{
  int numSets;
  int mergelevel;
  int lgeoloc;
  int lregion;
  int lprojtype;
  int lmetadata;
  dset_obj_t obj[MAX_DSETS];
};
}  // namespace

#ifdef HAVE_LIBHDF5
static void
print_filter(hid_t dset_id, char *varname)
{
  hid_t plist;
  unsigned int flags;
  int idx;
  unsigned int cd_values;
  int nfilter;
  size_t cd_nelmts = 1;
  size_t pnamelen = 64;
  char pname[64];

  // get filter
  plist = H5Dget_create_plist(dset_id);
  nfilter = H5Pget_nfilters(plist);

  for (idx = 0; idx < nfilter; idx++)
  {
    H5Pget_filter(plist, idx, &flags, &cd_nelmts, &cd_values, pnamelen, pname);
    cdo_print("Dataset %s: filter %d =  %s", varname, idx + 1, pname);
  }

  H5Pclose(plist);
}

static void
get_grid_info(double c0, double re, int *nrxp, int *nryp, double *r0p, double *s0p, double *cp)
{
  constexpr double pi = M_PI;

  double git = 2. * pi * re * std::cos(pi / 6.) / c0;
  // number of longitude pixels
  int nrx = 2 * (int) std::lround(0.5 * git);

  // central index in longitude
  double r0 = nrx / 2 + 0.5;

  // resolution in km
  double c = 2. * pi * re * std::cos(30. * pi / 180.) / nrx;

  double phi = pi / 2.;
  double s90 = re / c * std::sin(phi) / std::cos(30. * pi / 180.);

  int nry = (int) std::floor(s90);
  // central index in latitude
  double s0 = nry + 0.5;
  // number of latitude pixels
  nry = 2 * nry;

  *nrxp = nrx;
  *nryp = nry;
  *r0p = r0;
  *s0p = s0;
  *cp = c;
}

static double
det_lon_atovs(double r, double r0, double lts, double c, double re)
{
  constexpr double pi = M_PI;

  double xla = (r - r0) * c / re / std::cos(lts * pi / 180.); /* longitude */
  xla = 180. * xla / pi;

  return xla;
}

static double
det_lat_atovs(double s, double s0, double lts, double c, double re)
{
  constexpr double pi = M_PI;

  double siphi = (s - s0) * c * std::cos(lts * pi / 180.) / re;
  double phi = 180. * std::asin(siphi) / pi; /* latitude */

  return phi;
}

static int
defLonLatGrid(int nx, int ny, double c0, double lts, double re)
{
  int gridID;
  int nrx, nry;
  double c;
  double r0, s0;

  get_grid_info(c0, re, &nrx, &nry, &r0, &s0, &c);

  if (nx != nrx || ny != nry)
  {
    printf("nrx=%d nry=%d\n", nrx, nry);
    return -1;
  }

  Varray<double> xvals(nx), yvals(ny);
  Varray<double> xbounds(nx * 2), ybounds(ny * 2);

  for (int i = 0; i < nx; ++i)
  {
    double r = i + 1;
    auto xla = det_lon_atovs(r, r0, lts, c, re);
    xvals[i] = xla;
    xla = det_lon_atovs(r - 0.5, r0, lts, c, re);
    xbounds[2 * i] = xla;
    xla = det_lon_atovs(r + 0.5, r0, lts, c, re);
    xbounds[2 * i + 1] = xla;
    /* printf("xla[%d]=%g\n", i, xla); */
  }

  for (int i = 0; i < ny; ++i)
  {
    double s = (nry - i - 1) + 1;
    auto phi = det_lat_atovs(s, s0, lts, c, re);
    yvals[i] = phi;
    phi = det_lat_atovs(s - 0.5, s0, lts, c, re);
    ybounds[2 * i] = phi;
    phi = det_lat_atovs(s + 0.5, s0, lts, c, re);
    ybounds[2 * i + 1] = phi;
    /* printf("phi[%d]=%g\n", i, phi); */
  }

  gridID = gridCreate(GRID_LONLAT, nx * ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals.data());
  gridDefYvals(gridID, yvals.data());
  /*
  gridDefXbounds(gridID, xbounds.data());
  gridDefYbounds(gridID, ybounds.data());
  */

  return gridID;
}

static int
defSinusoidalGrid(int nx, int ny, double xmin, double ymax, double dx, double dy, double p1, double p2, double p3, double p4)
{
  (void) (p1);  // unused
  (void) (p2);  // unused
  (void) (p3);  // unused
  (void) (p4);  // unused
  Varray<double> xvals(nx), yvals(ny);

  for (int i = 0; i < nx; ++i) xvals[i] = xmin + i * dx + dx / 2;
  for (int i = 0; i < ny; ++i) yvals[i] = ymax - i * dy - dy / 2;

  auto gridID = gridCreate(GRID_PROJECTION, nx * ny);

  grid_def_params_sinu(gridID);

  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals.data());
  gridDefYvals(gridID, yvals.data());

  return gridID;
}

static int
defLaeaGrid(int nx, int ny, double xmin, double ymax, double dx, double dy, double a, double lon0, double lat0)
{
  Varray<double> xvals(nx), yvals(ny);

  for (int i = 0; i < nx; ++i) xvals[i] = xmin + i * dx + dx / 2;
  for (int i = 0; i < ny; ++i) yvals[i] = ymax - i * dy - dy / 2;

  auto gridID = gridCreate(GRID_PROJECTION, nx * ny);

  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals.data());
  gridDefYvals(gridID, yvals.data());

  grid_def_params_laea(gridID, a, lon0, lat0);

  return gridID;
}

static int
scan_pcs_def(char *pcs_def, char proj[128], double *a, double *lon0, double *lat0)
{
  char *pcs[64];
  int npcs = 0;
  int len;
  int nfound = 0;

  std::strcpy(proj, "unknown");
  *a = 1;
  *lon0 = 0;
  *lat0 = 0;

  pcs[npcs++] = &pcs_def[0];
  len = (int) std::strlen(pcs_def);
  for (int i = 0; i < len; ++i)
    if (pcs_def[i] == ',' && npcs < 64)
    {
      pcs_def[i] = 0;
      pcs[npcs++] = &pcs_def[i + 1];
    }

  for (int i = 0; i < npcs; ++i)
  {
    if (std::memcmp(pcs[i], "proj=", 5) == 0)
    {
      pcs[i] += 5;
      std::strcpy(proj, pcs[i]);
      nfound++;
    }
    else if (std::memcmp(pcs[i], "a=", 2) == 0)
    {
      pcs[i] += 2;
      *a = atof(pcs[i]);
      nfound++;
    }
    else if (std::memcmp(pcs[i], "lon_0=", 6) == 0)
    {
      pcs[i] += 6;
      *lon0 = atof(pcs[i]);
      nfound++;
    }
    else if (std::memcmp(pcs[i], "lat_0=", 6) == 0)
    {
      pcs[i] += 6;
      *lat0 = atof(pcs[i]);
      nfound++;
    }
  }

  return nfound;
}

static int
read_geolocation(hid_t loc_id, int nx, int ny, int lprojtype)
{
  int gridID = -1;
  hid_t grp_id;
  hid_t proj_id, region_id;
  hid_t proj_tid, region_tid;
  hid_t str_tid, fltarr_tid;
  hid_t ptype_id;
  hsize_t dims;
  int xsize, ysize;
  struct proj_t
  {
    char name[64] = { 0 };
    char ellipsoid[64] = { 0 };
    float parameter[10] = { 0.0 };
  };
  struct region_t
  {
    float xmin = 0.0;
    float xmax = 0.0;
    float ymin = 0.0;
    float ymax = 0.0;
    float dx = 0.0;
    float dy = 0.0;
  };

  proj_t proj;
  region_t region;
  char *projection_name = nullptr;

  if (Options::cdoVerbose) cdo_print("Read geolocation:");

  if (lprojtype)
  {
    ptype_id = H5Topen(loc_id, "ProjType");
    if (ptype_id >= 0)
    {
      projection_name = H5Tget_member_name(ptype_id, 0);
      H5Tclose(ptype_id);
    }
  }

  str_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_tid, 64);
  dims = 10;
  fltarr_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims, nullptr);

  proj_tid = H5Tcreate(H5T_COMPOUND, sizeof(proj_t));
  if (projection_name)
    H5Tinsert(proj_tid, projection_name, HOFFSET(proj_t, name), str_tid);
  else
    H5Tinsert(proj_tid, "Projection name", HOFFSET(proj_t, name), str_tid);
  H5Tinsert(proj_tid, "Reference ellipsoid", HOFFSET(proj_t, ellipsoid), str_tid);
  H5Tinsert(proj_tid, "Projection parameter", HOFFSET(proj_t, parameter), fltarr_tid);

  // if (projection_name) std::free(projection_name); failed with: pointer being freed was not allocated

  grp_id = H5Gopen(loc_id, "Geolocation");

  proj_id = H5Dopen(grp_id, "Projection");
  if (proj_id < 0) proj_id = H5Dopen(grp_id, "projection");
  /*
  {
    hid_t tid;
    int nmem;
    int im;

    tid = H5Dget_type(proj_id);
    nmem = H5Tget_nmembers(tid);
    for ( im = 0; im < nmem; ++im )
      {
        printf("%d %s\n", im, H5Tget_member_name(tid, im));
      }
  }
  */
  if (proj_id >= 0) H5Dread(proj_id, proj_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &proj);

  H5Dclose(proj_id);
  H5Tclose(proj_tid);
  H5Tclose(str_tid);
  H5Tclose(fltarr_tid);

  if (Options::cdoVerbose)
    cdo_print("  Projection: name=%s\n\t\t\tellipsoid=%s\n\t\t\tparameter=%g %g %g %g %g %g", proj.name, proj.ellipsoid,
              proj.parameter[0], proj.parameter[1], proj.parameter[2], proj.parameter[3], proj.parameter[4], proj.parameter[5]);

  region_tid = H5Tcreate(H5T_COMPOUND, sizeof(region_t));
  H5Tinsert(region_tid, "xmin", HOFFSET(region_t, xmin), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "xmax", HOFFSET(region_t, xmax), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "ymin", HOFFSET(region_t, ymin), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "ymax", HOFFSET(region_t, ymax), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "dx", HOFFSET(region_t, dx), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "dy", HOFFSET(region_t, dy), H5T_NATIVE_FLOAT);

  region_id = H5Dopen(grp_id, "Region");
  if (region_id < 0) region_id = H5Dopen(grp_id, "region");

  if (region_id >= 0) H5Dread(region_id, region_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &region);

  H5Dclose(region_id);
  H5Tclose(region_tid);

  if (region.xmin > region.xmax)
  {
    double xmin = region.xmin;
    region.xmin = region.xmax;
    region.xmax = xmin;
    if (Options::cdoVerbose) cdo_print("  Swap xmin/xmax");
  }

  if (Options::cdoVerbose)
    cdo_print("  Region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g", region.xmin, region.xmax, region.ymin, region.ymax,
              region.dx, region.dy);

  H5Gclose(grp_id);

  /* check region */
  xsize = (int) std::lround((region.xmax - region.xmin) / region.dx);
  ysize = (int) std::lround((region.ymax - region.ymin) / region.dy);

  if (Options::cdoVerbose) cdo_print("  Size: xsize=%d  ysize=%d", xsize, ysize);

  /* some CM-SAF files have incorrect entries for some metadata. */
  /* these are corrected in the following sections. */
  /* in case of questions on this, contact frank.kaspar@dwd.de */
  if (cdo_cmpstr(proj.ellipsoid, "WSG-84")) std::strcpy(proj.ellipsoid, "WGS-84");

  if ((int) region.xmin == -8887500 && (int) region.xmax == -8887500 && (int) region.ymin == 8887500 && (int) region.ymax == 8887500
      && (int) region.dx == 15000 && (int) region.dy == 15000)
  {
    region.xmax = 8887500.0;
    region.ymin = -8887500.0;
    if (Options::cdoVerbose)
      cdo_print("  Corrected region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g", region.xmin, region.xmax, region.ymin,
                region.ymax, region.dx, region.dy);

    xsize = (int) std::lround((region.xmax - region.xmin) / region.dx);
    ysize = (int) std::lround((region.ymax - region.ymin) / region.dy);
    if (Options::cdoVerbose) cdo_print("  Corrected size: xsize=%d  ysize=%d", xsize, ysize);
  }

  if (nx == 298 && ny == 371 && (int) region.xmin == -6709222 && (int) region.xmax == 6709222 && (int) region.ymin == -6664078
      && (int) region.ymax == 9984898 && (int) region.dx == 45000 && (int) region.dy == 45000)
  {
    region.xmin = -6705000;
    region.xmax = 6705000;
    region.ymin = -6705000;
    region.ymax = 9990000;
    cdo_print("  Corrected region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g", region.xmin, region.xmax, region.ymin, region.ymax,
              region.dx, region.dy);

    xsize = (int) std::lround((region.xmax - region.xmin) / region.dx);
    ysize = (int) std::lround((region.ymax - region.ymin) / region.dy);
    if (Options::cdoVerbose) cdo_print("  Corrected size: xsize=%d  ysize=%d", xsize, ysize);
  }

  if (!cdo_cmpstr(proj.name, "sinusoidal")
      && ((nx == xsize && ny == ysize && (int) region.xmin == -8887500 && (int) region.xmax == 8887500
           && (int) region.ymin == -8887500 && (int) region.ymax == 8887500 && (int) region.dx == 15000 && (int) region.dy == 15000)
          || (nx == xsize && ny == ysize && (int) region.xmin == -5827500 && (int) region.xmax == 5827500
              && (int) region.ymin == 3307500 && (int) region.ymax == 8887500 && (int) region.dx == 15000
              && (int) region.dy == 15000)
          || (nx == xsize && ny == ysize && (int) region.xmin == -5827500 && (int) region.xmax == 5827500
              && (int) region.ymin == 3307500 && (int) region.ymax == 8887500 && (int) region.dx == 45000
              && (int) region.dy == 45000)
          || (nx == xsize && ny == ysize && (int) region.xmin == -5827500 && (int) region.xmax == 5827500
              && (int) region.ymin == 3307500 && (int) region.ymax == 8887500 && (int) region.dx == 3000 && (int) region.dy == 3000)
          || (nx == 298 && ny == 371 && (int) region.xmin == -6709222 && (int) region.xmax == 6709222
              && (int) region.ymin == -6664078 && (int) region.ymax == 9984898 && (int) region.dx == 45000
              && (int) region.dy == 45000)
          || (nx == xsize && ny == ysize && (int) region.xmin == -6705000 && (int) region.xmax == 6705000
              && (int) region.ymin == -6705000 && (int) region.ymax == 9990000 && (int) region.dx == 45000
              && (int) region.dy == 45000)))
  {
    if (Options::cdoVerbose) cdo_print("Replacing incorrect projection parameters for sinusoidal products:");
    std::strcpy(proj.ellipsoid, "WGS-84");
    std::strcpy(proj.name, "sinusoidal");
    proj.parameter[0] = 0.0;
    proj.parameter[1] = 0.0;
    proj.parameter[2] = 0.0;
    proj.parameter[3] = 0.0;
    proj.parameter[4] = -99.99;
    proj.parameter[5] = -99.99;
    if (Options::cdoVerbose)
      cdo_print("proj1 = %g, proj2 = %g, proj3 = %g, proj4 = %g,", proj.parameter[0], proj.parameter[1], proj.parameter[2],
                proj.parameter[3]);
  }

  if (nx == xsize && ny == ysize && cdo_cmpstr(proj.name, "sinusoidal") && cdo_cmpstr(proj.ellipsoid, "WGS-84"))
  {
    gridID = defSinusoidalGrid(nx, ny, region.xmin, region.ymax, region.dx, region.dy, proj.parameter[0], proj.parameter[1],
                               proj.parameter[2], proj.parameter[3]);
  }
  /* modification by Frank Kaspar */
  else if (nx == xsize && ny == ysize && cdo_cmpstr(proj.name, "Lambert Azimuthal Equal Area")
           && std::memcmp(proj.ellipsoid, "Sphere", 6) == 0)
  {
    double a = (proj.parameter[4] < 0) ? 6370997.0 : proj.parameter[4];
    gridID = defLaeaGrid(nx, ny, region.xmin, region.ymax, region.dx, region.dy, a, proj.parameter[2], proj.parameter[3]);
  }
  else if (std::memcmp(proj.name, "Cylindrical Equal Area", 22) == 0 && std::memcmp(proj.ellipsoid, "Sphere", 6) == 0)
  {
    double c0 = 0.001 * std::sqrt(proj.parameter[5]); /* nominal spatial resolution */
    double lts = proj.parameter[3];
    double re = proj.parameter[4] / 1000; /* Earth radius [km]*/
    if (Options::cdoVerbose) cdo_print("  c0 = %g, lts = %g, re = %g", c0, lts, re);
    gridID = defLonLatGrid(nx, ny, c0, lts, re);
  }
  else if (nx == 386 && ny == 162)
  {
    double c0 = 90; /* nominal spatial resolution */
    double lts = 30;
    double re = 6371.228; /* Earth radius [km]*/
    if (Options::cdoVerbose) cdo_print("  c0 = %g, lts = %g, re = %g", c0, lts, re);
    gridID = defLonLatGrid(nx, ny, c0, lts, re);
  }

  return gridID;
}

static int
read_region(hid_t loc_id, int nx, int ny)
{
  int gridID = -1;
  struct region_t
  {
    double area_extent[4];
    int xsize;
    int ysize;
    float xscale;
    float yscale;
    float lat_0;
    float lon_0;
    float lat_ts;
    char id[128];
    char name[128];
    char pcs_id[128];
    char pcs_def[128];
  };
  region_t region;
  char proj[128];
  double a, lon0, lat0;

  if (Options::cdoVerbose) cdo_print("Read region:");

  /*
   * Create a data type for region
   */
  hid_t region_tid = H5Tcreate(H5T_COMPOUND, sizeof(region_t));
  hsize_t dims = 4;
  hid_t fltarr_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &dims, nullptr);
  hid_t str64_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str64_tid, 128);
  hid_t str128_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str128_tid, 128);

  H5Tinsert(region_tid, "area_extent", HOFFSET(region_t, area_extent), fltarr_tid);
  H5Tinsert(region_tid, "xsize", HOFFSET(region_t, xsize), H5T_NATIVE_INT);
  H5Tinsert(region_tid, "ysize", HOFFSET(region_t, ysize), H5T_NATIVE_INT);
  H5Tinsert(region_tid, "xscale", HOFFSET(region_t, xscale), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "yscale", HOFFSET(region_t, yscale), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lat_0", HOFFSET(region_t, lat_0), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lon_0", HOFFSET(region_t, lon_0), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lat_ts", HOFFSET(region_t, lat_ts), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "id", HOFFSET(region_t, id), str64_tid);
  H5Tinsert(region_tid, "name", HOFFSET(region_t, name), str64_tid);
  H5Tinsert(region_tid, "pcs_id", HOFFSET(region_t, pcs_id), str64_tid);
  H5Tinsert(region_tid, "pcs_def", HOFFSET(region_t, pcs_def), str128_tid);

  hid_t grp_id = H5Gopen(loc_id, "/");

  hid_t region_id = H5Dopen(grp_id, "region");
  /*
  {
    hid_t tid;
    int nmem;
    int im;

    tid = H5Dget_type(proj_id);
    nmem = H5Tget_nmembers(tid);
    for ( im = 0; im < nmem; ++im )
      {
        printf("%d %s\n", im, H5Tget_member_name(tid, im));
      }
  }
  */
  herr_t status = H5Dread(region_id, region_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &region);
  (void) (status);  // unused

  if (Options::cdoVerbose)
  {
    printf("area_extent[0] = %g\n", region.area_extent[0]);
    printf("area_extent[1] = %g\n", region.area_extent[1]);
    printf("area_extent[2] = %g\n", region.area_extent[2]);
    printf("area_extent[3] = %g\n", region.area_extent[3]);
    printf("xsize = %d\n", region.xsize);
    printf("ysize = %d\n", region.ysize);
    printf("xscale = %g\n", region.xscale);
    printf("yscale = %g\n", region.yscale);
    printf("lat_0 = %g\n", region.lat_0);
    printf("lon_0 = %g\n", region.lon_0);
    printf("lat_ts = %g\n", region.lat_ts);
    printf("id = %s\n", region.id);
    printf("name = %s\n", region.name);
    printf("pcs_id = %s\n", region.pcs_id);
    printf("pcs_def = %s\n", region.pcs_def);
  }

  H5Dclose(region_id);
  H5Tclose(region_tid);
  H5Tclose(str64_tid);
  H5Tclose(str128_tid);
  H5Tclose(fltarr_tid);

  H5Gclose(grp_id);

  /* check region */

  int nfound = scan_pcs_def(region.pcs_def, proj, &a, &lon0, &lat0);

  if (Options::cdoVerbose)
  {
    printf("proj = %s\n", proj);
    printf("a    = %g\n", a);
    printf("lon0 = %g\n", lon0);
    printf("lat0 = %g\n", lat0);
  }

  double xmin = region.area_extent[0];
  double ymin = region.area_extent[1];
  double xmax = region.area_extent[2];
  double ymax = region.area_extent[3];

  /*
  xsize = (int) std::lround((region.xmax-region.xmin)/region.dx);
  ysize = (int) std::lround((region.ymax-region.ymin)/region.dy);

  if ( Options::cdoVerbose ) cdo_print("  Size: xsize=%d  ysize=%d", xsize, ysize);
  */

  if (nfound == 4 && nx == region.xsize && ny == region.ysize && cdo_cmpstr(proj, "laea"))
  {
    double dx = (xmax - xmin) / nx;
    double dy = (ymax - ymin) / ny;
    gridID = defLaeaGrid(nx, ny, xmin, ymax, dx, dy, a, lon0, lat0);
  }

  return gridID;
}

static void
read_dataset(hid_t loc_id, const char *name, void *opdata)
{
  hid_t dataspace;
  hsize_t dims_out[9];  // dataset dimensions
  herr_t status;        // Generic return value
  hid_t attr, atype, atype_mem;
  int iattr;
  float fattr;
  double dattr;
  char attname[CDI_MAX_NAME];
  H5T_class_t atype_class;
  size_t atype_size;
  int rank;
  int nx = 0, ny = 0, nz = 0, nt = 0;
  size_t gridsize;
  double addoffset = 0.0, scalefactor = 1.0, missval = cdiInqMissval();
  bool haveAddoffset = false, haveScalefactor = false, haveMissvals = false;
  int nset;
  int ftype = 0;
  int len;
  int dtype = CDI_DATATYPE_FLT32;
  char attstring[4096];  // Buffer to read string attribute back
  char varname[CDI_MAX_NAME];
  size_t numMissVals;
  int num_attrs;

  attstring[0] = 0;
  std::strcpy(varname, name);

  hid_t dset_id = H5Dopen(loc_id, varname);

  hid_t type_id = H5Dget_type(dset_id); /* get datatype*/

  H5T_class_t type_class = H5Tget_class(type_id);
  if (type_class < 0) { cdo_abort(" Invalid datatype for %s", varname); }
  /*
  else {
    if(type_class == H5T_INTEGER)  puts("   Datatype is 'H5T_NATIVE_INTEGER'.\n");
    if(type_class == H5T_FLOAT)    puts("   Datatype is 'H5T_NATIVE_FLOAT'.\n");
    if(type_class == H5T_STRING)   puts("   Datatype is 'H5T_NATIVE_STRING'.\n");
    if(type_class == H5T_BITFIELD) puts("   Datatype is 'H5T_NATIVE_BITFIELD'.\n");
    if(type_class == H5T_OPAQUE)   puts("   Datatype is 'H5T_NATIVE_OPAQUE'.\n");
    if(type_class == H5T_COMPOUND) puts("   Datatype is 'H5T_NATIVE_COMPOUND'.\n");
  }
  */
  hid_t native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  if (H5Tequal(native_type, H5T_NATIVE_SCHAR) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_INT8;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_UCHAR) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_UINT8;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_SHORT) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_INT16;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_USHORT) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_UINT16;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_INT) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_INT32;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_UINT) > 0)
  {
    ftype = 0;
    dtype = CDI_DATATYPE_UINT32;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_FLOAT) > 0)
  {
    ftype = 1;
    dtype = CDI_DATATYPE_FLT32;
  }
  else if (H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0)
  {
    ftype = 1;
    dtype = CDI_DATATYPE_FLT64;
  }
  else
  {
    cdo_warning("Dataset %s skipped, unsupported native datatype!", varname);
    goto RETURN;
  }
  H5Tclose(native_type);

  dataspace = H5Dget_space(dset_id);  // dataspace handle
  rank = H5Sget_simple_extent_ndims(dataspace);
  status = H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);

  if (rank == 2)
  {
    nx = dims_out[1];
    ny = dims_out[0];
    nz = 1;
    nt = 1;
  }
  else if (rank == 3)
  {
    nx = dims_out[2];
    ny = dims_out[1];
    nz = 1;
    nt = dims_out[0];
  }
  else
  {
    cdo_warning("Dataset %s skipped, unsupported rank (=%d)!", varname, rank);
    goto RETURN;
  }

  len = (int) std::strlen(varname);
  if (len > 0 && ((datasets_t *) opdata)->mergelevel)
    if (std::isdigit(varname[len - 1]) && std::memcmp(varname, "Data", 4) != 0)
    {
      if (nt > 1) cdo_abort("Combination of nlevel > 1 and ntime > 1 not implemented!");

      nz = atoi(&varname[len - 1]);
      varname[len - 1] = 0;
    }

  gridsize = nx * ny;

  if (nz == 1)
    nset = ((datasets_t *) opdata)->numSets;
  else
  {
    for (nset = 0; nset < ((datasets_t *) opdata)->numSets; ++nset)
    {
      if (cdo_cmpstr(varname, ((datasets_t *) opdata)->obj[nset].name)) break;
    }

    if (nset >= ((datasets_t *) opdata)->numSets) cdo_abort("3D var %s not found!", varname);
  }

  if (nset < MAX_DSETS)
  {
    if (Options::cdoVerbose) print_filter(dset_id, varname);

    num_attrs = H5Aget_num_attrs(dset_id);
    for (int i = 0; i < num_attrs; ++i)
    {
      attr = H5Aopen_idx(dset_id, i);
      atype = H5Aget_type(attr);
      H5Aget_name(attr, sizeof(attname), attname);

      if (cdo_cmpstr(attname, "CLASS") || cdo_cmpstr(attname, "IMAGE_VERSION") || cdo_cmpstr(attname, "PALETTE")) continue;

      atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
      atype_size = H5Tget_size(atype);
      atype_class = H5Tget_class(atype);

      len = std::strlen(attname);
      for (int k = 0; k < len; ++k) attname[k] = std::tolower(attname[k]);

      if (cdo_cmpstr(attname, "intercept") || cdo_cmpstr(attname, "offset"))
      {
        if (atype_class == H5T_FLOAT)
        {
          if (atype_size == 4)
          {
            status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
            if (status >= 0)
            {
              addoffset = fattr;
              haveAddoffset = true;
            }
          }
          else
          {
            status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
            if (status >= 0)
            {
              addoffset = dattr;
              haveAddoffset = true;
            }
          }

          if (!haveAddoffset) cdo_warning("Reading of float attribute %s failed!", attname);
        }
        else if (atype_class == H5T_INTEGER)
        {
          status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
          if (status >= 0)
          {
            addoffset = iattr;
            haveAddoffset = true;
          }
          else
            cdo_warning("Reading of integer attribute %s failed!", attname);
        }
        else
          cdo_warning("Attribute %s has unsupported data type!", attname);
      }
      else if (cdo_cmpstr(attname, "gain") || cdo_cmpstr(attname, "scaling_factor"))
      {
        if (atype_class == H5T_FLOAT)
        {
          if (atype_size == 4)
          {
            status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
            if (status >= 0)
            {
              scalefactor = fattr;
              haveScalefactor = true;
            }
          }
          else
          {
            status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
            if (status >= 0)
            {
              scalefactor = dattr;
              haveScalefactor = true;
            }
          }

          if (!haveScalefactor) cdo_warning("Reading of float attribute %s failed!", attname);
        }
        else if (atype_class == H5T_INTEGER)
        {
          status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
          if (status >= 0)
          {
            scalefactor = iattr;
            haveScalefactor = true;
          }
          else
            cdo_warning("Reading of integer attribute %s failed!", attname);
        }
        else
          cdo_warning("Attribute %s has unsupported data type!", attname);
      }
      else if (std::strncmp(attname, "no_data", 7) == 0 || std::strncmp(attname, "nodata", 6) == 0)
      {
        if (atype_class == H5T_FLOAT)
        {
          if (atype_size == 4)
          {
            status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
            if (status >= 0)
            {
              missval = fattr;
              haveMissvals = true;
            }
          }
          else
          {
            status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
            if (status >= 0)
            {
              missval = dattr;
              haveMissvals = true;
            }
          }

          if (!haveMissvals) cdo_warning("Reading of float attribute %s failed!", attname);
        }
        else if (atype_class == H5T_INTEGER)
        {
          status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
          if (status >= 0)
          {
            missval = iattr;
            haveMissvals = true;
          }
          else
            cdo_warning("Reading of integer attribute %s failed!", attname);
        }
        else
          cdo_warning("Attribute %s has unsupported data type!", attname);
      }
      else if (cdo_cmpstr(attname, "description"))
      {
        H5Aread(attr, atype_mem, attstring);
        if (((datasets_t *) opdata)->obj[nset].description) std::free(((datasets_t *) opdata)->obj[nset].description);
        ((datasets_t *) opdata)->obj[nset].description = strdup(attstring);
      }
      else if (cdo_cmpstr(attname, "title"))
      {
        H5Aread(attr, atype_mem, attstring);
        if (((datasets_t *) opdata)->obj[nset].title) std::free(((datasets_t *) opdata)->obj[nset].title);
        ((datasets_t *) opdata)->obj[nset].title = strdup(attstring);
      }
      else if (cdo_cmpstr(attname, "time"))
      {
        H5Aread(attr, atype_mem, attstring);
        if (((datasets_t *) opdata)->obj[nset].time) std::free(((datasets_t *) opdata)->obj[nset].time);
        ((datasets_t *) opdata)->obj[nset].time = strdup(attstring);
      }
      else if (cdo_cmpstr(attname, "unit"))
      {
        H5Aread(attr, atype_mem, attstring);
        ((datasets_t *) opdata)->obj[nset].units = strdup(attstring);
      }

      H5Tclose(atype_mem);
      H5Aclose(attr);
      H5Tclose(atype);
    }

    int offset = gridsize * (nz - 1);
    double *array = ((datasets_t *) opdata)->obj[nset].array;
    array = (double *) std::realloc(array, gridsize * nz * nt * sizeof(double));
    ((datasets_t *) opdata)->obj[nset].array = array;
    array = array + offset;

    if (ftype)
    {
      if (dtype == CDI_DATATYPE_FLT32)
      {
        std::vector<float> farray(gridsize * nt);
        status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, farray.data());
        if (status < 0) cdo_abort("Reading of NATIVE_FLOAT variable %s failed!", varname);
        for (size_t i = 0; i < gridsize * nt; ++i) array[i] = farray[i];
      }
      else
      {
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
        if (status < 0) cdo_abort("Reading of NATIVE_DOUBLE variable %s failed!", varname);
      }
    }
    else
    {
      Varray<int> iarray(gridsize * nt);
      status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray.data());
      if (status < 0) cdo_abort("Reading of NATIVE_INT variable %s failed!", varname);
      for (size_t i = 0; i < gridsize * nt; ++i) array[i] = iarray[i];
    }

    ((datasets_t *) opdata)->obj[nset].name = strdup(varname);
    ((datasets_t *) opdata)->obj[nset].nx = nx;
    ((datasets_t *) opdata)->obj[nset].ny = ny;
    ((datasets_t *) opdata)->obj[nset].nz = nz;
    ((datasets_t *) opdata)->obj[nset].nt = nt;
    ((datasets_t *) opdata)->obj[nset].gridsize = gridsize;

    if (nz > 1)
    {
      if (((datasets_t *) opdata)->obj[nset].dtype != dtype) cdo_warning("Data type changes over levels!");

      if (haveAddoffset && fp_is_not_equal(((datasets_t *) opdata)->obj[nset].addoffset, addoffset))
        cdo_warning("Offset changes over levels!");

      if (haveScalefactor && fp_is_not_equal(((datasets_t *) opdata)->obj[nset].scalefactor, scalefactor))
        cdo_warning("Scalefactor changes over levels!");

      if (haveMissvals && fp_is_not_equal(((datasets_t *) opdata)->obj[nset].missval, missval))
        cdo_warning("Missing value changes over levels!");
    }

    if (nz == 1) ((datasets_t *) opdata)->numSets++;

    std::vector<bool> mask(gridsize * nt, false);

    numMissVals = 0;

    auto mm = varray_min_max(gridsize * nt, array);

    if (Options::cdoVerbose)
      cdo_print("Dataset %s: missval = %g  addoffset = %g  scalefactor = %g", varname, missval, addoffset, scalefactor);

    if (Options::cdoVerbose)
      cdo_print("Dataset %s: dtype = %d  minval = %g  maxval = %g  missval = %g", varname, dtype, mm.min, mm.max, missval);

    if (dtype == CDI_DATATYPE_UINT8)
    {
      if (mm.min >= 0 && mm.max <= 127) dtype = CDI_DATATYPE_INT8;
    }
    else if (dtype == CDI_DATATYPE_UINT16)
    {
      if (mm.min >= 0 && mm.max <= 32767) dtype = CDI_DATATYPE_INT16;
    }

    haveAddoffset = is_not_equal(addoffset, 0.0);
    haveScalefactor = is_not_equal(scalefactor, 1.0);

    if (haveAddoffset || haveScalefactor)
    {
      for (size_t i = 0; i < gridsize * nt; ++i)
        if (fp_is_not_equal(array[i], missval))
        {
          mask[i] = false;

          if (haveScalefactor) array[i] *= scalefactor;
          if (haveAddoffset) array[i] += addoffset;
        }
        else
        {
          numMissVals++;
          mask[i] = true;
        }
    }

    double minval = 1e35;
    double maxval = -1e35;
    for (size_t i = 0; i < gridsize * nt; ++i)
      if (mask[i] == false)
      {
        minval = std::min(minval, array[i]);
        maxval = std::max(maxval, array[i]);
      }

    if (Options::cdoVerbose)
      cdo_print("Dataset %s: dtype = %d  minval = %g  maxval = %g  missval = %g", varname, dtype, minval, maxval, missval);

    if (numMissVals)
    {
      if (!(missval < minval || missval > maxval))
      {
        if (fp_is_equal(missval, 255.) && dtype == CDI_DATATYPE_UINT8)
        {
          missval = -255;
          dtype = CDI_DATATYPE_INT16;
          cdo_print("Dataset %s: changed missval to %g and datatype to INT16!", varname, missval);

          for (size_t i = 0; i < gridsize * nt; ++i)
            if (mask[i]) array[i] = missval;
        }
        else
          cdo_warning(" Missing value is inside the range of valid values!\n"
                      "\tDataset %s,  Missval: %g,  Range: %g - %g",
                      varname, missval, minval, maxval);
      }
    }

    ((datasets_t *) opdata)->obj[nset].dtype = dtype;
    ((datasets_t *) opdata)->obj[nset].haveAddoffset = haveAddoffset;
    ((datasets_t *) opdata)->obj[nset].haveScalefactor = haveScalefactor;
    ((datasets_t *) opdata)->obj[nset].haveMissvals = haveMissvals;
    ((datasets_t *) opdata)->obj[nset].addoffset = addoffset;
    ((datasets_t *) opdata)->obj[nset].scalefactor = scalefactor;
    ((datasets_t *) opdata)->obj[nset].missval = missval;
  }
  else
  {
    cdo_warning("Too many datasets (MAX = %d)!", MAX_DSETS);
    goto RETURN;
  }

  H5Sclose(dataspace);

RETURN:

  H5Dclose(dset_id);
  H5Tclose(type_id);
}

static herr_t
obj_info(hid_t loc_id, const char *name, void *opdata)
{
  H5G_stat_t statbuf;
  H5Gget_objinfo(loc_id, name, false, &statbuf);

  H5G_obj_t obj_type = statbuf.type;

  switch (obj_type)
  {
    case H5G_GROUP:
      if (Options::cdoVerbose) cdo_print(" Object with name %s is a group", name);
      if (cdo_cmpstr(name, "Data"))
      {
        ((datasets_t *) opdata)->mergelevel = true;
        H5Giterate(loc_id, name, nullptr, obj_info, opdata);
      }
      else if (cdo_cmpstr(name, "Geolocation")) { ((datasets_t *) opdata)->lgeoloc = true; }
      else if (cdo_cmpstr(name, "Metadata")) { ((datasets_t *) opdata)->lmetadata = true; }
      break;
    case H5G_DATASET:
      if (Options::cdoVerbose) cdo_print(" Object with name %s is a dataset", name);
      if (std::strstr(name, "PALETTE"))
      {
        if (Options::cdoVerbose) cdo_print("   Skip dataset: %s", name);
      }
      // else if (std::strstr(name, "egion"))
      else if (cdo_cmpstr(name, "region")) { ((datasets_t *) opdata)->lregion = true; }
      else
      {
        if (Options::cdoVerbose) cdo_print("   Read dataset: %s", name);
        read_dataset(loc_id, name, opdata);
      }
      break;
    case H5G_TYPE:
      if (Options::cdoVerbose) cdo_print(" Object with name %s is a named datatype", name);
      if (cdo_cmpstr(name, "ProjType")) { ((datasets_t *) opdata)->lprojtype = true; }
      break;
    default: cdo_abort(" Unable to identify an object %s", name); break;
  }

  return 0;
}

static void
get_global_att(hid_t file_id, const char *obj_path, int vlistID)
{
  hid_t attr, atype, atype_mem, obj_id, grp_id = -1;
  char attname[CDI_MAX_NAME];
  H5T_class_t type_class;
  int attint;
  double attflt;
  int pos;
  int num_attrs;
  char attstring[4096]; /* Buffer to read string attribute back */

  attstring[0] = 0;

  obj_id = H5Gopen(file_id, obj_path);

  num_attrs = H5Aget_num_attrs(obj_id);

  for (int i = 0; i < num_attrs; ++i)
  {
    attr = H5Aopen_idx(obj_id, i);
    atype = H5Aget_type(attr);
    H5Aget_name(attr, sizeof(attname), attname);

    /* remove illegal characters */
    for (pos = 0; pos < (int) std::strlen(attname); ++pos)
      if (attname[pos] == '&') attname[pos] = '_';

    atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
    type_class = H5Tget_class(atype);
    if (type_class == H5T_STRING)
    {
      H5Aread(attr, atype_mem, attstring);
      cdiDefAttTxt(vlistID, CDI_GLOBAL, attname, (int) std::strlen(attstring), attstring);
    }
    else if (type_class == H5T_INTEGER)
    {
      H5Aread(attr, H5T_NATIVE_INT, &attint);
      cdiDefAttInt(vlistID, CDI_GLOBAL, attname, CDI_DATATYPE_INT32, 1, &attint);
    }
    else if (type_class == H5T_FLOAT)
    {
      H5Aread(attr, H5T_NATIVE_DOUBLE, &attflt);
      cdiDefAttFlt(vlistID, CDI_GLOBAL, attname, CDI_DATATYPE_FLT64, 1, &attflt);
    }
    H5Tclose(atype_mem);
    H5Aclose(attr);
    H5Tclose(atype);
  }

  if (grp_id >= 0) H5Gclose(grp_id);
}

static int
get_vdate(int vlistID)
{
  int64_t vdate = 0;
  int natts;
  int len, type;
  char name[CDI_MAX_NAME];
  char attstr[CDI_MAX_NAME];

  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);

  for (int i = 0; i < natts; ++i)
  {
    cdiInqAtt(vlistID, CDI_GLOBAL, i, name, &type, &len);
    if (type == CDI_DATATYPE_TXT)
    {
      if (cdo_cmpstr(name, "DateAndTime") || cdo_cmpstr(name, "Date_Time"))
      {
        cdiInqAttTxt(vlistID, CDI_GLOBAL, name, CDI_MAX_NAME, attstr);
        if (len > 8) len = 8;
        attstr[len] = 0;
        vdate = atoi(attstr);
        if (vdate < 999999) vdate = vdate * 100 + 1;
      }
    }
  }

  return vdate;
}

static void
dsets_init(datasets_t *dsets)
{
  dsets->numSets = 0;
  dsets->mergelevel = 0;
  dsets->lgeoloc = 0;
  dsets->lregion = 0;
  dsets->lprojtype = 0;
  dsets->lmetadata = 0;

  for (int i = 0; i < MAX_DSETS; ++i)
  {
    dsets->obj[i].nx = 0;
    dsets->obj[i].ny = 0;
    dsets->obj[i].nz = 0;
    dsets->obj[i].name = nullptr;
    dsets->obj[i].description = nullptr;
    dsets->obj[i].units = nullptr;
    dsets->obj[i].title = nullptr;
    dsets->obj[i].time = nullptr;
    dsets->obj[i].dtype = CdoDefault::DataType;
    dsets->obj[i].haveScalefactor = false;
    dsets->obj[i].haveAddoffset = false;
    dsets->obj[i].haveMissvals = false;
    dsets->obj[i].missval = cdiInqMissval();
    dsets->obj[i].array = nullptr;
  }
}
#endif
class Importcmsaf : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Importcmsaf",
    .operators = { { "import_cmsaf", ImportcmsafHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Importcmsaf> registration = RegisterEntry<Importcmsaf>();
#ifdef HAVE_LIBHDF5
  int gridID = -1, zaxisID;
  size_t numMissVals;
  int ivar;
  int varID, levelID, tsID;
  hid_t file_id;
  datasets_t dsets;
#endif

  CdoStreamID streamID;
  int taxisID;
  int vlistID;

public:
  void
  init() override
  {

    operator_check_argc(0);

    if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = CDI_FILETYPE_NC;

#ifdef HAVE_LIBHDF5
    dsets_init(&dsets);

    // Open an existing file
    file_id = H5Fopen(cdo_get_stream_name(0), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) cdo_abort("H5Fopen failed on %s", cdo_get_stream_name(0));

    // cmsaf_type = get_cmsaf_type(file_id);

    H5Giterate(file_id, "/", nullptr, obj_info, (void *) &dsets);

    if (dsets.numSets == 0) cdo_abort("No dataset found!");
#endif
  }

  void
  run() override
  {

#ifdef HAVE_LIBHDF5
    int gridsize = dsets.obj[0].gridsize;
    int nx = dsets.obj[0].nx;
    int ny = dsets.obj[0].ny;
    int nz = dsets.obj[0].nz;
    int nt = dsets.obj[0].nt;

    for (ivar = 0; ivar < dsets.numSets; ++ivar)
      if (dsets.obj[ivar].nt > 1)
      {
        nt = dsets.obj[ivar].nt;
        break;
      }

    Varray<int> vtimes;
    if (nt > 1)
    {
      vtimes.resize(nt);

      for (int i = 0; i < nt; ++i) vtimes[i] = i * 10000 + 45 * 100;

      if (dsets.obj[ivar].time)
      {
        char *pline = dsets.obj[ivar].time;
        for (int i = 0; i < nt; ++i)
        {
          int itime = ((int) std::strtol(pline, &pline, 10)) * 100;
          if (itime < 0 || itime > 240000)
          {
            cdo_warning("Wrong time string!");
            break;
          }
          vtimes[i] = itime;
        }
      }
    }

    if (Options::cdoVerbose)
      for (ivar = 0; ivar < dsets.numSets; ++ivar)
        cdo_print(" Var %d %-20s %dx%d nlev = %d nts = %d", ivar, dsets.obj[ivar].name, nx, ny, nz, dsets.obj[ivar].nt);

    for (ivar = 1; ivar < dsets.numSets; ++ivar)
    {
      if (nx != dsets.obj[0].nx || ny != dsets.obj[0].ny) cdo_abort("Gridsize must not change!");
      if (nz != dsets.obj[0].nz) cdo_abort("Number of levels must not change!");
    }

    if (dsets.lgeoloc) { gridID = read_geolocation(file_id, nx, ny, dsets.lprojtype); }
    else if (dsets.lregion) { gridID = read_region(file_id, nx, ny); }

    if (gridID == -1)
    {
      gridID = gridCreate(GRID_GENERIC, gridsize);
      gridDefXsize(gridID, nx);
      gridDefYsize(gridID, ny);
    }

    if (nz == 1)
      zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
    else
    {
      Varray<double> levels(nz);
      for (int i = 0; i < nz; ++i) levels[i] = i + 1;
      zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
      zaxisDefLevels(zaxisID, levels.data());
    }

    vlistID = vlistCreate();

    taxisID = cdo_taxis_create((nt > 1) ? TAXIS_RELATIVE : TAXIS_ABSOLUTE);
    taxisDefCalendar(taxisID, CALENDAR_STANDARD);
    vlistDefTaxis(vlistID, taxisID);

    for (ivar = 0; ivar < dsets.numSets; ++ivar)
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, dsets.obj[ivar].name);
      if (dsets.obj[ivar].description) cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, dsets.obj[ivar].description);
      if (dsets.obj[ivar].units) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, dsets.obj[ivar].units);
      if (dsets.obj[ivar].title)
        cdiDefAttTxt(vlistID, varID, "title", (int) std::strlen(dsets.obj[ivar].title), dsets.obj[ivar].title);

      vlistDefVarDatatype(vlistID, varID, dsets.obj[ivar].dtype);
      if (dsets.obj[ivar].haveMissvals) vlistDefVarMissval(vlistID, varID, dsets.obj[ivar].missval);
      if (dsets.obj[ivar].haveScalefactor) cdiDefKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, dsets.obj[ivar].scalefactor);
      if (dsets.obj[ivar].haveAddoffset) cdiDefKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, dsets.obj[ivar].addoffset);
    }

    get_global_att(file_id, "/", vlistID);
    if (dsets.lmetadata) get_global_att(file_id, "Metadata", vlistID);

    int vdate = get_vdate(vlistID);
    if (vdate == 0) vdate = 10101;

    streamID = cdo_open_write(1);

    cdo_def_vlist(streamID, vlistID);

    for (tsID = 0; tsID < nt; ++tsID)
    {
      int vtime = (vtimes.empty()) ? 0 : vtimes[tsID];
      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(vdate);
      vDateTime.time = cdiTime_set(vtime);
      taxisDefVdatetime(taxisID, vDateTime);
      cdo_def_timestep(streamID, tsID);

      for (ivar = 0; ivar < dsets.numSets; ++ivar)
      {
        varID = ivar;

        if (tsID > 0 && dsets.obj[ivar].nt == 1) continue;

        gridsize = dsets.obj[ivar].gridsize;
        double missval = dsets.obj[ivar].missval;

        for (levelID = 0; levelID < nz; ++levelID)
        {
          int offset = gridsize * levelID;
          if (nz == 1) offset = gridsize * tsID;
          double *array = dsets.obj[ivar].array + offset;

          auto mm = varray_min_max_mv(gridsize, array, missval);
          numMissVals = gridsize - mm.n;

          if (Options::cdoVerbose)
            cdo_print(" Write var %d,  level %d, numMissVals %zu, missval %g, minval %g, maxval %g", varID, levelID, numMissVals,
                      missval, mm.min, mm.max);

          cdo_def_field(streamID, varID, levelID);
          cdo_write_field(streamID, array, numMissVals);
        }
      }
    }
#else
    cdo_abort("HDF5 support not compiled in!");
#endif
  }

  void
  close() override
  {
#ifdef HAVE_LIBHDF5
    // Close file
    (void) H5Fclose(file_id);

    process_def_var_num(vlistNvars(vlistID));

    cdo_stream_close(streamID);
    vlistDestroy(vlistID);

    for (ivar = 0; ivar < dsets.numSets; ++ivar) std::free(dsets.obj[ivar].array);

#endif
  }
};
