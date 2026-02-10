/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#include "cdo_options.h"
#include "process_int.h"
#include "griddes.h"
#include <mpim_grid.h>
#include "remap_vars.h"
#include "knndata.h"

#ifdef HAVE_LIBNETCDF
static void
nce(int istat)
{
  // This routine provides a simple interface to NetCDF error message routine.
  if (istat != NC_NOERR) cdo_abort(nc_strerror(istat));
}

static size_t
cdf_read_dimlen(int ncfileid, const char *dimname)
{
  size_t dimlen = 0;
  int ncdimid;
  auto status = nc_inq_dimid(ncfileid, dimname, &ncdimid);
  if (status == NC_NOERR) nce(nc_inq_dimlen(ncfileid, ncdimid, &dimlen));
  return dimlen;
}

static void
cdf_read_att_text(int ncfileid, int ncvarid, const char *name, std::string &text)
{
  size_t attlen;
  nc_type atttype;
  nce(nc_inq_atttype(ncfileid, ncvarid, name, &atttype));
  nce(nc_inq_attlen(ncfileid, ncvarid, name, &attlen));

  if (atttype == NC_CHAR)
  {
    char *attbuf = new char[attlen + 1];

    nce(nc_get_att_text(ncfileid, ncvarid, name, attbuf));
    attbuf[attlen] = 0;
    text = attbuf;

    delete[] attbuf;
  }
}

static void
cdf_read_var_int(int ncfileid, const char *name, int *array)
{
  int ncvarid;
  nce(nc_inq_varid(ncfileid, name, &ncvarid));
  nce(nc_get_var_int(ncfileid, ncvarid, array));
}

static void
cdf_read_var_size(int ncfileid, const char *name, size_t len, size_t *array)
{
  if (len < 0x7FFFFC00)  // 2GB
  {
    std::vector<int> iarray(len);
    cdf_read_var_int(ncfileid, name, iarray.data());
    for (size_t i = 0; i < len; ++i) array[i] = (size_t) iarray[i];
  }
  else
  {
    int ncvarid;
    nce(nc_inq_varid(ncfileid, name, &ncvarid));
    nce(nc_get_var_ulonglong(ncfileid, ncvarid, (unsigned long long *) array));
  }
}

static void
cdf_read_var_double(int ncfileid, const char *name, double *array)
{
  int ncvarid;
  nce(nc_inq_varid(ncfileid, name, &ncvarid));
  nce(nc_get_var_double(ncfileid, ncvarid, array));
}

static void
cdf_read_coordinate_radian(int ncfileid, const char *name, Varray<double> &varray)
{
  int ncvarid;
  nce(nc_inq_varid(ncfileid, name, &ncvarid));
  nce(nc_get_var_double(ncfileid, ncvarid, varray.data()));

  std::string grid_units;
  cdf_read_att_text(ncfileid, ncvarid, "units", grid_units);
  if (string_to_LonLatUnits(grid_units, name) == LonLatUnits::Deg) grid_to_radian(varray);
}

namespace
{
struct RemapAttributes
{
  std::string map_name;
  std::string history;
  std::string cdo_version;
};
}  // namespace

struct RemapGridW
{
  std::string name;
  int rank = 0;               // rank of the grid
  size_t ncells = 0;          // total points on the grid
  size_t ncorners = 0;        // number of corners for each grid cell
  size_t dims[2] = { 0, 0 };  // size of grid dimension

  std::vector<int> mask;  // flag which cells participate

  Varray<double> cell_center_lon;  // lon/lat coordinates for
  Varray<double> cell_center_lat;  // each grid center in radians
  Varray<double> cell_corner_lon;  // lon/lat coordinates for
  Varray<double> cell_corner_lat;  // each grid corner in radians

  Varray<double> cell_area;  // total area of each grid cell
  Varray<double> cell_frac;  // fractional area of grid cells participating in remapping

  void
  set_num_cells(size_t ncells_)
  {
    ncells = ncells_;

    mask.resize(ncells);

    cell_center_lon.resize(ncells);
    cell_center_lat.resize(ncells);

    cell_frac.resize(ncells, 0.0);
  }

  void
  set_num_corners(size_t ncorners_)
  {
    ncorners = ncorners_;

    cell_corner_lon.resize(ncorners * ncells, 0.0);
    cell_corner_lat.resize(ncorners * ncells, 0.0);
  }
};

static void
remapGridWAlloc(RemapMethod mapType, RemapGridW &grid)
{
  if (mapType == RemapMethod::CONSERV) { grid.cell_area.resize(grid.ncells, 0.0); }
}

struct RemapVarsW
{
  RemapMethod mapType{ RemapMethod::UNDEF };  // identifier for remapping method
  NormOpt normOpt{ NormOpt::NONE };           // option for normalization (conserv only)
  size_t numLinks = 0;                        // number of links for remapping
  size_t numWeights = 0;                      // number of weights used in remapping

  Varray<size_t> srcCellIndices;  // source grid indices for each link
  Varray<size_t> tgtCellIndices;  // target grid indices for each link
  Varray<double> weights;         // map weights for each link [numLinks*numWeights]
};

static void
remapVarsWInit(RemapMethod mapType, int remapOrder, RemapVarsW &rv)
{
  // Determine the number of weights
  rv.numWeights = (mapType == RemapMethod::BICUBIC) ? 4 : 1;
  if (mapType == RemapMethod::CONSERV && remapOrder == 2) rv.numWeights = 3;
}

static NormOpt
remapGetNormOpt(std::string const &normOptStr)
{
  NormOpt normOpt(NormOpt::NONE);
  // clang-format off
  if      (normOptStr == "none")     normOpt = NormOpt::NONE;
  else if (normOptStr == "fracarea") normOpt = NormOpt::FRACAREA;
  else if (normOptStr == "destarea") normOpt = NormOpt::DESTAREA;
  else
    {
      cdo_print("normalize_opt = %s", normOptStr);
      cdo_abort("Invalid normalization option");
    }
  // clang-format on

  if (Options::cdoVerbose) cdo_print("normalize_opt = %s", normOptStr);

  return normOpt;
}

RemapSwitches get_maptype(int ncfileid);
KnnParams read_knn_params(int ncfileid);

static void
read_remapgrid_scrip(int ncfileid, std::string const &prefix, bool lgridarea, RemapGridW &grid)
{
  // Read all variables of the grid
  cdf_read_var_size(ncfileid, (prefix + "_dims").c_str(), 2, grid.dims);
  cdf_read_var_int(ncfileid, (prefix + "_imask").c_str(), grid.mask.data());

  cdf_read_coordinate_radian(ncfileid, (prefix + "_center_lat").c_str(), grid.cell_center_lat);
  cdf_read_coordinate_radian(ncfileid, (prefix + "_center_lon").c_str(), grid.cell_center_lon);

  if (grid.ncorners)
  {
    cdf_read_coordinate_radian(ncfileid, (prefix + "_corner_lat").c_str(), grid.cell_corner_lat);
    cdf_read_coordinate_radian(ncfileid, (prefix + "_corner_lon").c_str(), grid.cell_corner_lon);
  }

  if (lgridarea) cdf_read_var_double(ncfileid, (prefix + "_area").c_str(), grid.cell_area.data());
  cdf_read_var_double(ncfileid, (prefix + "_frac").c_str(), grid.cell_frac.data());
}

static void
read_remapweights_scrip(int ncfileid, RemapVarsW &rv)
{
  cdf_read_var_size(ncfileid, "src_address", rv.numLinks, rv.srcCellIndices.data());
  cdf_read_var_size(ncfileid, "tgt_address", rv.numLinks, rv.tgtCellIndices.data());

  for (size_t i = 0; i < rv.numLinks; ++i) rv.srcCellIndices[i]--;
  for (size_t i = 0; i < rv.numLinks; ++i) rv.tgtCellIndices[i]--;

  cdf_read_var_double(ncfileid, "remap_matrix", rv.weights.data());
}

static RemapSwitches
readRemapFileScrip(std::string const &remapFile, RemapAttributes &remapAtts, RemapGridW &srcGrid, RemapGridW &tgtGrid,
                   RemapVarsW &rv)
{
  // The routine reads a NetCDF file to extract remapping info in SCRIP format

  // Open file and read some global information
  auto ncfileid = cdo_cdf_openread(remapFile.c_str());

  cdf_read_att_text(ncfileid, NC_GLOBAL, "history", remapAtts.history);

  // Map name
  cdf_read_att_text(ncfileid, NC_GLOBAL, "title", remapAtts.map_name);

  if (Options::cdoVerbose)
  {
    cdo_print("Reading remapping: %s", remapAtts.map_name);
    cdo_print("From file: %s", remapFile);
  }

  // Map Type
  auto remapSwitches = get_maptype(ncfileid);
  if (remapSwitches.mapType == RemapMethod::KNN && remapSwitches.numNeighbors == -1) auto knnParams = read_knn_params(ncfileid);
  auto lgridarea = (remapSwitches.mapType == RemapMethod::CONSERV);

  remapVarsWInit(remapSwitches.mapType, remapSwitches.remapOrder, rv);

  rv.mapType = remapSwitches.mapType;

  // Normalization option
  std::string normalizeOptStr;  // character string for normalization option
  cdf_read_att_text(ncfileid, NC_GLOBAL, "normalization", normalizeOptStr);
  rv.normOpt = remapGetNormOpt(normalizeOptStr);

  // File convention
  std::string convention;  // character string for output convention
  cdf_read_att_text(ncfileid, NC_GLOBAL, "conventions", convention);
  if (convention != "SCRIP")
  {
    cdo_print("convention = %s", convention);
    if (convention == "NCAR-CSM")
      cdo_abort("Unsupported file convention: %s!", convention);
    else
      cdo_abort("Unknown file convention!");
  }

  // Read some additional global attributes

  // Source and destination grid names
  cdf_read_att_text(ncfileid, NC_GLOBAL, "source_grid", srcGrid.name);
  cdf_read_att_text(ncfileid, NC_GLOBAL, "dest_grid", tgtGrid.name);

  if (Options::cdoVerbose) cdo_print("Remapping between: %s and %s", srcGrid.name, tgtGrid.name);

  // Read dimension information
  srcGrid.set_num_cells(cdf_read_dimlen(ncfileid, "src_grid_size"));
  tgtGrid.set_num_cells(cdf_read_dimlen(ncfileid, "dst_grid_size"));

  srcGrid.set_num_corners(cdf_read_dimlen(ncfileid, "src_grid_corners"));
  tgtGrid.set_num_corners(cdf_read_dimlen(ncfileid, "dst_grid_corners"));

  srcGrid.rank = cdf_read_dimlen(ncfileid, "src_grid_rank");
  tgtGrid.rank = cdf_read_dimlen(ncfileid, "dst_grid_rank");

  remapGridWAlloc(rv.mapType, srcGrid);
  remapGridWAlloc(rv.mapType, tgtGrid);

  rv.numLinks = cdf_read_dimlen(ncfileid, "numLinks");
  //  if ( rv.numLinks == 0 ) cdo_abort("Number of remap links is 0, no remap weights found!");
  rv.numWeights = cdf_read_dimlen(ncfileid, "num_wgts");

  // Allocate address and weight arrays
  if (rv.numLinks > 0)
  {
    rv.srcCellIndices.resize(rv.numLinks);
    rv.tgtCellIndices.resize(rv.numLinks);
    rv.weights.resize(rv.numWeights * rv.numLinks);
  }

  read_remapgrid_scrip(ncfileid, "src_grid", lgridarea, srcGrid);
  read_remapgrid_scrip(ncfileid, "dst_grid", lgridarea, tgtGrid);

  if (rv.numLinks > 0) read_remapweights_scrip(ncfileid, rv);

  // Close input file
  cdo_cdf_close(ncfileid);

  return remapSwitches;
}  // readRemapFileScrip

static int
cdf_def_dim(int ncfileid, const char *name, size_t len)
{
  int ncdimid = -1;
  nce(nc_def_dim(ncfileid, name, len, &ncdimid));
  return ncdimid;
}

static int
cdf_def_var(int ncfileid, const char *name, nc_type xtype, int ndims, const int *dimidsp)
{
  int ncvarid = -1;
  nce(nc_def_var(ncfileid, name, xtype, ndims, dimidsp, &ncvarid));
  return ncvarid;
}

static void
cdfWriteAttText(int ncfileid, int ncvarid, const char *name, std::string const &text)
{
  if (text.size()) nce(nc_put_att_text(ncfileid, ncvarid, name, text.size(), text.c_str()));
}

static void
cdfWriteVarInt(int ncfileid, int ncvarid, int *array)
{
  nce(nc_put_var_int(ncfileid, ncvarid, array));
}

static void
cdfWriteVarDouble(int ncfileid, int ncvarid, double *array)
{
  nce(nc_put_var_double(ncfileid, ncvarid, array));
}

static void
cdfWriteVarSize(int ncfileid, int ncvarid, nc_type sizetype, size_t len, size_t *array)
{
  if (len == 0) return;

  if (sizetype == NC_INT)
  {
    std::vector<int> iarray(len);
    for (size_t i = 0; i < len; ++i) iarray[i] = (int) array[i];
    nce(nc_put_var_int(ncfileid, ncvarid, iarray.data()));
  }
  else { nce(nc_put_var_ulonglong(ncfileid, ncvarid, (unsigned long long *) array)); }
}

static void
checkRemapFilesize(const RemapGridW &srcGrid, const RemapGridW &tgtGrid, const RemapVarsW &rv, int &writemode, nc_type &sizetype)
{
  size_t nlinks = rv.numLinks;
  size_t nele1 = 4 * 8 + 4 + srcGrid.ncorners * 2 * 8;
  size_t nele2 = 4 * 8 + 4 + tgtGrid.ncorners * 2 * 8;
  size_t filesize = srcGrid.ncells * nele1 + tgtGrid.ncells * nele2 + nlinks * (4 + 4 + rv.numWeights * 8);

  if (Options::cdoVerbose)
  {
    cdo_print("Number of remap links:       %zu", nlinks);
    cdo_print("Filesize for remap weights: ~%zu", filesize);
  }

  if (filesize > 0x7FFFFC00)  // 2**31 - 1024 (<2GB)
  {
    size_t maxlinks = 0x3FFFFFFF;  // 1GB
    auto gridsizemax = (srcGrid.ncells > tgtGrid.ncells) ? srcGrid.ncells : tgtGrid.ncells;
    if (nlinks > maxlinks || filesize > 8 * maxlinks || gridsizemax > 0x7FFFFC00)
    {
      if (Options::cdoVerbose) cdo_print("Store weights and links to NetCDF4!");
      writemode |= NC_NETCDF4;
      if (gridsizemax > 0x7FFFFC00)
        sizetype = NC_UINT64;
      else
        writemode |= NC_CLASSIC_MODEL;
    }
    else
    {
#ifdef NC_64BIT_OFFSET
      writemode |= NC_64BIT_OFFSET;
      if (Options::cdoVerbose) cdo_print("Store weights and links to NetCDF2!");
#else
      cdo_print("Filesize for remap weights maybe too large!");
#endif
    }
  }
}

struct CDFgrid
{
  int dims_id;
  int cntrlat_id;
  int cntrlon_id;
  int crnrlat_id;
  int crnrlon_id;
  int imask_id;
  int area_id;
  int frac_id;
};

static CDFgrid
define_remapgrid_scrip(int ncfileid, std::string const &prefix, nc_type sizetype, bool lgridarea, const RemapGridW &grid)
{
  CDFgrid cdfGrid;

  // Define grid size dimension
  int nc_size_id = cdf_def_dim(ncfileid, (prefix + "_size").c_str(), grid.ncells);
  // Define grid corner dimension
  int nc_corn_id = grid.ncorners ? cdf_def_dim(ncfileid, (prefix + "_corners").c_str(), grid.ncorners) : -1;
  // Define grid rank dimension
  int nc_rank_id = cdf_def_dim(ncfileid, (prefix + "_rank").c_str(), grid.rank);

  // Define grid dimensions
  cdfGrid.dims_id = cdf_def_var(ncfileid, (prefix + "_dims").c_str(), sizetype, 1, &nc_rank_id);
  // Define grid center latitude array
  cdfGrid.cntrlat_id = cdf_def_var(ncfileid, (prefix + "_center_lat").c_str(), NC_DOUBLE, 1, &nc_size_id);
  // Define grid center longitude array
  cdfGrid.cntrlon_id = cdf_def_var(ncfileid, (prefix + "_center_lon").c_str(), NC_DOUBLE, 1, &nc_size_id);

  // Define grid corner lat/lon arrays
  int nc_dims2_id[2] = { nc_size_id, nc_corn_id };
  cdfGrid.crnrlat_id = grid.ncorners ? cdf_def_var(ncfileid, (prefix + "_corner_lat").c_str(), NC_DOUBLE, 2, nc_dims2_id) : -1;
  cdfGrid.crnrlon_id = grid.ncorners ? cdf_def_var(ncfileid, (prefix + "_corner_lon").c_str(), NC_DOUBLE, 2, nc_dims2_id) : -1;

  // Define units for all coordinate arrays
  std::string gridunits = "radians";
  cdfWriteAttText(ncfileid, cdfGrid.cntrlat_id, "units", gridunits);
  cdfWriteAttText(ncfileid, cdfGrid.cntrlon_id, "units", gridunits);
  if (grid.ncorners)
  {
    cdfWriteAttText(ncfileid, cdfGrid.crnrlat_id, "units", gridunits);
    cdfWriteAttText(ncfileid, cdfGrid.crnrlon_id, "units", gridunits);
  }

  // Define grid mask
  cdfGrid.imask_id = cdf_def_var(ncfileid, (prefix + "_imask").c_str(), NC_INT, 1, &nc_size_id);
  cdfWriteAttText(ncfileid, cdfGrid.imask_id, "units", "unitless");

  // Define grid area array
  cdfGrid.area_id = -1;  // id for area of source grid cells
  if (lgridarea)
  {
    cdfGrid.area_id = cdf_def_var(ncfileid, (prefix + "_area").c_str(), NC_DOUBLE, 1, &nc_size_id);
    cdfWriteAttText(ncfileid, cdfGrid.area_id, "units", "square radians");
  }

  // Define grid fraction array
  cdfGrid.frac_id = cdf_def_var(ncfileid, (prefix + "_frac").c_str(), NC_DOUBLE, 1, &nc_size_id);
  cdfWriteAttText(ncfileid, cdfGrid.frac_id, "units", "unitless");

  return cdfGrid;
}

static void
write_remapgrid_scrip(int ncfileid, const CDFgrid &cdfGrid, bool lgridarea, RemapGridW &grid)
{
  int dims[2] = { (int) grid.dims[0], (int) grid.dims[1] };
  cdfWriteVarInt(ncfileid, cdfGrid.dims_id, dims);

  cdfWriteVarInt(ncfileid, cdfGrid.imask_id, grid.mask.data());

  if (grid.cell_center_lat.size()) cdfWriteVarDouble(ncfileid, cdfGrid.cntrlat_id, grid.cell_center_lat.data());
  if (grid.cell_center_lon.size()) cdfWriteVarDouble(ncfileid, cdfGrid.cntrlon_id, grid.cell_center_lon.data());

  if (grid.ncorners)
  {
    cdfWriteVarDouble(ncfileid, cdfGrid.crnrlat_id, grid.cell_corner_lat.data());
    cdfWriteVarDouble(ncfileid, cdfGrid.crnrlon_id, grid.cell_corner_lon.data());
  }

  if (lgridarea) cdfWriteVarDouble(ncfileid, cdfGrid.area_id, grid.cell_area.data());

  cdfWriteVarDouble(ncfileid, cdfGrid.frac_id, grid.cell_frac.data());
}

static std::string
remap_set_mapmethod(const RemapSwitches &remapSwitches)
{
  auto submapLAF = (remapSwitches.submapType == SubmapType::LAF);

  // clang-format off
  std::string mapMethod;
  switch (remapSwitches.mapType)
    {
    case RemapMethod::CONSERV:
      mapMethod = submapLAF ? "Largest area fraction" : "Conservative remapping";
      break;
    case RemapMethod::BILINEAR: mapMethod = "Bilinear remapping"; break;
    case RemapMethod::BICUBIC:  mapMethod = "Bicubic remapping"; break;
    case RemapMethod::KNN:
      mapMethod = (remapSwitches.numNeighbors == 1) ? "Nearest neighbor" : "Distance weighted avg of nearest neighbors";
      break;
    default:
      mapMethod = "unknown";
    }
  // clang-format on

  return mapMethod;
}

static void
write_remapfile_scrip(std::string const &remapFile, const RemapAttributes &remapAtts, const RemapSwitches &remapSwitches,
                      RemapGridW &srcGrid, RemapGridW &tgtGrid, RemapVarsW &rv)
{
  // Writes remap data to a NetCDF file using SCRIP conventions

#ifdef HAVE_LIBNETCDF

  int nc_dims2_id[2];  // NetCDF ids for 2d array dims

  auto lgridarea = (remapSwitches.mapType == RemapMethod::CONSERV);

  // if ( rv.numLinks == 0 ) cdo_abort("Number of remap links is 0, no remap weights found!");

  int writemode = NC_CLOBBER;
  auto sizetype = NC_INT;
  checkRemapFilesize(srcGrid, tgtGrid, rv, writemode, sizetype);

  // Create NetCDF file for mapping and define some global attributes
  int ncfileid = -1;
  nce(nc_create(remapFile.c_str(), writemode, &ncfileid));

  // Map name
  cdfWriteAttText(ncfileid, NC_GLOBAL, "title", remapAtts.map_name);

  // Normalization option
  std::string normalizeOptStr;  // character string for normalization option
  // clang-format off
  switch (rv.normOpt)
    {
    case NormOpt::NONE:     normalizeOptStr = "none";     break;
    case NormOpt::FRACAREA: normalizeOptStr = "fracarea"; break;
    case NormOpt::DESTAREA: normalizeOptStr = "destarea"; break;
    default: normalizeOptStr = "unknown";
    }
  // clang-format on
  cdfWriteAttText(ncfileid, NC_GLOBAL, "normalization", normalizeOptStr);

  // Map method
  auto mapMethod = remap_set_mapmethod(remapSwitches);
  cdfWriteAttText(ncfileid, NC_GLOBAL, "map_method", mapMethod);

  // Remap order
  if (remapSwitches.mapType == RemapMethod::CONSERV && remapSwitches.submapType == SubmapType::NONE)
    nce(nc_put_att_int(ncfileid, NC_GLOBAL, "remap_order", NC_INT, 1L, &remapSwitches.remapOrder));

  // File convention
  cdfWriteAttText(ncfileid, NC_GLOBAL, "conventions", "SCRIP");

  // Source and destination grid names
  cdfWriteAttText(ncfileid, NC_GLOBAL, "source_grid", srcGrid.name);
  cdfWriteAttText(ncfileid, NC_GLOBAL, "dest_grid", tgtGrid.name);

  // History
  cdfWriteAttText(ncfileid, NC_GLOBAL, "history", remapAtts.history);

  cdfWriteAttText(ncfileid, NC_GLOBAL, "CDO", remapAtts.cdo_version);

  // Define grids
  auto cdfSrcGrid = define_remapgrid_scrip(ncfileid, "src_grid", sizetype, lgridarea, srcGrid);
  auto cdfTgtGrid = define_remapgrid_scrip(ncfileid, "dst_grid", sizetype, lgridarea, tgtGrid);

  // Define map size dimensions
  auto nc_numlinks_id = cdf_def_dim(ncfileid, "numLinks", rv.numLinks);
  auto nc_numwgts_id = cdf_def_dim(ncfileid, "num_wgts", rv.numWeights);

  // Define mapping arrays
  auto nc_srcadd_id = cdf_def_var(ncfileid, "src_address", sizetype, 1, &nc_numlinks_id);
  auto nc_dstadd_id = cdf_def_var(ncfileid, "tgt_address", sizetype, 1, &nc_numlinks_id);

  nc_dims2_id[0] = nc_numlinks_id;
  nc_dims2_id[1] = nc_numwgts_id;
  int nc_rmpmatrix_id = cdf_def_var(ncfileid, "remap_matrix", NC_DOUBLE, 2, nc_dims2_id);

  // End definition stage

  nce(nc_enddef(ncfileid));

  // Write mapping data

  write_remapgrid_scrip(ncfileid, cdfSrcGrid, lgridarea, srcGrid);
  write_remapgrid_scrip(ncfileid, cdfTgtGrid, lgridarea, tgtGrid);

  for (size_t i = 0; i < rv.numLinks; ++i) rv.srcCellIndices[i]++;
  for (size_t i = 0; i < rv.numLinks; ++i) rv.tgtCellIndices[i]++;

  cdfWriteVarSize(ncfileid, nc_srcadd_id, sizetype, rv.numLinks, rv.srcCellIndices.data());
  cdfWriteVarSize(ncfileid, nc_dstadd_id, sizetype, rv.numLinks, rv.tgtCellIndices.data());
  cdfWriteVarDouble(ncfileid, nc_rmpmatrix_id, rv.weights.data());

  nce(nc_close(ncfileid));

#else
  cdo_abort("NetCDF support not compiled in!");
#endif
}

static void
check_areas(size_t n_a, Varray<double> const &area_a, Varray<double> const &area_b, size_t n_s, Varray<size_t> const &col,
            Varray<size_t> const &row, Varray<double> const &S)
{
  if (area_a.size() == 0 || area_b.size() == 0) return;

  Varray<double> sum(n_a, 0.0);

  // sum weighted ratio of true areas
  //      ’a’ is source; ’b’ is destination
  for (size_t i = 0; i < n_s; ++i)  // loop over all elements of S (the weights)
    sum[col[i]] = sum[col[i]] + S[i] * area_a[col[i]] / area_b[row[i]];

  // check that sums are equal to 1 (within tolerance of 1.e-6)
  for (size_t i = 0; i < n_a; ++i)  // loop over all source cells
  {
    if (std::fabs(sum[i] - 1.0) > 1.e-6)
    {
      cdo_print("srcCellIndex=%zu areaSum=%g", i + 1, sum[i]);
      break;
    }
  }
}

static void
verify_weights(std::string const &remapFile)
{
  RemapGridW srcGrid, tgtGrid;
  RemapVarsW remapVars;
  RemapAttributes remapAtts;

  (void) readRemapFileScrip(remapFile, remapAtts, srcGrid, tgtGrid, remapVars);

  check_areas(srcGrid.ncells, srcGrid.cell_area, tgtGrid.cell_area, remapVars.numLinks, remapVars.srcCellIndices,
              remapVars.tgtCellIndices, remapVars.weights);
}

static void
write_remap_scrip(std::string const &remapFileIn, std::string const &remapFileOut)
{
  RemapGridW srcGrid, tgtGrid;
  RemapVarsW remapVars;
  RemapAttributes remapAtts;

  auto remapSwitches = readRemapFileScrip(remapFileIn, remapAtts, srcGrid, tgtGrid, remapVars);

  /*
  char history[1024] = "date and time";
  time_t date_and_time_in_sec = time(NULL);
  if (date_and_time_in_sec != -1)
    {
      struct tm *date_and_time = localtime(&date_and_time_in_sec);
      (void) str::strftime(history, 1024, "%d %b %Y : ", date_and_time);
      std::strcat(history, cdo::command_line());
    }
    */
  if (Options::VersionInfo) remapAtts.cdo_version = cdo_comment();
  write_remapfile_scrip(remapFileOut, remapAtts, remapSwitches, srcGrid, tgtGrid, remapVars);
}
#endif

class Verifyweights : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Verifyweights",
    .operators = { { "verifyweights", 0, 1, "remap filename" }, { "writeremapscrip", 0, 2, "input and output remap filename" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<Verifyweights> registration = RegisterEntry<Verifyweights>(module);

public:
  void
  init() override
  {
#ifndef HAVE_LIBNETCDF
    cdo_abort("NetCDF support not compiled in!");
#else
    auto VERIFYWEIGHTS = module.get_id("verifyweights");
    auto WRITEREMAPSCRIP = module.get_id("writeremapscrip");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));
    auto numFiles = cdo_operator_f2(operatorID);
    operator_check_argc(numFiles);

    auto remapFileIn = cdo_operator_argv(0);

    if (operatorID == VERIFYWEIGHTS)
      verify_weights(remapFileIn);
    else if (operatorID == WRITEREMAPSCRIP)
      write_remap_scrip(remapFileIn, cdo_operator_argv(1));
#endif
  }

  void
  run() override
  {
  }

  void
  close() override
  {
  }
};
