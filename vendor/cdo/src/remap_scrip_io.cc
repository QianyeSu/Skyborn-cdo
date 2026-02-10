/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "knndata.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#include <ctime>

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_default_values.h"
#include "cdo_output.h"
#include "griddes.h"
#include <mpim_grid.h>
#include "remap.h"
#include "commandline.h"

constexpr size_t IndexLimit = 0x7FFFFC00;  // 2**31 - 1024 (<2GB)
static auto filterAvail{ false };

enum struct Compress
{
  NONE,
  ZIP,
  ZSTD
};

#ifdef HAVE_LIBNETCDF
static void
nce(int istat)
{
  // This routine provides a simple interface to NetCDF error message routine.
  if (istat != NC_NOERR) cdo_abort("%s", nc_strerror(istat));
}

static void
write_array_int64(int ncId, int ncVarId, nc_type xtype, size_t arraySize, size_t *array)
{
  if (arraySize == 0) return;

  if (xtype == NC_NAT) {}
#ifdef NC_UINT64
  else if (xtype == NC_UINT64) { nce(nc_put_var_ulonglong(ncId, ncVarId, (unsigned long long *) array)); }
#endif
  else
  {
    std::vector<int> arrayInt(arraySize);
    for (size_t i = 0; i < arraySize; ++i) arrayInt[i] = (int) array[i];
    nce(nc_put_var_int(ncId, ncVarId, arrayInt.data()));
  }
}

static void
read_array_int64(int ncId, int ncVarId, size_t arraySize, size_t *array)
{
  nc_type xtype = NC_NAT;
  nce(nc_inq_vartype(ncId, ncVarId, &xtype));

  if (xtype == NC_NAT) {}
#ifdef NC_UINT64
  else if (xtype == NC_UINT64) { nce(nc_get_var_ulonglong(ncId, ncVarId, (unsigned long long *) array)); }
#endif
  else
  {
    std::vector<int> arrayInt(arraySize);
    nce(nc_get_var_int(ncId, ncVarId, arrayInt.data()));
    for (size_t i = 0; i < arraySize; ++i) array[i] = (size_t) arrayInt[i];
  }
}
#endif

#ifdef HAVE_LIBNETCDF
static void
define_compression_zip(int ncId, int ncVarId)
{
  int deflateLevel = Options::cdoCompLevel;
  if (deflateLevel < 1 || deflateLevel > 9) deflateLevel = 1;

  int shuffle = 0, deflate = 1;
  nce(nc_def_var_deflate(ncId, ncVarId, shuffle, deflate, deflateLevel));
}
#endif

#ifdef HAVE_LIBNETCDF
static int
define_var(Compress compress, int ncId, const char *name, nc_type xtype, int ndims, const int *dimidsp)
{
  int ncVarId = -1;
  nce(nc_def_var(ncId, name, xtype, ndims, dimidsp, &ncVarId));
  // clang-format off
  if      (compress == Compress::ZIP)  define_compression_zip(ncId, ncVarId);
  else if (Options::filterSpec.size() && filterAvail) cdf_def_var_filter(ncId, ncVarId, Options::filterSpec.c_str());
  // clang-format on
  return ncVarId;
}
#endif

#ifdef HAVE_LIBNETCDF
static int
set_write_mode(bool need_src_cell_corners, size_t srcGridSize, size_t srcGridNC, bool need_tgt_cell_corners, size_t tgtGridSize,
               size_t tgtGridNC, size_t numLinks, size_t numWeights, nc_type &largeSizetype)
{
  int writemode = NC_CLOBBER;

  size_t nlinks = numLinks;
  size_t nele1 = 4 * 8 + 4;
  size_t nele2 = 4 * 8 + 4;
  if (need_src_cell_corners) nele1 += srcGridNC * 2 * 8;
  if (need_tgt_cell_corners) nele2 += tgtGridNC * 2 * 8;
  size_t filesize = srcGridSize * nele1 + tgtGridSize * nele2 + nlinks * (4 + 4 + numWeights * 8);

  if (Options::cdoVerbose)
  {
    cdo_print("Number of remap links:       %zu", nlinks);
    cdo_print("Filesize for remap weights: ~%zu", filesize);
  }

  if (filesize > IndexLimit)
  {
    constexpr size_t maxlinks = 0x3FFFFFFF;  // 1GB
    size_t gridsizeMax = std::max(srcGridSize, tgtGridSize);
    if (nlinks > maxlinks || filesize > 8 * maxlinks || gridsizeMax > IndexLimit)
    {
      if (Options::cdoVerbose) cdo_print("Store weights and links to NetCDF4!");
      writemode |= NC_NETCDF4;
      if (gridsizeMax > IndexLimit)
        largeSizetype = NC_UINT64;
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

  return writemode;
}

static std::string
set_map_method(const RemapSwitches &remapSwitches, bool &needGridarea)
{
  std::string mapMethod = "unknown";

  switch (remapSwitches.mapType)
  {
    case RemapMethod::CONSERV:
      needGridarea = true;
      mapMethod = (remapSwitches.submapType == SubmapType::LAF) ? "Largest area fraction" : "Conservative remapping";
      break;
    case RemapMethod::BILINEAR: mapMethod = "Bilinear remapping"; break;
    case RemapMethod::BICUBIC: mapMethod = "Bicubic remapping"; break;
    case RemapMethod::KNN:
      if (remapSwitches.numNeighbors == -1) { mapMethod = "k-nearest neighbor"; }
      else { mapMethod = (remapSwitches.numNeighbors == 1) ? "Nearest neighbor" : "Distance weighted avg of nearest neighbors"; }
      break;
    case RemapMethod::UNDEF: break;
  }

  return mapMethod;
}

static void
put_att_text(int ncId, int ncVarId, const char *name, const char *text)
{
  nce(nc_put_att_text(ncId, ncVarId, name, std::strlen(text), text));
}

static int
def_dim(int ncId, const char *name, size_t len)
{
  int dimId = -1;
  nce(nc_def_dim(ncId, name, len, &dimId));
  return dimId;
}

static int
inq_varid(int ncId, const char *name)
{
  int varId = -1;
  auto status = nc_inq_varid(ncId, name, &varId);
  if (status != NC_NOERR) cdo_warning("var name not found: %s", name);
  nce(status);
  return varId;
}

static size_t
inq_dimlen(int ncId, const char *name)
{
  size_t dimlen = 0;
  int dimId = -1;
  auto status = nc_inq_dimid(ncId, name, &dimId);
  if (status != NC_NOERR) cdo_warning("dim name not found: %s", name);
  nce(status);
  nce(nc_inq_dimlen(ncId, dimId, &dimlen));
  return dimlen;
}

#endif

void
remap_write_data_scrip(std::string const &weightsfile, KnnParams const &knnParams, const RemapSwitches &remapSwitches,
                       RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv)
{
  // Writes remap data to a NetCDF file using SCRIP conventions

#ifdef HAVE_LIBNETCDF

  const char *normalizeOpt = "unknown";
  switch (rv.normOpt)
  {
    case NormOpt::NONE: normalizeOpt = "none"; break;
    case NormOpt::FRACAREA: normalizeOpt = "fracarea"; break;
    case NormOpt::DESTAREA: normalizeOpt = "destarea"; break;
  }

  // if (rv.numLinks == 0) cdo_abort("Number of remap links is 0, no remap weights found!");

  auto largeSizetype = NC_INT;
  auto writemode = set_write_mode(srcGrid.needCellCorners, srcGrid.size, srcGrid.numCorners, tgtGrid.needCellCorners, tgtGrid.size,
                                  tgtGrid.numCorners, rv.numLinks, rv.numWeights, largeSizetype);

  auto srcSizetype = (srcGrid.size > IndexLimit) ? largeSizetype : NC_INT;
  auto dstSizetype = (tgtGrid.size > IndexLimit) ? largeSizetype : NC_INT;

  auto compress{ Compress::NONE };
  if (CdoDefault::FileType == CDI_FILETYPE_NC4 || CdoDefault::FileType == CDI_FILETYPE_NC4C)
  {
    filterAvail = true;
    writemode = NC_CLOBBER | NC_NETCDF4;
    if (CdoDefault::FileType == CDI_FILETYPE_NC4C && largeSizetype != NC_UINT64) writemode |= NC_CLASSIC_MODEL;
    if (Options::cdoCompType == CDI_COMPRESS_ZIP) compress = Compress::ZIP;
  }

  // Create NetCDF file for mapping and define some global attributes
  int ncId = -1;
  nce(nc_create(weightsfile.c_str(), writemode, &ncId));

  // Map name
  const char *mapName = "CDO remapping";
  put_att_text(ncId, NC_GLOBAL, "title", mapName);

  // Normalization option
  put_att_text(ncId, NC_GLOBAL, "normalization", normalizeOpt);

  // Map method
  auto needGridarea = false;
  auto mapMethod = set_map_method(remapSwitches, needGridarea);
  put_att_text(ncId, NC_GLOBAL, "map_method", mapMethod.c_str());

  auto storeKnnParams = (remapSwitches.mapType == RemapMethod::KNN && mapMethod.starts_with("k-nearest"));

  if (remapSwitches.mapType == RemapMethod::KNN && remapSwitches.numNeighbors > 1 && remapSwitches.numNeighbors != 4)
    nce(nc_put_att_int(ncId, NC_GLOBAL, "num_neighbors", NC_INT, 1L, &remapSwitches.numNeighbors));

  // Remap order
  if (remapSwitches.mapType == RemapMethod::CONSERV && remapSwitches.submapType == SubmapType::NONE)
    nce(nc_put_att_int(ncId, NC_GLOBAL, "remap_order", NC_INT, 1L, &remapSwitches.remapOrder));

  // File convention
  put_att_text(ncId, NC_GLOBAL, "conventions", "SCRIP");

  // Source and destination grid names
  put_att_text(ncId, NC_GLOBAL, "source_grid", srcGrid.name.c_str());
  put_att_text(ncId, NC_GLOBAL, "dest_grid", tgtGrid.name.c_str());

  // History
  auto dateAndTimeInSec = std::time(NULL);
  if (dateAndTimeInSec != -1)
  {
    char history[1024] = "date and time";
    struct tm *dateAndTime = std::localtime(&dateAndTimeInSec);
    (void) std::strftime(history, 1024, "%d %b %Y : ", dateAndTime);
    std::strcat(history, cdo::command_line());
    put_att_text(ncId, NC_GLOBAL, "history", history);
  }

  if (Options::VersionInfo) put_att_text(ncId, NC_GLOBAL, "CDO", cdo_comment());

  // Prepare NetCDF dimension info

  // Define grid size dimensions
  auto nc_srcgrdsize_id = def_dim(ncId, "src_grid_size", srcGrid.size);
  auto nc_dstgrdsize_id = def_dim(ncId, "dst_grid_size", tgtGrid.size);

  // Define grid corner dimension
  auto nc_srcgrdcorn_id = (srcGrid.needCellCorners) ? def_dim(ncId, "src_grid_corners", srcGrid.numCorners) : -1;
  auto nc_dstgrdcorn_id = (tgtGrid.needCellCorners) ? def_dim(ncId, "dst_grid_corners", tgtGrid.numCorners) : -1;

  // Define grid rank dimension
  auto nc_srcgrdrank_id = def_dim(ncId, "src_grid_rank", srcGrid.rank);
  auto nc_dstgrdrank_id = def_dim(ncId, "dst_grid_rank", tgtGrid.rank);

  // Define map size dimensions
  auto nc_numlinks_id = def_dim(ncId, "num_links", rv.numLinks);
  auto nc_numwgts_id = def_dim(ncId, "num_wgts", rv.numWeights);

  // Define grid dimensions

  auto srcDimsXtype = (srcGrid.dims[0] > IndexLimit) ? largeSizetype : NC_INT;
  auto dstDimsXtype = (tgtGrid.dims[0] > IndexLimit) ? largeSizetype : NC_INT;

  int nc_srcgrddims_id = -1, nc_dstgrddims_id = -1;
  nce(nc_def_var(ncId, "src_grid_dims", srcDimsXtype, 1, &nc_srcgrdrank_id, &nc_srcgrddims_id));
  nce(nc_def_var(ncId, "dst_grid_dims", dstDimsXtype, 1, &nc_dstgrdrank_id, &nc_dstgrddims_id));

  // Define all arrays for NetCDF descriptors

  // Define grid center latitude array
  auto nc_srcgrdcntrlat_id = define_var(compress, ncId, "src_grid_center_lat", NC_DOUBLE, 1, &nc_srcgrdsize_id);
  auto nc_dstgrdcntrlat_id = define_var(compress, ncId, "dst_grid_center_lat", NC_DOUBLE, 1, &nc_dstgrdsize_id);

  // Define grid center longitude array
  auto nc_srcgrdcntrlon_id = define_var(compress, ncId, "src_grid_center_lon", NC_DOUBLE, 1, &nc_srcgrdsize_id);
  auto nc_dstgrdcntrlon_id = define_var(compress, ncId, "dst_grid_center_lon", NC_DOUBLE, 1, &nc_dstgrdsize_id);

  // Define grid corner lat/lon arrays

  int nc_dims2_id[2];  // NetCDF ids for 2d array dims
  nc_dims2_id[0] = nc_srcgrdsize_id;
  nc_dims2_id[1] = nc_srcgrdcorn_id;

  int nc_srcgrdcrnrlat_id = -1, nc_srcgrdcrnrlon_id = -1;
  if (srcGrid.needCellCorners)
  {
    nc_srcgrdcrnrlat_id = define_var(compress, ncId, "src_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id);
    nc_srcgrdcrnrlon_id = define_var(compress, ncId, "src_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id);
  }

  nc_dims2_id[0] = nc_dstgrdsize_id;
  nc_dims2_id[1] = nc_dstgrdcorn_id;

  int nc_dstgrdcrnrlat_id = -1, nc_dstgrdcrnrlon_id = -1;
  if (tgtGrid.needCellCorners)
  {
    nc_dstgrdcrnrlat_id = define_var(compress, ncId, "dst_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id);
    nc_dstgrdcrnrlon_id = define_var(compress, ncId, "dst_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id);
  }

  // Define units for all coordinate arrays

  const char *srcGridUnits = "radians";
  const char *tgtGridUnits = "radians";
  put_att_text(ncId, nc_srcgrdcntrlat_id, "units", srcGridUnits);
  put_att_text(ncId, nc_dstgrdcntrlat_id, "units", tgtGridUnits);
  put_att_text(ncId, nc_srcgrdcntrlon_id, "units", srcGridUnits);
  put_att_text(ncId, nc_dstgrdcntrlon_id, "units", tgtGridUnits);
  if (srcGrid.needCellCorners) put_att_text(ncId, nc_srcgrdcrnrlat_id, "units", srcGridUnits);
  if (srcGrid.needCellCorners) put_att_text(ncId, nc_srcgrdcrnrlon_id, "units", srcGridUnits);
  if (tgtGrid.needCellCorners) put_att_text(ncId, nc_dstgrdcrnrlat_id, "units", tgtGridUnits);
  if (tgtGrid.needCellCorners) put_att_text(ncId, nc_dstgrdcrnrlon_id, "units", tgtGridUnits);

  // Define grid mask

  auto nc_srcgrdimask_id = define_var(compress, ncId, "src_grid_imask", NC_INT, 1, &nc_srcgrdsize_id);
  put_att_text(ncId, nc_srcgrdimask_id, "units", "unitless");

  auto nc_dstgrdimask_id = define_var(compress, ncId, "dst_grid_imask", NC_INT, 1, &nc_dstgrdsize_id);
  put_att_text(ncId, nc_dstgrdimask_id, "units", "unitless");

  // Define grid area arrays

  int nc_srcgrdarea_id = -1, nc_dstgrdarea_id = -1;
  if (needGridarea)
  {
    nc_srcgrdarea_id = define_var(compress, ncId, "src_grid_area", NC_DOUBLE, 1, &nc_srcgrdsize_id);
    put_att_text(ncId, nc_srcgrdarea_id, "units", "square radians");

    nc_dstgrdarea_id = define_var(compress, ncId, "dst_grid_area", NC_DOUBLE, 1, &nc_dstgrdsize_id);
    put_att_text(ncId, nc_dstgrdarea_id, "units", "square radians");
  }

  // Define grid fraction arrays

  auto nc_srcgrdfrac_id = define_var(compress, ncId, "src_grid_frac", NC_DOUBLE, 1, &nc_srcgrdsize_id);
  put_att_text(ncId, nc_srcgrdfrac_id, "units", "unitless");

  auto nc_dstgrdfrac_id = define_var(compress, ncId, "dst_grid_frac", NC_DOUBLE, 1, &nc_dstgrdsize_id);
  put_att_text(ncId, nc_dstgrdfrac_id, "units", "unitless");

  // Define mapping arrays

  auto nc_srcadd_id = define_var(compress, ncId, "src_address", srcSizetype, 1, &nc_numlinks_id);
  auto nc_dstadd_id = define_var(compress, ncId, "dst_address", dstSizetype, 1, &nc_numlinks_id);

  nc_dims2_id[0] = nc_numlinks_id;
  nc_dims2_id[1] = nc_numwgts_id;

  auto nc_rmpmatrix_id = define_var(compress, ncId, "remap_matrix", NC_DOUBLE, 2, nc_dims2_id);

  if (storeKnnParams)
  {
    auto knn_params_id = define_var(Compress::NONE, ncId, "knn_params", NC_INT, 0, NULL);
    int k = knnParams.k;
    int kMin = knnParams.kMin;
    auto weightingMethodStr = weightingMethod_to_string(knnParams.weighted);
    nc_put_att_int(ncId, knn_params_id, "k", NC_INT, 1, &k);
    if (kMin > 0) { nc_put_att_int(ncId, knn_params_id, "kmin", NC_INT, 1, &kMin); }
    put_att_text(ncId, knn_params_id, "weighted", weightingMethodStr.c_str());
    if (knnParams.weighted == WeightingMethod::gaussWeighted)
    {
      nc_put_att_double(ncId, knn_params_id, "gauss_scale", NC_DOUBLE, 1, &knnParams.gaussScale);
    }
  }

  // End definition stage

  nce(nc_enddef(ncId));

  // Write mapping data

  write_array_int64(ncId, nc_srcgrddims_id, srcDimsXtype, 2, &srcGrid.dims[0]);
  write_array_int64(ncId, nc_dstgrddims_id, dstDimsXtype, 2, &tgtGrid.dims[0]);

  nce(nc_put_var_schar(ncId, nc_srcgrdimask_id, &srcGrid.mask[0]));
  nce(nc_put_var_schar(ncId, nc_dstgrdimask_id, &tgtGrid.mask[0]));

  if (!srcGrid.centerLats.empty()) nce(nc_put_var_double(ncId, nc_srcgrdcntrlat_id, srcGrid.centerLats.data()));
  if (!srcGrid.centerLons.empty()) nce(nc_put_var_double(ncId, nc_srcgrdcntrlon_id, srcGrid.centerLons.data()));

  if (srcGrid.needCellCorners)
  {
    nce(nc_put_var_double(ncId, nc_srcgrdcrnrlat_id, srcGrid.cornerLats.data()));
    nce(nc_put_var_double(ncId, nc_srcgrdcrnrlon_id, srcGrid.cornerLons.data()));
  }

  if (!tgtGrid.centerLats.empty()) nce(nc_put_var_double(ncId, nc_dstgrdcntrlat_id, tgtGrid.centerLats.data()));
  if (!tgtGrid.centerLons.empty()) nce(nc_put_var_double(ncId, nc_dstgrdcntrlon_id, tgtGrid.centerLons.data()));

  if (tgtGrid.needCellCorners)
  {
    nce(nc_put_var_double(ncId, nc_dstgrdcrnrlat_id, tgtGrid.cornerLats.data()));
    nce(nc_put_var_double(ncId, nc_dstgrdcrnrlon_id, tgtGrid.cornerLons.data()));
  }

  if (needGridarea) nce(nc_put_var_double(ncId, nc_srcgrdarea_id, &srcGrid.cellArea[0]));

  nce(nc_put_var_double(ncId, nc_srcgrdfrac_id, &srcGrid.cellFrac[0]));

  /*
  if (luse_cell_area)
    nce(nc_put_var_double(ncId, nc_dstgrdarea_id, tgtGrid.cell_area_in));
  else
  */
  if (needGridarea) nce(nc_put_var_double(ncId, nc_dstgrdarea_id, &tgtGrid.cellArea[0]));

  nce(nc_put_var_double(ncId, nc_dstgrdfrac_id, &tgtGrid.cellFrac[0]));

  if (rv.numLinks > 0)
  {
    for (size_t i = 0; i < rv.numLinks; ++i) rv.srcCellIndices[i]++;
    for (size_t i = 0; i < rv.numLinks; ++i) rv.tgtCellIndices[i]++;

    write_array_int64(ncId, nc_srcadd_id, srcSizetype, rv.numLinks, &rv.srcCellIndices[0]);
    write_array_int64(ncId, nc_dstadd_id, dstSizetype, rv.numLinks, &rv.tgtCellIndices[0]);

    nce(nc_put_var_double(ncId, nc_rmpmatrix_id, &rv.weights[0]));
  }

  nce(nc_close(ncId));

#else
  cdo_abort("NetCDF support not compiled in!");
#endif

}  // remap_write_data_scrip

/*****************************************************************************/

#ifdef HAVE_LIBNETCDF
static std::string
get_text_attribute(int ncId, int att_id, const char *att_name)
{
  char cstr[1024];
  nce(nc_get_att_text(ncId, att_id, att_name, cstr));
  size_t attlen;
  nce(nc_inq_attlen(ncId, att_id, att_name, &attlen));
  cstr[attlen] = 0;

  return std::string(cstr);
}

KnnParams
read_knn_params(int ncId)
{
  KnnParams knnParams;
  auto nc_knn_params_id = inq_varid(ncId, "knn_params");
  int k = 0;
  int status = nc_get_att_int(ncId, nc_knn_params_id, "k", &k);
  if (status == NC_NOERR && k > 0) knnParams.k = k;
  int kMin = 0;
  status = nc_get_att_int(ncId, nc_knn_params_id, "kmin", &kMin);
  if (status == NC_NOERR && kMin > 0) knnParams.kMin = kMin;
  auto weightingMethodStr = get_text_attribute(ncId, nc_knn_params_id, "weighted");
  knnParams.weighted = string_to_weightingMethod(weightingMethodStr);
  if (knnParams.weighted == WeightingMethod::gaussWeighted)
  {
    double gaussScale = 1.0;
    status = nc_get_att_double(ncId, nc_knn_params_id, "gauss_scale", &gaussScale);
    if (status == NC_NOERR) knnParams.gaussScale = gaussScale;
  }
  return knnParams;
}

RemapSwitches
get_maptype(int ncId)
{
  RemapSwitches remapSwitches;
  remapSwitches.remapOrder = 1;

  // Map method
  auto mapMethod = get_text_attribute(ncId, NC_GLOBAL, "map_method");

  if (mapMethod.starts_with("Conservative"))
  {
    remapSwitches.mapType = RemapMethod::CONSERV;
    int iatt;
    int status = nc_get_att_int(ncId, NC_GLOBAL, "remap_order", &iatt);
    if (status == NC_NOERR) remapSwitches.remapOrder = iatt;
  }
  else if (mapMethod.starts_with("Bilinear")) { remapSwitches.mapType = RemapMethod::BILINEAR; }
  else if (mapMethod.starts_with("Bicubic")) { remapSwitches.mapType = RemapMethod::BICUBIC; }
  else if (mapMethod.starts_with("k-nearest"))
  {
    remapSwitches.mapType = RemapMethod::KNN;
    remapSwitches.numNeighbors = -1;
  }
  else if (mapMethod.starts_with("Nearest"))
  {
    remapSwitches.mapType = RemapMethod::KNN;
    remapSwitches.numNeighbors = 1;
  }
  else if (mapMethod.starts_with("Distance"))
  {
    int numNeighbors = 4;
    int iatt;
    int status = nc_get_att_int(ncId, NC_GLOBAL, "num_neighbors", &iatt);
    if (status == NC_NOERR && iatt > 0) numNeighbors = iatt;

    remapSwitches.mapType = RemapMethod::KNN;
    remapSwitches.numNeighbors = numNeighbors;
  }
  else if (mapMethod.starts_with("Largest"))
  {
    remapSwitches.mapType = RemapMethod::CONSERV;
    remapSwitches.submapType = SubmapType::LAF;
  }
  else
  {
    cdo_print("mapType='%s'", mapMethod);
    cdo_abort("Invalid Map Type");
  }

  if (Options::cdoVerbose) cdo_print("mapType='%s'", mapMethod);

  return remapSwitches;
}
#endif

RemapSwitches
remap_read_data_scrip(std::string const &weightsfile, int gridID1, int gridID2, RemapGrid &srcGrid, RemapGrid &tgtGrid,
                      RemapVars &rv, KnnParams &knnParams)
{
  RemapSwitches remapSwitches;

  // The routine reads a NetCDF file to extract remapping info in SCRIP format

#ifdef HAVE_LIBNETCDF

  constexpr bool readGridCoordinates = false;
  int status;
  size_t dimlen;

  // Open file and read some global information

  auto ncId = cdo_cdf_openread(weightsfile.c_str());

  // Map name
  auto mapName = get_text_attribute(ncId, NC_GLOBAL, "title");

  if (Options::cdoVerbose)
  {
    cdo_print("Reading remap weights: %s", mapName);
    cdo_print("From file: %s", weightsfile);
  }

  // Map Tyoe
  remapSwitches = get_maptype(ncId);
  if (remapSwitches.mapType == RemapMethod::KNN && remapSwitches.numNeighbors == -1) knnParams = read_knn_params(ncId);

  auto needGridarea = false;

  remap_vars_init(remapSwitches.mapType, remapSwitches.remapOrder, rv);

  rv.mapType = remapSwitches.mapType;
  rv.numLinksPerValue = -1;

  // Normalization option
  auto normalizeOpt = get_text_attribute(ncId, NC_GLOBAL, "normalization");

  // clang-format off
  if      (normalizeOpt == "none")     rv.normOpt = NormOpt::NONE;
  else if (normalizeOpt == "fracarea") rv.normOpt = NormOpt::FRACAREA;
  else if (normalizeOpt == "destarea") rv.normOpt = NormOpt::DESTAREA;
  else
    {
      cdo_print("normalize_opt = %s", normalizeOpt);
      cdo_abort("Invalid normalization option");
    }
  // clang-format on

  if (Options::cdoVerbose) cdo_print("normalize_opt = %s", normalizeOpt);

  // File convention
  auto convention = get_text_attribute(ncId, NC_GLOBAL, "conventions");

  if (convention != "SCRIP")
  {
    cdo_print("convention = %s", convention);
    cdo_abort("%s file convention!", (convention == "NCAR-CSM") ? "Unsupported" : "Unknown");
  }

  // Read some additional global attributes

  // Source and destination grid names

  srcGrid.name = get_text_attribute(ncId, NC_GLOBAL, "source_grid");
  tgtGrid.name = get_text_attribute(ncId, NC_GLOBAL, "dest_grid");

  if (Options::cdoVerbose) cdo_print("Remapping between: %s and %s", srcGrid.name, tgtGrid.name);

  // Read dimension information
  srcGrid.size = inq_dimlen(ncId, "src_grid_size");
  // if (srcGrid.size != gridInqSize(gridID1)) cdo_abort("Source grids have different size!");

  tgtGrid.size = inq_dimlen(ncId, "dst_grid_size");

  // if (tgtGrid.size != gridInqSize(gridID2)) cdo_abort("Target grids have different size!");

  int nc_srcgrdcorn_id;
  status = nc_inq_dimid(ncId, "src_grid_corners", &nc_srcgrdcorn_id);
  if (status == NC_NOERR)
  {
    nce(nc_inq_dimlen(ncId, nc_srcgrdcorn_id, &dimlen));
    srcGrid.numCorners = dimlen;
    srcGrid.useCellCorners = true;
    srcGrid.needCellCorners = true;
  }

  int nc_dstgrdcorn_id;
  status = nc_inq_dimid(ncId, "dst_grid_corners", &nc_dstgrdcorn_id);
  if (status == NC_NOERR)
  {
    nce(nc_inq_dimlen(ncId, nc_dstgrdcorn_id, &dimlen));
    tgtGrid.numCorners = dimlen;
    tgtGrid.useCellCorners = true;
    tgtGrid.needCellCorners = true;
  }

  srcGrid.rank = inq_dimlen(ncId, "src_grid_rank");
  tgtGrid.rank = inq_dimlen(ncId, "dst_grid_rank");

  rv.numLinks = inq_dimlen(ncId, "num_links");
  // if (rv.numLinks == 0) cdo_abort("Number of remap links is 0, no remap weights found!");

  rv.numWeights = inq_dimlen(ncId, "num_wgts");

  srcGrid.gridID = gridID1;
  tgtGrid.gridID = gridID2;

  int gridID1_gme_c = -1;
  if (gridInqType(gridID1) == GRID_GME)
  {
    srcGrid.nvgp = gridInqSize(gridID1);
    gridID1_gme_c = gridToUnstructured(gridID1, NeedCorners::Yes);
  }

  if (readGridCoordinates)
  {
    remap_grid_alloc(rv.mapType, srcGrid);
    remap_grid_alloc(rv.mapType, tgtGrid);

    if (gridInqType(gridID1) == GRID_GME) gridInqMaskGME(gridID1_gme_c, &srcGrid.vgpm[0]);
  }
  else
  {
    srcGrid.mask.resize(srcGrid.size);
    tgtGrid.mask.resize(tgtGrid.size);
    if (remapSwitches.mapType == RemapMethod::CONSERV) tgtGrid.cellFrac.resize(tgtGrid.size);
  }

  // Allocate address and weight arrays for mapping 1
  if (rv.numLinks > 0)
  {
    rv.maxLinks = rv.numLinks;
    rv.srcCellIndices.resize(rv.numLinks);
    rv.tgtCellIndices.resize(rv.numLinks);
    rv.weights.resize(rv.numWeights * rv.numLinks);
  }

  // Get variable ids

  auto nc_srcgrddims_id = inq_varid(ncId, "src_grid_dims");
  auto nc_srcgrdimask_id = inq_varid(ncId, "src_grid_imask");
  auto nc_srcgrdcntrlat_id = inq_varid(ncId, "src_grid_center_lat");
  auto nc_srcgrdcntrlon_id = inq_varid(ncId, "src_grid_center_lon");

  auto nc_srcgrdcrnrlat_id = (srcGrid.numCorners) ? inq_varid(ncId, "src_grid_corner_lat") : -1;
  auto nc_srcgrdcrnrlon_id = (srcGrid.numCorners) ? inq_varid(ncId, "src_grid_corner_lon") : -1;

  auto nc_srcgrdarea_id = (needGridarea) ? inq_varid(ncId, "src_grid_area") : -1;
  auto nc_srcgrdfrac_id = inq_varid(ncId, "src_grid_frac");

  auto nc_dstgrddims_id = inq_varid(ncId, "dst_grid_dims");
  auto nc_dstgrdimask_id = inq_varid(ncId, "dst_grid_imask");
  auto nc_dstgrdcntrlat_id = inq_varid(ncId, "dst_grid_center_lat");
  auto nc_dstgrdcntrlon_id = inq_varid(ncId, "dst_grid_center_lon");

  auto nc_dstgrdcrnrlat_id = (tgtGrid.numCorners) ? inq_varid(ncId, "dst_grid_corner_lat") : -1;
  auto nc_dstgrdcrnrlon_id = (tgtGrid.numCorners) ? inq_varid(ncId, "dst_grid_corner_lon") : -1;

  auto nc_dstgrdarea_id = (needGridarea) ? inq_varid(ncId, "dst_grid_area") : -1;
  auto nc_dstgrdfrac_id = inq_varid(ncId, "dst_grid_frac");

  auto nc_srcadd_id = inq_varid(ncId, "src_address");
  auto nc_dstadd_id = inq_varid(ncId, "dst_address");
  auto nc_rmpmatrix_id = inq_varid(ncId, "remap_matrix");

  // Read all variables

  read_array_int64(ncId, nc_srcgrddims_id, 2, &srcGrid.dims[0]);

  nce(nc_get_var_schar(ncId, nc_srcgrdimask_id, &srcGrid.mask[0]));

  if (readGridCoordinates)
  {
    nce(nc_get_var_double(ncId, nc_srcgrdcntrlat_id, srcGrid.centerLats.data()));
    nce(nc_get_var_double(ncId, nc_srcgrdcntrlon_id, srcGrid.centerLons.data()));

    auto srcGridUnits = get_text_attribute(ncId, nc_srcgrdcntrlat_id, "units");
    if (string_to_LonLatUnits(srcGridUnits, "source grid center lon") == LonLatUnits::Deg) grid_to_radian(srcGrid.centerLons);
    if (string_to_LonLatUnits(srcGridUnits, "source grid center lat") == LonLatUnits::Deg) grid_to_radian(srcGrid.centerLats);

    if (srcGrid.numCorners)
    {
      nce(nc_get_var_double(ncId, nc_srcgrdcrnrlat_id, srcGrid.cornerLats.data()));
      nce(nc_get_var_double(ncId, nc_srcgrdcrnrlon_id, srcGrid.cornerLons.data()));

      srcGridUnits = get_text_attribute(ncId, nc_srcgrdcrnrlat_id, "units");
      if (string_to_LonLatUnits(srcGridUnits, "source grid corner lon") == LonLatUnits::Deg) grid_to_radian(srcGrid.cornerLons);
      if (string_to_LonLatUnits(srcGridUnits, "source grid corner lat") == LonLatUnits::Deg) grid_to_radian(srcGrid.cornerLats);
    }

    if (needGridarea) nce(nc_get_var_double(ncId, nc_srcgrdarea_id, &srcGrid.cellArea[0]));

    nce(nc_get_var_double(ncId, nc_srcgrdfrac_id, &srcGrid.cellFrac[0]));
  }

  read_array_int64(ncId, nc_dstgrddims_id, 2, &tgtGrid.dims[0]);

  nce(nc_get_var_schar(ncId, nc_dstgrdimask_id, &tgtGrid.mask[0]));

  if (readGridCoordinates)
  {
    nce(nc_get_var_double(ncId, nc_dstgrdcntrlat_id, tgtGrid.centerLats.data()));
    nce(nc_get_var_double(ncId, nc_dstgrdcntrlon_id, tgtGrid.centerLons.data()));

    auto tgtGridUnits = get_text_attribute(ncId, nc_dstgrdcntrlat_id, "units");
    if (string_to_LonLatUnits(tgtGridUnits, "target grid center lon") == LonLatUnits::Deg) grid_to_radian(tgtGrid.centerLons);
    if (string_to_LonLatUnits(tgtGridUnits, "target grid center lat") == LonLatUnits::Deg) grid_to_radian(tgtGrid.centerLats);

    if (tgtGrid.numCorners)
    {
      nce(nc_get_var_double(ncId, nc_dstgrdcrnrlat_id, tgtGrid.cornerLats.data()));
      nce(nc_get_var_double(ncId, nc_dstgrdcrnrlon_id, tgtGrid.cornerLons.data()));

      tgtGridUnits = get_text_attribute(ncId, nc_dstgrdcrnrlat_id, "units");
      if (string_to_LonLatUnits(tgtGridUnits, "target grid corner lon") == LonLatUnits::Deg) grid_to_radian(tgtGrid.cornerLons);
      if (string_to_LonLatUnits(tgtGridUnits, "target grid corner lat") == LonLatUnits::Deg) grid_to_radian(tgtGrid.cornerLats);
    }

    if (needGridarea) nce(nc_get_var_double(ncId, nc_dstgrdarea_id, &tgtGrid.cellArea[0]));
  }

  if (remapSwitches.mapType == RemapMethod::CONSERV) nce(nc_get_var_double(ncId, nc_dstgrdfrac_id, &tgtGrid.cellFrac[0]));

  if (rv.numLinks > 0)
  {
    read_array_int64(ncId, nc_srcadd_id, rv.numLinks, &rv.srcCellIndices[0]);
    read_array_int64(ncId, nc_dstadd_id, rv.numLinks, &rv.tgtCellIndices[0]);

    for (size_t i = 0; i < rv.numLinks; ++i) rv.srcCellIndices[i]--;
    for (size_t i = 0; i < rv.numLinks; ++i) rv.tgtCellIndices[i]--;

    nce(nc_get_var_double(ncId, nc_rmpmatrix_id, &rv.weights[0]));
  }

  // Close input file
  cdo_cdf_close(ncId);

#else
  cdo_abort("NetCDF support not compiled in!");
#endif

  return remapSwitches;
}  // remap_read_data_scrip
