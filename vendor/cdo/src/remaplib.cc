/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_timer.h"
#include "cdo_output.h"
#include "cdo_omp.h"
#include <mpim_grid.h>
#include "gridreference.h"
#include "remap.h"

static bool RemapGenWeights = true;
static bool RemapWriteRemap = false;

constexpr size_t DefaultMaxIteration = 100;
size_t RemapMaxIteration = DefaultMaxIteration;  // Max iteration count for i, j iteration

static void
print_healpix_warning_and_abort(void)
{
  cdo_warning("HEALPix cell edges are not located on great circles. This leads to inaccuracies in the calculation of "
              "cell areas.");
  cdo_abort("Use the CDO option --force to perform this operation.");
}

int
remap_check_mask_indices(const size_t (&indices)[4], Vmask const &mask)
{
  int searchResult = 1;
  if (mask.size() > 0)
  {
    for (int i = 0; i < 4; ++i)
    {
      if (mask[indices[i]] == 0) { searchResult = 0; }
    }
  }
  return searchResult;
}

void
remap_set_int(int remapvar, int value)
{
  // clang-format off
  if      (remapvar == REMAP_WRITE_REMAP)   RemapWriteRemap = (value > 0);
  else if (remapvar == REMAP_GENWEIGHTS)    RemapGenWeights = (value > 0);
  else if (remapvar == REMAP_MAX_ITER)      RemapMaxIteration = value;
  else    cdo_abort("Unsupported remap variable (%d)!", remapvar);
  // clang-format on
}

PointLonLat
remapgrid_get_lonlat(const RemapGrid *grid, size_t index)
{
  double lon{}, lat{};

  if (grid->type == RemapGridType::Reg2D)
  {
    auto nx = grid->dims[0];
    auto iy = index / nx;
    auto ix = index - iy * nx;
    lat = grid->centerLatsReg2d[iy];
    lon = grid->centerLonsReg2d[ix];
    if (lon < 0) lon += PI2;
  }
  else if (grid->type == RemapGridType::HealPix)
  {
    hp_index_to_lonlat(grid->hpParams, index, &lon, &lat);
    // if (lon < 0) lon += PI2;
  }
  else
  {
    lat = grid->centerLats[index];
    lon = grid->centerLons[index];
  }

  return PointLonLat(lon, lat);
}

void
check_lon_range(const char *txt, size_t nlons, Varray<double> &lons)
{
  assert(!lons.empty());

  if (txt)
  {
    double minval = 1.e36;
    double maxval = -1.e36;
#ifdef _OPENMP
#pragma omp parallel for if (nlons > cdoMinLoopSize) default(shared) schedule(static) reduction(min : minval) \
    reduction(max : maxval)
#endif
    for (size_t i = 0; i < nlons; ++i)
    {
      minval = std::min(minval, lons[i]);
      maxval = std::max(maxval, lons[i]);
    }
    if (minval < -PI2 || maxval > 2 * PI2)
      cdo_warning("%s grid cell center longitudes out of range (min=%.3g/max=%.3g)!", txt, RAD2DEG * minval, RAD2DEG * maxval);
  }

#ifdef _OPENMP
#pragma omp parallel for simd if (nlons > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < nlons; ++i)
  {
    // remove missing values
    if (lons[i] < -PI2) lons[i] = 0;
    if (lons[i] > 2 * PI2) lons[i] = PI2;

    if (lons[i] > PI2) lons[i] -= PI2;
    if (lons[i] < 0.0) lons[i] += PI2;
  }
}

void
check_lat_range(const char *txt, size_t nlats, Varray<double> &lats)
{
  assert(!lats.empty());

  if (txt)
  {
    double minval = 1.e36;
    double maxval = -1.e36;
#ifdef _OPENMP
#pragma omp parallel for if (nlats > cdoMinLoopSize) default(shared) schedule(static) reduction(min : minval) \
    reduction(max : maxval)
#endif
    for (size_t i = 0; i < nlats; ++i)
    {
      minval = std::min(minval, lats[i]);
      maxval = std::max(maxval, lats[i]);
    }
    if (minval < -(PIH + 0.0001) || maxval > (PIH + 0.0001))
      cdo_warning("%s grid cell center latitudes out of range (min=%.3g/max=%.3g)!", txt, RAD2DEG * minval, RAD2DEG * maxval);
  }

#ifdef _OPENMP
#pragma omp parallel for simd if (nlats > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < nlats; ++i)
  {
    if (lats[i] > PIH) lats[i] = PIH;
    if (lats[i] < -PIH) lats[i] = -PIH;
  }
}

static void
grid_check_lat_borders_rad(size_t n, Varray<double> &ybounds)
{
  constexpr double YMAX = PIH;
  constexpr double YLIM = 88 * DEG2RAD;
  auto lrev = (ybounds[0] > ybounds[n - 1]);
  if (lrev)
  {
    if (ybounds[0] > ybounds[1])
    {
      if (ybounds[0] > YLIM) ybounds[0] = YMAX;
      if (ybounds[n - 1] < -YLIM) ybounds[n - 1] = -YMAX;
    }
    else
    {
      if (ybounds[1] > YLIM) ybounds[1] = YMAX;
      if (ybounds[n - 2] < -YLIM) ybounds[n - 2] = -YMAX;
    }
  }
  else
  {
    if (ybounds[0] < ybounds[1])
    {
      if (ybounds[0] < -YLIM) ybounds[0] = -YMAX;
      if (ybounds[n - 1] > YLIM) ybounds[n - 1] = YMAX;
    }
    else
    {
      if (ybounds[1] < -YLIM) ybounds[1] = -YMAX;
      if (ybounds[n - 2] > YLIM) ybounds[n - 2] = YMAX;
    }
  }
}

static void
convert_bounds_reg2d(size_t n, Varray<double> const &boundsIn, Varray<double> &boundsOut)
{
  auto lrev = (boundsIn[0] > boundsIn[2 * n - 1]);
  if (boundsIn[0] > boundsIn[1]) lrev = !lrev;
  if (lrev)
  {
    boundsOut[0] = boundsIn[1];
    for (size_t i = 0; i < n; ++i) boundsOut[i + 1] = boundsIn[2 * i];
  }
  else
  {
    boundsOut[0] = boundsIn[0];
    for (size_t i = 0; i < n; ++i) boundsOut[i + 1] = boundsIn[2 * i + 1];
  }
}

static void
remap_define_reg2d(int gridID, RemapGrid &grid, bool conservMapping, const char *txt)
{
  auto nx = grid.dims[0];
  auto ny = grid.dims[1];
  auto nxp1 = nx + 1;
  auto nyp1 = ny + 1;

  auto nxm = nx;
  if (grid.isCyclic) nxm++;

  if (grid.size != nx * ny) cdo_abort("Internal error, wrong dimensions!");

  grid.centerLonsReg2d.resize(nxm);
  grid.centerLatsReg2d.resize(ny);

  grid.centerLonsReg2d[0] = 0.0;
  grid.centerLatsReg2d[0] = 0.0;
  gridInqXvals(gridID, grid.centerLonsReg2d.data());
  gridInqYvals(gridID, grid.centerLatsReg2d.data());

  static bool doCheck = true;
  if (doCheck)
  {
    doCheck = false;
    check_longitude_range(grid.centerLonsReg2d, "center", cdo_grid_get_units(gridID, CDI_XAXIS, "grid center lon"));
    check_latitude_range(grid.centerLatsReg2d, "center", cdo_grid_get_units(gridID, CDI_YAXIS, "grid center lat"));
  }

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, grid.centerLonsReg2d, "grid reg2d center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, grid.centerLatsReg2d, "grid reg2d center lat");

  if (grid.centerLonsReg2d[nx - 1] < grid.centerLonsReg2d[0])
    for (size_t i = 1; i < nx; ++i)
      if (grid.centerLonsReg2d[i] < grid.centerLonsReg2d[i - 1]) grid.centerLonsReg2d[i] += PI2;

  if (grid.isCyclic) grid.centerLonsReg2d[nx] = grid.centerLonsReg2d[0] + PI2;

  grid.cornerLonsReg2d.resize(nxp1);
  grid.cornerLatsReg2d.resize(nyp1);

  if (gridInqXbounds(gridID, nullptr))
  {
    Varray<double> xbounds(2 * nx);
    gridInqXbounds(gridID, xbounds.data());
    convert_bounds_reg2d(nx, xbounds, grid.cornerLonsReg2d);
    cdo_grid_to_radian(gridID, CDI_XAXIS, grid.cornerLonsReg2d, "grid reg2d corner lon");
  }
  else
  {
    if (conservMapping && nx == 1) cdo_abort("Longitude bounds of %s grid missing!", txt);
    grid_gen_corners(nx, grid.centerLonsReg2d, grid.cornerLonsReg2d);
  }

  if (gridInqYbounds(gridID, nullptr))
  {
    Varray<double> ybounds(2 * ny);
    gridInqYbounds(gridID, ybounds.data());
    convert_bounds_reg2d(ny, ybounds, grid.cornerLatsReg2d);
    cdo_grid_to_radian(gridID, CDI_YAXIS, grid.cornerLatsReg2d, "grid reg2d corner lat");
  }
  else
  {
    if (conservMapping && ny == 1) cdo_abort("Latitude bounds of %s grid missing!", txt);
    grid_gen_corners(ny, grid.centerLatsReg2d, grid.cornerLatsReg2d);
    grid_check_lat_borders_rad(ny + 1, grid.cornerLatsReg2d);
  }
}

static void
init_mask(int gridID, RemapGrid &grid)
{
  auto len = grid.size;
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < len; ++i) grid.mask[i] = 1;

  if (gridInqMask(gridID, nullptr))
  {
    std::vector<int> mask(len);
    gridInqMask(gridID, &mask[0]);
    for (size_t i = 0; i < len; ++i)
      if (mask[i] == 0) grid.mask[i] = 0;
  }
}

static void
remap_define_grid(RemapMethod mapType, int gridID, RemapGrid &grid, char const *txt)
{
  bool destroyGrid = false;
  int gridID_gme = -1;

  auto gridType = gridInqType(grid.gridID);
  auto isHealpixGrid = (grid.type == RemapGridType::HealPix);

  if (isHealpixGrid) { grid.hpParams = cdo::get_healpix_params(gridID); }

  if (gridType != GRID_UNSTRUCTURED && gridType != GRID_CURVILINEAR && !isHealpixGrid)
  {
    if (gridType == GRID_GME)
    {
      gridID_gme = gridToUnstructured(grid.gridID, NeedCorners::Yes);
      grid.nvgp = gridInqSize(gridID_gme);
      gridID = gridDuplicate(gridID_gme);
      gridCompress(gridID);
      grid.useCellCorners = true;
    }
    else if (gridType == GRID_GAUSSIAN_REDUCED || is_healpix_grid(gridID))
    {
      auto needCorners = grid.needCellCorners ? NeedCorners::Yes : NeedCorners::No;
      if (Options::force == false && is_healpix_grid(gridID) && needCorners == NeedCorners::Yes) print_healpix_warning_and_abort();
      destroyGrid = true;
      gridID = gridToUnstructured(grid.gridID, needCorners);
    }
    else if (RemapWriteRemap || grid.type != RemapGridType::Reg2D)
    {
      destroyGrid = true;
      gridID = gridToCurvilinear(grid.gridID, NeedCorners::Yes);
    }
  }

  grid.name = is_healpix_grid(grid.gridID) ? "healpix" : gridNamePtr(gridType);
  grid.size = gridInqSize(gridID);

  grid.dims[0] = isHealpixGrid ? grid.size : gridInqXsize(gridID);
  grid.dims[1] = gridInqYsize(gridID);
  if (gridType != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_UNSTRUCTURED && !isHealpixGrid)
  {
    if (grid.dims[0] == 0) cdo_abort("%s grid without longitude coordinates!", gridNamePtr(gridType));
    if (grid.dims[1] == 0) cdo_abort("%s grid without latitude coordinates!", gridNamePtr(gridType));
  }

  grid.isCyclic = (gridIsCircular(gridID) > 0);

  grid.rank = (gridInqType(gridID) == GRID_UNSTRUCTURED || isHealpixGrid) ? 1 : 2;

  grid.numCorners = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  remap_grid_alloc(mapType, grid);

  // Initialize logical mask
  init_mask(gridID, grid);

  if (!RemapWriteRemap && grid.type == RemapGridType::Reg2D) return;

  if (isHealpixGrid)
  {
    if (RemapWriteRemap)
    {
      if (grid.centerLons.size() == 0 || grid.centerLats.size() == 0)
        cdo_abort("Internal problem - lonlat coordinates not allocated!");

      for (size_t i = 0; i < grid.size; ++i) hp_index_to_lonlat(grid.hpParams, i, &grid.centerLons[i], &grid.centerLats[i]);
    }
  }
  else
  {
    if (!gridHasCoordinates(gridID)) cdo_abort("%s grid cell center coordinates missing!", txt);

    gridInqXvals(gridID, grid.centerLons.data());
    gridInqYvals(gridID, grid.centerLats.data());

    if (grid.needCellCorners)
    {
      if (!gridHasBounds(gridID)) cdo_abort("%s grid cell corner coordinates missing!", txt);

      gridInqXbounds(gridID, grid.cornerLons.data());
      gridInqYbounds(gridID, grid.cornerLats.data());
    }

    if (gridInqType(grid.gridID) == GRID_GME) gridInqMaskGME(gridID_gme, &grid.vgpm[0]);

    // Convert lat/lon units if required
    cdo_grid_to_radian(gridID, CDI_XAXIS, grid.centerLons, "grid center lon");
    cdo_grid_to_radian(gridID, CDI_YAXIS, grid.centerLats, "grid center lat");
    if (grid.numCorners && grid.needCellCorners)
    {
      cdo_grid_to_radian(gridID, CDI_XAXIS, grid.cornerLons, "grid corner lon");
      cdo_grid_to_radian(gridID, CDI_YAXIS, grid.cornerLats, "grid corner lat");
    }

    // Convert longitudes to 0,2pi interval
    check_lon_range(txt, grid.size, grid.centerLons);
    if (grid.numCorners && grid.needCellCorners) check_lon_range(nullptr, grid.numCorners * grid.size, grid.cornerLons);

    // Make sure input latitude range is within the machine values for +/- pi/2
    check_lat_range(txt, grid.size, grid.centerLats);
    if (grid.numCorners && grid.needCellCorners) check_lat_range(nullptr, grid.numCorners * grid.size, grid.cornerLats);
  }

  if (destroyGrid) gridDestroy(gridID);
}

void
remap_grid_alloc(RemapMethod mapType, RemapGrid &grid)
{
  if (grid.nvgp) grid.vgpm.resize(grid.nvgp);

  // only needed for srcGrid and remap_gen_weights
  grid.mask.resize(grid.size);

  if (RemapWriteRemap || (grid.type != RemapGridType::Reg2D && grid.type != RemapGridType::HealPix))
  {
    grid.centerLons.resize(grid.size);
    grid.centerLats.resize(grid.size);
  }

  auto needCellarea = (mapType == RemapMethod::CONSERV);
  if (needCellarea) grid.cellArea.resize(grid.size, 0.0);

  if (RemapGenWeights || mapType == RemapMethod::CONSERV) { grid.cellFrac.resize(grid.size, 0.0); }

  if (grid.needCellCorners && grid.numCorners > 0)
  {
    auto nalloc = grid.numCorners * grid.size;
    grid.cornerLons.resize(nalloc, 0);
    grid.cornerLats.resize(nalloc, 0);
  }
}

void
remap_grid_free(RemapGrid &grid, bool removeMask)
{
  varray_free(grid.vgpm);
  if (removeMask) varray_free(grid.mask);

  varray_free(grid.centerLatsReg2d);
  varray_free(grid.centerLonsReg2d);
  varray_free(grid.cornerLatsReg2d);
  varray_free(grid.cornerLonsReg2d);

  varray_free(grid.centerLats);
  varray_free(grid.centerLons);
  varray_free(grid.cornerLats);
  varray_free(grid.cornerLons);

  varray_free(grid.cellArea);
  varray_free(grid.cellFrac);

  if (grid.tmpgridID != -1)
  {
    gridDestroy(grid.tmpgridID);
    grid.tmpgridID = -1;
  }
}

static void
check_for_convex_cells(const RemapGrid &tgtGrid)
{
  if (tgtGrid.type == RemapGridType::Reg2D) return;

  auto numCorners = tgtGrid.numCorners;
  if (numCorners <= 4) return;

  auto numCells = tgtGrid.size;
  if (numCells > 1000) numCells = 1000;

  for (size_t i = 0; i < numCells; ++i) {}
}

void
remap_search_init(RemapMethod mapType, RemapSearch &search, RemapGrid &srcGrid, RemapGrid &tgtGrid)
{
  search.srcGrid = &srcGrid;
  search.tgtGrid = &tgtGrid;

  auto usePointsearch = (mapType == RemapMethod::KNN);
  if (srcGrid.type != RemapGridType::Reg2D)
  {
    if (srcGrid.type != RemapGridType::HealPix) usePointsearch |= (mapType == RemapMethod::BILINEAR);
    usePointsearch |= (mapType == RemapMethod::BICUBIC);
  }

  auto useCellsearch = (mapType == RemapMethod::CONSERV);

  std::string searchMethodStr;
  cdo::timer timer;

  if (usePointsearch)
  {
    searchMethodStr = "Point search";
    if (srcGrid.doExtrapolate) search.gps.enable_extrapolation();
    grid_pointsearch_create(search.gps, srcGrid);
  }
  else if (useCellsearch)
  {
    searchMethodStr = "Cell search";
    grid_cellsearch_create(search.gcs, srcGrid);
    // check_for_convex_cells(tgtGrid);
  }
  else if (srcGrid.type != RemapGridType::HealPix && srcGrid.type != RemapGridType::Reg2D)
  {
    cdo_abort("remap_search_init: internal error, search not initialized!");
  }

  if (Options::cdoVerbose && searchMethodStr.size()) cdo_print("%s created: %.2f seconds", searchMethodStr, timer.elapsed());
}

void
remap_init_grids(RemapMethod mapType, bool doExtrapolate, int gridID1, RemapGrid &srcGrid, int gridID2, RemapGrid &tgtGrid)
{
  auto reg2d_srcGridID = gridID1;
  auto reg2d_tgtGridID = gridID2;

  if (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC || mapType == RemapMethod::KNN
      || mapType == RemapMethod::CONSERV)
  {
    if (is_reg2d_grid(gridID1)) { srcGrid.type = RemapGridType::Reg2D; }
    else if (is_global_healpix_grid(gridID1) && mapType == RemapMethod::BILINEAR) { srcGrid.type = RemapGridType::HealPix; }
    else if (is_global_healpix_grid(gridID1) && mapType == RemapMethod::KNN && !RemapGenWeights)
    {
      srcGrid.type = RemapGridType::HealPix;
    }
  }

  if (srcGrid.type == RemapGridType::Reg2D)
  {
    if (is_reg2d_grid(gridID2) && mapType == RemapMethod::CONSERV) { tgtGrid.type = RemapGridType::Reg2D; }
    // else srcGrid.type = -1;
  }

  if (!RemapGenWeights && is_reg2d_grid(gridID2) && tgtGrid.type != RemapGridType::Reg2D)
  {
    if (mapType == RemapMethod::KNN) { tgtGrid.type = RemapGridType::Reg2D; }
    if (mapType == RemapMethod::BILINEAR && (srcGrid.type == RemapGridType::Reg2D || srcGrid.type == RemapGridType::HealPix))
    {
      tgtGrid.type = RemapGridType::Reg2D;
    }
  }

  if (!RemapGenWeights && is_healpix_grid(gridID2))
  {
    if (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::KNN) { tgtGrid.type = RemapGridType::HealPix; }
  }

  srcGrid.doExtrapolate = doExtrapolate;

  if (mapType == RemapMethod::CONSERV)
  {
    if (srcGrid.type != RemapGridType::Reg2D && srcGrid.type != RemapGridType::HealPix)
    {
      srcGrid.useCellCorners = true;
      srcGrid.needCellCorners = true;
    }

    if (tgtGrid.type != RemapGridType::Reg2D)
    {
      tgtGrid.useCellCorners = true;
      tgtGrid.needCellCorners = true;
    }
  }

  srcGrid.gridID = gridID1;
  tgtGrid.gridID = gridID2;

  if (gridInqType(gridID1) == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID1))
  {
    auto reference = dereferenceGrid(gridID1);
    if (reference.isValid) { srcGrid.gridID = gridID1 = reference.gridID; }
    if (reference.notFound) { cdo_abort("Reference to source grid not found!"); }
  }

  if (gridInqType(gridID2) == GRID_UNSTRUCTURED && !gridHasCoordinates(gridID2))
  {
    auto reference = dereferenceGrid(gridID2);
    if (reference.isValid) { tgtGrid.gridID = gridID2 = reference.gridID; }
    if (reference.notFound) { cdo_abort("Reference to target grid not found!"); }
  }

  auto sgridID = srcGrid.gridID;
  if (gridInqSize(sgridID) > 1 && gridProjIsSupported(sgridID) && srcGrid.type != RemapGridType::HealPix)
  {
    auto needCorners = srcGrid.needCellCorners ? NeedCorners::Yes : NeedCorners::No;
    if (is_healpix_grid(sgridID))
    {
      if (Options::force == false && needCorners == NeedCorners::Yes) { print_healpix_warning_and_abort(); }
      gridID1 = gridToUnstructured(srcGrid.gridID, needCorners);
    }
    else { gridID1 = gridToCurvilinear(srcGrid.gridID, needCorners); }
    srcGrid.gridID = gridID1;
    srcGrid.tmpgridID = srcGrid.gridID;
  }

  // if (srcGrid.type != RemapGridType::Reg2D)
  remap_define_grid(mapType, gridID1, srcGrid, "Source");
  remap_define_grid(mapType, gridID2, tgtGrid, "Target");

  auto conservMapping = (mapType == RemapMethod::CONSERV);
  if (srcGrid.type == RemapGridType::Reg2D) remap_define_reg2d(reg2d_srcGridID, srcGrid, conservMapping, "source");
  if (tgtGrid.type == RemapGridType::Reg2D) remap_define_reg2d(reg2d_tgtGridID, tgtGrid, conservMapping, "target");
}

/*****************************************************************************/

void
remap_check_area(size_t gridSize, Varray<double> const &cell_area, const char *name)
{
  for (size_t i = 0; i < gridSize; ++i)
  {
    if (cell_area[i] < -0.01) { cdo_print("%s grid area error: %zu %g", name, i, cell_area[i]); }
  }
}

template <typename T>
void
remap_set_mask(Varray<T> const &v, size_t n, size_t numMissVals, double mv, Vmask &mask)
{
  mask.resize(n);

  if (numMissVals)
  {
    T missval = mv;
    if (std::isnan(missval))
    {
#ifdef _OPENMP
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < n; ++i) { mask[i] = fp_is_not_equal(v[i], missval); }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < n; ++i) { mask[i] = !is_equal(v[i], missval); }
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n; ++i) { mask[i] = 1; }
  }
}

// Explicit instantiation
template void remap_set_mask(Varray<float> const &v, size_t gridsize, size_t numMissVals, double mv, Vmask &mask);
template void remap_set_mask(Varray<double> const &v, size_t gridsize, size_t numMissVals, double mv, Vmask &mask);

void
remap_set_mask(Field const &field1, size_t gridsize, size_t numMissVals, double missval, Vmask &imask)
{
  auto func = [&](auto const &v) { remap_set_mask(v, gridsize, numMissVals, missval, imask); };
  field_operation(func, field1);
}
