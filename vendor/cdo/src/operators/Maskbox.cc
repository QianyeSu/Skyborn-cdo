/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Maskbox    masklonlatbox   Mask lon/lat box
      Maskbox    maskindexbox    Mask index box
      Maskbox    maskregion      Mask regions
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "selboxinfo.h"
#include "region.h"

static void
maskbox(Vmask &mask, int gridID, const SelboxInfo &selboxInfo)
{
  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  for (long ilat = 0; ilat < nlat; ilat++)
    for (long ilon = 0; ilon < nlon; ilon++)
      if ((lat1 <= ilat && ilat <= lat2 && ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))))
        mask[nlon * ilat + ilon] = false;
}

void getlonlatparams(int argc_offset, double &xlon1, double &xlon2, double &xlat1, double &xlat2);

static void
maskbox_cell(Vmask &mask, int gridID)
{
  double xlon1 = 0, xlon2 = 0, xlat1 = 0, xlat2 = 0;
  getlonlatparams(0, xlon1, xlon2, xlat1, xlat2);

  auto gridID0 = gridID;
  gridID = generate_full_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  auto gridsize = gridInqSize(gridID);

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

  if (xlon1 > xlon2) cdo_abort("The second longitude have to be greater than the first one!");

  if (xlat1 > xlat2) std::swap(xlat1, xlat2);

  for (size_t i = 0; i < gridsize; ++i)
  {
    mask[i] = true;

    auto xval = xvals[i];
    auto yval = yvals[i];
    if (yval >= xlat1 && yval <= xlat2)
    {
      if (((xval >= xlon1 && xval <= xlon2) || (xval - 360 >= xlon1 && xval - 360 <= xlon2)
           || (xval + 360 >= xlon1 && xval + 360 <= xlon2)))
      {
        mask[i] = false;
      }
    }
  }

  if (gridID0 != gridID) gridDestroy(gridID);
}

static inline bool
is_point_inside(double xval, double yval, double xi, double xj, double yi, double yj)
{
  return (((yval >= yi && yval < yj) || (yval > yj && yval <= yi)) && (xval < ((xj - xi) * (yval - yi) / (yj - yi) + xi)));
}

static bool
point_is_inside(double xval, double yval, size_t n, const double *xcoords, const double *ycoords)
{
  auto c = false;

  for (size_t i = 0, j = n - 1; i < n; j = i++)
  {
    if (is_point_inside(xval, yval, xcoords[i], xcoords[j], ycoords[i], ycoords[j])) c = !c;
  }

  return c;
}

static bool
point_is_inside(double xval, double yval, double xmin, double xmax, const double *xcoords, const double *ycoords, size_t nofcoords)
{
  auto c = false;

  // clang-format off
  if      (xval >= xmin && xval <= xmax)
    c = point_is_inside(xval,         yval, nofcoords, xcoords, ycoords);
  else if (xval > 180.0 && xval - 360.0 >= xmin && xval - 360.0 <= xmax)
    c = point_is_inside(xval - 360.0, yval, nofcoords, xcoords, ycoords);
  else if (xval <   0.0 && xval + 360.0 >= xmin && xval + 360.0 <= xmax)
    c = point_is_inside(xval + 360.0, yval, nofcoords, xcoords, ycoords);
  // clang-format on

  return c;
}

static void
mask_region_regular(Vmask &mask, size_t nlon, size_t nlat, Varray<double> const &xvals, Varray<double> const &yvals,
                    const double *xcoords, const double *ycoords, size_t segmentSize)
{
  auto xmm = varray_min_max(segmentSize, xcoords);
  auto ymm = varray_min_max(segmentSize, ycoords);

  auto gridsize = nlon * nlat;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto ilat = i / nlon;
    auto yval = yvals[ilat];
    if (yval > ymm.min && yval < ymm.max)
    {
      if (point_is_inside(xvals[i - ilat * nlon], yval, xmm.min, xmm.max, xcoords, ycoords, segmentSize)) mask[i] = false;
    }
  }
}

static void
mask_region_cell(Vmask &mask, size_t gridsize, Varray<double> const &xvals, Varray<double> const &yvals, const double *xcoords,
                 const double *ycoords, size_t segmentSize)
{
  auto xmm = varray_min_max(segmentSize, xcoords);
  auto ymm = varray_min_max(segmentSize, ycoords);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto yval = yvals[i];
    if (yval >= ymm.min && yval <= ymm.max)
    {
      if (point_is_inside(xvals[i], yval, xmm.min, xmm.max, xcoords, ycoords, segmentSize)) mask[i] = false;
    }
  }
}

static int
get_gridID(int vlistID1, bool operIndexBox)
{
  std::vector<int> gridsFound;

  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    if (gridInqSize(gridID1) == 1) continue;

    auto gridtype = gridInqType(gridID1);
    auto projtype = gridInqProjType(gridID1);

    auto isReg2dGeoGrid = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR);
    auto projHasGeoCoords = (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_RLL);

    if (isReg2dGeoGrid || projHasGeoCoords || (operIndexBox && (gridtype == GRID_GENERIC || gridtype == GRID_PROJECTION))
        || (!operIndexBox && (gridtype == GRID_UNSTRUCTURED || is_healpix_grid(gridID1))))
    {
      gridsFound.push_back(gridID1);
    }
    else
    {
      if (gridInqSize(gridID1) > 2) cdo_warning("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
  }

  if (gridsFound.size() == 0) cdo_abort("No processable grid found!");
  if (gridsFound.size() > 1) cdo_abort("Too many different grids!");

  auto gridID = gridsFound[0];
  return gridID;
}

static std::vector<bool>
get_processVars(VarList const &varList, int gridID)
{
  auto numVars = varList.numVars();

  std::vector<bool> processVars(numVars, false);

  int varID;
  for (varID = 0; varID < numVars; ++varID)
    if (gridID == varList.vars[varID].gridID) processVars[varID] = true;

  for (varID = 0; varID < numVars; ++varID)
    if (processVars[varID]) break;

  if (varID >= numVars) cdo_abort("No processable variable found!");

  return processVars;
}

class Maskbox : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Maskbox",
    .operators = { { "masklonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude", MaskboxHelp },
                   { "maskindexbox", 0, 0, "index of first and last longitude and index of first and last latitude", MaskboxHelp },
                   { "maskregion", 0, 0, "DCW region or the path to region file", MaskregionHelp },
                   { "maskcircle", 0, 0, "Longitude, latitude of the center and radius of the circle", MaskregionHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Maskbox>();

  int MASKLONLATBOX{}, MASKINDEXBOX{}, MASKREGION{}, MASKCIRCLE{};

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  std::vector<bool> processVars{};
  Vmask mask{};
  size_t gridsize{};

  VarList varList1{};

public:
  void
  init() override
  {
    MASKLONLATBOX = module.get_id("masklonlatbox");
    MASKINDEXBOX = module.get_id("maskindexbox");
    MASKREGION = module.get_id("maskregion");
    MASKCIRCLE = module.get_id("maskcircle");

    auto operatorID = cdo_operator_id();
    auto operIndexBox = (operatorID == MASKINDEXBOX);

    operator_input_arg(cdo_operator_enter(operatorID));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto gridID = get_gridID(vlistID1, operIndexBox);

    processVars = get_processVars(varList1, gridID);

    operator_input_arg(cdo_operator_enter(operatorID));

    gridsize = gridInqSize(gridID);
    mask.resize(gridsize, true);

    auto gridtype = gridInqType(gridID);

    if (operatorID == MASKLONLATBOX)
    {
      if (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || is_healpix_grid(gridID))
        maskbox_cell(mask, gridID);
      else
        maskbox(mask, gridID, gen_lonlat_selbox(0, gridID));
    }
    else if (operatorID == MASKINDEXBOX) { maskbox(mask, gridID, gen_index_selbox(0, gridID)); }
    else if (operatorID == MASKREGION)
    {
      auto nlon = gridInqXsize(gridID);
      auto nlat = gridInqYsize(gridID);
      auto fullGrid = (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || is_healpix_grid(gridID));
      Varray<double> xvals(fullGrid ? gridsize : nlon), yvals(fullGrid ? gridsize : nlat);

      auto gridID0 = gridID;
      if (fullGrid)
      {
        gridID = generate_full_grid(gridID);
        if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");
      }

      gridInqXvals(gridID, xvals.data());
      gridInqYvals(gridID, yvals.data());

      // Convert lat/lon units if required
      cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
      cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

      auto numFiles = cdo_operator_argc();
      for (int i = 0; i < numFiles; ++i)
      {
        Regions regions;
        auto param = cdo_operator_argv(i);
        if (param.starts_with("dcw:"))
          read_regions_from_dcw(param.c_str() + 4, regions);
        else
          read_regions_from_file(param, regions);

        for (size_t k = 0; k < regions.numSegments; ++k)
        {
          auto segmentSize = regions.segmentSize[k];
          if (segmentSize < 3) continue;
          auto offset = regions.segmentOffset[k];
          auto xcoords = &regions.x[offset];
          auto ycoords = &regions.y[offset];
          if (fullGrid)
            mask_region_cell(mask, gridsize, xvals, yvals, xcoords, ycoords, segmentSize);
          else
            mask_region_regular(mask, nlon, nlat, xvals, yvals, xcoords, ycoords, segmentSize);
        }
      }

      if (gridID0 != gridID) gridDestroy(gridID);
    }
    else if (operatorID == MASKCIRCLE)
    {
      /*
      CirclePoint cpoint;

      selcircle_get_parameter(cpoint);
      if (cpoint.radius < 0.0 || cpoint.radius > 180.0) cdo_abort("radius=%g out of bounds (0-180 deg)!", cpoint.radius);

      if (varList1.gridsizeMax() < cpoint.maxpoints) cpoint.maxpoints = varList1.gridsizeMax();
      if (Options::cdoVerbose)
        cdo_print("lon = %g, lat = %g, radius = %gdeg(%gkm)", cpoint.lon, cpoint.lat, cpoint.radius,
                  radiusDegToKm(cpoint.radius));

      cpoint.radius *= DEG2RAD;
      cpoint.lon *= DEG2RAD;
      cpoint.lat *= DEG2RAD;
      */
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;

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

        if (processVars[varID])
        {
          auto const &var = varList1.vars[varID];
          field.init(var);
          cdo_read_field(streamID1, field);

          auto func = [&](auto &v)
          {
            for (size_t i = 0; i < gridsize; ++i)
              if (mask[i]) v[i] = var.missval;
          };
          field_operation(func, field);

          field_num_mv(field);
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, field);
        }
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
