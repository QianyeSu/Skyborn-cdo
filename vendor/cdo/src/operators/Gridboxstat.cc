/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Gridboxstat    gridboxrange        Gridbox range
      Gridboxstat    gridboxmin          Gridbox minimum
      Gridboxstat    gridboxmax          Gridbox maximum
      Gridboxstat    gridboxsum          Gridbox sum
      Gridboxstat    gridboxmean         Gridbox mean
      Gridboxstat    gridboxavg          Gridbox average
      Gridboxstat    gridboxstd          Gridbox standard deviation
      Gridboxstat    gridboxstd1         Gridbox standard deviation [Normalize by (n-1)]
      Gridboxstat    gridboxvar          Gridbox variance
      Gridboxstat    gridboxvar1         Gridbox variance [Normalize by (n-1)]
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "cdo_options.h"
#include "cdo_omp.h"
#include "field_functions.h"

static void
gen_boxgrid_reg2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);
  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

  {
    Varray<double> xvals1(nlon1), yvals1(nlat1);
    Varray<double> xvals2(nlon2), yvals2(nlat2);
    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    check_longitude_range(xvals1, "center", LonLatUnits::Deg);
    check_latitude_range(yvals1, "center", LonLatUnits::Deg);

    size_t j = 0;
    for (size_t i = 0; i < nlon1; i += xinc)
    {
      auto i1 = i + (xinc - 1);
      if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
      xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2.0;
      j++;
    }

    j = 0;
    for (size_t i = 0; i < nlat1; i += yinc)
    {
      auto i1 = i + (yinc - 1);
      if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
      yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2.0;
      j++;
    }

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }

  if (gridHasBounds(gridID1))
  {
    Varray<double> grid1_corner_lon(2 * nlon1), grid1_corner_lat(2 * nlat1);
    Varray<double> grid2_corner_lon(2 * nlon2), grid2_corner_lat(2 * nlat2);
    gridInqXbounds(gridID1, grid1_corner_lon.data());
    gridInqYbounds(gridID1, grid1_corner_lat.data());

    check_longitude_range(grid1_corner_lon, "corner", LonLatUnits::Deg);
    check_latitude_range(grid1_corner_lat, "corner", LonLatUnits::Deg);

    size_t j = 0;
    for (size_t i = 0; i < nlon1; i += xinc)
    {
      auto i1 = i + (xinc - 1);
      if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
      grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
      grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
      j++;
    }

    j = 0;
    for (size_t i = 0; i < nlat1; i += yinc)
    {
      auto i1 = i + (yinc - 1);
      if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
      grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
      grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
      j++;
    }

    gridDefNvertex(gridID2, 2);
    gridDefXbounds(gridID2, grid2_corner_lon.data());
    gridDefYbounds(gridID2, grid2_corner_lat.data());
  }
}

static void
gen_boxgrid_curv2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);
  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);
  auto gridsize1 = gridInqSize(gridID1);
  auto gridsize2 = gridInqSize(gridID2);

  auto circular = gridIsCircular(gridID1);
  double xvals2_0 = 0.0;

  Varray<double> xvals1(gridsize1), yvals1(gridsize1);
  Varray<double> xvals2(gridsize2), yvals2(gridsize2);
  gridInqXvals(gridID1, xvals1.data());
  gridInqYvals(gridID1, yvals1.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID1, CDI_XAXIS, xvals1, "grid center lon");
  cdo_grid_to_degree(gridID1, CDI_YAXIS, yvals1, "grid center lat");

  check_longitude_range(xvals1, "center", LonLatUnits::Deg);
  check_latitude_range(yvals1, "center", LonLatUnits::Deg);

  // Process grid2 bounds
  double area_norm = xinc * yinc;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t ig = 0; ig < nlat2 * nlon2; ++ig)
  {
    auto y2 = ig / nlon2;
    auto x2 = ig - y2 * nlon2;
    auto g2_add = (y2 * nlon2 + x2);

    for (size_t y1 = y2 * yinc; y1 < yinc * (y2 + 1); y1++)
    {
      auto use_y1 = (y1 >= nlat1) ? nlat1 - 1 : y1;
      for (size_t x1 = x2 * xinc; x1 < xinc * (x2 + 1); x1++)
      {
        auto use_x1 = x1;
        if (x1 >= nlon1)
        {
          if (circular && use_y1 == y1)
            use_y1 -= 1;
          else
            use_x1 = nlon1 - 1;
        }

        auto g1_add = (use_y1 * nlon1) + use_x1;
        auto xval1 = xvals1[g1_add];
        auto yval1 = yvals1[g1_add];

        constexpr double xrangeLim = 359.99;
        if (y1 == y2 * yinc && x1 == x2 * xinc)
        {
          xvals2_0 = xval1;
          xvals2[g2_add] = xval1 / area_norm;
          yvals2[g2_add] = yval1 / area_norm;
        }
        else if (std::fabs(xval1 - xvals2_0) > xrangeLim)
        {
          if ((xval1 - xvals2_0) > xrangeLim)
            xvals2[g2_add] += (xval1 - 360.0) / area_norm;
          else if ((xval1 - xvals2_0) < -xrangeLim)
            xvals2[g2_add] += (xval1 + 360.0) / area_norm;
          yvals2[g2_add] += yval1 / area_norm;
        }
        else
        {
          xvals2[g2_add] += xval1 / area_norm;
          yvals2[g2_add] += yval1 / area_norm;
        }
      }  // x1
    }    // y1

    //  while ( xvals2[g2_add] >  180. ) xvals2[g2_add] -= 360.;
    //  while ( xvals2[g2_add] < -180. ) xvals2[g2_add] += 360.;
  }

  gridDefXvals(gridID2, xvals2.data());
  gridDefYvals(gridID2, yvals2.data());
}

static int
gen_boxgrid(int gridID1, size_t xinc, size_t yinc)
{
  if (xinc < 1 || yinc < 1) cdo_abort("xinc and yinc must not be smaller than 1!");

  int gridID2 = -1;
  auto gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR || gridtype == GRID_GENERIC)
  {
    auto nlon1 = gridInqXsize(gridID1);
    auto nlat1 = gridInqYsize(gridID1);
    if (xinc > nlon1 || yinc > nlat1) cdo_abort("xinc and/or yinc exceeds gridsize!");

    auto nlon2 = nlon1 / xinc;
    auto nlat2 = nlat1 / yinc;
    if (nlon1 % xinc) nlon2++;
    if (nlat1 % yinc) nlat2++;
    auto gridsize2 = nlon2 * nlat2;

    gridID2 = gridCreate(gridtype, gridsize2);
    gridDefXsize(gridID2, nlon2);
    gridDefYsize(gridID2, nlat2);
  }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT) { gen_boxgrid_reg2D(gridID1, xinc, yinc, gridID2); }
  else if (gridtype == GRID_GENERIC) {}
  else if (gridtype == GRID_CURVILINEAR) { gen_boxgrid_curv2D(gridID1, xinc, yinc, gridID2); }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

static void
gridbox_stat(Field const &field1, Field &field2, size_t xinc, size_t yinc, int statfunc)
{
  auto useWeight = !field1.weightv.empty();

  auto boxsize = xinc * yinc;
  FieldVector fields(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    fields[i].resize(boxsize);
    if (useWeight) fields[i].weightv.resize(boxsize);
    fields[i].missval = field1.missval;
  }

  auto gridID1 = field1.grid;
  auto gridID2 = field2.grid;

  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t ig = 0; ig < nlat2 * nlon2; ++ig)
  {
    auto ompthID = cdo_omp_get_thread_num();

    auto &field = fields[ompthID];
    field.numMissVals = 0;

    auto ilat = ig / nlon2;
    auto ilon = ig - ilat * nlon2;

    size_t isize = 0;
    for (size_t j = 0; j < yinc; ++j)
    {
      auto jj = ilat * yinc + j;
      if (jj >= nlat1) break;
      for (size_t i = 0; i < xinc; ++i)
      {
        auto ii = ilon * xinc + i;
        auto index = jj * nlon1 + ii;
        if (ii >= nlon1) break;
        field.vec_d[isize] = field1[index];
        if (fp_is_equal(field.vec_d[isize], field.missval)) field.numMissVals++;
        if (useWeight) field.weightv[isize] = field1.weightv[index];
        isize++;
      }
    }

    field.size = isize;
    auto func = [&](auto &v) { v[ig] = field_function(field, statfunc); };
    field_operation(func, field2);
  }

  field_num_mv(field2);
}

class Gridboxstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Gridboxstat",
    .operators = { { "gridboxrange", FieldFunc_Range, 0, GridboxstatHelp },
                   { "gridboxmin", FieldFunc_Min, 0, GridboxstatHelp },
                   { "gridboxmax", FieldFunc_Max, 0, GridboxstatHelp },
                   { "gridboxsum", FieldFunc_Sum, 0, GridboxstatHelp },
                   { "gridboxmean", FieldFunc_Meanw, 1, GridboxstatHelp },
                   { "gridboxavg", FieldFunc_Avgw, 1, GridboxstatHelp },
                   { "gridboxstd", FieldFunc_Stdw, 1, GridboxstatHelp },
                   { "gridboxstd1", FieldFunc_Std1w, 1, GridboxstatHelp },
                   { "gridboxvar", FieldFunc_Varw, 1, GridboxstatHelp },
                   { "gridboxvar1", FieldFunc_Var1w, 1, GridboxstatHelp },
                   { "gridboxskew", FieldFunc_Skew, 0, GridboxstatHelp },
                   { "gridboxkurt", FieldFunc_Kurt, 0, GridboxstatHelp },
                   { "gridboxmedian", FieldFunc_Median, 0, GridboxstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Gridboxstat> registration = RegisterEntry<Gridboxstat>();

  int lastgrid = -1;
  bool wstatus = false;

  int operfunc{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  VarList varList2{};

  int xinc{};
  int yinc{};

  bool needWeights{};
  int gridID1{ CDI_UNDEFID };

public:
  void
  init() override
  {
    operator_input_arg("xinc, yinc");
    operator_check_argc(2);
    xinc = parameter_to_int(cdo_operator_argv(0));
    yinc = parameter_to_int(cdo_operator_argv(1));

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);
    needWeights = (cdo_operator_f2(operatorID) != 0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids > 1) cdo_abort("Too many different grids!");

    gridID1 = vlistGrid(vlistID1, 0);

    auto gridID2 = gen_boxgrid(gridID1, xinc, yinc);
    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    Field field1, field2;
    auto gridsize1 = gridInqSize(gridID1);
    if (needWeights) field1.weightv.resize(gridsize1);

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
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        if (needWeights && field1.grid != lastgrid)
        {
          lastgrid = field1.grid;
          wstatus = gridcell_weights(field1.grid, field1.weightv);
        }
        if (wstatus != 0 && tsID == 0 && levelID == 0)
          cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!",
                      varList1.vars[varID].name);

        field2.init(varList2.vars[varID]);
        gridbox_stat(field1, field2, xinc, yinc, operfunc);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field2);
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
