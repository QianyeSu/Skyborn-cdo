/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intgrid    intgridbil      Bilinear grid interpolation
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "interpol.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "matrix_view.h"
#include "remapknn.h"

template <typename T1, typename T2>
static void
thinout(Varray<T1> const &varray1, Varray<T2> &varray2, int gridID1, int gridID2, size_t xinc, size_t yinc)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

  MatrixView<const T1> xfield1(varray1.data(), nlat1, nlon1);
  MatrixView<T2> xfield2(varray2.data(), nlat2, nlon2);

  size_t olat = 0;
  for (size_t ilat = 0; ilat < nlat1; ilat += yinc)
  {
    size_t olon = 0;
    for (size_t ilon = 0; ilon < nlon1; ilon += xinc)
    {
      xfield2[olat][olon] = xfield1[ilat][ilon];
      olon++;
    }
    olat++;
  }
}

static void
thinout(Field const &field1, Field &field2, size_t xinc, size_t yinc)
{
  auto func = [&](auto const &v1, auto &v2, auto grid1, auto grid2) { thinout(v1, v2, grid1, grid2, xinc, yinc); };
  field_operation2(func, field1, field2, field1.grid, field2.grid);
  field_num_mv(field2);
}

static int
gen_thinout_grid(int gridID1, size_t xinc, size_t yinc)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1 / xinc;
  auto nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  auto gridsize2 = nlon2 * nlat2;

  auto gridtype = gridInqType(gridID1);
  auto gridID2 = gridCreate(gridtype, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  grid_copy_names(gridID1, gridID2);

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
  {
    Varray<double> xvals1(nlon1), yvals1(nlat1);
    Varray<double> xvals2(nlon2), yvals2(nlat2);
    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    size_t olat = 0;
    for (size_t ilat = 0; ilat < nlat1; ilat += yinc) yvals2[olat++] = yvals1[ilat];

    size_t olon = 0;
    for (size_t ilon = 0; ilon < nlon1; ilon += xinc) xvals2[olon++] = xvals1[ilon];

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }
  else if (gridtype == GRID_CURVILINEAR)
  {
    Varray<double> xvals1(nlon1 * nlat1), yvals1(nlon1 * nlat1);
    Varray<double> xvals2(nlon2 * nlat2), yvals2(nlon2 * nlat2);
    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    thinout(xvals1, xvals2, gridID1, gridID2, xinc, yinc);
    thinout(yvals1, yvals2, gridID1, gridID2, xinc, yinc);

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

static int
gen_boxavg_grid(int gridID1, size_t xinc, size_t yinc)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1 / xinc;
  auto nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  auto gridsize2 = nlon2 * nlat2;

  auto gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  auto gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
  {
    Varray<double> xvals1(nlon1), yvals1(nlat1);
    Varray<double> xvals2(nlon2), yvals2(nlat2);
    gridInqXvals(gridID1, &xvals1[0]);
    gridInqYvals(gridID1, &yvals1[0]);

    for (size_t i = 0, j = 0; i < nlon1; i += xinc)
    {
      auto i1 = i + (xinc - 1);
      if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
      xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2;
      j++;
    }
    for (size_t i = 0, j = 0; i < nlat1; i += yinc)
    {
      auto i1 = i + (yinc - 1);
      if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
      yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2;
      j++;
    }

    Varray<double> grid2_corner_lon, grid2_corner_lat;
    if (gridHasBounds(gridID1))
    {
      Varray<double> grid1_corner_lon, grid1_corner_lat;
      grid1_corner_lon.resize(2 * nlon1);
      grid1_corner_lat.resize(2 * nlat1);
      grid2_corner_lon.resize(2 * nlon2);
      grid2_corner_lat.resize(2 * nlat2);
      gridInqXbounds(gridID1, &grid1_corner_lon[0]);
      gridInqYbounds(gridID1, &grid1_corner_lat[0]);

      for (size_t i = 0, j = 0; i < nlon1; i += xinc)
      {
        auto i1 = i + (xinc - 1);
        if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
        grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
        grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
        j++;
      }
      for (size_t i = 0, j = 0; i < nlat1; i += yinc)
      {
        auto i1 = i + (yinc - 1);
        if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
        grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
        grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
        j++;
      }
    }

    gridDefXvals(gridID2, &xvals2[0]);
    gridDefYvals(gridID2, &yvals2[0]);

    if (!grid2_corner_lon.empty() && !grid2_corner_lat.empty())
    {
      gridDefNvertex(gridID2, 2);
      gridDefXbounds(gridID2, &grid2_corner_lon[0]);
      gridDefYbounds(gridID2, &grid2_corner_lat[0]);
    }
  }
  else { cdo_abort("Unsupported grid: %s", gridNamePtr(gridtype)); }

  return gridID2;
}

template <typename T1, typename T2>
static void
boxavg(Varray<T1> const &varray1, Varray<T2> &varray2, int gridID1, int gridID2, size_t xinc, size_t yinc)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

  MatrixView<const T1> xfield1(varray1.data(), nlat1, nlon1);
  MatrixView<T2> xfield2(varray2.data(), nlat2, nlon2);

  for (size_t ilat = 0; ilat < nlat2; ilat++)
    for (size_t ilon = 0; ilon < nlon2; ilon++)
    {
      double xsum = 0.0;
      size_t in = 0;
      for (size_t j = 0; j < yinc; ++j)
      {
        auto jj = ilat * yinc + j;
        if (jj >= nlat1) break;
        for (size_t i = 0; i < xinc; ++i)
        {
          auto ii = ilon * xinc + i;
          if (ii >= nlon1) break;
          in++;
          xsum += xfield1[jj][ii];
        }
      }

      xfield2[ilat][ilon] = xsum / in;
    }
}

static void
boxavg(Field const &field1, Field &field2, size_t xinc, size_t yinc)
{
  auto func = [&](auto const &v1, auto &v2, auto grid1, auto grid2) { boxavg(v1, v2, grid1, grid2, xinc, yinc); };
  field_operation2(func, field1, field2, field1.grid, field2.grid);
  field_num_mv(field2);
}

class Intgrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intgrid",
    .operators = { { "intgridbil" }, { "intgriddis" }, { "intgridnn" }, { "intgridknn" }, { "boxavg" }, { "thinout" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Intgrid> registration = RegisterEntry<Intgrid>(module);

private:
  int INTGRID_BIL{}, INTGRID_DIS{}, INTGRID_NN{}, INTGRID_KNN{}, BOXAVG{}, THINOUT{};
  int gridID2 = -1;
  int xinc = 1, yinc = 1;

  int operatorID{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool useFeldMem{};

  VarList varList1{};
  VarList varList2{};

  KnnParams knnParams{};

public:
  void
  init() override
  {
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }

    INTGRID_BIL = module.get_id("intgridbil");
    INTGRID_DIS = module.get_id("intgriddis");
    INTGRID_NN = module.get_id("intgridnn");
    INTGRID_KNN = module.get_id("intgridknn");
    BOXAVG = module.get_id("boxavg");
    THINOUT = module.get_id("thinout");

    operatorID = cdo_operator_id();

    useFeldMem = (operatorID == INTGRID_BIL || operatorID == THINOUT || operatorID == BOXAVG);
    useFeldMem = true;

    if (operatorID == INTGRID_BIL || operatorID == INTGRID_DIS || operatorID == INTGRID_NN)
    {
      operator_input_arg("grid description file or name");
      gridID2 = cdo_define_grid(cdo_operator_argv(0));
      if (operatorID == INTGRID_NN)
      {
        INTGRID_KNN = INTGRID_NN;
        knnParams.k = 1;
        knnParams.kMin = 1;
        knnParams.extrapolate = true;
        knnParams.weighted = WeightingMethod::distanceWeighted;
      }
      if (operatorID == INTGRID_DIS)
      {
        INTGRID_KNN = INTGRID_DIS;
        knnParams.k = 4;
        knnParams.kMin = 1;
        knnParams.extrapolate = true;
        knnParams.weighted = WeightingMethod::distanceWeighted;
      }
    }
    else if (operatorID == INTGRID_KNN)
    {
      auto remapParams = remapknn_get_parameter();
      if (Options::cdoVerbose) remapknn_print_parameter(remapParams);
      if (remapParams.gridString.empty()) cdo_abort("target grid parameter missing!");
      gridID2 = cdo_define_grid(remapParams.gridString);
      knnParams = remapParams.knnParams;
      if (knnParams.kMin == 0) knnParams.kMin = knnParams.k;
      remapknn_verify_parameter(knnParams);
      if (Options::cdoVerbose) print_knn_parameter(knnParams, "KNN parameter: ");
    }
    else if (operatorID == THINOUT || operatorID == BOXAVG)
    {
      operator_input_arg("xinc, yinc");
      operator_check_argc(2);
      xinc = parameter_to_int(cdo_operator_argv(0));
      yinc = parameter_to_int(cdo_operator_argv(1));
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);

      if (operatorID == BOXAVG)
      {
        if (index == 0)
        {
          if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN) cdo_abort("%s data unsupported!", gridNamePtr(gridtype));

          gridID2 = gen_boxavg_grid(gridID1, xinc, yinc);
        }
        else
          cdo_abort("Too many different grids!");
      }
      if (operatorID == THINOUT)
      {
        if (index == 0)
        {
          if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN && gridtype != GRID_CURVILINEAR)
            cdo_abort("%s data unsupported!", gridNamePtr(gridtype));

          gridID2 = gen_thinout_grid(gridID1, xinc, yinc);
        }
        else
          cdo_abort("Too many different grids!");
      }
      else if (operatorID == INTGRID_BIL)
      {
        auto ldistgen = (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2));
        if (!ldistgen && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN)
          cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
      }
      else if (operatorID == INTGRID_KNN)
      {
        auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID1));
        if (!gridProjIsSupported(gridID1) && !hasProjParams && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
            && gridtype != GRID_GME && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED && gridtype != GRID_HEALPIX)
          cdo_abort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
      }

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field1;
    if (!useFeldMem) field1.resize(varList1.gridsizeMax());

    auto gridsize = gridInqSize(gridID2);
    Field field2;
    if (!useFeldMem) field2.resize(gridsize);

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
        if (useFeldMem)
        {
          field1.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field1);

          field2.init(varList2.vars[varID]);
        }
        else
        {
          cdo_read_field(streamID1, field1);

          field1.grid = varList1.vars[varID].gridID;
          field1.missval = varList1.vars[varID].missval;
          field2.grid = gridID2;
          field2.missval = field1.missval;
          field2.numMissVals = 0;
        }

        // clang-format off
        if      (operatorID == INTGRID_BIL)  intgrid_bil(field1, field2);
        else if (operatorID == INTGRID_KNN)  intgrid_knn(knnParams, field1, field2);
        else if (operatorID == BOXAVG)       boxavg(field1, field2, xinc, yinc);
        else if (operatorID == THINOUT)      thinout(field1, field2, xinc, yinc);
        // clang-format on

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
