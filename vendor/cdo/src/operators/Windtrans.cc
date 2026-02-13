/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  This file is a extension of CDO.

  Author: M. Koutek, 2012, KNMI (NL)

*/

/*
   This module contains the following operators:

      Windtrans       uvDestag   : destagger U and/or V wind (in place)
                      rotuvNorth : rotate grid-relative wind(u,v) to North_pole-relative
*/

#include <cdi.h>

#include <mpim_grid.h>
#include "grid_define.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"

#define MAXARG 10

// NOTE: Variable with codes (3[3,4],105,10) will get from CDO typically a name:
// "10u", resp. "10v" Mostly u & v at level type 105 is not staggered, but test it to be sure..

#define CheckVarIsU(varID, varname, code)   \
  {                                         \
    VarIsU = 0;                             \
    if (lvar)                               \
    {                                       \
      if (varname == chvars[0]) VarIsU = 1; \
    }                                       \
    else if (code == chcodes[0])            \
      VarIsU = 1;                           \
  }
#define CheckVarIsV(varID, varname, code)   \
  {                                         \
    VarIsV = 0;                             \
    if (lvar)                               \
    {                                       \
      if (varname == chvars[1]) VarIsV = 1; \
    }                                       \
    else if (code == chcodes[1])            \
      VarIsV = 1;                           \
  }
#define CheckUVisStaggered(varID1, varID2, zaxisID1, zaxisID2) \
  {                                                            \
    VarsUVareStaggered = 1;                                    \
    auto gridID1tmp = varList1.vars[varID1].gridID;            \
    auto gridID2tmp = varList1.vars[varID2].gridID;            \
    if (zaxisID1 != zaxisID2)                                  \
      VarsUVareStaggered = 0;                                  \
    else if (gridID1tmp == gridID2tmp)                         \
      VarsUVareStaggered = 0;                                  \
  }

static void
destaggerUorV(double *fu, double *fuOut, int klev, long nlat, long nlon, long UorV, long offset)
{
  long lat0, lon0;
  double u0, u1;
  double u_dstg;
  long next;
  /* This does de-staggering (-0.5; -0.5) */

  long idx = offset;
  if (UorV == 0)
  {
    next = -1;  // U-wind
    lat0 = 0;
    lon0 = 1;
  }
  else
  {
    next = -nlon;  // V-wind: length of the row (2d array)
    lat0 = 1;
    lon0 = 0;
  }

  if (cdoDebugExt >= 20)
    cdo_print("destaggerUorV(): (nlon=%d, nlat=%d);               [lat0=%ld, lon0=%ld, next=%ld];  (default order destaggering)",
              nlon, nlat, lat0, lon0, next);

  for (int lev = 0; lev < klev; lev++)
  {
    if (lat0)  // the first row we let as it is
      for (int lon = 0; lon < nlon; lon++)
      {
        u0 = fu[idx];
        fuOut[idx] = u0;
        idx++;
      }
    for (int lat = lat0; lat < nlat; lat++)
    {
      if (lon0)  // the first column we let as it is
      {
        u0 = fu[idx];
        fuOut[idx] = u0;
        idx++;
      }
      for (int lon = lon0; lon < nlon; lon++)
      {
        u0 = fu[idx];
        u1 = fu[idx + next];
        u_dstg = 0.5 * (u0 + u1);
        fuOut[idx] = u_dstg;
        idx++;
      }
    }
  }
}

static void
destaggerUorV_positiveOrder(double *fu, double *fuOut, int klev, long nlat, long nlon, int UorV, long offset)
{
  long latE, lonE;
  double u0, u1;
  double u_dstg;
  long next;
  /* This does de-staggering (+0.5; +0.5) */

  long idx = offset;
  if (UorV == 0)  // U-wind: length of the row (2d array)
  {
    next = nlon;
    latE = nlat - 1;
    lonE = nlon;
  }
  else  // V-wind:
  {
    next = 1;
    latE = nlat;
    lonE = nlon - 1;
  }

  if (cdoDebugExt >= 20)
    cdo_print("destaggerUorV(): (nlon=%ld, nlat=%ld);               [latE=%d, lonE=%d, next=%d];  (positive order destaggering)",
              nlon, nlat, latE, lonE, next);

  for (int lev = 0; lev < klev; lev++)
  {
    for (long lat = 0; lat < latE; lat++)
    {
      for (long lon = 0; lon < lonE; lon++)
      {
        u0 = fu[idx];
        u1 = fu[idx + next];
        u_dstg = 0.5 * (u0 + u1);
        fuOut[idx] = u_dstg;
        idx++;
      }
      if (lonE < nlon)  // the last column we let as it is
      {
        u0 = fu[idx];
        fuOut[idx] = u0;
        idx++;
      }
    }

    if (latE < nlat)  // the last row we let as it is
      for (long lon = 0; lon < nlon; lon++)
      {
        u0 = fu[idx];
        fuOut[idx] = u0;
        idx++;
      }
  }
}

static void *
DestaggerUV()
{
  int varID1 = CDI_UNDEFID, varID2 = CDI_UNDEFID;
  int zaxisID1 = CDI_UNDEFID, zaxisID2 = CDI_UNDEFID;
  int varID1stg = CDI_UNDEFID, varID2stg = CDI_UNDEFID;
  size_t gridsize;
  int chcodes[MAXARG];
  const char *chvars[MAXARG];
  double destagGridOffsets[MAXARG];
  int gridID1 = CDI_UNDEFID, gridID2 = CDI_UNDEFID;
  int gridID0 = CDI_UNDEFID;
  int gridID;
  bool lcopy = false;
  int UorV;
  size_t nlon = 0, nlat = 0;
  Varray<double> ivar, ovar;
  double dxU = 0, dyU = 0, dxV = 0, dyV = 0;

  chcodes[0] = 0;
  chcodes[1] = 0;
  chvars[0] = 0;
  chvars[1] = 0;

  // Note: Already initialized by the caller! Don't call again:
  // cdo_initialize(process);

  operator_input_arg("Pair of u and v in the staggered system:\n\
    Usage: uvDestag,u,v -or- uvDestag,33,34 -or- uvDestag,u,v,-0.5,-0.5 -or- uvDestag,33,34,-0.5,-0.5\n \
    Destaggered grid offsets <,-/+0.5,-/+0.5> are optional.\n           \
    If file contains grid with temperature (name='t' or code=11) then grid_temp will be used for destaggered wind.");

  Debug(cdoDebugExt, "UVDESTAG (destaggering) requested)..");

  int nch = cdo_operator_argc();
  if (nch < 2) cdo_abort("Number of input arguments < 2; At least 2 arguments needed: uvDestag,33,34<,-0.5,-0.5> optional");
  if (nch >= MAXARG) cdo_abort("Number of input arguments >= %d", MAXARG);

  bool lvar = false;
  if (std::isdigit(cdo_operator_argv(0)[0]))
  {
    lvar = false;  // We have a list of codes
    for (int i = 0; i < 2; ++i) chcodes[i] = parameter_to_int(cdo_operator_argv(i));
  }
  else
  {
    lvar = true;  // We have a list of variables
    for (int i = 0; i < 2; ++i) chvars[i] = cdo_operator_argv(i).c_str();
  }

  destagGridOffsets[0] = -0.5;
  destagGridOffsets[1] = -0.5;

  if (nch > 2)
  {
    for (int i = 2; i < (2 + 2); ++i) destagGridOffsets[i - 2] = parameter_to_double(cdo_operator_argv(i));
  }

  Debug(cdoDebugExt, "destagGridOffsets = (%01.1f,%01.1f)", destagGridOffsets[0], destagGridOffsets[1]);

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  VarList varList1(vlistID1);
  VarList varList2(vlistID2);

  // Find the first occurance of staggered U and V variables (for example codes 33,34). We assume that one (grib-)file
  // contains only 1 sort of horizontal grids and at most 2 staggered grids.
  /*
    [Hirlam/hip_work] >cdo sinfo LAMH_D11_201302150000_00000_GB_only_UV
    File format: GRIB
    -1 : Institut Source   Param       Ttype   Dtype  Gridsize Num  Levels Num
     1 : KNMI     unknown  11.1        instant  P11    399300   1      19   1
     2 : KNMI     unknown  33.1        instant  P13    399300   1       6   2  ** non-staggered U  (leveltype=105)
     3 : KNMI     unknown  34.1        instant  P12    399300   1       6   2  ** non-staggered V  (leveltype=105)
     4 : KNMI     unknown  33.1        instant  P14    399300   2      60   3  ** STAGGERED U  (leveltype=109)
     5 : KNMI     unknown  34.1        instant  P14    399300   3      60   3  ** STAGGERED V  (leveltype=109)
     6 : KNMI     unknown  11.1        instant  P11    399300   1      60   3
     7 : KNMI     unknown  11.1        instant  P10    399300   1       1   4
     8 : KNMI     unknown  11.1        instant  P11    399300   1      11   5
     9 : KNMI     unknown  33.1        instant  P14    399300   2      11   5  ** STAGGERED U  (leveltype=100)
    10 : KNMI     unknown  34.1        instant  P14    399300   3      11   5  ** STAGGERED V  (leveltype=100)
   Horizontal grids :
     1 : lonlat   > size      : dim = 399300  nlon = 726  nlat = 550
                    rlon      : first = -30.2  last = 42.3  inc = 0.1  degrees
                    rlat      : first = -30.8  last = 24.1  inc = 0.1  degrees
                    northpole : lon = -195  lat = 30
     2 : lonlat       > size  : dim = 399300  nlon = 726  nlat = 550
                    rlon      : first = -30.15  last = 42.35  inc = 0.1  degrees
                    rlat      : first = -30.8  last = 24.1  inc = 0.1  degrees
                    northpole : lon = -195  lat = 30
     3 : lonlat       > size  : dim = 399300  nlon = 726  nlat = 550
                    rlon      : first = -30.2  last = 42.3  inc = 0.1  degrees
                    rlat      : first = -30.75  last = 24.15  inc = 0.1  degrees
                    northpole : lon = -195  lat = 30
   Vertical grids :
     1 : height                 m : 0 2 801 802 803 804 805 901 902 903 904 905 951 952 953 954 955 998 999
     2 : height                 m : 10 801 802 803 804 805
     3 : hybrid             level : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
                                   19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43
                                   44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
     4 : meansea            level : 0
     5 : pressure              Pa : 5000 10000 20000 25000 30000 40000 50000 70000 85000 92500 100000
  */
  int VarIsU, VarIsV;
  int VarsUVareStaggered;

  // Search for staggered u and v wind:
  int numVars = varList1.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var = varList1.vars[varID];

    int pnum, pcat, pdis;
    cdiDecodeParam(var.param, &pnum, &pcat, &pdis);
    auto code = pnum;
    auto zaxisID = var.zaxisID;
    auto ltype = zaxis_to_ltype(zaxisID);
    auto nlevs = zaxisInqSize(zaxisID);
    auto gridIDx = var.gridID;
    if (cdoDebugExt >= 20)
      cdo_print("Var.id [%4d] with grib code:%3d and has name: %6s; level "
                "type: %3d; number of levels: %3d; gridID: %d; zaxisID: %d",
                varID, code, var.name, ltype, nlevs, gridIDx, zaxisID);

    CheckVarIsU(varID, var.name, code);
    CheckVarIsV(varID, var.name, code);
    if (VarIsU)
    {
      varID1 = varID;
      zaxisID1 = zaxisID;
    }
    else if (VarIsV)
    {
      varID2 = varID;
      zaxisID2 = zaxisID;
    }

    if ((varID1 != CDI_UNDEFID) && (varID2 != CDI_UNDEFID))
    {
      CheckUVisStaggered(varID1, varID2, zaxisID1, zaxisID2);
      if (VarsUVareStaggered)
      {
        gridID1 = varList1.vars[varID1].gridID;
        gridID2 = varList1.vars[varID2].gridID;
        if (cdoDebugExt)
          cdo_print("Found STAGGERED U & V: varID1=%d (gridID1=%d), varID2=%d (gridID2=%d)", varID1, gridID2, varID2, gridID1);
        varID1stg = varID1;
        varID2stg = varID2;
        if (gridID0 != CDI_UNDEFID)
        {
          vlistChangeVarGrid(vlistID2, varID1stg, gridID0);  // set the variable onto the non-staggered grid
          vlistChangeVarGrid(vlistID2, varID2stg, gridID0);  // set the variable onto the non-staggered grid
        }
        // Allow a next level-type UV-pair to be found;
        // NOTE: There may be separate CDO staggerd variables for (33/34; 109; *) and (33/34; 100; *)
        varID1 = varID2 = CDI_UNDEFID;
      }
    }
    // search for a reference (non-staggered) grid
    // We take temperature for example as the new (horizontal) grid for
    // de-staggered uv If there will be no temperature field we will define grid the grid
    if (lvar)  // We have a list of variables
    {
      if (var.name == "t") gridID0 = var.gridID;
    }
    else
    {
      if (code == 11) gridID0 = var.gridID;
    }
  }  // end of for ( varID = 0; varID < nvars; ..

  if (gridID0 >= 0) Debug(cdoDebugExt, "Found DESTAGGERED grid for U, V: gridID0=%d", gridID0);

  if (varID1stg == -1 && varID2stg == -1)
  {
    cdo_print("NOTE: We did not find any staggered U,V wind components. Performing file-copy.");
    lcopy = true;
    gridID0 = gridID1;
  }

  auto numGrids = vlistNumGrids(vlistID1);

  double xincU = 0.0, yincU = 0.0, xincV = 0.0, yincV = 0.0;

  if (!lcopy)
  {
    if (gridInqXsize(gridID1) != gridInqXsize(gridID2)) cdo_abort("Xsize(gridID1) != Xsize(gridID2)");
    if (gridInqYsize(gridID1) != gridInqYsize(gridID2)) cdo_abort("Ysize(gridID1) != Ysize(gridID2)");

    nlon = gridInqXsize(gridID1);
    nlat = gridInqYsize(gridID1);

    xincU = gridInqXinc(gridID1);
    yincU = gridInqXinc(gridID1);
    xincV = gridInqXinc(gridID2);
    yincV = gridInqXinc(gridID2);

    auto gridtype1 = gridInqType(gridID1);
    auto gridtype2 = gridInqType(gridID2);

    if (!(gridtype1 == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL && gridtype2 == GRID_PROJECTION
          && gridInqProjType(gridID2) == CDI_PROJ_RLL))
    {
      cdo_print("U - wind: Grid nr. %d is gridtype: %d (%s)", gridID1, gridtype1, gridNamePtr(gridtype1));
      cdo_print("V - wind: Grid nr. %d is gridtype: %d (%s)", gridID2, gridtype2, gridNamePtr(gridtype2));
      cdo_abort("Destaggering supports only grid type = 'lonlat' (GRID_LONLAT).");
    }

    // define output grid
    if (gridID0 == -1)
    {
      if (cdoDebugExt)
        cdo_print("Calling define_destagered_grid( destagGridOffsets = (%01.1f,%01.1f) )", destagGridOffsets[0],
                  destagGridOffsets[1]);
      gridID0 = cdo_define_destagered_grid(gridID1, gridID2, destagGridOffsets);
    }

    if (gridID0 == -1) cdo_abort("Cannot define DESTAGGERED grid for U, V.");

    if (cdoDebugExt >= 10) cdo_print_griddes(gridID0, 1);

    double xfirst_R = gridInqXval(gridID0, 0);  // reference grid for non-staggered fields
                                                // (default: search for temperature; otherwise: create a new grid)
    double yfirst_R = gridInqYval(gridID0, 0);
    double xfirst_U = gridInqXval(gridID1, 0);  // grid of u-wind
    double yfirst_U = gridInqYval(gridID1, 0);
    double xfirst_V = gridInqXval(gridID2, 0);  // grid of v-wind
    double yfirst_V = gridInqYval(gridID2, 0);

    dxU = -xfirst_U + xfirst_R;
    dyU = -yfirst_U + yfirst_R;
    dxV = -xfirst_V + xfirst_R;
    dyV = -yfirst_V + yfirst_R;

    if (cdoDebugExt)
    {
      cdo_print("Grid info: (xfirst_R = %3.2f; yfirst_R = %3.2f); (xfirst_U "
                "= %3.2f; yfirst_U = %3.2f); (xfirst_V = %3.2f; yfirst_V = %3.2f);",
                xfirst_R, yfirst_R, xfirst_U, yfirst_U, xfirst_V, yfirst_V);
      cdo_print("Grid info: (dxU; dyU) = (%3.2f; %3.2f); (dxV; dyV) = (%3.2f; %3.2f) ", dxU, dyU, dxV, dyV);
      cdo_print("Grid info: nlon=%zu, nlat=%zu ", nlon, nlat);
    }

    if (cdoDebugExt)
    {
      auto order = (dxU < 0) ? "default" : "positive";
      cdo_print("About to perform destaggering (U-wind): (%3.2f; %3.2f) - %s order ", dxU, dyV, order);
    }

    if (cdoDebugExt)
    {
      auto order = (dyV < 0) ? "default" : "positive";
      cdo_print("About to perform destaggering (V-wind): (%3.2f; %3.2f) - %s order ", dxU, dyV, order);
    }

    for (int index = 0; index < numGrids; ++index)
    {
      gridID = vlistGrid(vlistID1, index);
      if (cdoDebugExt >= 10)
        cdo_print("Grid nr. %d is gridtype: %d (%s)", index, gridInqType(gridID), gridNamePtr(gridInqType(gridID)));
    }

    if (gridID0 == -1)
    {
      Debug(cdoDebugExt, "Last trial to find a reference grid for destaggered wind.");
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];

        gridID = var.gridID;
        Debug(cdoDebugExt, "Var.id %d has grid nr:%d", varID, gridID);
        if ((varID != varID1stg) && (varID != varID2stg))
        {
          // this will the new (horizontal) grid for de-staggered
          gridID0 = var.gridID;
        }
      }
    }
    if (gridID0 == -1) cdo_abort("Referencial horizontal grid not found!");

    if (gridInqSize(gridID2) != gridInqSize(gridID1))
      cdo_abort("gridSize of U-wind != gridSize of V-wind!  This should not happen!");

    if (cdoDebugExt)
      cdo_print("Allocating memory for maximum gridsize (for input) = %ld [%4.3f MB]", varList1.gridsizeMax(),
                varList1.gridsizeMax() * sizeof(double) / (1024.0 * 1024));
    ivar.resize(varList1.gridsizeMax());  // storage for other fields than

    gridsize = gridInqSize(gridID1);  // actual size of U-wind should be same as V-wind
    if (cdoDebugExt)
      cdo_print("Allocating memory for gridsize (destaggered output)= %ld; nlon=%zu, nlat=%zu", gridsize, nlon, nlat);
    ovar.resize(gridsize);
  }  // end of  if (!lcopy)

  auto streamID2 = cdo_open_write(1);

  if (varID1stg != CDI_UNDEFID && varID2stg != CDI_UNDEFID)
  {
    vlistChangeVarGrid(vlistID2, varID1stg, gridID0);  // set the variable onto the non-staggered grid
    vlistChangeVarGrid(vlistID2, varID2stg, gridID0);  // set the variable onto the non-staggered grid
  }

  cdo_def_vlist(streamID2, vlistID2);  // from this point the stream is using a different vlistID !!!!!

  vlistID2 = cdo_stream_inq_vlist(streamID2);  // refresh it

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
    if (numFields == 0) break;

    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsID);

    if (!lcopy && cdoDebugExt)
    {
      cdo_print("Processing timestep: %d", tsID);
      cdo_print("Starting destaggering. Total fields to be processed: %05d", numFields);
    }

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);

      auto const &var = varList1.vars[varID];

      int pnum, pcat, pdis;
      cdiDecodeParam(var.param, &pnum, &pcat, &pdis);
      auto code = pnum;
      auto zaxisID = var.zaxisID;
      auto ltype = zaxis_to_ltype(zaxisID);
      auto level = (int) zaxisInqLevel(zaxisID, levelID);

      if (!lcopy)
      {
        // NOTE: There may be separate CDO staggerd variables for (33/34; 109; *) and (33/34; 100; *)
        //       You cannot use this way of U/V variable detection!
        // VarIsU = (varID == varID1stg); ** DON'T USE
        // VarIsV = (varID == varID2stg); ** DON'T USE
        CheckVarIsU(varID, var.name, code);
        CheckVarIsV(varID, var.name, code);

        UorV = -1;  // -1: not U, neither V
        if (VarIsU)
          UorV = 0;  // 0: U-wind; 1: V-wind
        else if (VarIsV)
          UorV = 1;

        if (UorV >= 0)
        {
          gridID = var.gridID;
          if (UorV == 0)
          {
            if (gridID == gridID1)  // Has this variable the staggered grid?
            {
              if (cdoDebugExt >= 10)
                cdo_print("Destaggering U-wind record: %05d (timestep:%d); Var.id [%4d]; (code=%3d; "
                          "ltype=%3d; level=%4d; levelID=%3d); GridID %d => %d  *** <===",
                          fieldID, tsID, varID, code, ltype, level, levelID, var.gridID, varList2.vars[varID].gridID);
            }
            else
              UorV = -1;  // this U is not staggered, just copy the record..
          }
          if (UorV == 1)
          {
            if (gridID == gridID2)  // Has this variable the staggered grid?
            {
              if (cdoDebugExt >= 10)
                cdo_print("Destaggering V-wind record: %05d (timestep:%d); Var.id [%4d]; (code=%3d; "
                          "ltype=%3d; level=%4d; levelID=%3d); GridID %d => %d  *** <===",
                          fieldID, tsID, varID, code, ltype, level, levelID, var.gridID, varList2.vars[varID].gridID);
            }
            else
              UorV = -1;  // this V is not staggered, just copy the record..
          }
        }  // end of: if (UorV>=0)

        if (UorV >= 0)  // re-check again since it could mean that current record with U or V is not staggered
        {
          size_t numMissVals;
          cdo_read_field(streamID1, ivar.data(), &numMissVals);

          // void destaggerUorV(double *fu, double *fuOut, int klev, int nlat, int nlon, int UorV, long int offset);
          // We handle one level at the time; klev=1;offset=0.
          if ((dxU < 0.0) && (dyV < 0.0))
            // UorV = 0: U-wind; 1: V-wind
            destaggerUorV(ivar.data(), ovar.data(), 1, nlat, nlon, UorV, 0);
          else if ((dyU > 0.0) && (dxV > 0.0))
            // UorV = 0: U-wind; 1: V-wind
            destaggerUorV_positiveOrder(ivar.data(), ovar.data(), 1, nlat, nlon, UorV, 0);
          else
            cdo_abort("Unsupported destaggering grid offset: (dxU; dyU) = (%3.2f; %3.2f); (dxV; dyV) = (%3.2f; %3.2f) "
                      "where: xincU=%3.2f, yincU=%3.2f, xincV=%3.2f, yincV=%3.2f",
                      dxU, dyU, dxV, dyV, xincU, yincU, xincV, yincV);

          // Typical Hirlam LAMH_D11 situation with destaggering on (-0.5,-0.5)
          // cdo uvDestag: Grid info: (xfirst_R = -30.20; yfirst_R = -30.80);
          //                          (xfirst_U = -30.15; yfirst_U = -30.80);
          //                          (xfirst_V = -30.20; yfirst_V = -30.75);
          // cdo uvDestag: Grid info: (dxU; dyU) = (-0.05; 0.00); (dxV; dyV) = (0.00; -0.05)

          // Less typical would be to choose destaggering on (+0.5,+0.5)
          // cdo uvDestag: Grid info: (xfirst_R = -30.15; yfirst_R = -30.75);
          //                          (xfirst_U = -30.15; yfirst_U = -30.80);
          //                          (xfirst_V = -30.20; yfirst_V = -30.75);
          // cdo uvDestag: Grid info: (dxU; dyU) = (0.00; 0.05); (dxV; dyV) = (0.05; 0.00)
          if (cdoDebugExt >= 20) cdo_print("Setting GRID id from: %d => to: %d", var.gridID, varList2.vars[varID].gridID);

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, ovar.data(), numMissVals);
        }
        else
        {  // copy the record to the output unchanged...
          cdo_def_field(streamID2, varID, levelID);
          if (cdoDebugExt >= 20)
            cdo_print("Stream-copy data record:    %05d (timestep:%d); "
                      "Var.id [%4d]; (code=%3d; ltype=%3d; level=%4d; levelID=%3d)",
                      fieldID, tsID, varID, code, ltype, level, levelID);
          cdo_copy_field(streamID1, streamID2);
        }
      }  // end of: if (!lcopy)
      else
      {  // copy the record to the output unchanged...
        cdo_def_field(streamID2, varID, levelID);
        if (cdoDebugExt >= 20)
          cdo_print("Stream-copy data record:    %05d (timestep:%d); "
                    "Var.id [%4d]; (code=%3d; ltype=%3d; level=%4d; levelID=%3d)",
                    fieldID, tsID, varID, code, ltype, level, levelID);
        cdo_copy_field(streamID1, streamID2);
      }

    }  // end of for ( fieldID = ...

    tsID++;
  }  // end of while ( (numFields ...

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);

  return 0;
}

static void
rot_uv_north(int gridID, double *us, double *vs)
{
  // Function transforms grid-relative UV into north-pole relative UV in the (spherical) lat-lon coordinate space.
  // Input:  u,v : grid relative (u,v); NOT staggered
  //
  // This function needs that the gridpoint coordinates have been transformed
  // from modelspace into LatLon space using function gridToCurvilinear().
  //
  // xvals[], yvals[] contains the Lon-Lat coordinates.

  if (gridInqType(gridID) != GRID_CURVILINEAR)
    cdo_abort("%s(gridname=%s) grid must be GRID_CURVILINEAR!", __func__, gridNamePtr(gridInqType(gridID)));
  // this should never happen

  static double *rotationMatrixArray = nullptr;
  double lon_pnt0, lat_pnt0;
  double lon_pntEast = 0, lat_pntEast = 0;
  double lon_pntNorth, lat_pntNorth;
  double dLatEast = 0, dLonEast = 0;
  double dLatNorth, dLonNorth;
  double xpntEastSph, ypntEastSph, zpntEastSph = 0;
  double xpntNorthSph, ypntNorthSph, zpntNorthSph;
  double xpntNorthSphRot, ypntNorthSphRot, zpntNorthSphRot;
  double xpnt0Sph, ypnt0Sph, zpnt0Sph;
  double xnormSph, ynormSph, znormSph;
  double xncross, yncross, zncross;
  double vecAngle;
  size_t i, j;

  // The following "correction for the grid-step direction" from funtion
  // project_uv_latlon() cannot be used here! This cannot be done in the
  // geographic lon/lan (on the globe).
  // int signLon=( (xvals[1] - xvals[0]) < 0 )?-1:1;
  // int signLat=( (yvals[1] - yvals[0]) < 0 )?-1:1;
  /*
      Tested for scanning mode 64 LAMH_D11:
      -------------------------------------
      iScansNegatively = 0;
      jScansPositively = 1;
      jPointsAreConsecutive = 0;
      #-READ ONLY- alternativeRowScanning = 0;
      jDirectionIncrementInDegrees = 0.1;
      iDirectionIncrementInDegrees = 0.1;

      !!! When scanning mode 00 the lons & lats arrays: xvals[idx] & yvals[idx]; are upside-down !!!

        iScansNegatively = 0;
        jScansPositively = 0;
        jPointsAreConsecutive = 0;
        #-READ ONLY- alternativeRowScanning = 0;
        jDirectionIncrementInDegrees = 0.1;  // this should be -0.1 ??
        iDirectionIncrementInDegrees = 0.1;
  */

  auto nx = gridInqXsize(gridID);
  auto ny = gridInqYsize(gridID);

  Varray<double> xvals(nx * ny), yvals(nx * ny);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  int scanningMode = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_SCANNINGMODE, &scanningMode);
  auto jScansPositively = (scanningMode == 64);
  if (scanningMode == 64)
  {
    if (cdoDebugExt > 1) cdo_print("NOTICE: Processing data with scanning mode(%d); gridID=%d", scanningMode, gridID);
  }
  else if (scanningMode == 0)
  {
    if (cdoDebugExt > 1) cdo_print("NOTICE: Processing data with scanning mode(%d); gridID=%d", scanningMode, gridID);
  }
  else
  {
    cdo_abort("\n***\n***\n WARNING! Unsupported data scanning mode(%d); gridID=%d; For this operation we support only: 64,00 \n"
              "RESULT will be probably incorrect!\n***\n***",
              scanningMode, gridID);
  }

  if (cdoDebugExt)
    cdo_print("%s(gridname=%s) .. processing grid with UV [nx*ny] (%zu * %zu)", __func__, gridNamePtr(gridInqType(gridID)), nx, ny);

  if (gridInqSize(gridID) != (nx * ny)) cdo_abort("Incorrect gridsize (%zu) != nx*ny (%zu * %zu)", gridInqSize(gridID), nx, ny);
  // this should never happen

#define OPTrotuvNorth 1  // ACTIVATE SPEED - OPTIMIZATION

#define radians(aDeg) (DEG2RAD * aDeg)
#define NormVector(vec0, vec1, vec2)                                    \
  {                                                                     \
    double vecLen = std::sqrt(vec0 * vec0 + vec1 * vec1 + vec2 * vec2); \
    vec0 = vec0 / vecLen;                                               \
    vec1 = vec1 / vecLen;                                               \
    vec2 = vec2 / vecLen;                                               \
  }

#define CrossProd(vecx0, vecx1, vecx2, vecy0, vecy1, vecy2, vecz0, vecz1, vecz2) \
  {                                                                              \
    vecz0 = vecx1 * vecy2 - vecy1 * vecx2;                                       \
    vecz1 = vecx2 * vecy0 - vecy2 * vecx0;                                       \
    vecz2 = vecx0 * vecy1 - vecy0 * vecx1;                                       \
  }

  if (rotationMatrixArray == nullptr)
  {
    if (cdoDebugExt > 0) cdo_print("About to compute rotationMatrixArray for the whole grid [%zu x %zu]", nx, ny);

    rotationMatrixArray = new double[4 * nx * ny];
    for (j = 0; j < ny; ++j)
      for (i = 0; i < nx; ++i)
      {
        auto idx = j * nx + i;
        lon_pnt0 = xvals[idx];
        lat_pnt0 = yvals[idx];
#ifndef OPTrotuvNorth
        // For speed - optimization not used and not needed. Kept only for clarity.
        if ((i + 1) < nx)
        {
          lon_pntEast = xvals[idx + 1];  // longitude of grid point towards east
          lat_pntEast = yvals[idx + 1];  //  latitude of grid point towards east
        }
        else
        {
          lon_pntEast = lon_pnt0 + 1.0 * (lon_pnt0 - xvals[idx - 1]);  // at the grid border define extended gridpoint
          lat_pntEast = lat_pnt0 + 1.0 * (lat_pnt0 - yvals[idx - 1]);
        }
#endif
        if ((j + 1) < ny)
        {
          lon_pntNorth = xvals[idx + nx];  // longitude of grid point towards north
          lat_pntNorth = yvals[idx + nx];  //  latitude of grid point towards north
        }
        else
        {
          lon_pntNorth = lon_pnt0 + 1.0 * (lon_pnt0 - xvals[idx - nx]);  // at the grid border define extended gridpoint
          lat_pntNorth = lat_pnt0 + 1.0 * (lat_pnt0 - yvals[idx - nx]);
        }

        // (lon_pntNorth, lat_pntNorth)
        //     ^
        //     |       (lon_pntCenter, lat_pntCenter)   center of the
        //     cell-diagonal
        //     |
        // (lon_pnt0,lat_pnt0) ----> (lon_pntEast,lat_pntEast)

        // lon_pntCenter = 0.5*(lon_pntNorth + lon_pntEast);
        // lat_pntCenter = 0.5*(lat_pntNorth + lat_pntEast);
        // lon_pnt0 -= lon_pntCenter; lon_pntEast -= lon_pntCenter;
        // lon_pntNorth -= lon_pntCenter;  lat_pnt0 -= lat_pntCenter;
        // lat_pntEast -= lat_pntCenter; lat_pntNorth -= lat_pntCenter;

        // This is the local coordinate system of a grid cell where we have
        // (u,v) at location (xpnt0,ypnt0).

        // The local coordinate system is now centered around
        // (lon_pnt0,lat_pnt0) The vector towards north pole (UP-VECTOR) at
        // this location will be (0,1,0) The tangent plane at this location
        // is XY wil a normal (0, 0, 1)

        // Nummerical approach using projection onto a unit sphere
        lon_pnt0 = radians(lon_pnt0);
        lat_pnt0 = radians(lat_pnt0);
        xpnt0Sph = std::cos(lat_pnt0) * std::cos(lon_pnt0);
        ypnt0Sph = std::cos(lat_pnt0) * std::sin(lon_pnt0);  // # Get // [lon_pnt0,lat_pnt0] on the unit sphere.
        zpnt0Sph = std::sin(lat_pnt0);                       // # Only XY plane is needed.
        dLonNorth = radians(lon_pntNorth);
        dLatNorth = radians(lat_pntNorth);
        xpntNorthSph = std::cos(dLatNorth) * std::cos(dLonNorth);
        ypntNorthSph = std::cos(dLatNorth) * std::sin(dLonNorth);  // # Get [dLonNorth,dLatNorth] on the unit sphere.
        zpntNorthSph = std::sin(dLatNorth);                        // # Only XY plane is needed.
        xpntNorthSph -= xpnt0Sph, ypntNorthSph -= ypnt0Sph;
        zpntNorthSph -= zpnt0Sph;
        NormVector(xpntNorthSph, ypntNorthSph, zpntNorthSph);  // vecy

#ifndef OPTrotuvNorth
        // For speed - optimization not used and not needed. Kept only for clarity.
        dLonEast = radians(lon_pntEast);
        dLatEast = radians(lat_pntEast);
        xpntEastSph = std::cos(dLatEast) * std::cos(dLonEast);
        ypntEastSph = std::cos(dLatEast) * std::sin(dLonEast);  // # Get [dLonEast,dLatEast] on the unit sphere.
        zpntEastSph = std::sin(dLatEast);                       // # Only XY plane is needed.
        xpntEastSph -= xpnt0Sph;
        ypntEastSph -= ypnt0Sph;
        zpntEastSph -= zpnt0Sph;                            // make vectors from points
        NormVector(xpntEastSph, ypntEastSph, zpntEastSph);  // vecx
        // vecz = CrossProd(vecx,vecy)
        CrossProd(xpntEastSph, ypntEastSph, zpntEastSph, xpntNorthSph, ypntNorthSph, zpntNorthSph, xnormSph, ynormSph,
                  znormSph);  // vec z
#else
        xnormSph = xpnt0Sph;
        ynormSph = ypnt0Sph;
        znormSph = zpnt0Sph;
        NormVector(xnormSph, ynormSph, znormSph);  // vec z ... normal vector to the sphere at pnt0 (lon_pnt0,lat_pnt0)
#endif  // #ifdef OPTrotuvNorth

        // # vecUP = (0.0,0.0,1.0) .. up-vector in a global coordinate system
        // #   ^^ up-vector & vector-towards-north-pole & shere-normal-vector
        //  belong to ONE plane #      that crosses north pole and point pnt0
        //(lon_pnt0,lat_pnt0) # Project vecUP onto plane XY, where
        //  plane-normal is vecz # vecnProjXY = vecUP - D*vecz;   D=
        //  a*x1+b*y1+c*z1;  vecz=(a,b,c); vecUP = (x1,y1,z1)=(0,0,1) #
        //  D= vecz[2]*1; # vecyRot = NormVector( (0.0 - vecz[2]*vecz[0],0.0
        //- vecz[2]*vecz[1], 1.0  - vecz[2]*vecz[2]) )

        // double Dist =  xnormSph * 0.0 +  ynormSph * 0.0 + znormSph * 1.0;
        // // Left out for optimization
        xpntNorthSphRot = -znormSph * xnormSph;       // xpntNorthSphRot = 0.0 - Dist*xnormSph;
        ypntNorthSphRot = -znormSph * ynormSph;       // ypntNorthSphRot = 0.0 - Dist*ynormSph;
        zpntNorthSphRot = 1.0 - znormSph * znormSph;  // zpntNorthSphRot = 1.0 - Dist*znormSph;
        NormVector(xpntNorthSphRot, ypntNorthSphRot, zpntNorthSphRot);

        // This would create in 3D the rotated Easting vector; but we don't
        // need it in this routine. Left out to optimize the computation..
        // CrossProd( xpntNorthSphRot, ypntNorthSphRot, zpntNorthSphRot,
        // xnormSph, ynormSph, znormSph,
        //            xpntEastSph,  ypntEastSphRot,  zpntEastSphRot );
        //            //vecxRot = CrossProd(vecy,vecz)

        vecAngle = std::acos((xpntNorthSph * xpntNorthSphRot + ypntNorthSph * ypntNorthSphRot + zpntNorthSph * zpntNorthSphRot));
        // Determine the sign of the angle
        CrossProd(xpntNorthSphRot, ypntNorthSphRot, zpntNorthSphRot, xpntNorthSph, ypntNorthSph, zpntNorthSph, xncross, yncross,
                  zncross);
        if ((xncross * xnormSph + yncross * ynormSph + zncross * znormSph) > 0.0)  // dotProduct
          vecAngle *= -1.0;

        if (!jScansPositively) vecAngle += radians(180.0);  // this is needed in scanning mode 00

        xpntNorthSph = std::sin(vecAngle);  // Rotate the point/vector (0,1) around Z-axis with vecAngle
        ypntNorthSph = std::cos(vecAngle);
        xpntEastSph = ypntNorthSph;  // Rotate the same point/vector around Z-axis with 90 degrees
        ypntEastSph = -xpntNorthSph;

        // zpntNorthSph = 0; zpntEastSph = 0;  // not needed in 2D

        // 1) Build the rotation matrix and put the axes-base vectors into the matrix
        auto VJaa = xpntEastSph;
        auto VJab = xpntNorthSph;
        auto VJba = ypntEastSph;
        auto VJbb = ypntNorthSph;

        auto idx4 = 4 * idx;
        // Caching the rotation matrix for later usage ..
        rotationMatrixArray[idx4++] = VJaa;
        rotationMatrixArray[idx4++] = VJab;
        rotationMatrixArray[idx4++] = VJba;
        rotationMatrixArray[idx4++] = VJbb;

        if (cdoDebugExt >= 20)
          if (((i < 3) && (j < 3)) || ((i > (nx - 3)) && (j > (ny - 3))))
          {
            cdo_print("grid point [%zu,%zu] with latlon[%3.6f,%3.6f]; (lon_pntNorth, lat_pntNorth) = [%3.6f,%3.6f]; "
                      "dLonNorth=%3.6f; dLatNorth=%3.6f (Northing grid relative) ",
                      i, j, lon_pnt0, lat_pnt0, lon_pntNorth, lat_pntNorth, RAD2DEG * dLonNorth, RAD2DEG * dLatNorth);
            cdo_print("grid point [%zu,%zu] with latlon[%3.6f,%3.6f]; (lon_pntEast,lat_pntEast    )= [%3.6f,%3.6f]; "
                      "dLonEast =%3.6f; dLatEast =%3.6f (Easting grid relative ) ",
                      i, j, lon_pnt0, lat_pnt0, lon_pntEast, lat_pntEast, RAD2DEG * dLonEast, RAD2DEG * dLatEast);
            // cdo_print("(xpntNorthSph, ypntNorthSph)= [%3.6f,%3.6f];
            // (xpntEastSph,ypntEastSph) = [%3.6f,%3.6f];",
            //         xpntNorthSph, ypntNorthSph, xpntEastSph,ypntEastSph
            //         );
            // vecAngle = RAD2DEG * std::acos( (xpntEastSph*xpntNorthSph +
            // ypntEastSph*ypntNorthSph + zpntEastSph*zpntNorthSph) );
            // vecAngle = RAD2DEG * std::acos( (xpntEastSph*xpntNorthSph +
            // ypntEastSph*ypntNorthSph) );
            cdo_print("(xpntNorthSph, ypntNorthSph, zpntNorthSph)= [%3.6f,%3.6f,%3.6f]; (xpntEastSph,ypntEastSph, "
                      "zpntEastSph) = [%3.6f,%3.6f,%3.6f]; vecAngle= %3.6f",
                      xpntNorthSph, ypntNorthSph, zpntNorthSph, xpntEastSph, ypntEastSph, zpntEastSph, vecAngle);
            cdo_print("rotation matrix for grid point [%zu,%zu] with latlon[%3.6f,%3.6f]: (VJaa, VJab, VJba, VJbb) = "
                      "(%3.6f,%3.6f,%3.6f,%3.6f)",
                      i, j, lon_pnt0, lat_pnt0, VJaa, VJab, VJba, VJbb);
          }

      }  // end of for ( i = 0; i < nx; i++ )
  }  // end of if (rotationMatrixArray== nullptr)

  // Take the rotation matrix from the cache
  for (j = 0; j < ny; ++j)
    for (i = 0; i < nx; ++i)
    {
      auto idx = (j * nx + i);
      auto idx4 = 4 * idx;
      auto VJaa = rotationMatrixArray[idx4++];
      auto VJab = rotationMatrixArray[idx4++];
      auto VJba = rotationMatrixArray[idx4++];
      auto VJbb = rotationMatrixArray[idx4++];

      // 2) Transform the UV vector with jacobian matrix
      auto u = us[idx];
      auto v = vs[idx];
      // u = 6.0;  v = 0.0; // test: 6 m/s along the easting direction of the grid
      auto magnitude = hypot(u, v);  // old vector magnitude in the model space
      //(uu) =   (VJaa VJab) * ( u )
      //(vv)     (VJba VJbb)   ( v )
      auto uu = VJaa * u + VJab * v;
      auto vv = VJba * u + VJbb * v;
      //(uu) =   (VJaa VJab VJac) * ( u )
      //(vv)     (VJba VJbb VJbc)   ( v )
      //(ww)     (VJba VJbb VJcc)   ( w )

      // 3) Apply scaling of the vector so that the vector keeps the original length (model space)
      auto newMagnitude = hypot(uu, vv);
      us[idx] = uu * magnitude / newMagnitude;
      vs[idx] = vv * magnitude / newMagnitude;
    }

  Debug(cdoDebugExt >= 20, "%s(gridname=%s) finished.", __func__, gridNamePtr(gridInqType(gridID)));
}

static void
rot_uv_back_mode64(int gridID, double *us, double *vs)
{
  // This function is partially based on rot_uv_back() of the CDO rotuv
  // operator. This routine expects the data to be in scanning-mode 64. This
  // routine gives comparable (not numerically same) results as rot_uv_north().
  // rot_uv_back_mode64() is significantly slower than rot_uv_north().

  int scanningMode = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_SCANNINGMODE, &scanningMode);
  if (scanningMode == 64)
  {
    if (cdoDebugExt > 1) cdo_print("NOTICE: Processing data with scanning mode(%d); gridID=%d", scanningMode, gridID);
  }
  else
  {
    cdo_print("\n***\n***\n WARNING! Unsupported data scanning mode(%d); gridID=%d; For this operation we support only: 64\n"
              "RESULT will be probably incorrect!\n***\n***",
              scanningMode, gridID);
  }

  double xpole = 0, ypole = 0, angle = 0;
  if (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL)
    gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);

  Varray<double> xvals(nlon), yvals(nlat);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, 1, &xpole, "xpole");
  cdo_grid_to_degree(gridID, CDI_YAXIS, 1, &ypole, "ypole");
  cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

  double u, v;
  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++)
    {
      auto i = ilat * nlon + ilon;

      auto xval = lamrot_to_lam(yvals[ilat], xvals[ilon], ypole, xpole, angle);
      auto yval = phirot_to_phi(yvals[ilat], xvals[ilon], ypole, angle);

      usvs_to_uv(us[i], vs[i], yval, xval, ypole, xpole, &u, &v);

      us[i] = u;
      vs[i] = v;
    }
}

static void
project_uv_latlon(int gridID, double *us, double *vs)
{
  // The function project_uv_latlon() is WARPING the UV-vector from
  // grid-relative to North-pole-relative definition in LAT-LON projection,
  // flatten into 2D space. The resulting vectors do NOT exit in a SPHERICAL
  // LAT-LON space BUT in a strongly WARPED (2D) LAT-LON projection space.
  // Imagine what happens to the grid-points and grid-cells close to the poles.
  // And what affect it has on the "grid-cell relative UV-vectors" in
  // NON-rectangular grid-cells. This routine is VERY fast as it uses JACOBIANs.
  // Depending on shape distorsion of each cell we see the "flow" vectors to
  // adapt. The resulting UV-vectors can by directly plot with (2D) plotting
  // tools in LAT-LON. Even (2D) streamlines can be used on this way transformed
  // vectors. This function expectes scaning-mode 64.

  if (gridInqType(gridID) != GRID_CURVILINEAR)
    cdo_abort("%s(gridname=%s) transformation grid must be GRID_CURVILINEAR!", __func__, gridNamePtr(gridInqType(gridID)));
  // this should never happen

  auto nx = gridInqXsize(gridID);
  auto ny = gridInqYsize(gridID);

  Varray<double> xvals(nx * ny), yvals(nx * ny);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  int signLon = ((xvals[1] - xvals[0]) < 0) ? -1 : 1;
  int signLat = ((yvals[1] - yvals[0]) < 0) ? -1 : 1;

  if (cdoDebugExt)
    cdo_print("%s(gridname=%s) .. processing grid with UV [nx*ny] (%zu * %zu)", __func__, gridNamePtr(gridInqType(gridID)), nx, ny);

  if (gridInqSize(gridID) != (nx * ny)) cdo_abort("Incorrect gridsize (%zu) != nx*ny (%zu * %zu)", gridInqSize(gridID), nx, ny);
  // this should never happen

  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      auto idx = j * nx + i;
      auto xpnt0 = xvals[idx];
      auto ypnt0 = yvals[idx];
      double xpntEast, ypntEast;
      double xpntNorth, ypntNorth;
      if ((i + 1) < nx)
      {
        xpntEast = xvals[idx + 1];  // longitude of grid point towards east
        ypntEast = yvals[idx + 1];  //  latitude of grid point towards east
      }
      else
      {
        xpntEast = xpnt0 + 0.01 * (xpnt0 - xvals[idx - 1]);  // at the grid border define extended gridpoint
        ypntEast = ypnt0 + 0.01 * (ypnt0 - yvals[idx - 1]);
      }
      if ((j + 1) < ny)
      {
        xpntNorth = xvals[idx + nx];  // longitude of grid point towards north
        ypntNorth = yvals[idx + nx];  //  latitude of grid point towards north
      }
      else
      {
        xpntNorth = xpnt0 + 0.01 * (xpnt0 - xvals[idx - nx]);  // at the grid border define extended gridpoint
        ypntNorth = ypnt0 + 0.01 * (ypnt0 - yvals[idx - nx]);
      }
      /*
          This is the local coordinate system of a grid cell where we have
         (u,v) at location (xpnt0,ypnt0). Gridpoint coordinates have been
         transformed from modelspace into LatLon space using function gridToCurvilinear().

          (xpntNorth, ypntNorth)
              ^
              |
              |
          (xpnt0,ypnt0) ----> (xpntEast,ypntEast)

          modelXLon=modelX+deltaX;    xpntEast
          modelYLon=modelY;           ypntEast
          modelXLat=modelX;           xpntNorth
          modelYLat=modelY+deltaY;    ypntNorth

          distLon = hypot(modelXLon-lo, modelYLon-la);
          distLat = hypot(modelXLat-lo, modelYLat-la);
      */
      ///  Basically transforms grid-relative UV into north pole relative UV
      ///  (~ lat-lon space) u,v : grid relative (u,v) not staggered get(u); get(v)
      // 1) Build the jacobian matrix
      auto distLon = hypot(xpntEast - xpnt0, ypntEast - ypnt0);
      auto distLat = hypot(xpntNorth - xpnt0, ypntNorth - ypnt0);
      auto VJaa = signLon * (xpntEast - xpnt0) / distLon;
      auto VJab = signLon * (xpntNorth - xpnt0) / distLat;
      auto VJba = signLat * (ypntEast - ypnt0) / distLon;
      auto VJbb = signLat * (ypntNorth - ypnt0) / distLat;
      if (cdoDebugExt >= 20)
        if (((i < 3) && (j < 3)) || ((i > (nx - 3)) && (j > (ny - 3))))
          cdo_print("Jacobian for grid point [%zu,%zu] with latlon[%3.6f,%3.6f]: (VJaa, VJab, VJba, VJbb) = "
                    "(%3.6f,%3.6f,%3.6f,%3.6f)",
                    i, j, xpnt0, ypnt0, VJaa, VJab, VJba, VJbb);
      // 2) Transform the UV vector with jacobian matrix
      auto u = us[idx];
      auto v = vs[idx];
      // u = 6.0;  v = 0.0; // test: 6 m/s along the easting direction of the grid
      auto magnitude = hypot(u, v);  // old vector magnitude in the model space
      auto uu = VJaa * u + VJab * v;
      auto vv = VJba * u + VJbb * v;
      // 3) Apply scaling of the vector so that the vector keeps the original length (model space)
      auto newMagnitude = hypot(uu, vv);
      us[idx] = uu * magnitude / newMagnitude;
      vs[idx] = vv * magnitude / newMagnitude;
    }

  Debug(cdoDebugExt >= 20, "%s(gridname=%s) finished.", __func__, gridNamePtr(gridInqType(gridID)));
}

class Windtrans : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Windtrans",
    .operators = { { "uvDestag", WindtransHelp },
                   { "rotuvN", WindtransHelp },
                   { "rotuvNorth", WindtransHelp },
                   { "projuvLatLon", WindtransHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Windtrans> registration = RegisterEntry<Windtrans>();
  int UVDESTAG{}, ROTUVNORTH{}, ROTUVN{}, PROJUVLATLON{};
  int operatorID{};

  void *
  TransformUV()
  {
    int varID1, varID2, nlevel1, nlevel2;
    size_t gridsize = 0;
    int code, gridID;
    int ltype, level, nlevs, zaxisID;
    int pnum, pcat, pdis;
    int chcodes[MAXARG];
    const char *chvars[MAXARG];
    int gridIDcurvl = -1;
    int gridIDlastused = -1;

    // Note: Already initialized by the caller! Don't call again: cdo_initialize(process);

    operator_input_arg("Pairs of u and v in the rotated system;\n usage: rotuvNorth,u,v  -or- rotuvNorth,33,34");

    int nch = cdo_operator_argc();
    if (nch != 2) cdo_abort("Number of input arguments != 2");
    if (nch >= MAXARG) cdo_abort("Number of input arguments >= %d", MAXARG);

    bool lvar = false;  // We have a list of codes
    int len = (int) cdo_operator_argv(0).size();
    int ix = (cdo_operator_argv(0)[0] == '-') ? 1 : 0;
    for (int i = ix; i < len; ++i)
      if (!std::isdigit(cdo_operator_argv(0)[i]))
      {
        lvar = true;  // We have a list of variables
        break;
      }

    if (lvar)
    {
      for (int i = 0; i < nch; ++i) chvars[i] = cdo_operator_argv(i).c_str();
    }
    else
    {
      for (int i = 0; i < nch; ++i) chcodes[i] = parameter_to_int(cdo_operator_argv(i));
    }

    auto streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    VarList varList1(vlistID1);
    VarList varList2(vlistID2);

    auto numVars = varList1.numVars();
    auto maxFields = varList1.maxFields();

    std::vector<FieldInfo> fieldInfoList(maxFields);

    std::vector<std::vector<size_t>> varnumMissVals(numVars);
    Varray3D<double> vardata(numVars);

    // U & V are NOT grid relative
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(vlistID2, varID, CDI_KEY_UVRELATIVETOGRID, 0);

    bool lfound[MAXARG];
    for (int i = 0; i < nch; ++i) lfound[i] = false;

    if (lvar)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int i = 0; i < nch; ++i)
          if (varList1.vars[varID].name == chvars[i]) lfound[i] = true;
      }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Variable %s not found!", chvars[i]);
    }
    else
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int i = 0; i < nch; ++i)
          if (varList2.vars[varID].code == chcodes[i]) lfound[i] = true;
      }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Code %d not found!", chcodes[i]);
    }

    int VarIsU, VarIsV;

    // NOTE: Variable with codes (3[3,4],105,10) will get from CDO typically a name: "10u", resp. "10v"

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];

      cdiDecodeParam(var.param, &pnum, &pcat, &pdis);
      code = pnum;
      zaxisID = var.zaxisID;
      ltype = zaxis_to_ltype(zaxisID);
      nlevs = zaxisInqSize(zaxisID);

      gridID = var.gridID;
      if (cdoDebugExt >= 20)
        cdo_print(
            "Var.id [%4d] with grib code:%3d and has name: %6s; level type: %3d; number of levels: %3d; gridID: %d; zaxisID: %d",
            varID, code, var.name, ltype, nlevs, gridID, zaxisID);

      if (!(gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL))
        cdo_abort("Only rotated lon/lat grids supported!");

      CheckVarIsU(varID, var.name, code);
      CheckVarIsV(varID, var.name, code);
      if (VarIsU || VarIsV)
      {
        gridsize = gridInqSize(gridID);
        if (cdoDebugExt)
          cdo_print("Allocating memory for variableID %4d (code=%3d): gridsize(%zu)*nlevels(%d) = %zu [%4.3f MB]", varID,
                    varList2.vars[varID].code, gridsize, nlevs, gridsize * nlevs,
                    gridsize * nlevs * sizeof(double) / (1024.0 * 1024));
        varnumMissVals[varID].resize(nlevs);
        vardata[varID].resize(nlevs);
        for (int levelID = 0; levelID < nlevs; ++levelID) vardata[varID][levelID].resize(gridsize);
      }
    }

    Debug(cdoDebugExt, "Neccessary memory has been allocated.");

    auto taxisID1 = vlistInqTaxis(vlistID1);
    auto taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);          // from this point the stream is using a different vlistID !!!!!
    vlistID2 = cdo_stream_inq_vlist(streamID2);  // refresh it

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      if (cdoDebugExt) cdo_print("About to read U & V data to memory. Other data will be stream-copied to the output file.");

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        code = varList1.vars[varID].code;
        zaxisID = varList1.vars[varID].zaxisID;
        ltype = zaxis_to_ltype(zaxisID);
        level = zaxisInqLevel(zaxisID, levelID);

        if (vardata[varID].empty())
        {                                     // This means that it is not eighter U neither V.
          fieldInfoList[fieldID].varID = -1;  // We will NOT record/store this field in memory
          fieldInfoList[fieldID].levelID = -1;
          // We will stream-copy this data
          cdo_def_field(streamID2, varID, levelID);
          // if ( cdoDebugExt>10 ) cdo_print("Copying data record.. %05d (timestep:%05d)", fieldID, tsID);
          if (cdoDebugExt >= 20)
            cdo_print("Stream-copy data record:    %05d (timestep:%d); Var.id [%4d]; (code=%3d; ltype=%3d; level=%4d; "
                      "levelID=%3d)",
                      fieldID, tsID, varID, code, ltype, level, levelID);
          cdo_copy_field(streamID1, streamID2);  // cannot do this ! We have to set the flag uvGridRelative = 0
        }
        else
        {
          fieldInfoList[fieldID].varID = varID;
          fieldInfoList[fieldID].levelID = levelID;
          if (cdoDebugExt >= 10)
            cdo_print("Memmory-read data record:   %05d (timestep:%d); Var.id [%4d]; (code=%3d; ltype=%3d; level=%4d; "
                      "levelID=%3d)",
                      fieldID, tsID, varID, code, ltype, level, levelID);
          cdo_read_field(streamID1, vardata[varID][levelID].data(), &varnumMissVals[varID][levelID]);
          if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
        }
      }  // end of for ( fieldID = 0; fieldID < numFields; ..

      Debug(cdoDebugExt, "All neccessary U & V data are in memory. About to transform the windvectors...");

      int code1, zaxisID1, ltype1;
      int code2, zaxisID2, ltype2;
      Debug(cdoDebugExt, "Looping over %d variables to look for U-wind..", numVars);

      // find u-variables:
      for (varID1 = 0; varID1 < numVars; ++varID1)
        if (vardata[varID1].empty())
        {
          // if ( cdoDebugExt ) cdo_print("Checking U-wind: vardata[%d]== nullptr ",varID1);
        }
        else  // This means that it is U or V.
        {
          code1 = varList2.vars[varID1].code;
          zaxisID1 = varList2.vars[varID1].zaxisID;
          ltype1 = zaxis_to_ltype(zaxisID1);
          nlevel1 = zaxisInqSize(zaxisID1);
          CheckVarIsV(varID1, varList2.vars[varID1].name, code1);
          if (VarIsV) continue;
          if (cdoDebugExt >= 20)
            cdo_print("Checking U-wind: Var.id [%4d] with grib code:%3d; name: %6s; level type: %3d; number of levels: %3d; "
                      "zaxisID: %d",
                      varID1, code1, varList2.vars[varID1].name, ltype1, nlevel1, zaxisID1);
          CheckVarIsU(varID1, varList2.vars[varID1].name, code1);
          if (!VarIsU) continue;
          if (cdoDebugExt >= 10)
            cdo_print("** FOUND U-wind; Var.id [%4d] with grib code:%3d; name: %6s; level type: %3d; number of levels: %3d; "
                      "zaxisID: %d",
                      varID1, code1, varList2.vars[varID1].name, ltype1, nlevel1, zaxisID1);
          auto usvarID = varID1;
          // find corresponding v-variable to u-variable:
          for (varID2 = 0; varID2 < numVars; ++varID2)
            if (vardata[varID2].empty())
            {
              // if ( cdoDebugExt ) cdo_print("Checking V-wind: vardata[%d]== nullptr ",varID1);
            }
            else  // This means that it is U or V.
            {
              auto const &var2 = varList2.vars[varID2];
              code2 = var2.code;
              zaxisID2 = var2.zaxisID;
              ltype2 = zaxis_to_ltype(zaxisID2);
              nlevel2 = var2.nlevels;
              CheckVarIsU(varID2, var2.name, code2);
              if (VarIsU) continue;
              if (cdoDebugExt >= 20)
                cdo_print("Checking V-wind: Var.id [%4d] with grib code:%3d; name: %6s; level type: %3d; number of "
                          "levels: %3d; zaxisID: %d",
                          varID2, code2, var2.name, ltype2, nlevel2, zaxisID2);
              CheckVarIsV(varID2, var2.name, code2);
              if (!VarIsV) continue;
              if (!((ltype1 == ltype2) && (nlevel1 == nlevel2) && (zaxisID1 == zaxisID2))) continue;
              if (cdoDebugExt >= 10)
                cdo_print("** FOUND V-wind; Var.id [%4d] with grib code:%3d; name: %6s; level type: %3d; number of "
                          "levels: %3d; zaxisID: %d",
                          varID2, code2, var2.name, ltype1, nlevel2, zaxisID2);
              auto vsvarID = varID2;
              if (cdoDebugExt >= 20)
                cdo_print("Using code %d [%d](u) and code %d [%d](v)", varList1.vars[varID1].code, code1,
                          varList1.vars[varID2].code, code2);
              gridID = varList1.vars[varID1].gridID;
              if (operatorID != ROTUVN)  // ROTUVN operator does not need creation of gridIDcurvl ...
              {
                if ((gridIDcurvl != -1) && (gridIDlastused != gridID))
                  cdo_abort("The gridID (%d) used just previously for uv-wind tranformation is not same this "
                            "time(%d)!",
                            gridIDlastused, gridID);

                if (gridIDcurvl == -1)
                {
                  if (cdoDebugExt) cdo_print("Building LAT-LON grid for the direction to the North. (First time only).");
                  gridIDlastused = gridID;
                  // Compute 2D array with latlons only once. We expect that all horizontal grids for UV are same.
                  // NOTE: At this stage U and V cannot be staggered!
                  gridIDcurvl = gridToCurvilinear(gridID, NeedCorners::Yes);
                  if (cdoDebugExt)
                    cdo_print("Transformed rotated-latLon grid (id:%d) to curvilinear (id:%d) with true lat-lon "
                              "coordinates.",
                              gridID, gridIDcurvl);
                  // Grid definition with id: "gridIDcurvl" contains latlons of every gridpoint..
                  // For details see: ./libcdi/src/cdi.h; Setgridtype to GRID_CURVILINEAR

                  if (gridIDcurvl == -1) cdo_abort("Creation of curvilinear grid definition failed!");

                  if (gridInqType(gridIDcurvl) != GRID_CURVILINEAR)
                  {
                    gridDestroy(gridIDcurvl);
                    cdo_abort("Creation of curvilinear grid definition failed: type != GRID_CURVILINEAR");
                  }
                  if (cdoDebugExt)
                  {
                    double xpole = 0, ypole = 0, angle = 0;
                    if (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL)
                      gridInqParamRLL(gridID, &xpole, &ypole, &angle);

                    cdo_print("GRID_PROJECTION(id: %d) && CDI_PROJ_RLL:", gridID);
                    cdo_print("grid Xsize   %zu, grid Ysize   %zu", gridInqXsize(gridID), gridInqYsize(gridID));
                    cdo_print("grid Xfirst  %4.3f, grid Yfirst  %4.3f", gridInqXval(gridID, 0), gridInqYval(gridID, 0));
                    cdo_print("grid Xinc   %4.3f, grid Yinc   %4.3f", gridInqXinc(gridID), gridInqYinc(gridID));
                    cdo_print("grid Xpole   %4.3f, grid Ypole   %4.3f", xpole, ypole);
                    cdo_print("GRID_CURVILINEAR (id: %d):", gridIDcurvl);
                    cdo_print("grid Xsize   %zu, grid Ysize   %zu", gridInqXsize(gridIDcurvl), gridInqYsize(gridIDcurvl));
                    cdo_print("grid Xfirst  %4.3f, grid Yfirst  %4.3f", gridInqXval(gridIDcurvl, 0), gridInqYval(gridIDcurvl, 0));
                    cdo_print("grid Xlast   %4.3f, grid Ylast   %4.3f", gridInqXval(gridIDcurvl, gridInqSize(gridIDcurvl) - 1),
                              gridInqYval(gridIDcurvl, gridInqSize(gridIDcurvl) - 1));
                    if (cdoDebugExt >= 20)
                    {
                      printf("Xvals (size=%zu):\n", gridInqSize(gridIDcurvl));
                      size_t ii;
                      for (ii = 0; ii < 10; ++ii) printf("%4.3f ", gridInqXval(gridIDcurvl, ii));
                      printf("\n...\n");
                      for (ii = gridInqSize(gridIDcurvl) - 10; ii < gridInqSize(gridIDcurvl); ++ii)
                        printf("%4.3f ", gridInqXval(gridIDcurvl, ii));
                      printf("\n");
                      printf("Yvals (size=%zu):\n", gridInqSize(gridIDcurvl));
                      for (ii = 0; ii < 10; ++ii) printf("%4.3f ", gridInqYval(gridIDcurvl, ii));
                      printf("\n...\n");
                      for (ii = gridInqSize(gridIDcurvl) - 10; ii < gridInqSize(gridIDcurvl); ++ii)
                        printf("%4.3f ", gridInqYval(gridIDcurvl, ii));
                      printf("\n");
                    }
                  }  // end of if (cdoDebugExt)
                  Debug(cdoDebugExt, "LAT-LON grid created.");
                }  // end of if (gridIDcurvl==-1)
              }  // end of if (operatorID != ROTUVN)

              int uRelativeToGrid = -1, vRelativeToGrid = -1;
              cdiInqKeyInt(vlistID1, varID1, CDI_KEY_UVRELATIVETOGRID, &uRelativeToGrid);
              cdiInqKeyInt(vlistID1, varID2, CDI_KEY_UVRELATIVETOGRID, &vRelativeToGrid);

              if (uRelativeToGrid != 1)
              {
                cdo_warning("U component is NOT relative to grid. No transformation to north-pole takes place!");
              }
              else if (vRelativeToGrid != 1)
              {
                cdo_warning("V component is NOT relative to grid. No transformation to north-pole takes place!");
              }
              else
              {
                for (int levelID = 0; levelID < nlevel1; ++levelID)
                {
                  if (cdoDebugExt)
                    cdo_print("RotuvNorth(): processing  level type: %d; level %d (out of [0:%d])", ltype1, levelID, nlevel1 - 1);
                  auto usvar = vardata[usvarID][levelID].data();
                  auto vsvar = vardata[vsvarID][levelID].data();
                  if (operatorID == ROTUVNORTH)
                  {
                    rot_uv_north(gridIDcurvl, usvar, vsvar);
                    // rot_uv_north(gridIDlastused, usvar, vsvar);
                    // transform "in-place" the uv from grid relative into north-pole related
                  }
                  else if (operatorID == PROJUVLATLON)
                    project_uv_latlon(gridIDcurvl, usvar, vsvar);
                  else if (operatorID == ROTUVN)
                    rot_uv_back_mode64(gridID, usvar, vsvar);
                }
              }

              Debug(cdoDebugExt, "Finished processing level type: %d", ltype1);
              break;
            }  // end  for varID2
        }  // end  for varID1

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto varID = fieldInfoList[fieldID].varID;
        auto levelID = fieldInfoList[fieldID].levelID;
        if (varID != -1)
        {
          code = varList1.vars[varID].code;
          zaxisID = varList1.vars[varID].zaxisID;
          ltype = zaxis_to_ltype(zaxisID);
          level = zaxisInqLevel(zaxisID, levelID);
          if (cdoDebugExt >= 10)
            cdo_print("Write modified data record: %05d (timestep:%d); Var.id [%4d]; (code=%3d; ltype=%3d; level=%4d; "
                      "levelID=%3d)",
                      fieldID, tsID, varID, code, ltype, level, levelID);

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, vardata[varID][levelID].data(), varnumMissVals[varID][levelID]);
        }
      }

      tsID++;
    }  // end of while (true)

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (gridIDcurvl != -1) gridDestroy(gridIDcurvl);  // at the end must Free the allocated curvilinear grid definition...

    return 0;
  }

public:
  void
  init() override
  {
    UVDESTAG = module.get_id("uvDestag");
    ROTUVNORTH = module.get_id("rotuvNorth");
    ROTUVN = module.get_id("rotuvN");
    // Cylindrical Equidistant projection
    PROJUVLATLON = module.get_id("projuvLatLon");

    operatorID = cdo_operator_id();
  }

  void
  run() override
  {

    // will be calling: rot_uv_north(int gridID, double *us, double *vs);
    if (operatorID == ROTUVNORTH) TransformUV();

    // will be calling: rot_uv_back_mode64(int gridID, double *us, double *vs);
    if (operatorID == ROTUVN) TransformUV();

    // will be calling: project_uv_latlon(int gridID, double *us, double *vs)
    if (operatorID == PROJUVLATLON) TransformUV();

    if (operatorID == UVDESTAG) DestaggerUV();

    cdo_abort("Unexpected operatorID %d", operatorID);
  }

  void
  close() override
  {
    // TODO: split the TransformUV and DestaggerUV into seperate classes and fix their close
  }
};
