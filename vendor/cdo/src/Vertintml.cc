/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertint    ml2pl           Model to pressure level interpolation
      Vertint    ml2hl           Model to height level interpolation
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "constants.h"
#include "const.h"
#include "param_conversion.h"
#include "vertint_util.h"

/*
static void
field_copy_div_2d_to_3d(MemType memType, size_t gridsize, int nlevels, Field const &field2d, Field3D &field3d)
{
  if (memType == MemType::Float)
    for (size_t i = 0; i < gridsize; ++i) field3d.vec_f[gridsize * nlevels + i] = field2d.vec_f[i] / PlanetGrav;
  else
    for (size_t i = 0; i < gridsize; ++i) field3d.vec_d[gridsize * nlevels + i] = field2d.vec_d[i] / PlanetGrav;
}
*/
template <typename T>
static void
field_check_sgeopot(Varray<T> const &sgeopot, const T *gheightAtSurface)
{
  auto len = sgeopot.size();
  double sumDiff = 0;
  for (size_t i = 0; i < len; ++i) sumDiff += std::fabs((sgeopot[i] / PlanetGrav) - gheightAtSurface[i]);

  constexpr double lim = 0.1;  // 10cm per grid point
  if ((sumDiff / len) > lim) cdo_warning("Bottom level of gheight differ from sgeopot/%g (diff=%g)!", PlanetGrav, sumDiff / len);
  // printf("sumDiff %g %g\n", sumDiff, sumDiff / len);
}

static void
field_check_sgeopot(MemType memType, size_t gridsize, int nlevels, Field const &field2d, const Field3D &field3d)
{
  if (memType == MemType::Float)
    field_check_sgeopot(field2d.vec_f, &field3d.vec_f[gridsize * (nlevels - 1)]);
  else
    field_check_sgeopot(field2d.vec_d, &field3d.vec_d[gridsize * (nlevels - 1)]);
}

static void
vct_to_hybrid_pressure(MemType memType, Field3D &pressure_FL, Field3D &pressure_HL, Varray<double> const &vct, Field const &ps,
                       long numHybridLevels, long ngp)
{
  if (memType == MemType::Float)
    vct_to_hybrid_pressure(pressure_FL.vec_f.data(), pressure_HL.vec_f.data(), vct, ps.vec_f.data(), numHybridLevels, ngp);
  else
    vct_to_hybrid_pressure(pressure_FL.vec_d.data(), pressure_HL.vec_d.data(), vct, ps.vec_d.data(), numHybridLevels, ngp);
}

static void
invert_vct(Varray<double> &vct)
{
  Varray<double> vctbuf = vct;
  auto vctSize = vct.size();
  for (size_t i = 0; i < vctSize / 2; ++i)
  {
    vct[vctSize / 2 - 1 - i] = vctbuf[i];
    vct[vctSize - 1 - i] = vctbuf[i + vctSize / 2];
  }
}

template <class ForwardIt>
static void
check_range(double rmin, double rmax, ForwardIt first, ForwardIt last, std::string const &name)
{
  auto [pmin, pmax] = std::minmax_element(first, last);
  if (std::fabs(*pmax - *pmin) <= 1.0e-9 || (*pmin < rmin) || (*pmax > rmax))
    cdo_warning("%s out of range (%g - %g) min=%g max=%g", name, rmin, rmax, *pmin, *pmax);
}

static void
check_vct(Varray<double> const &vct, int numHalfLevels)
{
  check_range(0.0, 50000.0, vct.begin(), vct.begin() + numHalfLevels, "vct A");
  check_range(0.0, 1.0, vct.begin() + numHalfLevels, vct.end(), "vct B");
}

static void
check_range_ps(int stepNum, Field const &psProg)
{
  auto mm = field_min_max(psProg);
  if (mm.min < MIN_PS || mm.max > MAX_PS)
    cdo_warning("Surface pressure out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static void
check_range_sgeopot(int stepNum, Field const &sgeopot)
{
  auto mm = field_min_max(sgeopot);
  if (mm.min < MIN_FIS || mm.max > MAX_FIS)
    cdo_warning("Surface geopotential out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
  if (sgeopot.gridsize > 1 && mm.min >= 0.0 && mm.max <= 9000.0 && is_not_equal(mm.min, mm.max))
    cdo_warning("Surface geopotential has an unexpected range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static bool
zaxis_is_hybrid(int zaxisType)
{
  return (zaxisType == ZAXIS_HYBRID || zaxisType == ZAXIS_HYBRID_HALF);
}

static void
change_hybrid_zaxis(int vlistID1, int vlistID2, int vctSize, const double *vct, int zaxisID2, int numFullLevels, int numHalfLevels)
{
  auto numZaxes = vlistNumZaxis(vlistID1);
  for (int iz = 0; iz < numZaxes; ++iz)
  {
    auto zaxisID = vlistZaxis(vlistID1, iz);
    auto nlevels = zaxisInqSize(zaxisID);
    auto zaxisType = zaxisInqType(zaxisID);

    if (zaxis_is_hybrid(zaxisType) && (nlevels == numHalfLevels || nlevels == numFullLevels))
    {
      auto vctSize2 = zaxisInqVctSize(zaxisID);
      if (vctSize2 == vctSize && std::memcmp(vct, zaxisInqVctPtr(zaxisID), vctSize * sizeof(double)) == 0)
        vlistChangeZaxisIndex(vlistID2, iz, zaxisID2);
    }
  }
}

static void
print_vars_found(const VarIDs &varIDs, const CdoVars &vars)
{
  cdo_print("Found:");
  // clang-format off
  if (-1 != varIDs.taID)      cdo_print("  %s -> %s", var_stdname(air_temperature), vars[varIDs.taID].name);
  if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), vars[varIDs.psID].name);
  if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), vars[varIDs.lnpsID].name);
  if (-1 != varIDs.sgeopotID) cdo_print("  %s -> %s", var_stdname(surface_geopotential), vars[varIDs.sgeopotID].name);
  if (-1 != varIDs.geopotID)  cdo_print("  %s -> %s", var_stdname(geopotential), vars[varIDs.geopotID].name);
  if (-1 != varIDs.gheightID) cdo_print("  %s -> %s", var_stdname(geopotential_height), vars[varIDs.gheightID].name);
  // clang-format on
}

static void
pressure_level_interpolation(Varray<double> &pressureLevels, bool useHeightLevel, bool extrapolate)
{
  int numPL = pressureLevels.size();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList1(vlistID1);
  varList_set_unique_memtype(varList1);
  auto memType = varList1.vars[0].memType;

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxisID_PL = zaxisCreate(useHeightLevel ? ZAXIS_ALTITUDE : ZAXIS_PRESSURE, numPL);
  zaxisDefLevels(zaxisID_PL, pressureLevels.data());

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
  int vctSize = vct.size();

  // check VCT
  if (zaxisID_ML != -1) check_vct(vct, numHalfLevels);

  change_hybrid_zaxis(vlistID1, vlistID2, vctSize, vct.data(), zaxisID_PL, numFullLevels, numHalfLevels);

  VarList varList2(vlistID2);
  varList_set_memtype(varList2, memType);

  int psvarID = -1;
  auto vctIsInverted = false;
  if (vctSize && vctSize % 2 == 0)
  {
    psvarID = varList_get_psvarid(varList1, zaxisID_ML);

    int i = vctSize / 2 + 1;
    for (; i < vctSize; ++i)
      if (vct[i] > vct[i - 1]) break;
    if (i == vctSize) vctIsInverted = true;
  }

  if (Options::cdoVerbose) cdo_print("vctIsInverted = %d", static_cast<int>(vctIsInverted));

  if (vctIsInverted) invert_vct(vct);

  auto numVars = varList1.numVars();

  std::vector<bool> processVars(numVars), interpVars(numVars);
  Varray2D<size_t> varnumMissVals(numVars);
  Field3DVector vardata1(numVars), vardata2(numVars);

  auto maxLevels = std::max(numHalfLevels, numPL);

  Varray<size_t> pnumMissVals;
  if (!extrapolate) pnumMissVals.resize(numPL);

  // check levels
  if (zaxisID_ML != -1)
  {
    if (zaxisInqSize(zaxisID_ML) != numHybridLevels) cdo_abort("Internal error, wrong number of hybrid level!");
  }

  Field3D pressure_FL, pressure_HL;
  if (zaxisID_ML != -1 && gridsize > 0)
  {
    CdoVar var3Df, var3Dh;
    var3Df.gridsize = gridsize;
    var3Df.nlevels = numFullLevels;
    var3Df.memType = memType;
    pressure_FL.init(var3Df);

    var3Dh.gridsize = gridsize;
    var3Dh.nlevels = numHalfLevels;
    var3Dh.memType = memType;
    pressure_HL.init(var3Dh);
  }
  else
    cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!");

  if (useHeightLevel)
  {
    std::vector<double> phlev(numPL);
    height_to_pressure(pressureLevels.data(), phlev.data(), numPL);

    if (Options::cdoVerbose)
      for (int i = 0; i < numPL; ++i) cdo_print("level=%d  height=%g  pressure=%g", i + 1, pressureLevels[i], phlev[i]);

    pressureLevels = phlev;
  }

  auto varIDs = varList_search_varIDs(varList1, numFullLevels);

  // vertical_interp_Z() is implemented for gheight on model half levels only
  if (-1 != varIDs.gheightID && varList1.vars[varIDs.gheightID].nlevels == numFullLevels)
  {
    if (Options::cdoVerbose)
      cdo_print("%s(%s) on model full levels found!", var_stdname(geopotential_height), varList1.vars[varIDs.gheightID].name);
    varIDs.gheightID = -1;
  }

  if (Options::cdoVerbose) print_vars_found(varIDs, varList1.vars);

  for (int varID = 0; varID < numVars; ++varID)
  {
    auto &var1 = varList1.vars[varID];
    auto nlevels = var1.nlevels;

    if (var1.gridType == GRID_SPECTRAL && zaxis_is_hybrid(var1.zaxisType)) cdo_abort("Spectral data on model level unsupported!");

    if (var1.gridType == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

    // if (varID == varIDs.gheightID) var1.nlevels = nlevels + 1;
    vardata1[varID].init(var1);
    // if (varID == varIDs.gheightID) var1.nlevels = nlevels;

    // interpVars[varID] = (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && nlevels == numHybridLevels);
    interpVars[varID]
        = (var1.zaxisID == zaxisID_ML
           || (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && (nlevels == numHalfLevels || nlevels == numFullLevels)));

    if (interpVars[varID])
    {
      varnumMissVals[varID].resize(maxLevels, 0);
      vardata2[varID].init(varList2.vars[varID]);
    }
    else
    {
      varnumMissVals[varID].resize(nlevels);
      if (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && nlevels > 1)
        cdo_warning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)", varID + 1, var1.name, nlevels);
    }
  }

  auto needVertIndexHalf{ false };
  for (int varID = 0; varID < numVars; ++varID)
  {
    if (interpVars[varID])
    {
      auto const &var1 = varList1.vars[varID];
      if (var1.nlevels == numHalfLevels && varID != varIDs.gheightID) needVertIndexHalf = true;
    }
  }

  std::vector<int> vertIndex_FL;
  std::vector<int> vertIndex_HL;
  if (zaxisID_ML != -1 && gridsize > 0)
  {
    auto num3D = gridsize * numPL;
    if (num3D > std::numeric_limits<int>::max()) cdo_abort("gridSize*numLevels=%zu exceeds the limits of integer.", num3D);
    vertIndex_FL.resize(num3D);
    if (needVertIndexHalf) vertIndex_HL.resize(num3D);
  }

  if (zaxisID_ML != -1 && varIDs.gheightID != -1 && varIDs.taID == -1)
    cdo_abort("%s not found, needed for vertical interpolation of %s!", var_stdname(air_temperature),
              var_stdname(geopotential_height));

  auto presID = (psvarID != -1) ? psvarID : varIDs.lnpsID;

  if (zaxisID_ML != -1 && presID == -1)
  {
    if (varIDs.psID == -1) cdo_abort("%s not found!", var_stdname(surface_air_pressure));
    presID = varIDs.psID;
  }

  if (Options::cdoVerbose && presID != -1)
  {
    if (presID == varIDs.lnpsID)
      cdo_print("Using LOG(%s) from %s", var_stdname(surface_air_pressure), varList1.vars[presID].name);
    else
      cdo_print("Using %s from %s", var_stdname(surface_air_pressure), varList1.vars[presID].name);
  }

  Field psProg;
  if (zaxisID_ML != -1 && presID != -1) psProg.init(varList1.vars[presID]);

  auto sgeopotNeeded = (varIDs.taID != -1 || varIDs.gheightID != -1);

  Field sgeopot;
  if (zaxisID_ML != -1 && sgeopotNeeded)
  {
    sgeopot.init(varList1.vars[presID]);
    if (varIDs.sgeopotID == -1)
    {
      if (extrapolate)
      {
        if (varIDs.geopotID == -1)
          cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
        else
          cdo_print("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));
      }
      field_fill(sgeopot, 0.0);
    }
  }

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
    if (numFields == 0) break;

    for (int varID = 0; varID < numVars; ++varID) processVars[varID] = false;

    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsID);

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto const &var = varList1.vars[varID];

      if (vctIsInverted && zaxisID_ML != -1 && var.zaxisID == zaxisID_ML) levelID = var.nlevels - 1 - levelID;

      cdo_read_field(streamID1, vardata1[varID], levelID, &varnumMissVals[varID][levelID]);

      processVars[varID] = true;
    }

    if (zaxisID_ML != -1)
    {
      if (sgeopotNeeded)
      {
        if (varIDs.sgeopotID != -1)
          field_copy(vardata1[varIDs.sgeopotID], sgeopot);
        else if (varIDs.geopotID != -1)
          field_copy(vardata1[varIDs.geopotID], numFullLevels - 1, sgeopot);

        // check range of surface geopot
        if (extrapolate && (varIDs.sgeopotID != -1 || varIDs.geopotID != -1)) check_range_sgeopot(tsID + 1, sgeopot);
      }

      if (presID == varIDs.lnpsID)
        field_transform(vardata1[varIDs.lnpsID], psProg, unary_op_exp);
      else if (presID != -1)
        field_copy(vardata1[presID], psProg);

      // check range of psProg
      check_range_ps(tsID + 1, psProg);

      vct_to_hybrid_pressure(memType, pressure_FL, pressure_HL, vct, psProg, numFullLevels, gridsize);

      gen_vert_index(vertIndex_FL, pressureLevels, pressure_FL, gridsize);
      if (!extrapolate) gen_vert_index_mv(vertIndex_FL, pressureLevels, gridsize, psProg, pnumMissVals);

      if (needVertIndexHalf)
      {
        gen_vert_index(vertIndex_HL, pressureLevels, pressure_HL, gridsize);
        if (!extrapolate) gen_vert_index_mv(vertIndex_HL, pressureLevels, gridsize, psProg, pnumMissVals);
      }
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (processVars[varID])
      {
        auto const &var = varList1.vars[varID];

        if (tsID > 0 && var.isConstant) continue;

        if (interpVars[varID])
        {
          if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
            cdo_abort("Number of hybrid level differ from full/half level (param=%s)!", var.name);

          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
          }

          if (varID == varIDs.taID)
          {
            if (var.nlevels == numHalfLevels) cdo_abort("Temperature on half level unsupported!");

            vertical_interp_T(var.nlevels, pressure_FL, pressure_HL, vardata1[varID], vardata2[varID], sgeopot, vertIndex_FL,
                              pressureLevels, gridsize);
          }
          else if (varID == varIDs.gheightID)
          {
            // field_copy_div_2d_to_3d(memType, gridsize, var.nlevels, sgeopot, vardata1[varID]);
            field_check_sgeopot(memType, gridsize, var.nlevels, sgeopot, vardata1[varID]);

            vertical_interp_Z(numFullLevels, pressure_FL, pressure_HL, vardata1[varID], vardata2[varID], vardata1[varIDs.taID],
                              sgeopot, vertIndex_FL, pressureLevels, gridsize);
          }
          else
          {
            auto const &levels3D = (var.nlevels == numFullLevels) ? pressure_FL : pressure_HL;
            auto const &vertIndex3D = (var.nlevels == numFullLevels) ? vertIndex_FL : vertIndex_HL;
            vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex3D, pressureLevels, gridsize);
          }

          if (!extrapolate) varray_copy(numPL, pnumMissVals, varnumMissVals[varID]);
        }

        for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
        {
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                          varnumMissVals[varID][levelID]);
        }
      }
    }

    tsID++;
  }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}

// #define ENABLE_HEIGHT_LEVEL_INTERPOLATION

#ifdef ENABLE_HEIGHT_LEVEL_INTERPOLATION
// obsolete, use intlevel!!!
static void
height_level_interpolation(Varray<double> &heightLevels, bool extrapolate)
{
  int numHeightLevels = heightLevels.size();

  auto streamID1 = cdo_open_read(0);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = vlistDuplicate(vlistID1);

  VarList varList1(vlistID1);
  varList_setUniqueMemtype(varList1);
  auto memType = varList1.vars[0].memType;

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto gridsize = vlist_check_gridsize(vlistID1);

  auto zaxisID_HL = zaxisCreate(ZAXIS_HEIGHT, numHeightLevels);
  zaxisDefLevels(zaxisID_HL, heightLevels.data());

  int zaxisID_ML = -1;
  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;
  auto vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
  int vctSize = vct.size();

  // check VCT
  if (zaxisID_ML != -1) check_vct(vct, numHalfLevels);

  change_hybrid_zaxis(vlistID1, vlistID2, vctSize, vct.data(), zaxisID_HL, numFullLevels, numHalfLevels);

  VarList varList2(vlistID2);
  varList_setMemtype(varList2, memType);
  /*
  auto vctIsInverted = false;
  if (vctSize && vctSize % 2 == 0)
    {
      (void) cdoVars_get_psvarid(varList1.vars, zaxisID_ML);

      int i;
      for (i = vctSize / 2 + 1; i < vctSize; ++i)
        if (vct[i] > vct[i - 1]) break;
      if (i == vctSize) vctIsInverted = true;
    }

  if (Options::cdoVerbose) cdo_print("vctIsInverted = %d", static_cast<int>(vctIsInverted));

  if (vctIsInverted) invert_vct(vct);
  */
  auto numVars = varList1.numVars();

  std::vector<bool> processVars(numVars), interpVars(numVars);
  Varray2D<size_t> varnumMissVals(numVars);
  Field3DVector vardata1(numVars), vardata2(numVars);

  auto maxLevels = std::max(numHalfLevels, numHeightLevels);

  Varray<size_t> pnumMissVals;
  if (!extrapolate) pnumMissVals.resize(numHeightLevels);

  // check levels
  if (zaxisID_ML != -1)
  {
    auto nlev = zaxisInqSize(zaxisID_ML);
    if (nlev != numHybridLevels) cdo_abort("Internal error, wrong number of hybrid level!");
  }

  std::vector<int> vertIndex;
  if (zaxisID_ML != -1 && gridsize > 0) { vertIndex.resize(gridsize * numHeightLevels); }
  else { cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!"); }

  VarIDs varIDs = search_varIDs(varList1.vars, vlistID1, numFullLevels);

  if (Options::cdoVerbose) print_vars_found(varIDs, varList1.vars);

  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var1 = varList1.vars[varID];
    auto nlevels = var1.nlevels;

    if (var1.gridType == GRID_SPECTRAL && zaxis_is_hybrid(var1.zaxisType)) cdo_abort("Spectral data on model level unsupported!");
    if (var1.gridType == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

    vardata1[varID].init(var1);

    // interpVars[varID] = (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && nlevels == numHybridLevels);
    interpVars[varID]
        = (var1.zaxisID == zaxisID_ML
           || (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && (nlevels == numHalfLevels || nlevels == numFullLevels)));

    if (interpVars[varID])
    {
      varnumMissVals[varID].resize(maxLevels, 0);
      vardata2[varID].init(varList2.vars[varID]);
    }
    else
    {
      varnumMissVals[varID].resize(nlevels);
      if (zaxis_is_hybrid(var1.zaxisType) && zaxisID_ML != -1 && nlevels > 1)
        cdo_warning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)", varID + 1, var1.name, nlevels);
    }
  }

  if (zaxisID_ML != -1 && varIDs.gheightID == -1) cdo_abort("%s not found!", var_stdname(geopotential_height));
  /*
  auto sgeopotNeeded = (!extrapolate && varIDs.gheightID != -1);

  Field sgeopot;
  if (zaxisID_ML != -1 && sgeopotNeeded)
    {
      if (varIDs.sgeopotID == -1)
        {
          CdoVar var2D;
          var2D.gridsize = gridsize;
          var2D.nlevels = 1;
          var2D.memType = memType;
          sgeopot.init(var2D);

          if (extrapolate) cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
        }
      else { sgeopot.init(varList1.vars[varIDs.sgeopotID]); }

      field_fill(sgeopot, 0.0);
    }
  */
  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  Field heightBottom;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
    if (numFields == 0) break;

    for (int varID = 0; varID < numVars; ++varID) processVars[varID] = false;

    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsID);

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      cdo_read_field(streamID1, vardata1[varID], levelID, &varnumMissVals[varID][levelID]);

      processVars[varID] = true;
    }

    if (zaxisID_ML != -1)
    {
      auto lreverse = true;
      gen_vert_index(vertIndex, heightLevels, vardata1[varIDs.gheightID], gridsize, lreverse);
      if (!extrapolate)
      {
        heightBottom.init(varList1.vars[varIDs.gheightID]);
        field_copy(vardata1[varIDs.gheightID], numFullLevels - 1, heightBottom);
        gen_vert_index_mv(vertIndex, heightLevels, gridsize, heightBottom, pnumMissVals, lreverse);
      }
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (processVars[varID])
      {
        auto const &var = varList1.vars[varID];
        if (tsID > 0 && var.isConstant) continue;

        if (interpVars[varID])
        {
          if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
            cdo_abort("Number of hybrid level differ from full/half level (param=%s)!", var.name);

          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
          }

          auto const &levels3D = vardata1[varIDs.gheightID];
          vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex, heightLevels, gridsize);

          if (!extrapolate) varray_copy(numHeightLevels, pnumMissVals, varnumMissVals[varID]);
        }

        for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
        {
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                          varnumMissVals[varID][levelID]);
        }
      }
    }

    tsID++;
  }

  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}
#endif

class Vertintml : public Process
{
  enum
  {
    func_pl,
    func_hl
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertintml",
    // clang-format off
    .operators = { { "ml2pl", func_pl, 0, "pressure levels in pascal", VertintmlHelp },
                   { "ml2hl", func_hl, 0, "height levels in meter", VertintmlHelp },
                   { "ml2plx", func_pl, 0, "pressure levels in pascal", VertintmlHelp },
                   { "ml2hlx", func_hl, 0, "height levels in meter", VertintmlHelp },
#ifdef ENABLE_HEIGHT_LEVEL_INTERPOLATION
                   { "ml2height", func_hl, 1, "height levels in meter", VertintmlHelp },
                   { "ml2heightx", func_hl, 1, "height levels in meter", VertintmlHelp }
#endif
                 },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertintml> registration = RegisterEntry<Vertintml>(module);

  bool doHeightInterpolation{ false };
  bool useHeightLevel{ false };
  bool extrapolate{ false };

  Varray<double> levels;

public:
  void
  init() override
  {
    auto ML2PLX = module.get_id("ml2plx");
    auto ML2HLX = module.get_id("ml2hlx");
#ifdef ENABLE_HEIGHT_LEVEL_INTERPOLATION
    auto ML2HEIGHTX = module.get_id("ml2heightx");
#endif

    auto operatorID = cdo_operator_id();
    useHeightLevel = (cdo_operator_f1(operatorID) == func_hl);
    doHeightInterpolation = (cdo_operator_f2(operatorID) == 1);

#ifdef ENABLE_HEIGHT_LEVEL_INTERPOLATION
    extrapolate = (operatorID == ML2PLX || operatorID == ML2HLX || operatorID == ML2HEIGHTX);
#else
    extrapolate = (operatorID == ML2PLX || operatorID == ML2HLX);
#endif
    if (extrapolate == false) extrapolate = getenv_extrapolate();

    operator_input_arg(cdo_operator_enter(operatorID));

    if (cdo_operator_argc() == 1 && cdo_operator_argv(0) == "default")
    {
      if (useHeightLevel)
        levels = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
      else
        levels
            = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
    }
    else { levels = cdo_argv_to_fltarr(cdo_get_oper_argv()); }
  }

  void
  run() override
  {
#ifdef ENABLE_HEIGHT_LEVEL_INTERPOLATION
    if (doHeightInterpolation)
      height_level_interpolation(levels, extrapolate);
    else
#endif
      pressure_level_interpolation(levels, useHeightLevel, extrapolate);
  }

  void
  close() override
  {
  }
};
