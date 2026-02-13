/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Derivepar     gheight           Geopotential height on full-levels
      Derivepar     gheight_half      Geopotential height on half-levels
      Derivepar     sealevelpressure  Sea level pressure
*/

// variables derived from ECHAM or ERA data

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "const.h"
#include "cdo_zaxis.h"
#include "cdo_options.h"

void geopot_height_half(double *gheight, const double *ta_fl, const double *hus_fl, const double *p_hl, long ngp, long nlev);
void geopot_height_full(double *gheight, const double *ta_fl, const double *hus_fl, const double *p_hl, long ngp, long nlev);

static void
check_range_var2d(int stepNum, Varray<double> const &var2d, double rMin, double rMax, const char *varname)
{
  auto mm = varray_min_max(var2d);
  if (mm.min < rMin || mm.max > rMax)
    cdo_warning("%s out of range (min=%g max=%g) [timestep:%d]!", varname, mm.min, mm.max, stepNum);
}

static void
check_range_var3d(int stepNum, int nlevels, size_t gridsize, Varray<double> const &var3d, double rMin, double rMax,
                  const char *varname)
{
  static auto printWarning = true;
  if (printWarning)
  {
    double minVal = 1.e33, maxVal = -1.e33;
    for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      auto mm = varray_min_max(gridsize, &var3d[gridsize * levelID]);
      minVal = std::min(minVal, mm.min);
      maxVal = std::max(maxVal, mm.max);
    }

    if (minVal < rMin || maxVal > rMax)
    {
      printWarning = false;
      cdo_warning("%s out of range (min=%g max=%g) [timestep:%d]!", varname, minVal, maxVal, stepNum);
    }
  }
}

static void
compute_rho(int numFullLevels, size_t gridsize, Varray<double> const &pfull, Varray<double> const &ta, Varray<double> const &hus,
            Varray<double> &rho)
{
  //-- Specific gas constant for dry air
  constexpr double R_L = 287.085;  //[J/(kg*K)

  for (int levelID = 0; levelID < numFullLevels; ++levelID)
  {
    auto offset = levelID * gridsize;
    for (size_t i = 0; i < gridsize; ++i)
    {
      auto tv = ta[offset + i] * (1. + 0.6078 * hus[offset + i]);
      rho[offset + i] = pfull[offset + i] / (R_L * tv);
    }
  }
}

class Derivepar : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Derivepar",
    // clang-format off
    .operators = { { "gheight", DeriveparHelp },
                   { "gheight_half", DeriveparHelp },
                   { "air_density", DeriveparHelp },
                   { "sealevelpressure", DeriveparHelp } },
    // clang-format on
    .aliases = { { "gheight_full", "gheight" }, { "gheighthalf", "gheight_half" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Derivepar>();

private:
  int GHEIGHT_FULL{}, GHEIGHT_HALF{}, AIR_DENSITY{}, SEALEVELPRESSURE{};
  int surfaceID = -1;
  int presID = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  size_t gridsize{};
  VarIDs varIDs{};
  int zaxisID_ML = -1;
  int numHybridLevels = 0;
  int numFullLevels = 0;
  int numHalfLevels = 0;

  Varray<double> array;
  Varray<double> sgeopot;
  Varray<double> ps;
  Varray<double> ta;
  Varray<double> hus;
  Varray<double> gheight;
  Varray<double> halfPress;
  Varray<double> fullPress;
  Varray<double> rho;
  Varray<double> sealevelpressure;
  Varray<double> vct;

  VarList varList1;

  int operatorID{};

public:
  void
  init() override
  {
    GHEIGHT_FULL = module.get_id("gheight");
    GHEIGHT_HALF = module.get_id("gheight_half");
    AIR_DENSITY = module.get_id("air_density");
    SEALEVELPRESSURE = module.get_id("sealevelpressure");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    auto gridID0 = vlistGrid(vlistID1, 0);
    if (gridInqType(gridID0) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

    gridsize = vlist_check_gridsize(vlistID1);

    vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);
    int vctSize = vct.size();

    if (Options::cdoVerbose)
    {
      auto numAB = vctSize / 2;
      for (int i = 0; i < 10; ++i) cdo_print("vct: %5d %25.17f %25.17f", i, vct[i], vct[numAB + i]);
      cdo_print("vct:   ...");
      for (int i = numAB - 10; i < numAB; ++i) cdo_print("vct: %5d %25.17f %25.17f", i, vct[i], vct[numAB + i]);
    }

    if (zaxisID_ML == -1) cdo_abort("No 3D variable with hybrid sigma pressure coordinate found!");

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);

    auto numVars = varList1.numVars();
    auto const &vars1 = varList1.vars;

    varIDs = varList_search_varIDs(varList1, numFullLevels);

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.husID)     cdo_print("  %s -> %s", var_stdname(specific_humidity), vars1[varIDs.husID].name);
      if (-1 != varIDs.taID)      cdo_print("  %s -> %s", var_stdname(air_temperature), vars1[varIDs.taID].name);
      if (-1 != varIDs.psID)      cdo_print("  %s -> %s", var_stdname(surface_air_pressure), vars1[varIDs.psID].name);
      if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), vars1[varIDs.lnpsID].name);
      if (-1 != varIDs.sgeopotID) cdo_print("  %s -> %s", var_stdname(surface_geopotential), vars1[varIDs.sgeopotID].name);
      if (-1 != varIDs.geopotID)  cdo_print("  %s -> %s", var_stdname(geopotential), vars1[varIDs.geopotID].name);
      if (-1 != varIDs.gheightID) cdo_print("  %s -> %s", var_stdname(geopotential_height), vars1[varIDs.gheightID].name);
      // clang-format on
    }

    if (varIDs.lnpsID != -1 && varIDs.lnpsID2 != -1)
      cdo_abort("Found LOG(%s) twice: lsp and lnps!", var_stdname(surface_air_pressure));

    if (varIDs.taID == -1) cdo_abort("%s not found!", var_stdname(air_temperature));

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (operatorID == SEALEVELPRESSURE) varIDs.husID = -1;

      auto const &var = vars1[varID];
      if (var.gridType == GRID_SPECTRAL && var.zaxisType == ZAXIS_HYBRID) cdo_abort("Spectral data on model level unsupported!");
      if (var.gridType == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");
    }

    array.resize(gridsize);
    sgeopot.resize(gridsize);
    ps.resize(gridsize);
    ta.resize(gridsize * numFullLevels);
    halfPress.resize(gridsize * (numFullLevels + 1));

    auto shumidity = var_stdname(specific_humidity);
    if (operatorID == GHEIGHT_FULL || operatorID == GHEIGHT_HALF)
    {
      if (varIDs.husID == -1) { cdo_warning("%s not found - using algorithm without %s!", shumidity, shumidity); }
      else { hus.resize(gridsize * numFullLevels); }

      gheight.resize(gridsize * (numFullLevels + 1));
    }
    else if (operatorID == AIR_DENSITY)
    {
      if (varIDs.husID == -1) { cdo_abort("%s not found !", shumidity); }

      hus.resize(gridsize * numFullLevels);
      fullPress.resize(gridsize * numFullLevels);
      rho.resize(gridsize * numFullLevels);
    }
    else if (operatorID == SEALEVELPRESSURE)
    {
      fullPress.resize(gridsize * numFullLevels);

      surfaceID = zaxis_from_name("surface");
      sealevelpressure.resize(gridsize);
    }

    if (operatorID != AIR_DENSITY && zaxisID_ML != -1 && varIDs.sgeopotID == -1)
    {
      if (varIDs.geopotID == -1)
        cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
      else
        cdo_print("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));

      std::ranges::fill(sgeopot, 0.0);
    }

    presID = varIDs.lnpsID;
    if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      if (varIDs.psID == -1) { cdo_abort("%s not found!", var_stdname(surface_air_pressure)); }
      else { presID = varIDs.psID; }
    }

    if (Options::cdoVerbose)
    {
      if (presID == varIDs.lnpsID) { cdo_print("using LOG(%s)", var_stdname(surface_air_pressure)); }
      else { cdo_print("using %s", var_stdname(surface_air_pressure)); }
    }

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    {
      int var_id = -1;
      int varID = -1;

      if (operatorID == GHEIGHT_FULL)
      {
        var_id = geopotential_height;
        varID = vlistDefVar(vlistID2, gridID0, zaxisID_ML, TIME_VARYING);
      }
      else if (operatorID == GHEIGHT_HALF)
      {
        auto zaxisID_ML_Half = zaxisCreate(ZAXIS_HYBRID_HALF, numHalfLevels);
        zaxisDefVct(zaxisID_ML_Half, 2 * numHalfLevels, vct.data());
        Varray<double> levs(numHalfLevels);
        for (int i = 0; i < numHalfLevels; ++i) levs[i] = i + 1;
        zaxisDefLevels(zaxisID_ML_Half, levs.data());
        var_id = geopotential_height;
        varID = vlistDefVar(vlistID2, gridID0, zaxisID_ML_Half, TIME_VARYING);
      }
      else if (operatorID == AIR_DENSITY)
      {
        var_id = air_density;
        varID = vlistDefVar(vlistID2, gridID0, zaxisID_ML, TIME_VARYING);
      }
      else if (operatorID == SEALEVELPRESSURE)
      {
        var_id = air_pressure_at_sea_level;
        varID = vlistDefVar(vlistID2, gridID0, surfaceID, TIME_VARYING);
      }
      else
        cdo_abort("Internal problem, invalid operatorID: %d!", operatorID);

      vlistDefVarParam(vlistID2, varID, cdiEncodeParam(var_echamcode(var_id), 128, 255));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, var_name(var_id));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, var_stdname(var_id));
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, var_units(var_id));
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
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
        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);

        auto offset = gridsize * levelID;

        if (zaxisID_ML != -1)
        {
          if (varID == varIDs.sgeopotID) { varray_copy(gridsize, array, sgeopot); }
          else if (varID == varIDs.geopotID && varIDs.sgeopotID == -1 && (levelID + 1) == numFullLevels)
          {
            varray_copy(gridsize, array, sgeopot);
          }
          else if (varID == presID)
          {
            if (varIDs.lnpsID != -1)
              for (size_t i = 0; i < gridsize; ++i) ps[i] = std::exp(array[i]);
            else if (varIDs.psID != -1) { varray_copy(gridsize, array, ps); }
          }
          else if (varID == varIDs.taID) { array_copy(gridsize, array.data(), &ta[offset]); }
          else if (varID == varIDs.husID) { array_copy(gridsize, array.data(), &hus[offset]); }
        }
      }

      if (zaxisID_ML != -1)
      {
        // check range of psProg
        check_range_var2d(tsID + 1, ps, MIN_PS, MAX_PS, "Surface pressure");
        // check range of surface geopot
        check_range_var2d(tsID + 1, sgeopot, MIN_FIS, MAX_FIS, "Orography");
      }

      check_range_var3d(tsID + 1, varList1.vars[varIDs.taID].nlevels, gridsize, ta, MIN_T, MAX_T, "Temperature");

      if (varIDs.husID != -1)
        check_range_var3d(tsID + 1, varList1.vars[varIDs.husID].nlevels, gridsize, hus, -0.1, MAX_Q, "Humidity");

      if (operatorID == GHEIGHT_FULL || operatorID == GHEIGHT_HALF)
      {
        vct_to_hybrid_pressure((double *) nullptr, halfPress.data(), vct, ps.data(), numFullLevels, gridsize);
        array_copy(gridsize, sgeopot.data(), gheight.data() + gridsize * numFullLevels);
        if (operatorID == GHEIGHT_FULL)
          geopot_height_full(gheight.data(), ta.data(), hus.data(), halfPress.data(), gridsize, numFullLevels);
        else
          geopot_height_half(gheight.data(), ta.data(), hus.data(), halfPress.data(), gridsize, numFullLevels);

        int varID = 0;
        auto numLevels = (operatorID == GHEIGHT_FULL) ? numFullLevels : numHalfLevels;
        for (int levelID = 0; levelID < numLevels; ++levelID)
        {
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, gheight.data() + levelID * gridsize, 0);
        }
      }
      else if (operatorID == AIR_DENSITY)
      {
        vct_to_hybrid_pressure(fullPress.data(), halfPress.data(), vct, ps.data(), numFullLevels, gridsize);

        int varID = 0;
        auto numLevels = numFullLevels;
        compute_rho(numFullLevels, gridsize, fullPress, ta, hus, rho);
        for (int levelID = 0; levelID < numLevels; ++levelID)
        {
          auto offset = levelID * gridsize;
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, rho.data() + offset, 0);
        }
      }
      else if (operatorID == SEALEVELPRESSURE)
      {
        vct_to_hybrid_pressure(fullPress.data(), halfPress.data(), vct, ps.data(), numFullLevels, gridsize);

        extrapolate_P(sealevelpressure.data(), &halfPress[gridsize * numFullLevels], &fullPress[gridsize * (numFullLevels - 1)],
                      sgeopot.data(), &ta[gridsize * (numFullLevels - 1)], gridsize);

        cdo_def_field(streamID2, 0, 0);
        cdo_write_field(streamID2, sealevelpressure.data(), 0);
      }
      else
        cdo_abort("Internal error");

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
