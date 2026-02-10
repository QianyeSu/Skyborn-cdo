/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pressure    pressure             Pressure on model full-levels
      Pressure    pressure_half        Pressure on model half-levels
      Pressure    delta_pressure       Difference of two pressure half-levels
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "const.h"

class Pressure : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Pressure",
    // clang-format off
    .operators = { { "pressure", PressureHelp },
                   { "pressure_half", PressureHelp },
                   { "delta_pressure", PressureHelp } },
    // clang-format on
    .aliases = { { "pressure_full", "pressure" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Pressure> registration = RegisterEntry<Pressure>(module);

private:
  int PRESSURE_FULL{}, PRESSURE_HALF{}, DELTA_PRESSURE{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int pvarID{};

  int operatorID{};

  size_t gridsize{};
  int zaxisID_ML = -1;

  VarIDs varIDs;
  Varray<double> psProg;
  Varray<double> vct;
  Varray<double> array;
  Varray<double> deltaPressure;
  Varray<double> pressureFull;
  Varray<double> pressureHalf;

  int numHybridLevels = 0, numFullLevels = 0, numHalfLevels = 0;

public:
  void
  init() override
  {
    PRESSURE_FULL = module.get_id("pressure");
    PRESSURE_HALF = module.get_id("pressure_half");
    DELTA_PRESSURE = module.get_id("delta_pressure");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1(vlistID1);
    // varList_setUniqueMemtype(varList1);
    // auto memType = varList1.vars[0].memType;

    gridsize = vlist_check_gridsize(vlistID1);

    vct = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels, numHalfLevels);

    auto hasVars3D = (zaxisID_ML != -1 && gridsize > 0);
    if (!hasVars3D) cdo_abort("No 3D variable with hybrid sigma pressure coordinate found!");

    psProg = Varray<double>(gridsize);
    deltaPressure = Varray<double>(gridsize * numFullLevels);
    pressureFull = Varray<double>(gridsize * numFullLevels);
    pressureHalf = Varray<double>(gridsize * numHalfLevels);

    int zaxisID_PL = -1;
    int zaxisType = -1;
    if (operatorID == PRESSURE_FULL || operatorID == DELTA_PRESSURE)
    {
      if (numFullLevels == numHybridLevels)
        zaxisID_PL = zaxisID_ML;
      else
      {
        if (Options::cdoVerbose) cdo_print("Creating ZAXIS_HYBRID .. (numFullLevels=%d)", numFullLevels);
        zaxisType = ZAXIS_HYBRID;
        numHybridLevels = numFullLevels;
      }
    }
    else
    {
      if (numHalfLevels == numHybridLevels)
        zaxisID_PL = zaxisID_ML;
      else
      {
        if (Options::cdoVerbose) cdo_print("Creating ZAXIS_HYBRID_HALF .. (numHalfLevels=%d)", numHalfLevels);
        zaxisType = ZAXIS_HYBRID_HALF;
        numHybridLevels = numHalfLevels;
      }
    }

    if (zaxisID_PL == -1)
    {
      zaxisID_PL = zaxisCreate(zaxisType, numHybridLevels);

      Varray<double> level(numHalfLevels);
      for (int l = 0; l < numHalfLevels; ++l) level[l] = l + 1;
      zaxisDefLevels(zaxisID_PL, level.data());
      zaxisDefVct(zaxisID_PL, 2 * numHalfLevels, vct.data());
    }

    varIDs = varList_search_varIDs(varList1, numFullLevels);

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != varIDs.psID)   cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1.vars[varIDs.psID].name);
      if (-1 != varIDs.lnpsID) cdo_print("  LOG(%s) -> %s", var_stdname(surface_air_pressure), varList1.vars[varIDs.lnpsID].name);
      // clang-format on
    }

    pvarID = varIDs.lnpsID;
    if (zaxisID_ML != -1 && varIDs.lnpsID != -1)
    {
      auto gridType = varList1.vars[varIDs.lnpsID].gridType;
      if (gridType == GRID_SPECTRAL)
      {
        varIDs.lnpsID = -1;
        cdo_warning("Spectral LOG(%s) not supported - using %s!", var_stdname(surface_air_pressure),
                    var_stdname(surface_air_pressure));
      }
    }

    if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      pvarID = varIDs.psID;
      if (varIDs.psID == -1) cdo_abort("%s not found!", var_stdname(surface_air_pressure));
    }

    auto gridID = varList1.vars[pvarID].gridID;
    if (gridInqType(gridID) == GRID_SPECTRAL)
      cdo_abort("%s on spectral representation not supported!", var_stdname(surface_air_pressure));

    array = Varray<double>(gridsize);

    auto vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    auto ovarID = vlistDefVar(vlistID2, gridID, zaxisID_PL, TIME_VARYING);
    vlistDefVarParam(vlistID2, ovarID, cdiEncodeParam(1, 255, 255));
    cdiDefKeyString(vlistID2, ovarID, CDI_KEY_NAME, "pressure");
    cdiDefKeyString(vlistID2, ovarID, CDI_KEY_STDNAME, "air_pressure");
    cdiDefKeyString(vlistID2, ovarID, CDI_KEY_UNITS, "Pa");

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
        if (varID == pvarID)
        {
          size_t numMissVals;
          cdo_read_field(streamID1, array.data(), &numMissVals);
          if (numMissVals) cdo_abort("Missing valus unsupported!");
        }
      }

      if (zaxisID_ML != -1)
      {
        if (varIDs.lnpsID != -1)
          for (size_t i = 0; i < gridsize; ++i) psProg[i] = std::exp(array[i]);
        else if (varIDs.psID != -1)
          varray_copy(gridsize, array, psProg);

        // check range of psProg
        auto mm = varray_min_max(psProg);
        if (mm.min < MIN_PS || mm.max > MAX_PS) cdo_warning("Surface pressure out of range (min=%g max=%g)!", mm.min, mm.max);

        vct_to_hybrid_pressure(pressureFull.data(), pressureHalf.data(), vct, psProg.data(), numFullLevels, gridsize);
      }

      double *pout = nullptr;
      int nlevels = 0;
      if (operatorID == PRESSURE_FULL)
      {
        nlevels = numFullLevels;
        pout = pressureFull.data();
      }
      else if (operatorID == DELTA_PRESSURE)
      {
        nlevels = numFullLevels;
        for (int k = 0; k < numFullLevels; ++k)
          for (size_t i = 0; i < gridsize; ++i)
            deltaPressure[k * gridsize + i] = pressureHalf[(k + 1) * gridsize + i] - pressureHalf[k * gridsize + i];

        pout = deltaPressure.data();
      }
      else if (operatorID == PRESSURE_HALF)
      {
        nlevels = numHalfLevels;
        pout = pressureHalf.data();
      }

      int varID = 0;
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, &pout[levelID * gridsize], 0);
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
