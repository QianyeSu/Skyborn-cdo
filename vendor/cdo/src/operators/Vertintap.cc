/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     ap2pl           Model air pressure level to pressure level interpolation
*/

#include "cdo_options.h"
#include "cdo_output.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "field_vinterp.h"
#include "stdnametable.h"
#include "util_string.h"
#include "const.h"
#include "param_conversion.h"
#include "vertint_util.h"

static void
check_range_ps(int stepNum, Field const &psProg)
{
  auto mm = field_min_max(psProg);
  if (mm.min < MIN_PS || mm.max > MAX_PS)
    cdo_warning("Surface pressure out of range (min=%g max=%g) [timestep:%d]!", mm.min, mm.max, stepNum);
}

static bool
is_height_axis(int zaxisID)
{
  auto isHeight = false;
  if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
  {
    // auto units = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
    auto stdname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_STDNAME);
    // if (stdname == "height" && units.empty()) isHeight = true;
    if (stdname == "height") isHeight = true;
  }
  return isHeight;
}

template <typename T>
static void
calc_half_press(size_t gridsize, int numFullLevels, Varray<T> const &fullPress, int numHalfLevels, Varray<T> &halfPress)
{
  for (size_t i = 0; i < gridsize; ++i) halfPress[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (int k = 1; k < numFullLevels; ++k)
  {
    auto fullPress_km1 = &fullPress[(k - 1) * gridsize];
    auto fullPress_k = &fullPress[k * gridsize];
    auto halfPress_k = &halfPress[k * gridsize];
    for (size_t i = 0; i < gridsize; ++i) halfPress_k[i] = 0.5 * (fullPress_km1[i] + fullPress_k[i]);
  }
  for (size_t i = 0; i < gridsize; ++i)
    halfPress[(numHalfLevels - 1) * gridsize + i] = fullPress[(numFullLevels - 1) * gridsize + i];
}

static void
calc_half_press(const Field3D &fullPress, Field3D &halfPress)
{
  if (fullPress.memType == MemType::Float)
    calc_half_press(fullPress.gridsize, fullPress.nlevels, fullPress.vec_f, halfPress.nlevels, halfPress.vec_f);
  else
    calc_half_press(fullPress.gridsize, fullPress.nlevels, fullPress.vec_d, halfPress.nlevels, halfPress.vec_d);
}

class Vertintap : public Process
{
  enum
  {
    func_pl,
    func_hl
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertintap",
    // clang-format off
    .operators = { { "ap2pl", func_pl, 0, "pressure levels in pascal", VertintapHelp },
                   { "ap2plx", func_pl, 0, "pressure levels in pascal" },
                   { "ap2hl", func_hl, 0, "height levels in meter" },
                   { "ap2hlx", func_hl, 0, "height levels in meter" } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertintap> registration = RegisterEntry<Vertintap>();

  int AP2PLX{}, AP2HLX{};
  int airPressID_FL = -1, airPressID_HL = -1, deltaPressID = -1;
  int psID = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  size_t gridsize{};

  int numPL{};
  int numVars{};
  int zaxisID_FL{};
  int zaxisID_HL{};
  int numFullLevels{};
  int numHalfLevels{};

  bool extrapolate{};

  VarList varList1{};
  VarList varList2{};

  Varray<double> pressureLevels{};

  std::vector<bool> processVars, interpVars;
  Varray2D<size_t> varnumMissVals;
  Field3DVector vardata1, vardata2;

  Varray<size_t> numMiss_FL, numMiss_HL;
  std::vector<int> vertIndex_FL, vertIndex_HL;
  Field psProg;
  Field3D fullPress, halfPress;

  CdoVar var3Dfull{}, var3Dhalf{};

public:
  void
  init() override
  {
    AP2PLX = module.get_id("ap2plx");
    AP2HLX = module.get_id("ap2hlx");

    auto operatorID = cdo_operator_id();
    auto useHeightLevel = (cdo_operator_f1(operatorID) == func_hl);

    extrapolate = (operatorID == AP2PLX || operatorID == AP2HLX);
    if (extrapolate == false) extrapolate = getenv_extrapolate();

    operator_input_arg(cdo_operator_enter(operatorID));

    if (cdo_operator_argc() == 1 && cdo_operator_argv(0) == "default")
    {
      if (useHeightLevel)
        pressureLevels = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000 };
      else
        pressureLevels
            = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
    }
    else { pressureLevels = cdo_argv_to_fltarr(cdo_get_oper_argv()); }

    numPL = pressureLevels.size();

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    gridsize = vlist_check_gridsize(vlistID1);

    auto zaxistype = useHeightLevel ? ZAXIS_HEIGHT : ZAXIS_PRESSURE;
    auto zaxisIDp = zaxisCreate(zaxistype, numPL);
    zaxisDefLevels(zaxisIDp, pressureLevels.data());

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);
    auto memtype = varList1.vars[0].memType;

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto stdname = string_to_lower(varList1.vars[varID].stdname);

      // clang-format off
      if      (stdname == var_stdname(surface_air_pressure))                psID = varID;
      else if (stdname == var_stdname(pressure_thickness))                  deltaPressID = varID;
      else if (stdname == var_stdname(air_pressure) && airPressID_FL == -1) airPressID_FL = varID;
      else if (stdname == var_stdname(air_pressure) && airPressID_HL == -1) airPressID_HL = varID;
      // clang-format on
    }

    if (-1 != airPressID_FL && -1 != airPressID_HL)
    {
      if (varList1.vars[airPressID_FL].nlevels == varList1.vars[airPressID_HL].nlevels)
        cdo_abort("Found two %s variables (%s/%s) with the same number of levels."
                  " Select one of them before using this operator!",
                  var_stdname(air_pressure), varList1.vars[airPressID_FL].name, varList1.vars[airPressID_HL].name);

      if (varList1.vars[airPressID_FL].nlevels == (varList1.vars[airPressID_HL].nlevels + 1))
        std::swap(airPressID_FL, airPressID_HL);

      if ((varList1.vars[airPressID_FL].nlevels + 1) != varList1.vars[airPressID_HL].nlevels)
        cdo_abort("Unexpected number of % levels in %s and %s. Number of half levels must be the number of full levels plus 1!",
                  var_stdname(air_pressure), varList1.vars[airPressID_FL].name, varList1.vars[airPressID_HL].name);
    }

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != psID)          cdo_print("  %s -> %s", var_stdname(surface_air_pressure), varList1.vars[psID].name);
      if (-1 != deltaPressID)  cdo_print("  %s -> %s", var_stdname(pressure_thickness), varList1.vars[deltaPressID].name);
      if (-1 != airPressID_FL) cdo_print("  %s (full) -> %s", var_stdname(air_pressure), varList1.vars[airPressID_FL].name);
      if (-1 != airPressID_HL) cdo_print("  %s (half) -> %s", var_stdname(air_pressure), varList1.vars[airPressID_HL].name);
      // clang-format on
    }

    if (-1 == airPressID_FL) cdo_abort("%s not found!", var_stdname(air_pressure));

    zaxisID_FL = (-1 == airPressID_FL) ? -1 : varList1.vars[airPressID_FL].zaxisID;
    zaxisID_HL = (-1 == airPressID_HL) ? -1 : varList1.vars[airPressID_HL].zaxisID;
    numFullLevels = (-1 == zaxisID_FL) ? 0 : varList1.vars[airPressID_FL].nlevels;
    numHalfLevels = (-1 == zaxisID_HL) ? numFullLevels + 1 : varList1.vars[airPressID_HL].nlevels;

    auto numZaxes = vlistNumZaxis(vlistID1);
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisID_FL || zaxisID == zaxisID_HL
          || (is_height_axis(zaxisID) && (nlevels == numHalfLevels || nlevels == numFullLevels)))
        vlistChangeZaxis(vlistID2, zaxisID, zaxisIDp);
    }

    varList2 = VarList(vlistID2);
    varList_set_memtype(varList2, memtype);

    processVars.resize(numVars);
    interpVars.resize(numVars);
    varnumMissVals.resize(numVars);
    vardata1.resize(numVars);
    vardata2.resize(numVars);

    auto maxLevels = std::max(std::max(numFullLevels, numHalfLevels), numPL);

    if (!extrapolate) numMiss_FL.resize(numPL);
    if (!extrapolate) numMiss_HL.resize(numPL);

    vertIndex_FL.resize(gridsize * numPL);
    vertIndex_HL.resize(gridsize * numPL);

    var3Dfull.gridsize = gridsize;
    var3Dfull.nlevels = numFullLevels;
    var3Dfull.memType = memtype;
    fullPress.init(var3Dfull);

    var3Dhalf.gridsize = gridsize;
    var3Dhalf.nlevels = numHalfLevels;
    var3Dhalf.memType = memtype;
    halfPress.init(var3Dhalf);

    if (useHeightLevel)
    {
      Varray<double> phlev(numPL);
      height_to_pressure(pressureLevels.data(), phlev.data(), numPL);

      if (Options::cdoVerbose)
        for (int i = 0; i < numPL; ++i) cdo_print("level=%d   height=%g   pressure=%g", i + 1, pressureLevels[i], phlev[i]);

      pressureLevels = phlev;
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto isHeightAxis = is_height_axis(var.zaxisID);

      if (gridInqType(var.gridID) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      vardata1[varID].init(var);

      interpVars[varID] = (var.zaxisID == zaxisID_FL
                           || (isHeightAxis && zaxisID_FL != -1 && (var.nlevels == numHalfLevels || var.nlevels == numFullLevels)));

      if (interpVars[varID])
      {
        varnumMissVals[varID].resize(maxLevels, 0);
        vardata2[varID].init(varList2.vars[varID]);
      }
      else
      {
        if (isHeightAxis && zaxisID_FL != -1 && var.nlevels > 1)
          cdo_warning("Parameter %d has wrong number of levels, skipped! (name=%s nlevel=%d)", varID + 1, var.name, var.nlevels);

        varnumMissVals[varID].resize(var.nlevels);
      }
    }

    if (zaxisID_FL != -1 && psID == -1)
    {
      if (deltaPressID != -1)
        cdo_warning("Surface pressure not found - set to vertical sum of %s!", var_stdname(pressure_thickness));
      else
        cdo_warning("Surface pressure not found - set to lower bound of %s!", var_stdname(air_pressure));
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (interpVars[varID] && varList1.vars[varID].isConstant) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

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

      for (int varID = 0; varID < numVars; ++varID)
      {
        processVars[varID] = false;
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID) varnumMissVals[varID][levelID] = 0;
      }

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, vardata1[varID], levelID, &varnumMissVals[varID][levelID]);
        processVars[varID] = true;
      }

      for (int varID = 0; varID < numVars; ++varID)
        if (interpVars[varID]) processVars[varID] = true;

      if (zaxisID_FL != -1)
      {
        if (tsID == 1 && varList1.vars[airPressID_FL].timeType == TIME_CONSTANT)
          cdo_warning("%s does not vary in time!", var_stdname(air_pressure));

        if (psID != -1)
        {
          psProg.init(varList1.vars[psID]);
          field_copy(vardata1[psID], psProg);
        }
        else if (deltaPressID != -1)
        {
          psProg.init(varList1.vars[deltaPressID]);
          field_fill(psProg, 0);
          for (int k = 0; k < numFullLevels; ++k) field_add(psProg, vardata1[deltaPressID], k);
        }
        else
        {
          psProg.init(varList1.vars[airPressID_FL]);
          field_copy(vardata1[airPressID_FL], numFullLevels - 1, psProg);
        }

        // check range of psProg
        check_range_ps(tsID + 1, psProg);

        field_copy(vardata1[airPressID_FL], fullPress);

        if (-1 != zaxisID_HL)
          field_copy(vardata1[airPressID_HL], halfPress);
        else
          calc_half_press(fullPress, halfPress);

        gen_vert_index(vertIndex_FL, pressureLevels, fullPress, gridsize);
        if (!extrapolate) gen_vert_index_mv(vertIndex_FL, pressureLevels, gridsize, psProg, numMiss_FL);

        gen_vert_index(vertIndex_HL, pressureLevels, halfPress, gridsize);
        if (!extrapolate) gen_vert_index_mv(vertIndex_HL, pressureLevels, gridsize, psProg, numMiss_HL);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        if (processVars[varID])
        {
          auto const &var = varList1.vars[varID];

          if (tsID > 0 && !interpVars[varID] && var.isConstant) continue;

          if (interpVars[varID])
          {
            if (var.nlevels != numFullLevels && var.nlevels != numHalfLevels)
              cdo_abort("Number of generalized height level differ from full/half level (param=%s)!", var.name);

            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
            }

            auto const &levels3D = (var.nlevels == numFullLevels) ? fullPress : halfPress;
            auto const &vertIndex3D = (var.nlevels == numFullLevels) ? vertIndex_FL : vertIndex_HL;
            vertical_interp_X(levels3D, vardata1[varID], vardata2[varID], vertIndex3D, pressureLevels, gridsize);

            if (!extrapolate)
            {
              auto const &numMiss = (var.nlevels == numFullLevels) ? numMiss_FL : numMiss_HL;
              varray_copy(numPL, numMiss, varnumMissVals[varID]);
            }
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
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
