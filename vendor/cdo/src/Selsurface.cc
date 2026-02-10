/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "interpol.h"
#include "field_functions.h"

template <typename T>
static void
isosurface_kernel(double isoval, size_t numMissVals, Varray<double> const &levels, int nlevels, size_t gridsize, double mv,
                  const Varray<const T *> &data3D, Varray<T> &data2D)
{
  T missval = mv;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    data2D[i] = missval;

    for (int k = 0; k < (nlevels - 1); ++k)
    {
      double val1 = data3D[k][i];
      double val2 = data3D[k + 1][i];

      if (numMissVals)
      {
        auto hasMissvals1 = fp_is_equal(val1, missval);
        auto hasMissvals2 = fp_is_equal(val2, missval);
        if (hasMissvals1 && hasMissvals2) continue;
        if (hasMissvals1 && is_equal(isoval, val2)) data2D[i] = levels[k + 1];
        if (hasMissvals2 && is_equal(isoval, val1)) data2D[i] = levels[k];
        if (hasMissvals1 || hasMissvals2) continue;
      }

      if ((isoval >= val1 && isoval <= val2) || (isoval >= val2 && isoval <= val1))
      {
        data2D[i] = is_equal(val1, val2) ? levels[k] : intlin(isoval, levels[k], val1, levels[k + 1], val2);
        break;
      }
    }
  }
}

static void
isosurface(double isoval, int nlevels, Varray<double> const &levels, FieldVector const &field3D, Field &field2D)
{
  auto gridsize = gridInqSize(field3D[0].grid);
  auto missval = field3D[0].missval;

  auto numMissVals = field3D[0].numMissVals;
  for (int k = 1; k < nlevels; ++k) numMissVals += field3D[k].numMissVals;
  /*
  auto func = [&](auto &v)
  {
    Varray<std::decay_t<decltype(v[0])> const *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_f.data();
    isosurface_kernel(isoval, numMissVals, levels, nlevels, gridsize, missval, data3D, v);
  };
  field_operation(func, field2D);
  */
  if (field3D[0].memType == MemType::Float)
  {
    Varray<const float *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_f.data();
    isosurface_kernel(isoval, numMissVals, levels, nlevels, gridsize, missval, data3D, field2D.vec_f);
  }
  else
  {
    Varray<const double *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_d.data();
    isosurface_kernel(isoval, numMissVals, levels, nlevels, gridsize, missval, data3D, field2D.vec_d);
  }

  field_num_mv(field2D);
}

template <typename T>
static void
layer_value_min_kernel(int nlevels, size_t gridsize, double mv, const Varray<const T *> &data3D, Varray<T> &data2D)
{
  T missval = mv;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    data2D[i] = missval;

    for (int k = 0; k < nlevels; ++k)
    {
      auto val = data3D[k][i];
      if (fp_is_not_equal(val, missval))
      {
        data2D[i] = val;
        break;
      }
    }
  }
}

static void
layer_value_min(int nlevels, FieldVector const &field3D, Field &field2D)
{
  auto gridsize = gridInqSize(field3D[0].grid);
  auto missval = field3D[0].missval;

  if (field3D[0].memType == MemType::Float)
  {
    Varray<const float *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_f.data();
    layer_value_min_kernel(nlevels, gridsize, missval, data3D, field2D.vec_f);
  }
  else
  {
    Varray<const double *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_d.data();
    layer_value_min_kernel(nlevels, gridsize, missval, data3D, field2D.vec_d);
  }

  field_num_mv(field2D);
}

template <typename T>
static void
layer_value_max_kernel(int nlevels, size_t gridsize, double mv, const Varray<const T *> &data3D, Varray<T> &data2D)
{
  T missval = mv;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    data2D[i] = missval;

    for (int k = nlevels - 1; k >= 0; --k)
    {
      auto val = data3D[k][i];
      if (fp_is_not_equal(val, missval))
      {
        data2D[i] = val;
        break;
      }
    }
  }
}

static void
layer_value_max(int nlevels, FieldVector const &field3D, Field &field2D)
{
  auto gridsize = gridInqSize(field3D[0].grid);
  auto missval = field3D[0].missval;

  if (field3D[0].memType == MemType::Float)
  {
    Varray<const float *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_f.data();
    layer_value_max_kernel(nlevels, gridsize, missval, data3D, field2D.vec_f);
  }
  else
  {
    Varray<const double *> data3D(nlevels);
    for (int k = 0; k < nlevels; ++k) data3D[k] = field3D[k].vec_d.data();
    layer_value_max_kernel(nlevels, gridsize, missval, data3D, field2D.vec_d);
  }

  field_num_mv(field2D);
}

class Selsurface : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selsurface",
    // clang-format off
    .operators = { { "isosurface", SelsurfaceHelp },
                   { "bottomvalue", SelsurfaceHelp },
                   { "topvalue", SelsurfaceHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Selsurface> registration = RegisterEntry<Selsurface>(module);

private:
  int ISOSURFACE{}, BOTTOMVALUE{}, TOPVALUE{};
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int zaxisID1 = -1;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int numVars{};
  std::vector<bool> isVar3D;
  std::vector<bool> foundVar;
  int operatorID{};

  double isoval = 0.0;

  VarList varList1;
  Field field2;
  FieldVector2D varsData1;

  int numLevels{};
  Varray<double> levels;

  bool isPositive{};
  bool isReverse{};

public:
  void
  init() override
  {
    ISOSURFACE = module.get_id("isosurface");
    BOTTOMVALUE = module.get_id("bottomvalue");
    TOPVALUE = module.get_id("topvalue");

    operatorID = cdo_operator_id();

    if (operatorID == ISOSURFACE)
    {
      operator_input_arg("isoval");
      operator_check_argc(1);
      isoval = parameter_to_double(cdo_operator_argv(0));
    }

    if (Options::cdoVerbose) cdo_print("Isoval: %g", isoval);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numZaxes = varList1.numZaxes();
    for (int i = 0; i < numZaxes; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevels_zaxis = zaxisInqSize(zaxisID);
      if (zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF)
        if (nlevels_zaxis > 1)
        {
          zaxisID1 = zaxisID;
          break;
        }
    }
    if (zaxisID1 == -1) cdo_abort("No processable variable found!");

    numLevels = zaxisInqSize(zaxisID1);
    levels = Varray<double>(numLevels);
    cdo_zaxis_inq_levels(zaxisID1, levels.data());

    isPositive = !(levels[0] < 0.0 && levels[numLevels - 1] < 0.0);
    isReverse = (levels[0] > levels[numLevels - 1]);

    auto zaxisIDsfc = zaxis_from_name("surface");
    for (int i = 0; i < numZaxes; ++i)
      if (zaxisID1 == vlistZaxis(vlistID1, i)) vlistChangeZaxisIndex(vlistID2, i, zaxisIDsfc);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    numVars = varList1.numVars();
    isVar3D = std::vector<bool>(numVars);
    foundVar = std::vector<bool>(numVars);

    field2D_init(varsData1, varList1);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      isVar3D[varID] = (var.zaxisID == zaxisID1);
    }
  }

  void
  run() override
  {
    auto bottom_value_func = isReverse ? layer_value_max : layer_value_min;
    auto top_value_func = isReverse ? layer_value_min : layer_value_max;
    if (isPositive && positive_is_down(zaxisID1)) std::swap(bottom_value_func, top_value_func);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID) foundVar[varID] = false;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];
        auto &field1 = varsData1[varID][levelID];
        field1.init(var);
        cdo_read_field(streamID1, field1);
        foundVar[varID] = true;
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        if (foundVar[varID])
        {
          auto const &var = varList1.vars[varID];
          if (isVar3D[varID])
          {
            field2.init(var);
            // clang-format off
            if      (operatorID == ISOSURFACE)  isosurface(isoval, numLevels, levels, varsData1[varID], field2);
            else if (operatorID == BOTTOMVALUE) bottom_value_func(numLevels, varsData1[varID], field2);
            else if (operatorID == TOPVALUE)    top_value_func(numLevels, varsData1[varID], field2);
            // clang-format on

            cdo_def_field(streamID2, varID, 0);
            cdo_write_field(streamID2, field2);
          }
          else
          {
            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              cdo_def_field(streamID2, varID, levelID);
              cdo_write_field(streamID2, varsData1[varID][levelID]);
            }
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

    vlistDestroy(vlistID2);
  }
};
