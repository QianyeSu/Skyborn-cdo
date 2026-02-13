/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intlevel   intlevel        Linear level interpolation
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_omp.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "pmlist.h"
#include "param_conversion.h"

template <typename T>
static double
vert_interp_lev_kernel(float w1, float w2, T var1L1, T var1L2, T missval)
{
  if (fp_is_equal(var1L1, missval)) w1 = 0.0f;
  if (fp_is_equal(var1L2, missval)) w2 = 0.0f;

  // clang-format off
  if      (is_equal(w1, 0.0f) && is_equal(w2, 0.0f)) return missval;
  else if (is_equal(w1, 0.0f)) return (w2 >= 0.5f) ? var1L2 : missval;
  else if (is_equal(w2, 0.0f)) return (w1 >= 0.5f) ? var1L1 : missval;
  else                         return var1L1 * (double)w1 + var1L2 * (double)w2;
  // clang-format on
}

constexpr int BottomLevel = 32000;
constexpr int TopLevel = 32001;

static void
restore_index_and_weights(int nlev1, int idx, float wgt, int &idx1, int &idx2, float &wgt1, float &wgt2)
{
  if (idx == BottomLevel)
  {
    idx1 = 0;
    idx2 = 0;
    wgt1 = 0.0f;
    wgt2 = wgt;
  }
  else if (idx == TopLevel)
  {
    idx1 = nlev1 - 1;
    idx2 = nlev1 - 1;
    wgt1 = wgt;
    wgt2 = 0.0f;
  }
  else
  {
    idx1 = (idx < 0) ? -idx : idx;
    idx2 = (idx < 0) ? idx1 - 1 : idx1 + 1;
    wgt1 = wgt;
    wgt2 = 1.0f - wgt;
  }
  // printf("%d %d %g %g\n", idx1, idx2, wgt1, wgt2);
}

//  1D vertical interpolation
template <typename T1, typename T2>
static void
vert_interp_lev(size_t gridsize, int nlev1, double mv, Varray<T1> const &vardata1, Varray<T2> &vardata2, int nlev2,
                Varray<int> const &lev_idx, Varray<float> const &lev_wgt)
{
  T1 missval = mv;

  for (int ilev = 0; ilev < nlev2; ++ilev)
  {
    auto idx = lev_idx[ilev];
    auto wgt = lev_wgt[ilev];
    int idx1, idx2;
    float wgt1, wgt2;
    restore_index_and_weights(nlev1, idx, wgt, idx1, idx2, wgt1, wgt2);

    // upper/lower values from input field
    auto var1L1 = &vardata1[gridsize * idx1];
    auto var1L2 = &vardata1[gridsize * idx2];

    auto var2 = &vardata2[gridsize * ilev];

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (size_t i = 0; i < gridsize; ++i) { var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1[i], var1L2[i], missval); }
  }
}

static void
vert_interp_lev(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                Varray<int> const &lev_idx, Varray<float> const &lev_wgt)
{
  auto func = [&](auto &v1, auto &v2) { vert_interp_lev(gridsize, nlev1, missval, v1, v2, nlev2, lev_idx, lev_wgt); };
  field_operation2(func, field1, field2);
}

//  3D vertical interpolation
template <typename T1, typename T2>
void
vert_interp_lev3d(size_t gridsize, int nlev1, double mv, Varray<T1> const &vardata1, Varray<T2> &vardata2, int nlev2,
                  Varray<int> const &lev_idx, Varray<float> const &lev_wgt)
{
  T1 missval = mv;

  for (int ilev = 0; ilev < nlev2; ilev++)
  {
    auto offset = ilev * gridsize;
    auto var2 = &vardata2[offset];

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (size_t i = 0; i < gridsize; ++i)
    {
      auto idx = lev_idx[offset + i];
      auto wgt = lev_wgt[offset + i];
      int idx1, idx2;
      float wgt1, wgt2;
      restore_index_and_weights(nlev1, idx, wgt, idx1, idx2, wgt1, wgt2);

      // upper/lower values from input field
      auto var1L1 = vardata1[idx1 * gridsize + i];
      auto var1L2 = vardata1[idx2 * gridsize + i];

      var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1, var1L2, missval);
    }
  }
}

void
vert_interp_lev3d(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                  Varray<int> const &lev_idx, Varray<float> const &lev_wgt)
{
  auto func = [&](auto &v1, auto &v2) { vert_interp_lev3d(gridsize, nlev1, missval, v1, v2, nlev2, lev_idx, lev_wgt); };
  field_operation2(func, field1, field2);
}

void
vert_gen_weights(int extrapolate, int nlev1, Varray<double> const &lev1, int nlev2, Varray<double> const &lev2,
                 Varray<int> &lev_idx, Varray<float> &lev_wgt)
{
  for (int i2 = 0; i2 < nlev2; ++i2)
  {
    int idx1 = 0, idx2 = 0;
    double val2 = 0.0;

    // Because 2 levels were added to the source vertical coordinate (one on top, one at the bottom), its loop starts at 1
    int i1;
    for (i1 = 1; i1 < nlev1; ++i1)
    {
      auto lev1_isUp = (lev1[i1 - 1] < lev1[i1]);
      idx1 = lev1_isUp ? i1 - 1 : i1;
      idx2 = lev1_isUp ? i1 : i1 - 1;
      auto val1 = lev1[idx1];
      val2 = lev1[idx2];
      if (lev2[i2] > val1 && lev2[i2] <= val2) break;
    }

    if (i1 == nlev1) cdo_abort("Level %g not found!", lev2[i2]);

    if (i1 - 1 == 0)  // destination levels is not covert by the first two input z levels
    {
      lev_idx[i2] = BottomLevel;
      lev_wgt[i2] = static_cast<float>(extrapolate || is_equal(lev2[i2], val2));
    }
    else if (i1 == nlev1 - 1)  // destination level is beyond the last value of the input z field
    {
      lev_idx[i2] = TopLevel;
      lev_wgt[i2] = static_cast<float>(extrapolate || is_equal(lev2[i2], val2));
    }
    else  // target z values has two bounday values in input z field
    {
      lev_idx[i2] = idx1 - 1;
      if (idx1 > idx2) lev_idx[i2] = -lev_idx[i2];
      lev_wgt[i2] = (lev1[idx2] - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
    }
    // printf("%d %g %d %d %g %g %d %g\n", i2, lev2[i2], idx1, idx2, lev1[idx1], lev1[idx2], lev_idx[i2], lev_wgt[i2]);
  }
}

bool
levelDirUp(int nlev, const double *const lev)
{
  auto lup = (nlev > 1 && lev[1] > lev[0]);
  if (lup)
  {
    for (int k = 1; k < nlev - 1; ++k)
      if (lev[k + 1] <= lev[k]) return false;
  }

  return lup;
}

bool
levelDirDown(int nlev, const double *const lev)
{
  auto ldown = (nlev > 1 && lev[1] < lev[0]);
  if (ldown)
  {
    for (int k = 1; k < nlev - 1; ++k)
      if (lev[k + 1] >= lev[k]) return false;
  }

  return ldown;
}

template <typename T>
static void
vert_gen_weights3d1d(bool extrapolate, size_t gridsize, int nlev1, Varray<T> const &xlev1, int nlev2, Varray<double> const &lev2,
                     Varray<int> &xlev_idx, Varray<float> &xlev_wgt)
{
  auto nthreads = Threading::ompNumMaxThreads;
  Varray2D<double> lev1p2(nthreads, Varray<double>(nlev1 + 2));
  Varray2D<float> lev_wgt(nthreads, Varray<float>(nlev2));
  Varray2D<int> lev_idx(nthreads, Varray<int>(nlev2));

  // Check monotony of vertical levels
  for (int k = 0; k < nlev1; ++k) lev1p2[0][k] = xlev1[k * gridsize];
  auto lup = levelDirUp(nlev1, lev1p2[0].data());
  auto ldown = levelDirDown(nlev1, lev1p2[0].data());
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  double level_0 = lup ? -1.e33 : 1.e33;
  double level_N = lup ? 1.e33 : -1.e33;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto ompthID = cdo_omp_get_thread_num();

    lev1p2[ompthID][0] = level_0;
    lev1p2[ompthID][nlev1 + 1] = level_N;
    for (int k = 0; k < nlev1; ++k) lev1p2[ompthID][k + 1] = xlev1[k * gridsize + i];

    vert_gen_weights(extrapolate, nlev1 + 2, lev1p2[ompthID], nlev2, lev2, lev_idx[ompthID], lev_wgt[ompthID]);

    for (int k = 0; k < nlev2; ++k) xlev_idx[k * gridsize + i] = lev_idx[ompthID][k];
    for (int k = 0; k < nlev2; ++k) xlev_wgt[k * gridsize + i] = lev_wgt[ompthID][k];
  }
}

static void
vert_gen_weights3d1d(bool extrapolate, size_t gridsize, int nlev1, Field3D &field1, int nlev2, Varray<double> const &lev2,
                     Varray<int> &lev_idx, Varray<float> &lev_wgt)
{
  auto func = [&](auto &v) { vert_gen_weights3d1d(extrapolate, gridsize, nlev1, v, nlev2, lev2, lev_idx, lev_wgt); };
  field_operation(func, field1);
}

static int
create_zaxis_from_zvar(Varray<double> const &levels, int vlistID, int varID)
{
  auto nlevels = levels.size();
  auto zaxisID = zaxisCreate(ZAXIS_GENERIC, nlevels);

  std::string name = "zlev";
  cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name.c_str());
  auto longname = cdo::inq_var_longname(vlistID, varID);
  if (longname.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname.c_str());
  auto units = cdo::inq_var_units(vlistID, varID);
  if (units.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units.c_str());

  zaxisDefLevels(zaxisID, levels.data());

  return zaxisID;
}

static int
create_zaxis_from_zaxis(Varray<double> const &levels, int zaxisID1)
{
  auto nlevels = levels.size();
  auto zaxisID2 = zaxisCreate(zaxisInqType(zaxisID1), nlevels);

  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_NAME, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_LONGNAME, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_UNITS, zaxisID2);
  cdiCopyKey(zaxisID1, CDI_GLOBAL, CDI_KEY_DATATYPE, zaxisID2);

  zaxisDefLevels(zaxisID2, levels.data());

  return zaxisID2;
}

namespace
{
struct Parameter
{
  Varray<double> levels;
  std::string zdescription;
  std::string zvarname;
  bool extrapolate{ false };
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      // if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &values = kv.values;
      auto const &value = kv.values[0];
      // int numValues = kv.nvalues;
      // if (numValues == 1 && value.empty()) numValues = 0;

      if (key == "level")
      {
        if (!params.zdescription.empty()) cdo_abort("Parameter level and zdescription can't be mixed!");
        params.levels = cdo_argv_to_fltarr(values);
      }
      else if (key == "zdescription")
      {
        if (!params.levels.empty()) cdo_abort("Parameter zdescription and level can't be mixed!");
        params.zdescription = value;
      }
      else if (key == "zvarname") { params.zvarname = value; }
      else if (key == "extrapolate") { params.extrapolate = parameter_to_bool(value); }
      else { cdo_abort("Invalid parameter key >%s<!", key); }
    }
  }

  return params;
}

static void
handle_zvar(size_t &zvarGridsize, size_t &wisize, std::string const &zvarname, VarList const &varList1, int &zvarID,
            bool &zvarIsVarying, Varray<double> const &lev2, int &nlev1, int const &nlev2, int const &vlistID1, int &zaxisID2,
            int &zaxisID1)
{
  for (auto const &var : varList1.vars)
  {
    if (zvarname == var.name)
    {
      zvarID = var.ID;
      break;
    }
  }

  auto const &zvar = varList1.vars[zvarID];
  if (zvarID == CDI_UNDEFID) cdo_abort("Variable %s not found!", zvarname);
  zvarIsVarying = (zvar.timeType == TIME_VARYING);
  zvarGridsize = zvar.gridsize;
  nlev1 = zvar.nlevels;

  if (zaxisID2 == CDI_UNDEFID) zaxisID2 = create_zaxis_from_zvar(lev2, vlistID1, zvarID);

  wisize = zvarGridsize * nlev2;

  auto numZaxes = vlistNumZaxis(vlistID1);
  int i = 0;
  for (; i < numZaxes; ++i)
  {
    auto zaxisID = vlistZaxis(vlistID1, i);
    auto nlevels = zaxisInqSize(zaxisID);
    if (nlevels == nlev1)
    {
      zaxisID1 = zaxisID;
      break;
    }
  }
  if (i == numZaxes) cdo_abort("No processable variable found!");
}

static void
handle_empty_zvar(size_t &wisize, Varray<double> &lev1, Varray<double> const &lev2, int &nlev1, int const &nlev2,
                  int const &vlistID1, int &zaxisID2, int &zaxisID1)
{
  int numLevels = 0;
  auto numZaxes = vlistNumZaxis(vlistID1);
  int i = 0;
  for (; i < numZaxes; ++i)
  {
    auto zaxisID = vlistZaxis(vlistID1, i);
    numLevels = zaxisInqSize(zaxisID);
    // if (zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF)
    if (numLevels > 1)
    {
      zaxisID1 = zaxisID;
      break;
    }
  }
  if (i == numZaxes) cdo_abort("No processable variable found!");

  if (zaxisID2 == CDI_UNDEFID) zaxisID2 = create_zaxis_from_zaxis(lev2, zaxisID1);

  nlev1 = numLevels;
  lev1.resize(nlev1 + 2);
  cdo_zaxis_inq_levels(zaxisID1, &lev1[1]);

  auto lup = levelDirUp(nlev1, &lev1[1]);
  auto ldown = levelDirDown(nlev1, &lev1[1]);
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  lev1[0] = lup ? -1.e33 : 1.e33;
  lev1[nlev1 + 1] = lup ? 1.e33 : -1.e33;

  if (Options::cdoVerbose)
    for (i = 0; i < nlev1 + 2; ++i) cdo_print("level_in %d: %g", i, lev1[i]);

  wisize = nlev2;
}

class Intlevel : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intlevel",
    // clang-format off
    .operators = { { "intlevel", IntlevelHelp },
                   { "intlevelx", IntlevelHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Intlevel> registration = RegisterEntry<Intlevel>();

  int zaxisID1 = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  Varray<double> lev2;
  int zaxisID2 = CDI_UNDEFID;

  int numLevels1 = 0;
  bool zvarIsVarying = false;

  Varray<double> lev1;
  int numLevels2{};

  size_t zvarGridsize = 0;
  int zvarID = CDI_UNDEFID;

  MemType memType{};

  Varray<int> lev_idx;
  Varray<float> lev_wgt;

  VarList varList1;
  VarList varList2;

  Parameter params{};

public:
  void
  init() override
  {
    operator_input_arg("level|zdescription[, zvarname, extrapolate]");

    zaxisID2 = CDI_UNDEFID;
    auto const &argList = cdo_get_oper_argv();

    if (std::isdigit((int) argList[0][0])) { params.levels = cdo_argv_to_fltarr(argList); }
    else
    {
      params = get_parameter();
      if (!params.zdescription.empty())
      {
        auto zfilename = params.zdescription;
        auto fobj = c_fopen(zfilename, "r");
        if (!fobj.get()) cdo_abort("Open failed on %s", zfilename);
        zaxisID2 = zaxis_from_file(fobj.get(), zfilename.c_str());
        if (zaxisID2 == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zfilename);
        auto nlevels = zaxisInqSize(zaxisID2);
        params.levels.resize(nlevels);
        zaxisInqLevels(zaxisID2, params.levels.data());
      }
    }

    lev2 = params.levels;
    numLevels2 = lev2.size();

    if (Options::cdoVerbose)
      for (int i = 0; i < numLevels2; ++i) cdo_print("level out %d: %g", i, lev2[i]);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);
    memType = varList1.vars[0].memType;

    // Find z-variable
    size_t wisize = 0;

    if (params.zvarname.empty()) { handle_empty_zvar(wisize, lev1, lev2, numLevels1, numLevels2, vlistID1, zaxisID2, zaxisID1); }
    else
    {
      handle_zvar(zvarGridsize, wisize, params.zvarname, varList1, zvarID, zvarIsVarying, lev2, numLevels1, numLevels2, vlistID1,
                  zaxisID2, zaxisID1);
    }

    auto numZaxes = vlistNumZaxis(vlistID1);
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto numLevels = zaxisInqSize(zaxisID);
      if (zaxisID == zaxisID1 || numLevels == numLevels1) vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
    }

    varList2 = VarList(vlistID2);
    varList_set_memtype(varList2, memType);

    lev_idx.resize(wisize);
    lev_wgt.resize(wisize);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto numVars = varList1.numVars();
    std::vector<bool> processVars(numVars);
    std::vector<bool> interpVars(numVars, false);
    std::vector<std::vector<size_t>> varnumMissVals(numVars);
    Field3DVector vardata1(numVars);
    Field3DVector vardata2(numVars);

    int maxLevels = std::max(numLevels1, numLevels2);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      vardata1[varID].init(var1);
      interpVars[varID] = (var1.zaxisID == zaxisID1 || var1.nlevels == numLevels1);

      if (interpVars[varID])
      {
        varnumMissVals[varID].resize(maxLevels, 0);
        vardata2[varID].init(varList2.vars[varID]);
      }
      else { varnumMissVals[varID].resize(var1.nlevels); }
    }

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

      if (tsID == 0 || zvarIsVarying)
      {
        if (!params.zvarname.empty())
          vert_gen_weights3d1d(params.extrapolate, zvarGridsize, numLevels1, vardata1[zvarID], numLevels2, lev2, lev_idx, lev_wgt);
        else
          vert_gen_weights(params.extrapolate, numLevels1 + 2, lev1, numLevels2, lev2, lev_idx, lev_wgt);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        if (processVars[varID] && interpVars[varID])
        {
          auto const &var1 = varList1.vars[varID];

          auto missval = var1.missval;
          auto gridsize = var1.gridsize;

          if (!params.zvarname.empty())
            vert_interp_lev3d(gridsize, numLevels1, missval, vardata1[varID], vardata2[varID], numLevels2, lev_idx, lev_wgt);
          else
            vert_interp_lev(gridsize, numLevels1, missval, vardata1[varID], vardata2[varID], numLevels2, lev_idx, lev_wgt);

          for (int levelID = 0; levelID < numLevels2; ++levelID)
          {
            auto offset = gridsize * levelID;
            auto func = [&](auto const &v) { varnumMissVals[varID][levelID] = array_num_mv(gridsize, &v[offset], missval); };
            field_operation(func, vardata2[varID]);
          }
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        if (processVars[varID])
        {
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
