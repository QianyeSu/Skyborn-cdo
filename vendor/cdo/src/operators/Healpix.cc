/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "workerthread.h"
#include "process_int.h"
#include "cdo_math.h"
#include "cdo_options.h"
#include "cdo_omp.h"
#include "grid_healpix.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "mpim_grid.h"
#include "progress.h"

enum struct Stat
{
  Mean = 1,
  Avg = 2
};

namespace
{
struct Parameter
{
  int fact{ 1 };
  int nsideIn{ 0 };
  int nsideOut{ 0 };
  int zoom{ 0 };
  HpOrder orderIn{ HpOrder::Undef };
  HpOrder orderOut{ HpOrder::Undef };
  Stat stat{ Stat::Mean };
  double power{ 0.0 };
  bool doDegrade{ true };
};
}  // namespace

template <typename T>
static std::pair<double, size_t>
stat_mean_mv_kernel(T const *const v, size_t n, T missval)
{
  double sum = 0.0;
  size_t nOut = 0;
  for (size_t i = 0; i < n; ++i)
    if (fp_is_not_equal(v[i], missval))
    {
      sum += v[i];
      nOut++;
    }

  return std::make_pair(sum, nOut);
}

template <typename T>
static T
stat_avg_mv(T const *const v, size_t n, T missval, double scale)
{
  auto [sum, nOut] = stat_mean_mv_kernel(v, n, missval);
  return (nOut == n) ? (sum / nOut) * scale : missval;
}

template <typename T>
static T
stat_mean_mv(T const *const v, size_t n, T missval, double scale)
{
  auto [sum, nOut] = stat_mean_mv_kernel(v, n, missval);
  return (nOut > 0) ? (sum / nOut) * scale : missval;
}

template <typename T>
static T
stat_mean(T const *const v, size_t n)
{
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i) { sum += v[i]; }
  return sum / n;
}

static double
get_scalefactor(Parameter const &params)
{
  double nsideQuot = static_cast<double>(params.nsideIn) / params.nsideOut;
  return (std::fabs(params.power) > 0.0) ? std::pow(nsideQuot, -params.power) : 1.0;
}

template <typename T1, typename T2>
static void
degrade(Varray<T1> const &v1, size_t n2, Varray<T2> &v2, bool hasMissvals, double mv, Parameter const &params)
{
  auto scale = get_scalefactor(params);
  size_t nvals = params.fact * params.fact;
  size_t nx = n2 * nvals;

  if (hasMissvals)
  {
    T1 missval = mv;
    if (params.stat == Stat::Mean)
    {
#ifdef _OPENMP
#pragma omp parallel for if (nx > 4 * cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < n2; ++i) { v2[i] = stat_mean_mv(&v1[i * nvals], nvals, missval, scale); }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for if (nx > 4 * cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < n2; ++i) { v2[i] = stat_avg_mv(&v1[i * nvals], nvals, missval, scale); }
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (nx > 4 * cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n2; ++i) { v2[i] = stat_mean(&v1[i * nvals], nvals) * scale; }
  }
}

static void
hp_degrade(Field const &field1, Field &field2, Parameter const &params)
{
  auto hasMissvals = (field1.numMissVals > 0);
  auto func = [&](auto const &v1, auto &v2, auto n2, double mv) { degrade(v1, n2, v2, hasMissvals, mv, params); };
  field_operation2(func, field1, field2, field2.gridsize, field1.missval);
  if (hasMissvals) field_num_mv(field2);
}

template <typename T1, typename T2>
static void
upgrade(size_t n1, Varray<T1> const &v1, Varray<T2> &v2, bool hasMissvals, double mv, Parameter const &params)
{
  auto scale = get_scalefactor(params);
  size_t nvals = params.fact * params.fact;
  if (hasMissvals)
  {
    T1 missval = mv;
    for (size_t i = 0; i < n1; ++i)
    {
      for (size_t k = 0; k < nvals; ++k) { v2[i * nvals + k] = fp_is_equal(v1[i], missval) ? missval : v1[i] * scale; }
    }
  }
  else
  {
    for (size_t i = 0; i < n1; ++i)
    {
      for (size_t k = 0; k < nvals; ++k) { v2[i * nvals + k] = v1[i] * scale; }
    }
  }
}

static void
hp_upgrade(Field const &field1, Field &field2, Parameter const &params)
{
  auto hasMissvals = (field1.numMissVals > 0);
  auto func = [&](auto const &v1, auto &v2, auto n1, double mv) { upgrade(n1, v1, v2, hasMissvals, mv, params); };
  field_operation2(func, field1, field2, field1.gridsize, field1.missval);
  if (hasMissvals) field_num_mv(field2);
}

template <typename T>
static void
ring_to_nested(int nside, size_t gridsize, Varray<T> &v)
{
  Varray<T> vtmp = v;
  hp_ring_to_nested(nside, gridsize, vtmp.data(), v.data());
}

static void
ring_to_nested(Field &field, int nside)
{
  auto func = [&](auto &v, auto gridsize) { ring_to_nested(nside, gridsize, v); };
  field_operation(func, field, field.gridsize);
}

template <typename T>
static void
nested_to_ring(int nside, size_t gridsize, Varray<T> &v)
{
  Varray<T> vtmp = v;
  hp_nested_to_ring(nside, gridsize, vtmp.data(), v.data());
}

static void
nested_to_ring(Field &field, int nside)
{
  auto func = [&](auto &v, auto gridsize) { nested_to_ring(nside, gridsize, v); };
  field_operation(func, field, field.gridsize);
}

static Stat
set_stat(std::string const &statString)
{
  if (statString == "mean") { return Stat::Mean; }
  if (statString == "avg") { return Stat::Avg; }

  cdo_abort("Parameter value stat=%s unsupported!", statString);
  return Stat::Mean;
}

static Parameter
get_parameter(void)
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
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "nside") params.nsideOut = parameter_to_int(value);
      else if (key == "zoom")  params.zoom     = parameter_to_int(value);
      else if (key == "order") params.orderOut = hp_get_order(parameter_to_word(value));
      else if (key == "fact")  params.fact     = parameter_to_int(value);
      else if (key == "stat")  params.stat     = set_stat(parameter_to_word(value));
      else if (key == "power") params.power    = parameter_to_double(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
verify_parameter(Parameter const &params)
{
  if (params.fact > 1 && params.nsideOut > 0) cdo_abort("Parameter 'fact' can't be combined with 'nside'!");
  if (params.fact > 1 && params.zoom > 0) cdo_abort("Parameter 'fact' can't be combined with 'zoom'!");
  if (params.zoom > 0 && params.nsideOut > 0) cdo_abort("Parameter 'zoom' can't be combined with 'nside'!");
}

static int
define_healpix_grid(size_t gridsize, int nside, HpOrder order)
{
  int refinementLevel = static_cast<int>(std::log2(nside));
  auto orderString = (order == HpOrder::Ring) ? "ring" : "nested";

  auto gridID = gridCreate(GRID_HEALPIX, gridsize);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_DIMNAME, "cell");
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "crs");
  const std::string gridmapName = "healpix";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gridmapName.c_str());
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) gridmapName.size(), gridmapName.c_str());
  cdiDefAttInt(gridID, CDI_GLOBAL, "refinement_level", CDI_DATATYPE_INT32, 1, &refinementLevel);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "indexing_scheme", (int) std::strlen(orderString), orderString);

  if (refinementLevel < 14) cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_INT32);

  return gridID;
}

static int
define_healpix_proj(size_t gridsize, int nside, HpOrder order)
{
  auto orderString = (order == HpOrder::Ring) ? "ring" : "nested";
  auto projection = "healpix";
  auto gridID = gridCreate(GRID_PROJECTION, gridsize);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_DIMNAME, "cells");
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, projection);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, projection);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) std::strlen(projection), projection);
  cdiDefAttInt(gridID, CDI_GLOBAL, "healpix_nside", CDI_DATATYPE_INT32, 1, &nside);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "healpix_order", (int) std::strlen(orderString), orderString);

  return gridID;
}

static int
hp_define_grid(int gridID1, Parameter &params)
{
  auto gridType = gridInqType(gridID1);

  auto hpParams = cdo::get_healpix_params(gridID1);
  auto nside = hpParams.nside();
  auto order = hpParams.order();
  params.nsideIn = nside;
  params.orderIn = order;

  if (!cdo::is_power_of_two(params.nsideIn)) cdo_abort("Input healpix: nside must be a power of two!");

  if (params.nsideOut == 0)
  {
    auto fact = params.fact;
    params.nsideOut = (fact > 1) ? (params.doDegrade ? nside / fact : nside * fact) : nside;
  }
  else
  {
    if (params.doDegrade)
    {
      if (params.nsideOut > params.nsideIn)
        cdo_abort("Parameter nside=%d must be less than input nside=%d!", params.nsideOut, params.nsideIn);
      params.fact = params.nsideIn / params.nsideOut;
    }
    else
    {
      if (params.nsideOut < params.nsideIn)
        cdo_abort("Parameter nside=%d must be greater than input nside=%d!", params.nsideOut, params.nsideIn);
      params.fact = params.nsideOut / params.nsideIn;
    }
  }

  if (!cdo::is_power_of_two(params.nsideOut)) cdo_abort("Parameter nside must be a power of two!");

  if (params.orderOut == HpOrder::Undef) params.orderOut = params.orderIn;

  auto nsideOut = static_cast<size_t>(params.nsideOut);
  auto gridsize = 12 * nsideOut * nsideOut;
  auto define_healpix_func = (gridType == GRID_HEALPIX) ? define_healpix_grid : define_healpix_proj;
  return define_healpix_func(gridsize, params.nsideOut, params.orderOut);
}

class Healpix : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Healpix",
    .operators = { { "hpupgrade", HealpixHelp }, { "hpdegrade", HealpixHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Healpix> registration = RegisterEntry<Healpix>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  bool doDegrade{};

  Parameter params{};

  VarList varList1{};
  VarList varList2{};

public:
  void
  init() override
  {
    auto HPDEGRADE = module.get_id("hpdegrade");

    auto operatorID = cdo_operator_id();
    doDegrade = (operatorID == HPDEGRADE);

    params = get_parameter();
    params.doDegrade = doDegrade;
    verify_parameter(params);
    if (params.zoom > 0) params.nsideOut = std::lround(std::pow(2, params.zoom));

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids > 1) cdo_abort("Too many different grids!");

    auto gridID = vlistGrid(vlistID1, 0);
    if (!is_healpix_grid(gridID)) cdo_abort("Input grid is not healpix!");

    vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto gridID2 = hp_define_grid(gridID, params);
    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  static void
  process_and_write_data(Field &field1, Field &field2, CdoStreamID streamID2, int varID, int levelID, Parameter const &params)
  {
    auto hp_func = params.doDegrade ? hp_degrade : hp_upgrade;

    if (params.orderIn == HpOrder::Ring) ring_to_nested(field1, params.nsideIn);

    hp_func(field1, field2, params);

    if (params.orderOut == HpOrder::Ring) nested_to_ring(field2, params.nsideOut);

    cdo_def_field(streamID2, varID, levelID);
    cdo_write_field(streamID2, field2);
  }

  void
  run() override
  {
    auto runAsync = (doDegrade && Options::CDO_Async_Read > 0);
    auto workerThread = runAsync ? std::make_unique<WorkerThread>() : nullptr;
    auto numTasks = runAsync ? 2 : 1;

    Field fieldVector1[2];
    Field field2;

    auto numSteps1 = varList1.numSteps();
    cdo::Progress progress(get_id());

    int numSets = 0;
    int tsID1 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID1);
      if (numFields == 0) break;

      if (numSteps1 > 1) progress.update((tsID1 + 1.0) / numSteps1);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field1 = fieldVector1[numSets % numTasks];
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        if (runAsync && numSets > 0) { workerThread->wait(); }
        numSets++;

        field2.init(varList2.vars[varID]);

        std::function<void()> process_task
            = std::bind(process_and_write_data, std::ref(field1), std::ref(field2), streamID2, varID, levelID, std::cref(params));

        runAsync ? workerThread->doAsync(process_task) : process_task();
      }

      if (runAsync) { workerThread->wait(); }

      tsID1++;
    }

    if (runAsync) { workerThread->wait(); }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
