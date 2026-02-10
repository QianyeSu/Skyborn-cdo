/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Remapweights  genbil          Generate bilinear interpolation weights
      Remapweights  genbic          Generate bicubic interpolation weights
      Remapweights  genknn          Generate k-nearest neighbor weights
      Remapweights  gennn           Generate nearest neighbor weights
      Remapweights  gendis          Generate distance-weighted averaging weights
      Remapweights  gencon          Generate YAC first order conservative remap weights
      Remapweights  genlaf          Generate largest area fraction weights
*/

#include <algorithm>
#include <thread>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "remap_utils.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "cdo_options.h"
#include "util_string.h"
#include "pmlist.h"

int get_remapIndex(int numRemaps, std::vector<RemapType> &remapList, int gridID, size_t numMissVals, bool useMask,
                   Vmask const &imask);
void pack_gme_vgpm(Varray<int> const &vgpm, Vmask &imask);

namespace
{
struct RemapweightsParams
{
  std::string gridString;
  KnnParams knnParams;
};
}  // namespace

static RemapweightsParams
get_parameter_knn()
{
  RemapweightsParams params;

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
      if      (key == "k")           params.knnParams.k = parameter_to_int(value);
      else if (key == "kmin")        params.knnParams.kMin = parameter_to_int(value);
      else if (key == "weighted")    params.knnParams.weighted = string_to_weightingMethod(parameter_to_word(value));
      else if (key == "gauss_scale") params.knnParams.gaussScale = parameter_to_double(value);
      else if (key == "extrapolate") params.knnParams.extrapolate = parameter_to_bool(value);
      else if (key == "grid")        params.gridString = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
print_parameter(const RemapweightsParams &params)
{
  std::stringstream outbuffer;
  outbuffer << "grid=" << params.gridString;
  outbuffer << ", k=" << params.knnParams.k;
  outbuffer << ", kmin=" << params.knnParams.kMin;
  outbuffer << ", weighted=" << weightingMethod_to_string(params.knnParams.weighted);
  outbuffer << ", gauss_scale=" << params.knnParams.gaussScale;
  outbuffer << ", extrapolate=" << params.knnParams.extrapolate;

  cdo_verbose("%s", outbuffer.str());
}

static void
get_parameter_map3d(int offset, int &neighbors, bool &map3D, std::string &grid)
{
  auto numArgs = cdo_operator_argc() - offset;
  if (numArgs)
  {
    auto argsList = cdo_get_oper_argv();
    argsList.erase(argsList.begin(), argsList.begin() + offset);

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argsList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "grid")                    grid = parameter_to_word(value);
      else if (key == "neighbors")               neighbors = parameter_to_int(value);
      else if (key == "map3D" || key == "map3d") map3D = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

static void
remap_write_weights(std::string const &remapWeightsFile, KnnParams const &knnParams, const RemapSwitches &remapSwitches,
                    RemapType &remap)
{
  remap_write_data_scrip(remapWeightsFile, knnParams, remapSwitches, remap.srcGrid, remap.tgtGrid, remap.vars);

  constexpr auto removeMask{ false };
  remap_vars_free(remap.vars);
  remap_grid_free(remap.srcGrid, removeMask);
  remap_grid_free(remap.tgtGrid);
}

class Remapweights : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Remapweights",
    .operators = { { "genbil", GENBIL, 0, RemapbilHelp },
                   { "genbic", GENBIC, 0, RemapbicHelp },
                   { "genknn", GENKNN, 0, {} },
                   { "gennn", GENNN, 0, RemapnnHelp },
                   { "gendis", GENDIS, 0, RemapdisHelp },
                   { "gencon", GENCON, 0, RemapconHelp },
                   { "genycon2test", GENYCON2, 0, nullptr },
                   { "genlaf", GENLAF, 0, RemaplafHelp } },
    .aliases = { { "genycon", "gencon" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Remapweights> registration = RegisterEntry<Remapweights>(module);

  KnnParams knnParams{};
  RemapSwitches remapSwitches{};
  int numRemaps{ 0 };
  int numNeighbors{ 0 };
  bool map3D{ false };

  CdoStreamID streamID1{};
  int vlistID1{ CDI_UNDEFID };

  int gridID2{};

  bool useMask{};
  bool extrapolateIsSet{};
  bool remapExtrapolate{};
  bool needGradients{};

  int operfunc{};
  int maxRemaps{};

  VarList varList1{};
  std::vector<bool> remapGrids{};
  std::vector<RemapType> remapList{};

  RemapDefaults remapDefaults{};
  RemapMethod mapType{};
  int remapOrder{};

  NormOpt normOpt{};

  bool remap_genweights{ true };

  Field field1{};

public:
  std::string
  get_parameter()
  {
    std::string targetGridName;

    if (operfunc == GENKNN)
    {
      auto remapParams = get_parameter_knn();
      if (Options::cdoVerbose) print_parameter(remapParams);
      if (remapParams.gridString.empty()) cdo_abort("grid parameter missing!");
      targetGridName = remapParams.gridString;
      knnParams = remapParams.knnParams;
      if (knnParams.kMin == 0) knnParams.kMin = knnParams.k;
    }
    else
    {
      operator_input_arg("grid description file or name");
      targetGridName = cdo_operator_argv(0);
      int offset = targetGridName.starts_with("grid=") ? 0 : 1;
      if (cdo_operator_argc() > offset)
      {
        int numNeighborsParam = 0;
        get_parameter_map3d(offset, numNeighborsParam, map3D, targetGridName);
        if (map3D) remapDefaults.genMultiWeights = 1;
        if (operfunc == GENDIS)
        {
          if (numNeighborsParam < 0) cdo_abort("Number of nearest neighbors out of range (>0)!");
          if (numNeighborsParam > 0) numNeighbors = numNeighborsParam;
        }
      }
      else { operator_check_argc(1); }
    }

    return targetGridName;
  }

  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    auto writeRemapWeightsOnly = true;

    remap_set_int(REMAP_WRITE_REMAP, writeRemapWeightsOnly);

    remapDefaults = remap_get_params();
    extrapolateIsSet = (remapDefaults.extrapolate != -1);
    remapExtrapolate = extrapolateIsSet ? (bool) remapDefaults.extrapolate : remap_func_is_dist(operfunc);

    if (Options::cdoVerbose) cdo_print("Extrapolation %s!", remapExtrapolate ? "enabled" : "disabled");

    auto targetGridName = get_parameter();

    gridID2 = cdo_define_grid(targetGridName);
    if (gridInqType(gridID2) == GRID_GENERIC) cdo_abort("Unsupported target grid type (generic)!");

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    varList1 = VarList(vlistID1);

    auto findOnlyFirst = true;
    remapGrids = remap_set_grids(vlistID1, varList1, findOnlyFirst);
    auto numRemapGrids = std::ranges::count_if(remapGrids, [](auto flag) { return (flag == true); });
    if (numRemapGrids == 0) cdo_abort("No remappable grid found!");

    maxRemaps = remapDefaults.maxRemaps;
    if (maxRemaps == -1) maxRemaps = remap_get_max_maps(vlistID1);
    if (maxRemaps < 1) cdo_abort("maxRemaps out of range (>0)!");

    remapList.resize(maxRemaps);

    remapSwitches = remap_operfunc_to_maptype(operfunc);
    if (remapSwitches.mapType == RemapMethod::KNN)
    {
      if (numNeighbors) remapSwitches.numNeighbors = numNeighbors;
      if (remapSwitches.numNeighbors != -1)
      {
        knnParams.k = remapSwitches.numNeighbors;
        knnParams.kMin = 1;
        knnParams.extrapolate = remapExtrapolate;
        knnParams.weighted = WeightingMethod::distanceWeighted;
      }
    }

    mapType = remapSwitches.mapType;
    remapOrder = remapSwitches.remapOrder;
    useMask = !(!remap_genweights
                && (mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC || mapType == RemapMethod::KNN
                    || mapType == RemapMethod::CONSERV));

    remap_set_int(REMAP_GENWEIGHTS, (int) remap_genweights);

    normOpt = NormOpt(NormOpt::NONE);
    if (mapType == RemapMethod::CONSERV) normOpt = remap_get_normOpt();

    needGradients = (mapType == RemapMethod::BICUBIC);
    if (mapType == RemapMethod::CONSERV && remapOrder == 2)
    {
      if (Options::cdoVerbose) cdo_print("Second order remapping");
      needGradients = true;
    }
  }

  void
  run() override
  {
    std::thread writeWorker;
    Vmask imask;

    int tsID = 0;
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID);

    int gridIDout = -1;
    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto &var = varList1.vars[varID];
      field1.init(var);
      cdo_read_field(streamID1, field1);
      auto numMissVals1 = useMask ? field1.numMissVals : 0;

      auto gridIndex = vlistGridIndex(vlistID1, var.gridID);
      if (remapGrids[gridIndex])
      {
        if (numRemaps == 0) { gridIDout = var.gridID; }
        else if (gridIDout != var.gridID) { continue; }

        if (mapType != RemapMethod::CONSERV && var.gridType == GRID_GME)
          cdo_abort("Only conservative remapping is available to remap between GME grids!");

        if (gridIsCircular(var.gridID) && !extrapolateIsSet) remapExtrapolate = true;

        remap_set_mask(field1, var.gridsize, numMissVals1, var.missval, imask);

        int remapIndex = get_remapIndex(numRemaps, remapList, var.gridID, numMissVals1, useMask, imask);
        if (remapIndex >= 0) continue;

        if (numRemaps >= maxRemaps) break;

        remapIndex = numRemaps;
        numRemaps++;

        auto &remap = remapList[remapIndex];

        //  remap.srcGrid.luse_cell_area = false;
        //  remap.tgtGrid.luse_cell_area = false;

        remap.vars.normOpt = normOpt;

        if ((mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC)
            && (var.gridType == GRID_GME || var.gridType == GRID_UNSTRUCTURED))
          cdo_abort("Bilinear/bicubic interpolation doesn't support unstructured source grids!");

        // Initialize grid information for both grids
        remap_init_grids(mapType, remapExtrapolate, var.gridID, remap.srcGrid, gridID2, remap.tgtGrid);
        remap_search_init(mapType, remap.search, remap.srcGrid, remap.tgtGrid);

        remap.gridID = var.gridID;
        remap.numMissVals = numMissVals1;

        if (var.gridType == GRID_GME) { pack_gme_vgpm(remap.srcGrid.vgpm, imask); }

        varray_copy(remap.srcGrid.size, imask, remap.srcGrid.mask);

        if (mapType == RemapMethod::CONSERV)
        {
          std::ranges::fill(remap.srcGrid.cellArea, 0.0);
          std::ranges::fill(remap.srcGrid.cellFrac, 0.0);
          std::ranges::fill(remap.tgtGrid.cellArea, 0.0);
        }
        std::ranges::fill(remap.tgtGrid.cellFrac, 0.0);

        // initialize some remapping variables
        remap_vars_init(mapType, remapOrder, remap.vars);

        remap_print_info(operfunc, remap_genweights, remap.srcGrid, remap.tgtGrid, numMissVals1, knnParams);

        if (needGradients && remap.srcGrid.rank != 2 && remapOrder == 2)
        {
          cdo_abort("Second order remapping is not available for unstructured grids!");
        }

        remap_gen_weights(remapSwitches.mapType, knnParams, remap);

        std::string outFile = cdo_get_stream_name(1);
        if (remapDefaults.genMultiWeights) { outFile += string_format("%05d", numRemaps) + ".nc"; }

        // remap_write_weights(outFile, remapSwitches, remap);

        if (numRemaps > 1) writeWorker.join();

        writeWorker = std::thread(remap_write_weights, outFile, knnParams, remapSwitches, std::ref(remap));

        if (!remapDefaults.genMultiWeights) break;
      }
    }

    writeWorker.join();

    for (int remapIndex = 0; remapIndex < numRemaps; remapIndex++)
    {
      auto &remap = remapList[remapIndex];
      remap_grid_free(remap.srcGrid);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
  }
};
