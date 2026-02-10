/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Interpolate remap           Grid remapping
      Interpolate remapbil        Bilinear interpolation
      Interpolate remapbic        Bicubic interpolation
      Interpolate remapknn        K-nearest neighbor remapping
      Interpolate remapnn         Nearest-neighbor remapping
      Interpolate remapdis        Distance-weighted averaging
      Interpolate remapcon        YAC first order conservative remapping
      Interpolate remaplaf        Largest area fraction remapping
*/

#include <algorithm>

#include <cdi.h>

#include "grid_convert.h"
#include "process_int.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "remap_utils.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "cdo_options.h"
#include "util_wildcards.h"
#include "remapknn.h"

int
get_remapIndex(int numRemaps, std::vector<RemapType> &remaps, int gridID, size_t numMissVals, bool useMask, Vmask const &imask)
{
  int remapIndex = -1;
  for (remapIndex = numRemaps - 1; remapIndex >= 0; remapIndex--)
  {
    auto &remap = remaps[remapIndex];
    if (gridID == remap.gridID)
    {
      if ((useMask && numMissVals == remap.numMissVals && imask == remap.srcGrid.mask) || !useMask)
      {
        remap.nused++;
        break;
      }
    }
  }

  return remapIndex;
}

void
pack_gme_vgpm(Varray<int> const &vgpm, Vmask &imask)
{
  auto n = imask.size();
  assert(n <= vgpm.size());
  for (size_t i = 0, j = 0; i < n; ++i)
    if (vgpm[i]) imask[j++] = imask[i];
}

static int
maptype_to_operfunc(const RemapSwitches &remapSwitches)
{
  int operfunc = -1;

  if (remapSwitches.mapType == RemapMethod::CONSERV)
    operfunc = (remapSwitches.submapType == SubmapType::LAF) ? REMAPLAF : ((remapSwitches.remapOrder == 2) ? REMAPYCON2 : REMAPCON);
  else if (remapSwitches.mapType == RemapMethod::BILINEAR)
    operfunc = REMAPBIL;
  else if (remapSwitches.mapType == RemapMethod::BICUBIC)
    operfunc = REMAPBIC;
  else if (remapSwitches.mapType == RemapMethod::KNN)
  {
    if (remapSwitches.numNeighbors == -1)
      operfunc = REMAPKNN;
    else
      operfunc = (remapSwitches.numNeighbors == 1) ? REMAPNN : REMAPDIS;
  }
  else
    cdo_abort("Unsupported mapping method (mapType = %d)", remapSwitches.mapType);

  return operfunc;
}

static void
scale_gridbox_area(Field const &field1, Field &field2, Varray<double> const &grid2_area)
{
  auto gridsize1 = field1.size;
  auto gridsize2 = field2.size;
  double v1sum;
  auto func = [&](auto const &v1, auto &v2)
  {
    v1sum = varray_sum(gridsize1, v1);
    auto v2sum = varray_sum(gridsize2, grid2_area);
    for (size_t i = 0; i < gridsize2; ++i) v2[i] = grid2_area[i] / v2sum * v1sum;
  };
  field_operation2(func, field1, field2);
  static auto hasGridboxInfo = true;
  if (hasGridboxInfo)
  {
    cdo_print("gridbox_area replaced and scaled to %g", v1sum);
    hasGridboxInfo = false;
  }
}

template <typename T>
void gme_grid_restore(T *p, int ni, int nd);

template <typename T>
static void
restore_gme_grid(int gridID, size_t srvGridSize, Varray<T> &v, size_t tgtGridSize, Varray<int> const &vgpm)
{
  int nd, ni, ni2, ni3;
  gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);

  auto n = tgtGridSize;
  for (size_t i = srvGridSize; i > 0; i--)
    if (vgpm[i - 1]) v[i - 1] = v[--n];

  gme_grid_restore(v.data(), ni, nd);
}

static void
restore_gme_grid(Field &field, const RemapGrid &tgtGrid)
{
  auto func = [&](auto &v) { restore_gme_grid(field.grid, field.gridsize, v, tgtGrid.size, tgtGrid.vgpm); };
  field_operation(func, field);
}

template <typename T>
static void
store_gme_grid(size_t gridsize, Varray<T> &v, Varray<int> const &vgpm)
{
  for (size_t i = 0, j = 0; i < gridsize; ++i)
    if (vgpm[i]) v[j++] = v[i];
}

static void
store_gme_grid(Field &field, Varray<int> const &vgpm)
{
  auto func = [&](auto &v) { store_gme_grid(field.gridsize, v, vgpm); };
  field_operation(func, field);
}

template <typename T>
static void
remap_normalize_field(NormOpt normOpt, size_t gridsize, Varray<T> &array, double mv, const RemapGrid &tgtGrid)
{
  T missval = mv;
  // used to check the result of remapcon

  if (normOpt == NormOpt::NONE)
  {
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_not_equal(array[i], missval))
      {
        auto gridError = tgtGrid.cellFrac[i] * tgtGrid.cellArea[i];
        array[i] = (std::fabs(gridError) > 0.0) ? (array[i] / gridError) : missval;
      }
    }
  }
  else if (normOpt == NormOpt::DESTAREA)
  {
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_not_equal(array[i], missval))
      {
        array[i] = (std::fabs(tgtGrid.cellFrac[i]) > 0.0) ? (array[i] / tgtGrid.cellFrac[i]) : missval;
      }
    }
  }
}

static void
remap_normalize_field(NormOpt normOpt, Field &field, const RemapGrid &tgtGrid)
{
  auto func = [&](auto &v) { remap_normalize_field(normOpt, field.gridsize, v, field.missval, tgtGrid); };
  field_operation(func, field);
}

template <typename T>
static void
remap_set_fracmin(double fracMin, size_t gridsize, Varray<T> &array, double mv, const RemapGrid *tgtGrid)
{
  if (fracMin > 0.0)
  {
    T missval = mv;
    for (size_t i = 0; i < gridsize; ++i)
      if (tgtGrid->cellFrac[i] < fracMin) array[i] = missval;
  }
}

static void
remap_set_fracmin(double fracMin, Field &field, const RemapGrid *tgtGrid)
{
  auto func = [&](auto &v) { remap_set_fracmin(fracMin, field.gridsize, v, field.missval, tgtGrid); };
  field_operation(func, field);
}

static void
remap_field(RemapMethod mapType, KnnParams const &knnParams, RemapType &remap, Field const &field1, Field &field2)
{
  // clang-format off
  if      (mapType == RemapMethod::BILINEAR) remap_bilinear(remap.search, field1, field2);
  else if (mapType == RemapMethod::BICUBIC)  remap_bicubic(remap.search, field1, field2);
  else if (mapType == RemapMethod::KNN)      remap_knn(knnParams, remap.search, field1, field2);
  else if (mapType == RemapMethod::CONSERV)  remap_conserv(remap.vars.normOpt, remap.search, field1, field2);
  // clang-format on
}

static RemapSwitches
remap_read_weights(std::string const &remapWeightsFile, int gridID1, int gridID2, RemapType &remap0, bool extrapolateIsSet,
                   bool &remapExtrapolate, KnnParams &knnParams)
{
  auto remapSwitches
      = remap_read_data_scrip(remapWeightsFile, gridID1, gridID2, remap0.srcGrid, remap0.tgtGrid, remap0.vars, knnParams);

  auto gridsize = remap0.srcGrid.size;
  remap0.gridID = gridID1;
  remap0.gridsize = gridInqSize(gridID1);

  if (remapSwitches.mapType == RemapMethod::KNN && !extrapolateIsSet) remapExtrapolate = true;
  if (gridIsCircular(gridID1) && !extrapolateIsSet) remapExtrapolate = true;

  if (remapSwitches.mapType == RemapMethod::KNN)
  {
    if (remapSwitches.numNeighbors == -1)
    {
      if (knnParams.kMin == 0) knnParams.kMin = knnParams.k;
    }
    else
    {
      knnParams.k = remapSwitches.numNeighbors;
      knnParams.kMin = 1;
      knnParams.extrapolate = remapExtrapolate;
      knnParams.weighted = WeightingMethod::distanceWeighted;
    }
  }

  if (gridInqType(gridID1) == GRID_GME) gridsize = remap0.srcGrid.nvgp;

  if (gridsize != remap0.gridsize) cdo_abort("Size of source grid and weights from %s differ!", remapWeightsFile);

  if (gridInqType(gridID1) == GRID_GME) gridsize = remap0.srcGrid.size;

  for (size_t i = 0; i < gridsize; ++i)
    if (remap0.srcGrid.mask[i] == false) remap0.numMissVals++;

  auto gridsize2 = gridInqSize(gridID2);
  if (gridInqType(gridID2) == GRID_GME)
  {
    remap0.tgtGrid.nvgp = gridInqSize(gridID2);
    remap0.tgtGrid.vgpm.resize(gridInqSize(gridID2));
    auto gridID2_gme = gridToUnstructured(gridID2, NeedCorners::Yes);
    gridInqMaskGME(gridID2_gme, &remap0.tgtGrid.vgpm[0]);
    gridDestroy(gridID2_gme);
    size_t isize = 0;
    for (size_t i = 0; i < gridsize2; ++i)
      if (remap0.tgtGrid.vgpm[i]) isize++;
    gridsize2 = isize;
  }
  // printf("grid2 %zu %d %zu\n", gridsize2, remap0.tgtGrid.nvgp, remap0.tgtGrid.size);
  if (remap0.tgtGrid.size != gridsize2) cdo_abort("Size of target grid and weights from %s differ!", remapWeightsFile);

  return remapSwitches;
}

double get_planet_radius_in_meter(int gridID);

static void
print_node_info(int gridID2, RemapVars const &rv, const RemapGrid &srcGrid)
{
  static auto printNodeInfo = true;
  if (printNodeInfo == false) return;
  printNodeInfo = false;

  double lon2{}, lat2{};
  if (gridInqSize(gridID2) == 1)
  {
    gridInqXvals(gridID2, &lon2);
    gridInqYvals(gridID2, &lat2);
  }

  if (rv.numLinks == 1)
  {
    auto srcCellIndex = rv.srcCellIndices[0];
    auto pointLL = remapgrid_get_lonlat(&srcGrid, srcCellIndex);
    auto lon1 = pointLL.lon();
    auto lat1 = pointLL.lat();
    if (lon1 > PI) lon1 -= PI2;
    auto distance = orthodrome(lon1, lat1, DEG2RAD * lon2, DEG2RAD * lat2) * get_planet_radius_in_meter(srcGrid.gridID) / 1000;

    cdo_print("Target Point: lon=%g/lat=%g  Source Point: index=%zu lon=%g/lat=%g distance=%.3fkm", lon2, lat2, srcCellIndex + 1,
              RAD2DEG * lon1, RAD2DEG * lat1, distance);
  }
  else if (rv.numLinks == 0) { cdo_print("Target Point: lon=%g/lat=%g  Source Point: not found", lon2, lat2); }
}

RemapknnParams
remapknn_get_parameter()
{
  RemapknnParams params;
  auto &knnParams = params.knnParams;

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
      if      (key == "k")           knnParams.k = parameter_to_int(value);
      else if (key == "kmin")        knnParams.kMin = parameter_to_int(value);
      else if (key == "weighted")    knnParams.weighted = string_to_weightingMethod(parameter_to_word(value));
      else if (key == "gauss_scale") knnParams.gaussScale = parameter_to_double(value);
      else if (key == "rbf_scale")   knnParams.rbfScale = parameter_to_double(value);
      else if (key == "extrapolate") knnParams.extrapolate = parameter_to_bool(value);
      else if (key == "grid")        params.gridString = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

void
remapknn_verify_parameter(KnnParams const &knnParams)
{
  if (knnParams.weighted == WeightingMethod::gaussWeighted)
  {
    if (knnParams.k <= 1) cdo_abort("Parameter k must be greater than 1!");
    if (knnParams.kMin <= 1) cdo_abort("Parameter kmin must be greater than 1!");
  }
}

void
print_knn_parameter(KnnParams const &knnParams, std::string const &prefix)
{
  std::stringstream outbuffer;
  outbuffer << prefix;
  outbuffer << "k=" << knnParams.k;
  outbuffer << ", kmin=" << knnParams.kMin;
  outbuffer << ", weighted=" << weightingMethod_to_string(knnParams.weighted);
  if (knnParams.weighted == WeightingMethod::gaussWeighted) outbuffer << ", gauss_scale=" << knnParams.gaussScale;
  if (knnParams.weighted == WeightingMethod::rbf) outbuffer << ", rbf_scale=" << knnParams.rbfScale;
  outbuffer << ", extrapolate=" << knnParams.extrapolate;

  cdo_verbose("%s", outbuffer.str());
}

void
remapknn_print_parameter(RemapknnParams const &params)
{
  print_knn_parameter(params.knnParams, "grid=" + params.gridString + ", ");
}

class Remapgrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Remapgrid",
    // clang-format off
    .operators = { { "remap", REMAP, 0, RemapHelp },
                   { "remapbil", REMAPBIL, 0, RemapbilHelp },
                   { "remapbic", REMAPBIC, 0, RemapbicHelp },
                   { "remapknn", REMAPKNN, 0, {} },
                   { "remapnn", REMAPNN, 0, RemapnnHelp },
                   { "remapdis", REMAPDIS, 0, RemapdisHelp },
                   { "remapcon", REMAPCON, 1, RemapconHelp },
                   { "remapycon2test", REMAPYCON2, 1, nullptr},
                   { "remaplaf", REMAPLAF, 1, RemaplafHelp },
                   { "remapavgtest", REMAPAVG, 0, nullptr} },
    // clang-format on
    .aliases = { { "remapycon", "remapcon" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Remapgrid> registration = RegisterEntry<Remapgrid>(module);

  KnnParams knnParams{};
  RemapSwitches remapSwitches{};
  bool remap_genweights = Options::REMAP_genweights;
  int numRemaps = 0;
  int numNeighbors = 0;
  std::vector<std::string> remapWeightsFiles{};

  CdoStreamID streamID1{};
  int taxisID1{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };

  CdoStreamID streamID2{};
  int taxisID2{ CDI_UNDEFID };
  int gridID2{};

  bool useMask{};
  bool extrapolateIsSet{};
  bool remapExtrapolate{};
  bool needGradients{};
  bool doRemap{};

  int operfunc{};
  int maxRemaps{};

  std::vector<bool> remapGrids{};
  std::vector<RemapType> remapList{};

  RemapGradients gradients{};
  RemapDefaults remapDefaults{};
  RemapMethod mapType{};
  int remapOrder{};

  NormOpt normOpt{};

  VarList varList1{};
  VarList varList2{};

public:
  std::string
  get_parameter()
  {
    std::string targetGridName;

    if (doRemap)
    {
      operator_input_arg("grid description file or name, remap weights file (SCRIP NetCDF)");
      operator_check_argc(2);
      remapWeightsFiles.push_back(cdo_operator_argv(1));
      targetGridName = cdo_operator_argv(0);
    }
    else
    {
      if (operfunc == REMAPKNN)
      {
        auto remapParams = remapknn_get_parameter();
        if (Options::cdoVerbose) remapknn_print_parameter(remapParams);
        if (remapParams.gridString.empty()) cdo_abort("target grid parameter missing!");
        targetGridName = remapParams.gridString;
        knnParams = remapParams.knnParams;
        if (knnParams.kMin == 0) knnParams.kMin = knnParams.k;
        remapknn_verify_parameter(knnParams);
        if (Options::cdoVerbose) print_knn_parameter(knnParams, "KNN parameter: ");
        extrapolateIsSet = true;
        remapExtrapolate = knnParams.extrapolate;
      }

      else
      {
        operator_input_arg("grid description file or name");
        if (operfunc == REMAPDIS && cdo_operator_argc() == 2)
        {
          auto inum = parameter_to_int(cdo_operator_argv(1));
          if (inum < 1) cdo_abort("Number of nearest neighbors out of range (>0)!");
          numNeighbors = inum;
        }
        else { operator_check_argc(1); }
        targetGridName = cdo_operator_argv(0);
      }
    }

    return targetGridName;
  }

  void
  set_knn_parameter()
  {
    knnParams.k = remapSwitches.numNeighbors;
    knnParams.kMin = 1;
    knnParams.extrapolate = remapExtrapolate;
    knnParams.weighted = WeightingMethod::distanceWeighted;
  }

  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto needCellCorners = cdo_operator_f2(operatorID);
    if (not needCellCorners && this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }

    auto writeRemapWeightsOnly = false;
    doRemap = (operfunc == REMAP);

    remap_set_int(REMAP_WRITE_REMAP, writeRemapWeightsOnly);

    remapDefaults = remap_get_params();
    extrapolateIsSet = (remapDefaults.extrapolate != -1);
    remapExtrapolate = extrapolateIsSet ? (bool) remapDefaults.extrapolate : remap_func_is_dist(operfunc);

    auto targetGridName = get_parameter();

    if (Options::cdoVerbose) cdo_print("Extrapolation %s!", remapExtrapolate ? "enabled" : "disabled");

    gridID2 = cdo_define_grid(targetGridName);
    if (gridInqType(gridID2) == GRID_GENERIC) cdo_abort("Unsupported target grid type (generic)!");

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    remapGrids = remap_set_grids(vlistID1, varList1);
    auto numRemapGrids = std::ranges::count_if(remapGrids, [](auto flag) { return (flag == true); });
    if (numRemapGrids == 0) cdo_abort("No remappable grid found!");

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
      if (remapGrids[index]) vlistChangeGridIndex(vlistID2, index, gridID2);

    maxRemaps = remapDefaults.maxRemaps;
    if (maxRemaps == -1) maxRemaps = remap_get_max_maps(vlistID1);
    if (maxRemaps < 1) cdo_abort("maxRemaps out of range (>0)!");

    remapList.resize(maxRemaps);

    if (doRemap) remap_genweights = true;

    if (doRemap)
    {
      int gridIndex1;
      for (gridIndex1 = 0; gridIndex1 < numGrids; ++gridIndex1)
        if (remapGrids[gridIndex1]) break;
      auto gridID1 = vlistGrid(vlistID1, gridIndex1);

      remapWeightsFiles = expand_path_names(remapWeightsFiles);
      for (auto &remapWeightsFile : remapWeightsFiles)
      {
        // printf("read remapWeightsFile %d: %s\n", numRemaps+1, remapWeightsFile.c_str());
        auto &remap = remapList[numRemaps];
        remapSwitches
            = remap_read_weights(remapWeightsFile, gridID1, gridID2, remap, extrapolateIsSet, remapExtrapolate, knnParams);
        if (numRemaps == 0) { operfunc = maptype_to_operfunc(remapSwitches); }
        else if (operfunc != maptype_to_operfunc(remapSwitches)) { cdo_abort("Remapping method changed in input weights files!"); }
        numRemaps++;
      }
    }
    else
    {
      remapSwitches = remap_operfunc_to_maptype(operfunc);
      if (remapSwitches.mapType == RemapMethod::KNN)
      {
        if (numNeighbors) remapSwitches.numNeighbors = numNeighbors;
        if (remapSwitches.numNeighbors != -1) set_knn_parameter();
      }
    }

    mapType = remapSwitches.mapType;
    remapOrder = remapSwitches.remapOrder;

    varList2 = VarList(vlistID2);

    // if (!remap_genweights && mapType == RemapMethod::CONSERV) remap_genweights = true;
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

    if (needGradients) gradients.init(varList1.gridsizeMax());

    if (remap_genweights)
    {
      // remap() gives rounding errors on target arrays with single precision floats
      for (auto &var : varList2.vars) var.memType = MemType::Double;
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  int
  new_remapIndex()
  {
    if (numRemaps < maxRemaps)
    {
      auto remapIndex = numRemaps;
      numRemaps++;
      return remapIndex;
    }

    auto remapIndex = maxRemaps - 1;
    auto &remapLast = remapList[remapIndex];
    remap_vars_free(remapLast.vars);
    remap_grid_free(remapLast.srcGrid);
    remap_grid_free(remapLast.tgtGrid);
    remapList.pop_back();
    remapList.push_back({});
    return remapIndex;
  }

  void
  grid_search_init(RemapType &remap, const CdoVar &var)
  {
    if (gridIsCircular(var.gridID) && !extrapolateIsSet) remapExtrapolate = true;

    //  remap.srcGrid.luse_cell_area = false;
    //  remap.tgtGrid.luse_cell_area = false;

    remap.vars.normOpt = normOpt;

    if ((mapType == RemapMethod::BILINEAR || mapType == RemapMethod::BICUBIC)
        && (var.gridType == GRID_GME || var.gridType == GRID_UNSTRUCTURED))
      cdo_abort("Bilinear/bicubic interpolation doesn't support unstructured source grids!");

    // Initialize grid information for both grids
    remap_init_grids(mapType, remapExtrapolate, var.gridID, remap.srcGrid, gridID2, remap.tgtGrid);
    remap_search_init(mapType, remap.search, remap.srcGrid, remap.tgtGrid);
  }

  void
  remap_init(RemapType &remap, const CdoVar &var, size_t numMissVals1, Vmask &imask)
  {
    if (remap.gridID != var.gridID) { grid_search_init(remap, var); }

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
    // std::ranges::fill(remap.tgtGrid.cellFrac, 0.0); failed with pgc++ 23.9.0
    for (auto &cf : remap.tgtGrid.cellFrac) cf = 0.0;

    // initialize some remapping variables
    remap_vars_init(mapType, remapOrder, remap.vars);

    remap_print_info(operfunc, remap_genweights, remap.srcGrid, remap.tgtGrid, numMissVals1, knnParams);

    if (remap_genweights)
    {
      if (needGradients)
      {
        if (remap.srcGrid.rank != 2 && remapOrder == 2)
          cdo_abort("Second order remapping is not available for unstructured grids!");
      }

      remap_gen_weights(remapSwitches.mapType, knnParams, remap);
    }
  }

  void
  run() override
  {
    Field field1, field2;
    Vmask imask;

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
        auto &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID1, field1);
        auto numMissVals1 = useMask ? field1.numMissVals : 0;

        field2.init(varList2.vars[varID]);

        auto gridIndex = vlistGridIndex(vlistID1, var.gridID);

        if (!remapGrids[gridIndex]) { field_copy(field1, field2); }
        else
        {
          if (mapType != RemapMethod::CONSERV && var.gridType == GRID_GME)
            cdo_abort("Only conservative remapping is available to remap between GME grids!");

          if (gridIsCircular(var.gridID) && !extrapolateIsSet) remapExtrapolate = true;

          remap_set_mask(field1, var.gridsize, numMissVals1, var.missval, imask);

          int remapIndex = get_remapIndex(numRemaps, remapList, var.gridID, numMissVals1, useMask, imask);
          if (Options::cdoVerbose && remapIndex >= 0) cdo_print("Using remap %d", remapIndex);
          if (remapIndex < 0)
          {
            remapIndex = new_remapIndex();
            remap_init(remapList[remapIndex], var, numMissVals1, imask);
          }

          auto &remap = remapList[remapIndex];
          if (var.gridType == GRID_GME) store_gme_grid(field1, remapList[remapIndex].srcGrid.vgpm);

          auto gridsize2 = gridInqSize(gridID2);

          if (remap_genweights)
          {
            remap.nused++;

            if (needGradients)
            {
              if (remap.srcGrid.rank != 2 && remapOrder == 2)
                cdo_abort("Second order remapping is not available for unstructured grids!");

              remap::gradients(field1, remap.srcGrid, gradients);
            }

            if (Options::cdoVerbose && operfunc == REMAPNN && gridsize2 == 1)
            {
              print_node_info(gridID2, remap.vars, remap.srcGrid);
            }

            if (operfunc == REMAPLAF)
              remap_laf(field2, var.missval, gridsize2, remap.vars, field1);
            else if (operfunc == REMAPAVG)
              remap_avg(field2, var.missval, gridsize2, remap.vars, field1);
            else
              remap_field(field2, var.missval, gridsize2, remap.vars, field1, gradients);
          }
          else { remap_field(remapSwitches.mapType, knnParams, remap, field1, field2); }

          if (operfunc == REMAPCON || operfunc == REMAPYCON2)
          {
            // used only to check the result of remapcon
            if (0) remap_normalize_field(remap.vars.normOpt, field2, remap.tgtGrid);
            remap_set_fracmin(remapDefaults.fracMin, field2, &remap.tgtGrid);
            if (var.name == "gridbox_area") scale_gridbox_area(field1, field2, remap.tgtGrid.cellArea);
          }

          // calculate some statistics
          if (Options::cdoVerbose) remap::stat(remapOrder, remap.srcGrid, remap.tgtGrid, remap.vars, field1, field2);

          if (gridInqType(gridID2) == GRID_GME) restore_gme_grid(field2, remap.tgtGrid);

          field2.numMissVals = field_num_mv(field2);
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field2);
      }

      tsID++;
    }

    if (doRemap && remap_genweights && remapList[0].nused == 0)
      remap_print_warning(remapWeightsFiles[0], operfunc, remapList[0].srcGrid, remapList[0].numMissVals);

    for (int remapIndex = 0; remapIndex < numRemaps; remapIndex++)
    {
      auto &remap = remapList[remapIndex];
      remap_vars_free(remap.vars);
      remap_grid_free(remap.srcGrid);
      remap_grid_free(remap.tgtGrid);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
