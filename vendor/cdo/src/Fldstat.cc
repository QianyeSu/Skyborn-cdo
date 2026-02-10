/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     fldmin          Field minimum
     fldmax          Field maximum
     fldrange        Field range
     fldsum          Field sum
     fldint          Field integral
     fldmean         Field mean
     fldavg          Field average
     fldstd          Field standard deviation
     fldstd1         Field standard deviation (n-1)
     fldvar          Field variance
     fldvar1         Field variance (n-1)
     fldskew         Field skewness
     fldkurt         Field kurtosis
     fldmedian       Field median
     fldcount        Field count
     fldpctl         Field percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "pmlist.h"
#include "cdo_zaxis.h"
#include "printinfo.h"
#include "progress.h"
#include "field_functions.h"

void gridcell_areas(int gridID, Varray<double> &array);

static void
print_location_LL(double value, size_t index, int gridID, double level, bool isReg2d, size_t nlon, bool isHealpix,
                  HpParams const &hpParams, char const *funcName, std::string const &varName, CdiDateTime vDateTime)
{
  double xval{}, yval{};
  if (isHealpix)
  {
    hp_index_to_lonlat(hpParams, index, &xval, &yval);
    xval = rad_to_deg(xval);
    yval = rad_to_deg(yval);
  }
  else
  {
    auto j = index / nlon;
    auto i = index - j * nlon;
    xval = gridInqXval(gridID, isReg2d ? i : index);
    yval = gridInqYval(gridID, isReg2d ? j : index);
  }

  static auto printHeader{ true };
  if (printHeader)
  {
    fprintf(stdout, "  Date       Time          Name   Level      Cell       Lon       Lat       %s\n", funcName);
    printHeader = false;
  }

  fprintf(stdout, "%s %s %10s %7g %9zu %9.7g %9.7g %12.5g\n", date_to_string(vDateTime.date).c_str(),
          time_to_string(vDateTime.time).c_str(), varName.c_str(), level, index + 1, xval, yval, value);
}

template <typename T>
static void
print_location(int operfunc, CdoVar const &var, int levelID, double sglval, Varray<T> const &v, CdiDateTime vDateTime)
{
  auto funcName = (operfunc == FieldFunc_Min) ? "Minval" : "Maxval";

  auto isHealpix = is_healpix_grid(var.gridID);
  auto isReg2d = (var.gridType == GRID_GAUSSIAN || var.gridType == GRID_LONLAT);
  if (isHealpix || isReg2d || var.gridType == GRID_CURVILINEAR || var.gridType == GRID_UNSTRUCTURED)
  {
    T value = sglval;
    size_t ij = 0;
    for (; ij < var.gridsize; ++ij)
    {
      if (fp_is_equal(v[ij], value)) { break; }
    }
    if (ij == var.gridsize) cdo_abort("Internal error: %s not found", funcName);

    HpParams hpParams{};
    if (isHealpix) hpParams = cdo::get_healpix_params(var.gridID);
    auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
    auto nlon = gridInqXsize(var.gridID);
    print_location_LL(sglval, ij, var.gridID, level, isReg2d, nlon, isHealpix, hpParams, funcName, var.name, vDateTime);
  }
}

static void
print_location(int operfunc, const CdoVar &var, int levelID, double sglval, Field const &field, CdiDateTime vDateTime)
{
  auto func = [&](auto const &v) { print_location(operfunc, var, levelID, sglval, v, vDateTime); };
  field_operation(func, field);
}

template <typename T>
static void
field_mul_weights(Varray<T> &v1, Varray<double> const &v2, size_t numMissVals, double mv)
{
  T missval = mv;
  assert(v1.size() > 0);
  assert(v2.size() == v1.size());

  auto gridSize = v1.size();
  if (numMissVals)
  {
    for (size_t i = 0; i < gridSize; ++i)
      if (fp_is_not_equal(v1[i], missval)) { v1[i] *= v2[i]; }
  }
  else
  {
    for (size_t i = 0; i < gridSize; ++i) { v1[i] *= v2[i]; }
  }
}

static void
field_mul_weights(Field &field)
{
  auto func = [&](auto &v, auto const &w, auto numMissVals, double missval) { field_mul_weights(v, w, numMissVals, missval); };
  field_operation(func, field, field.weightv, field.numMissVals, field.missval);
}

static void
remove_global_grid_attr(int vlistID)
{
  cdiDelAtt(vlistID, CDI_GLOBAL, "ICON_grid_file_uri");
  cdiDelAtt(vlistID, CDI_GLOBAL, "number_of_grid_used");
  cdiDelAtt(vlistID, CDI_GLOBAL, "uuidOfHGrid");
}

static int
gen_target_gridpoint(int gridID1)
{
  int gridID2 = -1;

  auto gridType = gridInqType(gridID1);
  if (gridType == GRID_UNSTRUCTURED)
  {
    gridID2 = gridCreate(gridType, 1);
    grid_copy_names(gridID1, gridID2);
  }
  else if (gridType == GRID_GENERIC)
  {
    gridID2 = gridCreate(GRID_GENERIC, 1);
    grid_copy_names(gridID1, gridID2);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
  }
  else
  {
    gridID2 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
  }

  auto value = 0.0;
  gridDefXvals(gridID2, &value);
  gridDefYvals(gridID2, &value);

  return gridID2;
}

static void
print_weights_warning(int numGrids, std::string const &varname)
{
  if (numGrids == 1)
    cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
  else
    cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", varname);
}

namespace
{
struct Parameter
{
  double pctlNumber{};  // percentile number
  bool useWeights{ true };
  bool verbose{ false };
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params{};
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

      if (key == "weights") { params.useWeights = parameter_to_bool(value); }
      else if (key == "verbose") { params.verbose = parameter_to_bool(value); }
      else if (key == "pn") { params.pctlNumber = parameter_to_double(value); }
      else { cdo_abort("Invalid parameter key >%s<!", key); }
    }
  }

  return params;
}

static void
print_parameter(Parameter const &params)
{
  cdo_print("parameter: pn=%g  weights=%d  verbose=%d", params.pctlNumber, params.useWeights, params.verbose);
}

static int
get_gridcell_weights(Field &field, bool useWeights, bool doPrintWarning, int numGrids, std::string const &varName)
{
  auto gridSize = field.size;
  field.weightv.resize(gridSize);
  if (!useWeights)
  {
    cdo_print("Using constant grid cell area weights!");
    std::ranges::fill(field.weightv, 1.0);
  }

  field.weightv[0] = 1;
  if (useWeights && field.size > 1)
  {
    auto wstatus = (gridcell_weights(field.grid, field.weightv) != 0);
    if (wstatus && doPrintWarning) print_weights_warning(numGrids, varName);
  }

  return field.grid;
}

static int
get_gridcell_areas(Field &field)
{
  auto gridSize = field.size;
  field.weightv.resize(gridSize);

  gridcell_areas(field.grid, field.weightv);

  return field.grid;
}

class Fldstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Fldstat",
    // clang-format off
    .operators = { { "fldrange", FieldFunc_Range, 0, FldstatHelp },
                   { "fldmin", FieldFunc_Min, 0, FldstatHelp },
                   { "fldmax", FieldFunc_Max, 0, FldstatHelp },
                   { "fldsum", FieldFunc_Sum, 0, FldstatHelp },
                   { "fldint", FieldFunc_Sum, 0, FldstatHelp },
                   { "fldmean", FieldFunc_Meanw, 1, FldstatHelp },
                   { "fldavg", FieldFunc_Avgw, 1, FldstatHelp },
                   { "fldstd", FieldFunc_Stdw, 1, FldstatHelp },
                   { "fldstd1", FieldFunc_Std1w, 1, FldstatHelp },
                   { "fldvar", FieldFunc_Varw, 1, FldstatHelp },
                   { "fldvar1", FieldFunc_Var1w, 1, FldstatHelp },
                   { "fldskew", FieldFunc_Skew, 0, FldstatHelp },
                   { "fldkurt", FieldFunc_Kurt, 0, FldstatHelp },
                   { "fldmedian", FieldFunc_Median, 0, FldstatHelp },
                   { "fldcount", FieldFunc_Count, 0, FldstatHelp },
                   { "fldpctl", FieldFunc_Pctl, 0, FldstatHelp } },
    // clang-format on
    .aliases = { { "globavg", "fldavg" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Fldstat> registration = RegisterEntry<Fldstat>(module);

  int FLDINT{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  Parameter params{};

  bool isMinMaxFunc{};
  bool needWeights{};
  bool needCellarea{};

  int operfunc{};
  int numGrids{};

  VarList varList1{};

public:
  void
  init() override
  {
    FLDINT = module.get_id("fldint");

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    isMinMaxFunc = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);
    needWeights = (cdo_operator_f2(operatorID) != 0);
    needCellarea = (operatorID == FLDINT);

    auto loadGrid = (needWeights || needCellarea || (isMinMaxFunc && Options::cdoVerbose));
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (not loadGrid && this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (not loadGrid && this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    auto isIntArg = (operfunc == FieldFunc_Pctl && cdo_operator_argc() == 1 && std::isdigit(cdo_get_oper_argv()[0][0]));
    if (isIntArg) { params.pctlNumber = parameter_to_double(cdo_operator_argv(0)); }
    else { params = get_parameter(); }

    if (operfunc == FieldFunc_Pctl) { operator_check_argc(1); }

    if (Options::cdoVerbose) { print_parameter(params); }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    if (!isMinMaxFunc) vlist_unpack(vlistID2);
    remove_global_grid_attr(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    numGrids = vlistNumGrids(vlistID1);

    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridID2 = gen_target_gridpoint(gridID1);
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    bool printLocation = ((Options::cdoVerbose || params.verbose) && isMinMaxFunc);
    Field field;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int lastgrid = -1;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      auto vDateTime = taxisInqVdatetime(taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto fstatus = ((tsID + (fieldID + 1.0) / numFields) / numSteps);
        if (numSteps > 1) progress.update(fstatus);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        auto doPrintWarning = (tsID == 0 && levelID == 0);
        if (needWeights && field.grid != lastgrid)
          lastgrid = get_gridcell_weights(field, params.useWeights, doPrintWarning, numGrids, var.name);
        else if (needCellarea && field.grid != lastgrid)
          lastgrid = get_gridcell_areas(field);

        if (needCellarea) field_mul_weights(field);

        auto singleValue = (operfunc == FieldFunc_Pctl) ? field_pctl(field, params.pctlNumber) : field_function(field, operfunc);

        if (printLocation) { print_location(operfunc, var, levelID, singleValue, field, vDateTime); }

        size_t numMissVals = fp_is_equal(singleValue, field.missval);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, &singleValue, numMissVals);
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
