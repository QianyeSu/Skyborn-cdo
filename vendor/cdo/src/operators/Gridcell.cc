/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Gridcell   gridarea        Grid cell area in m^2
      Gridcell   gridweights     Grid cell weights
      Gridcell   gridmask        Grid mask
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include <mpim_grid.h>
#include "constants.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "util_string.h"

double gridGetPlanetRadius(int gridID);

double
planet_radius_env(void)
{
  static double planetRadiusEnv{ 0 };
  static auto readFromEnv = true;
  if (readFromEnv)
  {
    readFromEnv = false;
    auto envString = getenv_string("PLANET_RADIUS");
    if (envString.size())
    {
      auto fval = radius_str_to_meter(envString);
      if (is_not_equal(fval, 0)) planetRadiusEnv = fval;
    }
  }

  return planetRadiusEnv;
}

double
get_planet_radius_in_meter(int gridID)
{
  auto planetRadiusInMeter{ PlanetRadiusDefault };

  if (is_not_equal(planet_radius_env(), 0))
  {
    planetRadiusInMeter = planet_radius_env();
    cdo_print("Using planet radius from env.var. PLANET_RADIUS: %.8gm", planetRadiusInMeter);
  }
  else if (gridGetPlanetRadius(gridID) > 1)
  {
    planetRadiusInMeter = gridGetPlanetRadius(gridID);
    cdo_print("Using planet radius from grid description: %.8gm", planetRadiusInMeter);
  }
  else { cdo_print("Using default planet radius: %.8gm", planetRadiusInMeter); }

  return planetRadiusInMeter;
}

static void
gridcell_areas(int gridID, Varray<double> &array, double planetRadiusInMeter)
{
  auto gridType = gridInqType(gridID);

  if (gridProjIsSupported(gridID) || gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_GME
      || gridType == GRID_CURVILINEAR || gridType == GRID_UNSTRUCTURED || gridType == GRID_GAUSSIAN_REDUCED
      || gridType == GRID_HEALPIX)
  {
    if (gridHasArea(gridID))
    {
      cdo_print("Using existing grid cell area!");
      gridInqArea(gridID, array.data());
    }
    else
    {
      auto status = gridGenArea(gridID, array);
      // clang-format off
      if      (status == 1) cdo_abort("%s: Cell corner coordinates missing!", __func__);
      else if (status == 2) cdo_abort("%s: Can't compute grid cell area for this grid!", __func__);
      // clang-format on

      if (planetRadiusInMeter <= 0) planetRadiusInMeter = get_planet_radius_in_meter(gridID);
      auto ngp = gridInqSize(gridID);
      for (size_t i = 0; i < ngp; ++i) array[i] *= planetRadiusInMeter * planetRadiusInMeter;
    }
  }
  else { cdo_abort("%s: Unsupported grid type: %s", __func__, gridNamePtr(gridType)); }
}

void
gridcell_areas(int gridID, Varray<double> &array)
{
  gridcell_areas(gridID, array, 0.0);
}

static void
grid_dx(int gridID, Varray<double> &array, double planetRadiusInMeter)
{
  if (planetRadiusInMeter <= 0) planetRadiusInMeter = get_planet_radius_in_meter(gridID);

  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xv(gridsize), yv(gridsize);
  gridInqXvals(gridID, xv.data());
  gridInqYvals(gridID, yv.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xv, "grid longitudes");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yv, "grid latitudes");

  for (size_t j = 0; j < ysize; ++j)
  {
    auto joff = j * xsize;
    for (size_t i = 0; i < xsize; ++i)
    {
      double len1, len2;
      if (i == 0)
      {
        len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joff + i + 1], yv[joff + i + 1]);
        len1 = len2;
      }
      else if (i == (xsize - 1))
      {
        len1 = orthodrome(xv[joff + i - 1], yv[joff + i - 1], xv[joff + i], yv[joff + i]);
        len2 = len1;
      }
      else
      {
        len1 = orthodrome(xv[joff + i - 1], yv[joff + i - 1], xv[joff + i], yv[joff + i]);
        len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joff + i + 1], yv[joff + i + 1]);
      }

      array[joff + i] = 0.5 * (len1 + len2) * planetRadiusInMeter;
    }
  }
}

static void
grid_dy(int gridID, Varray<double> &array, double planetRadiusInMeter)
{
  if (planetRadiusInMeter <= 0) planetRadiusInMeter = get_planet_radius_in_meter(gridID);

  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xv(gridsize), yv(gridsize);
  gridInqXvals(gridID, xv.data());
  gridInqYvals(gridID, yv.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xv, "grid longitudes");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yv, "grid latitudes");

  for (size_t i = 0; i < xsize; ++i)
  {
    for (size_t j = 0; j < ysize; ++j)
    {
      auto joff = j * xsize;
      auto joffp1 = (j + 1) * xsize;
      auto joffm1 = (j - 1) * xsize;
      double len1, len2;
      if (j == 0)
      {
        len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joffp1 + i], yv[joffp1 + i]);
        len1 = len2;
      }
      else if (j == (ysize - 1))
      {
        len1 = orthodrome(xv[joffm1 + i], yv[joffm1 + i], xv[joff + i], yv[joff + i]);
        len2 = len1;
      }
      else
      {
        len1 = orthodrome(xv[joffm1 + i], yv[joffm1 + i], xv[joff + i], yv[joff + i]);
        len2 = orthodrome(xv[joff + i], yv[joff + i], xv[joffp1 + i], yv[joffp1 + i]);
      }

      array[joff + i] = 0.5 * (len1 + len2) * planetRadiusInMeter;
    }
  }
}

static void
get_parameter(double &radiusInMeter)
{
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
      if (key == "radius") radiusInMeter = radius_str_to_meter(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Gridcell : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Gridcell",
    .operators = { { "gridarea", 1, 0, GridcellHelp },
                   { "gridweights", GridcellHelp },
                   { "gridmask", GridcellHelp },
                   { "griddx", 1, 0, GridcellHelp },
                   { "griddy", 1, 0, GridcellHelp },
                   { "gridcellidx", GridcellHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Gridcell> registration = RegisterEntry<Gridcell>();

  int GRIDAREA{}, GRIDWEIGHTS{}, GRIDMASK{}, GRIDDX{}, GRIDDY{}, GRIDCELLIDX{};

  CdoStreamID streamID1{};

  int operatorID{};
  int gridID{};
  int vlistID2{ CDI_UNDEFID };

  size_t gridsize{};

  double planetRadiusInMeter{ 0 };

  Varray<double> array{};

public:
  void
  define_var(int varID)
  {
    if (operatorID == GRIDAREA)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cell_area");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, "area");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "area of grid cell");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m2");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
    else if (operatorID == GRIDWEIGHTS)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cell_weights");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
    else if (operatorID == GRIDMASK)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "grid_mask");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_UINT8);
    }
    else if (operatorID == GRIDDX)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "dx");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "delta x");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m");
    }
    else if (operatorID == GRIDDY)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "dy");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "delta y");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m");
    }
    else if (operatorID == GRIDCELLIDX)
    {
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "gridcellidx");
      cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "grid cell index");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_INT32);
    }
  }

  void
  init() override
  {
    GRIDAREA = module.get_id("gridarea");
    GRIDWEIGHTS = module.get_id("gridweights");
    GRIDMASK = module.get_id("gridmask");
    GRIDDX = module.get_id("griddx");
    GRIDDY = module.get_id("griddy");
    GRIDCELLIDX = module.get_id("gridcellidx");

    operatorID = cdo_operator_id();

    auto needRadius = (cdo_operator_f1(operatorID) > 0);

    if (needRadius && cdo_operator_argc() == 1)
    {
      double radiusInMeter{ 0 };
      get_parameter(radiusInMeter);
      if (is_not_equal(radiusInMeter, 0))
      {
        planetRadiusInMeter = radiusInMeter;
        cdo_print("Using user defined planet radius: %.8gm", planetRadiusInMeter);
      }
    }
    else { operator_check_argc(0); }

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    auto numGrids = vlistNumGrids(vlistID1);

    if (numGrids > 1) cdo_warning("Found more than 1 grid, using the first one!");

    gridID = vlistGrid(vlistID1, 0);
    auto zaxisID = zaxis_from_name("surface");

    vlistID2 = vlistCreate();
    auto varID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_CONSTANT);
    vlistDefNtsteps(vlistID2, 0);

    define_var(varID);

    auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID);

    gridsize = gridInqSize(gridID);
    array.resize(gridsize);
  }

  void
  run() override
  {
    if (operatorID == GRIDAREA) { gridcell_areas(gridID, array, planetRadiusInMeter); }
    else if (operatorID == GRIDWEIGHTS)
    {
      auto status = gridcell_weights(gridID, array);
      if (status != 0) cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
    }
    else if (operatorID == GRIDMASK)
    {
      std::vector<int> mask(gridsize, 1);
      if (gridInqMask(gridID, nullptr)) gridInqMask(gridID, &mask[0]);

      for (size_t i = 0; i < gridsize; ++i) array[i] = mask[i];
    }
    else if (operatorID == GRIDCELLIDX)
    {
      for (size_t i = 0; i < gridsize; ++i) array[i] = i + 1;
    }
    else if (operatorID == GRIDDX || operatorID == GRIDDY)
    {
      auto gridtype = gridInqType(gridID);
      if (gridProjIsSupported(gridID) || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR)
      {
        if (gridtype != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, NeedCorners::Yes);

        (operatorID == GRIDDX) ? grid_dx(gridID, array, planetRadiusInMeter) : grid_dy(gridID, array, planetRadiusInMeter);
      }
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

    auto streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_timestep(streamID2, 0);
    cdo_def_field(streamID2, 0, 0);
    cdo_write_field(streamID2, array.data(), 0);
    cdo_stream_close(streamID2);
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
  }
};
