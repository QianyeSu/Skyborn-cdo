/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Müller

*/

/*
   This module contains the following operators:

      Vargen     const           Create a constant field
      Vargen     random          Field with random values
      Vargen     stdatm          Field values for pressure and temperature for the standard atmosphere
*/

#include <cstdlib>
#include <cassert>
#include <cdi.h>

#include "cdo_options.h"
#include "cdo_data.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "grid_healpix.h"
#include "griddes.h"
#include "stdnametable.h"
#include "param_conversion.h"

static int
random_init(int operatorID)
{
  unsigned int seed = Options::Random_Seed;
  operator_input_arg(cdo_operator_enter(operatorID));
  if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
  if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
  auto gridID = cdo_define_grid(cdo_operator_argv(0));
  if (cdo_operator_argc() == 2)
  {
    auto idum = parameter_to_int(cdo_operator_argv(1));
    if (idum >= 0 && idum < 0x7FFFFFFF) seed = idum;
  }
  std::srand(seed);
  return gridID;
}

static void
conv_generic_grid(int gridID, size_t gridsize, Varray<double> &xvals2D, Varray<double> &yvals2D)
{
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  assert(gridsize == xsize * ysize);

  Varray<double> xcoord(xsize), ycoord(ysize);
  gridInqXvals(gridID, &xcoord[0]);
  gridInqYvals(gridID, &ycoord[0]);

  auto xrange = varray_range(xsize, xcoord);
  auto yrange = varray_range(ysize, ycoord);

  for (size_t j = 0; j < ysize; ++j)
    for (size_t i = 0; i < xsize; ++i)
    {
      xvals2D[j * xsize + i] = xcoord[i] * M_PI / xrange;
      yvals2D[j * xsize + i] = ycoord[j] * M_PI / yrange;
    }
}

static int
generate_full_point_grid_radian(int gridID, Varray<double> &xvals, Varray<double> &yvals)
{
  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Target cell center coordinates missing!");

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yvals, "grid center lat");

  return gridID;
}

static size_t
calc_index_ii(size_t nx, double xval)
{
  if (xval >= 180.0) xval -= 360.0;
  if (xval < -180.0) xval += 360.0;
  size_t ii = (xval + 180.0) * 2.0;
  if (ii >= nx) ii = nx - 1;
  return ii;
}

static size_t
calc_index_jj(size_t ny, double yval)
{
  size_t jj = (yval + 90.0) * 2.0;
  if (jj >= ny) jj = ny - 1;
  return jj;
}

static void
remap_nn_reg2d_to_reg2d(size_t nx, size_t ny, Varray<float> const &data, int gridID, Varray<float> &array)
{
  auto gridtype = gridInqType(gridID);
  if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN) cdo_abort("Internal error, wrong grid type!");

  auto nxvals = gridInqXsize(gridID);
  auto nyvals = gridInqYsize(gridID);
  Varray<double> xvals(nxvals), yvals(nyvals);

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t j = 0; j < nyvals; ++j)
  {
    auto jj = calc_index_jj(ny, yvals[j]);
    for (size_t i = 0; i < nxvals; ++i)
    {
      auto ii = calc_index_ii(nx, xvals[i]);
      array[j * nxvals + i] = data[jj * nx + ii];
    }
  }
}

static void
remap_nn_reg2d_to_nonreg2d(size_t nx, size_t ny, Varray<float> const &data, int gridID, Varray<float> &array)
{
  auto gridsize = gridInqSize(gridID);
  Varray<double> xvals(gridsize), yvals(gridsize);

  auto gridID2 = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID2)) cdo_abort("Target cell center coordinates missing!");

  gridInqXvals(gridID2, xvals.data());
  gridInqYvals(gridID2, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID2, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID2, CDI_YAXIS, yvals, "grid center lat");

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto jj = calc_index_jj(ny, yvals[i]);
    auto ii = calc_index_ii(nx, xvals[i]);
    array[i] = data[jj * nx + ii];
  }

  if (gridID != gridID2) gridDestroy(gridID2);
}

static void
remap_nn_reg2d_to_healpix(size_t nx, size_t ny, Varray<float> const &data, int gridID, Varray<float> &array)
{
  auto gridID2 = gridID;
  auto gridsize = gridInqSize(gridID2);
  auto hpParams = cdo::get_healpix_params(gridID2);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    double xval, yval;
    hp_index_to_lonlat(hpParams, i, &xval, &yval);
    auto jj = calc_index_jj(ny, yval * RAD2DEG);
    auto ii = calc_index_ii(nx, xval * RAD2DEG);
    array[i] = data[jj * nx + ii];
  }
}

static void
remap_nn_reg2d(size_t nx, size_t ny, Varray<float> const &data, int gridID, Varray<float> &array)
{
  auto gridtype = gridInqType(gridID);
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)
    remap_nn_reg2d_to_reg2d(nx, ny, data, gridID, array);
  else if (is_healpix_grid(gridID))
    remap_nn_reg2d_to_healpix(nx, ny, data, gridID, array);
  else
    remap_nn_reg2d_to_nonreg2d(nx, ny, data, gridID, array);
}

static int
define_point_grid()
{
  auto gridID = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID, 1);
  gridDefYsize(gridID, 1);
  double value = 0.0;
  gridDefXvals(gridID, &value);
  gridDefYvals(gridID, &value);

  return gridID;
}

static int
define_zaxis(bool lstdatm, int nlevels, double *levels)
{
  int zaxisID = -1;

  if (lstdatm)
  {
    zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevels);
    zaxisDefLevels(zaxisID, levels);
    cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, "level");
    cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, "Level");
    cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, "m");
  }
  else { zaxisID = zaxis_from_name("surface"); }

  return zaxisID;
}

static void
define_pressure_attributes(int vlistID, int varID)
{
  vlistDefVarParam(vlistID, varID, cdiEncodeParam(1, 255, 255));
  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "P");
  cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, "air_pressure");
  cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, "pressure");
  cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "hPa");
}

static void
define_temperature_attributes(int vlistID, int varID)
{
  vlistDefVarParam(vlistID, varID, cdiEncodeParam(130, 128, 255));
  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "T");
  cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, var_stdname(air_temperature));
  cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, "temperature");
  cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "K");
}

class Vargen : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vargen",
    // clang-format off
    .operators = { { "random", 0, 0, "grid description file or name, <seed>", VargenHelp },
                   { "const", 0, 0, "constant value, grid description file or name", VargenHelp },
                   { "sincos", 0, 0, "grid description file or name", VargenHelp },
                   { "coshill", 0, 0, "grid description file or name", VargenHelp },
                   { "testfield", 0, 0, "grid description file or name", VargenHelp },
                   { "seq", 0, 0, "start,end,<increment>", VargenHelp },
                   { "topo", VargenHelp },
                   { "temp", VargenHelp },
                   { "mask", VargenHelp },
                   { "stdatm", 0, 0, "height levels[m]", VargenHelp } },
    // clang-format on
    .aliases = { { "for", "seq" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 1, NoRestriction },
  };
  inline static RegisterEntry<Vargen> registration = RegisterEntry<Vargen>();

  int RANDOM, SINCOS, COSHILL, TESTFIELD, CONST, SEQ, TOPO, TEMP, MASK, STDATM;

private:
  static constexpr size_t nlat = 360, nlon = 720;
  double lon[nlon]{}, lat[nlat]{};
  int gridID = -1, gridIDdata = -1;
  double rstart = 0.0, rstop = 0.0, rinc = 0.0;
  double rconst = 0.0;
  std::vector<double> levels{};
  int numSteps{};

  int taxisID{};
  CdoStreamID streamID{};
  int varID2{};
  int vlistID{};
  int operatorID{};

public:
  void
  init() override
  {
    int nlevels = 1;

    RANDOM = module.get_id("random");
    SINCOS = module.get_id("sincos");
    COSHILL = module.get_id("coshill");
    // not used todo: make unavailable for non developers
    TESTFIELD = module.get_id("testfield");
    CONST = module.get_id("const");
    SEQ = module.get_id("seq");
    TOPO = module.get_id("topo");
    TEMP = module.get_id("temp");
    MASK = module.get_id("mask");
    STDATM = module.get_id("stdatm");

    operatorID = cdo_operator_id();

    if (operatorID == RANDOM) { gridID = random_init(operatorID); }
    else if (operatorID == SINCOS || operatorID == COSHILL || operatorID == TESTFIELD)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      operator_check_argc(1);
      gridID = cdo_define_grid(cdo_operator_argv(0));
    }
    else if (operatorID == CONST)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      operator_check_argc(2);
      rconst = parameter_to_double(cdo_operator_argv(0));
      gridID = cdo_define_grid(cdo_operator_argv(1));
    }
    else if (operatorID == TOPO || operatorID == TEMP || operatorID == MASK)
    {
      gridIDdata = gridCreate(GRID_LONLAT, nlon * nlat);
      gridDefXsize(gridIDdata, nlon);
      gridDefYsize(gridIDdata, nlat);

      for (size_t i = 0; i < nlon; ++i) lon[i] = -179.75 + i * 0.5;
      for (size_t i = 0; i < nlat; ++i) lat[i] = -89.75 + i * 0.5;

      gridDefXvals(gridIDdata, lon);
      gridDefYvals(gridIDdata, lat);

      gridID = gridIDdata;

      if (cdo_operator_argc() == 1) gridID = cdo_define_grid(cdo_operator_argv(0));
      if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");
    }
    else if (operatorID == SEQ)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      if (cdo_operator_argc() < 2) cdo_abort("Too few arguments!");
      if (cdo_operator_argc() > 3) cdo_abort("Too many arguments!");

      rstart = parameter_to_double(cdo_operator_argv(0));
      rstop = parameter_to_double(cdo_operator_argv(1));
      rinc = (cdo_operator_argc() == 3) ? parameter_to_double(cdo_operator_argv(2)) : 1;
      if (fp_is_equal(rinc, 0.0)) cdo_abort("Increment is zero!");

      gridID = define_point_grid();
    }
    else if (operatorID == STDATM)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      levels = cdo_argv_to_fltarr(cdo_get_oper_argv());
      nlevels = levels.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nlevels; ++i) printf("levels %d: %g\n", i, levels[i]);

      gridID = define_point_grid();
    }

    auto zaxisID = define_zaxis(operatorID == STDATM, nlevels, levels.data());

    vlistID = vlistCreate();

    auto timetype = (operatorID == SEQ) ? TIME_VARYING : TIME_CONSTANT;

    auto varID = vlistDefVar(vlistID, gridID, zaxisID, timetype);
    /*
       For the standard atmosphere two output variables are generated: pressure and temperature.
       The first (varID) is pressure, second (varID2) is temperature. Add an additional variable for the standard atmosphere.
    */
    varID2 = (operatorID == STDATM) ? vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT) : -1;

    if (operatorID == MASK) vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT8);

    if (operatorID == STDATM)
    {
      define_pressure_attributes(vlistID, varID);
      define_temperature_attributes(vlistID, varID2);
    }
    else
    {
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, cdo_operator_name(operatorID));
      if (operatorID == TOPO) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "m");
      if (operatorID == TEMP) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, "K");
    }

    taxisID = cdo_taxis_create(TAXIS_RELATIVE);
    vlistDefTaxis(vlistID, taxisID);

    if (operatorID != SEQ) vlistDefNtsteps(vlistID, 1);

    streamID = cdo_open_write(0);

    cdo_def_vlist(streamID, vlistID);

    numSteps = (operatorID == SEQ) ? 1.001 + ((rstop - rstart) / rinc) : 1;
    if (operatorID != SEQ) vlistDefNtsteps(vlistID, 0);
  }

  void
  run() override
  {
    VarList varList(vlistID);

    auto julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

    size_t gridsize = gridInqSize(gridID);
    Varray<float> array(gridsize);

    for (int tsID = 0; tsID < numSteps; ++tsID)
    {
      auto rval = rstart + rinc * tsID;
      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(julday_to_date(CALENDAR_PROLEPTIC, julday + tsID));
      taxisDefVdatetime(taxisID, vDateTime);
      cdo_def_timestep(streamID, tsID);

      // this should either be 1 or 2, two for atmosphere
      for (int varID = 0; varID < varList.numVars(); ++varID)
      {
        for (int levelID = 0; levelID < varList.vars[varID].nlevels; ++levelID)
        {
          cdo_def_field(streamID, varID, levelID);

          if (operatorID == RANDOM) { cdo::fill_random(array); }
          else if (operatorID == SINCOS || operatorID == COSHILL || operatorID == TESTFIELD)
          {
            Varray<double> xvals(gridsize), yvals(gridsize);

            if (grid_is_distance_generic(gridID)) { conv_generic_grid(gridID, gridsize, xvals, yvals); }
            else { gridID = generate_full_point_grid_radian(gridID, xvals, yvals); }

            // clang-format off
            if      (operatorID == SINCOS)    cdo::fill_sincos(array, xvals, yvals);
            else if (operatorID == COSHILL)   cdo::fill_coshill(array, xvals, yvals);
            else if (operatorID == TESTFIELD) cdo::fill_testfield(array, xvals, yvals);
            // clang-format on
          }
          else if (operatorID == CONST) { std::ranges::fill(array, rconst); }
          else if (operatorID == TOPO || operatorID == TEMP || operatorID == MASK)
          {
            Varray<float> data;

            // clang-format off
            if      (operatorID == TOPO) data = cdo::unpack_data(cdo::topoData);
            else if (operatorID == TEMP) data = cdo::unpack_data(cdo::tempData);
            else if (operatorID == MASK) data = cdo::unpack_data(cdo::maskData);
            // clang-format on

            if (gridID != gridIDdata && gridIDdata != -1) { remap_nn_reg2d(nlon, nlat, data, gridID, array); }
            else { array = data; }
          }
          else if (operatorID == SEQ) { array[0] = rval; }
          else if (operatorID == STDATM)
          {
            array[0] = (varID == varID2) ? cdo::std_atm_temperatur(levels[levelID]) : cdo::std_atm_pressure(levels[levelID]);
          }

          cdo_write_field_f(streamID, array.data(), 0);
        }
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);

    vlistDestroy(vlistID);
  }
};
