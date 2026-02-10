/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Brockmann Consult

  Author: Ralf Quast

*/

#include <assert.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <cdi.h>

#include "cdo_cdi_wrapper.h"
#include "varray.h"

#define TO_KELVIN(x) ((x) + 273.15)
#define MISSVAL -9.0E33

static int
equals(double expected, double actual, double eps)
{
  return (int) (std::fabs(expected - actual) < eps);
}
/*
static
double humidityIndex(double t, double p, double r, double missval)
{
  static double tmin = 26.0;
  static double rmin = 40.0;

  if ( t < tmin || r < rmin )
    return missval;

  return t + (5.0 / 9.0) * ((0.01 * r * p * 6.112 * std::pow(10.0, (7.5 * t) / (237.7 t))) - 10.0);
}
*/
// reads a file containing data for a single grid point
static void
read_file(const char path[], Varray2D<double> &vars, int numVars, int numSteps)
{
  auto streamID = streamOpenRead(path);
  if (streamID < 0)
  {
    fprintf(stderr, "%s\n", cdiStringError(streamID));
    exit(EXIT_FAILURE);
  }

  auto vlistID = streamInqVlist(streamID);

  assert(vlistNtsteps(vlistID) == numSteps);
  assert(cdo_vlist_gridsizemax(vlistID) == 1);
  assert(vlistNvars(vlistID) == numVars);

  // auto taxisID = vlistInqTaxis(vlistID);

  for (int tsID = 0; tsID < numSteps; ++tsID)
  {
    auto numFields = streamInqTimestep(streamID, tsID);
    assert(numFields == numVars);

    // auto vDateTime = taxisInqVdatetime(taxisID);

    size_t numMissVals;
    for (int varID = 0; varID < numVars; ++varID) streamReadVar(streamID, varID, &vars[varID][tsID], &numMissVals);
  }

  streamClose(streamID);
}

// write file containing data for a single grid point
static void
write_file(const char path[], const double array[], int length)
{
  double lons[] = { 0.0 };
  double lats[] = { 0.0 };

  auto gridID = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID, 1);
  gridDefYsize(gridID, 1);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);

  auto zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  auto vlistID = vlistCreate();

  auto varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "test_values");
  vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_FLT64);
  vlistDefVarMissval(vlistID, varID, MISSVAL);

  auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  auto streamID = streamOpenWrite(path, CDI_FILETYPE_SRV);
  if (streamID < 0)
  {
    fprintf(stderr, "%s\n", cdiStringError(streamID));
    exit(EXIT_FAILURE);
  }

  streamDefVlist(streamID, vlistID);

  for (int tsID = 0; tsID < length; ++tsID)
  {
    int vdate = (tsID < 6) ? (20060625 + tsID) : (20060701 + tsID - 6);
    int vtime = 235900;
    CdiDateTime vDateTime{};
    vDateTime.date = cdiDate_set(vdate);
    vDateTime.time = cdiTime_set(vtime);
    taxisDefVdatetime(taxisID, vDateTime);
    streamDefTimestep(streamID, tsID);

    auto value = array[tsID];
    size_t numMissVals = fp_is_equal(value, MISSVAL) ? 1 : 0;

    streamWriteVar(streamID, varID, &value, numMissVals);
  }

  streamClose(streamID);

  vlistDestroy(vlistID);
}

// gets the path of the CDO binary executable
static std::string
get_cdo_path()
{
  char *cdoPath = getenv("CDO_PATH");
  if (cdoPath == nullptr)
  {
    struct stat filestat;
    int status;
    status = stat("$HOME/bin/cdo", &filestat);
    if (status == 0)
      return "$HOME/bin/cdo";
    else
    {
      fprintf(stderr, "cdo binary not found! Use CDO_PATH to set the path to cdo.\n");
      exit(-1);
    }
  }

  return cdoPath;
}

// submits a CDO command
static int
submit_cdo_command(const char *argument)
{
  auto cdoPath = get_cdo_path();
  auto cdoCommand = cdoPath + " -b 64 " + argument;
  printf("%s\n", cdoCommand.c_str());

  return system(cdoCommand.c_str());
}

static void
testEcaFd()
{
  constexpr double array[] = { MISSVAL, MISSVAL, TO_KELVIN(1.0), TO_KELVIN(1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0) };

  int nvars = 1;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in.srv", array, 2);

  submit_cdo_command("eca_fd in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));

  write_file("in.srv", array, 3);

  submit_cdo_command("eca_fd in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 0.0));

  write_file("in.srv", array, 5);

  submit_cdo_command("eca_fd in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.));

  write_file("in.srv", array, 6);

  submit_cdo_command("eca_fd in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 2.));
}

static void
testEcaSu()
{
  constexpr double array[] = { MISSVAL, MISSVAL, TO_KELVIN(26.0), TO_KELVIN(24.0), TO_KELVIN(26.0), TO_KELVIN(24.0) };

  int nvars = 1;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in.srv", array, 2);

  submit_cdo_command("eca_su in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));

  write_file("in.srv", array, 6);

  submit_cdo_command("eca_su in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 2.));

  submit_cdo_command("eca_su,20.0 in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 4.));

  submit_cdo_command("eca_su,30.0 in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 0.));
}

static void
testFdns()
{
  constexpr double array1[] = { MISSVAL, TO_KELVIN(1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0) };
  constexpr double array2[] = { 0.0, 0.0, 1.0, 0.0, MISSVAL };

  int nvars = 1;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in1.srv", array1, 1);
  write_file("in2.srv", array2, 1);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));

  write_file("in1.srv", array1, 2);
  write_file("in2.srv", array2, 2);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 0.));

  write_file("in1.srv", array1, 3);
  write_file("in2.srv", array2, 3);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 0.));

  write_file("in1.srv", array1, 4);
  write_file("in2.srv", array2, 4);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.));

  write_file("in1.srv", array1, 5);
  write_file("in2.srv", array2, 5);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.));

  write_file("in1.srv", array1 + 4, 1);
  write_file("in2.srv", array2 + 4, 1);

  submit_cdo_command("fdns in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));
}

static void
testEcaGsl()
{
  constexpr double array1[]
      = { TO_KELVIN(6.0),  TO_KELVIN(6.0),  TO_KELVIN(6.0),  TO_KELVIN(6.0),  TO_KELVIN(6.0),  TO_KELVIN(6.0),  TO_KELVIN(6.0),
          TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0), TO_KELVIN(-1.0) };
  constexpr double array2[] = { 0.5 };

  int nvars = 2;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in1.srv", array1, 14);
  write_file("in2.srv", array2, 1);

  submit_cdo_command("eca_gsl in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 7.0));
  assert(fp_is_equal(vars[1][0], 181.0));

  //  submit_cdo_command("eca_gsl,6,5.0,0.6 in1.srv in2.srv out.srv");
  //  read_file("out.srv", vars, nvars, nts);
  //  assert(fp_is_equal(vars[0][0], MISSVAL));
  //  assert(fp_is_equal(vars[1][0], MISSVAL));

  write_file("in1.srv", array1, 7);

  submit_cdo_command("eca_gsl in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.0));
  assert(fp_is_equal(vars[1][0], 181.0));

  write_file("in1.srv", array1, 8);

  submit_cdo_command("eca_gsl in1.srv in2.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 2.0));
  assert(fp_is_equal(vars[1][0], 181.0));

  //  write_file("in1.srv", array1, 4);

  //  submit_cdo_command("eca_gsl in1.srv in2.srv out.srv");
  //  read_file("out.srv", vars, nvars, nts);
  //  assert(fp_is_equal(vars[0][0], MISSVAL));
  //  assert(fp_is_equal(vars[1][0], MISSVAL));
}
/*
static
void testHi()
{
  constexpr double array1[] = {MISSVAL, 70.0, 36.0, 46.0, 56.0};
  constexpr double array2[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  constexpr double array3[] = {0.4, 0.4, 0.3, 0.5, 0.6};

  int nvars = 1;
  int nts   = 5;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in1.srv", array1, 5);
  write_file("in2.srv", array2, 5);
  write_file("in3.srv", array3, 5);

  submit_cdo_command("hi in1.srv in2.srv in3.srv out.srv");

  read_file("out.srv", vars, nvars, nts);
  assert(equals(vars[0][0], humidityIndex(array1[0], array2[0], array3[0], MISSVAL), 1.0e-5));
  assert(equals(vars[0][1], humidityIndex(array1[1], array2[1], array3[1], MISSVAL), 1.0e-5));
  assert(equals(vars[0][2], humidityIndex(array1[2], array2[2], array3[2], MISSVAL), 1.0e-5));
  assert(equals(vars[0][3], humidityIndex(array1[3], array2[3], array3[3], MISSVAL), 1.0e-5));
  assert(equals(vars[0][4], humidityIndex(array1[4], array2[4], array3[4], MISSVAL), 1.0e-5));
}
*/
static void
testTimcount()
{
  constexpr double array[] = { MISSVAL, MISSVAL, TO_KELVIN(1.0), MISSVAL, TO_KELVIN(1.0), TO_KELVIN(1.0) };

  // number of output variables and time steps
  int nvars = 1;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in.srv", array, 2);

  submit_cdo_command("timcount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));

  write_file("in.srv", array, 3);

  submit_cdo_command("timcount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.));

  write_file("in.srv", array, 5);

  submit_cdo_command("timcount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 2.));

  write_file("in.srv", array, 6);

  submit_cdo_command("timcount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 3.));
}
/*
static
void testSeascount()
{
  constexpr double array[] = {MISSVAL, MISSVAL, TO_KELVIN(1.0), MISSVAL,
    TO_KELVIN(1.0), TO_KELVIN(1.0)};

  // number of output variables and time steps
  int nvars = 1;
  int nts   = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in.srv", array, 2);

  submit_cdo_command("seascount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], MISSVAL));

  write_file("in.srv", array, 3);

  submit_cdo_command("seascount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 1.));

  write_file("in.srv", array, 5);

  submit_cdo_command("seascount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 2.));

  write_file("in.srv", array, 6);

  submit_cdo_command("seascount in.srv out.srv");
  read_file("out.srv", vars, nvars, nts);
  assert(fp_is_equal(vars[0][0], 3.));
}
*/
static void
testWct()
{
  constexpr double array1[] = { -3.1102 };
  constexpr double array2[] = { 1.9787 };

  int nvars = 1;
  int nts = 1;

  Varray2D<double> vars(nvars, Varray<double>(nts));

  write_file("in1.srv", array1, 1);
  write_file("in2.srv", array2, 1);

  submit_cdo_command("wct in1.srv in2.srv out.srv");

  read_file("out.srv", vars, nvars, nts);
  assert(equals(vars[0][0], -6.34597, 1.0e-5));
}

int
main(void)
{
  testEcaFd();
  testEcaSu();
  testEcaGsl();

  testFdns();
  // testHi();
  testTimcount();
  testWct();

  return 0;
}
