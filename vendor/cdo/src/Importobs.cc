/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <fstream>

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "varray.h"
#include <mpim_grid.h>
#include "griddes.h"

static void
init_vars(int vlistID, int gridID, int zaxisID, int nvars)
{
  const int code[] = { 11, 17, 33, 34, 1, 2 /*, 3*/ };
  const char *name[] = { "temp", "depoint", "u", "v", "height", "pressure" /*, "station"*/ };
  const char *units[] = { "Celsius", "", "m/s", "m/s", "m", "hPa" /*, ""*/ };

  for (int i = 0; i < nvars; ++i)
  {
    auto varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
    vlistDefVarParam(vlistID, varID, cdiEncodeParam(code[i], 255, 255));
    cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
    cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
    vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_FLT32);
  }
}

static void
init_data(VarList const &varList, Varray2D<double> &data)
{
  for (auto const &var : varList.vars)
  {
    for (size_t i = 0; i < var.gridsize; ++i) data[var.ID][i] = var.missval;
  }
}

static void
write_data(CdoStreamID streamID, VarList const &varList, Varray2D<double> &data)
{
  for (auto const &var : varList.vars)
  {
    auto numMissVals = varray_num_mv(var.gridsize, data[var.ID], var.missval);
    cdo_def_field(streamID, var.ID, 0);
    cdo_write_field(streamID, data[var.ID].data(), numMissVals);
  }
}

static int
get_date(const char *name)
{
  const char *pname = strchr(name, '_');
  int date = pname ? atoi(pname + 1) : 0;
  return date;
}

class Importobs : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Importobs",
    .operators = { { "import_obs", 0, 0, "grid description file or name" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Importobs>(module);

  int numVars = 6;
  int vtime = 0;
  char dummy[32], station[32], datetime[32];
  float lat{}, lon{}, height1{}, pressure{}, height2{}, value{};
  double latmin = 90, latmax = -90, lonmin = 360, lonmax = -360;
  int code{};

  CdoStreamID streamID{};

  int taxisID{};
  int vlistID{};
  int zaxisID{};
  int gridID{};

  int vdate{};

  Varray2D<double> data;

  Varray<double> xvals;
  size_t xsize{};
  Varray<double> yvals;
  size_t ysize{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    gridID = cdo_define_grid(cdo_operator_argv(0));

    if (gridInqType(gridID) != GRID_LONLAT) cdo_abort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID)));

    auto gridsize = gridInqSize(gridID);
    xsize = gridInqXsize(gridID);
    ysize = gridInqYsize(gridID);

    xvals = Varray<double>(gridsize);
    yvals = Varray<double>(gridsize);
    gridInqXvals(gridID, xvals.data());
    gridInqYvals(gridID, yvals.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
    cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

    vdate = get_date(cdo_get_stream_name(0));
    if (vdate <= 999999) vdate = vdate * 100 + 1;

    streamID = cdo_open_write(1);
    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
    vlistID = vlistCreate();

    taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID, taxisID);
    data.resize(numVars);
    for (int i = 0; i < numVars; ++i) data[i].resize(gridsize);

    init_vars(vlistID, gridID, zaxisID, numVars);

    cdo_def_vlist(streamID, vlistID);
  }

  void
  run() override
  {
    VarList varList(vlistID);

    std::ifstream file(cdo_get_stream_name(0));
    if (!file.is_open()) cdo_abort("Open failed on: %s\n", cdo_get_stream_name(0));

    int vdate0 = 0;
    int vtime0 = 0;
    // ntime = 0;
    int tsID = 0;
    std::string line;
    while (std::getline(file, line))
    {
      std::sscanf(line.c_str(), "%s %s %s %g %g %g %d %g %g %g", dummy, station, datetime, &lat, &lon, &height1, &code, &pressure,
                  &height2, &value);
      long vdate_l;
      std::sscanf(datetime, "%ld_%d", &vdate_l, &vtime);
      vdate = vdate_l;

      if (vdate != vdate0 || vtime != vtime0)
      {
        if (tsID > 0) write_data(streamID, varList, data);

        vdate0 = vdate;
        vtime0 = vtime;
        // printf("%s %d %d %g %g %g %d %g %g %g\n", station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);
        CdiDateTime vDateTime{};
        vDateTime.date = cdiDate_set(vdate);
        vDateTime.time = cdiTime_set(vtime);
        taxisDefVdatetime(taxisID, vDateTime);
        cdo_def_timestep(streamID, tsID);

        init_data(varList, data);

        tsID++;
      }

      if (lon < lonmin) lonmin = lon;
      if (lon > lonmax) lonmax = lon;
      if (lat < latmin) latmin = lat;
      if (lat > latmax) latmax = lat;

      auto dy = yvals[1] - yvals[0];
      size_t j;
      for (j = 0; j < ysize; ++j)
        if (lat >= (yvals[j] - dy / 2) && lat < (yvals[j] + dy / 2)) break;

      auto dx = xvals[1] - xvals[0];
      if (lon < (xvals[0] - dx / 2) && lon < 0) lon += 360;
      size_t i;
      for (i = 0; i < xsize; ++i)
        if (lon >= (xvals[i] - dx / 2) && lon < (xvals[i] + dx / 2)) break;

      long index = -1;
      if (code == 11) index = 0;
      if (code == 17) index = 1;
      if (code == 33) index = 2;
      if (code == 34) index = 3;

      // printf("%d %d %d %g %g %g %g\n", i, j, index, dx, dy, lon, lat);
      if (i < xsize && j < ysize && index >= 0)
      {
        char *pstation = station;
        while (std::isalpha(*pstation)) pstation++;
        // printf("station %s %d\n", pstation, atoi(pstation));
        data[index][j * xsize + i] = value;
        data[4][j * xsize + i] = height1;
        data[5][j * xsize + i] = pressure;
        // data[    6][j*xsize+i] = atoi(pstation);
      }

      /*
        printf("%s %d %d %g %g %g %d %g %g %g\n", station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);
      */
    }

    file.close();

    write_data(streamID, varList, data);

    if (Options::cdoVerbose) printf("lonmin=%g, lonmax=%g, latmin=%g, latmax=%g\n", lonmin, lonmax, latmin, latmax);

    process_def_var_num(vlistNvars(vlistID));
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
    vlistDestroy(vlistID);
  }
};
