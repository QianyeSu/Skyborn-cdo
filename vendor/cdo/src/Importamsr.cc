/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

#define NLON 1440
#define NLAT 720
#define MAX_VARS 6

static void
init_amsr_day(int vlistID, int gridID, int zaxisID, int numVars)
{
  /*
    Version-5 RSS AMSR-E or AMSR-J daily files

    filename  with path in form satname_yyyymmdd_v5.gz
    where satname  = name of satellite (amsre or amsr)
             yyyy	= year
               mm	= month
               dd	= day of month

    1:time	time of measurement in fractional hours GMT
    2:sst     	sea surface temperature in deg Celcius
    3:wind	10m surface wind in meters/second
    4:vapor	columnar water vapor in millimeters
    5:cloud	cloud liquid water in millimeters
    6:rain   	rain rate in millimeters/hour
  */
  const char *name[] = { "hours", "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "h", "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  constexpr double xscale[] = { 0.1, 0.15, 0.2, 0.3, 0.01, 0.1 };
  constexpr double xminval[] = { 0., -3., 0., 0., 0., 0. };

  for (int i = 0; i < numVars; ++i)
  {
    int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
    cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
    cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
    vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
    vlistDefVarMissval(vlistID, varID, 254);
    cdiDefKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, xscale[i]);
    cdiDefKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, xminval[i]);
  }
}

static void
init_amsr_averaged(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 AMSR-E or AMSR-J time-averaged files including:
          3-day		(average of 3 days ending on file date)
          weekly	(average of 7 days ending on Saturday of file date)
          monthly	(average of all days in month)


     filename
           format of file names are:
                3-day	satname_yyyymmddv5_d3d.gz
                weekly	satname_yyyymmddv5.gz
                monthly	satname_yyyymmv5.gz

        where	satname	=name of satellite (amsre or amsr)
                        yyyy	=year
                        mm	=month
                        dd	=day of month

    1:sst       sea surface temperature in deg Celcius
    2:wind      10m surface wind in meters/second
    3:vapor	columnar water vapor in millimeters
    4:cloud     cloud liquid water in millimeters
    5:rain	rain rate in millimeters/hour
  */
  const char *name[] = { "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  constexpr double xscale[] = { 0.15, 0.2, 0.3, 0.01, 0.1 };
  constexpr double xminval[] = { -3., 0., 0., 0., 0. };

  for (int i = 0; i < nvars; ++i)
  {
    // varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);
    int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
    cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name[i]);
    cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units[i]);
    vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
    vlistDefVarMissval(vlistID, varID, 254);
    cdiDefKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, xscale[i]);
    cdiDefKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, xminval[i]);
  }
}

static void
read_amsr(FILE *fp, VarList const &varList, Varray2D<double> &data, size_t *numMissVals)
{
  std::vector<unsigned char> amsr_data;

  for (auto const &var : varList.vars)
  {
    amsr_data.resize(var.gridsize);
    size_t size = std::fread(amsr_data.data(), 1, var.gridsize, fp);
    if (size != var.gridsize) cdo_abort("Read error!");

    double xminval = 0.0, xscale = 1.0;
    cdiInqKeyFloat(varList.vlistID, var.ID, CDI_KEY_SCALEFACTOR, &xscale);
    cdiInqKeyFloat(varList.vlistID, var.ID, CDI_KEY_ADDOFFSET, &xminval);

    numMissVals[var.ID] = 0;
    for (size_t i = 0; i < var.gridsize; ++i)
    {
      if (amsr_data[i] <= 250) { data[var.ID][i] = amsr_data[i] * xscale + xminval; }
      else
      {
        data[var.ID][i] = var.missval;
        numMissVals[var.ID]++;
      }
    }
  }
}

static void
write_data(CdoStreamID streamID, int nvars, Varray2D<double> &data, size_t *numMissVals)
{
  for (int varID = 0; varID < nvars; ++varID)
  {
    cdo_def_field(streamID, varID, 0);
    cdo_write_field(streamID, data[varID].data(), numMissVals[varID]);
  }
}

static int
get_date(const char *name)
{
  const char *pname = strchr(name, '_');
  int date = pname ? atoi(pname + 1) : 0;
  return date;
}

class Importamsr : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Importamsr",
    .operators = { { "import_amsr", ImportamsrHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Importamsr> registration = RegisterEntry<Importamsr>(module);

  int tsID{};
  int numVars{};
  int vtime = 0;
  double xvals[NLON], yvals[NLAT];
  Varray2D<double> data;
  size_t numMissVals[MAX_VARS];

  CdoStreamID streamID{};

  int vlistID{};
  int gridID{};
  int zaxisID{};
  int taxisID{};

  size_t gridsize{};
  size_t fsize{};

  FILE *fp{};
  int vdate{};

  void
  A()
  {
    numVars = 6;
    for (int i = 0; i < numVars; ++i) data[i].resize(gridsize);

    init_amsr_day(vlistID, gridID, zaxisID, numVars);
    VarList varList(vlistID);

    cdo_def_vlist(streamID, vlistID);

    vtime = 13000;  // 1:30:00
    for (tsID = 0; tsID < 2; ++tsID)
    {
      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(vdate);
      vDateTime.time = cdiTime_set(vtime);
      taxisDefVdatetime(taxisID, vDateTime);
      vtime += 120000;  // 13:30:00
      cdo_def_timestep(streamID, tsID);

      read_amsr(fp, varList, data, numMissVals);

      write_data(streamID, numVars, data, numMissVals);
    }
  }

  void
  B()
  {
    numVars = 5;
    for (int i = 0; i < numVars; ++i) data[i].resize(gridsize);

    init_amsr_averaged(vlistID, gridID, zaxisID, numVars);
    VarList varList(vlistID);

    // vlistDefNtsteps(vlistID, 0);
    cdo_def_vlist(streamID, vlistID);

    CdiDateTime vDateTime{};
    vDateTime.date = cdiDate_set(vdate);
    vDateTime.time = cdiTime_set(vtime);
    taxisDefVdatetime(taxisID, vDateTime);
    tsID = 0;
    cdo_def_timestep(streamID, tsID);

    read_amsr(fp, varList, data, numMissVals);

    write_data(streamID, numVars, data, numMissVals);
  }

public:
  void
  init() override
  {
    operator_check_argc(0);

    fp = std::fopen(cdo_get_stream_name(0), "r");
    if (fp == nullptr)
    {
      perror(cdo_get_stream_name(0));
      exit(EXIT_FAILURE);
    }

    fseek(fp, 0L, SEEK_END);
    fsize = (size_t) ftell(fp);
    fseek(fp, 0L, SEEK_SET);

    vdate = get_date(cdo_get_stream_name(0));
    if (vdate <= 999999) vdate = vdate * 100 + 1;

    streamID = cdo_open_write(1);

    /*
      Longitude  is 0.25*xdim-0.125    degrees east
      Latitude   is 0.25*ydim-90.125
    */
    gridsize = NLON * NLAT;
    gridID = gridCreate(GRID_LONLAT, gridsize);
    gridDefXsize(gridID, NLON);
    gridDefYsize(gridID, NLAT);
    for (int i = 0; i < NLON; ++i) xvals[i] = 0.25 * (i + 1) - 0.125;
    for (int i = 0; i < NLAT; ++i) yvals[i] = 0.25 * (i + 1) - 90.125;
    gridDefXvals(gridID, xvals);
    gridDefYvals(gridID, yvals);

    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

    vlistID = vlistCreate();
    taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID, taxisID);

    data = Varray2D<double>(MAX_VARS);
  }

  void
  run() override
  {
    if (fsize == 12441600) { A(); }
    else if (fsize == 5184000) { B(); }
    else { cdo_abort("Unexpected file size for AMSR data!"); }
    process_def_var_num(vlistNvars(vlistID));
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
    vlistDestroy(vlistID);

    std::fclose(fp);
  }
};
