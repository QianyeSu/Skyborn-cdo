/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include <algorithm>

#include "cdo_options.h"
#include "cdo_omp.h"
#include "cdo_timer.h"
#include "cdo_zaxis.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "util_files.h"

static double
coshill(double x, double y)
{
  return 2.0 - std::cos(std::acos(std::cos(x) * std::cos(y)) / 1.2);
}

static double
testfield(double xval, double yval, double start, double shift)
{
  double xyz[3];
  gcLLtoXYZ(xval + shift, yval + 0.5 * shift, xyz);
  auto x = xyz[0];
  auto y = xyz[1];
  auto z = xyz[2];
  return start + 1.0 + std::pow(x, 8.0) + std::exp(2.0 * y * y * y) + std::exp(2.0 * x * x) + 10.0 * x * y * z;
}

static void
print_stat(const char *sinfo, MemType memtype, int datatype, int filetype, off_t nvalues, double dataSize, double fileSize,
           double tw)
{
  nvalues /= 1000000;
  dataSize /= 1024. * 1024. * 1024.;

  double rout = (tw > 0) ? nvalues / tw : -1;
  cdo_print("%s Wrote %.1f GB of %d bit floats to %s %s, %.1f MVal/s", sinfo, dataSize, (memtype == MemType::Float) ? 32 : 64,
            cdo::datatype_to_cstr(datatype), cdo::filetype_to_cstr(filetype), rout);

  fileSize /= 1024. * 1024. * 1024.;

  rout = (tw > 0) ? 1024 * fileSize / tw : -1;
  cdo_print("%s Wrote %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, fileSize, tw, rout);
}

static int
create_zaxis(int numLevels)
{
  int zaxisID = -1;

  if (numLevels == 1) { zaxisID = zaxis_from_name("surface"); }
  else
  {
    Varray<double> levels(numLevels);
    for (int i = 0; i < numLevels; ++i) levels[i] = 100 * i;
    zaxisID = zaxisCreate(ZAXIS_HEIGHT, numLevels);
    zaxisDefLevels(zaxisID, &levels[0]);
  }

  return zaxisID;
}

namespace
{
struct Parameter
{
  int nruns = 1;
  int nvars = 10;
  int nlevs = 0;
  int nsteps = 30;
  std::string grid = "global_.2";
  bool varySteps = false;
};
}  // namespace

static Parameter
get_parameter(void)
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    // kvlist.name = cdo_module_name();
    kvlist.name = "CDIwrite";
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "nruns")      params.nruns = parameter_to_int(value);
      else if (key == "nvars")      params.nvars = parameter_to_int(value);
      else if (key == "nlevs")      params.nlevs = parameter_to_int(value);
      else if (key == "nsteps")     params.nsteps = parameter_to_int(value);
      else if (key == "grid")       params.grid = parameter_to_word(value);
      else if (key == "varysteps")  params.varySteps = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
verify_parameter(Parameter &params)
{
  params.nruns = std::clamp(params.nruns, 0, 9999);
  params.nvars = std::max(params.nvars, 1);
  params.nlevs = std::clamp(params.nlevs, 1, 255);
  params.nsteps = std::max(params.nsteps, 1);
}

class CDIwrite : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CDIwrite",
    .operators = { { "cdiwrite" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 1, NoRestriction },
  };
  inline static RegisterEntry<CDIwrite> registration = RegisterEntry<CDIwrite>();

  MemType memtype = Options::CDO_Memtype;
  int filetype = -1, datatype = -1;
  char sinfo[64] = { 0 };
  off_t nvalues = 0;
  double fileSize = 0, dataSize = 0;
  double runTimeSum = 0.0;
  Parameter params;

  int vlistID{};
  int taxisID{};
  size_t gridsize{};

  Varray3D<double> vars;
  Varray<double> array;
  Varray<double> xvals, yvals;
  std::vector<float> farray;

public:
  void
  init() override
  {
    if (Options::cdoVerbose) cdo_print("parameter: nruns/nvars/nlevs/nsteps/grid/varysteps");

    params = get_parameter();
    verify_parameter(params);

    auto gridID = cdo_define_grid(params.grid);
    gridsize = gridInqSize(gridID);
    auto zaxisID = create_zaxis(params.nlevs);

    if (Options::cdoVerbose)
    {
      cdo_print("nruns     : %d", params.nruns);
      cdo_print("nvars     : %d", params.nvars);
      cdo_print("nlevs     : %d", params.nlevs);
      cdo_print("nsteps    : %d", params.nsteps);
      cdo_print("gridsize  : %zu", gridsize);
      cdo_print("varysteps : %d", params.varySteps);
    }

    array.resize(gridsize);
    xvals.resize(gridsize);
    yvals.resize(gridsize);

    auto gridID2 = generate_full_point_grid(gridID);
    if (!gridHasCoordinates(gridID2)) cdo_abort("Target cell center coordinates missing!");

    gridInqXvals(gridID2, &xvals[0]);
    gridInqYvals(gridID2, &yvals[0]);

    // Convert lat/lon units if required
    cdo_grid_to_radian(gridID2, CDI_XAXIS, xvals, "grid center lon");
    cdo_grid_to_radian(gridID2, CDI_YAXIS, yvals, "grid center lat");

    for (size_t i = 0; i < gridsize; ++i) array[i] = coshill(xvals[i], yvals[i]);

    vars.resize(params.nvars);
    for (int varID = 0; varID < params.nvars; ++varID)
    {
      vars[varID].resize(params.nlevs);
      for (int levelID = 0; levelID < params.nlevs; ++levelID)
      {
        auto &var = vars[varID][levelID];
        var.resize(gridsize);
        for (size_t i = 0; i < gridsize; ++i) { var[i] = varID + array[i] * (levelID + 1); }
      }
    }

    if (memtype == MemType::Float) farray.resize(gridsize);

    vlistID = vlistCreate();

    for (int i = 0; i < params.nvars; ++i)
    {
      auto varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(varID + 1, 255, 255));
    }

    taxisID = cdo_taxis_create(TAXIS_RELATIVE);
    vlistDefTaxis(vlistID, taxisID);

    vlistDefNtsteps(vlistID, params.nsteps);
  }

  void
  run() override
  {
    for (int irun = 0; irun < params.nruns; ++irun)
    {
      cdo::timer runTimer;
      dataSize = 0;
      nvalues = 0;

      auto streamID = cdo_open_write(0);
      cdo_def_vlist(streamID, vlistID);

      filetype = cdo_inq_filetype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);
      if (datatype == CDI_UNDEFID) datatype = CDI_DATATYPE_FLT32;

      auto julday = date_to_julday(CALENDAR_PROLEPTIC, 19870101);

      for (int tsID = 0; tsID < params.nsteps; ++tsID)
      {
        cdo::timer stepTimer;
        CdiDateTime vDateTime{};
        vDateTime.date = cdiDate_set(julday_to_date(CALENDAR_PROLEPTIC, julday + tsID));
        taxisDefVdatetime(taxisID, vDateTime);
        cdo_def_timestep(streamID, tsID);

        if (params.varySteps)
        {
          for (int varID = 0; varID < params.nvars; ++varID)
          {
            vars[varID].resize(params.nlevs);
            for (int levelID = 0; levelID < params.nlevs; ++levelID)
            {
              auto &var = vars[varID][levelID];
              var.resize(gridsize);
#ifdef _OPENMP
#pragma omp parallel for if (gridsize > cdoMinLoopSize) default(shared) schedule(static)
#endif
              for (size_t i = 0; i < gridsize; ++i)
              {
                var[i] = varID + testfield(xvals[i], yvals[i], 0.1 * tsID, 0.1 * tsID) * (levelID + 1);
              }
            }
          }
        }

        for (int varID = 0; varID < params.nvars; ++varID)
        {
          for (int levelID = 0; levelID < params.nlevs; ++levelID)
          {
            auto &var = vars[varID][levelID];
            nvalues += gridsize;
            cdo_def_field(streamID, varID, levelID);
            if (memtype == MemType::Float)
            {
              for (size_t i = 0; i < gridsize; ++i) { farray[i] = var[i]; }
              cdo_write_field_f(streamID, &farray[0], 0);
              dataSize += gridsize * 4;
            }
            else
            {
              cdo_write_field(streamID, &var[0], 0);
              dataSize += gridsize * 8;
            }
          }
        }

        if (Options::cdoVerbose) { cdo_print("Timestep %d: %.3f seconds", tsID + 1, stepTimer.elapsed()); }
      }

      cdo_stream_close(streamID);

      auto runTime = runTimer.elapsed();
      runTimeSum += runTime;

      fileSize = (double) FileUtils::size(cdo_get_stream_name(0));

      if (params.nruns > 1) std::snprintf(sinfo, sizeof(sinfo), "(run %d)", irun + 1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, dataSize, fileSize, runTime);
    }

    if (params.nruns > 1) print_stat("(mean)", memtype, datatype, filetype, nvalues, dataSize, fileSize, runTimeSum / params.nruns);
  }

  void
  close() override
  {
    vlistDestroy(vlistID);
  }
};
