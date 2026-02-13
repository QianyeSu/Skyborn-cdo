/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
      Output     outputtab       Table output
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "printinfo.h"
#include "cdo_zaxis.h"

static void
outputarr(int dig, size_t gridsize, Varray<double> const &array)
{
  for (size_t i = 0; i < gridsize; ++i) { std::fprintf(stdout, "  arr[%zu] = %.*g;\n", i, dig, array[i]); }
}

static void
outputsp(size_t gridsize, Varray<double> const &array, long ntr)
{
  auto mm = varray_min_max(gridsize, array);
  if (/* T11 */ mm.min >= -1 && mm.max <= 12)
  {
    auto spc = array.data();
    for (long m = 0; m <= ntr; ++m)
    {
      for (long n = m; n <= ntr; ++n)
      {
        std::fprintf(stdout, "%3d", (int) *spc++);
        std::fprintf(stdout, "%3d", (int) *spc++);
      }
      std::fprintf(stdout, "\n");
    }
  }
}

static void
output(size_t gridsize, Varray<double> const &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (nout == 6)
    {
      nout = 0;
      std::fprintf(stdout, "\n");
    }
    std::fprintf(stdout, " %12.6g", array[i]);
    nout++;
  }
  std::fprintf(stdout, "\n");
}

static void
print_xyz(std::FILE *fp, double x, double y, double z)
{
  std::fprintf(fp, "%g %g %g\n", x, y, z);
}

static void
outputxyz(size_t gridsize, Varray<double> const &array, double missval, size_t nlon, size_t nlat, Varray<double> const &lon,
          Varray<double> const &lat)
{
  double fmin = 0.0;
  for (size_t i = 0; i < gridsize; ++i)
    if (fp_is_not_equal(array[i], missval))
    {
      if (array[i] < fmin) fmin = array[i];
      std::fprintf(stdout, "%g\t%g\t%g\t%g\n", lon[i], lat[i], array[i], array[i]);
    }
  auto fname = "frontplane.xyz";
  auto fobj = c_fopen(fname, "w");
  if (fobj.get() == nullptr) cdo_abort("Open failed on %s", fname);
  // first front plane
  auto dx = (lon[1] - lon[0]);
  auto x0 = lon[0] - dx / 2;
  auto y0 = lat[0] - dx / 2;
  auto z0 = fmin;
  std::fprintf(fobj.get(), ">\n");
  for (size_t i = 0; i < nlon; ++i)
  {
    print_xyz(fobj.get(), x0, y0, z0);
    z0 = array[i];
    print_xyz(fobj.get(), x0, y0, z0);
    x0 = x0 + dx;
    print_xyz(fobj.get(), x0, y0, z0);
  }
  z0 = fmin;
  print_xyz(fobj.get(), x0, y0, z0);
  x0 = lon[0] - dx / 2;
  print_xyz(fobj.get(), x0, y0, z0);

  // second front plane
  x0 = lon[0] - dx / 2;
  y0 = lat[0] - dx / 2;
  z0 = fmin;
  std::fprintf(fobj.get(), ">\n");
  for (size_t i = 0; i < nlat; ++i)
  {
    print_xyz(fobj.get(), x0, y0, z0);
    z0 = array[i * nlon];
    print_xyz(fobj.get(), x0, y0, z0);
    y0 += dx;
    print_xyz(fobj.get(), x0, y0, z0);
  }
  z0 = fmin;
  print_xyz(fobj.get(), x0, y0, z0);
  y0 = lat[0] - dx / 2;
  print_xyz(fobj.get(), x0, y0, z0);
}

static void
read_xy_coordinates(bool hasRegxyCoordinates, int gridID0, Varray<double> &grid_xvals, Varray<double> &grid_yvals)
{
  auto gridsize = gridInqSize(gridID0);
  auto xsize = gridInqXsize(gridID0);
  auto ysize = gridInqYsize(gridID0);

  if (hasRegxyCoordinates)
  {
    grid_xvals.resize(xsize);
    grid_yvals.resize(ysize);
    gridInqXvals(gridID0, grid_xvals.data());
    gridInqYvals(gridID0, grid_yvals.data());
  }
  else
  {
    auto gridIDx = generate_full_point_grid(gridID0);
    if (!gridHasCoordinates(gridIDx)) cdo_abort("Cell center coordinates missing!");

    grid_xvals.resize(gridsize);
    grid_yvals.resize(gridsize);
    gridInqXvals(gridIDx, grid_xvals.data());
    gridInqYvals(gridIDx, grid_yvals.data());

    if (gridIDx != gridID0) gridDestroy(gridIDx);
  }
}

static void
read_lonlat_coordinates(int gridID0, Varray<double> &grid_center_lon, Varray<double> &grid_center_lat)
{
  auto gridsize = gridInqSize(gridID0);
  auto gridIDx = generate_full_point_grid(gridID0);
  if (!gridHasCoordinates(gridIDx)) cdo_abort("Cell center coordinates missing!");

  grid_center_lon.resize(gridsize);
  grid_center_lat.resize(gridsize);
  gridInqXvals(gridIDx, grid_center_lon.data());
  gridInqYvals(gridIDx, grid_center_lat.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridIDx, CDI_XAXIS, grid_center_lon, "grid center lon");
  cdo_grid_to_degree(gridIDx, CDI_YAXIS, grid_center_lat, "grid center lat");

  if (gridIDx != gridID0) gridDestroy(gridIDx);
}

static void
outputint(size_t gridsize, Varray<double> const &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (nout == 8)
    {
      nout = 0;
      std::fprintf(stdout, "\n");
    }
    std::fprintf(stdout, " %8d", (int) array[i]);
    nout++;
  }
  std::fprintf(stdout, "\n");
}

static void
outputf(int nelem, std::string const &format, size_t gridsize, Varray<double> const &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (nout == nelem)
    {
      nout = 0;
      std::fprintf(stdout, "\n");
    }
    std::fprintf(stdout, format.c_str(), array[i]);
    nout++;
  }
  std::fprintf(stdout, "\n");
}

class Output : public Process
{
  enum
  {
    knohead,
    kvalue,
    kparam,
    kcode,
    kname,
    kx,
    ky,
    klon,
    klat,
    klev,
    kbin,
    kxind,
    kyind,
    ktimestep,
    kdate,
    ktime,
    kyear,
    kmonth,
    kday
  };
  struct KeyLenEntry
  {
    std::string key;
    int idx;
    int len;
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Output",
    .operators = { { "output", 0, 1, OutputHelp },
                   { "outputint", OutputHelp },
                   { "outputsrv", OutputHelp },
                   { "outputext", OutputHelp },
                   { "outputf", OutputHelp },
                   { "outputts", OutputHelp },
                   { "outputfld", OutputHelp },
                   { "outputarr", OutputHelp },
                   { "outputxyz", OutputHelp },
                   { "outputtab", OutputtabHelp } },
    .aliases = { { "outputkey", "outputtab" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { -1, 0, NoRestriction },
  };
  inline static RegisterEntry<Output> registration = RegisterEntry<Output>();

  int OUTPUT{}, OUTPUTINT{}, OUTPUTSRV{}, OUTPUTEXT{}, OUTPUTF{}, OUTPUTTS{}, OUTPUTFLD{}, OUTPUTARR{}, OUTPUTXYZ{}, OUTPUTTAB{};
  size_t numMissVals{};
  int nelem = 1;
  int index{};
  std::string format{};
  char paramstr[32]{};
  int year{}, month{}, day{};
  std::vector<int> keyIndices;
  bool opercplx{};
  int operatorID{};

  // clang-format off
  std::vector<KeyLenEntry> keyMap = {
    { "nohead",   knohead,    0 },
    { "value",    kvalue,     8 },
    { "param",    kparam,    11 },
    { "code",     kcode,      4 },
    { "name",     kname,      8 },
    { "x",        kx,         6 },
    { "y",        ky,         6 },
    { "lon",      klon,       6 },
    { "lat",      klat,       6 },
    { "lev",      klev,       6 },
    { "bin",      kbin,       6 },
    { "xind",     kxind,      4 },
    { "yind",     kyind,      4 },
    { "timestep", ktimestep,  6 },
    { "date",     kdate,     10 },
    { "time",     ktime,      8 },
    { "year",     kyear,      5 },
    { "month",    kmonth,     2 },
    { "day",      kday,       2 },
  };
  // clang-format on

public:
  void
  init() override
  {
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }

    OUTPUT = module.get_id("output");
    OUTPUTINT = module.get_id("outputint");
    OUTPUTSRV = module.get_id("outputsrv");
    OUTPUTEXT = module.get_id("outputext");
    OUTPUTF = module.get_id("outputf");
    OUTPUTTS = module.get_id("outputts");
    OUTPUTFLD = module.get_id("outputfld");
    OUTPUTARR = module.get_id("outputarr");
    OUTPUTXYZ = module.get_id("outputxyz");
    OUTPUTTAB = module.get_id("outputtab");

    (void) (OUTPUT);  // unused

    operatorID = cdo_operator_id();
    opercplx = (cdo_operator_f2(operatorID) == 1);

    if (operatorID == OUTPUTF)
    {
      operator_input_arg("format and number of elements [optional]");

      if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

      format = cdo_operator_argv(0);
      if (cdo_operator_argc() == 2) nelem = parameter_to_int(cdo_operator_argv(1));
    }
    else if (operatorID == OUTPUTTAB)
    {
      auto lhead = true;

      operator_input_arg("keys to print");

      auto numArgs = cdo_operator_argc();
      auto const &argList = cdo_get_oper_argv();

      if (Options::cdoVerbose)
        for (int i = 0; i < numArgs; ++i) cdo_print("key%d=%s", i + 1, argList[i]);

      keyIndices.reserve(numArgs);
      for (int i = 0; i < numArgs; ++i)
      {
        auto sz = argList[i].find_first_of(':');
        auto currentName = (sz == std::string::npos) ? argList[i] : argList[i].substr(0, sz);
        auto formatLen = (sz == std::string::npos) ? -1 : std::stoi(argList[i].substr(sz + 1));
        size_t k;
        for (k = 0; k < keyMap.size(); ++k)
        {
          auto const &key = keyMap[k].key;
          if (key == currentName)
          {
            if (keyMap[k].idx == knohead)
              lhead = false;
            else
            {
              keyIndices.push_back(k);
              if (formatLen != -1) keyMap[k].len = formatLen;
            }
            break;
          }
        }

        if (k == keyMap.size()) cdo_abort("Key >%s< unsupported!", currentName);
      }

      if (Options::cdoVerbose)
        for (auto ki : keyIndices) cdo_print("idx=%d/%d  len=%d  name=%s", ki, keyMap[ki].idx, keyMap[ki].len, keyMap[ki].key);

      if (lhead)
      {
        std::fprintf(stdout, "#");
        for (auto ki : keyIndices) std::fprintf(stdout, "%*s ", keyMap[ki].len, keyMap[ki].key.c_str());
        std::fprintf(stdout, "\n");
      }
    }
    else { operator_check_argc(0); }
  }

  void
  run() override
  {
    for (int fileIdx = 0; fileIdx < cdo_stream_cnt(); fileIdx++)
    {
      auto streamID = cdo_open_read(fileIdx);
      auto vlistID = cdo_stream_inq_vlist(streamID);

      VarList varList(vlistID);

      auto numGrids = vlistNumGrids(vlistID);
      int numDiffGrids = 0;
      for (index = 1; index < numGrids; ++index)
        if (vlistGrid(vlistID, 0) != vlistGrid(vlistID, index)) numDiffGrids++;

      if (numDiffGrids > 0) cdo_abort("Too many different grids!");

      auto gridID0 = vlistGrid(vlistID, 0);
      auto gridtype = gridInqType(gridID0);
      auto hasRegxyCoordinates = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION);
      auto gridsize = gridInqSize(gridID0);
      auto xsize = gridInqXsize(gridID0);
      size_t nwpv = (vlistNumber(vlistID) == CDI_COMP) ? 2 : 1;
      if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");
      auto gridsizemax = nwpv * gridsize;
      Varray<double> array(gridsizemax);
      Varray<double> grid_center_lon, grid_center_lat;
      Varray<double> grid_xvals, grid_yvals;

      if (operatorID == OUTPUTTAB) read_xy_coordinates(hasRegxyCoordinates, gridID0, grid_xvals, grid_yvals);

      if (operatorID == OUTPUTFLD || operatorID == OUTPUTXYZ || operatorID == OUTPUTTAB)
        read_lonlat_coordinates(gridID0, grid_center_lon, grid_center_lat);

      int tsID = 0;
      auto taxisID = vlistInqTaxis(vlistID);
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID);
        auto vDateStr = date_to_string(vDateTime.date);
        auto vTimeStr = time_to_string(vDateTime.time);

        cdiDate_decode(vDateTime.date, &year, &month, &day);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID);
          auto const &var = varList.vars[varID];

          auto code = var.code;
          auto gridID = var.gridID;
          auto dig = (var.dataType == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;
          gridsize = var.nwpv * var.gridsize;
          auto nlon = gridInqXsize(gridID);
          auto nlat = gridInqYsize(gridID);
          auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
          auto missval = var.missval;

          cdiParamToString(var.param, paramstr, sizeof(paramstr));

          if (nlon * nlat != gridsize)
          {
            nlon = gridsize;
            nlat = 1;
          }

          cdo_read_field(streamID, array.data(), &numMissVals);

          auto vdate = cdiDate_get(vDateTime.date);
          auto vtime = cdiTime_get(vDateTime.time);

          if (operatorID == OUTPUTSRV)
            std::fprintf(stdout, "%4d %8g %8ld %4d %8zu %8zu %d %d\n", code, level, (long) vdate, vtime, nlon, nlat, 0, 0);

          if (operatorID == OUTPUTEXT) std::fprintf(stdout, "%8ld %4d %8g %8zu\n", (long) vdate, code, level, gridsize);

          if (operatorID == OUTPUTINT) { outputint(gridsize, array); }
          else if (operatorID == OUTPUTF) { outputf(nelem, format, gridsize, array); }
          else if (operatorID == OUTPUTTS)
          {
            if (gridsize > 1) cdo_abort("operator works only with one gridpoint!");

            std::fprintf(stdout, "%s %s %.*g\n", vDateStr.c_str(), vTimeStr.c_str(), dig, array[0]);
          }
          else if (operatorID == OUTPUTFLD)
          {
            int hour, minute, second, ms;
            cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
            double xdate = vdate - (vdate / 100) * 100 + (hour * 3600 + minute * 60 + second) / 86400.0;
            for (size_t i = 0; i < gridsize; ++i)
              if (fp_is_not_equal(array[i], missval))
                std::fprintf(stdout, "%g\t%g\t%g\t%.*g\n", xdate, grid_center_lat[i], grid_center_lon[i], dig, array[i]);
          }
          else if (operatorID == OUTPUTTAB)
          {
            auto is2dGrid = (gridtype == GRID_CURVILINEAR || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION);
            for (size_t i = 0; i < gridsize; ++i)
            {
              auto yind = i;
              auto xind = i;
              if (is2dGrid)
              {
                yind /= xsize;
                xind -= yind * xsize;
              }
              auto lon = grid_center_lon[i];
              auto lat = grid_center_lat[i];

              auto xval = hasRegxyCoordinates ? grid_xvals[xind] : grid_xvals[i];
              auto yval = hasRegxyCoordinates ? grid_yvals[yind] : grid_yvals[i];

              for (auto ki : keyIndices)
              {
                auto len = keyMap[ki].len;
                // clang-format off
                switch (keyMap[ki].idx)
                {
                  case kvalue:    std::fprintf(stdout, "%*.*g ", len, dig, array[i]); break;
                  case kx:        std::fprintf(stdout, "%*.*g ", len, dig, xval); break;
                  case ky:        std::fprintf(stdout, "%*.*g ", len, dig, yval); break;
                  case klon:      std::fprintf(stdout, "%*.*g ", len, dig, lon); break;
                  case klat:      std::fprintf(stdout, "%*.*g ", len, dig, lat); break;
                  case klev:      std::fprintf(stdout, "%*.*g ", len, dig, level); break;
                  case kbin:      std::fprintf(stdout, "%*.*g ", len, dig, level); break;
                  case kparam:    std::fprintf(stdout, "%*s ", len, paramstr); break;
                  case kcode:     std::fprintf(stdout, "%*d ", len, code); break;
                  case kname:     std::fprintf(stdout, "%*s ", len, var.name.c_str()); break;
                  case kxind:     std::fprintf(stdout, "%*zu ", len, xind + 1); break;
                  case kyind:     std::fprintf(stdout, "%*zu ", len, yind + 1); break;
                  case ktimestep: std::fprintf(stdout, "%*d ", len, tsID + 1); break;
                  case kdate:     std::fprintf(stdout, "%*s ", len, vDateStr.c_str()); break;
                  case ktime:     std::fprintf(stdout, "%*s ", len, vTimeStr.c_str()); break;
                  case kyear:     std::fprintf(stdout, "%*d ", len, year); break;
                  case kmonth:    std::fprintf(stdout, "%*d ", len, month); break;
                  case kday:      std::fprintf(stdout, "%*d ", len, day); break;
                }
                // clang-format on
              }
              std::fprintf(stdout, "\n");
            }
          }
          else if (operatorID == OUTPUTXYZ)
          {
            if (tsID == 0 && fieldID == 0) outputxyz(gridsize, array, missval, nlon, nlat, grid_center_lon, grid_center_lat);
          }
          else if (operatorID == OUTPUTARR) { outputarr(dig, gridsize, array); }
          else
          {
            if (gridInqType(gridID) == GRID_SPECTRAL && gridsize <= 156)
              outputsp(gridsize, array, gridInqTrunc(gridID));
            else
              output(gridsize, array);
          }
        }

        tsID++;
      }

      cdo_stream_close(streamID);
    }
  }

  void
  close() override
  {
  }
};
