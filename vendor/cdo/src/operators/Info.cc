/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Info       info            Dataset information
      Info       map             Dataset information and simple map
*/

#include <cdi.h>
#include <numbers>

#include "workerthread.h"
#include "cdo_options.h"
#include "cdo_math.h"
#include "process_int.h"
#include "mpmo_color.h"
#include "varray.h"
#include "printinfo.h"
#include "field_functions.h"
#include "cdo_zaxis.h"

namespace
{
struct InfoStat
{
  double min{ std::numeric_limits<double>::max() };
  double max{ -std::numeric_limits<double>::max() };
  double sum{ 0.0 };
  double sumi{ 0.0 };
  size_t numVals{ 0 };
  size_t numMissVals{ 0 };
  int numLevels{ 0 };
};
}  // namespace

static void
field_min_max_sum(Field const &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  auto func = [&](auto const &v, auto n) { return varray_min_max_sum(v, n, mms); };
  mms = field_operation(func, field, field.size);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
}

static size_t
field_min_max_sum_mv(Field const &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  auto func = [&](auto const &v, auto n, double mv) { return varray_min_max_sum_mv(v, n, mms, mv); };
  mms = field_operation(func, field, field.size, field.missval);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
  return mms.n;
}

static void
print_grid_index(int nlon, int nlat, int i)
{
  int index = (nlat < 10) ? 2 : (nlat < 100) ? 3 : (nlat < i) ? 4 : 5;

  std::stringstream s;
  s << std::string(index, ' ');
  for (int ilon = 0; ilon < nlon; ilon++) s << ((ilon + 1) / i) % 10;

  printf("%s\n", s.str().c_str());
}

static void
compute_level(double min, double max, double (&level)[10])
{
  auto step = (max - min) / 10.0;

  if (is_not_equal(step, 0.0))
  {
    auto a = std::pow(10.0, std::floor(std::log(step) / std::numbers::ln10));
    auto b = step / a;

    // clang-format off
    if      (b > 5) b = 0.5 * std::ceil(b / 0.5);
    else if (b > 2) b = 0.2 * std::ceil(b / 0.2);
    else if (b > 1) b = 0.1 * std::ceil(b / 0.1);
    else            b = 1;
    // clang-format on

    step = b * a;

    if (min < 0.0 && max > 0.0)
    {
      int min_n = (int) std::floor(10.0 * (-min) / (max - min) - 0.5);
      int max_n = (int) std::ceil(10.0 * (-min) / (max - min) - 0.5);
      min_n = std::max(min_n, 0);
      max_n = std::max(max_n, 1);
      level[min_n] = 0;
      for (int i = min_n - 1; i >= 0; i--) level[i] = level[i + 1] - step;
      for (int i = max_n; i < 9; ++i) level[i] = level[i - 1] + step;
    }
    else
    {
      level[0] = step * std::ceil(min / step + 0.5);
      for (int i = 1; i < 9; ++i) level[i] = level[i - 1] + step;
    }
  }
  else
    for (int i = 0; i < 9; ++i) level[i] = min;
}

static unsigned char
val_to_char(double x, double missval, double min, double max, const double (&level)[10])
{
  unsigned char c;
  if (fp_is_equal(x, missval)) { c = '.'; }
  else if (fp_is_equal(x, min) && fp_is_not_equal(min, max)) { c = 'm'; }
  else if (fp_is_equal(x, max) && fp_is_not_equal(min, max)) { c = 'M'; }
  else if (fp_is_equal(x, 0.0)) { c = '*'; }
  else if (x < 0)
  {
    c = '9';
    for (int i = 0; i < 9; ++i)
      if (level[i] > x)
      {
        c = i + '0';
        break;
      }
  }
  else
  {
    c = '0';
    for (int i = 8; i >= 0; i--)
      if (level[i] < x)
      {
        c = i + 1 + '0';
        break;
      }
  }

  return c;
}

static std::pair<TextMode, TextColor>
char_to_mode_and_color(unsigned char c, int &blinkMin, int &blinkMax)
{
  TextMode mode(MODELESS);
  TextColor color(BLACK);
  switch (c)
  {
      // clang-format off
    case '0': mode = BRIGHT  ; color = BLUE   ; break;
    case '1': mode = MODELESS; color = BLUE   ; break;
    case '2': mode = BRIGHT  ; color = CYAN   ; break;
    case '3': mode = MODELESS; color = CYAN   ; break;
    case '4': mode = MODELESS; color = GREEN  ; break;
    case '5': mode = MODELESS; color = YELLOW ; break;
    case '6': mode = MODELESS; color = RED    ; break;
    case '7': mode = BRIGHT  ; color = RED    ; break;
    case '8': mode = MODELESS; color = MAGENTA; break;
    case '9': mode = BRIGHT  ; color = MAGENTA; break;
    // clang-format on
    case 'm':
      (blinkMax == 1) ? mode = BLINK : mode = MODELESS, color = BLACK;
      if (blinkMax) blinkMax = 0;
      break;
    case 'M':
      (blinkMin == 1) ? mode = BLINK : mode = MODELESS, color = BLACK;
      if (blinkMin) blinkMin = 0;
      break;
  }

  return std::make_pair(mode, color);
}

template <typename T>
static void
print_map(int nlon, int nlat, Varray<T> const &varray, double mv, double min, double max)
{
  auto missval = static_cast<T>(mv);
  // source code from PINGO
  double level[10] = {};
  compute_level(min, max, level);

  printf("\n");

  for (int i = 1; i <= 4; ++i)
  {
    int current = 10000 / std::pow(10, i);
    if (nlon >= current) print_grid_index(nlon, nlat, current);
  }
  printf("\n");

  int blinkMin = 1, blinkMax = 1;
  for (int ilat = 0; ilat < nlat; ilat++)
  {
    printf("%0*d ", (nlat < 10) ? 1 : (nlat < 100) ? 2 : (nlat < 1000) ? 3 : 4, ilat + 1);
    for (int ilon = 0; ilon < nlon; ilon++)
    {
      double x = varray[ilat * nlon + ilon];
      auto c = val_to_char(x, missval, min, max, level);
      auto [mode, color] = char_to_mode_and_color(c, blinkMin, blinkMax);
      set_text_color(stdout, mode, color);
      putchar(c);
      reset_text_color(stdout);
    }
    printf(" %0*d\n", (nlat < 10) ? 1 : (nlat < 100) ? 2 : (nlat < 1000) ? 3 : 4, ilat + 1);
  }
  printf("\n");

  for (int i = 1; i <= 4; ++i)
  {
    int current = 10000 / std::pow(10, i);
    if (nlon >= current) print_grid_index(nlon, nlat, current);
  }
  printf("\n");

  for (int i = 0; i < 10; ++i)
  {
    printf("%d=%c%+9.3e,%+9.3e%c%s", i, '[', (i == 0) ? min : level[i - 1], (i == 9) ? max : level[i], ']',
           (i != 2 && i != 5 && i != 8) ? "  " : "");

    if (i == 2 || i == 5 || i == 8) printf("\n");
  }

  printf("*=0  .=miss  m=min=%+9.3e  M=max=%+9.3e\n", min, max);
  printf("\n");
}

static void
print_map(int nlon, int nlat, Field const &field, double min, double max)
{
  auto func = [&](auto const &v, double mv) { print_map(nlon, nlat, v, mv, min, max); };
  field_operation(func, field, field.missval);
}

template <typename T>
static size_t
complex_sum(Varray<T> const &v, double mv, size_t gridsize, double &sumr, double &sumi)
{
  T missval = static_cast<T>(mv);
  size_t n = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (fp_is_not_equal(v[i * 2], missval) && fp_is_not_equal(v[i * 2 + 1], missval))
    {
      sumr += v[i * 2];
      sumi += v[i * 2 + 1];
      n++;
    }
  }

  return n;
}

static size_t
field_complex_sum(Field const &field, double &sumr, double &sumi)
{
  auto func = [&](auto const &v, double mv, size_t gridsize) { return complex_sum(v, mv, gridsize, sumr, sumi); };
  return field_operation(func, field, field.missval, field.gridsize);
}

static void
infostat_init(InfoStat &infoStat)
{
  infoStat.numVals = 0;
  infoStat.numMissVals = 0;
  infoStat.numLevels = 0;
  infoStat.min = std::numeric_limits<double>::max();
  infoStat.max = -std::numeric_limits<double>::max();
  infoStat.sum = 0.0;
  infoStat.sumi = 0.0;
}

static void
print_header(int fileIndex, bool lvinfo, int operfunc)
{
  auto e = (operfunc == Func_Name) ? "Parameter name" : ((operfunc == Func_Code) ? "Code number" : "Parameter ID");

  set_text_color(stdout, BRIGHT);
  if (fileIndex)
    std::fprintf(stdout, "%6d :       Date     Time   %s Gridsize    Miss :     Minimum        Mean     Maximum : %s\n", fileIndex,
                 lvinfo ? "Nlevs" : "Level", e);
  else
    std::fprintf(stdout, "       :       Date     Time   %s Gridsize    Miss :     Minimum        Mean     Maximum : %s\n",
                 lvinfo ? "Nlevs" : "Level", e);
  reset_text_color(stdout);
}

static void
print_xheader(int fileIndex)
{
  auto e = "Parameter name";

  set_text_color(stdout, BRIGHT);
  if (fileIndex)
    std::fprintf(stdout, "%6d : NumSteps NumLevels  Gridsize   NumMiss :     Minimum        Mean     Maximum : %s\n", fileIndex, e);
  reset_text_color(stdout);
}

static void
compute_stat_real(Field const &field, InfoStat &infoStat, size_t &imiss, size_t gridsize)
{
  if (infoStat.numMissVals)
  {
    auto numVals = field_min_max_sum_mv(field, infoStat.min, infoStat.max, infoStat.sum);
    imiss = gridsize - numVals;
    infoStat.numVals += numVals;
  }
  else if (gridsize == 1)
  {
    infoStat.sum = (infoStat.numVals == 0) ? field[0] : infoStat.sum + field[0];
    infoStat.min = (infoStat.numVals == 0) ? field[0] : std::min(infoStat.min, field[0]);
    infoStat.max = (infoStat.numVals == 0) ? field[0] : std::max(infoStat.max, field[0]);
    infoStat.numVals += 1;
  }
  else
  {
    field_min_max_sum(field, infoStat.min, infoStat.max, infoStat.sum);
    infoStat.numVals += gridsize;
  }
}

static void
compute_stat_comp(Field const &field, InfoStat &infoStat, size_t &imiss, size_t gridsize)
{
  auto numVals = field_complex_sum(field, infoStat.sum, infoStat.sumi);
  imiss = gridsize - numVals;
  infoStat.numVals += numVals;
}

static void
print_stat_real(const InfoStat &infoStat)
{
  if (infoStat.numVals == 0)
    std::fprintf(stdout, "                     nan            ");
  else if (infoStat.numVals == 1)
    std::fprintf(stdout, "            %#12.5g            ", infoStat.sum);
  else
    std::fprintf(stdout, "%#12.5g%#12.5g%#12.5g", infoStat.min, infoStat.sum / infoStat.numVals, infoStat.max);
}

static void
print_stat_comp(const InfoStat &infoStat)
{
  auto arrmean_r = (infoStat.numVals > 0) ? infoStat.sum / infoStat.numVals : 0.0;
  auto arrmean_i = (infoStat.numVals > 0) ? infoStat.sumi / infoStat.numVals : 0.0;
  std::fprintf(stdout, "   -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
}

static void
print_xinfo(int numSteps, CdoVar const &var, InfoStat const &infoStat)
{
  std::fprintf(stdout, "%6d : ", var.ID + 1);
  set_text_color(stdout, GREEN);
  std::fprintf(stdout, "%8d   %7d %9zu %9zu ", numSteps, var.nlevels, var.gridsize, infoStat.numMissVals);
  reset_text_color(stdout);

  std::fprintf(stdout, ":");

  set_text_color(stdout, BLUE);
  (var.nwpv == CDI_REAL) ? print_stat_real(infoStat) : print_stat_comp(infoStat);
  reset_text_color(stdout);

  std::fprintf(stdout, " : ");

  // set_text_color(stdout, GREEN);
  std::fprintf(stdout, "%-14s", var.name.c_str());
  // reset_text_color(stdout);

  std::fprintf(stdout, "\n");
}

static void
print_info(int setNum, int levelID, CdiDateTime const &vDateTime, CdoVar const &var, int operfunc, bool lvinfo,
           InfoStat const &infoStat)
{
  char paramstr[32];
  cdiParamToString(var.param, paramstr, sizeof(paramstr));

  std::fprintf(stdout, "%6d :", setNum);

  auto vdateString = date_to_string(vDateTime.date);
  auto vtimeString = time_to_string(vDateTime.time);

  set_text_color(stdout, MAGENTA);
  std::fprintf(stdout, "%s %s ", vdateString.c_str(), vtimeString.c_str());
  reset_text_color(stdout);

  set_text_color(stdout, GREEN);
  if (lvinfo)
    std::fprintf(stdout, "%7d ", var.nlevels);
  else
    std::fprintf(stdout, "%7g ", cdo_zaxis_inq_level(var.zaxisID, levelID));

  std::fprintf(stdout, "%8zu %7zu ", var.gridsize, infoStat.numMissVals);
  reset_text_color(stdout);

  std::fprintf(stdout, ":");

  set_text_color(stdout, BLUE);
  (var.nwpv == CDI_REAL) ? print_stat_real(infoStat) : print_stat_comp(infoStat);
  reset_text_color(stdout);

  std::fprintf(stdout, " : ");

  // set_text_color(stdout, GREEN);
  // clang-format off
  if      (operfunc == Func_Name) std::fprintf(stdout, "%-14s", var.name.c_str());
  else if (operfunc == Func_Code) std::fprintf(stdout, "%4d   ", var.code);
  else                            std::fprintf(stdout, "%-14s", paramstr);
  // clang-format on
  // reset_text_color(stdout);

  std::fprintf(stdout, "\n");
}

static void
info(Field &field, int setNum, int streamIndex, int levelID, CdiDateTime vDateTime, CdoVar &var, int operfunc, bool printMap,
     bool lvinfo, bool lcinfo, InfoStat &infoStat)
{
  if (printMap) print_header(-(streamIndex + 1), lvinfo, operfunc);

  auto numMissVals = field.numMissVals;
  auto loutput = (not lvinfo and not lcinfo);

  if (loutput) infostat_init(infoStat);

  infoStat.numMissVals += numMissVals;
  infoStat.numLevels += 1;
  if (not lcinfo and (var.nlevels == infoStat.numLevels)) loutput = true;

  size_t numNANs = (Options::fast || std::isnan(field.missval)) ? 0 : field_num_NANs(field);
  var.counter += numNANs;
  if (numNANs && field.numMissVals == 0)
  {
    field.missval = cdo::NaN();
    infoStat.numMissVals += numNANs;
  }

  size_t imiss = 0;
  (var.nwpv == CDI_REAL) ? compute_stat_real(field, infoStat, imiss, var.gridsize)
                         : compute_stat_comp(field, infoStat, imiss, var.gridsize);

  if (loutput) print_info(setNum, levelID, vDateTime, var, operfunc, lvinfo, infoStat);

  if (imiss != numMissVals && numMissVals) cdo_warning("Found %zu of %zu missing values (%s)!", imiss, numMissVals, var.name);

  if (printMap)
  {
    auto gridID = var.gridID;
    auto gridtype = var.gridType;
    auto nlon = gridInqXsize(gridID);
    auto nlat = gridInqYsize(gridID);

    if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR
        || (gridtype == GRID_GENERIC && nlon * nlat == var.gridsize && nlon < 2048))
    {
      print_map(nlon, nlat, field, infoStat.min, infoStat.max);
    }
  }
}

class Info : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Info",
    // clang-format off
    .operators = { { "info", Func_Param, 0, InfoHelp },
                   { "infop", Func_Param, 0, InfoHelp },
                   { "infon", Func_Name, 0, InfoHelp },
                   { "infoc", Func_Code, 0, InfoHelp },
                   { "vinfo", Func_Name, 0, InfoHelp },
                   { "cinfo", Func_Name, 0, InfoHelp },
                   { "map", Func_Param, 0, InfoHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { -1, 0, NoRestriction },
  };
  inline static RegisterEntry<Info> registration = RegisterEntry<Info>();

  int operfunc{};

  bool printMap{};
  bool lvinfo{};
  bool lcinfo{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    auto VINFO = module.get_id("vinfo");
    auto CINFO = module.get_id("cinfo");
    auto MAP = module.get_id("map");

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    printMap = (operatorID == MAP);
    lvinfo = (operatorID == VINFO);
    lcinfo = (operatorID == CINFO);
  }

  void
  run() override
  {
    int numSets = 0;
    auto numStreams = cdo_stream_cnt();

    for (int streamIndex = 0; streamIndex < numStreams; streamIndex++)
    {
      auto streamID = cdo_open_read(streamIndex);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      auto taxisID = vlistInqTaxis(vlistID);

      VarList varList(vlistID);
      auto numVars = varList.numVars();
      if (numVars == 0) continue;

      auto runAsync = (Options::CDO_Async_Read > 0);
      auto workerThread = runAsync ? std::make_unique<WorkerThread>() : nullptr;
      auto numTasks = runAsync ? 2 : 1;

      Field fieldVector[2];
      std::vector<InfoStat> infoStatList(numVars);

      if (lcinfo)
        for (auto &infoStat : infoStatList) infostat_init(infoStat);

      if (lcinfo) { print_xheader(-(streamIndex + 1)); }
      else if (not printMap) { print_header(-(streamIndex + 1), lvinfo, operfunc); }

      numSets = 0;
      int tsID = 0;
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields == 0) break;

        auto vDateTime = taxisInqVdatetime(taxisID);

        if (lvinfo)
        {
          if (numFields == 1 && runAsync && numSets > 0) { workerThread->wait(); }
          for (auto &infoStat : infoStatList) infostat_init(infoStat);
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID);
          auto &var = varList.vars[varID];
          auto &field = fieldVector[numSets % numTasks];
          field.init(var);
          cdo_read_field(streamID, field);

          if (runAsync && numSets > 0) { workerThread->wait(); }

          numSets = lvinfo ? varID + 1 : numSets + 1;

          std::function<void()> info_task
              = std::bind(info, std::ref(field), numSets, streamIndex, levelID, vDateTime, std::ref(var), operfunc, printMap,
                          lvinfo, lcinfo, std::ref(infoStatList[varID]));

          runAsync ? workerThread->doAsync(info_task) : info_task();
        }

        tsID++;
      }

      if (runAsync) { workerThread->wait(); }

      cdo_stream_close(streamID);

      if (lcinfo)
        for (auto const &var : varList.vars) print_xinfo(tsID, var, infoStatList[var.ID]);

      for (auto const &var : varList.vars)
      {
        if (var.counter > 0)
        {
          cdo_warning("%s contains %zu NaNs which are not treated as missing values. "
                      "This can lead to incorrect CDO results in all other arithmetic functions!",
                      var.name, var.counter);
        }
      }
    }

    if (numSets > 36 && !printMap && !lcinfo) print_header(0, lvinfo, operfunc);
  }

  void
  close() override
  {
  }
};
