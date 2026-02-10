/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
/* ============================================================= */
/*                                                               */
/* postprocessing program for ECHAM data and ECMWF analysis data */
/*                                                               */
/* Luis     Kornblueh   - MPI    Hamburg                         */
/* Uwe      Schulzweida - MPI    Hamburg                         */
/* Arno     Hellbach    - DKRZ   Hamburg                         */
/* Edilbert Kirk        - MI Uni Hamburg                         */
/* Michael  Ponater     - DLR    Oberpfaffenhofen                */
/*                                                               */
/* ============================================================= */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cdi.h>

#include "cdo_default_values.h"

#define streamOpenWrite cdo_open_write
#define streamDefVlist cdo_def_vlist
#define streamDefTimestep cdo_def_timestep

#include "afterburner.h"
#include "compare.h"
#include "cdo_options.h"
#include "cdi_lockedIO.h"
#include "gaussian_latitudes.h"

int scan_par_obsolete(char *namelist, const char *name, int def);
int scan_par(int verbose, char *namelist, const char *name, int def);
void scan_code(char *namelist, struct Variable *vars, int maxCodes, int *numCodes);
int scan_time(int verbose, char *namelist, int *hours, int max_hours);
void scan_darray(char *namelist, const char *name, double *values, int maxValues, int *numValues);

char zaxistypename[CDI_MAX_NAME];

struct ReadArgs
{
  int lana;
  int numFields;
  int numFieldsNext;
  struct Variable *vars;
  AfterControl *globs;
};

static bool lstdout = true;
static int Source = 0;
static int ofiletype = -1;
static int DataType = -1;

static char *filename;
static std::vector<const char *> ifiles;
static char *ifile = nullptr;

static int ofileidx = 0;

static int specGridID = -1;
static int gaussGridID = -1;
static int iVertID = -1;
static int oVertID = -1;
static int oVertID_half = -1;

static bool Lhybrid2pressure = false;

static int TsID;

#define TIMESTEP_INTERVAL -1
#define MONTHLY_INTERVAL 0
#define DAILY_INTERVAL 1
#define UNLIM_INTERVAL 2

#define MaxHours 24
static int nrqh;
static int hours[MaxHours + 1];

static std::vector<double> LevelFound;

static inline bool
filetype_is_netcdf(int filetype)
{
  return filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4
         || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR;
}

static void
lprintf(FILE *fp)
{
  int num = 67;
  int cval = '-';

  fprintf(fp, " ");
  for (int inum = 0; inum < num; inum++) fprintf(fp, "%c", cval);
  fprintf(fp, "\n");
}

static void
FreeMean(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].mean)
    {
      std::free(vars[code].mean);
      vars[code].mean = nullptr;
    }
}

static void
after_PostProcess(AfterControl &globs)
{
  if (globs.EndOfInterval)
  {
    if (lstdout)
    {
      if (globs.OutputInterval == DAILY_INTERVAL)
        fprintf(stdout, " Processed Day %2d  Month %2d  Year %04d", globs.OldDate.dy, globs.OldDate.mo, globs.OldDate.yr);
      else if (globs.OutputInterval == MONTHLY_INTERVAL)
        fprintf(stdout, " Processed Month %2d  Year %04d", globs.OldDate.mo, globs.OldDate.yr);
      else if (globs.OutputInterval == UNLIM_INTERVAL)
        fprintf(stdout, " Processed range from %6.4d-%2.2d-%2.2d to %6.4d-%2.2d-%2.2d", globs.StartDate.yr, globs.StartDate.mo,
                globs.StartDate.dy, globs.OldDate.yr, globs.OldDate.mo, globs.OldDate.dy);

      if (globs.Mean)
        fprintf(stdout, "  (Mean of %3d Terms)\n", globs.MeanCount);
      else
        fprintf(stdout, "   Terms %3d\n", globs.MeanCount);
    }

    globs.EndOfInterval = false;
    globs.MeanCount = 0;
  }
}

/* ================= */
/* switch input file */
/* ================= */
static void
after_SwitchFile(AfterControl *globs)
{
  bool echam4 = false;
  char y3, y2, y1, y0;
  char m1, m0;
  char d1, d0;

  streamClose(globs->istreamID);

  if (globs->Multi > 0)
  {
    int i = std::strlen(ifile);
    if (i < 10)
    {
      fprintf(stderr, " Not a valid filename: %s \n", ifile);
      exit(1);
    }

    if (ifile[i - 3] == '.')
    {
      echam4 = true;
      y3 = ifile[i - 9];
      y2 = ifile[i - 8];
      y1 = ifile[i - 7];
      y0 = ifile[i - 6];
      m1 = ifile[i - 5];
      m0 = ifile[i - 4];
      d1 = ifile[i - 2];
      d0 = ifile[i - 1];
    }
    else
    {
      y3 = ifile[i - 6];
      y2 = ifile[i - 5];
      y1 = ifile[i - 4];
      y0 = ifile[i - 3];
      m1 = ifile[i - 2];
      m0 = ifile[i - 1];
      d1 = '0';
      d0 = '1';
    }

    for (int n = 0; n < globs->DayIn; ++n)
    {
      if (d0 == '9')
      {
        d0 = '0';
        d1++;
      }
      else
        d0++;
      if (d1 == '3' && d0 > '0')
      {
        d1 = '0';
        d0 = '1';
        if (m1 == '0')
        {
          if (m0 == '9')
          {
            m0 = '0';
            m1 = '1';
          }
          else
            m0++;
        }
        else
        {
          if (m0 < '2')
            m0++;
          else
          {
            m1 = '0';
            m0 = '1';
            y0++;
            if (y0 > '9')
            {
              y0 = '0';
              y1++;
            }
            if (y1 > '9')
            {
              y1 = (char) '0';
              if (std::isdigit((int) y2))
                y2++;
              else
                y2 = '1';
              if (y2 > '9')
              {
                y2 = (char) '0';
                if (std::isdigit((int) y3))
                  y3++;
                else
                  y3 = '1';
              }
            }
          }
        }
      }
    }

    if (echam4)
    {
      ifile[i - 9] = y3;
      ifile[i - 8] = y2;
      ifile[i - 7] = y1;
      ifile[i - 6] = y0;
      ifile[i - 5] = m1;
      ifile[i - 4] = m0;
      ifile[i - 2] = d1;
      ifile[i - 1] = d0;
    }
    else
    {
      ifile[i - 6] = y3;
      ifile[i - 5] = y2;
      ifile[i - 4] = y1;
      ifile[i - 3] = y0;
      ifile[i - 2] = m1;
      ifile[i - 1] = m0;
    }

    globs->Multi--;
  }

  if (globs->Nfiles > 0) ifile = (char *) ifiles[--globs->Nfiles];

  fprintf(stderr, " Continuation file: %s\n", ifile);

  globs->istreamID = stream_open_read_locked(ifile);

  globs->ivlistID = streamInqVlist(globs->istreamID);
  globs->taxisID = vlistInqTaxis(globs->ivlistID);
}

static CdiDateTime
after_getDateTime(const struct Date &datetime)
{
  CdiDateTime cdiDateTime{};
  cdiDateTime.date = cdiDate_encode(datetime.yr, datetime.mo, datetime.dy);
  cdiDateTime.time = cdiTime_encode(datetime.hr, datetime.mn, 0, 0);
  return cdiDateTime;
}

static void
after_setDateTime(struct Date *datetime, CdiDateTime const &cdiDateTime)
{
  int sec, ms;
  cdiDate_decode(cdiDateTime.date, &datetime->yr, &datetime->mo, &datetime->dy);
  cdiTime_decode(cdiDateTime.time, &datetime->hr, &datetime->mn, &sec, &ms);
}

static void
after_printProcessStatus(int tsID)
{
  static bool counthead = false;

  if (tsID == -1)
  {
    if (cdo::stdoutIsTerminal)
    {
      fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      fflush(stdout);
    }

    counthead = false;
  }
  else
  {
    if (!counthead)
    {
      if (cdo::stdoutIsTerminal) fprintf(stdout, " Process timestep :       ");

      counthead = true;
    }

    if (cdo::stdoutIsTerminal)
    {
      fprintf(stdout, "\b\b\b\b\b\b%6d", tsID);
      fflush(stdout);
    }
  }
}

static int
after_setNextDate(AfterControl *globs)
{
  int numFields = 0;

  bool righttime = false;
  while (true)
  {
    numFields = streamInqTimestep(globs->istreamID, TsID);
    if (numFields == 0 && (globs->Multi > 0 || globs->Nfiles > 0))
    {
      if (lstdout) after_printProcessStatus(-1);

      after_SwitchFile(globs);

      if (globs->istreamID >= 0)
      {
        TsID = 0;
        numFields = streamInqTimestep(globs->istreamID, TsID);
      }
    }
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(globs->taxisID);
    after_setDateTime(&globs->NextDate, vDateTime);

    for (int i = 0; i < nrqh; ++i)
      if (hours[i] < 0 || hours[i] == globs->NextDate.hr)
      {
        righttime = true;
        break;
      }

    if (righttime) break;

    TsID += 1;
  }

  return numFields;
}

static void
after_readTimestep(ReadArgs &readArgs)
{
  int varID, gridID, zaxisID, levelID, timeType;
  size_t numMissVals;

  auto numFields = readArgs.numFields;
  auto analysisData = readArgs.lana;
  auto vars = readArgs.vars;
  auto globs = readArgs.globs;

  for (int code = 0; code < MaxCodes; ++code) vars[code].numMissVals0 = 0;

  while (numFields--)
  {
    streamInqField(globs->istreamID, &varID, &levelID);

    auto code = vlistInqVarCode(globs->ivlistID, varID);
    if (code <= 0 || code >= MaxCodes) continue;

    // Skip fields containing unneeded codes
    if (!vars[code].needed0) continue;

    vlistInqVar(globs->ivlistID, varID, &gridID, &zaxisID, &timeType);

    auto leveltype = zaxisInqType(zaxisID);

    // Skip fields with unselected levels
    int levelOffset = -1;
    // if ( vars[code].ozaxisID != vars[code].izaxisID && ! Lhybrid2pressure )
    if ((vars[code].ozaxisID != vars[code].izaxisID) && (leveltype == ZAXIS_PRESSURE))
    {
      auto level = (int) zaxisInqLevel(zaxisID, levelID);
      for (int i = 0; i < globs->NumLevelRequest; ++i)
      {
        if (is_equal(globs->LevelRequest[i], level))
        {
          levelOffset = i;
          break;
        }
      }

      if (levelOffset < 0) continue;

      zaxisID = vars[code].ozaxisID;
      levelID = levelOffset;
    }

    if (globs->Debug)
    {
      fprintf(stderr, "T%d", globs->Truncation);

      fprintf(stderr, "  Code %3d   Level%6d   %6.4d-%2.2d-%2.2d  %2.2d:%2.2d:00\n", code, (int) zaxisInqLevel(zaxisID, levelID),
              globs->OldDate.yr, globs->OldDate.mo, globs->OldDate.dy, globs->OldDate.hr, globs->OldDate.mn);
    }

    if (analysisData)
    {
      streamReadField(globs->istreamID, globs->varray.data(), &numMissVals);
      after_AnalysisAddField(globs, vars, code, gridID, zaxisID, levelID, numMissVals);
    }
    else
    {
      double *dataptr = after_get_dataptr(vars, code, gridID, zaxisID, levelID);
      streamReadField(globs->istreamID, dataptr, &numMissVals);
      after_EchamAddField(globs, vars, code, gridID, zaxisID, levelID, numMissVals);
    }

    if (iVertID != -1 && oVertID != -1 && (vars[code].izaxisID == iVertID)) vars[code].ozaxisID = oVertID;
  }

  TsID++;

  // printf("%3d  date = %d  time = %04d\n", TsID, vdate, vtime);

  readArgs.numFieldsNext = after_setNextDate(globs);
}

static void
after_defineNextTimestep(const AfterControl &globs)
{
  static int otsID = 0;
  taxisDefVdatetime(globs.taxisID2, after_getDateTime(globs.OldDate));

  if (globs.Mean != 2)
  {
    if (otsID == 0)
    {
      auto nvars = vlistNvars(globs.ovlistID);
      if (nvars == 0) afterAbort("No variable selected!");
      vlistDefTaxis(globs.ovlistID, globs.taxisID2);
      streamDefVlist(globs.ostreamID, globs.ovlistID);
    }
    taxisDefNumavg(globs.taxisID2, globs.MeanCount + 1);
    streamDefTimestep(globs.ostreamID, otsID);
  }

  otsID++;
}

static void
after_setEndOfInterval(AfterControl &globs, int numFields)
{
  if (numFields == 0)
  {
    globs.EndOfInterval = true;
    return;
  }

  if (globs.OutputInterval == DAILY_INTERVAL)
    globs.EndOfInterval = (globs.NewDate.dy != globs.OldDate.dy);
  else if (globs.OutputInterval == MONTHLY_INTERVAL)
    globs.EndOfInterval = (globs.NewDate.mo != globs.OldDate.mo);
  else if (globs.OutputInterval == UNLIM_INTERVAL)
    globs.EndOfInterval = false;
  else
    afterAbort("output interval %d not implemented!", globs.OutputInterval);
}

static void
after_moveTimestep(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code) vars[code].numMissVals = vars[code].numMissVals0;

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].hybrid0)
    {
      vars[code].hybrid = vars[code].hybrid0;
      vars[code].hybrid0 = nullptr;
    }

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].spectral0)
    {
      vars[code].spectral = vars[code].spectral0;
      vars[code].spectral0 = nullptr;
    }

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].grid0)
    {
      vars[code].grid = vars[code].grid0;
      vars[code].grid0 = nullptr;
    }
}

static void
after_check_content(struct Variable *vars, int timestep)
{
  for (int code = 0; code < 272; ++code)
  {
    // if ( code == GEOPOTENTIAL ) continue;
    if (code == SLP) continue;
    if (code == GEOPOTHEIGHT) continue;
    if (code == STREAM) continue;
    if (code == VELOPOT) continue;
    if (code == U_WIND) continue;
    if (code == V_WIND) continue;
    if (code == OMEGA) continue;
    if (code == RHUMIDITY) continue;
    if (code == LOW_CLOUD) continue;
    if (code == MID_CLOUD) continue;
    if (code == HIH_CLOUD) continue;
    if (code == PS) continue;
    if (code == HUMIDITY)
    {
      if (vars[code].needed && !vars[code].selected && vars[code].spectral == nullptr && vars[code].hybrid == nullptr)
      {
        static auto printWarning = true;
        if (printWarning) cdo_warning("No humidity in data file, set to zero !");
        printWarning = false;
        vars[code].needed = false;
      }
    }
    else
    {
      if (vars[code].needed && !vars[code].comp && vars[code].spectral == nullptr && vars[code].hybrid == nullptr)
      {
        afterAbort("Code  %3d not found at timestep %d!", code, timestep);
      }
    }
  }
}

static void
after_control(AfterControl &globs, struct Variable *vars)
{
  int numFields = 0;
  ReadArgs readArgs;

  if (afterReadAsync) afterWorkerThread = new WorkerThread;

  for (int code = 0; code < MaxCodes; ++code) vars[code].needed0 = vars[code].needed;

  TsID = 0;

  bool righttime = false;
  while (true)
  {
    numFields = streamInqTimestep(globs.istreamID, TsID);
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(globs.taxisID);
    after_setDateTime(&globs.StartDate, vDateTime);
    after_setDateTime(&globs.NewDate, vDateTime);

    for (int i = 0; i < nrqh; ++i)
      if (hours[i] < 0 || hours[i] == globs.NewDate.hr)
      {
        righttime = true;
        break;
      }

    if (righttime) break;

    TsID++;
  }

  auto taxis_is_relative = (taxisInqType(globs.taxisID) == TAXIS_RELATIVE);
  auto rDateTime = taxis_is_relative ? taxisInqRdatetime(globs.taxisID) : after_getDateTime(globs.StartDate);

  if (filetype_is_netcdf(ofiletype))
  {
    taxisDefCalendar(globs.taxisID2, CALENDAR_PROLEPTIC);
    taxisDefType(globs.taxisID2, TAXIS_RELATIVE);
    taxisDefTunit(globs.taxisID2, TUNIT_DAY);
    taxisDefRdatetime(globs.taxisID2, rDateTime);
  }

  globs.OldDate = globs.NewDate;

  bool tsFirst = true;

  while (numFields > 0)
  {
    readArgs.numFields = numFields;
    readArgs.lana = globs.AnalysisData;
    readArgs.vars = vars;
    readArgs.globs = &globs;
    std::function<void()> readTimestep_func = std::bind(after_readTimestep, std::ref(readArgs));

    if (tsFirst || !afterReadAsync)
    {
      afterReadAsync ? afterWorkerThread->doAsync(readTimestep_func) : readTimestep_func();

      if (tsFirst && globs.Type > 0) after_legini_setup(globs, vars);

      if (afterReadAsync) { afterWorkerThread->wait(); }
      tsFirst = false;
    }
    else { afterWorkerThread->wait(); }

    numFields = readArgs.numFieldsNext;

    globs.MeanCount0 = globs.MeanCount;
    globs.NewDate = globs.NextDate;

    after_moveTimestep(vars);

    if (numFields && afterReadAsync) { afterWorkerThread->doAsync(readTimestep_func); }

    after_setEndOfInterval(globs, numFields);

    if (lstdout) after_printProcessStatus(TsID);

    if (lstdout && globs.EndOfInterval) after_printProcessStatus(-1);

    if (globs.Mean == 0 || globs.EndOfInterval)
    {
      if (!globs.AnalysisData) after_check_content(vars, globs.TermCount + 1);
      after_defineNextTimestep(globs);
    }

    globs.AnalysisData ? after_processPL(globs, vars) : after_processML(globs, vars);

    after_PostProcess(globs);

    if (numFields)
    {
      if (globs.AnalysisData)
        after_AnalysisDependencies(vars, MaxCodes);
      else
        after_EchamDependencies(vars, MaxCodes, globs.Type, Source);
    }

    globs.OldDate = globs.NewDate;
  }

  if (afterWorkerThread) delete afterWorkerThread;
}

static void
after_setLevel(AfterControl &globs)
{
  int found;
  bool removeLevel[MaxLevel];
  double level;
  bool checkLevel = true;
  // default pressure level
  long plevelDefault[]
      = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
  // default height level
  long hlevelDefault[] = { 0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000 };

  int numplevelDefault = sizeof(plevelDefault) / sizeof(plevelDefault[0]);
  int numhlevelDefault = sizeof(hlevelDefault) / sizeof(hlevelDefault[0]);

  int iLevelType = (iVertID != -1) ? zaxisInqType(iVertID) : -1;
  if (iLevelType == ZAXIS_HYBRID && globs.Type > 20) Lhybrid2pressure = true;

  if (globs.Verbose) lprintf(stdout);

  if (globs.NumLevelRequest == 0)
  {
    if (iVertID == -1)
    {
      if (globs.Verbose) fprintf(stdout, " No level detected\n");
    }
    else
    {
      if (Lhybrid2pressure)
      {
        if (globs.unitsel == 0)
        {
          if (globs.Verbose) fprintf(stdout, " Default pressure level selected:\n");
          globs.NumLevelRequest = numplevelDefault;
          for (int l = 0; l < globs.NumLevelRequest; ++l) globs.LevelRequest[l] = plevelDefault[l];
          oVertID = zaxisCreate(ZAXIS_PRESSURE, globs.NumLevelRequest);
          zaxisDefLevels(oVertID, globs.LevelRequest);
        }
        else
        {
          if (globs.Verbose) fprintf(stdout, " Default height level selected:\n");
          globs.NumLevelRequest = numhlevelDefault;
          for (int l = 0; l < globs.NumLevelRequest; ++l) globs.LevelRequest[l] = hlevelDefault[l];
          oVertID = zaxisCreate(ZAXIS_HEIGHT, globs.NumLevelRequest);
          zaxisDefLevels(oVertID, globs.LevelRequest);
        }
      }
      else
      {
        if (globs.Verbose)
        {
          if (iLevelType == ZAXIS_HYBRID)
            fprintf(stdout, " All detected hybrid level selected:\n");
          else
            fprintf(stdout, " All detected pressure level selected:\n");
        }
        globs.NumLevelRequest = globs.NumLevelFound;
        for (int l = 0; l < globs.NumLevelRequest; ++l) globs.LevelRequest[l] = LevelFound[l];
        int maxLev = globs.NumLevelFound;
        if (maxLev == (globs.numHalfLevels - 1)) globs.LevelRequest[maxLev] = globs.numHalfLevels;
        oVertID = iVertID;
      }
    }
    checkLevel = false;
  }
  else
  {
    if (iVertID == -1)
    {
      if (globs.Verbose) fprintf(stdout, " No level detected\n");
      checkLevel = false;
    }
    else if (globs.NumLevelRequest == 1 && is_equal(globs.LevelRequest[0], 0))
    {
      if (globs.Verbose) fprintf(stdout, " No level selected\n");
      globs.NumLevelRequest = 0;
      checkLevel = false;
    }
    else if (globs.Verbose)
    {
      if (Lhybrid2pressure) { fprintf(stdout, " Selected %s level:\n", (globs.unitsel == 0) ? "pressure" : "height"); }
      else
      {
        if (iLevelType == ZAXIS_HYBRID)
          fprintf(stdout, " Selected hybrid level:\n");
        else
          fprintf(stdout, " Selected %s level:\n", (globs.unitsel == 0) ? "pressure" : "height");
      }
    }
  }

  if (globs.Verbose && iVertID != -1)
    for (int l = 0; l < globs.NumLevelRequest; ++l) fprintf(stdout, "  Level %2d = %13.4f\n", l + 1, globs.LevelRequest[l]);

  if (checkLevel)
  {
    for (int k = 0; k < globs.NumLevelRequest; ++k) removeLevel[k] = false;
    for (int k = 0; k < globs.NumLevelRequest; ++k)
    {
      level = globs.LevelRequest[k];
      for (int l = k + 1; l < globs.NumLevelRequest; ++l)
        if (removeLevel[l] == false && is_equal(level, globs.LevelRequest[l]))
        {
          if (globs.Verbose) fprintf(stdout, "  Level %2d = %13.4f double request\n", l + 1, globs.LevelRequest[l]);
          removeLevel[l] = true;
        }
    }

    {
      int l = 0;
      for (int k = 0; k < globs.NumLevelRequest; ++k)
        if (removeLevel[k] == false) globs.LevelRequest[l++] = globs.LevelRequest[k];

      globs.NumLevelRequest = l;
    }

    if (globs.AnalysisData || globs.Type < 30)
    {
      for (int k = 0; k < globs.NumLevelRequest; ++k) removeLevel[k] = false;
      for (int k = 0; k < globs.NumLevelRequest; ++k)
      {
        level = globs.LevelRequest[k];
        found = false;
        for (int l = 0; l < globs.NumLevelFound; ++l)
          if (is_equal(level, LevelFound[l])) found = true;

        if (!found)
        {
          fprintf(stdout, "  Level %2d = %14.4f not in input\n", k + 1, globs.LevelRequest[k]);
          removeLevel[k] = true;
        }
      }

      {
        int l = 0;
        for (int k = 0; k < globs.NumLevelRequest; ++k)
          if (removeLevel[k] == false) globs.LevelRequest[l++] = globs.LevelRequest[k];

        if (l != globs.NumLevelRequest)
        {
          if (globs.Verbose) lprintf(stdout);
          afterAbort("Inconsistent or invalid level list!");
        }

        globs.NumLevelRequest = l;
      }
    }
  }

  if (globs.Verbose) lprintf(stdout);
}

static void
after_defineLevel(const AfterControl &globs, struct Variable *vars)
{
  // hybrid, pressure, height

  switch (globs.Type)
  {
    case 0:
    case 10:
    case 11:
    case 20:
    {
      if (iVertID == -1) break;

      int iVertType = zaxisInqType(iVertID);
      if (iVertType == ZAXIS_HYBRID || iVertType == ZAXIS_HYBRID_HALF)
      {
        if (oVertID == -1)
        {
          if (globs.NumLevelRequest > globs.NumLevelFound) afterAbort("Too many level requested");

          if (globs.NumLevelFound == globs.NumLevelRequest)
          {
            int i;
            for (i = 0; i < globs.NumLevelRequest; ++i)
              if (is_not_equal(globs.LevelRequest[i], LevelFound[i])) break;

            if (i == globs.NumLevelRequest) oVertID = iVertID;
          }

          if (oVertID == -1 && globs.NumLevelRequest > 0)
          {
            oVertID = zaxisCreate(ZAXIS_HYBRID, globs.NumLevelRequest);
            zaxisDefLevels(oVertID, globs.LevelRequest);
            zaxisDefVct(oVertID, globs.nvct, globs.vct.data());
          }
        }

        if (vars[GEOPOTHEIGHT].selected && oVertID_half == -1)
        {
          if (oVertID_half == -1 && globs.NumLevelRequest > 0)
          {
            oVertID_half = zaxisCreate(ZAXIS_HYBRID_HALF, globs.NumLevelRequest + 1);
            zaxisDefLevels(oVertID_half, globs.LevelRequest);
            zaxisDefVct(oVertID_half, globs.nvct, globs.vct.data());
          }
        }

        for (int code = 0; code < MaxCodes; ++code)
        {
          if (vars[code].selected)
          {
            int zaxisID = vars[code].izaxisID;
            if (zaxisID != -1)
            {
              int zaxisType = zaxisInqType(zaxisID);
              int numLevels = zaxisInqSize(zaxisID);
              if ((zaxisType == ZAXIS_HYBRID || zaxisType == ZAXIS_HYBRID_HALF) && numLevels >= globs.NumLevelRequest)
              {
                vars[code].ozaxisID = (code == GEOPOTHEIGHT) ? oVertID_half : oVertID;
              }
            }
          }
        }
      }
      else
      {
        zaxisName(zaxisInqType(iVertID), zaxistypename);
        afterAbort("%s level data unsupported for TYPE %d", zaxistypename, globs.Type);
      }
      break;
    }
    case 30:
    case 40:
    case 41:
    case 50:
    case 60:
    case 61:
    case 70:
    {
      if (iVertID == -1) break;

      if (oVertID == -1)
      {
        oVertID = zaxisCreate((globs.unitsel == 0) ? ZAXIS_PRESSURE : ZAXIS_HEIGHT, globs.NumLevelRequest);
        zaxisDefLevels(oVertID, globs.LevelRequest);
      }

      for (int code = 0; code < MaxCodes; ++code)
      {
        if (vars[code].selected)
        {
          int zaxisID = vars[code].izaxisID;
          if (zaxisID != -1)
          {
            int zaxisType = zaxisInqType(zaxisID);
            int numLevels = zaxisInqSize(zaxisID);
            if (zaxisType == zaxisInqType(iVertID) && (numLevels == globs.NumLevel || numLevels == globs.NumLevel + 1)
                && numLevels > 1)
              vars[code].ozaxisID = oVertID;
          }
        }
      }

      break;
    }
    default: afterAbort("TYPE %d unsupported", globs.Type);
  }
}

static void
after_defineGrid(const AfterControl &globs, struct Variable *vars)
{
  int ogridID = -1;

  // spectral, fourier, gauss, zonal mean

  switch (globs.Type)
  {
    case 0:
    case 50:
    {
      if (specGridID == -1)
      {
        if (globs.DimSP == 0) afterAbort("dim spectral undefined");
        if (globs.Truncation == 0) afterAbort("truncation undefined");

        specGridID = gridCreate(GRID_SPECTRAL, globs.DimSP);
        gridDefTrunc(specGridID, globs.Truncation);
      }

      ogridID = specGridID;
      break;
    }
    case 20:
    case 30:
    case 70:
    {
      if (gaussGridID == -1)
      {
        if (globs.Longitudes == 0) afterAbort("number of longitudes undefined");
        if (globs.Latitudes == 0) afterAbort("number of latitudes undefined");

        gaussGridID = gridCreate(GRID_GAUSSIAN, globs.Longitudes * globs.Latitudes);
        gridDefXsize(gaussGridID, globs.Longitudes);
        gridDefYsize(gaussGridID, globs.Latitudes);

        size_t nlon = globs.Longitudes;
        std::vector<double> lons(nlon);
        for (size_t i = 0; i < nlon; ++i) lons[i] = i * (360.0 / nlon);
        gridDefXvals(gaussGridID, lons.data());

        size_t nlat = globs.Latitudes;
        std::vector<double> lats(nlat), latw(nlat);
        gaussian_latitudes(nlat, lats.data(), latw.data());
        for (size_t j = 0; j < nlat; ++j) lats[j] = 180. / M_PI * std::asin(lats[j]);
        gridDefYvals(gaussGridID, lats.data());
      }

      ogridID = gaussGridID;
      break;
    }
    case 10:
    case 40:
    case 60:
    {
      if (globs.Fouriers == 0) afterAbort("number of fourier coefficients undefined");
      if (globs.Latitudes == 0) afterAbort("number of latitudes undefined");

      ogridID = gridCreate(GRID_FOURIER, globs.Fouriers * globs.Latitudes);
      gridDefXsize(ogridID, globs.Latitudes);
      gridDefYsize(ogridID, globs.Fouriers);
      break;
    }
    case 11:
    case 41:
    case 61:
    {
      if (globs.Latitudes == 0) afterAbort("Number of latitudes undefined");

      ogridID = gridCreate(GRID_GAUSSIAN, globs.Latitudes);
      gridDefXsize(ogridID, 1);
      gridDefYsize(ogridID, globs.Latitudes);
      break;
    }
    default: afterAbort("TYPE %d unsupported", globs.Type);
  }

  if (ogridID != -1)
    for (int code = 0; code < MaxCodes; ++code)
    {
      if (vars[code].selected) vars[code].ogridID = ogridID;
    }

  if (ogridID == -1) afterAbort("out grid undefined");
}

static void
after_setCodes(const AfterControl &globs, struct Variable *vars, int maxCodes, int numCodes)
{
  if (globs.Verbose) lprintf(stdout);

  if (numCodes == 0)
  {
    if (globs.Verbose) fprintf(stdout, " All detected codes selected:\n");

    for (int code = 0; code < maxCodes; ++code)
      if (vars[code].detected) vars[code].selected = 1;
  }
  else if (globs.Verbose)
    fprintf(stdout, " Selected codes:\n");

  if (globs.Verbose)
  {
    fprintf(stdout, "  Table Code Name              Longname\n");
    fprintf(stdout, "  ----- ---- ----              --------\n");
  }

  for (int code = 0; code < maxCodes; ++code)
    if (vars[code].selected)
    {
      char name[CDI_MAX_NAME];
      name[0] = 0;
      char longname[CDI_MAX_NAME];
      longname[0] = 0;
      int tableID;
      int table = 0;
      auto varID = vars[code].ivarID;

      if (varID == CDI_UNDEFID)
      {
        auto modelID = vlistInqVarModel(globs.ivlistID, 0);
        table = 128;
        tableID = tableInq(modelID, table, nullptr);

        vars[code].tableID = tableID;
      }
      else
      {
        tableID = vlistInqVarTable(globs.ivlistID, varID);
        table = tableInqNum(tableID);
        vlistInqVarName(globs.ivlistID, varID, name);
        vlistInqVarLongname(globs.ivlistID, varID, longname);
      }

      if (!name[0]) tableInqEntry(tableID, code, -1, name, longname, nullptr);

      if (globs.Verbose)
      {
        fprintf(stdout, " %5d", table);
        fprintf(stdout, " %4d", code);
        if (!name[0])
          fprintf(stdout, "  var%d", code);
        else
        {
          fprintf(stdout, "  %-16s", name);
          if (longname[0]) fprintf(stdout, "  %s", longname);
        }
        fprintf(stdout, "\n");
      }
    }
}

static void
after_checkNamelist(const AfterControl &globs)
{
  if (globs.Mean && globs.Type < 20) afterAbort("Mean is only available for TYPE >= 20!");

  if (globs.extrapolate == false && globs.Type >= 30)
  {
    if (globs.Type > 30) afterAbort("EXTRAPOLATE = 0 is only available for TYPE = 30!");
    if (globs.Mean) afterAbort("EXTRAPOLATE = 0 is only available with MEAN = 0!");
  }
}

static void
after_parini(AfterControl &globs, struct Variable *vars)
{
  char namelist[65536];

  if (cdo::stdinIsTerminal)
  {
    fprintf(stderr, "Default namelist: \n");
    fprintf(stderr, "  TYPE=0, CODE=-1, LEVEL=-1, INTERVAL=0, MEAN=0, EXTRAPOLATE=1\n");
    fprintf(stdout, "Enter namelist parameter:\n");
  }
  else
  {
    std::fseek(stdin, 0L, SEEK_END);
    long length = ftell(stdin);
    if (length == 0L) fprintf(stderr, "\n stdin not connected\n");
    std::fseek(stdin, 0L, SEEK_SET);
  }

  int i = 1;
  namelist[0] = ' ';
  int c = getchar();
  while ((c != EOF) && i < (int) (sizeof(namelist) - 1))
  {
    if ((c >= '0' && c <= '9') || (c == '-' || c == '.'))
      namelist[i++] = c;
    else if (c >= 'a' && c <= 'z')
      namelist[i++] = c;
    else if (c >= 'A' && c <= 'Z')
      namelist[i++] = std::tolower(c);
    else
      c = ' ';

    if (c == ' ' && namelist[i - 1] != ' ') namelist[i++] = c;
    c = getchar();
  }
  namelist[i] = 0;

  globs.Debug = scan_par(globs.Verbose, namelist, "debug", 0);

  if (globs.Debug)
  {
    lprintf(stderr);
    fprintf(stderr, "  Length of namelist:%4d bytes\n", (int) std::strlen(namelist));

    for (i = 0; i < (int) std::strlen(namelist); i += 60) fprintf(stderr, "  namelist[%02d]=%-60.60s\n", i, namelist + i);
    lprintf(stderr);
  }

  if (globs.Verbose)
  {
    lprintf(stdout);
    fprintf(stdout, " Namelist:\n");
  }

  globs.Type = scan_par(globs.Verbose, namelist, "type", 0);
  globs.Multi = scan_par(globs.Verbose, namelist, "multi", 0);
  globs.Mean = scan_par(globs.Verbose, namelist, "mean", 0);
  globs.OutputInterval = scan_par(globs.Verbose, namelist, "interval", MONTHLY_INTERVAL);

  if (globs.Mean >= 2) afterAbort("Namelist parameter MEAN=%d out of bounds (0:1)", globs.Mean);

  auto fileFormat = scan_par(globs.Verbose, namelist, "format", -1);
  auto gribFormat = scan_par_obsolete(namelist, "grib", 0);
  auto cdfFormat = scan_par_obsolete(namelist, "netcdf", 0);

  if (gribFormat && cdfFormat) afterAbort("GRIB or NetCDF?");

  switch (fileFormat)
  {
    case -1: ofiletype = -1; break;
    case 0: ofiletype = CDI_FILETYPE_SRV; break;
    case 1: ofiletype = CDI_FILETYPE_GRB; break;
    case 2: ofiletype = CDI_FILETYPE_NC; break;
    case 3: ofiletype = CDI_FILETYPE_EXT; break;
    case 4: ofiletype = CDI_FILETYPE_NC2; break;
    case 5: ofiletype = CDI_FILETYPE_NC5; break;
    case 6: ofiletype = CDI_FILETYPE_NC4; break;
    default: afterAbort("unknown file format %d", fileFormat);
  }

  if (gribFormat) ofiletype = CDI_FILETYPE_GRB;
  if (cdfFormat) ofiletype = CDI_FILETYPE_NC;

  int precision = scan_par(globs.Verbose, namelist, "precision", 0);
  if (precision) switch (precision)
    {
      case 8: DataType = CDI_DATATYPE_PACK8; break;
      case 16: DataType = CDI_DATATYPE_PACK16; break;
      case 24: DataType = CDI_DATATYPE_PACK24; break;
      case 32: DataType = CDI_DATATYPE_FLT32; break;
      case 64: DataType = CDI_DATATYPE_FLT64; break;
      default: afterAbort("unsupported data precision %d", precision);
    }

  globs.unitsel = scan_par(globs.Verbose, namelist, "unitsel", 0);
  globs.DayIn = scan_par(globs.Verbose, namelist, "dayinc", 30);
  globs.extrapolate = (bool) scan_par(globs.Verbose, namelist, "extrapolate", 1);
  globs.szip = (bool) scan_par(globs.Verbose, namelist, "szip", 0);

  if (globs.Multi) --globs.Multi;

  nrqh = scan_time(globs.Verbose, namelist, hours, MaxHours);
  scan_code(namelist, vars, MaxCodes, &globs.NumCodesRequest);

  scan_darray(namelist, "level", globs.LevelRequest, MaxLevel, &globs.NumLevelRequest);
  if (globs.NumLevelRequest == 1)
    if (is_equal(globs.LevelRequest[0], -1)) globs.NumLevelRequest = 0;

  if (globs.Verbose) lprintf(stdout);

  after_checkNamelist(globs);
}

static void
after_dimcalc(AfterControl &globs)
{
  if (globs.AnalysisData) globs.NumLevel = globs.NumLevelRequest;

  if (globs.Latitudes == 0)
  {
    globs.Latitudes = 2 * ((globs.Truncation * 3 + 3) / 4);
    if (globs.Truncation == 30) globs.Latitudes = 48;
  }

  if (globs.Longitudes == 0)
  {
    globs.Longitudes = globs.Latitudes * 2;
    if (globs.Truncation == 62) globs.Longitudes = 192;
  }

  globs.Waves = globs.Truncation + 1;
  globs.Fouriers = globs.Waves * 2;
  globs.DimSP = ((long) globs.Truncation + 1) * ((long) globs.Truncation + 2);
  globs.DimFC = (long) globs.Latitudes * (long) globs.Fouriers;
  globs.DimGP = (long) globs.Latitudes * (long) globs.Longitudes;
  globs.Dim3GP = globs.NumLevel * globs.DimGP;
  globs.Dim3FC = globs.NumLevel * globs.DimFC;
  globs.Dim3SP = globs.NumLevel * globs.DimSP;
  globs.HalfLevels = globs.NumLevel + 1;
  globs.DimSP_half = globs.DimSP / 2;

  if (globs.AnalysisData) fprintf(stdout, " Found Ana or Re-Ana Data\n");

  if (globs.Verbose)
  {
    fprintf(stdout, " Dimensions:\n");
    fprintf(stdout, "  Truncation        = %4d\n", globs.Truncation);
    fprintf(stdout, "  Levels            = %4d\n", globs.NumLevel);
    fprintf(stdout, "  Latitudes         = %4d\n", globs.Latitudes);
    fprintf(stdout, "  Longitudes        = %4d\n", globs.Longitudes);
    lprintf(stdout);
  }
}

/* ----------------------------------------------------------- */
/* Extract basic dimension information                         */
/* ----------------------------------------------------------- */
static void
after_precntl(AfterControl &globs, struct Variable *vars)
{
  int vertfound = 0;
  int nhzaxis = 0;
  int gridSizeMax = 0;

  auto numVars = vlistNvars(globs.ivlistID);
  auto numGrids = vlistNumGrids(globs.ivlistID);
  auto numZaxes = vlistNumZaxis(globs.ivlistID);
  auto numSteps = vlistNtsteps(globs.ivlistID);

  if (globs.Debug)
  {
    fprintf(stderr, "numVars    = %d\n", numVars);
    fprintf(stderr, "numGrids   = %d\n", numGrids);
    fprintf(stderr, "numZaxes   = %d\n", numZaxes);
    fprintf(stderr, "numSteps   = %d\n", numSteps);
  }

  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(globs.ivlistID, index);
    auto gridtype = gridInqType(gridID);
    int datasize = gridInqSize(gridID);

    if (datasize > gridSizeMax) gridSizeMax = datasize;

    if (gridtype == GRID_SPECTRAL && globs.Truncation == 0)
    {
      specGridID = gridID;
      globs.Truncation = gridInqTrunc(gridID);
    }
    else if (gridtype == GRID_GAUSSIAN && globs.Latitudes == 0)
    {
      gaussGridID = gridID;
      globs.Longitudes = gridInqXsize(gridID);
      globs.Latitudes = gridInqYsize(gridID);
    }
  }

  if (globs.Truncation == 0 && globs.Latitudes == 0) afterAbort("Unsupported file structure (no spectral or Gaussian data found)!");

  if (globs.Truncation == 0)
  {
    if (globs.Latitudes)
    {
      switch (globs.Latitudes)
      {
        case 512: globs.Truncation = 511; break;
        case 320: globs.Truncation = 213; break;
        case 192: globs.Truncation = 127; break;
        case 160: globs.Truncation = 106; break;
        case 128: globs.Truncation = 85; break;
        case 96: globs.Truncation = 63; break;
        case 94: globs.Truncation = 62; break;
        case 64: globs.Truncation = 42; break;
        case 48: globs.Truncation = 31; break;
        case 32: globs.Truncation = 21; break;
        case 1: break;
        default: fprintf(stderr, "%d Gaussian latitudes not supported.\n", globs.Latitudes);
      }
    }
  }

  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(globs.ivlistID, index);
    auto leveltype = zaxisInqType(zaxisID);
    auto numlevel = zaxisInqSize(zaxisID);
    if (numlevel > 1)
    {
      if (leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_HYBRID_HALF || leveltype == ZAXIS_PRESSURE)
      {
        if ((leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_HYBRID_HALF) && globs.nvct == 0)
        {
          nhzaxis++;
          int nvct = zaxisInqVctSize(zaxisID);
          if (numlevel != (nvct / 2 - 1) && numlevel != (nvct / 2))
          {
            if (nvct == 0)
            {
              if (numlevel != 191) cdo_warning("VCT missing for hybrid level data with %d levels!", numlevel);
            }
            else { cdo_warning("Skip %d hybrid level data with %d levels!", (nvct / 2 - 1), numlevel); }
            continue;
          }
        }
        else if ((leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_HYBRID_HALF) && globs.nvct == zaxisInqVctSize(zaxisID))
        {
          continue;
        }

        if (iVertID != -1) cdo_warning("More than %d different vertical grid structure found!", vertfound);

        vertfound++;

        if (iVertID != -1) continue;

        iVertID = zaxisID;
        globs.NumLevelFound = numlevel;
        LevelFound.resize(globs.NumLevelFound);
        for (int l = 0; l < globs.NumLevelFound; ++l) LevelFound[l] = (int) zaxisInqLevel(zaxisID, l);

        if (leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_HYBRID_HALF)
        {
          if (globs.nvct == 0)
          {
            if (zaxisInqVctSize(zaxisID))
            {
              globs.nvct = zaxisInqVctSize(zaxisID);
              globs.numHalfLevels = globs.nvct / 2;
              if ((int) globs.vct.size() != globs.nvct)
              {
                globs.vct.resize(globs.nvct);
                zaxisInqVct(zaxisID, globs.vct.data());
              }
            }
            else { afterAbort("VCT not defined in inputfile!"); }
          }

          if (numlevel == globs.numHalfLevels) { globs.NumLevelFound--; }
          else if (numlevel != (globs.numHalfLevels - 1))
          {
            afterAbort("Number of hybrid levels %d does not match VCT levels %d", numlevel, globs.numHalfLevels - 1);
          }

          if (globs.Debug)
            for (int i = 0; i < globs.numHalfLevels; ++i)
              fprintf(stderr, " vct: %4d %10.4f %10.4f\n", i, globs.vct[i], globs.vct[i + globs.numHalfLevels]);
        }

        if (leveltype == ZAXIS_PRESSURE) globs.AnalysisData = true;
      }
    }
  }

  if (nhzaxis > 0 && globs.nvct == 0) afterAbort("VCT missing!");

  globs.NumLevel = globs.NumLevelFound;

  if (specGridID != -1) globs.Spectral = true;
  if (gaussGridID != -1) globs.Gaussian = true;

  if (globs.Debug) fprintf(stderr, "   T = %3d   L = %2d\n", globs.Truncation, globs.NumLevelFound);

  if (globs.Debug) fprintf(stderr, " CODE CHECK\n");

  if (globs.Verbose)
  {
    auto instID = vlistInqVarInstitut(globs.ivlistID, 0);
    auto modelID = vlistInqVarModel(globs.ivlistID, 0);

    lprintf(stdout);
    fprintf(stdout, " Institute : ");
    if (instID == CDI_UNDEFID)
      fprintf(stdout, "unknown\n");
    else
    {
      if (institutInqLongnamePtr(instID))
        fprintf(stdout, "%s\n", institutInqLongnamePtr(instID));
      else
        fprintf(stdout, "name unknown\n");
    }

    fprintf(stdout, " Source    : ");
    if (modelID == CDI_UNDEFID)
      fprintf(stdout, "unknown\n");
    else
    {
      if (modelInqNamePtr(modelID))
      {
        if (std::strncmp(modelInqNamePtr(modelID), "ECHAM5", 6) == 0) Source = ECHAM5_Source;
        fprintf(stdout, "%s\n", modelInqNamePtr(modelID));
      }
      else
        fprintf(stdout, "name unknown\n");
    }
  }

  for (int varID = 0; varID < numVars; ++varID)
  {
    int gridID, zaxisID, timeType;
    vlistInqVar(globs.ivlistID, varID, &gridID, &zaxisID, &timeType);
    auto code = vlistInqVarCode(globs.ivlistID, varID);
    if (code <= 0 || code >= MaxCodes)
    {
      cdo_warning("Code number %d out of range, variable ignored!", code);
      continue;
    }
    auto gridtype = gridInqType(gridID);
    auto numlevel = zaxisInqSize(zaxisID);
    auto leveltype = zaxisInqType(zaxisID);

    vars[code].ivarID = varID;
    vars[code].igridID = gridID;
    vars[code].ogridID = gridID;
    vars[code].izaxisID = zaxisID;
    vars[code].ozaxisID = zaxisID;

    vars[code].detected = true;

    if (globs.Debug)
      fprintf(stderr, "Code %3d  Levels = %3d  LevelType = %3d  GridType = %3d\n", code, numlevel, leveltype, gridtype);
  }

  if (globs.Debug) fprintf(stderr, "FieldDim = %d\n", gridSizeMax);

  globs.varray.resize(gridSizeMax);

  if (globs.Debug)
    for (int code = 0; code < MaxCodes; ++code)
    {
      if (vars[code].detected) fprintf(stderr, " Detected Code %3d with %3d level\n", code, zaxisInqSize(vars[code].izaxisID));
    }
}

/*
 * -----------------------------------------------------------
 * Define output variables
 * -----------------------------------------------------------
 */
static void
after_postcntl(const AfterControl &globs, struct Variable *vars)
{
  char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  int datatype;

  if (globs.Debug) lprintf(stdout);
  if (globs.Debug)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].detected)
      {
        auto gridID = vars[code].igridID;
        auto zaxisID = vars[code].izaxisID;
        zaxisName(zaxisInqType(zaxisID), zaxistypename);
        fprintf(stderr, " Detected Code %3d  grid %-8s size %5zu  level %2d %-8s\n", code, gridNamePtr(gridInqType(gridID)),
                gridInqSize(gridID), zaxisInqSize(zaxisID), zaxistypename);
      }

  if (globs.Debug) lprintf(stdout);
  if (globs.Debug)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed) { fprintf(stderr, "   Needed Code %3d\n", code); }

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].selected)
    {
      name[0] = 0;
      longname[0] = 0;
      units[0] = 0;
      auto ivarID = vars[code].ivarID;
      auto ogridID = vars[code].ogridID;
      auto ozaxisID = vars[code].ozaxisID;

      if (ogridID == -1)
      {
        /*
        cdo_warning( "undefined grid for code %d", code);
        */
        continue;
      }
      if (ozaxisID == -1)
      {
        /*
        cdo_warning( "undefined level for code %d", code);
        */
        continue;
      }

      auto instID = vlistInqVarInstitut(globs.ivlistID, ivarID);
      auto modelID = vlistInqVarModel(globs.ivlistID, ivarID);
      auto tableID = vlistInqVarTable(globs.ivlistID, ivarID);

      vars[code].missval = vlistInqVarMissval(globs.ivlistID, ivarID);
      vars[code].samp = nullptr;

      if (DataType != -1)
        datatype = DataType;
      else
        datatype = vlistInqVarDatatype(globs.ivlistID, ivarID);

      if (vars[code].comp) { tableID = vars[code].tableID; }
      else
      {
        vlistInqVarName(globs.ivlistID, ivarID, name);
        vlistInqVarLongname(globs.ivlistID, ivarID, longname);
        vlistInqVarUnits(globs.ivlistID, ivarID, units);
      }

      if (globs.Mean != 2)
      {
        vlistDefTaxis(globs.ovlistID, globs.taxisID2);
        auto ovarID = vlistDefVar(globs.ovlistID, ogridID, ozaxisID, TIME_VARYING);
        if (globs.Mean) vlistDefVarTsteptype(globs.ovlistID, ovarID, TSTEP_AVG);
        vlistDefVarCode(globs.ovlistID, ovarID, code);
        vars[code].ovarID = ovarID;
        vlistDefVarInstitut(globs.ovlistID, ovarID, instID);
        vlistDefVarModel(globs.ovlistID, ovarID, modelID);
        vlistDefVarTable(globs.ovlistID, ovarID, tableID);
        if (name[0]) cdiDefKeyString(globs.ovlistID, ovarID, CDI_KEY_NAME, name);
        if (longname[0]) cdiDefKeyString(globs.ovlistID, ovarID, CDI_KEY_LONGNAME, longname);
        if (units[0]) cdiDefKeyString(globs.ovlistID, ovarID, CDI_KEY_UNITS, units);
        vlistDefVarDatatype(globs.ovlistID, ovarID, datatype);
        vlistDefVarMissval(globs.ovlistID, ovarID, vars[code].missval);
      }

      if (globs.Mean >= 2)
      {
        vlistDefTaxis(globs.ovlistID2, globs.taxisID2);
        auto ovarID2 = vlistDefVar(globs.ovlistID2, ogridID, ozaxisID, TIME_VARYING);
        if (globs.Mean) vlistDefVarTsteptype(globs.ovlistID2, ovarID2, TSTEP_AVG);
        vlistDefVarCode(globs.ovlistID2, ovarID2, code);
        vars[code].ovarID2 = ovarID2;
        vlistDefVarInstitut(globs.ovlistID2, ovarID2, instID);
        vlistDefVarModel(globs.ovlistID2, ovarID2, modelID);
        vlistDefVarTable(globs.ovlistID2, ovarID2, tableID);
        if (name[0]) cdiDefKeyString(globs.ovlistID2, ovarID2, CDI_KEY_NAME, name);
        if (longname[0]) cdiDefKeyString(globs.ovlistID2, ovarID2, CDI_KEY_LONGNAME, longname);
        if (units[0]) cdiDefKeyString(globs.ovlistID2, ovarID2, CDI_KEY_UNITS, units);
        vlistDefVarDatatype(globs.ovlistID2, ovarID2, datatype);
        vlistDefVarMissval(globs.ovlistID2, ovarID2, vars[code].missval);
      }
    }

  if (globs.Debug) lprintf(stdout);
  if (globs.Debug)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        auto gridID = vars[code].ogridID;
        auto zaxisID = vars[code].ozaxisID;
        zaxisName(zaxisInqType(zaxisID), zaxistypename);
        fprintf(stderr, " Selected Code %3d  grid %-8s size %5zu  level %2d %-8s\n", code, gridNamePtr(gridInqType(gridID)),
                gridInqSize(gridID), zaxisInqSize(zaxisID), zaxistypename);
      }
}

static void
after_readVct(AfterControl &globs, const char *vctfile)
{
  char line[1024];
  int nlines = 0;

  auto fp = std::fopen(vctfile, "r");
  if (fp == nullptr) cdo_sys_error("Open failed on %s", vctfile);

  while (fgets(line, 1023, fp))
  {
    if (line[0] == '#' || line[0] == '\0') continue;
    nlines++;
  }

  globs.nvct = nlines * 2;
  globs.numHalfLevels = globs.nvct / 2;
  globs.vct.resize(globs.nvct);

  std::rewind(fp);

  int i = 0;
  while (fgets(line, 1023, fp))
  {
    if (line[0] == '#' || line[0] == '\0') continue;
    int n;
    double va, vb;
    std::sscanf(line, "%d %lg %lg", &n, &va, &vb);
    globs.vct[i] = va;
    globs.vct[i + globs.numHalfLevels] = vb;
    i++;
  }
  fprintf(stdout, "  Read VCT with %d hybrid levels from file %s\n", globs.numHalfLevels - 1, vctfile);

  std::fclose(fp);
}

static void
after_variable_init(struct Variable *vars)
{
  memset(vars, 0, sizeof(struct Variable));

  vars->ivarID = -1;
  vars->ovarID = -1;
  vars->ovarID2 = -1;
  vars->izaxisID = -1;
  vars->ozaxisID = -1;
  vars->igridID = -1;
  vars->ogridID = -1;
  vars->tableID = -1;
}

static void
after_processing(AfterControl &globs, struct Variable *vars)
{
  globs.istreamID = stream_open_read_locked(ifile);

  if (ofiletype == -1) ofiletype = streamInqFiletype(globs.istreamID);

  globs.ivlistID = streamInqVlist(globs.istreamID);
  globs.taxisID = vlistInqTaxis(globs.ivlistID);
  globs.taxisID2 = taxisDuplicate(globs.taxisID);

  if (globs.Mean != 2)
  {
    globs.ostreamID = cdo_open_write(ofileidx, ofiletype);

    if (globs.szip) cdo_def_comp_type(globs.ostreamID, CDI_COMPRESS_SZIP);

    globs.ovlistID = vlistCreate();
  }

  /* ---------------- */
  /*  pre-processing  */
  /* ---------------- */
  after_precntl(globs, vars);

  /* ----------------- */
  /*  initializations  */
  /* ----------------- */
  after_setCodes(globs, vars, MaxCodes, globs.NumCodesRequest);

  if (globs.unitsel == 2)
    for (int i = 0; i < globs.NumLevelRequest; ++i) globs.LevelRequest[i] = globs.LevelRequest[i] * 1000;

  if (!globs.AnalysisData)
    for (int i = 0; i < globs.NumLevelRequest; ++i)
    {
      if ((globs.LevelRequest[i] >= 65535) && globs.unitsel && ofiletype == CDI_FILETYPE_GRB)
      {
        fprintf(stderr, "\n Level %9.2f out of range (max=65535)!\n", globs.LevelRequest[i]);
        exit(1);
      }

      if (!globs.unitsel && globs.Type >= 20 && globs.NumLevelRequest > 1 && is_equal(globs.LevelRequest[i], 0))
      {
        fprintf(stderr, "\n Level %9.2f illegal for Type %d\n", globs.LevelRequest[i], globs.Type);
        exit(1);
      }
    }

  after_setLevel(globs);

  after_dimcalc(globs);

  globs.rcoslat.resize(globs.Latitudes);
  globs.coslat.resize(globs.Latitudes);
  globs.derivationFactor.resize(globs.Latitudes);

  if (globs.Type < 50 && globs.AnalysisData)
  {
    fprintf(stderr, " ::::::::::::::::::::::::::::::::::::::::::::::\n");
    fprintf(stderr, " -> Type < 50 is not appropriate for Analysis.\n");
    fprintf(stderr, " -> Please check wether you can use Type >= 50.\n");
    fprintf(stderr, " -> Premature Exit. Sorry.\n");
    exit(1);
  }

  if (globs.Type == 10 || globs.Type == 40 || globs.Type == 60)
  {
    if (ofiletype == CDI_FILETYPE_GRB)
      afterAbort("Can't write fourier coefficients to GRIB!");
    else if (filetype_is_netcdf(ofiletype))
      afterAbort("Can't write fourier coefficients to NetCDF!");
  }

  filename = strrchr(ifile, '/');
  if (filename == 0)
    filename = ifile;
  else
    filename++;

  if (globs.Type >= 30 && globs.Type < 50
      && (vars[DIVERGENCE].selected || vars[VELOPOT].selected || vars[VORTICITY].selected || vars[STREAM].selected
          || globs.AnalysisData))
  {
    if (globs.Type == 30) globs.Type = 70;
    if (globs.Type == 40) globs.Type = 60;
    if (globs.Type == 41) globs.Type = 61;

    if (globs.AnalysisData)
      fprintf(stderr, "\n TYPE changed to %d (for analysis data)\n", globs.Type);
    else
      fprintf(stderr, "\n TYPE changed to %d (with code %d, %d, %d or %d)\n", globs.Type, DIVERGENCE, VELOPOT, VORTICITY, STREAM);
  }

  if (globs.AnalysisData)
    after_AnalysisDependencies(vars, MaxCodes);
  else
  {
    after_EchamDependencies(vars, MaxCodes, globs.Type, Source);
    vars[GEOPOTENTIAL].needed |= globs.Type >= 30 || vars[SLP].comp || vars[GEOPOTHEIGHT].comp;
  }

  // if ( vars[U_WIND].needed || vars[V_WIND].needed )
  if (vars[U_WIND].comp || vars[V_WIND].comp)
  {
    globs.dv2uv_f1.resize(globs.DimSP_half);
    globs.dv2uv_f2.resize(globs.DimSP_half);
    geninx(globs.Truncation, globs.dv2uv_f1.data(), globs.dv2uv_f2.data());
  }

  /* --------- */
  /*  Control  */
  /* --------- */

  after_defineLevel(globs, vars);

  after_defineGrid(globs, vars);

  after_postcntl(globs, vars);  // define output variables

  after_control(globs, vars);

  if (globs.ostreamID != CDO_STREAM_UNDEF) cdo_stream_close(globs.ostreamID);

  process_def_var_num(vlistNvars(globs.ivlistID));

  streamClose(globs.istreamID);
}

class Afterburner : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Afterburner",
    .operators = { { "after", AfterburnerHelp } },
    .aliases = { { "afterburner", "after" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Afterburner> registration = RegisterEntry<Afterburner>(module);

  AfterControl globs = {};
  struct Variable vars[MaxCodes + 5];

public:
  void
  init() override
  {
    lstdout = !Options::silentMode;

    globs.Verbose = Options::cdoVerbose;

    if (cdo_operator_argc() == 1) after_readVct(globs, cdo_operator_argv(0).c_str());

    for (int code = 0; code < MaxCodes + 5; ++code) after_variable_init(&vars[code]);

    after_parini(globs, vars);  // read namelist parameter

    if (CdoDefault::FileType != CDI_UNDEFID) ofiletype = CdoDefault::FileType;

    auto streamCnt = cdo_stream_cnt();
    auto numFiles = streamCnt - 1;

    ofileidx = numFiles;

    ifile = strdup(cdo_get_stream_name(0));

    globs.Nfiles = numFiles - 1;
    if (globs.Nfiles > 0)
    {
      if (globs.Multi > 0) afterAbort("Namelist parameter MULTI works only with one inputfile");

      ifiles.resize(globs.Nfiles);
      for (int i = 0; i < globs.Nfiles; ++i) ifiles[i] = cdo_get_stream_name(--numFiles);
      for (int i = 0; i < globs.Nfiles; ++i) printf("files %d %s\n", i + 1, ifiles[i]);
    }
  }

  void
  run() override
  {
    after_processing(globs, vars);
  }

  void
  close() override
  {
    FreeMean(vars);
  }
};
