/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Tinfo      tinfo           Time information
*/

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "util_string.h"

constexpr int MaxGaps = 64;
constexpr int MaxNTSM = 128;
constexpr int LimNTSM = 1024;

int iunits[] = { 1, 60, 3600, 86400, 1, 12 };

static void
print_bounds(int taxisID, int calendar)
{
  CdiDateTime vDateTime0{}, vDateTime1{};
  taxisInqVdatetimeBounds(taxisID, &vDateTime0, &vDateTime1);

  fprintf(stdout, " %s %s", date_to_string(vDateTime0.date).c_str(), time_to_string(vDateTime0.time).c_str());
  fprintf(stdout, " %s %s", date_to_string(vDateTime1.date).c_str(), time_to_string(vDateTime1.time).c_str());

  auto julianDate0 = julianDate_encode(calendar, vDateTime0);
  auto julianDate1 = julianDate_encode(calendar, vDateTime1);
  auto jdelta = julianDate_to_seconds(julianDate_sub(julianDate1, julianDate0));

  auto timeIncr = get_time_increment(jdelta, vDateTime0.date, vDateTime1.date);

  // fprintf(stdout, "  %g  %g  %g  %d", jdelta, jdelta/3600, std::fmod(jdelta,3600), timeIncr.period%3600);
  int len = fprintf(stdout, " %3ld %s%s", (long) timeIncr.period, time_units_cstr(timeIncr.units),
                    (std::abs(timeIncr.period) != 1) ? "s" : "");
  for (int i = 0; i < 11 - len; ++i) fprintf(stdout, " ");
}

static int
fill_gap(int ngaps, int (&ntsm)[MaxNTSM], int (&rangetsm)[MaxGaps][2], CdiDateTime (&vDateTimesM)[MaxGaps][MaxNTSM], int tsID,
         TimeIncrement timeIncr0, CdiDateTime vDateTime, CdiDateTime vDateTime0, int calendar, int day0, JulianDate julianDate,
         JulianDate julianDate0)
{
  int its = 0;
  int year, month, day;
  CdiDateTime nDateTime{};
  int64_t ijulinc = timeIncr0.period * iunits[(int) timeIncr0.units];

  if (ijulinc > 0 && ngaps < MaxGaps)
  {
    rangetsm[ngaps][0] = tsID;
    rangetsm[ngaps][1] = tsID + 1;

    if (timeIncr0.units == TimeUnits::MONTHS || timeIncr0.units == TimeUnits::YEARS)
    {
      its = 0;
      nDateTime = vDateTime0;
      while (true)
      {
        cdiDate_decode(nDateTime.date, &year, &month, &day);

        month += (int) ijulinc;
        adjust_month_and_year(month, year);

        if (day0 == 31) day = days_per_month(calendar, year, month);

        nDateTime.date = cdiDate_encode(year, month, day);
        if (cdiDate_get(nDateTime.date) >= cdiDate_get(vDateTime.date)) break;

        if (its < MaxNTSM)
          vDateTimesM[ngaps][its] = nDateTime;
        else if (its >= LimNTSM)
          break;

        its++;
      }
    }
    else
    {
      its = 0;
      julianDate0 = julianDate_add_seconds(julianDate0, ijulinc);
      while (julianDate_to_seconds(julianDate0) < julianDate_to_seconds(julianDate))
      {
        nDateTime = julianDate_decode(calendar, julianDate0);
        julianDate0 = julianDate_add_seconds(julianDate0, ijulinc);
        if (its < MaxNTSM)
          vDateTimesM[ngaps][its] = nDateTime;
        else if (its >= LimNTSM)
          break;

        its++;
      }
    }
    ntsm[ngaps] = its;
  }

  return its;
}

class Tinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Tinfo",
    .operators = { { "tinfo" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Tinfo> registration = RegisterEntry<Tinfo>(module);

private:
  CdoStreamID streamID{};
  int taxisID{};
  int numSteps{};
  int calendar{};

  CdiDateTime vDateTime{};
  CdiDateTime vDateTime0{};
  CdiDateTime vDateTimeFirst{};
  int tsID = 0, ntimeout{};
  int year0, month0, day0 = 0;
  int year, month, day;
  bool lforecast = false;
  TimeIncrement timeIncr{}, timeIncr0{};
  int ntsm[MaxNTSM];
  int rangetsm[MaxGaps][2];
  CdiDateTime vDateTimesM[MaxGaps][MaxNTSM]{};
  JulianDate julianDate{}, julianDate0{};
  double jdelta = 0, jdelta0 = 0;
  int arrow = 0;

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    operator_check_argc(0);

    streamID = cdo_open_read(0);
    auto vlistID = cdo_stream_inq_vlist(streamID);

    fprintf(stdout, "\n");

    taxisID = vlistInqTaxis(vlistID);
    numSteps = vlistNtsteps(vlistID);
    if (numSteps != 0)
    {
      if (numSteps == CDI_UNDEFID)
        fprintf(stdout, "   Time axis :  unlimited steps\n");
      else
        fprintf(stdout, "   Time axis :  %d step%s\n", numSteps, (numSteps == 1) ? "" : "s");

      if (taxisID != CDI_UNDEFID)
      {
        if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
        {
          auto rDateTime = taxisInqRdatetime(taxisID);
          fprintf(stdout, "     RefTime = %s %s", date_to_string(rDateTime.date).c_str(), time_to_string(rDateTime.time).c_str());

          auto unit = taxisInqTunit(taxisID);
          if (unit != CDI_UNDEFID) fprintf(stdout, "  Units = %s", tunit_to_cstr(unit));

          calendar = taxisInqCalendar(taxisID);
          if (calendar != CDI_UNDEFID) fprintf(stdout, "  Calendar = %s", calendar_to_cstr(calendar));

          if (taxisHasBounds(taxisID)) fprintf(stdout, "  Bounds = true");

          fprintf(stdout, "\n");

          if (taxisInqType(taxisID) == TAXIS_FORECAST)
          {
            auto fDateTime = taxisInqFdatetime(taxisID);
            fprintf(stdout, "     Forecast RefTime = %s", datetime_to_string(fDateTime).c_str());

            unit = taxisInqForecastTunit(taxisID);
            if (unit != CDI_UNDEFID) fprintf(stdout, "  Units = %s", tunit_to_cstr(unit));

            fprintf(stdout, "\n");

            lforecast = true;
          }
        }
      }

      calendar = taxisInqCalendar(taxisID);

      fprintf(stdout, "\n");
      fprintf(stdout, "         Verification Time              ");
      if (lforecast) fprintf(stdout, " Forecast Reference Time     ");
      if (taxisHasBounds(taxisID)) fprintf(stdout, " lower bound          upper bound");
      fprintf(stdout, "\n");

      fprintf(stdout, "Timestep YYYY-MM-DD hh:mm:ss   Increment");
      if (lforecast) fprintf(stdout, " YYYY-MM-DD hh:mm:ss   Period");
      if (taxisHasBounds(taxisID)) fprintf(stdout, " YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  Difference");
      fprintf(stdout, "\n");
    }
  }

  void
  run() override
  {
    if (numSteps == 0) return;

    int numGabs = 0;
    tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;

      vDateTime = taxisInqVdatetime(taxisID);

      cdiDate_decode(vDateTime.date, &year, &month, &day);

      fprintf(stdout, "%6d  %s %s", tsID + 1, date_to_string(vDateTime.date).c_str(), time_to_string(vDateTime.time).c_str());

      if (tsID)
      {
        cdiDate_decode(vDateTime0.date, &year0, &month0, &day0);

        julianDate0 = julianDate_encode(calendar, vDateTime0);
        julianDate = julianDate_encode(calendar, vDateTime);
        jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));

        timeIncr = get_time_increment(jdelta, vDateTime0.date, vDateTime.date);

        int len = fprintf(stdout, " %3ld %s%s", (long) timeIncr.period, time_units_cstr(timeIncr.units),
                          (std::abs(timeIncr.period) != 1) ? "s" : "");
        for (int i = 0; i < 11 - len; ++i) fprintf(stdout, " ");
      }
      else
      {
        vDateTimeFirst = vDateTime;
        fprintf(stdout, "   --------");
      }

      if (lforecast)
      {
        auto fDateTime = taxisInqFdatetime(taxisID);
        fprintf(stdout, " %s", datetime_to_string(fDateTime).c_str());

        auto fc_period = taxisInqForecastPeriod(taxisID);
        fprintf(stdout, " %7g", fc_period);
      }

      if (taxisHasBounds(taxisID)) print_bounds(taxisID, calendar);

      if (tsID > 1 && timeIncr != timeIncr0)
      {
        int its{};
        if (tsID == 2 && (jdelta0 > jdelta))
        {
          jdelta0 = jdelta;
          timeIncr0 = timeIncr;

          its = fill_gap(numGabs, ntsm, rangetsm, vDateTimesM, 1, timeIncr0, vDateTimeFirst, vDateTime, calendar, day, julianDate0,
                         julianDate_encode(calendar, vDateTimeFirst));

          arrow = '^';
        }
        else
        {
          its = fill_gap(numGabs, ntsm, rangetsm, vDateTimesM, tsID, timeIncr0, vDateTime, vDateTime0, calendar, day0, julianDate,
                         julianDate0);

          arrow = '<';

          if (its == 0 && timeIncr.period < 0)
          {
            its = -1;
            vDateTime = vDateTime0;
          }
        }

        if (its > 0)
        {
          numGabs++;
          if (Options::cdoVerbose)
            fprintf(stdout, "  %c--- Gap %d, missing %s%d timestep%s", arrow, numGabs, (its >= LimNTSM) ? "more than " : "", its,
                    (its != 1) ? "s" : "");
        }
        else if (its < 0)
        {
          if (Options::cdoVerbose) fprintf(stdout, "  %c--- Wrong date/time information, negative increment!", arrow);
        }
      }

      if (tsID == 1)
      {
        jdelta0 = jdelta;
        timeIncr0 = timeIncr;
      }

      fprintf(stdout, "\n");

      vDateTime0 = vDateTime;

      tsID++;
    }

    fprintf(stdout, "\n");

    fprintf(stdout, " Start date          : %s %s\n", date_to_string(vDateTimeFirst.date).c_str(),
            time_to_string(vDateTimeFirst.time).c_str());
    fprintf(stdout, " End date            : %s %s\n", date_to_string(vDateTime.date).c_str(),
            time_to_string(vDateTime.time).c_str());

    fprintf(stdout, " Increment           : %3ld %s%s\n", (long) timeIncr0.period, time_units_cstr(timeIncr0.units),
            (timeIncr0.period != 1) ? "s" : "");
    fprintf(stdout, " Number of timesteps : %d\n", tsID);
    fprintf(stdout, " Gaps identified     : %d\n", numGabs);

    if (Options::cdoVerbose && numGabs)
    {
      fprintf(stdout, "\nFound potentially %d gap%s in the time series", numGabs, (numGabs != 1) ? "s" : "");
      if (numGabs >= MaxGaps)
      {
        numGabs = MaxGaps;
        fprintf(stdout, ", here are the first %d", numGabs);
      }
      fprintf(stdout, ":\n");
      for (int igap = 0; igap < numGabs; ++igap)
      {
        fprintf(stdout, "  Gap %d between timestep %d and %d, missing %d timestep%s", igap + 1, rangetsm[igap][0],
                rangetsm[igap][1], ntsm[igap], (ntsm[igap] != 1) ? "s" : "");
        if (ntsm[igap] >= MaxNTSM)
        {
          ntsm[igap] = MaxNTSM;
          fprintf(stdout, ", here are the first %d", ntsm[igap]);
        }
        fprintf(stdout, ":\n");

        ntimeout = 0;
        for (int its = 0; its < ntsm[igap]; ++its)
        {
          if (ntimeout == 4)
          {
            ntimeout = 0;
            fprintf(stdout, "\n");
          }

          vDateTime = vDateTimesM[igap][its];
          fprintf(stdout, "  %s %s", date_to_string(vDateTime.date).c_str(), time_to_string(vDateTime.time).c_str());

          ntimeout++;
          tsID++;
        }
        fprintf(stdout, "\n");
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
