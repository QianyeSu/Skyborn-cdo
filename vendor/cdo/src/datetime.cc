/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "datetime.h"
#include "param_conversion.h"
#include "printinfo.h"
#include "util_string.h"

TimeStat CDO_Timestat_Date = TimeStat::UNDEF;
bool CDO_Ignore_Time_Bounds = false;
bool CDO_Use_Time_Bounds = false;

const char *
time_units_cstr(TimeUnits timeUnit)
{
  if (timeUnit == TimeUnits::SECONDS) return "second";
  if (timeUnit == TimeUnits::MINUTES) return "minute";
  if (timeUnit == TimeUnits::HOURS) return "hour";
  if (timeUnit == TimeUnits::DAYS) return "day";
  if (timeUnit == TimeUnits::MONTHS) return "month";
  if (timeUnit == TimeUnits::YEARS) return "year";

  return NULL;
}

void
set_timestat_date(std::string const &p_optarg)
{
  TimeStat timestatdate = TimeStat::UNDEF;

  // clang-format off
  if      (p_optarg == "first")   timestatdate = TimeStat::FIRST;
  else if (p_optarg == "last")    timestatdate = TimeStat::LAST;
  else if (p_optarg == "middle")  timestatdate = TimeStat::MEAN;
  else if (p_optarg == "midhigh") timestatdate = TimeStat::MIDHIGH;
  // clang-format on

  if (timestatdate == TimeStat::UNDEF) cdo_abort("option --%s: unsupported argument: %s", "timestat_date", p_optarg);

  CDO_Timestat_Date = timestatdate;
}

static void
get_timestat_date(TimeStat &tstatDate)
{
  auto envString = getenv_string("CDO_TIMESTAT_DATE");
  if (envString.empty()) envString = getenv_string("RUNSTAT_DATE");
  if (envString.size())
  {
    TimeStat envDate = TimeStat::UNDEF;
    auto envstrl = string_to_lower(envString);

    // clang-format off
    if      (envstrl == "first")   envDate = TimeStat::FIRST;
    else if (envstrl == "last")    envDate = TimeStat::LAST;
    else if (envstrl == "middle")  envDate = TimeStat::MEAN;
    else if (envstrl == "midhigh") envDate = TimeStat::MIDHIGH;
    // clang-format on

    if (envDate != TimeStat::UNDEF)
    {
      tstatDate = envDate;
      if (Options::cdoVerbose) cdo_print("Set CDO_TIMESTAT_DATE to %s", envString);
    }
  }
}

void
DateTimeList::init()
{
  static bool dateTimeInit{ false };
  if (!dateTimeInit && CDO_Timestat_Date == TimeStat::UNDEF) get_timestat_date(CDO_Timestat_Date);
  dateTimeInit = true;
}

CdiDateTime
DateTimeList::vDateTime(int tsID)
{
  if (tsID < 0 || (size_t) tsID >= m_size) cdo_abort("Internal error: tsID out of bounds!");

  return m_dtInfo[tsID].c;
}

void
DateTimeList::shift()
{
  for (size_t inp = 0; inp < m_size - 1; ++inp) m_dtInfo[inp] = m_dtInfo[inp + 1];
}

void
DateTimeList::taxis_inq_timestep(int taxisID, int tsID)
{
  auto nalloc = m_dtInfo.size();
  if ((size_t) tsID >= nalloc) { m_dtInfo.resize((nalloc >= 1024) ? nalloc + 512 : nalloc * 2); }

  if ((size_t) tsID >= m_size) m_size = (size_t) tsID + 1;

  m_dtInfo[tsID].v = taxisInqVdatetime(taxisID);
  m_dtInfo[tsID].c = m_dtInfo[tsID].v;

  if (tsID == 0)
  {
    if (m_hasBounds == -1) m_hasBounds = CDO_Ignore_Time_Bounds ? 0 : taxisHasBounds(taxisID);
    if (m_calendar == -1) m_calendar = taxisInqCalendar(taxisID);
  }

  if (m_hasBounds)
  {
    taxisInqVdatetimeBounds(taxisID, &(m_dtInfo[tsID].b[0]), &(m_dtInfo[tsID].b[1]));

    auto time = cdiTime_get(m_dtInfo[tsID].b[1].time);
    if (CDO_Use_Time_Bounds && time == 0 && cdiDateTime_isEQ(m_dtInfo[tsID].v, m_dtInfo[tsID].b[1]))
    {
      auto julianDate1 = julianDate_encode(m_calendar, m_dtInfo[tsID].b[0]);
      auto julianDate2 = julianDate_encode(m_calendar, m_dtInfo[tsID].b[1]);

      if (julianDate_to_seconds(julianDate1) < julianDate_to_seconds(julianDate2))
      {
        auto julianDate = julianDate_add_seconds(julianDate2, -1);
        m_dtInfo[tsID].c = julianDate_decode(m_calendar, julianDate);
      }
    }
  }
  else
  {
    cdiDateTime_init(&m_dtInfo[tsID].b[0]);
    cdiDateTime_init(&m_dtInfo[tsID].b[1]);
  }
}

void
DateTimeList::taxis_set_next_timestep(int taxisID)
{
  int tsID = m_size;
  this->taxis_inq_timestep(taxisID, tsID);
}

void
DateTimeList::taxis_def_timestep(int taxisID, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= m_size) cdo_abort("Internal error; tsID out of bounds!");

  taxisDefVdatetime(taxisID, m_dtInfo[tsID].v);
  if (m_hasBounds) taxisDefVdatetimeBounds(taxisID, m_dtInfo[tsID].b[0], m_dtInfo[tsID].b[1]);
}

void
DateTimeList::mean(int numSteps)
{
  if (numSteps % 2 == 0)
  {
#ifdef TEST_DTLIST_MEAN
    auto julianDate0 = julianDate_encode(m_calendar, dtInfo[0].v);

    double seconds = 0.0;
    for (int i = 1; i < numSteps; ++i)
    {
      auto julianDate = julianDate_encode(m_calendar, dtInfo[i].v);
      seconds += julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
    }

    auto julianDate = julianDate_add_seconds(julianDate0, std::lround(seconds / numSteps));
    timestat.v = julianDate_decode(calendar, julianDate);
#else
    auto julianDate1 = julianDate_encode(m_calendar, m_dtInfo[numSteps / 2 - 1].v);
    auto julianDate2 = julianDate_encode(m_calendar, m_dtInfo[numSteps / 2].v);

    auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
    auto julianDatem = julianDate_add_seconds(julianDate1, std::lround(seconds));
    m_timestat.v = julianDate_decode(m_calendar, julianDatem);
#endif
  }
  else { m_timestat.v = m_dtInfo[numSteps / 2].v; }
}

void
DateTimeList::midhigh(int numSteps)
{
  m_timestat.v = m_dtInfo[numSteps / 2].v;
}

void
DateTimeList::stat_taxis_def_timestep(int taxisID, int numSteps)
{
  if ((size_t) numSteps > m_size) cdo_abort("Internal error; unexpected numSteps=%d (limit=%ld)!", numSteps, m_size);

  if (CDO_Timestat_Date != TimeStat::UNDEF) m_stat = CDO_Timestat_Date;

  // clang-format off
  if      (m_stat == TimeStat::MEAN)    mean(numSteps);
  else if (m_stat == TimeStat::MIDHIGH) midhigh(numSteps);
  else if (m_stat == TimeStat::FIRST)   m_timestat.v = m_dtInfo[0].v;
  else if (m_stat == TimeStat::LAST)    m_timestat.v = m_dtInfo[numSteps - 1].v;
  else cdo_abort("Internal error; implementation missing for timestat=%d", (int)m_stat);
  // clang-format on

  if (m_hasBounds)
  {
    m_timestat.b[0] = m_dtInfo[0].b[0];
    m_timestat.b[1] = m_dtInfo[numSteps - 1].b[1];
  }
  else
  {
    m_timestat.b[0] = m_dtInfo[0].v;
    m_timestat.b[1] = m_dtInfo[numSteps - 1].v;
  }

  if (m_year) { m_timestat.v.date.year = m_year; }

  taxisDefVdatetime(taxisID, m_timestat.v);
  // if (m_hasBounds)
  {
    taxisDefVdatetimeBounds(taxisID, m_timestat.b[0], m_timestat.b[1]);
  }
}

void
DateTimeList::stat_taxis_def_timestep(int taxisID)
{
  int numSteps = m_size;
  this->stat_taxis_def_timestep(taxisID, numSteps);
}

CdiDateTime
datetime_avg(int calendar, int ndates, std::vector<CdiDateTime> const &cdiDateTimes)
{
  if (ndates % 2 == 0)
  {
    auto julianDate1 = julianDate_encode(calendar, cdiDateTimes[ndates / 2 - 1]);
    auto julianDate2 = julianDate_encode(calendar, cdiDateTimes[ndates / 2]);

    auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
    auto julianDatem = julianDate_add_seconds(julianDate1, std::lround(seconds));
    return julianDate_decode(calendar, julianDatem);
  }
  else { return cdiDateTimes[ndates / 2]; }
}

void
adjust_month_and_year(int &month, int &year)
{
  // clang-format off
  while (month > 12) { month -= 12; year++; }
  while (month <  1) { month += 12; year--; }
  // clang-format on
}

double
delta_time_step_0(int tsID, int calendar, CdiDateTime const &vDateTime, JulianDate &julianDate0, double &deltat1)
{
  double zj = 0.0;

  auto julianDate = julianDate_encode(calendar, vDateTime);

  if (tsID == 0) { julianDate0 = julianDate; }
  else
  {
    auto deltat = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
    if (tsID == 1) deltat1 = deltat;
    zj = deltat / deltat1;
  }

  return zj;
}

TimeIncrement
get_time_increment(double jdelta, CdiDate vDate0, CdiDate vDate1)
{
  auto seconds = (jdelta < 0.0) ? (int64_t) (jdelta - 0.5) : (int64_t) (jdelta + 0.5);

  int sign = 1;
  if (seconds < 0)
  {
    std::swap(vDate0, vDate1);
    seconds = -seconds;
    sign = -1;
  }

  int year0, month0, day0;
  cdiDate_decode(vDate0, &year0, &month0, &day0);
  int year1, month1, day1;
  cdiDate_decode(vDate1, &year1, &month1, &day1);

  auto deltay = year1 - year0;
  auto deltam = deltay * 12 + (month1 - month0);
  if (deltay == 0) deltay = 1;
  if (deltam == 0) deltam = 1;

  constexpr int secPerDay = 3600 * 24;
  TimeIncrement timeIncr;
  if (seconds >= (secPerDay * 30 * 12)) { timeIncr = { deltay, TimeUnits::YEARS }; }
  else if (seconds >= (secPerDay * 30) && seconds / (secPerDay * 30) < 12) { timeIncr = { deltam, TimeUnits::MONTHS }; }
  else if (seconds >= secPerDay && seconds / secPerDay < 32)
  {
    timeIncr = { seconds / secPerDay, TimeUnits::DAYS };
    if (timeIncr.period >= 27 && deltam == 1) timeIncr = { 1, TimeUnits::MONTHS };
  }
  else if (seconds >= 3600 && seconds % 3600 == 0) { timeIncr = { seconds / 3600, TimeUnits::HOURS }; }
  else if (seconds >= 60 && seconds % 60 == 0) { timeIncr = { seconds / 60, TimeUnits::MINUTES }; }
  else { timeIncr = { seconds, TimeUnits::SECONDS }; }

  timeIncr.period *= sign;

  return timeIncr;
}

void
check_time_increment(int tsID, int calendar, CdiDateTime const &vDateTime, CheckTimeIncr &checkTimeIncr)
{
  auto julianDate = julianDate_encode(calendar, vDateTime);

  if (tsID)
  {
    auto jdeltat = julianDate_to_seconds(julianDate_sub(julianDate, checkTimeIncr.julianDate0));
    auto timeIncr = get_time_increment(jdeltat, checkTimeIncr.vDate0, vDateTime.date);

    if (tsID == 1) checkTimeIncr.timeIncr = timeIncr;

    if (checkTimeIncr.printWarning
        && (timeIncr.period != checkTimeIncr.timeIncr.period || timeIncr.units != checkTimeIncr.timeIncr.units))
    {
      checkTimeIncr.printWarning = false;
      cdo_warning("Time increment in step %d (%lld%s) differs from step 1 (%lld%s)!"
                  " Set parameter equal=false for unequal time increments!",
                  tsID + 1, timeIncr.period, time_units_cstr(timeIncr.units), checkTimeIncr.timeIncr.period,
                  time_units_cstr(checkTimeIncr.timeIncr.units));
    }

    /*
    if (Options::cdoVerbose)
      std::fprintf(stdout, "Timestep: %d  increment: %3ld %s%s\n",
              tsID+1, (long) incrPeriod, tunits[(int)incrUnits], (std::abs(incrPeriod) != 1) ? "s" : "");
    */
  }

  checkTimeIncr.vDate0 = vDateTime.date;
  checkTimeIncr.julianDate0 = julianDate;
}

int
decode_month(CdiDate const &date)
{
  return date.month;
}

int
decode_month_and_day(CdiDate const &date)
{
  return (date.month * 100 + date.day);
}

int
decode_day_of_year(CdiDate const &date)
{
  if (date.day < 1 || date.day > 31) { return 0; }
  if (date.month < 1 || date.month > 12) { return 0; }

  return (date.month - 1) * 31 + date.day;
}

int
decode_hour_of_year(CdiDateTime const &cdiDateTime, int maxHours)
{
  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  int hourOfDay = 0;
  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
    hourOfDay = ((month - 1) * 31 + day - 1) * 25 + hour + 1;

  if (hourOfDay < 0 || hourOfDay >= maxHours)
    cdo_abort("Hour of year %d out of range (%s)!", hourOfDay, datetime_to_string(cdiDateTime));

  return hourOfDay;
}

int
decode_hour_of_day(CdiDateTime const &cdiDateTime, int maxHours)
{
  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  int hourOfDay = 0;
  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24) hourOfDay = hour + 1;

  if (hourOfDay < 0 || hourOfDay >= maxHours)
    cdo_abort("Hour of day %d out of range (%s)!", hourOfDay, datetime_to_string(cdiDateTime));

  return hourOfDay;
}

int
decode_minute_of_day(CdiDateTime const &cdiDateTime, int maxMinutes)
{
  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  int minuteOfDay = 0;
  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
  {
    if (minute >= 0 && minute < 60) minuteOfDay = hour * 60 + minute + 1;
  }

  if (minuteOfDay < 0 || minuteOfDay >= maxMinutes)
    cdo_abort("Minute of day %d out of range (%s)!", minuteOfDay, datetime_to_string(cdiDateTime));

  return minuteOfDay;
}

void
set_date_time(CdiDateTime &datetime1, CdiDateTime datetime2)
{
  if (datetime2.date.month == 12) datetime2.date.year -= 1;

  if (cdiDate_get(datetime2.date) > cdiDate_get(datetime1.date)) datetime1 = datetime2;
}

static void
get_timeunits(std::string const &unitsStr, int &incrPeriod, int &incrUnits, int &timeUnits)
{
  // clang-format off
  static std::unordered_map<std::string, std::pair<int, int>> unitsMap
      = { { "seconds", {     1, TUNIT_SECOND } },
          { "minutes", {    60, TUNIT_MINUTE } },
          { "hours",   {  3600, TUNIT_HOUR } },
          { "days",    { 86400, TUNIT_DAY } },
          { "months",  {     1, TUNIT_MONTH } },
          { "years",   {    12, TUNIT_YEAR } } };
  // clang-format on

  for (auto const &entry : unitsMap)
  {
    if (entry.first.starts_with(unitsStr))
    {
      incrUnits = entry.second.first;
      timeUnits = entry.second.second;
      // clang-format off
      if (timeUnits == TUNIT_HOUR)
      {
        if      (incrPeriod ==  3) { incrPeriod = 1; incrUnits = 10800; timeUnits = TUNIT_3HOURS;  }
        else if (incrPeriod ==  6) { incrPeriod = 1; incrUnits = 21600; timeUnits = TUNIT_6HOURS;  }
        else if (incrPeriod == 12) { incrPeriod = 1; incrUnits = 43200; timeUnits = TUNIT_12HOURS; }
      }
      // clang-format on

      return;
    }
  }

  cdo_abort("Time units >%s< unsupported!", unitsStr);
}

CdiDate
decode_datestring(std::string const &dateString)
{
  if (std::strchr(dateString.c_str() + 1, '-'))
  {
    int year = 1, month = 1, day = 1;
    std::sscanf(dateString.c_str(), "%d-%d-%d", &year, &month, &day);
    return cdiDate_encode(year, month, day);
  }
  else { return cdiDate_set(parameter_to_long(dateString)); }
}

CdiTime
decode_timestring(std::string const &timeString)
{
  if (std::strchr(timeString.c_str(), ':'))
  {
    int hour = 0, minute = 0;
    double fseconds = 0.0;
    std::sscanf(timeString.c_str(), "%d:%d:%lf", &hour, &minute, &fseconds);
    int second = (int) fseconds;
    int ms = (fseconds - second) * 1000;
    return cdiTime_encode(hour, minute, second, ms);
  }
  else { return cdiTime_set(parameter_to_int(timeString)); }
}

void
decode_timeunits(std::string const &timeUnitsString, int &incrPeriod, int &incrUnits, int &timeUnits)
{
  incrPeriod = 0;
  incrUnits = 0;
  timeUnits = 0;

  char *pUnits = nullptr;
  auto fperiod = std::strtod(timeUnitsString.c_str(), &pUnits);
  if (pUnits != timeUnitsString.c_str()) incrPeriod = std::lround(fperiod);

  if (pUnits) get_timeunits(pUnits, incrPeriod, incrUnits, timeUnits);
}

/**
 * Computes the day-of-year correspnding a given Gregorian date.
 *
 * @param date a Gregorian date in the form YYYYMMDD
 *
 * @return the day-of-year
 */
int
day_of_year(int calendar, int64_t date)
{
  int year = date / 10000;
  return (date_to_julday(calendar, date) - date_to_julday(calendar, cdiEncodeDate(year, 1, 1)) + 1);
}
