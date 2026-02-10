/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef DATETIME_H
#define DATETIME_H

#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>

#include "julian_date.h"

enum struct TimeUnits
{
  SECONDS = 0,
  MINUTES,
  HOURS,
  DAYS,
  MONTHS,
  YEARS
};

enum struct TimeStat
{
  UNDEF,
  FIRST,
  LAST,
  MEAN,
  MIDHIGH,
};

struct DateTimeInfo
{
  CdiDateTime c{};     // corrected verification time
  CdiDateTime v{};     // verification time
  CdiDateTime b[2]{};  // time bounds
};

class TimeIncrement
{
public:
  int64_t period = 0;
  TimeUnits units = TimeUnits::SECONDS;

  TimeIncrement() {}
  TimeIncrement(int64_t _period, TimeUnits _units) : period{ _period }, units{ _units } {}

  bool
  operator==(const TimeIncrement &timeIncr) const
  {
    return (period == timeIncr.period && units == timeIncr.units);
  }

  bool
  operator!=(const TimeIncrement &timeIncr) const
  {
    return (period != timeIncr.period || units != timeIncr.units);
  }
};

struct CheckTimeIncr
{
  JulianDate julianDate0{};
  TimeIncrement timeIncr{};
  CdiDate vDate0{};
  bool printWarning{ true };
};

// clang-format off
class  // DateTimeList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
DateTimeList
// clang-format on
{
public:
  // clang-format off
  DateTimeList() { init(); m_dtInfo.resize(2); }
  void set_year(int year) { m_year = year; }
  void set_stat(TimeStat stat) { m_stat = stat; }
  void set_calendar(int calendar) { m_calendar = calendar; }
  size_t size(){ return m_size; }
  std::vector<DateTimeInfo> const& info() { return m_dtInfo; }
  // clang-format on
  CdiDateTime vDateTime(int tsID);
  void shift();

  void taxis_set_next_timestep(int taxisID);
  void taxis_inq_timestep(int taxisID, int tsID);
  void taxis_def_timestep(int taxisID, int tsID);
  void stat_taxis_def_timestep(int taxisID, int numSteps);
  void stat_taxis_def_timestep(int taxisID);

private:
  size_t m_size{ 0 };
  int m_year{ 0 };
  int m_hasBounds{ -1 };
  int m_calendar{ -1 };
  TimeStat m_stat{ TimeStat::LAST };
  DateTimeInfo m_timestat{};
  std::vector<DateTimeInfo> m_dtInfo{};

  void init();
  void mean(int nsteps);
  void midhigh(int nsteps);
};

CdiDateTime datetime_avg(int calendar, int ndates, std::vector<CdiDateTime> const &cdiDateTimes);
void set_timestat_date(std::string const &p_optarg);

void adjust_month_and_year(int &month, int &year);

double delta_time_step_0(int tsID, int calendar, CdiDateTime const &vDateTime, JulianDate &juldate0, double &deltat1);

TimeIncrement get_time_increment(double jdelta, CdiDate vDate0, CdiDate vDate1);

void check_time_increment(int tsID, int calendar, CdiDateTime const &vDateTime, CheckTimeIncr &checkTimeIncr);

int decode_month(CdiDate const &date);
int decode_month_and_day(CdiDate const &date);

int decode_day_of_year(CdiDate const &date);
int decode_hour_of_year(CdiDateTime const &cdiDateTime, int maxHours);
int decode_hour_of_day(CdiDateTime const &cdiDateTime, int maxHours);
int decode_minute_of_day(CdiDateTime const &cdiDateTime, int maxMinutes);

void set_date_time(CdiDateTime &datetime1, CdiDateTime datetime2);

const char *time_units_cstr(TimeUnits timeUnit);

CdiDate decode_datestring(std::string const &dateString);
CdiTime decode_timestring(std::string const &timeString);
void decode_timeunits(std::string const &timeUnitsString, int &incrPeriod, int &incrUnits, int &timeUnits);

int day_of_year(int calendar, int64_t date);

#endif /* DATETIME_H */
