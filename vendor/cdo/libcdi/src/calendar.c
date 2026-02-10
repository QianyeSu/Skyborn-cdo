/* DO NOT REMOVE the config.h include file under any circumstances,
 * it's very much needed on some platforms */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
/* DO NOT REMOVE the above config.h include file under any
 * circumstances as long as it's the autoconf configuration header
 * used to build this package. When it's missing on some platforms,
 * some poor person has to do long, tedious debugging sessions, where
 * struct offsets almost imperceptibly change from one file to the
 * next to find out what happened */

#include <stdio.h>

#include "calendar.h"

static const int month_360[12] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };
static const int month_365[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
static const int month_366[12] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

static const int *
get_dayspermonth_array(int daysPerYear)
{
  // clang-format off
  return  (daysPerYear == 360) ? month_360 :
          (daysPerYear == 365) ? month_365 :
          (daysPerYear == 366) ? month_366 : NULL;
  // clang-format on
}

int
days_per_month(int calendar, int year, int month)
{
  int daysPerYear = calendar_dpy(calendar);
  const int *daysPerMonthArray = (daysPerYear == 360) ? month_360 : ((daysPerYear == 365) ? month_365 : month_366);

  int daysPerMonth = (month >= 1 && month <= 12) ? daysPerMonthArray[month - 1] : 0;

  if (daysPerYear == 0 && month == 2) daysPerMonth = ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) ? 29 : 28;

  return daysPerMonth;
}

int
days_per_year(int calendar, int year)
{
  int daysPerYear = calendar_dpy(calendar);
  if (daysPerYear == 0)
  {
    if (year == 1582 && (calendar == CALENDAR_STANDARD || calendar == CALENDAR_GREGORIAN))
      daysPerYear = 355;
    else if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)
      daysPerYear = 366;
    else
      daysPerYear = 365;
  }

  return daysPerYear;
}

void
decode_calday(int daysPerYear, int64_t days, int *year, int *month, int *day)
{
  *year = (int) ((days - 1) / daysPerYear);
  days -= (*year * daysPerYear);

  const int *daysPerMonthArray = get_dayspermonth_array(daysPerYear);

  int i = 0;
  if (daysPerMonthArray)
    for (i = 0; i < 12; i++)
    {
      if (days > daysPerMonthArray[i])
        days -= daysPerMonthArray[i];
      else
        break;
    }

  *month = i + 1;
  *day = (int) days;
}

int64_t
encode_calday(int daysPerYear, int year, int month, int day)
{
  int64_t rval = (int64_t) daysPerYear * year + day;

  const int *daysPerMonthArray = get_dayspermonth_array(daysPerYear);

  if (daysPerMonthArray)
    for (int i = 0; i < month - 1; i++) rval += daysPerMonthArray[i];

  return rval;
}
