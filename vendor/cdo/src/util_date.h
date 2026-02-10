#ifndef UTIL_DATE_H
#define UTIL_DATE_H

#include <cstring>
#include <cstdio>

#include "cdi_datetime.h"

enum
{
  CMP_HOUR = 4,
  CMP_DAY = 6,
  CMP_MONTH = 8,
  CMP_YEAR = 10,
  CMP_DATE = 31  // xxxxxxxxxxxxxxxxYYYYMMDDhhmmss allocate DATE_LEN+1 !!!!
};

#define SET_DATE(dtstr, date, time) (snprintf(dtstr, sizeof(dtstr), "%*ld%*d", CMP_DATE - 6, (long) date, 6, time))

inline bool
date_is_neq(CdiDateTime const &dateTime1, CdiDateTime const &dateTime2, int compareDate)
{
  char dateStr1[CMP_DATE + 1], dateStr2[CMP_DATE + 1];
  SET_DATE(dateStr1, cdiDate_get(dateTime1.date), cdiTime_get(dateTime1.time));
  SET_DATE(dateStr2, cdiDate_get(dateTime2.date), cdiTime_get(dateTime2.time));
  return (memcmp(dateStr1, dateStr2, CMP_DATE - compareDate) != 0);
}

#endif
