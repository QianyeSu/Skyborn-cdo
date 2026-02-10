/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

#include <algorithm>
#include <assert.h>

#include "cdi.h"
#include "julian_date.h"

#include "process_int.h"
#include "datetime.h"
#include "ecautil.h"

/**
 * Counts the number of nonumMissValsing values. The result of the operation
 * is computed according to the following rules:
 *
 * field1  field2  mode  result
 * a       b       0     a + 1
 * a       miss    0     a
 * miss    b       0     1
 * miss    miss    0     0
 *
 * a       b       1     a + 1
 * a       miss    1     0
 * miss    b       1     1
 * miss    miss    1     0
 *
 * a       b       n     b < n ? a : a + n
 * a       miss    n     a
 * miss    b       n     b < n ? 0 : b
 * miss    miss    n     0
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 * @param mode   the counting mode, must be an exact mathematical integer
 */
static void
count(Field &field1, Field const &field2, double mode)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals)
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (fp_is_equal(array2[i], missval2))
            {
              if (is_equal(mode, 1.0) || fp_is_equal(array1[i], missval1)) array1[i] = 0.0;
              continue;
            }

          if (fp_is_not_equal(array1[i], missval1))
            {
              if (is_equal(mode, 0.0) || is_equal(mode, 1.0))
                array1[i] += 1.0;
              else if (fp_is_equal(array2[i], mode) || array2[i] > mode)
                array1[i] += array2[i];
            }
          else
            {
              if (is_equal(mode, 0.0) || is_equal(mode, 1.0))
                array1[i] = 1.0;
              else if (array2[i] < mode)
                array1[i] = 0.0;
              else
                array1[i] = array2[i];
            }
        }

      field_num_mv(field1);
    }
  else
    {
      if (field2.numMissVals)
        {
          for (size_t i = 0; i < len; ++i)
            {
              if (fp_is_equal(array2[i], missval2))
                {
                  if (is_equal(mode, 1.0)) array1[i] = 0.0;
                }
              else if (is_equal(mode, 0.0) || is_equal(mode, 1.0))
                array1[i] += 1.0;
              else if (array2[i] >= mode)
                array1[i] += array2[i];
            }
        }
      else
        {
          for (size_t i = 0; i < len; ++i)
            {
              if (is_equal(mode, 0.0) || is_equal(mode, 1.0))
                array1[i] += 1.0;
              else if (array2[i] >= mode)
                array1[i] += array2[i];
            }
        }
    }
}

/**
 * Selects all field elements that compare to the corresponding
 * element of a reference field. The result of the operation is
 * computed according to the following rules:
 *
 * field1  field2  result
 * a       b       comp(a, b) ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1  the input field, also holds the result
 * @param field2  the reference field
 * @param compare the comparator
 */
static void
selcomp(Field &field1, Field const &field2, int (*compare)(double, double))
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
    {
      for (size_t i = 0; i < len; ++i)
        if (fp_is_equal(array1[i], missval1) || fp_is_equal(array2[i], missval2) || !compare(array1[i], array2[i]))
          array1[i] = missval1;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (!compare(array1[i], array2[i])) array1[i] = missval1;
    }

  field_num_mv(field1);
}

/**
 * Selects all field elements that compare to a certain reference
 * value. The result of the operation is computed according to the
 * following rules:
 *
 * field  c      result
 * a      c      comp(a, c) ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field   the input field, also holds the result
 * @param c       the refence value
 * @param compare the comparator
 */
static void
selcompc(Field &field, double c, int (*compare)(double, double))
{
  auto missval = field.missval;
  auto &array = field.vec_d;

  auto len = field.size;

  if (fp_is_equal(c, missval))
    {
      for (size_t i = 0; i < len; ++i) array[i] = missval;
    }
  else if (field.numMissVals)
    {
      for (size_t i = 0; i < len; ++i)
        if (fp_is_equal(array[i], missval) || !compare(array[i], c)) array[i] = missval;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (!compare(array[i], c)) array[i] = missval;
    }

  field_num_mv(field);
}

static int
le(double a, double b)
{
  return a <= b;
}

static int
lt(double a, double b)
{
  return a < b;
}

static int
ge(double a, double b)
{
  return a >= b;
}

static int
gt(double a, double b)
{
  return a > b;
}

static int
eq(double a, double b)
{
  return fp_is_equal(a, b);
}

static int
ne(double a, double b)
{
  return fp_is_not_equal(a, b);
}

void
vfarnum(Field &field1, Field const &field2)
{
  count(field1, field2, 0.0);
}

void
vfarnum2(Field &field1, Field const &field2)
{
  count(field1, field2, 1.0);
}

void
vfarnum3(Field &field1, Field const &field2, double n)
{
  count(field1, field2, n);
}

void
vfarsel(Field &field1, Field const &field2)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different gridsize (%s)", __func__);

  if (field2.numMissVals)
    {
      for (size_t i = 0; i < len; ++i)
        if (fp_is_equal(array2[i], missval2) || fp_is_equal(array2[i], 0.0)) array1[i] = missval1;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (is_equal(array2[i], 0.0)) array1[i] = missval1;
    }

  field_num_mv(field1);
}

void
vfarselle(Field &field1, Field const &field2)
{
  selcomp(field1, field2, le);
}

void
vfarsellt(Field &field1, Field const &field2)
{
  selcomp(field1, field2, lt);
}

void
vfarselge(Field &field1, Field const &field2)
{
  selcomp(field1, field2, ge);
}

void
vfarselgt(Field &field1, Field const &field2)
{
  selcomp(field1, field2, gt);
}

void
vfarseleq(Field &field1, Field const &field2)
{
  selcomp(field1, field2, eq);
}

void
vfarselne(Field &field1, Field const &field2)
{
  selcomp(field1, field2, ne);
}

void
vfarsellec(Field &field, double c)
{
  selcompc(field, c, le);
}

void
vfarselltc(Field &field, double c)
{
  selcompc(field, c, lt);
}

void
vfarselgec(Field &field, double c)
{
  selcompc(field, c, ge);
}

void
vfarseleqc(Field &field, double c)
{
  selcompc(field, c, eq);
}

void
vfarselnec(Field &field, double c)
{
  selcompc(field, c, ne);
}

void
vfarselgtc(Field &field, double c)
{
  selcompc(field, c, gt);
}

void
update_hist(FieldVector2D &field, int nlevels, size_t gridsize, std::vector<double> const &yvals, bool onlyNorth)
{
  for (int levelID = 0; levelID < nlevels; ++levelID)
    for (size_t i = 0; i < gridsize; ++i)
      if (onlyNorth)
        {
          if (yvals[i] >= 0.0) field[1][levelID].vec_d[i] = field[0][levelID].vec_d[i];
        }
      else
        field[1][levelID].vec_d[i] = field[0][levelID].vec_d[i];
}

void
define_mid_of_time(int frequency, int taxisID, int year, int month, int MaxMonths)
{
  CdiDateTime vDateTimeBound{};
  CdiDateTime vDateTimeBoundP1{};

  auto calendar = taxisInqCalendar(taxisID);

  if (frequency == 8)
    {
      vDateTimeBound.date = cdiDate_encode(year, month, 1);

      auto boundmonth = (month + 1 <= MaxMonths) ? month + 1 : 1;
      auto boundyear = (boundmonth != 1) ? year : year + 1;
      vDateTimeBoundP1.date = cdiDate_encode(boundyear, boundmonth, 1);
    }
  else
    {
      vDateTimeBound.date = cdiDate_encode(year, 1, 1);
      vDateTimeBoundP1.date = cdiDate_encode(year + 1, 1, 1);
    }

  auto julianDate1 = julianDate_encode(calendar, vDateTimeBound);
  auto julianDate2 = julianDate_encode(calendar, vDateTimeBoundP1);

  auto seconds = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / 2;
  auto julianDatem = julianDate_add_seconds(julianDate1, std::lround(seconds));
  auto vDateTime = julianDate_decode(calendar, julianDatem);

  taxisDefVdatetime(taxisID, vDateTime);
}

void
adjust_end_date(int nlevels, size_t gridsize, std::vector<double> const &yvals, double missval, int64_t ovdate,
                const FieldVector2D &startDateWithHist, FieldVector2D &endDateWithHist)
{
  auto ovdateSouth = std::min(cdiEncodeDate(ovdate / 10000, 6, 30), (int) ovdate);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      for (size_t i = 0; i < gridsize; ++i)
        {
          // start with southern sphere
          if (yvals[i] < 0)
            {
              if (fp_is_equal(startDateWithHist[1][levelID].vec_d[i], missval))
                {
                  endDateWithHist[0][levelID].vec_d[i] = missval;
                  continue;
                }
              if (fp_is_equal(endDateWithHist[0][levelID].vec_d[i], missval))
                {
                  endDateWithHist[0][levelID].vec_d[i] = ovdateSouth;
                }
            }
          else
            {
              if (fp_is_equal(startDateWithHist[0][levelID].vec_d[i], missval))
                {
                  endDateWithHist[0][levelID].vec_d[i] = missval;
                  continue;
                }

              if (fp_is_equal(endDateWithHist[0][levelID].vec_d[i], missval)) { endDateWithHist[0][levelID].vec_d[i] = ovdate; }
            }
        }
    }
}

void
compute_gsl(int nlevels, size_t gridsize, std::vector<double> &yvals, double missval, FieldVector2D &startDateWithHist,
            FieldVector2D &endDateWithHist, FieldVector &gslDuration, FieldVector &gslFirstDay, bool useCurrentYear)
{
  double firstDay, duration;

  if (!useCurrentYear)
    {
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          for (size_t i = 0; i < gridsize; ++i)
            {
              // start with southern sphere
              if (yvals[i] < 0.0)
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[0][levelID].vec_d[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[1][levelID].vec_d[i]));
                }
              else
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[1][levelID].vec_d[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[1][levelID].vec_d[i]));
                }

              if (fp_is_equal(startDateWithHist[1][levelID].vec_d[i], missval))
                firstDay = missval;
              else
                firstDay = (double) day_of_year(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[1][levelID].vec_d[i]);

              gslDuration[levelID].vec_d[i] = duration;
              gslFirstDay[levelID].vec_d[i] = firstDay;
            }
        }
    }
  else
    {
      // the current year can only have values for the northern hemisphere
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          for (size_t i = 0; i < gridsize; ++i)
            {
              // start with southern sphere
              if (yvals[i] < 0.0)
                {
                  gslDuration[levelID].vec_d[i] = missval;
                  gslFirstDay[levelID].vec_d[i] = missval;
                }
              else
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[0][levelID].vec_d[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[0][levelID].vec_d[i]));

                  if (fp_is_equal(startDateWithHist[0][levelID].vec_d[i], missval))
                    firstDay = missval;
                  else
                    firstDay = (double) day_of_year(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[0][levelID].vec_d[i]);

                  gslDuration[levelID].vec_d[i] = duration;
                  gslFirstDay[levelID].vec_d[i] = firstDay;
                }
            }
        }
    }

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      field_num_mv(gslDuration[levelID]);
      field_num_mv(gslFirstDay[levelID]);
    }
}

void
write_gsl_stream(CdoStreamID ostreamID, int otaxisID, int otsID, int ovarID1, int ovarID2, int ivlistID1, int first_var_id,
                 FieldVector &gslDuration, FieldVector &gslFirstDay, CdiDateTime const &vDateTime, int nlevels)
{
  (void) ivlistID1;
  (void) first_var_id;

  taxisDefVdatetime(otaxisID, vDateTime);
  cdo_def_timestep(ostreamID, otsID);

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      cdo_def_field(ostreamID, ovarID1, levelID);
      cdo_write_field(ostreamID, gslDuration[levelID]);
    }

  for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      cdo_def_field(ostreamID, ovarID2, levelID);
      cdo_write_field(ostreamID, gslFirstDay[levelID]);
    }
}
