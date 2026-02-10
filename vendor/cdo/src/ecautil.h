/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/
#ifndef ECAUTIL_H_
#define ECAUTIL_H_

#include "field.h"
#include "cdoStream.h"

/**
 * Computes the day-of-year correspnding a given Gregorian date.
 *
 * @param date a Gregorian date in the form YYYYMMDD
 *
 * @return the day-of-year
 */
unsigned long day_of_year(int date);

/**
 * Counts the number of nonumMissValsing values in a field. The result
 * of the operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    0
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */
void vfarnum(Field &field1, Field const &field2);

/**
 * Counts the number of consecutive nonumMissValsing values in a field.
 * The result of the operation is computed according to the following
 * rules:
 *
 * field1  field2  result
 * a       b       a + 1
 * a       miss    0
 * miss    b       1
 * miss    miss    0
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */
void vfarnum2(Field &field1, Field const &field2);

/**
 * Counts the number of values in series of at least n consecutive
 * nonumMissValsing values. The result of the operation is computed according
 * to the following rules:
 *
 * field1  field2  result
 * a       b       b < n ? a : b > n ? a + 1 : a + n
 * a       miss    a
 * miss    b       b < n ? 0 : b
 * miss    miss    0
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 * @param n      the number of consecutive values, must be an exact
 *               mathematical integer
 */
void vfarnum3(Field &field1, Field const &field2, double n);

/**
 * Selects field elements according to a given mask. The result of
 * the operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       b != 0.0 ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1  the input field, also holds the result
 * @param field2  the mask
 */
void vfarsel(Field &field1, Field const &field2);

/**
 * Selects all field elements that are less than or equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a <= b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarselle(Field &field1, Field const &field2);

/**
 * Selects all field elements that are less than the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a < b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarsellt(Field &field1, Field const &field2);

/**
 * Selects all field elements that are greater than or equal to
 * the corresponding element of a reference field. The result of
 * the operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a >= b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarselge(Field &field1, Field const &field2);

/**
 * Selects all field elements that are greater than the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a > b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarselgt(Field &field1, Field const &field2);

/**
 * Selects all field elements that are equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a == b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarseleq(Field &field1, Field const &field2);

/**
 * Selects all field elements that are not equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a != b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */
void vfarselne(Field &field1, Field const &field2);

/**
 * Selects all field elements that are less than or equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a <= c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarsellec(Field &field, double c);

/**
 * Selects all field elements that are less a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a < c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarselltc(Field &field, double c);

/**
 * Selects all field elements that are greater than or equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a >= c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarselgec(Field &field, double c);

/**
 * Selects all field elements that are greater than a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a > c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarselgtc(Field &field, double c);

/**
 * Selects all field elements that are equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a == c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarseleqc(Field &field, double c);

/**
 * Selects all field elements that are not equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 *
 * field  c      result
 * a      c      a != c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field the input field, also holds the result
 * @param c     the reference value
 */
void vfarselnec(Field &field, double c);

/**
 * reset the fields real values to the missval for all levels
 *
 * @param field     list of fields: 0 is index of the current values, 1 hold
 *                  the values of the previous year
 * @param nlevels   number of available levels
 * @param gridsize  number of grid points
 * @param yvals     list of latitudes
 * @param onlyNorth boolean for processing only the norther hemisphere
 */
void update_hist(FieldVector2D &field, int nlevels, size_t gridsize, std::vector<double> const &yvals, bool onlyNorth);

/*
 * Compute the Gsl and its starting day
 *
 * @param int nlevels
 * @param size_t gridsize
 * @param double *yvals = array of latitudes
 * @param int ysize = number of gridpoints in lat-direction
 * @param double missval
 * @param int ovdate = the last year, which has been fully processed
 * @param Field *startDate
 * @param Field *endDate
 * @param Field *startDateWithHist[2]
 * @param Field *endDateWithHist[2]
 * @param Field *gslDuration
 * @param Field *gslFirstDay
 * @param bool useCurrentYear = if true, only measurements of the current year
 *                             (index 0) are used for computation, i.e. that
 *                             gsl can only be computed for the northern
 *                             hemisphere (see definition of GSL: EcaGsl()
 */
void compute_gsl(int nlevels, size_t gridsize, std::vector<double> &yvals, double missval, FieldVector2D &startDateWithHist,
                 FieldVector2D &endDateWithHist, FieldVector &gslDuration, FieldVector &gslFirstDay, bool useCurrentYear);

/*
 * Adjust the endDates found in the current year:
 * if a starting date for gsl could be found, but no ending date, the end
 * should be the last day of the corresponding year for norther and June, 30th
 * for southern hemisphere
 */
void adjust_end_date(int nlevels, size_t gridsize, std::vector<double> const &yvals, double missval, int64_t ovdate,
                     const FieldVector2D &startDateWithHist, FieldVector2D &endDateWithHist);
/*
 * Calculates the mid of the year or month
 */
void define_mid_of_time(int frequency, int taxisID, int year, int month, int MaxMonths);
/*
 * Write GSL related fields to an output stream
 */
void write_gsl_stream(CdoStreamID ostreamID, int otaxisID, int otsID, int ovarID1, int ovarID2, int ivlistID1, int first_var_id,
                      FieldVector &gslDuration, FieldVector &gslFirstDay, CdiDateTime const &vDateTime, int nlevels);
#endif /*ECAUTIL_H_*/
