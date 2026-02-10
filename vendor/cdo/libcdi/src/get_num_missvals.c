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

#include "get_num_missvals.h"
#include "cdi_int.h"

size_t
get_num_missvalsSP(size_t size, float *data, float missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
  {
    for (size_t i = 0; i < size; i++)
      if (DBL_IS_EQUAL(data[i], missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < size; i++)
      if (IS_EQUAL(data[i], missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}

size_t
get_num_missvalsDP(size_t size, double *data, double missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
  {
    for (size_t i = 0; i < size; i++)
      if (DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float) missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < size; i++)
      if (IS_EQUAL(data[i], missval) || IS_EQUAL(data[i], (float) missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}

size_t
get_cplx_num_missvalsSP(size_t size, float *data, float missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
  {
    for (size_t i = 0; i < 2 * size; i += 2)
      if (DBL_IS_EQUAL(data[i], missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < 2 * size; i += 2)
      if (IS_EQUAL(data[i], missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}

size_t
get_cplx_num_missvalsDP(size_t size, double *data, double missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
  {
    for (size_t i = 0; i < 2 * size; i += 2)
      if (DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float) missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < 2 * size; i += 2)
      if (IS_EQUAL(data[i], missval) || IS_EQUAL(data[i], (float) missval))
      {
        data[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}
