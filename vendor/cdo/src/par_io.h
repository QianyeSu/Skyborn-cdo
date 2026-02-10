#ifndef PAR_IO_H
#define PAR_IO_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include <cstddef>
#include "cdoStream.h"

struct read_arg_t
{
  CdoStreamID streamID;
  int *varID, *levelID;
  size_t *numMissVals;
  double *array;
};

struct par_io_t
{
  int varID, levelID;
  size_t numMissVals;
  double *array;
  int array_size;
  int fieldID;
  int numFields;
  read_arg_t read_arg;
#ifdef HAVE_LIBPTHREAD
  pthread_t thrID;
  pthread_attr_t attr;
#endif
};

void par_read_field(CdoStreamID streamID, int *varID, int *levelID, double *array, size_t *numMissVals, par_io_t *parIO);

#endif /* PAR_IO_H */
