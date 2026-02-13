/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include "process_int.h"
#include "par_io.h"

void *
read_field(void *arg)
{
  auto *read_arg = (read_arg_t *) arg;
  auto streamID = read_arg->streamID;
  auto p_varID = read_arg->varID;
  auto p_levelID = read_arg->levelID;
  auto p_numMissVals = read_arg->numMissVals;
  auto p_array = read_arg->array;

  // std::fprintf(stderr, "read_field: streamID = %d\n", streamID);
  auto [varID, levelID] = cdo_inq_field(streamID);
  *p_varID = varID;
  *p_levelID = levelID;
  cdo_read_field(streamID, p_array, p_numMissVals);
  // std::fprintf(stderr, "read_field: varID %d levelID %d\n", *p_varID, *p_levelID);

  return nullptr;
}

void
par_read_field(CdoStreamID streamID, int *varID, int *levelID, double *array, size_t *numMissVals, par_io_t *parIO)
{
  bool lpario = false;
  int fieldID = 0, numFields = 0;
#ifdef HAVE_LIBPTHREAD
  pthread_t thrID = 0;
  // pthread_attr_t attr;
  int rval;

  if (parIO)
  {
    lpario = true;
    fieldID = parIO->fieldID;
    numFields = parIO->numFields;
    thrID = parIO->thrID;
  }
#endif

  if (fieldID == 0 || !lpario)
  {
    read_arg_t read_arg;
    read_arg.streamID = streamID;
    read_arg.varID = varID;
    read_arg.levelID = levelID;
    read_arg.numMissVals = numMissVals;
    read_arg.array = array;

    read_field(&read_arg);
  }
#ifdef HAVE_LIBPTHREAD
  else
  {
    // std::fprintf(stderr, "parIO1: %ld streamID %d %d %d\n", (long)thrID, streamID, fieldID, numFields);
    rval = pthread_join(thrID, nullptr);
    if (rval != 0) cdo_abort("pthread_join failed!");

    *varID = parIO->varID;
    *levelID = parIO->levelID;
    *numMissVals = parIO->numMissVals;
    // std::fprintf(stderr, "parIO2: %ld streamID %d %d %d\n", (long)thrID, streamID, *varID, *levelID);
    array_copy(parIO->array_size, parIO->array, array);
  }

  if (lpario && numFields > 1)
  {
    read_arg_t *read_arg = &(parIO->read_arg);
    if ((fieldID + 1) < numFields)
    {
      if (fieldID == 0)
      {
        pthread_attr_init(&parIO->attr);
        pthread_attr_setdetachstate(&parIO->attr, PTHREAD_CREATE_JOINABLE);
      }

      read_arg->streamID = streamID;
      read_arg->varID = &parIO->varID;
      read_arg->levelID = &parIO->levelID;
      read_arg->numMissVals = &parIO->numMissVals;
      read_arg->array = parIO->array;

      // std::fprintf(stderr, "pthread_create: streamID %d %d\n", read_arg->streamID,streamID);
      rval = pthread_create(&thrID, &parIO->attr, read_field, read_arg);
      if (rval != 0) cdo_abort("pthread_create failed!");

      // std::fprintf(stderr, "thrID = %ld\n", (long) thrID);
      parIO->thrID = thrID;
    }
    else
      pthread_attr_destroy(&parIO->attr);
  }
#endif
}
