#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "dmemory.h"

// the single image implementation
static int
cdiStreamReadVar(int streamID, int varID, int memType, void *data, size_t *numMissVals)
{
  // May fail if memType == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(numMissVals);

  stream_t *streamPtr = stream_to_pointer(streamID);
  const int fileType = streamPtr->filetype;

  *numMissVals = 0;

  if (memType == MEMTYPE_FLOAT && cdiFiletypeIsExse(fileType)) return 1;

  switch (cdiBaseFiletype(fileType))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_var(streamPtr, varID, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvReadVarDP(streamPtr, varID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extReadVarDP(streamPtr, varID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegReadVarDP(streamPtr, varID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_var(streamPtr, varID, memType, data, numMissVals); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(fileType));
  }

  return status;
}

/*
@Function  streamReadVar
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, double *data, SizeType *numMissVals)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  numMissVals     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVar(int streamID, int varID, double *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  cdiStreamReadVar(streamID, varID, MEMTYPE_DOUBLE, data, &numMiss);
  *numMissVals = (SizeType) numMiss;
}

/*
@Function  streamReadVarF
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, float *data, SizeType *numMissVals)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  numMissVals     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarF(int streamID, int varID, float *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  if (cdiStreamReadVar(streamID, varID, MEMTYPE_FLOAT, data, &numMiss))
  {
    // In case the file format does not support single precision reading,
    // we fall back to double precision reading, converting the data on the fly.
    size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
    elementCount *= (size_t) zaxisInqSize(vlistInqVarZaxis(streamInqVlist(streamID), varID));
    double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
    streamReadVar(streamID, varID, conversionBuffer, numMissVals);
    for (size_t i = elementCount; i--;) data[i] = (float) conversionBuffer[i];
    Free(conversionBuffer);
  }
  *numMissVals = (SizeType) numMiss;
}

static int
cdiStreamReadVarSlice(int streamID, int varID, int levelID, int memType, void *data, size_t *numMissVals)
{
  // May fail if memType == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(numMissVals);

  stream_t *streamPtr = stream_to_pointer(streamID);
  const int fileType = streamPtr->filetype;

  *numMissVals = 0;

  if (memType == MEMTYPE_FLOAT && cdiFiletypeIsExse(fileType)) return 1;

  switch (cdiBaseFiletype(fileType))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_var_slice(streamPtr, varID, levelID, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvReadVarSliceDP(streamPtr, varID, levelID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extReadVarSliceDP(streamPtr, varID, levelID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegReadVarSliceDP(streamPtr, varID, levelID, (double *) data, numMissVals); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_var_slice(streamPtr, varID, levelID, memType, data, numMissVals); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(fileType));
  }

  return status;
}

/*
@Function  streamReadVarSlice
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSlice(int streamID, int varID, int levelID, double *data, SizeType *numMissVals)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  numMissVals     Number of missing values.

@Description
The function streamReadVarSlice reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarSlice(int streamID, int varID, int levelID, double *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  if (cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, data, &numMiss))
  {
    Warning("Unexpected error returned from cdiStreamReadVarSlice()!");
    size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
    memset(data, 0, elementCount * sizeof(*data));
  }
  *numMissVals = (SizeType) numMiss;
}

/*
@Function  streamReadVarSliceF
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, SizeType *numMissVals)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  numMissVals     Number of missing values.

@Description
The function streamReadVarSliceF reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarSliceF(int streamID, int varID, int levelID, float *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  if (cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, data, &numMiss))
  {
    // In case the file format does not support single precision reading,
    // we fall back to double precision reading, converting the data on the fly.
    size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
    double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
    streamReadVarSlice(streamID, varID, levelID, conversionBuffer, numMissVals);
    for (size_t i = elementCount; i--;) data[i] = (float) conversionBuffer[i];
    Free(conversionBuffer);
  }
  *numMissVals = (SizeType) numMiss;
}

static void
stream_read_record(int streamID, int memType, void *data, size_t *numMissVals)
{
  check_parg(data);
  check_parg(numMissVals);

  stream_t *streamPtr = stream_to_pointer(streamID);

  if (streamPtr->lockIO) CDI_IO_LOCK();

  *numMissVals = 0;

  switch (cdiBaseFiletype(streamPtr->filetype))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_field(streamPtr, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srv_read_field(streamPtr, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: ext_read_field(streamPtr, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: ieg_read_field(streamPtr, memType, data, numMissVals); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_field(streamPtr, memType, data, numMissVals); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(streamPtr->filetype));
  }

  if (streamPtr->lockIO) CDI_IO_UNLOCK();
}

void
streamReadField(int streamID, double *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  stream_read_record(streamID, MEMTYPE_DOUBLE, (void *) data, &numMiss);
  *numMissVals = (SizeType) numMiss;
}

void
streamReadFieldF(int streamID, float *data, SizeType *numMissVals)
{
  size_t numMiss = 0;
  stream_read_record(streamID, MEMTYPE_FLOAT, (void *) data, &numMiss);
  *numMissVals = (SizeType) numMiss;
}

static int
cdiStreamReadVarSlicePart(int streamID, int varID, int levelID, int varType, int start, size_t size, int memType, void *data,
                          size_t *numMissVals)
{
  int status = 0;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(numMissVals);

  stream_t *streamPtr = stream_to_pointer(streamID);
  int fileType = streamPtr->filetype;

  *numMissVals = 0;

  // currently we only care for netcdf data
  switch (cdiBaseFiletype(fileType))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB:
    {
      grb_read_var_slice(streamPtr, varID, levelID, memType, data, numMissVals);
      break;
    }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF:
    {
      cdf_read_var_slice_part(streamPtr, varID, levelID, memType, varType, start, size, data, numMissVals);
      break;
    }
#endif
    default:
    {
      Error("%s support not compiled in!", strfiletype(fileType));
      status = 2;
      break;
    }
  }

  return status;
}

static void
cdiStreamReadVarPart(int streamID, int varID, int varType, int start, size_t size, int memType, void *data, size_t *numMissVals)
{
  (void) (varType);
  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(numMissVals);

  stream_t *streamPtr = stream_to_pointer(streamID);
  int fileType = streamPtr->filetype;

  *numMissVals = 0;

  // currently we only care for netcdf data
  switch (cdiBaseFiletype(fileType))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB:
    {
      grb_read_var(streamPtr, varID, memType, data, numMissVals);
      break;
    }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF:
    {
      cdf_read_var_part(streamPtr, varID, memType, varType, start, size, data, numMissVals);
      break;
    }
#endif
    default:
    {
      Error("%s support not compiled in!", strfiletype(fileType));
      break;
    }
  }
}

void
streamReadVarSlicePart(int streamID, int varID, int levelID, int varType, int start, SizeType size, void *data,
                       SizeType *numMissVals, int memMype)
{
  size_t numMiss = 0;
  if (cdiStreamReadVarSlicePart(streamID, varID, levelID, varType, start, (size_t) size, memMype, data, &numMiss))
  {
    Error("Unexpected error returned from cdiStreamReadVarSlicePart()!");
  }
  *numMissVals = (SizeType) numMiss;
}

void
streamReadVarPart(int streamID, int varID, int varType, int start, SizeType size, void *data, SizeType *numMissVals, int memType)
{
  size_t numMiss = 0;
  cdiStreamReadVarPart(streamID, varID, varType, start, (size_t) size, memType, data, &numMiss);
  *numMissVals = (SizeType) numMiss;
}
