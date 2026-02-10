/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cthread_debug.h"

#include <mutex>

static std::mutex streamOpenMutex;
static std::mutex streamMutex;

int
stream_open_read_locked(const char *const p_filename)
{
  open_lock();
  auto streamID = streamOpenRead(p_filename);
  open_unlock();
  if (streamID < 0) cdi_open_error(streamID, "Open failed on >%s<", p_filename);

  return streamID;
}

void
stream_close_locked(int p_fileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamClose(p_fileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_inq_field_locked(int p_fileID, int *const p_varID, int *const p_levelID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamInqField(p_fileID, p_varID, p_levelID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_def_field_locked(int p_fileID, int p_varID, int p_levelID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamDefField(p_fileID, p_varID, p_levelID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_read_field_float_locked(int p_fileID, float *const p_data, size_t *const p_numMissVals)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamReadFieldF(p_fileID, p_data, p_numMissVals);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_read_field_double_locked(int p_fileID, double *const p_data, size_t *const p_numMissVals)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamReadField(p_fileID, p_data, p_numMissVals);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_def_vlist_locked(int p_fileID, int p_vlistID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamDefVlist(p_fileID, p_vlistID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

int
stream_inq_vlist_locked(int p_fileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  auto vlistID = streamInqVlist(p_fileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);

  return vlistID;
}

void
stream_write_field_double_locked(int p_fileID, const double *const p_data, size_t p_numMissVals)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamWriteField(p_fileID, p_data, p_numMissVals);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

void
stream_write_field_float_locked(int p_fileID, const float *const p_data, size_t p_numMissVals)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamWriteFieldF(p_fileID, p_data, p_numMissVals);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
}

int
stream_inq_time_step_locked(int p_fileID, int p_tsID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  auto numFields = streamInqTimestep(p_fileID, p_tsID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);

  return numFields;
}

int
stream_def_time_step_locked(int p_fileID, int p_tsID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  auto success = streamDefTimestep(p_fileID, p_tsID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
  return success;
}

int
stream_copy_field_locked(int p_fileID, int p_targetFileID)
{
  if (Threading::cdoLockIO) cthread_mutex_lock(streamMutex);
  streamCopyField(p_fileID, p_targetFileID);
  if (Threading::cdoLockIO) cthread_mutex_unlock(streamMutex);
  return p_targetFileID;
}

void
vlist_copy_flag_locked(int p_vlistID2, int p_vlistID1)
{
  cthread_mutex_lock(streamMutex);
  vlistCopyFlag(p_vlistID2, p_vlistID1);
  cthread_mutex_unlock(streamMutex);
}

void
open_lock()
{
  cthread_mutex_lock(Threading::cdoLockIO ? streamMutex : streamOpenMutex);
}

void
open_unlock()
{
  cthread_mutex_unlock(Threading::cdoLockIO ? streamMutex : streamOpenMutex);
}

void
cdo_vlist_copy_flag(int vlistID2, int vlistID1)
{
  vlist_copy_flag_locked(vlistID2, vlistID1);
}
