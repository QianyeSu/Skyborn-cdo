/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#include "pipeStream.h"
#include "cdo_output.h"

#ifdef HAVE_LIBPTHREAD

#include "pthread_debug.h"

PipeStream::PipeStream(int p_processID)
{
  ispipe = true;
  m_pipe->pipe_set_name(p_processID, m_cdoStreamID);
  m_name = m_pipe->name;
}

int PipeStream::get_id(){
  return m_cdoStreamID;
}

int
PipeStream::open_read()
{
  rthreadID = pthread_self();
  isopen = true;
  return m_cdoStreamID;
}

int
PipeStream::open_write(int filetype)
{
  Debug(PIPE_STREAM, "pipe %s", m_pipe->name);

  wthreadID = pthread_self();
  m_filetype = filetype;
  isopen = true;

  return m_cdoStreamID;
}

int
PipeStream::open_append()
{
  cdo_warning("Operator does not support pipes!");
  return -1;
}

int
PipeStream::inq_vlist()
{
  // m_vlist is changed when the whlist was successfully defined by another pipe!!
  auto vlistID = m_pipe->pipe_inq_vlist(m_vlistID);
  if (vlistID == -1) cdo_abort("Couldn't read data from input stream %s!", m_name);
  return vlistID;
}

void
PipeStream::def_vlist(int p_vlistID)
{
  Debug(PIPE_STREAM, "%s pstreamID %d", m_pipe->name, m_cdoStreamID);
  auto vlistIDcp = vlistDuplicate(p_vlistID);
  m_pipe->pipe_def_vlist(m_vlistID, vlistIDcp);
}

void
PipeStream::inq_field(int *varID, int *levelID)
{
  m_pipe->pipe_inq_field(varID, levelID);
}

void
PipeStream::def_field(int varID, int levelID)
{
  m_pipe->pipe_def_field(varID, levelID);
}

void
PipeStream::read_field(float *p_data, size_t *p_numMissVals)
{
  m_nvals += m_pipe->pipe_read_field(m_vlistID, p_data, p_numMissVals);
}

void
PipeStream::read_field(double *p_data, size_t *p_numMissVals)
{
  m_nvals += m_pipe->pipe_read_field(m_vlistID, p_data, p_numMissVals);
}

void
PipeStream::read_field(Field *p_field, size_t *p_numMissVals)
{
  m_nvals += m_pipe->pipe_read_field(m_vlistID, p_field, p_numMissVals);
}

void
PipeStream::write_field(const float *p_data, size_t p_numMissVals)
{
  m_pipe->pipe_write_field(p_data, p_numMissVals);
}

void
PipeStream::write_field(const double *p_data, size_t p_numMissVals)
{
  m_pipe->pipe_write_field(p_data, p_numMissVals);
}

void
PipeStream::write_field(const Field *p_field, size_t p_numMissVals)
{
  m_pipe->pipe_write_field(p_field, p_numMissVals);
}

void
PipeStream::copy_field(CdoStreamID p_destination)
{
  (void) p_destination;
  // Not implemented for pipes
  // Cdi handles this. And in cdo we would have to decompress and recompress for copy operations
  // which is very resource intensive (also lossy compression)
  cdo_warning("Copy field not possible with piped streams!");
}

int
PipeStream::inq_timestep(int p_tsID)
{
  auto numFields = m_pipe->pipe_inq_timestep(p_tsID);
  m_tsID = p_tsID;
  Debug(PIPE_STREAM, "PipeStream: Current TsID: %d,  numFields: %d", m_tsID, numFields);
  return numFields;
}

void
PipeStream::def_timestep(int p_tsID)
{
  Debug(PIPE_STREAM, "%s pstreamID %d", m_pipe->name, m_cdoStreamID);
  m_pipe->pipe_def_timestep(m_vlistID, p_tsID);
}

int
PipeStream::inqFileType()
{
  return m_filetype;
}

int
PipeStream::inqByteorder()
{
  return m_filetype;
}

void
PipeStream::waitForPipe()
{
  m_pipe->close();
  std::unique_lock<std::mutex> lock(m_pipe->m_mutex);
  while (isopen)
    {
      Debug(PIPE_STREAM, "wait of read close");
      m_pipe->isClosed_cond.wait(lock);
    }
}

void
PipeStream::close()
{
  auto threadID = pthread_self();

  Debug(PIPE_STREAM, "streamID: %d thID: %ld rthID: %ld wthID: %ld", get_id(), threadID, rthreadID, wthreadID);

  if (pthread_equal(threadID, rthreadID))
    {
      isopen = false;
      m_pipe->close();
      pthread_join(wthreadID, nullptr);
    }
  else if (pthread_equal(threadID, wthreadID)) { waitForPipe(); }
  else { cdo_abort("Internal problem! Close pipe %s", m_name); }
}

size_t
PipeStream::getNvals()
{
  return m_nvals;
}

#endif
