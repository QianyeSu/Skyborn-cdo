/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef PIPESTREAM_H
#define PIPESTREAM_H

#include "cdoStream.h"
#include "pipe.h"

#ifdef HAVE_LIBPTHREAD

class FileStream;  // Predeclaration only for copy_field(...)

class PipeStream : public CdoStream
{
public:
  // Constructors
  explicit PipeStream(int p_processID);
  // ---

  // CdoStream Interface functions
  int open_read() override;
  int open_write(int p_filetype) override;
  int open_append() override;
  int get_id() override;

  int inq_vlist() override;
  void def_vlist(int p_vlistID) override;

  void inq_field(int *varID, int *levelID) override;
  void def_field(int varID, int levelID) override;

  void read_field(float *const p_data, size_t *numMissVals) override;
  void read_field(double *const p_data, size_t *numMissVals) override;
  void read_field(Field *const p_field, size_t *numMissVals) override;

  void write_field(const float *const p_data, size_t numMissVals) override;
  void write_field(const double *const p_data, size_t numMissVals) override;
  void write_field(const Field *const p_field, size_t numMissVals) override;

  void copy_field(CdoStreamID p_fileStream) override;

  int inq_timestep(int tsID) override;
  void def_timestep(int tsID) override;

  int inqFileType() override;
  int inqByteorder() override;

  void close() override;

  size_t getNvals() override;
  // ---

  // FileStreamOnly
  // ---

private:
  PipeStream() = delete;
  std::shared_ptr<pipe_t> m_pipe = std::make_shared<pipe_t>();
  pthread_t rthreadID;  // read  thread ID
  pthread_t wthreadID;  // write thread ID
  void waitForPipe();
};
#endif

#endif
