/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PIPE_H
#define PIPE_H

#include <mutex>
#include <condition_variable>

#include "field.h"

struct pipe_t
{
public:
  pipe_t() { pipe_init(); }

  void pipe_init();

  int pipe_inq_vlist(int &vlistID);
  void pipe_def_vlist(int &target_vlistID, int new_vlistID);

  int pipe_inq_timestep(int p_tsID);
  void pipe_def_timestep(int p_vlistID, int tsID);

  int pipe_inq_field(int *varID, int *levelID);
  void pipe_def_field(int p_varId, int p_levelID);

  void pipe_write_field(const double *const p_data, size_t p_numMissVals);
  void pipe_write_field(const float *const p_data, size_t p_numMissVals);
  void pipe_write_field(const Field *const p_data, size_t p_numMissVals);

  size_t pipe_read_field(int p_vlistID, double *data, size_t *numMissVals);
  size_t pipe_read_field(int p_vlistID, float *data, size_t *numMissVals);
  size_t pipe_read_field(int p_vlistID, Field *data, size_t *numMissVals);

  size_t pipe_read_pipe_field(double *data, int vlistID, size_t *p_numMissVals);
  size_t pipe_read_pipe_field(float *data, int vlistID, size_t *p_numMissVals);

  void pipe_set_name(int processID, int inputIDX);
  void close();

  bool EOP;
  bool usedata;
  bool hasdata;

  int varID, levelID;
  int fieldIDr, fieldIDw, tsIDr, tsIDw;

  size_t numMissVals;
  int numFields;

  bool dataIsFloat;
  double *data_d;
  float *data_f;

  std::mutex m_mutex;
  std::condition_variable tsDef_cond, tsInq_cond, vlistDef_cond, isClosed_cond;
  std::condition_variable recDef_cond, recInq_cond;
  std::condition_variable write_cond, read_cond;

  std::string name;

private:
  void wait_for_read();
};

#endif /* PIPE_H */
