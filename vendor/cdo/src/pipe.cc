/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#include <cdi.h>

#include "pipe.h"
#include "cdo_output.h"
#include "cthread_debug.h"

void
pipe_t::close()
{
  {
    std::scoped_lock lock(m_mutex);
    EOP = true;
    Debug(PIPE, "%s write closed", name);
  }

  cthread_cond_signal(tsDef_cond);
  cthread_cond_signal(tsInq_cond);
  cthread_cond_signal(recInq_cond);
  cthread_cond_signal(isClosed_cond);
}

void
pipe_t::pipe_init()
{
  EOP = false;

  fieldIDr = -1;
  fieldIDw = -1;
  tsIDr = -1;
  tsIDw = -1;

  varID = -1;
  levelID = -1;

  numFields = 0;
  numMissVals = 0;
  data_d = nullptr;
  data_f = nullptr;
  hasdata = false;
  usedata = true;
  dataIsFloat = false;
}

int
pipe_t::pipe_inq_timestep(int p_tsID)
{
  int numrecs = 0;

  {
    std::unique_lock<std::mutex> lock(m_mutex);
    usedata = false;
    fieldIDr = -1;
    if (p_tsID != tsIDr + 1)
    {
      if (!(p_tsID == tsIDr && tsIDr == tsIDw && fieldIDr == -1))
        cdo_abort("%s unexpected tsID %d %d %d", name, p_tsID, tsIDr + 1, tsIDw);
    }

    tsIDr = p_tsID;
    while (tsIDw != p_tsID)
    {
      if (EOP)
      {
        Debug(PIPE, "%s EOP", name.c_str());
        break;
      }
      if (hasdata)
      {
        Debug(PIPE, "%s has data", name.c_str());
        hasdata = false;
        data_d = nullptr;
        data_f = nullptr;
        dataIsFloat = false;
        read_cond.notify_all();
      }
      Debug(PIPE && !hasdata, "%s has no data", name);

      recInq_cond.notify_all(); /* o.k. ??? */

      Debug(PIPE, "%s wait for tsDef_cond", name.c_str());
      tsDef_cond.wait(lock);
    }

    numrecs = EOP ? 0 : numFields;
  }

  tsInq_cond.notify_all();

  return numrecs;
}

void
pipe_t::pipe_def_vlist(int &target_vlistID, int new_vlistID)
{
  {
    std::scoped_lock lock(m_mutex);
    target_vlistID = new_vlistID;
  }

  // lets the program know that the vlist is now defined
  vlistDef_cond.notify_all();
}

int
pipe_t::pipe_inq_vlist(int &p_vlistID)
{
  constexpr std::chrono::milliseconds timeOut(1000);
  constexpr int maxWaitCycles = 3600;

  Debug(PIPE, "Inquiring vlist for vlistID: %d", p_vlistID);
  std::chrono::milliseconds timeToWait(0);
  int nwaitcycles = 0;

  std::unique_lock<std::mutex> lock(m_mutex);
  while (p_vlistID == -1 && nwaitcycles < maxWaitCycles && !EOP)
  {
    timeToWait += timeOut;
    Debug(PIPE, "%s wait for vlistDef_cond", name.c_str());
    vlistDef_cond.wait_for(lock, timeToWait);
    nwaitcycles++;
  }

  return p_vlistID;
}

void
pipe_t::pipe_def_timestep(int p_vlistID, int p_tsID)
{
  {
    std::scoped_lock lock(m_mutex);
    fieldIDw = -1;
    tsIDw++;
    if (p_tsID != tsIDw) cdo_abort("unexpected p_tsID %d(%d) for %s", p_tsID, tsIDw, name);

    int numrecs = 0;
    if (p_tsID == 0) { numrecs = vlistNumFields(p_vlistID); }
    else
    {
      auto vlistID = p_vlistID;
      for (int i = 0; i < vlistNvars(vlistID); ++i)
      {
        if (vlistInqVarTimetype(vlistID, i) != TIME_CONSTANT) numrecs += zaxisInqSize(vlistInqVarZaxis(vlistID, i));
      }
      Debug(PIPE, " %s numrecs= %d nvars= %d ", name.c_str(), numrecs, vlistNvars(vlistID));
    }

    numFields = numrecs;
    Debug(PIPE, "%s numrecs %d p_tsID %d %d %d", name.c_str(), numrecs, p_tsID, tsIDw, tsIDr);
    if (numrecs == 0) EOP = true;
  }

  tsDef_cond.notify_all();
  // sleep(1);

  std::unique_lock<std::mutex> lock(m_mutex);
  while (tsIDr < p_tsID)
  {
    if (EOP)
    {
      Debug(PIPE, "EOP");
      break;
    }
    Debug(PIPE, " wait for tsInq_cond (p_tsID %d %d)", name.c_str(), p_tsID, tsIDr);
    tsInq_cond.wait(lock);
  }
}

int
pipe_t::pipe_inq_field(int *p_varID, int *p_levelID)
{
  bool condSignal = false;

  {
    std::scoped_lock lock(m_mutex);
    Debug(PIPE, "%s has no data %d %d ", name.c_str(), fieldIDr, " ", fieldIDw);
    if (hasdata || usedata)
    {
      hasdata = false;
      data_d = nullptr;
      data_f = nullptr;
      usedata = false;
      condSignal = true;
    }
  }

  if (condSignal) read_cond.notify_all();

  {
    std::unique_lock<std::mutex> lock(m_mutex);
    usedata = true;
    fieldIDr++;

    Debug(PIPE, "%s fieldID: r=%d w=%d", name.c_str(), fieldIDr, fieldIDw);

    while (fieldIDw != fieldIDr)
    {
      if (EOP)
      {
        Debug(PIPE, "EOP");
        break;
      }
      Debug(PIPE, "%s wait for recDef_cond", name);
      recDef_cond.wait(lock);
    }

    if (EOP)
    {
      *p_varID = -1;
      *p_levelID = -1;
    }
    else
    {
      *p_varID = varID;
      *p_levelID = levelID;
    }
  }

  recInq_cond.notify_all();

  return 0;
}

void
pipe_t::pipe_def_field(int p_varID, int p_levelID)
{
  bool condSignal = false;

  {
    std::scoped_lock lock(m_mutex);
    Debug(PIPE, "%s has data %d", name.c_str(), fieldIDr, fieldIDw);  //<- TODO: rethink positioning
    if (hasdata)
    {
      hasdata = false;
      data_d = nullptr;
      data_f = nullptr;
      condSignal = true;
    }
  }

  if (condSignal) read_cond.notify_all();

  {
    std::scoped_lock lock(m_mutex);
    usedata = true;
    fieldIDw++;
    varID = p_varID;
    levelID = p_levelID;
    Debug(PIPE, "%s fieldIDs: r=%d w=%d", name.c_str(), fieldIDr, fieldIDw);
  }

  recDef_cond.notify_all();

  std::unique_lock<std::mutex> lock(m_mutex);
  while (fieldIDr < fieldIDw)
  {
    if (tsIDw != tsIDr) break;
    if (EOP) break;
    Debug(PIPE, "%s wait for recInq_cond %d", name.c_str(), fieldIDr);
    recInq_cond.wait(lock);
  }
}

/***
 * copys data from a pipe to data
 *
 * @param data destination for the field data
 * @param pipe pipe that has the wanted data
 */
size_t
pipe_t::pipe_read_pipe_field(double *const p_data, int vlistID, size_t *const p_numMissVals)
{
  if (!p_data) cdo_abort("No data pointer for %s", name);

  auto datasize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  if (vlistNumber(vlistID) != CDI_REAL) datasize *= 2;

  if (dataIsFloat)
  {
    for (size_t i = 0; i < datasize; ++i) p_data[i] = (double) data_f[i];
  }
  else { std::memcpy(p_data, data_d, datasize * sizeof(double)); }

  *p_numMissVals = numMissVals;
  return datasize;
}

size_t
pipe_t::pipe_read_pipe_field(float *const p_data, int vlistID, size_t *const p_numMissVals)
{
  if (!p_data) cdo_abort("No data pointer for %s", name);

  auto datasize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  if (vlistNumber(vlistID) != CDI_REAL) datasize *= 2;

  if (dataIsFloat) { std::memcpy(p_data, data_f, datasize * sizeof(float)); }
  else
  {
    for (size_t i = 0; i < datasize; ++i) p_data[i] = (float) data_d[i];
  }

  *p_numMissVals = numMissVals;
  return datasize;
}

size_t
pipe_t::pipe_read_field(int p_vlistID, double *const p_data, size_t *const p_numMissVals)
{
  *p_numMissVals = 0;
  size_t nvals = 0;

  {
    std::unique_lock<std::mutex> lock(m_mutex);
    while (!hasdata)
    {
      Debug(PIPE, "%s wait for write_cond", name.c_str());
      write_cond.wait(lock);
    }

    if (hasdata) { nvals = pipe_read_pipe_field(p_data, p_vlistID, p_numMissVals); }
    else { cdo_abort("data type %d not implemented", hasdata); }

    Debug(PIPE, "%s read field %d", name.c_str(), fieldIDr);

    hasdata = false;
    data_d = nullptr;
  }

  read_cond.notify_all();

  return nvals;
}

size_t
pipe_t::pipe_read_field(int p_vlistID, float *const p_data, size_t *const p_numMissVals)
{
  *p_numMissVals = 0;
  size_t nvals = 0;

  {
    std::unique_lock<std::mutex> lock(m_mutex);
    while (!hasdata)
    {
      Debug(PIPE, "%s wait for write_cond", name.c_str());
      write_cond.wait(lock);
    }

    if (hasdata) { nvals = pipe_read_pipe_field(p_data, p_vlistID, p_numMissVals); }
    else { cdo_abort("data type %d not implemented", hasdata); }

    Debug(PIPE, "%s read field %d", name.c_str(), fieldIDr);

    hasdata = false;
    data_f = nullptr;
  }

  read_cond.notify_all();

  return nvals;
}

size_t
pipe_t::pipe_read_field(int p_vlistID, Field *const p_field, size_t *const p_numMissVals)
{
  return pipe_read_field(p_vlistID, p_field->vec_d.data(), p_numMissVals);
}

void
pipe_t::wait_for_read()
{
  write_cond.notify_all();

  Debug(PIPE, "%s write field %d", name, fieldIDw);

  std::unique_lock<std::mutex> lock(m_mutex);
  while (hasdata)
  {
    if (!usedata) break;
    if (fieldIDw != fieldIDr) break;

    if (EOP)
    {
      Debug(PIPE, "EOP");
      break;
    }
    Debug(PIPE, " wait for read_cond %s", name);
    read_cond.wait(lock);
  }
}

void
pipe_t::pipe_write_field(const double *const p_data, size_t p_numMissVals)
{
  {
    std::scoped_lock lock(m_mutex);
    hasdata = true;  // data pointer
    dataIsFloat = false;
    data_d = (double *) p_data;
    numMissVals = p_numMissVals;
  }

  wait_for_read();
}

void
pipe_t::pipe_write_field(const float *const p_data, size_t p_numMissVals)
{
  {
    std::scoped_lock lock(m_mutex);
    hasdata = true;  // data pointer
    dataIsFloat = true;
    data_f = (float *) p_data;
    numMissVals = p_numMissVals;
  }

  wait_for_read();
}

void
pipe_t::pipe_write_field(const Field *const p_field, size_t p_numMissVals)
{
  pipe_write_field(p_field->vec_d.data(), p_numMissVals);
}

void
pipe_t::pipe_set_name(int processID, int inputIDX)
{
  name = "(pipe" + std::to_string(processID + 1) + "." + std::to_string(inputIDX) + ")";
}
