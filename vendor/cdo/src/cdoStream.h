#ifndef CDOSTREAM_H
#define CDOSTREAM_H

#include <vector>
#include <memory>

#include "datarangelist.h"
#include "field.h"

class FileStream;

class CdoStream;

using CdiStreamID = int;
using CdoStreamID = std::shared_ptr<CdoStream>;
#define CDO_STREAM_UNDEF nullptr

class CdoStream
{
public:
  // Constructors
  virtual int open_read() = 0;
  virtual int open_write(int p_filetype) = 0;
  virtual int open_append() = 0;

  virtual int inq_vlist() = 0;
  virtual void def_vlist(int p_vlistID) = 0;

  virtual void inq_field(int *varID, int *levelID) = 0;
  virtual void def_field(int varID, int levelID) = 0;

  virtual void read_field(float *const p_data, size_t *numMissVals) = 0;
  virtual void read_field(double *const p_data, size_t *numMissVals) = 0;
  virtual void read_field(Field *const p_field, size_t *numMissVals) = 0;

  virtual void write_field(const float *const p_data, size_t numMissVals) = 0;
  virtual void write_field(const double *const p_data, size_t numMissVals) = 0;
  virtual void write_field(const Field *const p_field, size_t numMissVals) = 0;

  virtual void copy_field(CdoStreamID dest) = 0;

  virtual int inq_timestep(int tsID) = 0;
  virtual void def_timestep(int tsID) = 0;

  virtual int inqFileType() = 0;
  virtual int inqByteorder() = 0;

  virtual void close() = 0;

  virtual size_t getNvals() = 0;

  virtual int get_id() = 0;
  int getTsID();

  int m_cdoStreamID = -1;  // aka the id of the pstream
  int m_filetype = CDI_UNDEFID;
  size_t m_nvals = 0;
  bool isopen = false;

  std::string m_name;
  std::vector<Datarange> m_datarangelist;

  int m_vlistID = -1;
  int m_tsID = -1;
  int m_varID = -1;  // next varID defined with streamDefVar

  // to be removed or to be moved to FileStream! // some operators need some refactoring for these to be able to be moved
  int m_fileID = 0;
  bool ispipe = false;

protected:
  CdoStream();
  virtual ~CdoStream() = 0;

private:
};

#endif /* CDOSTREAM_H */
