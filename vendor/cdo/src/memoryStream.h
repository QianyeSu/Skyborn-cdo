#ifndef MEMORY_STREAM_H
#define MEMORY_STREAM_H

#include <string>
#include "fileStream.h"

class MemoryStream : public FileStream
{
public:
  int ncid;
  // Constructors
  MemoryStream(int p_ncid);
  MemoryStream(int p_ncid, int cdi_id);
  // ---

  // CdoStream Interface functions
  void close() override;
  int open_read() override;
  int open_write(int p_filetype = CDI_FILETYPE_NC4) override;
  int get_id() override;

  ~MemoryStream() {}

  int create_mem_output();
  //  int open_write(int p_filetype) override;
  //  int open_append() override;
  //
  //  int inq_vlist() override;
  //  void def_vlist(int p_vlistID) override;
  //
  //  void inq_record(int *varID, int *levelID) override;
  //  void defRecord(int varID, int levelID) override;
  //
  //  void read_record(float *const p_data, size_t *numMissVals) override;
  //  void read_record(double *const p_data, size_t *numMissVals) override;
  //  void read_record(Field *const p_field, size_t *numMissVals) override;
  //
  //  void write_record(const float *const p_data, size_t numMissVals) override;
  //  void write_record(const double *const p_data, size_t numMissVals) override;
  //  void write_record(const Field *const p_field, size_t numMissVals) override;
  //
  //  void copyRecord(CdoStreamID p_fileStream) override;
  //
  //  int inq_timestep(int tsID) override;
  //  void def_timestep(int tsID) override;
  //
  //  int inqFileType() override;
  //  int inqByteorder() override;
  //  void defDatarangeList(int p_vlistID);
  //
  //  void close() override;
  //
  //  size_t getNvals() override;

protected:
  static bool TimerEnabled;

private:
  void *ptr;
  MemoryStream() = delete;
  std::string m_filename;
};

#endif
