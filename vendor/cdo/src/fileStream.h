#ifndef FILESTREAM_H
#define FILESTREAM_H

#include <string>
#include "cdoStream.h"

class FileStream : public CdoStream
{
public:
  // Constructors
  explicit FileStream(std::string const &p_fileStream);
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
  void defDatarangeList(int p_vlistID);

  void close() override;

  size_t getNvals() override;
  // ---

  // FileStreamOnly
  int getFileID();

  static void
  enableTimers(bool p_enable)
  {
    FileStream::TimerEnabled = p_enable;
  }

  static bool
  timersEnabled()
  {
    return FileStream::TimerEnabled;
  }
  // ---

protected:
  static bool TimerEnabled;

private:
  std::string m_filename;
  void checkDatarange(int varID, double *array, size_t numMissVals);

protected:
  FileStream() = default;
};

#endif
