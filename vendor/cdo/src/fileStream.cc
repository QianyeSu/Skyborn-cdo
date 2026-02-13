/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_NC4HDF5_THREADSAFE */
#endif

#include <sys/stat.h> /* stat */
#include <cdi.h>

#include "fileStream.h"

#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_default_values.h"
#include "cdo_history.h"
#include "cdo_query.h"
#include "cdo_timer.h"
#include "commandline.h"

bool FileStream::TimerEnabled = false;
#ifndef HAVE_NC4HDF5_THREADSAFE
static auto inputFileTypeIsNetCDF4 = false;
#endif

FileStream::FileStream(std::string const &p_fileName) : m_filename(p_fileName)
{
  m_name = p_fileName;
  m_fileID = CDI_UNDEFID;
}

int
FileStream::get_id()
{
  return m_fileID;
}

int
FileStream::open_read()
{
  int fileID = -1;

  if (FileStream::timersEnabled()) cdo::readTimer.start();
  open_lock();

  if (m_filename.size() > 6 && m_filename.rfind("query:", 0) == 0)
  {
    CdiQuery *query = cdiQueryCreate();
    auto path = set_query_parameter(m_filename.substr(6), query);
    if (Options::cdoVerbose) cdiQueryPrint(query);

    fileID = streamOpenReadQuery(path.c_str(), query);
  }
  else if (Options::cdoQueryParameter.size() > 0)
  {
    CdiQuery *query = cdiQueryCreate();
    (void) set_query_parameter(Options::cdoQueryParameter, query);
    if (Options::cdoVerbose) cdiQueryPrint(query);

    fileID = streamOpenReadQuery(m_filename.c_str(), query);
  }
  else { fileID = streamOpenRead(m_filename.c_str()); }

  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_filename.c_str());
  isopen = true;

  m_filetype = streamInqFiletype(fileID);
  if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = m_filetype;
  m_fileID = fileID;

#ifndef HAVE_NC4HDF5_THREADSAFE
  if (m_filetype == CDI_FILETYPE_NC4 || m_filetype == CDI_FILETYPE_NC4C) inputFileTypeIsNetCDF4 = true;
#endif

  Debug(FILE_STREAM, "fileID: %d  path: %s", m_fileID, m_name);
  if (Options::numStreamWorker > 0) streamDefNumWorker(fileID, Options::numStreamWorker);

  open_unlock();
  if (FileStream::timersEnabled()) cdo::readTimer.stop();

  return fileID;
}

int
FileStream::open_write(int p_filetype)
{
  if (Options::cdoInteractive)
  {
    struct stat stbuf;
    auto rstatus = stat(m_name.c_str(), &stbuf);
    // If permanent file already exists, query user whether to overwrite or exit
    if (rstatus != -1) query_user_exit(m_name);
  }
  if (p_filetype == CDI_UNDEFID) p_filetype = CDI_FILETYPE_GRB;
  /*
#ifndef HAVE_NC4HDF5_THREADSAFE
  auto outputFileTypeIsNetCDF4 = (p_filetype == CDI_FILETYPE_NC4 || p_filetype == CDI_FILETYPE_NC4C);
  if (inputFileTypeIsNetCDF4 && outputFileTypeIsNetCDF4 && get_process_num() > 1 && Threading::cdoLockIO == false)
    {
      cdo_warning("Using a non-thread-safe NetCDF4/HDF5 library in a multi-threaded environment may lead to erroneous results!");
      cdo_warning("Use a thread-safe NetCDF4/HDF5 library or the CDO option -L to avoid such errors.");
    }
#endif
  */
  // TODO FIX THIS: if (FileStream::timersEnabled()) cdo::writeTimer.start();

  open_lock();
  auto fileID = streamOpenWrite(m_filename.c_str(), p_filetype);
  open_unlock();

  // TODO FIX THIS: if(FileStream::timersEnabled()) cdo::writeTimer.stop();
  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_name.c_str());
  isopen = true;

  if (CdoDefault::Byteorder != CDI_UNDEFID) streamDefByteorder(fileID, CdoDefault::Byteorder);

  set_compression(fileID, p_filetype);

  m_fileID = fileID;
  m_filetype = p_filetype;

  Debug(FILE_STREAM, "fileID: %d  path: %s", m_fileID, m_name);

  if (Options::PrintFilename) std::fprintf(stdout, "%s ", m_name.c_str());

  return m_cdoStreamID;
}

int
FileStream::open_append()
{
  if (FileStream::timersEnabled()) cdo::writeTimer.start();

  open_lock();
  auto fileID = streamOpenAppend(m_filename.c_str());
  open_unlock();

  if (FileStream::timersEnabled()) cdo::writeTimer.stop();

  if (fileID < 0) cdi_open_error(fileID, "Open failed on >%s<", m_filename.c_str());

  isopen = true;

  m_filetype = streamInqFiletype(fileID);
  set_compression(fileID, m_filetype);

  m_fileID = fileID;

  return m_fileID;
}

void
FileStream::def_vlist(int p_vlistID)
{
  if (m_filetype == CDI_FILETYPE_NCZARR)
  {
    auto maxSteps = vlistNtsteps(p_vlistID);
    if (maxSteps >= 0)
      streamDefMaxSteps(m_fileID, maxSteps);
    else
      cdo_warning("Unknown number of timesteps, use operator setmaxsteps to set max. number of timesteps!");
  }

  cdo_append_history(p_vlistID, cdo::command_line());

  if (CdoDefault::DataType != CDI_UNDEFID)
  {
    auto nvars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < nvars; ++varID) vlistDefVarDatatype(p_vlistID, varID, CdoDefault::DataType);
    if (CdoDefault::DataType == CDI_DATATYPE_FLT64 || CdoDefault::DataType == CDI_DATATYPE_FLT32)
    {
      for (int varID = 0; varID < nvars; ++varID)
      {
        double addoffset = 0.0, scalefactor = 1.0;
        auto haveAddoffset = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
        auto haveScalefactor = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
        if (haveAddoffset || haveScalefactor)
        {
          cdiDeleteKey(p_vlistID, varID, CDI_KEY_ADDOFFSET);
          cdiDeleteKey(p_vlistID, varID, CDI_KEY_SCALEFACTOR);
        }
      }
    }
  }

  if (Options::nsb > 0)
  {
    auto nvars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < nvars; ++varID) vlistDefVarNSB(p_vlistID, varID, Options::nsb);
  }

  if (Options::cdoChunkType != CDI_UNDEFID)
  {
    auto numVars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKTYPE, Options::cdoChunkType);
    auto numGrids = vlistNumGrids(p_vlistID);
    for (int index = 0; index < numGrids; ++index)
      cdiDefKeyInt(vlistGrid(p_vlistID, index), CDI_GLOBAL, CDI_KEY_CHUNKTYPE, Options::cdoChunkType);
  }

  if (Options::cdoChunkSize != CDI_UNDEFID)
  {
    auto numVars = vlistNvars(p_vlistID);
    constexpr int chunkKeys[] = { CDI_KEY_CHUNKSIZE_DIMT, CDI_KEY_CHUNKSIZE_DIMZ, CDI_KEY_CHUNKSIZE_DIMY, CDI_KEY_CHUNKSIZE_DIMX };
    for (int varID = 0; varID < numVars; ++varID)
    {
      for (auto chunkKey : chunkKeys)
      {
        if (cdo::inq_key_int(p_vlistID, varID, chunkKey) != 0) cdiDeleteKey(p_vlistID, varID, chunkKey);
      }
    }

    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE, Options::cdoChunkSize);
    auto numGrids = vlistNumGrids(p_vlistID);
    for (int index = 0; index < numGrids; ++index)
      cdiDefKeyInt(vlistGrid(p_vlistID, index), CDI_GLOBAL, CDI_KEY_CHUNKSIZE, Options::cdoChunkSize);
  }

  if (Options::cdoChunkSizeDimX != 0)
  {
    auto numVars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE_DIMX, Options::cdoChunkSizeDimX);
  }

  if (Options::cdoChunkSizeDimY != 0)
  {
    auto numVars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE_DIMY, Options::cdoChunkSizeDimY);
  }

  if (Options::cdoChunkSizeDimZ != 0)
  {
    auto numVars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE_DIMZ, Options::cdoChunkSizeDimZ);
    auto numZaxis = vlistNumZaxis(p_vlistID);
    for (int index = 0; index < numZaxis; ++index)
      cdiDefKeyInt(vlistZaxis(p_vlistID, index), CDI_GLOBAL, CDI_KEY_CHUNKSIZE_DIMZ, Options::cdoChunkSizeDimZ);
  }

  if (Options::cdoChunkSizeDimT != 0)
  {
    auto numVars = vlistNvars(p_vlistID);
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyInt(p_vlistID, varID, CDI_KEY_CHUNKSIZE_DIMT, Options::cdoChunkSizeDimT);
  }

  if (Options::CMOR_Mode)
  {
    cdo_def_tracking_id(p_vlistID, "tracking_id");
    cdo_def_creation_date(p_vlistID);
  }

  if (Options::VersionInfo) cdiDefAttTxt(p_vlistID, CDI_GLOBAL, "CDO", (int) std::strlen(cdo_comment()), cdo_comment());

#ifdef _OPENMP
  if (Threading::ompNumMaxThreads > 1)
    cdiDefAttInt(p_vlistID, CDI_GLOBAL, "cdo_openmp_thread_number", CDI_DATATYPE_INT32, 1, &Threading::ompNumMaxThreads);
#endif
  defDatarangeList(p_vlistID);

  if (FileStream::timersEnabled()) cdo::writeTimer.start();
  stream_def_vlist_locked(m_fileID, p_vlistID);
  if (FileStream::timersEnabled()) cdo::writeTimer.stop();
}

int
FileStream::inq_vlist()
{
  if (FileStream::timersEnabled()) cdo::readTimer.start();
  auto vlistID = stream_inq_vlist_locked(m_fileID);
  if (FileStream::timersEnabled()) cdo::readTimer.stop();
  if (vlistID == -1) cdo_abort("Couldn't read data from input fileID %d!", m_fileID);

  auto nsubtypes = vlistNsubtypes(vlistID);
  if (nsubtypes > 1) cdo_warning("Subtypes are unsupported, the processing results are possibly wrong!");

  if (CdoDefault::TaxisType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(vlistID), CdoDefault::TaxisType);

  m_vlistID = vlistID;
  return m_vlistID;
}

void
FileStream::inq_field(int *const varID, int *const levelID)
{
  if (FileStream::timersEnabled()) cdo::readTimer.start();
  stream_inq_field_locked(m_fileID, varID, levelID);
  if (FileStream::timersEnabled()) cdo::readTimer.stop();
  m_varID = *varID;
}

void
FileStream::def_field(int varID, int levelID)
{
  if (FileStream::timersEnabled()) cdo::writeTimer.start();
  stream_def_field_locked(m_fileID, varID, levelID);
  if (FileStream::timersEnabled()) cdo::writeTimer.stop();
  m_varID = varID;
}

void
FileStream::read_field(float *const p_data, size_t *const numMissVals)
{
  if (FileStream::timersEnabled()) cdo::readTimer.start();
  stream_read_field_float_locked(m_fileID, p_data, numMissVals);
  if (FileStream::timersEnabled()) cdo::readTimer.stop();
}

void
FileStream::read_field(double *const p_data, size_t *const numMissVals)
{
  if (FileStream::timersEnabled()) cdo::readTimer.start();
  stream_read_field_double_locked(m_fileID, p_data, numMissVals);
  if (FileStream::timersEnabled()) cdo::readTimer.stop();
}

void
FileStream::read_field(Field *const p_field, size_t *const numMissVals)
{
  read_field(p_field->vec_d.data(), numMissVals);
}

void
FileStream::write_field(const float *const p_data, size_t p_numMissVals)
{
  if (FileStream::timersEnabled()) cdo::writeTimer.start();

  auto varID = m_varID;
  if (varID < (int) m_datarangelist.size())
    if (m_datarangelist[varID].checkDatarange) m_datarangelist[varID].check_datarange(p_data, p_numMissVals);

  Debug(FILE_STREAM, "writing");
  stream_write_field_float_locked(m_fileID, p_data, p_numMissVals);

  if (FileStream::timersEnabled()) cdo::writeTimer.stop();
}

void
FileStream::write_field(const double *const p_data, size_t p_numMissVals)
{
  if (FileStream::timersEnabled()) cdo::writeTimer.start();

  auto varID = m_varID;
  if (varID < (int) m_datarangelist.size())
    if (m_datarangelist[varID].checkDatarange) m_datarangelist[varID].check_datarange(p_data, p_numMissVals);

  stream_write_field_double_locked(m_fileID, p_data, p_numMissVals);

  if (FileStream::timersEnabled()) cdo::writeTimer.stop();
}

void
FileStream::write_field(const Field *const p_field, size_t p_numMissVals)
{
  write_field(p_field->vec_d.data(), p_numMissVals);
}

void
FileStream::copy_field(CdoStreamID p_destination)
{
  FileStream *fStream = dynamic_cast<FileStream *>(p_destination.get());
  stream_copy_field_locked(m_fileID, fStream->getFileID());
}
/*
 * FileStream::inq_timestep(int p_tsID)
 * stets internal state of the cdi datastructure to work on the given timestep (p_tsID) and returns the number of fields that the
 * timestep contains.
 * Inquires and defines the time axis type if the timestep ID is 0 AND the taxis type is yet to be defined.
 * When the timestep inquiry was successfull m_tsID is set to the wanted p_tsID IF p_tsID != m_tsID
 * When only one process is running the timers are enabled.
 * -- last Documentation update(2019-06-14) --
 */
int
FileStream::inq_timestep(int p_tsID)
{
  if (FileStream::timersEnabled()) cdo::readTimer.start();
  auto numFields = stream_inq_time_step_locked(m_fileID, p_tsID);
  if (FileStream::timersEnabled()) cdo::readTimer.stop();

  if (p_tsID == 0 && CdoDefault::TaxisType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(m_vlistID), CdoDefault::TaxisType);

  if (numFields && p_tsID != m_tsID) m_tsID = p_tsID;
  Debug(FILE_STREAM, "Current TsID: %d,  numFields: %d", m_tsID, numFields);

  return numFields;
}

void
FileStream::def_timestep(int p_tsID)
{
  if (FileStream::timersEnabled()) cdo::writeTimer.start();
  // don't use sync -> very slow on GPFS
  //  if ( p_tsID > 0 ) streamSync(fileID);

  stream_def_time_step_locked(m_fileID, p_tsID);

  if (FileStream::timersEnabled()) cdo::writeTimer.stop();
}

int
FileStream::inqFileType()
{
  return streamInqFiletype(m_fileID);
}

int
FileStream::inqByteorder()
{
  return streamInqByteorder(m_fileID);
}

size_t
FileStream::getNvals()
{
  // set when the stream is closed
  // see: FileStream::close()
  return m_nvals;
}

int
FileStream::getFileID()
{
  return m_fileID;
}

void
FileStream::close()
{
  Debug(FILE_STREAM, "fileID: %d  path: %s", m_fileID, m_name);

  m_nvals = streamNvals(m_fileID);
  stream_close_locked(m_fileID);

  isopen = false;
  m_vlistID = -1;

  if (m_datarangelist.size())
  {
    m_datarangelist.clear();
    m_datarangelist.shrink_to_fit();
  }
}

void
FileStream::defDatarangeList(int p_vlistID)
{
  auto filetype = m_filetype;

  if (m_vlistID != -1) cdo_abort("Internal problem, vlist already defined!");

  if (m_datarangelist.size() != 0) cdo_abort("Internal problem, datarangelist already allocated!");

  auto nvars = vlistNvars(p_vlistID);
  assert(nvars > 0);

  m_datarangelist.resize(nvars);

  for (int varID = 0; varID < nvars; ++varID)
  {
    m_datarangelist[varID].gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
    m_datarangelist[varID].datatype = vlistInqVarDatatype(p_vlistID, varID);
    m_datarangelist[varID].missval = vlistInqVarMissval(p_vlistID, varID);

    double addoffset = 0.0, scalefactor = 1.0;
    auto haveAddoffset = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
    auto haveScalefactor = (cdiInqKeyFloat(p_vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
    if (haveAddoffset) m_datarangelist[varID].addoffset = addoffset;
    if (haveScalefactor) m_datarangelist[varID].scalefactor = scalefactor;

    m_datarangelist[varID].checkDatarange = false;

    auto datatype = m_datarangelist[varID].datatype;

    if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C
        || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR)
    {
      if (datatype == CDI_DATATYPE_UINT8
          && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
      {
        datatype = CDI_DATATYPE_INT16;
        m_datarangelist[varID].datatype = datatype;
      }

      if (datatype == CDI_DATATYPE_UINT16
          && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
      {
        datatype = CDI_DATATYPE_INT32;
        m_datarangelist[varID].datatype = datatype;
      }

      if (haveAddoffset || haveScalefactor)
      {
        if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
            || datatype == CDI_DATATYPE_UINT16)
          m_datarangelist[varID].checkDatarange = true;
      }
      else if (Options::CheckDatarange) { m_datarangelist[varID].checkDatarange = true; }
    }
  }

  m_vlistID = p_vlistID;  // used for -r/-a
}
