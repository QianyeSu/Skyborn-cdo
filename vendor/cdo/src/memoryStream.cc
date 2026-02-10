/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <sys/stat.h> /* stat */
#include <cdi.h>

#include "memoryStream.h"

#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_default_values.h"
#include "cdo_history.h"
#include "cdo_query.h"
#include "cdo_process.h"

#include "commandline.h"

MemoryStream::MemoryStream(int p_ncid) : FileStream()
{
  Debug(FILE_STREAM, "creating mem stream: p_ncid = %d", ncid);
  ncid = p_ncid;
  m_fileID = p_ncid;
}

MemoryStream::MemoryStream(int p_ncid, int cdi_id) : FileStream()
{
  Debug(FILE_STREAM, "creating mem stream: p_ncid = %d, cdi_id = %d", ncid, cdi_id);
  ncid = p_ncid;
  m_fileID = cdi_id;
}

int
MemoryStream::open_read()
{
  if (ncid < 0) cdi_open_error(ncid, "Open failed on >%s<", m_filename.c_str());
  isopen = true;

  m_filetype = streamInqFiletype(ncid);
  if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = m_filetype;
  m_fileID = ncid;

  Debug(FILE_STREAM, "Set number of worker to %d", Options::numStreamWorker);
  if (Options::numStreamWorker > 0) streamDefNumWorker(ncid, Options::numStreamWorker);

  open_unlock();

  return m_fileID;
}
int
MemoryStream::open_write(int p_filetype)
{
  Debug(FILE_STREAM, "Open write in memoryStream called");

  if (m_fileID < 0) cdi_open_error(m_fileID, "Open failed on >%s<", m_name.c_str());
  isopen = true;

  if (CdoDefault::Byteorder != CDI_UNDEFID)
  {
    std::cout << "attempting streamDefByteorder" << std::endl;
    streamDefByteorder(m_fileID, CdoDefault::Byteorder);
  }

  set_compression(m_fileID, CDI_FILETYPE_NC4);

  m_filetype = CDI_FILETYPE_NC4;
  Debug(FILE_STREAM, "finished open_write with filetype set to: %d", inqFileType());

  return m_fileID;
}

void
MemoryStream::close()
{
  Debug(FILE_STREAM, "%s fileID %d", m_name, m_fileID);

  m_nvals = streamNvals(m_fileID);

  streamCloseNCMem(m_fileID);

  isopen = false;
  m_vlistID = -1;

  if (m_datarangelist.size())
  {
    m_datarangelist.clear();
    m_datarangelist.shrink_to_fit();
  }
  Debug(FILE_STREAM, "called close function of MemoryStream, this will do nothing, destruction is handled by deconstructor");
}

int
MemoryStream::get_id()
{
  Debug(FILE_STREAM, "getting id from memory stream");
  return ncid;
}
