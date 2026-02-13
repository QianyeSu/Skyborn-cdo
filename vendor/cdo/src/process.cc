/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#ifdef HAVE_NETCDF
#include <netcdf.h>
#include <netcdf_mem.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>

#include "process.h"
#include "cdo_options.h"
#include "fileStream.h"
#include "pipeStream.h"

// temporary include for setting the local process
#include "process_int.h"

static int processNum = 0;

int
get_process_num()
{
  return processNum;
}

void
set_process_num(int p_num)
{
  processNum = p_num;
}

std::string
Process::replace_alias(std::string const &p_calledBy, const CdoModule &p_module)
{
  std::string originalName = p_calledBy;
  int aliasID = p_module.is_alias(p_calledBy);
  if (aliasID != -1) originalName = p_module.aliases[aliasID].original;
  return originalName;
}

Process::Process(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments,
                 const CdoModule &p_module)
    : arguments(p_module.arguments), m_ID(p_ID), m_module(p_module), m_oargv(p_arguments)
{
#ifdef HAVE_LIBPTHREAD
  threadID = pthread_self();
#endif

  operatorName = replace_alias(p_operatorName, p_module);
  def_prompt(operatorName);
}

int
Process::get_stream_cnt_in()
{
  return inputStreams.size();
}

int
Process::get_stream_cnt_out()
{
  return outputStreams.size();
}

void
Process::def_prompt(std::string const &name)
{
  if (m_ID == 0)
    std::snprintf(prompt, sizeof(prompt), "%s    %s", cdo::progname, name.c_str());
  else
    std::snprintf(prompt, sizeof(prompt), "%s(%d) %s", cdo::progname, m_ID, name.c_str());
}

const char *
Process::inq_prompt() const
{
  return prompt;
}

int
Process::get_operator_id()
{
  if (m_module.operators.size() <= 0) { cdo_abort("Operator not initialized!"); }

  Debug(PROCESS, "searching for %s", operatorName);
  for (size_t operID = 0; operID < m_module.operators.size(); operID++)
  {
    Debug(PROCESS, "comparing  %s and %s", operatorName, m_module.operators[operID].name);
    if (operatorName == m_module.operators[operID].name) return operID;
  }

  cdo_abort("Operator not callable by this name! Name is: %s", operatorName);

  return -1;
}

void
Process::add_file_in_stream(std::string const &file)
{
  inputStreams.push_back(std::make_shared<FileStream>(file));
  m_streamCnt++;
}

#include "memoryStream.h"
void
Process::add_mem_in_stream(int ncid)
{
  Debug(PROCESS, "adding memory stream with ncid: %d", ncid);
  inputStreams.push_back(std::make_shared<MemoryStream>(ncid));
  m_streamCnt++;
}

void
Process::add_mem_out_stream(int const &ncid)
{
  outputStreams.push_back(std::make_shared<MemoryStream>(ncid));
  m_streamCnt++;
}
void
Process::add_mem_out_stream(int const &ncid, int const &file_id)
{
  outputStreams.push_back(std::make_shared<MemoryStream>(ncid, file_id));
  m_streamCnt++;
}

void
Process::add_file_out_stream(std::string const &file)
{
  if (file[0] == '-') { cdo_abort("Missing output file. Found an operator instead of filename: %s", file); }
  outputStreams.push_back(std::make_shared<FileStream>(file));
  m_streamCnt++;
}

void
Process::add_child(const std::shared_ptr<Process> &childProcess)
{
  childProcesses.push_back(childProcess);
  add_pipe_in_stream();
}

void
Process::add_pipe_in_stream()
{
#ifdef HAVE_LIBPTHREAD
  inputStreams.push_back(std::make_shared<PipeStream>(m_ID));
  m_streamCnt++;
#else
  cdo_abort("Cannot use pipes, pthread support not compiled in!");
#endif
}

void
Process::add_parent(const std::shared_ptr<Process> &parentProcess)
{
  m_posInParent = parentProcess->inputStreams.size() - 1;
  outputStreams.push_back(parentProcess->inputStreams[m_posInParent]);
  m_streamCnt++;
}

void *
execute(void *process)
{
  Process *p = (Process *) process;

  p->cdo_initialize();

  p->init();
  p->run();
  p->close();

  p->cdo_finish();

  return nullptr;
}

pthread_t
Process::start_thread()
{
  Debug(PROCESS, "starting new thread for process %d", m_ID);
  pthread_attr_t attr;
  auto status = pthread_attr_init(&attr);
  if (status) cdo_sys_error("pthread_attr_init failed for '%s'", operatorName);
  status = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (status) cdo_sys_error("pthread_attr_setdetachstate failed for '%s'", operatorName);
  /*
    param.sched_priority = 0;
    status = pthread_attr_setschedparam(&attr, &param);
    if ( status ) cdo_sys_error("pthread_attr_setschedparam failed for '%s'", newarg+1);
  */
  /* status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED); */
  /* if ( status ) cdo_sys_error("pthread_attr_setinheritsched failed for '%s'", newarg+1); */

  int pthreadScope;
  pthread_attr_getscope(&attr, &pthreadScope);

  /* status = pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS); */
  /* if ( status ) cdo_sys_error("pthread_attr_setscope failed for '%s'", newarg+1); */
  /* If system scheduling scope is specified, then the thread is scheduled against all threads in the system */
  /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */

  size_t stacksize = 0;
  status = pthread_attr_getstacksize(&attr, &stacksize);
  if (stacksize < 2097152)
  {
    stacksize = 2097152;
    pthread_attr_setstacksize(&attr, stacksize);
  }

  pthread_t thrID;
  auto rval = pthread_create(&thrID, &attr, execute, this);
  if (rval != 0)
  {
    errno = rval;
    cdo_sys_error("pthread_create failed for '%s'", operatorName);
  }

  return thrID;
}

// local helper function

bool
Process::has_out_stream(CdoStreamID p_streamID)
{
  for (auto const &streamID : outputStreams)
  {
    if (streamID == p_streamID) return true;
  }
  return false;
}

bool
Process::has_in_stream(CdoStreamID p_streamID)
{
  for (auto const &streamID : inputStreams)
  {
    if (streamID == p_streamID) return true;
  }
  return false;
}

size_t
Process::inq_nvals()
{
  size_t nvals = 0;
  for (size_t i = 0, n = inputStreams.size(); i < n; ++i)
  {
    Debug(PROCESS, "Inquiring nvals from instream %s", inputStreams[i]->m_name);
    nvals += inputStreams[i]->getNvals();
  }
  return nvals;
}

bool
Process::has_no_pipes()
{
  return (childProcesses.size() == 0);
}

const char *
Process::get_out_stream_name()
{
  return outputStreams[0]->m_name.c_str();
}

size_t
Process::get_oper_argc()
{
  return m_oargv.size();
}

std::string
Process::get_argv(int p_idx)
{
  if (!(p_idx > (int) get_oper_argc() && p_idx > 0))
    cdo_abort("Process Argv not found. Idx: %d, Process argc: %d", p_idx, m_oargv.size());

  return m_oargv[p_idx];
}

std::tuple<int, int>
Process::create_output()
{
  int ncid = -1;

#ifdef HAVE_NETCDF
  if (auto retVal = nc_create_mem("test_name", 0, 4096, &ncid))
  {
    printf("Error: %s\n", nc_strerror(retVal));
    std::exit(1);
  }
  if (PROCESS) std::cout << "created ncid: " << ncid << std::endl;

  int streamID = streamOpenWriteNCMem(ncid);
#else
  int streamID = -100103;
#endif

  std::stringstream ss = std::stringstream();
  ss << "ERROR: could not open output stream to memory: errcode: ";
  ss << std::to_string(streamID);
  if (streamID < 0) { throw std::runtime_error(ss.str()); }
  return { ncid, streamID };
}

CdoStreamID
Process::open_write(std::string const &p_filename, int filetype)
{
  if (filetype == CDI_UNDEFID) filetype = cdo_filetype();

  if (write_mode == ProcessWriteMode::FILEIO) { add_file_out_stream(p_filename); }
  else if (write_mode == ProcessWriteMode::MEMORY)
  {
    auto [fid, sid] = create_output();
    add_mem_out_stream(fid, sid);
  }

  const auto pstreamID = outputStreams.back()->open_write(filetype);
  if (pstreamID == -1) cdo_abort("Could not create pstream for file: %s", p_filename);

  return outputStreams.back();
}

void
Process::set_obase(std::string const &obase)
{
  m_obase = obase;
}  // TODO into cc

std::string const &
Process::get_obase()
{
  return m_obase;
}

void
Process::close_streams()
{
  for (auto &s : inputStreams) { s->close(); }
  for (auto &s : outputStreams) { s->close(); }
}

int
Process::get_id()
{
  return m_ID;
}

void
Process::cdo_initialize()
{
#ifdef _OPENMP
  omp_set_num_threads(Threading::ompNumMaxThreads);  // Has to be called for every module (pthread)!
#endif
#ifdef HAVE_LIBPTHREAD
  threadID = pthread_self();
#endif
  Debug(PROCESS_INT, "Initializing process: %s (id: %d)", operatorName, get_id());

  set_local_process(this);
#ifdef HAVE_LIBPTHREAD
  Debug(PROCESS_INT, "process %d thread %ld", m_ID, pthread_self());
#endif
}

void
Process::cdo_finish(void)
{
  Debug(PROCESS_INT, "Finishing process: %d", get_id());

#ifdef HAVE_LIBPTHREAD
  Debug(PROCESS_INT, "process %d thread %ld", m_ID, pthread_self());
#endif
}
