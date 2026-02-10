/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef PROCESS_H
#define PROCESS_H

#include "cdoStream.h"
#include "cdo_module.h"

#include <vector>
#include <string>

enum struct ProcessStatus
{
  Ok = 0,
  UnlimitedIOCounts = -1,
  MissInput = -2,
  MissOutput = -3,
  TooManyStreams = -4,
  TooFewStreams = -5,
};

enum struct ProcessWriteMode
{
  FILEIO = 0,
  MEMORY = 1
};

void *execute(void *process);

class Process
{
private:
  ArgumentHandler arguments{};

public:
  // Member Variables

  ProcessWriteMode write_mode = ProcessWriteMode::FILEIO;
  int m_ID{};
  const CdoModule &m_module;
  int m_posInParent{};
  std::vector<std::shared_ptr<Process>> childProcesses{};
  std::vector<CdoStreamID> inputStreams{};
  std::vector<CdoStreamID> outputStreams{};

  int m_nvars = 0;
  int ntimesteps = 0;
  int m_streamCnt = 0;

  char prompt[64];

  std::string m_operatorCommand = "UNINITALIZED";
  std::string operatorName{};
  std::string m_obase{};
  std::vector<std::string> m_oargv{};

#ifdef HAVE_LIBPTHREAD
  pthread_t threadID{};
#endif

  /* Member Functions  */
  Process(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module);
  virtual ~Process() { Debug(PROCESS, "destruction of %s", operatorName); }

  pthread_t start_thread();
  virtual void init() = 0;
  virtual void run() = 0;
  virtual void close() = 0;

  std::tuple<int, int> create_output();

  CdoStreamID open_write(std::string const &p_filename, int filetype = CDI_UNDEFID);
  void
  parse_arguments()
  {
    arguments.parse(m_oargv);
  }

  template <typename T>
  void
  get_argument(std::string const &key, T &destination)
  {
    arguments.get(key, destination);
  }
  /**
   * returns the number of in streams this process currently has.
   **/
  int get_stream_cnt_in();
  /**
   * returns the number of out streams this process currently has.
   */
  int get_stream_cnt_out();

  /**
   * Adds a Process as child and creates and adds a new pipe stream.
   */
  void add_child(const std::shared_ptr<Process> &child_process);
  /**
   * Adds a Process as parent and adds the parents input stream to out streams.
   */
  void add_parent(const std::shared_ptr<Process> &parent_process);
  /**
   * Adds and creates a new file pstream to the in streams
   */
  void add_file_in_stream(std::string const &file);
  /**
   * Adds and creates a new file pstream to the out streams
   */
  void add_file_out_stream(std::string const &file);
  /**
   * Adds ands creates a memory stream to the out streams vector*/
  void add_mem_out_stream(int const &ncid);
  void add_mem_out_stream(int const &ncid, int const &file_id);
  /**
   * Adds a already created memory region to the process as input
   */
  void add_mem_in_stream(int ncid);
  /**
   * Adds and creates a new pipe pstream to the in streams
   */
  void add_pipe_in_stream();
  /**
   * Adds and creates a new file pstream to the out streams
   */
  void add_pipe_out_stream();
  /**
   * returns the operatorID of the currently in use operator
   */
  int get_operator_id();
  const char *inq_prompt() const;

  const char *get_out_stream_name();
  bool has_out_stream(CdoStreamID p_streamPtr);
  bool has_in_stream(CdoStreamID p_streamPtr);

  bool has_no_pipes();

  size_t inq_nvals();

  size_t get_oper_argc();
  std::string get_argv(int idx);

  std::string const &get_obase();
  void set_obase(std::string const &obase);

  int get_id();

  void close_streams();

  std::string replace_alias(std::string const &p_calledBy, const CdoModule &p_module);
  void def_prompt(std::string const &name);

  void cdo_initialize();
  void cdo_finish(void);
};

int get_process_num();
void set_process_num(int p_num);

#endif /* PROCESS_H */
