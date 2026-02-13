/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Oliver Heidmann
          Uwe Schulzweida

*/

#include "processManager.h"
#include "process.h"
#include "cdo_timer.h"
#include "cdo_output.h"
#include "cdo_options.h"
#include "factory.h"
#include "util_string.h"

#include <mutex>

constexpr int MAX_PROCESS = 65536;

static std::mutex processMutex;

static std::string parse_err_msg = "";

static const int IS_OBASE = -1;

extern "C" size_t getPeakRSS();
static std::string
get_max_memstring()
{
  std::stringstream memString;
  auto memMax = getPeakRSS();
  if (memMax)
  {
    constexpr std::array<const char *, 7> memUnitsList = { "B", "KB", "MB", "GB", "TB", "PB", "EB" };
    size_t memUnitsIdx = 0;
    for (; memMax > 9999 && memUnitsIdx < (memUnitsList.size() - 1); ++memUnitsIdx) { memMax /= 1024; }
    memString << " " << memMax << memUnitsList[memUnitsIdx];
  }

  return memString.str();
}

static void
print_benchmarks(double runTime, double readTime, double writeTime)
{
  auto memString = get_max_memstring();
  auto numberOfUsedThreads = get_process_num();
  if (Options::test)
  {
    auto in = std::lround(100 * readTime / runTime);
    auto out = std::lround(100 * writeTime / runTime);
    std::fprintf(stdout, " [%.2fs%s IO:%ld/%ld%% %dthread%s]", runTime, memString.c_str(), in, out, numberOfUsedThreads,
                 ADD_PLURAL(numberOfUsedThreads));
  }
  else
  {
    std::fprintf(stdout, " [%.2fs%s]", runTime, memString.c_str());
  }
}

static void
print_processed_values(Process *p_process, double runTime, double readTime, double writeTime)
{
  set_text_color(stdout, GREEN);
  std::fprintf(stdout, "%s: ", p_process->prompt);
  reset_text_color(stdout);

  auto nvals = p_process->inq_nvals();

  auto nvars = p_process->m_nvars;
  if (nvals > 0)
  {
    std::fprintf(stdout, "Processed %zu value%s from %d variable%s", nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
  }
  else if (nvars > 0) { std::fprintf(stdout, "Processed %d variable%s", nvars, ADD_PLURAL(nvars)); }

  auto ntimesteps = p_process->ntimesteps;
  if ((nvals || nvars) && ntimesteps > 0) std::fprintf(stdout, " over %d timestep%s", ntimesteps, ADD_PLURAL(ntimesteps));

  if (p_process->m_ID == 0) { print_benchmarks(runTime, readTime, writeTime); }

  // if (m_nvars > 0 || nvals > 0 || ntimesteps > 0 || m_ID == 0) std::fprintf(stdout, ".");
  std::fprintf(stdout, "\n");
}

void
ProcessManager::handle_child_construction(std::shared_ptr<Process> &parent, const std::shared_ptr<Node> &child)
{
  Debug(PROCESS, "handling child: %s", child->oper);
  if (child->type == Node::NodeType::OUT_MEM_BUFFER)
  {
    cdo_abort("%s", "Memory Out fild attempted to be attached as child, should only be parents");
  }
  else if (child->isInFile())
  {
    Debug(PROCESS, "Adding FILE in stream: %s", child->oper);
    parent->add_file_in_stream(child->oper);
  }
  else if (child->type == Node::NodeType::IN_MEM_BUFFER)
  {
    Debug(PROCESS, "Adding MEMORY in stream: %s", child->oper);
    parent->add_mem_in_stream(child->ncid);
  }
  else
  {
    auto c_ptr = build_node(child);
    parent->add_child(c_ptr);
    c_ptr->add_parent(parent);
  }
}

void
ProcessManager::buildProcessTree(std::vector<std::shared_ptr<Node>> roots)
{
  Debug(PROCESS, "Building process Tree");
  std::shared_ptr<Node> node = roots[0]->isOperator() ? roots[0] : roots[0]->children[0];

  std::shared_ptr<Process> first_process;
  try
  {
    first_process = create_process(node->oper, split_args(node->arguments));
  }
  catch (std::runtime_error &e)
  {
    cdo_abort("%s: %s", node->oper, e.what());
  }

  if (node->numOut() == IS_OBASE)
  {
    Debug(PROCESS, "Setting obase for %s", node->oper);
    first_process->set_obase(roots[0]->oper);
  }
  else if (node->numOut() > 0)
  {
    for (auto const &n : roots)
    {
      Debug(PROCESS, "adding out files to %s", node->oper);
      Debug(PROCESS, "node type = %d", n->type);
      if (n->type == Node::NodeType::OUTFILE) { first_process->add_file_out_stream(n->oper); }
    }
  }

  for (auto const &c : node->children) { handle_child_construction(first_process, c); }

  set_process_num(m_processes.size());
}

std::shared_ptr<Process>
ProcessManager::build_node(std::shared_ptr<Node> parent_node)
{
  Debug(PROCESS, "Building process for %s", parent_node->oper);
  auto parent_process = create_process(parent_node->oper, split_args(parent_node->arguments));
  for (auto &child_node : parent_node->children) { handle_child_construction(parent_process, child_node); }
  return parent_process;
}

void
ProcessManager::run_processes()
{
  for (auto &idProcessPair : m_processes)
  {
    if (idProcessPair.first)
    {
      /*TEMP*/
      if (!Options::silentMode && (cdo::stdoutIsTerminal || Options::cdoVerbose))
      {
        // MpMO::Print(Green("%s: ") + "Process started", idProcessPair.second->prompt);
        set_text_color(stdout, GREEN);
        std::fprintf(stdout, "%s: ", idProcessPair.second->prompt);
        reset_text_color(stdout);
        std::fprintf(stdout, "Process started\n");
      }
      m_threadIDs.push_back(idProcessPair.second->start_thread());
    }
  }
  m_threadIDs.push_back(pthread_self());
  // MpMO::PrintCerr(Green("%s: ") + "xProcess started", get_process_from_id(0).inq_prompt());
  Process *processZero = get_process_from_id(0).get();

  cdo::timer runTime;
  execute(processZero);

  if (Options::PrintFilename) std::fprintf(stdout, "\n");

  if (!Options::silentMode && (cdo::stdoutIsTerminal || Options::cdoVerbose))
    print_processed_values(processZero, runTime.elapsed(), cdo::readTimer.elapsed(), cdo::writeTimer.elapsed());
}

void
ProcessManager::kill_processes()
{
  for (auto threadID : m_threadIDs)
  {
    if (threadID != pthread_self())
    {
      pthread_cancel(threadID);
      Debug(PROCESS_MANAGER, "process killed: %ld", threadID);
    }
  }
}

void
ProcessManager::clear_processes()
{
  Debug(PROCESS_MANAGER, "Deleting Processes");
  m_processes.clear();
  m_numProcesses = 0;
  m_numProcessesActive = 0;
}

const std::shared_ptr<Process>
ProcessManager::create_process(std::string const &operatorName, std::vector<std::string> const &arguments)
{
  std::shared_ptr<Process> new_process;
  if ((m_numProcesses + 1) >= MAX_PROCESS) { cdo_abort("Limit of %d processes reached!", MAX_PROCESS); }
  auto processID = m_numProcesses++;

  auto it = Factory::find(operatorName, [&operatorName, &processID]()
                          { cdo_abort("Process %s (id:%d) could not be created", operatorName, processID); });

  auto constructor_function = Factory::get_constructor(it);
  auto success = m_processes.insert(std::make_pair(processID, constructor_function(processID, operatorName, arguments)));
  new_process = success.first->second;
  m_numProcessesActive++;

  return new_process;
}

int
ProcessManager::get_num_processes(void)
{
  std::scoped_lock lock(processMutex);
  int pnums = m_processes.size();
  return pnums;
}

int
ProcessManager::get_num_active_processes(void)
{
  std::scoped_lock lock(processMutex);
  int pnums = m_numProcessesActive;
  return pnums;
}

const std::shared_ptr<Process> &
ProcessManager::get_process_from_id(int p_processID)
{
  std::scoped_lock lock(processMutex);

  auto process = m_processes.find(p_processID);
  if (process == m_processes.end()) cdo_abort("Process with ID: %d not found", p_processID);

  return process->second;
}
