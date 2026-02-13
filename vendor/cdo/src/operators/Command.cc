/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:
*/

#include <cdi.h>

#include <climits>
#include <string>
#include <vector>
#include <iterator>

#include "cdo_options.h"
#include "cdo_timer.h"
#include "process_int.h"
#include "varray.h"

static bool mode_auto = false;
static bool mode_demo = false;

static VarList varList;

static CdiStreamID gl_streamID = 0;
static int gl_vlistID = 0;
static int gl_varID = 0;
static int gl_numVars = 0;
static int gl_levelID = 0;
static int gl_tsID1 = 0;
static int gl_tsID2 = 0;
static int gl_numSteps = -1;
static Varray<double> gl_data;

static int Done = 0;

static int com_help(std::string const &);
static int com_list(std::string const &);
static int com_quit(std::string const &);
static int com_stat(std::string const &);
static int com_set(std::string const &);
static int com_vars(std::string const &);

namespace
{
struct command_t
{
  int (*func)(std::string const &);  // Function to call to do the job.
  const std::string name;            // User printable name of the function.
  const std::string doc;             // Documentation for this function.
};
}  // namespace

static const std::vector<command_t> commands = { { com_help, "help", "Display this text" },
                                                 { com_help, "?", "Synonym for 'help'" },
                                                 { com_list, "list", "List files in DIR" },
                                                 { com_quit, "quit", "Quit using CDO" },
                                                 { com_stat, "stat", "Statistic for selected field" },
                                                 { com_set, "set", "set variables" },
                                                 { com_vars, "vars", "list variables" } };
static const int ncommands = commands.size();

// Return non-zero if ARG is a valid argument for CALLER, else print an error message and return zero.
/*
static int
valid_argument(char *caller, char *arg)
{
  if (!arg || !*arg)
    {
      std::fprintf(stderr, "%s: Argument required.\n", caller);
      return 0;
    }
  return 1;
}
*/
// Print out help for ARG, or for all of the commands if ARG is not present.
static int
com_help(std::string const &arg)
{
  int printed = 0;

  for (int i = 0; i < ncommands; ++i)
  {
    if (arg.empty() || (arg == commands[i].name))
    {
      printf("%s\t\t%s.\n", commands[i].name.c_str(), commands[i].doc.c_str());
      printed++;
    }
  }

  if (!printed)
  {
    printf("No commands match '%s'. Possibilties are:\n", arg.c_str());
    for (int i = 0; i < ncommands; ++i)
    {
      if (printed == 6)  // Print in six columns
      {
        printed = 0;
        printf("\n");
      }
      printf("%s\t", commands[i].name.c_str());
      printed++;
    }

    if (printed) printf("\n");
  }

  return 0;
}

// List the file(s) named in arg
static int
com_list(std::string const &arg)
{
  (void) (arg);

  return 0;
}

// The user wishes to quit using this program. Just set DONE non-zero.
static int
com_quit(std::string const &arg)
{
  (void) (arg);

  Done = 1;

  return 0;
}

static int
com_stat(std::string const &arg)
{
  (void) (arg);  // unused

  auto const &var = varList.vars[gl_varID];
  auto name = var.name.c_str();

  auto tsID2 = gl_tsID2;
  if (tsID2 == -1) tsID2 = (gl_numSteps != -1) ? gl_numSteps - 1 : INT_MAX - 1;

  for (int tsID = gl_tsID1; tsID <= tsID2; ++tsID)
  {
    auto numFields = streamInqTimestep(gl_streamID, tsID);
    if (numFields == 0)
    {
      if (gl_numSteps == -1) { gl_numSteps = tsID; }
      else { std::fprintf(stderr, "Timestep %d out of range!\n", tsID + 1); }
      break;
    }
    else
    {
      cdo::timer stepTimer;

      auto nlevels = var.nlevels;
      auto gridsize = var.gridsize;
      auto missval = var.missval;

      auto levelID = (nlevels > 1) ? gl_levelID : 0;

      size_t numMissVals;
      streamReadVarSlice(gl_streamID, gl_varID, levelID, gl_data.data(), &numMissVals);

      auto mmm = varray_min_max_mean_mv(gl_data, gridsize, missval);

      std::fprintf(stdout, "%s:  z=%d  t=%d  size=%zu numMissVals=%zu  min=%.5g mean=%.5g max=%.5g [%.2fs]\n", name, levelID + 1,
                   tsID + 1, gridsize, numMissVals, mmm.min, mmm.mean, mmm.max, stepTimer.elapsed());
    }
  }

  return 0;
}

static inline void
paramWarning(const char *name, const char *cstring, const char *endptr)
{
  std::fprintf(stdout, "%s parameter >%s< contains invalid character at position %d!\n", name, cstring,
               (int) (endptr - cstring + 1));
}

int
param2int(std::string const &str)
{
  char *endptr = nullptr;
  int ival = (int) std::strtol(str.c_str(), &endptr, 10);
  if (*endptr) paramWarning("Integer", str.c_str(), endptr);
  return ival;
}

static int
check_tsID(int tsID)
{
  if (gl_numSteps >= 0 && (tsID < -1 || tsID >= gl_numSteps))
  {
    std::fprintf(stdout, "t=%d out of range (max=%d)!\n", (tsID >= 0) ? tsID + 1 : tsID, gl_numSteps);
    return 1;
  }

  return 0;
}

static void
set_tsID(int tsID1, int tsID2)
{
  if (check_tsID(tsID1)) return;
  if (check_tsID(tsID2)) return;
  gl_tsID1 = tsID1;
  gl_tsID2 = tsID2;

  if (mode_auto) com_stat("");
}

static int
check_levelID(int levelID)
{
  auto nlevels = varList.vars[gl_varID].nlevels;
  if (levelID >= nlevels)
  {
    std::fprintf(stdout, "z=%d out of range (max=%d)!\n", levelID + 1, nlevels);
    return 1;
  }

  return 0;
}

static void
set_levelID(int levelID)
{
  if (check_levelID(levelID)) return;
  gl_levelID = levelID;

  if (mode_auto) com_stat("");
}

static int
check_varID(int varID)
{
  if (varID < 0 || varID >= (int) varList.numVars())
  {
    std::fprintf(stdout, "varID out of range (max=%d)!\n", (int) varList.numVars());
    return 1;
  }

  return 0;
}

static void
set_varID(int varID)
{
  if (check_varID(varID)) return;
  gl_varID = varID;

  if (mode_auto) com_stat("");
}

static void
set_t(std::vector<std::string> const &argv)
{
  auto argc = argv.size();
  if (argc == 1)
  {
    std::fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
    return;
  }
  else if (argc > 3)
  {
    std::fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
    return;
  }

  auto t1 = param2int(argv[1]);
  auto t2 = (argc == 3) ? param2int(argv[2]) : t1;
  set_tsID((t1 > 0) ? t1 - 1 : t1, (t2 > 0) ? t2 - 1 : t2);
}

static void
set_z(std::vector<std::string> const &argv)
{
  auto argc = argv.size();
  if (argc == 1)
  {
    std::fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
    return;
  }
  else if (argc > 2)
  {
    std::fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
    return;
  }

  auto z = param2int(argv[1]);
  set_levelID(z - 1);
}

static void
set_var(std::vector<std::string> const &argv)
{
  auto argc = argv.size();
  if (argc == 1)
  {
    std::fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
    return;
  }
  else if (argc > 2)
  {
    std::fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
    return;
  }

  auto &name = argv[1];
  auto lfound = false;
  for (int varID = 0; varID < gl_numVars; ++varID)
  {
    if (name == varList.vars[varID].name)
    {
      lfound = true;
      set_varID(varID);
      break;
    }
  }

  if (!lfound)
  {
    std::fprintf(stdout, "  set %s: Variable name <%s> not found!\n", argv[0].c_str(), name.c_str());
    return;
  }
}

static int
com_set(std::string const &arg)
{
  printf("com_set: >%s<\n", arg.c_str());
  if (arg.empty())
  {
    std::fprintf(stdout, "  command set: argument missing!\n");
    return -1;
  }

  std::istringstream iss(arg);
  std::vector<std::string> argv(std::istream_iterator<std::string>{ iss }, std::istream_iterator<std::string>());

  // for (int i = 0, n = argv.size(); i < n; ++i) printf(">%s<\n", argv[i].c_str());

  if (argv[0] == "t") { set_t(argv); }
  if (argv[0] == "z") { set_z(argv); }
  else if (argv[0] == "var") { set_var(argv); }

  return 0;
}

static int
com_vars(std::string const &arg)
{
  char paramstr[32];

  printf("com_vars: %s %d\n", arg.c_str(), gl_numVars);

  for (int varID = 0; varID < gl_numVars; ++varID)
  {
    auto const &var = varList.vars[varID];
    cdiParamToString(var.param, paramstr, sizeof(paramstr));

    std::fprintf(stdout, "varID=%3d, param=%s, name=%s, longname=\"%s\", units=\"%s\"\n", varID + 1, paramstr, var.name.c_str(),
                 var.longname.c_str(), var.units.c_str());
  }

  return 0;
}

/* Look up NAME as the name of a command, and return a pointer to that command. Return a nullptr pointer if NAME isn't a command
 * name.
 */
static const command_t *
find_command(std::string const &name)
{
  for (int i = 0; i < ncommands; ++i)
    if (name == commands[i].name) return &commands[i];

  return (command_t *) nullptr;
}

// Execute a command line.
static int
execute_line(std::string const &line)
{
  if (line.empty()) return 0;

  // Isolate the command word.
  int i = 0;
  while (line[i] && std::isspace(line[i])) i++;
  int pos = i;
  while (line[i] && !std::isspace(line[i])) i++;
  int count = i;

  auto word = line.substr(pos, count);

  const command_t *command = find_command(word);
  if (!command)
  {
    std::fprintf(stderr, "%s: No such command!\n", word.c_str());
    return -1;
  }
  // Get argument to command, if any.
  while (std::isspace(line[i])) i++;

  pos = i;
  auto args = line.substr(pos);

  // Call the function.
  return (*(command->func))(args);
}

std::string
trim(std::string const &str, std::string const &chars = "\t\n\v\f\r ")
{
  auto first = str.find_first_not_of(chars);
  if (std::string::npos == first) return str;
  auto last = str.find_last_not_of(chars);

  return str.substr(first, (last - first + 1));
}

extern "C" size_t getPeakRSS();

static std::string
peakRSS_string()
{
  char memstring[32] = { "" };
  size_t memmax = getPeakRSS();
  if (memmax)
  {
    size_t muindex = 0;
    static const char *mu[] = { "B", "KB", "MB", "GB", "TB", "PB" };
    static const size_t nmu = sizeof(mu) / sizeof(char *);
    while (memmax > 9999 && muindex < nmu - 1)
    {
      memmax /= 1024;
      muindex++;
    }
    std::snprintf(memstring, sizeof(memstring), "%zu%s", memmax, mu[muindex]);
  }

  return memstring;
}

static void
read_line(std::string const &p_prompt, std::string &line)
{
  fputs(p_prompt.c_str(), stdout);
  if (Options::cdoVerbose)
  {
    fputs(" [", stdout);
    fputs(peakRSS_string().c_str(), stdout);
    fputs("]", stdout);
  }
  fputs("> ", stdout);
  fflush(stdout);

  std::getline(std::cin, line);
}

static void
command_init()
{
  gl_vlistID = streamInqVlist(gl_streamID);

  auto taxisID = vlistInqTaxis(gl_vlistID);
  (void) (taxisID);  // unused

  varList = VarList(gl_vlistID);
  gl_numVars = varList.numVars();

  auto numSteps = varList.numSteps();
  if (numSteps == 0) numSteps = 1;
  gl_numSteps = numSteps;

  gl_data.resize(varList.gridsizeMax());

  set_varID(0);
}

static void
run_demo()
{
  gl_tsID1 = 0;
  gl_tsID2 = -1;
  for (int varID = 0; varID < gl_numVars; ++varID)
  {
    auto const &var = varList.vars[varID];
    auto nlevels = var.nlevels;
    for (int levelID = 0; levelID < nlevels; ++levelID)
    {
      set_varID(varID);
      set_levelID(levelID);
      com_stat(var.name);
    }
  }
}

class Command : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Command",
    .operators = { { "command" }, { "com" }, { "cmd" } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Command> registration = RegisterEntry<Command>();

public:
  void
  init() override
  {

    if (cdo_operator_argc() == 1)
    {
      if (cdo_operator_argv(0) == "auto") { mode_auto = true; }
      else if (cdo_operator_argv(0) == "demo") { mode_demo = true; }
      else { cdo_abort("Unsupported parameter: %s", cdo_operator_argv(0)); }
    }

    gl_streamID = streamOpenRead(cdo_get_stream_name(0));

    command_init();
  }

  void
  run() override
  {
    if (mode_demo) { run_demo(); }
    else
    {
      // Loop reading and executing lines until the user quits.
      const std::string custom_prompt = "cdo cmd";
      while (!Done)
      {
        std::string line;
        read_line(custom_prompt, line);
        execute_line(trim(line));
      }
    }
  }

  void
  close() override
  {
    streamClose(gl_streamID);
  }
};
