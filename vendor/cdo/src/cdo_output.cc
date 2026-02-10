/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <string>
#include <limits>
#include <bitset>

#include "cdo_output.h"

// Debug Switches
int cdoDebug = 0;
int cdoDebugExt = 0;  //  Debug level for the KNMI extensions
                      // Subsystem Debug Switches
unsigned PROCESS = 0;
unsigned PIPE = 0;
unsigned PIPE_STREAM = 0;
unsigned FILE_STREAM = 0;
unsigned PTHREAD = 0;
unsigned PROCESS_MANAGER = 0;
unsigned CDO_NODE = 0;
unsigned PARSER = 0;
unsigned PROCESS_INT = 0;
unsigned FACTORY = 0;
unsigned KVLIST = 0;
unsigned MODULE_INFO = 0;
unsigned ARGUMENTS = 0;

std::string debug_option_string = "DebugLevels:\n"
                                  "     0: off \n"
                                  "     1: all debugs messages enabled\n"
                                  "  Cdi:\n"
                                  "     2: cdi\n"
                                  "     3: memory\n"
                                  "     4: file\n"
                                  "     5: format\n"
                                  "  Cdo:\n"
                                  "     6: CDO\n"
                                  "     7: PipeStream\n"
                                  "     8: FileStream\n"
                                  "     9: Pipe\n"
                                  "    10: Pthread\n"
                                  "    11: Process\n"
                                  "    12: Process manager\n"
                                  "    13: CDO nodes\n"
                                  "    14: Parser\n"
                                  "    15: Process interface\n"
                                  "    16: Factory\n"
                                  "    17: KVList\n"
                                  "    18: Module Info\n"
                                  "    19: Arguments\n";

void
print_debug_options()
{
  std::cout << debug_option_string;
}

namespace cdo
{
void
parse_debug_arguments(std::vector<std::string> const &tokens, unsigned &cdoDebugLevel, unsigned &cdiDebugLevel)
{
  for (auto const &t : tokens)
    {
      if (t.substr(0, 4).compare("ext=") == 0)
        {
          cdoDebugExt = std::stoul(t.substr(5));
          continue;
        }
      unsigned int_token = std::stoul(t);
      switch (int_token)
        {
        case 0:
          cdoDebugLevel = 0;
          cdiDebugLevel = 0;
          break;
        case 1:
          cdoDebugLevel = std::numeric_limits<unsigned>::max();
          cdiDebugLevel = std::numeric_limits<unsigned>::max();
          break;
        case 2: cdiDebugLevel = std::numeric_limits<unsigned>::max(); break;
        case 6: cdoDebugLevel = std::numeric_limits<unsigned>::max(); break;

        default:
          if (int_token > 6) { cdoDebugLevel = (cdoDebugLevel | (1 << (int_token))); }
          else { cdiDebugLevel = (cdiDebugLevel | (1 << (int_token - 1))); }
        }
    }
}

// unused, only for debug reasons
void
print_debug_levels(const unsigned cdoDebugLevel, const unsigned cdiDebugLevel)
{
  // clang-format off
  std::string deb_level = "CDO Debug Levels:\n";
  deb_level   +="PIPE_STREAM:     "  + std::string(PIPE_STREAM ? "ON" : "OFF")     + "\n"
              + "FILE_STREAM:     "  + std::string(FILE_STREAM ? "ON" : "OFF")     + "\n"
              + "PIPE:            "  + std::string(PIPE ? "ON" : "OFF")            + "\n"
              + "PTHREAD:         "  + std::string(PTHREAD ? "ON" : "OFF")         + "\n"
              + "PROCESS:         "  + std::string(PROCESS ? "ON" : "OFF")         + "\n"
              + "PROCESS_MANAGER: "  + std::string(PROCESS_MANAGER ? "ON" : "OFF") + "\n"
              + "CDO_NODE:        "  + std::string(CDO_NODE ? "ON" : "OFF")        + "\n"
              + "PARSER:          "  + std::string(PARSER ? "ON" : "OFF")          + "\n"
              + "PROCESS_INT:     "  + std::string(PROCESS_INT ? "ON" : "OFF")     + "\n"
              + "FACTORY:         "  + std::string(FACTORY ? "ON" : "OFF")         + "\n"
              + "KVLIST:          "  + std::string(KVLIST ? "ON" : "OFF")          + "\n"
              + "MODULE_INFO:     "  + std::string(MODULE_INFO ? "ON" : "OFF")     + "\n"
              + "ARGUMENTS:       "  + std::string(ARGUMENTS ? "ON" : "OFF")       + "\n";
  std::cout << deb_level << std::endl;

  if (cdoDebugLevel  == 1)
    {
      std::bitset<32> cdo_dbg(cdoDebugLevel);
      std::cout << "CDO BITSET: " << cdo_dbg << '\n';
      std::bitset<32> cdi_dbg(cdiDebugLevel);
      std::cout << "CDI BITSET: " << cdi_dbg << '\n';
    }
  // clang-format on
}
void
set_debug(unsigned p_debug_level)
{
  // cdi dbg = 1 << 1
  // cdi mem = 1 << 2
  // cdi mem = 1 << 3
  // cdi format = 1 << 4
  if (p_debug_level & (1u << 6)) p_debug_level = (std::numeric_limits<unsigned>::max());
  if (p_debug_level & (1u << 7)) PIPE_STREAM = 1;
  if (p_debug_level & (1u << 8)) FILE_STREAM = 1;
#ifdef HAVE_LIBPTHREAD
  if (p_debug_level & (1u << 9)) PIPE = 1;
  if (p_debug_level & (1u << 10)) PTHREAD = 1;
#endif
  if (p_debug_level & (1u << 11)) PROCESS = 1;
  if (p_debug_level & (1u << 12)) PROCESS_MANAGER = 1;
  if (p_debug_level & (1u << 13)) CDO_NODE = 1;
  if (p_debug_level & (1u << 14)) PARSER = 1;
  if (p_debug_level & (1u << 15)) PROCESS_INT = 1;
  if (p_debug_level & (1u << 16)) FACTORY = 1;
  if (p_debug_level & (1u << 17)) KVLIST = 1;
  if (p_debug_level & (1u << 18)) MODULE_INFO = 1;
  if (p_debug_level & (1u << 19)) ARGUMENTS = 1;
  MpMO::DebugLevel = p_debug_level;
  cdoDebug = (p_debug_level > 0);
}

bool
dbg()
{
  return (cdoDebug > 0);
}

void
default_exit(std::string msg = "")
{
  (void) msg;
  exit(EXIT_FAILURE);
}

const char *
default_context()
{
  return "cdo init";
}

void (*exitProgram)(std::string) = default_exit;
const char *(*getContext)(void) = default_context;

void
set_exit_function(void (*func)(std::string str))
{
  exitProgram = func;
}

void
set_context_function(const char *(*func)(void) )
{
  getContext = func;
}
}  // namespace cdo

static std::string
getGRB2ErrStr(void)
{
  std::string errorString(
      "To create a CDO application with GRIB2 support use: ./configure --with-eccodes=<ecCodes root directory> ...");
  return errorString;
}

static std::string
getNCErrString(int filetype)
{
  const std::string ncv = (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C)
                              ? "4"
                              : ((filetype == CDI_FILETYPE_NC2) ? "2" : ((filetype == CDI_FILETYPE_NC5) ? "5" : ""));
#ifdef HAVE_LIBNETCDF
  std::string errorString = "CDO was build with a NetCDF version which doesn't support NetCDF" + ncv + " data !";
#else
  std::string errorString = "To create a CDO application with NetCDF" + ncv + " support use: ./configure --with-netcdf=<NetCDF"
                            + ncv + " root directory> ...";
#endif

  return errorString;
}

static std::string
checkForMissingLib(int filetype)
{
  std::string errStr;

  switch (filetype)
    {
    case CDI_FILETYPE_GRB: break;
    case CDI_FILETYPE_GRB2:
      {
        errStr = getGRB2ErrStr();
        break;
      }
    case CDI_FILETYPE_SRV: break;
    case CDI_FILETYPE_EXT: break;
    case CDI_FILETYPE_IEG: break;
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
    case CDI_FILETYPE_NCZARR:
      {
        errStr = getNCErrString(filetype);
        break;
      }
    default: break;
    }

  return errStr;
}

void
cdi_open_error(int cdiErrno, std::string const &format, const char *path)
{
  std::string context = cdo::getContext();
  MpMO::PrintCerr(Red("%s: ") + format + "\n" + std::string(context.size() + 2, ' ') + "%s", context, path,
                  cdiStringError(cdiErrno));

  if (cdiErrno == CDI_ELIBNAVAIL)
    {
      int byteorder;
      auto filetype = cdiGetFiletype(path, &byteorder);
      auto errStr = checkForMissingLib(filetype);
      if (errStr.size()) { MpMO::PrintCerr("%s\n", errStr); }
    }

  if (MpMO::exitOnError) cdo::exitProgram("cdi_open_error");
}

void
query_user_exit(std::string const &argument)
{
  // modified code from NCO
#define USR_RPL_MAX_LNG 10 /* Maximum length for user reply */
#define USR_RPL_MAX_NBR 10 /* Maximum number of chances for user to reply */
  char usr_rpl[USR_RPL_MAX_LNG];
  int usr_rpl_int;
  short nbr_itr = 0;
  size_t usr_rpl_lng = 0;

  // Initialize user reply string
  usr_rpl[0] = 'z';
  usr_rpl[1] = '\0';

  while (!(usr_rpl_lng == 1 && (*usr_rpl == 'o' || *usr_rpl == 'O' || *usr_rpl == 'e' || *usr_rpl == 'E')))
    {
      if (nbr_itr++ > USR_RPL_MAX_NBR)
        {
          (void) fprintf(stderr, "\n%s: ERROR %d failed attempts to obtain valid interactive input.\n", cdo::getContext(),
                         nbr_itr - 1);
          exit(EXIT_FAILURE);
        }

      if (nbr_itr > 1) (void) fprintf(stdout, "%s: ERROR Invalid response.\n", cdo::getContext());
      (void) fprintf(stdout, "%s: %s exists ---`e'xit, or `o'verwrite (delete existing file) (e/o)? ", cdo::getContext(),
                     argument.c_str());
      (void) fflush(stdout);
      if (fgets(usr_rpl, USR_RPL_MAX_LNG, stdin) == nullptr) continue;

      // Ensure last character in input string is \n and replace that with \0
      usr_rpl_lng = std::strlen(usr_rpl);
      if (usr_rpl_lng >= 1)
        if (usr_rpl[usr_rpl_lng - 1] == '\n')
          {
            usr_rpl[usr_rpl_lng - 1] = '\0';
            usr_rpl_lng--;
          }
    }

  // Ensure one case statement for each exit condition in preceding while loop
  usr_rpl_int = (int) usr_rpl[0];
  switch (usr_rpl_int)
    {
    case 'E':
    case 'e': exit(EXIT_SUCCESS); break;
    case 'O':
    case 'o': break;
    default: exit(EXIT_FAILURE); break;
    }
}

std::string
cdo_argv_to_string(std::vector<std::string> const &argv)
{
  std::string s_argv = "";
  for (auto const &x : argv) { s_argv += x + " "; }
  return s_argv;
}
