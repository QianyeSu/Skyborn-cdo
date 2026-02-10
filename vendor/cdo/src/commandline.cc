/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <vector>
#include <cstring>

static int gargc = 0;
static char **gargv;

static std::vector<char> CDO_CommandLine;

static void
initCommandLine(void)
{
  size_t maxlen = 1;
  for (int iarg = 0; iarg < gargc; iarg++) maxlen += std::strlen(gargv[iarg]) + 1;

  CDO_CommandLine.resize(maxlen);

  char *pargv;
  size_t offset = 0;
  for (int iarg = 0; iarg < gargc; iarg++)
    {
      if (iarg == 0)
        {
          pargv = strrchr(gargv[0], '/');
          if (pargv == 0)
            pargv = gargv[0];
          else
            pargv++;
        }
      else
        pargv = gargv[iarg];

      auto len = std::strlen(pargv);
      if (offset + len + 1 > maxlen) break;
      std::memcpy(&CDO_CommandLine[offset], pargv, len);
      offset += len;
      CDO_CommandLine[offset] = ' ';
      offset++;
    }

  CDO_CommandLine[offset - 1] = '\0';
}

namespace cdo
{

const char *
command_line(void)
{
  static bool init = false;

  if (!init)
    {
      initCommandLine();
      init = true;
    }

  return CDO_CommandLine.data();
}

void
set_command_line(int argc, char **argv)
{
  gargc = argc;
  gargv = argv;
}

};  // namespace cdo
