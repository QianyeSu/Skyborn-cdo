/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <vector>
#include <string>
#include <cstring>
#include <time.h>

#include <cdi.h>

#include "cdo_options.h"
#include "cdi_uuid.h"

static int numVlist;
static constexpr int maxVlist = 256;
static int vlistHistory[maxVlist];

static char strtime[32];
static char datetimestr[32];

static void
init_strtime()
{
  auto tp = time(nullptr);
  if (tp != -1)
    {
      auto const *ltime = localtime(&tp);
      (void) strftime(strtime, sizeof(strtime), "%a %b %d %H:%M:%S %Y: ", ltime);
      (void) strftime(datetimestr, sizeof(datetimestr), "%Y-%m-%dT%H:%M:%SZ", ltime);
    }
}

static const char *
get_strtimeptr()
{
  if (std::strlen(strtime) == 0) init_strtime();

  return strtime;
}

std::string
cdo_inq_history(int vlistID)
{
  static const char *historyAttrName = "history";

  std::vector<char> ghistory;
  auto atttype = cdiInqAttType(vlistID, CDI_GLOBAL, historyAttrName);
  if (atttype == CDI_DATATYPE_TXT)
    {
      auto ghistorysize = cdiInqAttLen(vlistID, CDI_GLOBAL, historyAttrName);
      if (ghistorysize > 0)
        {
          ghistory.resize(ghistorysize + 1);
          cdiInqAttTxt(vlistID, CDI_GLOBAL, historyAttrName, ghistorysize, ghistory.data());
          ghistory[ghistorysize] = 0;
        }
    }

  return ghistory.data();
}

void
cdo_append_history(int vlistID, const char *histstring)
{
  for (int i = 0; i < numVlist; ++i)
    if (vlistHistory[i] == vlistID) return;

  if (numVlist < maxVlist) vlistHistory[numVlist++] = vlistID;

  static const char *historyAttrName = "history";

  if (Options::CDO_Reset_History) cdiDelAtt(vlistID, CDI_GLOBAL, historyAttrName);

  if (!Options::CDO_Append_History) return;

  std::vector<char> ghistory;
  const auto atttype = cdiInqAttType(vlistID, CDI_GLOBAL, historyAttrName);
  if (atttype == CDI_DATATYPE_TXT)
    {
      auto ghistorysize = cdiInqAttLen(vlistID, CDI_GLOBAL, historyAttrName);
      if (ghistorysize < 0) ghistorysize = 0;
      if (ghistorysize > 0)
        {
          ghistory.resize(ghistorysize + 1);
          cdiInqAttTxt(vlistID, CDI_GLOBAL, historyAttrName, ghistorysize, ghistory.data());
          ghistory[ghistorysize] = 0;
        }
    }
  else if (atttype != -1) { return; }

  std::string history = get_strtimeptr();
  history += histstring;

  if (!ghistory.empty())
    {
      history += "\n";
      history += ghistory.data();
    }

  cdiDefAttTxt(vlistID, CDI_GLOBAL, historyAttrName, history.size(), history.c_str());
}

void
cdo_def_creation_date(int vlistID)
{
  if (std::strlen(datetimestr) == 0) init_strtime();
  cdiDefAttTxt(vlistID, CDI_GLOBAL, "creation_date", (int) std::strlen(datetimestr), datetimestr);
}

static void
get_uuid(char *uuidstr)
{
  unsigned char uuid[CDI_UUID_SIZE];
  cdiCreateUUID(uuid);
  cdiUUID2Str(uuid, uuidstr);
}

void
cdo_def_tracking_id(int vlistID, const char *uuid_attribute)
{
  char uuidstr[uuidNumHexChars + 1] = { 0 };
  get_uuid(uuidstr);
  cdiDefAttTxt(vlistID, CDI_GLOBAL, uuid_attribute, uuidNumHexChars, uuidstr);
}
