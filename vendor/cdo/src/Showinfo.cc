/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Showinfo   showparam       Show parameters
      Showinfo   showcode        Show code numbers
      Showinfo   showname        Show variable names
      Showinfo   showstdname     Show variable standard names
      Showinfo   showlevel       Show levels
      Showinfo   showyear        Show years
      Showinfo   showmon         Show months
      Showinfo   showdate        Show dates
      Showinfo   showtime        Show timesteps
      Showinfo   showltype       Show level types
      Showinfo   showformat      Show file format
*/

#include <cdi.h>

#include "process_int.h"
#include "printinfo.h"
#include "cdo_history.h"
#include "cdo_zaxis.h"

static void
print_newline(int nout, int maxOut)
{
  if (!Options::silentMode && !(nout % maxOut)) fprintf(stdout, "\n");
}

static void
print_newline_if_missing(int nout, int maxOut)
{
  if (Options::silentMode || (nout % maxOut)) fprintf(stdout, "\n");
}

static void
show_year(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto numSteps = vlistNtsteps(vlistID);
  if (numSteps == 0) return;

  constexpr int maxOut = 20;
  int nout = 0;
  int year0 = 0;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(taxisID);
    int year = vDateTime.date.year;

    if (tsID == 0 || year0 != year)
    {
      nout++;
      year0 = year;
      fprintf(stdout, " %4d", year0);
      print_newline(nout, maxOut);
    }

    tsID++;
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_mon(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto numSteps = vlistNtsteps(vlistID);
  if (numSteps == 0) return;

  constexpr int maxOut = 36;
  int nout = 0;
  int month0 = 0;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(taxisID);
    int month = vDateTime.date.month;

    if (tsID == 0 || month0 != month)
    {
      nout++;
      month0 = month;
      fprintf(stdout, " %2d", month0);
      print_newline(nout, maxOut);
    }

    tsID++;
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_date(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto numSteps = vlistNtsteps(vlistID);
  if (numSteps == 0) return;

  constexpr int maxOut = 12;
  int nout = 0;
  int64_t date0 = 0;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(taxisID);
    auto vdate = cdiDate_get(vDateTime.date);

    if (tsID == 0 || date0 != vdate)
    {
      nout++;
      date0 = vdate;
      fprintf(stdout, " %s", date_to_string(vDateTime.date).c_str());
      print_newline(nout, maxOut);
    }

    tsID++;
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_time(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto numSteps = vlistNtsteps(vlistID);
  if (numSteps == 0) return;

  constexpr int maxOut = 12;
  int nout = 0;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    auto vDateTime = taxisInqVdatetime(taxisID);
    nout++;
    fprintf(stdout, " %s", time_to_string(vDateTime.time).c_str());
    print_newline(nout, maxOut);

    tsID++;
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_timestamp(CdoStreamID streamID)
{
  auto vlistID = cdo_stream_inq_vlist(streamID);
  auto taxisID = vlistInqTaxis(vlistID);
  auto numSteps = vlistNtsteps(vlistID);
  if (numSteps == 0) return;

  constexpr int maxOut = 4;
  int nout = 0;

  int tsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    nout++;
    fprintf(stdout, " %s", datetime_to_string(taxisInqVdatetime(taxisID)).c_str());
    print_newline(nout, maxOut);

    tsID++;
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_code(VarList const &varList)
{
  constexpr int maxOut = 25;
  int nout = 0;

  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    nout++;
    auto const &var = varList.vars[varID];
    fprintf(stdout, " %d", var.code);
    print_newline(nout, maxOut);
  }
  if (Options::silentMode || (nout % maxOut)) fprintf(stdout, "\n");
}

static void
show_grid(VarList const &varList)
{
  fprintf(stdout, "# param nr | grid nr | z-axis nr:   /* Use in combination with operatores: griddes and zaxisdes */\n");
  auto vlistID = varList.vlistID;
  for (auto const &var : varList.vars)
  {
    fprintf(stdout, "      %3d     %3d      %3d\n", var.code, vlistGridIndex(vlistID, var.gridID) + 1,
            vlistZaxisIndex(vlistID, var.zaxisID) + 1);
  }
}

static void
show_unit(VarList const &varList)
{
  constexpr int maxOut = 10;
  int nout = 0;

  for (auto const &var : varList.vars)
  {
    nout++;
    if (var.units.size()) fprintf(stdout, " %s", var.units.c_str());
    print_newline(nout, maxOut);
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_param(VarList const &varList)
{
  constexpr int maxOut = 10;
  int nout = 0;

  char paramstr[32];
  for (auto const &var : varList.vars)
  {
    nout++;
    cdiParamToString(var.param, paramstr, sizeof(paramstr));

    fprintf(stdout, " %s", paramstr);
    print_newline(nout, maxOut);
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_name(VarList const &varList)
{
  int nout = 0;
  for (auto const &var : varList.vars)
  {
    nout++;
    if (nout > 1) fprintf(stdout, " ");
    fprintf(stdout, "%s", var.name.c_str());
  }
  fprintf(stdout, "\n");
}

static void
show_stdname(VarList const &varList)
{
  constexpr int maxOut = 1;
  int nout = 0;

  for (auto const &var : varList.vars)
  {
    nout++;
    fprintf(stdout, " %s", var.stdname.size() ? var.stdname.c_str() : "unknown");
    print_newline(nout, maxOut);
  }
  print_newline_if_missing(nout, maxOut);
}

static void
show_level(VarList const &varList)
{
  for (auto const &var : varList.vars)
  {
    for (int levelID = 0; levelID < var.nlevels; ++levelID) fprintf(stdout, " %.9g", cdo_zaxis_inq_level(var.zaxisID, levelID));
    fprintf(stdout, "\n");
  }
}

static void
show_ltype(int vlistID)
{
  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    auto ltype = zaxis_to_ltype(zaxisID);

    if (ltype != -1) fprintf(stdout, " %d", ltype);
  }
  fprintf(stdout, "\n");
}

static void
show_filter(int vlistID, VarList const &varList)
{
  char filterSpec[CDI_MAX_NAME];
  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto comptype = vlistInqVarCompType(vlistID, varID);
    if (comptype == CDI_COMPRESS_FILTER || comptype == CDI_COMPRESS_ZIP || comptype == CDI_COMPRESS_SZIP)
    {
      auto const &var = varList.vars[varID];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(vlistID, varID, CDI_KEY_FILTERSPEC_IN, filterSpec, &length);
      if (length > 0) fprintf(stdout, "%s=\"%s\"\n", var.name.c_str(), filterSpec);
    }
  }
}

static void
show_chunkspec(int vlistID, VarList const &varList)
{
  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var = varList.vars[varID];
    auto chunkSpecString = cdo::get_chunkspec_string(vlistID, varID);
    if (chunkSpecString.size() > 0) fprintf(stdout, "%s=\"%s\"\n", var.name.c_str(), chunkSpecString.c_str());
  }
}

static void
show_history(int vlistID)
{
  auto historyString = cdo_inq_history(vlistID);
  fprintf(stdout, "%s\n", historyString.c_str());
}

class Showinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Showinfo",
    .operators = { { "showyear", ShowinfoHelp },
                   { "showmon", ShowinfoHelp },
                   { "showdate", ShowinfoHelp },
                   { "showtime", ShowinfoHelp },
                   { "showtimestamp", ShowinfoHelp },
                   { "showcode", ShowinfoHelp },
                   { "showunit", ShowinfoHelp },
                   { "showparam", ShowinfoHelp },
                   { "showname", ShowinfoHelp },
                   { "showstdname", ShowinfoHelp },
                   { "showlevel", ShowinfoHelp },
                   { "showltype", ShowinfoHelp },
                   { "showformat", ShowinfoHelp },
                   { "showgrid", ShowinfoHelp },
                   { "showchunkspec", ShowinfoHelp },
                   { "showhistory" },
                   { "showfilter", ShowinfoHelp } },
    .aliases = { { "showvar", "showname" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Showinfo> registration = RegisterEntry<Showinfo>(module);

public:
  void
  init() override
  {
    cdiDefGlobal("COPY_CHUNKSPEC", true);

    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    auto SHOWYEAR = module.get_id("showyear");
    auto SHOWMON = module.get_id("showmon");
    auto SHOWDATE = module.get_id("showdate");
    auto SHOWTIME = module.get_id("showtime");
    auto SHOWTIMESTAMP = module.get_id("showtimestamp");
    auto SHOWCODE = module.get_id("showcode");
    auto SHOWUNIT = module.get_id("showunit");
    auto SHOWPARAM = module.get_id("showparam");
    auto SHOWNAME = module.get_id("showname");
    auto SHOWSTDNAME = module.get_id("showstdname");
    auto SHOWLEVEL = module.get_id("showlevel");
    auto SHOWLTYPE = module.get_id("showltype");
    auto SHOWFORMAT = module.get_id("showformat");
    auto SHOWGRID = module.get_id("showgrid");
    auto SHOWCHUNKSPEC = module.get_id("showchunkspec");
    auto SHOWFILTER = module.get_id("showfilter");
    auto SHOWHISTORY = module.get_id("showhistory");

    auto operatorID = cdo_operator_id();

    operator_check_argc(0);

    auto streamID = cdo_open_read(0);
    auto vlistID = cdo_stream_inq_vlist(streamID);
    VarList varList(vlistID);

    // clang-format off
    if      (operatorID == SHOWYEAR)      show_year(streamID);
    else if (operatorID == SHOWMON)       show_mon(streamID);
    else if (operatorID == SHOWDATE)      show_date(streamID);
    else if (operatorID == SHOWTIME)      show_time(streamID);
    else if (operatorID == SHOWTIMESTAMP) show_timestamp(streamID);
    else if (operatorID == SHOWCODE)      show_code(varList);
    else if (operatorID == SHOWGRID)      show_grid(varList);
    else if (operatorID == SHOWUNIT)      show_unit(varList);
    else if (operatorID == SHOWPARAM)     show_param(varList);
    else if (operatorID == SHOWNAME)      show_name(varList);
    else if (operatorID == SHOWSTDNAME)   show_stdname(varList);
    else if (operatorID == SHOWLEVEL)     show_level(varList);
    else if (operatorID == SHOWLTYPE)     show_ltype(vlistID);
    else if (operatorID == SHOWFORMAT)    print_filetype(streamID, vlistID);
    else if (operatorID == SHOWCHUNKSPEC) show_chunkspec(vlistID, varList);
    else if (operatorID == SHOWFILTER)    show_filter(vlistID, varList);
    else if (operatorID == SHOWHISTORY)   show_history(vlistID);
    // clang-format on

    cdo_stream_close(streamID);
  }

  void
  run() override
  {
  }

  void
  close() override
  {
  }
};
