/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/

#include <time.h>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_season.h"
#include "util_files.h"
#include "util_string.h"

constexpr int MaxStreams = 32;

static struct tm
datetime_to_tm(CdiDateTime vDateTime)
{
  struct tm stime;
  memset(&stime, 0, sizeof(struct tm));

  stime.tm_sec = vDateTime.time.second;
  stime.tm_min = vDateTime.time.minute;
  stime.tm_hour = vDateTime.time.hour;
  stime.tm_mday = vDateTime.date.day;
  stime.tm_mon = vDateTime.date.month - 1;
  stime.tm_year = vDateTime.date.year - 1900;

  return stime;
}

class Splittime : public Process
{
  enum
  {
    func_time,
    func_date
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Splittime",
    .operators = { { "splithour", func_time, 10000, SplittimeHelp },
                   { "splitday", func_date, 1, SplittimeHelp },
                   { "splitmon", func_date, 100, SplittimeHelp },
                   { "splitseas", func_date, 100, SplittimeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Splittime> registration = RegisterEntry<Splittime>();

private:
  int SPLITMON{}, SPLITSEAS{};
  CdoStreamID streamID1{};
  CdoStreamID streamIDs[MaxStreams]{};
  int tsIDs[MaxStreams]{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operfunc{};
  int operintval{};
  int operatorID{};

  std::string fileSuffix{};
  const char *format = nullptr;

  VarList varList1{};

  bool dataIsUnchanged{};

  void
  init_fields(FieldVector2D &fields)
  {
    auto numVars = varList1.numVars();
    fields.resize(numVars);
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      if (var1.isConstant)
      {
        fields[varID].resize(var1.nlevels);
        for (auto &field : fields[varID]) { field.init(var1); }
      }
    }
  }

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    SPLITMON = module.get_id("splitmon");
    SPLITSEAS = module.get_id("splitseas");

    operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    operintval = cdo_operator_f2(operatorID);

    if (operatorID == SPLITMON && cdo_operator_argc() == 1)
      format = cdo_operator_argv(0).c_str();
    else
      operator_check_argc(0);

    std::ranges::fill(streamIDs, CDO_STREAM_UNDEF);
    std::ranges::fill(tsIDs, 0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector2D fields{};
    auto haveConstVars = (varList1.numConstVars() > 0);
    if (haveConstVars) { init_fields(fields); }

    auto seasonNames = get_season_name();
    Field field;
    int stepIndex = 0;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (operfunc == func_date)
      {
        stepIndex = (cdiDate_get(vDateTime.date) / operintval) % 100;
        if (stepIndex < 0) stepIndex = -stepIndex;

        if (operatorID == SPLITSEAS) stepIndex = month_to_season(stepIndex);
      }
      else if (operfunc == func_time) { stepIndex = (cdiTime_get(vDateTime.time) / operintval) % 100; }

      if (stepIndex < 0 || stepIndex >= MaxStreams) cdo_abort("Step index out of range!");

      auto &streamID2 = streamIDs[stepIndex];
      if (streamID2 == CDO_STREAM_UNDEF)
      {
        auto fileName = cdo_get_obase();
        if (operatorID == SPLITSEAS) { fileName += string_format("%3s", seasonNames[stepIndex]); }
        else
        {
          char oformat[32];
          std::strcpy(oformat, "%02d");

          if (operatorID == SPLITMON && format)
          {
            char sbuf[32];
            auto stime = datetime_to_tm(vDateTime);
            auto slen = strftime(sbuf, sizeof(sbuf), format, &stime);
            if (slen) std::strcpy(oformat, sbuf);
          }

          fileName += string_format(oformat, stepIndex);
        }

        if (fileSuffix.size() > 0) fileName += fileSuffix;

        streamID2 = open_write(fileName);
        cdo_def_vlist(streamID2, vlistID2);
      }

      cdo_def_timestep(streamID2, tsIDs[stepIndex]);

      if (tsID > 0 && tsIDs[stepIndex] == 0 && haveConstVars)
      {
        for (int varID = 0; varID < varList1.numVars(); ++varID)
        {
          auto const &var = varList1.vars[varID];
          if (var.isConstant)
          {
            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              cdo_def_field(streamID2, varID, levelID);
              cdo_write_field(streamID2, fields[varID][levelID]);
            }
          }
        }
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        if (dataIsUnchanged && !(tsID == 0 && haveConstVars)) { cdo_copy_field(streamID1, streamID2); }
        else
        {
          auto const &var1 = varList1.vars[varID];
          field.init(var1);
          cdo_read_field(streamID1, field);
          cdo_write_field(streamID2, field);

          if (tsID == 0 && haveConstVars && var1.isConstant) { field_copy(field, fields[varID][levelID]); }
        }
      }

      tsIDs[stepIndex]++;
      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);

    for (auto &streamID : streamIDs)
    {
      if (streamID != CDO_STREAM_UNDEF) cdo_stream_close(streamID);
    }

    vlistDestroy(vlistID2);
  }
};
