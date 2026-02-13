/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splitdate   splitdate        Split into dates
*/

#include <cdi.h>

#include "process_int.h"
#include "util_files.h"
#include "util_string.h"

class Splitdate : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Splitdate",
    .operators = { { "splitdate", SplitdateHelp }, { "splitdatetime", SplitdateHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Splitdate> registration = RegisterEntry<Splitdate>();

private:
  int SPLITDATE{};
  CdoStreamID streamID1{};
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  bool splitDate{};
  bool dataIsUnchanged{};

  std::string fileSuffix{};

  VarList varList1{};

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

    SPLITDATE = module.get_id("splitdate");

    auto operatorID = cdo_operator_id();
    splitDate = (operatorID == SPLITDATE);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));
  }

  void
  run() override
  {
    FieldVector2D fields{};
    auto haveConstVars = (varList1.numConstVars() > 0);
    if (haveConstVars) { init_fields(fields); }

    int64_t vDate0 = -1;

    streamID2 = CDO_STREAM_UNDEF;
    Field field;
    int tsID2 = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (splitDate)
      {
        auto vDate = cdiDate_get(vDateTime.date);
        if (vDate != vDate0)
        {
          if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

          vDate0 = vDate;
          tsID2 = 0;

          auto const &date = vDateTime.date;
          auto fileName = cdo_get_obase() + string_format("%04d-%02d-%02d", date.year, date.month, date.day);
          if (fileSuffix.size() > 0) fileName += fileSuffix;

          streamID2 = open_write(fileName);
          cdo_def_vlist(streamID2, vlistID2);
        }
      }
      else
      {
        if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

        tsID2 = 0;

        auto const &date = vDateTime.date;
        auto const &time = vDateTime.time;
        auto fileName = cdo_get_obase()
                        + string_format("%04d-%02d-%02dT%02d:%02d:%02d", date.year, date.month, date.day, time.hour, time.minute,
                                        time.second);
        if (fileSuffix.size() > 0) fileName += fileSuffix;

        streamID2 = open_write(fileName);
        cdo_def_vlist(streamID2, vlistID2);
      }

      cdo_def_timestep(streamID2, tsID2);

      if (tsID > 0 && tsID2 == 0 && haveConstVars)
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

      tsID++;
      tsID2++;
    }
  }

  void
  close() override
  {
    if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
