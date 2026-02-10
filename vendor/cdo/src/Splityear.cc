/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Splityear  splityear       Split in years
     Splityear  splityearmon    Split in years and month
*/

#include <climits>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_files.h"
#include "util_string.h"

constexpr int MaxYears = 99999;

class Splityear : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Splityear",
    .operators = { { "splityear", SplittimeHelp }, { "splityearmon", SplittimeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Splityear> registration = RegisterEntry<Splityear>(module);

private:
  int SPLITYEAR{}, SPLITYEARMON{};
  int operatorID{};
  CdoStreamID streamID1{};
  Varray<int> cyear = Varray<int>(MaxYears, 0);
  std::string fileSuffix{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
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

    SPLITYEAR = module.get_id("splityear");
    SPLITYEARMON = module.get_id("splityearmon");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
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

    Field field;
    int ic = 0;
    int index1 = -INT_MAX;
    int index2;
    int year1 = -1, year2;
    int mon1 = -1, mon2;
    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      auto vDate = taxisInqVdatetime(taxisID1).date;
      int day;
      cdiDate_decode(vDate, &year2, &mon2, &day);

      if (operatorID == SPLITYEAR)
      {
        if (tsID == 0 || year1 != year2 || mon1 > mon2)
        {
          tsID2 = 0;

          ic = (year1 != year2) ? 0 : ic + 1;
          if (year2 >= 0 && year2 < MaxYears)
          {
            ic = cyear[year2];
            cyear[year2]++;
          }

          year1 = year2;

          if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

          auto fileName = cdo_get_obase() + string_format("%04d", year1);
          if (ic > 0) fileName += string_format("_%d", ic + 1);
          if (fileSuffix.size() > 0) fileName += fileSuffix;

          if (Options::cdoVerbose) cdo_print("create file %s", fileName);

          streamID2 = open_write(fileName);
          cdo_def_vlist(streamID2, vlistID2);
        }
        mon1 = mon2;
      }
      else if (operatorID == SPLITYEARMON)
      {
        index2 = year2 * 100 + mon2;

        if (tsID == 0 || index1 != index2)
        {
          tsID2 = 0;

          index1 = index2;

          if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

          auto fileName = cdo_get_obase() + string_format("%04d", index1);
          if (fileSuffix.size() > 0) fileName += fileSuffix;

          if (Options::cdoVerbose) cdo_print("create file %s", fileName);

          streamID2 = open_write(fileName);
          cdo_def_vlist(streamID2, vlistID2);
        }
      }

      cdo_def_timestep(streamID2, tsID2);

      if (tsID > 0 && tsID2 == 0 && haveConstVars)
      {
        auto numVars = varList1.numVars();
        for (int varID = 0; varID < numVars; ++varID)
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

      tsID2++;
      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
