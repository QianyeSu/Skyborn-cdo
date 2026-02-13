/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Splitsel   splitsel        Split time selection
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"
#include "util_string.h"

class Splitsel : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Splitsel",
    .operators = { { "splitsel", SplitselHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Splitsel> registration = RegisterEntry<Splitsel>();

private:
  CdoStreamID streamID1{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int numSets{};
  int numOffset{};
  int numSkip{};

  bool dataIsUnchanged{};

  std::string fileSuffix{};

  VarList varList1;

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

    // operator_input_arg("numSets <noffset <numSkip>>");

    auto nargc = cdo_operator_argc();
    if (nargc < 1) cdo_abort("Too few arguments! Need %d found %d.", 1, nargc);

    numSets = parameter_to_int(cdo_operator_argv(0));
    numOffset = (nargc > 1) ? parameter_to_int(cdo_operator_argv(1)) : 0;
    numSkip = (nargc > 2) ? parameter_to_int(cdo_operator_argv(2)) : 0;

    if (Options::cdoVerbose) cdo_print("numSets = %d, noffset = %d, numSkip = %d", numSets, numOffset, numSkip);

    if (numSets < 1) cdo_abort("numSets must be greater than 0!");
    if (numOffset < 0) cdo_abort("noffset must be greater or equal 0!");
    if (numSkip < 0)
    {
      if (-numSkip >= numSets) cdo_abort("Absolute value of negative numSkip must be less than numSets!");
      if (cdo_assert_files_only() == false) cdo_abort("Negative numSkip not allowed in combination with other operators!");
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    // int taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));
  }

  void
  run() override
  {
    FieldVector2D fields;
    auto haveConstVars = (varList1.numConstVars() > 0);
    if (haveConstVars) { init_fields(fields); }

    int index = 1;
    int numFields = 0;

    // offset
    int tsID = 0;
    for (; tsID < numOffset; ++tsID)
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0)
      {
        cdo_warning("noffset is larger than number of timesteps!");
        return;
      }

      if (tsID == 0 && haveConstVars)
        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);

          auto const &var = varList1.vars[varID];
          if (var.isConstant) { cdo_read_field(streamID1, fields[varID][levelID]); }
        }
    }

    Field field;
    while (true)
    {
      auto fileName = cdo_get_obase() + string_format("%06d", index);
      if (fileSuffix.size() > 0) fileName += fileSuffix;

      if (Options::cdoVerbose) cdo_print("create file %s", fileName);

      CdoStreamID streamID2 = CDO_STREAM_UNDEF;
      int tsID2 = 0;

      for (int set = 0; set < numSets; ++set)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        cdo_taxis_copy_timestep(taxisID2, taxisID1);

        if (streamID2 == CDO_STREAM_UNDEF)
        {
          streamID2 = open_write(fileName);
          cdo_def_vlist(streamID2, vlistID2);
        }

        cdo_def_timestep(streamID2, tsID2);

        if (tsID > 0 && tsID2 == 0 && haveConstVars)
        {
          int numVars = varList1.numVars();
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

        tsID++;
        tsID2++;
      }

      cdo_stream_close(streamID2);
      if (numFields == 0) break;

      if (cdo_stream_inq_timestep(streamID1, tsID) == 0) break;

      // skip
      for (int i = 0; i < numSkip; ++i)
        if (cdo_stream_inq_timestep(streamID1, tsID + i) == 0) break;

      tsID += numSkip;

      if (cdo_stream_inq_timestep(streamID1, tsID) == 0) break;

      index++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
