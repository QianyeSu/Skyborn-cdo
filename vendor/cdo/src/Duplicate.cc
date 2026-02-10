/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "field_functions.h"

class Duplicate : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Duplicate",
    .operators = { { "duplicate", DuplicateHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Duplicate> registration = RegisterEntry<Duplicate>(module);

private:
  FieldVector3D varsData;
  std::vector<CdiDateTime> vDateTimes;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };

  int numVars{};
  int numDuplicates{};

  VarList varList1;

public:
  void
  init() override
  {
    if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

    numDuplicates = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 2;
    if (Options::cdoVerbose) cdo_print("ndup = %d", numDuplicates);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    auto numSteps = varList1.numSteps();

    if (numSteps == 1 && varList1.numVaryingVars() == 0) numSteps = 0;

    if (numSteps == 0)
    {
      for (int varID = 0; varID < numVars; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
      numSteps = 1;
    }

    auto numStepsOut = (numSteps > 0) ? numDuplicates * numSteps : -1;
    vlistDefNtsteps(vlistID2, numStepsOut);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= vDateTimes.size()) vDateTimes.resize(vDateTimes.size() + NALLOC_INC);
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      field2D_init(varsData[tsID], varList1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field = varsData[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }

      tsID++;
    }

    auto nts = tsID;
    Field fieldOut;

    for (int idup = 0; idup < numDuplicates; idup++)
    {
      for (tsID = 0; tsID < nts; ++tsID)
      {
        taxisDefVdatetime(taxisID2, vDateTimes[tsID]);
        cdo_def_timestep(streamID2, idup * nts + tsID);

        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var1 = varList1.vars[varID];
          fieldOut.init(var1);
          for (int levelID = 0; levelID < var1.nlevels; ++levelID)
          {
            auto const &field = varsData[tsID][varID][levelID];
            if (field.hasData())
            {
              field_copy(field, fieldOut);
              cdo_def_field(streamID2, varID, levelID);
              cdo_write_field(streamID2, fieldOut);
            }
          }
        }
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
