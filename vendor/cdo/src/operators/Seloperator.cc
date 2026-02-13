/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "cdi_lockedIO.h"

class Seloperator : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seloperator",
    .operators = { { "seloperator" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Seloperator> registration = RegisterEntry<Seloperator>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

  bool dataIsUnchanged{};

public:
  void
  init() override
  {
    bool selfound = false;

    dataIsUnchanged = data_is_unchanged();

    operator_input_arg("code, ltype, level");

    auto scode = parameter_to_int(cdo_operator_argv(0));
    auto sltype = parameter_to_int(cdo_operator_argv(1));

    auto slevel = (cdo_operator_argc() == 3) ? parameter_to_double(cdo_operator_argv(2)) : 0.0;

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto code = vlistInqVarCode(vlistID1, varID);
      auto zaxisID = vlistInqVarZaxis(vlistID1, varID);
      auto nlevels = zaxisInqSize(zaxisID);
      auto ltype = zaxis_to_ltype(zaxisID);

      for (int levID = 0; levID < nlevels; levID++)
      {
        auto level = cdo_zaxis_inq_level(zaxisID, levID);
        auto sellevel = (cdo_operator_argc() == 3) ? is_equal(level, slevel) : true;
        auto selcode = (scode == -1 || scode == code);
        auto selltype = (sltype == -1 || sltype == ltype);

        if (selcode && selltype && sellevel)
        {
          vlistDefFlag(vlistID1, varID, levID, true);
          selfound = true;
        }
      }
    }

    if (!selfound) cdo_warning("Code %d, ltype %d, level %g not found!", scode, sltype, slevel);

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (vlistInqFlag(vlistID1, varID, levelID) == true)
        {
          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
          cdo_def_field(streamID2, varID2, levelID2);

          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            field.init(varList1.vars[varID]);
            cdo_read_field(streamID1, field);
            cdo_write_field(streamID2, field);
          }
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
