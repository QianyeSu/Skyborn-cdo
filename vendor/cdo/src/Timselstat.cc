/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timselstat    timselrange        Time selection range
      Timselstat    timselmin          Time selection minimum
      Timselstat    timselmax          Time selection maximum
      Timselstat    timselsum          Time selection sum
      Timselstat    timselmean         Time selection mean
      Timselstat    timselavg          Time selection average
      Timselstat    timselvar          Time selection variance
      Timselstat    timselvar1         Time selection variance [Normalize by (n-1)]
      Timselstat    timselstd          Time selection standard deviation
      Timselstat    timselstd1         Time selection standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_stepstat.h"
#include "process_int.h"
#include "param_conversion.h"
#include "datetime.h"
#include "progress.h"
#include "field_functions.h"

class Timselstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timselstat",
    .operators = { { "timselrange", FieldFunc_Range, 0, TimselstatHelp },
                   { "timselmin", FieldFunc_Min, 0, TimselstatHelp },
                   { "timselmax", FieldFunc_Max, 0, TimselstatHelp },
                   { "timselsum", FieldFunc_Sum, 0, TimselstatHelp },
                   { "timselmean", FieldFunc_Mean, 0, TimselstatHelp },
                   { "timselavg", FieldFunc_Avg, 0, TimselstatHelp },
                   { "timselvar", FieldFunc_Var, 0, TimselstatHelp },
                   { "timselvar1", FieldFunc_Var1, 0, TimselstatHelp },
                   { "timselstd", FieldFunc_Std, 0, TimselstatHelp },
                   { "timselstd1", FieldFunc_Std1, 0, TimselstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timselstat> registration = RegisterEntry<Timselstat>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int noffset{};
  int ndates{};
  int nskip{};

  int maxFields{};

  cdo::StepStat2D stepStat{};

  std::vector<FieldInfo> fieldInfoList{};
  DateTimeList dtlist{};
  VarList varList1{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    stepStat.init(operfunc);

    operator_input_arg("nsets <noffset <nskip>>");

    auto nargc = cdo_operator_argc();
    ndates = parameter_to_int(cdo_operator_argv(0));
    noffset = (nargc > 1) ? parameter_to_int(cdo_operator_argv(1)) : 0;
    nskip = (nargc > 2) ? parameter_to_int(cdo_operator_argv(2)) : 0;

    if (Options::cdoVerbose) cdo_print("nsets=%d, noffset=%d, nskip=%d", ndates, noffset, nskip);
    if (ndates < 1) cdo_abort("nsets must be greater than 0!");
    if (noffset < 0) cdo_abort("noffset must be greater equal 0!");
    if (nskip < -1) cdo_abort("nskip must be greater equal -1!");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!stepStat.lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxFields = varList1.maxFields();
    fieldInfoList = std::vector<FieldInfo>(maxFields);

    dtlist.set_stat(TimeStat::MEAN);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    stepStat.alloc(varList1, VARS_MEMTYPE);
  }

  void
  skip_noffset_steps(int &tsID)
  {
    for (tsID = 0; tsID < noffset; ++tsID)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
      }
    }

    if (tsID < noffset) { cdo_abort("noffset is larger than number of timesteps!"); }
  }

  void
  skip_nskip_steps(int &tsID, int &numFields)
  {
    for (int i = 0; i < nskip; ++i)
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;
      tsID++;
    }
  }

  void
  run() override
  {
    FieldVector2D lastData;  // used if nskip=-1
    Field field;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    if (nskip == -1)
    {
      int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
      field2D_init(lastData, varList1, FIELD_VEC | VARS_MEMTYPE);
    }

    int tsID;
    skip_noffset_steps(tsID);

    int otsID = 0;

    while (true)
    {
      int numSetsStart{ 0 };
      if (nskip == -1 && otsID > 0)
      {
        numSetsStart = 1;
        for (int varID = 0; varID < varList1.numVars(); ++varID)
        {
          for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
          {
            stepStat.add_field(lastData[varID][levelID], varID, levelID, 0);
          }
        }
      }

      int numFields{ 0 };
      int numSets{ 0 };
      for (numSets = numSetsStart; numSets < ndates; numSets++)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

        dtlist.taxis_inq_timestep(taxisID1, numSets);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          stepStat.add_field(field, varID, levelID, numSets);
          if (nskip == -1 && numSets == (ndates - 1)) { field_copy(field, lastData[varID][levelID]); }
        }

        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;
      if (numFields == 0 && numSets == 1 && nskip == -1) break;

      cdo::fields_process(fieldInfoList, varList1, stepStat, numSets);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo::write_out_stream(streamID2, fieldInfoList, varList1, stepStat, otsID);

      if (numSets < ndates)
      {
        cdo_warning("Last output step %d contains only %d of %d input time step%s!", otsID + 1, numSets, ndates,
                    (numSets == 1) ? "" : "s");
        break;
      }

      if (numFields == 0) break;
      otsID++;

      skip_nskip_steps(tsID, numFields);

      if (numFields == 0) break;
    }

    if (otsID == 0) cdo_abort("No output step found!");
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
