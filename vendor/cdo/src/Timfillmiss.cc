/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "datetime.h"
#include "cdo_omp.h"
#include "field_functions.h"
#include "pmlist.h"
#include "fill_1d.h"

static double
julianDate_to_double(int calendar, CdiDateTime const &dateTime1, CdiDateTime const &datetime0)
{
  return julianDate_to_seconds(julianDate_sub(julianDate_encode(calendar, dateTime1), julianDate_encode(calendar, datetime0)))
         / 86400.0;
}

class Timfillmiss : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timfillmiss",
    .operators = { { "timfillmiss", TimfillmissHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timfillmiss> registration = RegisterEntry<Timfillmiss>(module);

private:
  DateTimeList dtlist{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  FieldVector3D varsData{};

  int calendar{};
  int numVars{};

  FillMethod method{ FillMethod::Nearest };
  int limit{ 0 };
  int maxGaps{ 0 };

  Varray2D<double> dataValues2D{};

  Varray<double> timeValues{};
  int numSteps{};

  void
  get_parameter()
  {
    auto numArgs = cdo_operator_argc();
    if (numArgs)
    {
      auto const &argList = cdo_get_oper_argv();

      KVList kvlist;
      kvlist.name = cdo_module_name();
      if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      for (auto const &kv : kvlist)
      {
        auto const &key = kv.key;
        if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
        if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
        auto const &value = kv.values[0];

        // clang-format off
        if      (key == "method")   method = convert<FillMethod>(value);
        else if (key == "limit")    limit = parameter_to_int(value);
        else if (key == "max_gaps") maxGaps = parameter_to_int(value);
        else cdo_abort("Invalid parameter key >%s<!", key);
        // clang-format on
      }
    }
  }

public:
  void
  init() override
  {
    get_parameter();
    limit = std::max(limit, 0);
    maxGaps = std::max(maxGaps, 0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_write(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    cdo_def_vlist(streamID2, vlistID2);

    calendar = taxisInqCalendar(taxisID1);

    varList1 = VarList(vlistID1);
    numVars = varList1.numVars();
  }

  void
  step(int varID)
  {
    auto const &var = varList1.vars[varID];
    auto fieldMemType = var.memType;
    auto gridsize = var.gridsize;
    auto missval = var.missval;
    for (int levelID = 0; levelID < var.nlevels; ++levelID)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (size_t i = 0; i < gridsize; ++i)
      {
        auto ompthID = cdo_omp_get_thread_num();
        auto &dataValues = dataValues2D[ompthID];

        if (fieldMemType == MemType::Float)
          for (int t = 0; t < numSteps; ++t) dataValues[t] = varsData[t][varID][levelID].vec_f[i];
        else
          for (int t = 0; t < numSteps; ++t) dataValues[t] = varsData[t][varID][levelID].vec_d[i];

        // clang-format off
        if      (method == FillMethod::Nearest)  fill_1d_nearest(numSteps, timeValues, dataValues, missval, limit, maxGaps);
        else if (method == FillMethod::Linear)   fill_1d_linear(numSteps, timeValues, dataValues, missval, limit, maxGaps);
        else if (method == FillMethod::Forward)  fill_1d_forward(numSteps, dataValues, missval, limit, maxGaps);
        else if (method == FillMethod::Backward) fill_1d_backward(numSteps, dataValues, missval, limit, maxGaps);
        // clang-format on

        if (fieldMemType == MemType::Float)
          for (int t = 0; t < numSteps; ++t) varsData[t][varID][levelID].vec_f[i] = dataValues[t];
        else
          for (int t = 0; t < numSteps; ++t) varsData[t][varID][levelID].vec_d[i] = dataValues[t];
      }
    }
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
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

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

    numSteps = tsID;
    if (numSteps <= 1) cdo_abort("Number of time steps %d!", numSteps);

    timeValues = Varray<double>(numSteps);
    for (tsID = 0; tsID < numSteps; ++tsID)
    {
      timeValues[tsID] = julianDate_to_double(calendar, dtlist.vDateTime(tsID), dtlist.vDateTime(0));
    }

    dataValues2D = Varray2D<double>(Threading::ompNumMaxThreads);
    for (auto &dataValues : dataValues2D) dataValues.resize(numSteps);

    for (int varID = 0; varID < numVars; ++varID) { step(varID); }

    for (tsID = 0; tsID < numSteps; ++tsID)
    {
      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
        {
          auto &field = varsData[tsID][varID][levelID];
          if (field.hasData())
          {
            cdo_def_field(streamID2, varID, levelID);
            field_num_mv(field);
            cdo_write_field(streamID2, field);
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
