/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "cdo_omp.h"
#include "field_functions.h"
#include "pmlist.h"
#include "fill_1d.h"

class Vertfillmiss : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertfillmiss",
    .operators = { { "vertfillmiss", VertfillmissHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertfillmiss> registration = RegisterEntry<Vertfillmiss>(module);

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  FieldVector2D varsData;

  int calendar{};

  FillMethod method{ FillMethod::Nearest };
  int limit{ 0 };
  int maxGaps{ 0 };

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

    field2D_init(varsData, varList1);
  }

  void
  fillmiss(int varID)
  {
    auto const &var = varList1.vars[varID];
    auto fieldMemType = var.memType;
    auto gridsize = var.gridsize;
    auto numLevels = var.nlevels;
    auto missval = var.missval;

    Varray2D<double> dataValues2D(Threading::ompNumMaxThreads);
    for (auto &dataValues : dataValues2D) dataValues.resize(numLevels);

    Varray<double> levelValues(numLevels);
    zaxisInqLevels(var.zaxisID, levelValues.data());

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (size_t i = 0; i < gridsize; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &dataValues = dataValues2D[ompthID];

      if (fieldMemType == MemType::Float)
        for (int levelID = 0; levelID < numLevels; ++levelID) dataValues[levelID] = varsData[varID][levelID].vec_f[i];
      else
        for (int levelID = 0; levelID < numLevels; ++levelID) dataValues[levelID] = varsData[varID][levelID].vec_d[i];

      // clang-format off
      if      (method == FillMethod::Nearest)  fill_1d_nearest(numLevels, levelValues, dataValues, missval, limit, maxGaps);
      else if (method == FillMethod::Linear)   fill_1d_linear(numLevels, levelValues, dataValues, missval, limit, maxGaps);
      else if (method == FillMethod::Forward)  fill_1d_forward(numLevels, dataValues, missval, limit, maxGaps);
      else if (method == FillMethod::Backward) fill_1d_backward(numLevels, dataValues, missval, limit, maxGaps);
      // clang-format on

      if (fieldMemType == MemType::Float)
        for (int levelID = 0; levelID < numLevels; ++levelID) varsData[varID][levelID].vec_f[i] = dataValues[levelID];
      else
        for (int levelID = 0; levelID < numLevels; ++levelID) varsData[varID][levelID].vec_d[i] = dataValues[levelID];
    }
  }

  void
  run() override
  {
    auto numVars = varList1.numVars();
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
        auto &field = varsData[varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        int numLevels = varList1.vars[varID].nlevels;
        if (numLevels > 1)
        {
          size_t numMissVals = 0;
          for (int levelID = 0; levelID < numLevels; ++levelID) { numMissVals += varsData[varID][levelID].numMissVals; }
          if (numMissVals > 0) fillmiss(varID);
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
        {
          auto &field = varsData[varID][levelID];
          if (field.hasData())
          {
            cdo_def_field(streamID2, varID, levelID);
            field_num_mv(field);
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
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
