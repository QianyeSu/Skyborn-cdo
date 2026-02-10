/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pardup     pardup          Duplicate parameters
      Pardup     parmul          Multiply parameters
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"

class Pardup : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Pardup",
    .operators = { { "pardup" }, { "parmul" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Pardup> registration = RegisterEntry<Pardup>(module);

private:
  int PARDUP{}, PARMUL{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int numDuplicates = 0;

  VarList varList1{};

public:
  void
  init() override
  {
    PARDUP = module.get_id("pardup");
    PARMUL = module.get_id("parmul");

    auto operatorID = cdo_operator_id();

    if (operatorID == PARDUP) { numDuplicates = 2; }
    else if (operatorID == PARMUL)
    {
      operator_input_arg("number of multiply");
      numDuplicates = parameter_to_int(cdo_operator_argv(0));
    }
    else
      cdo_abort("operator not implemented!");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);
    auto numVars = varList1.numVars();

    for (int i = 1; i < numDuplicates; ++i)
    {
      vlistCat(vlistID2, vlistID1);
      for (int varID = 0; varID < numVars; ++varID)
        vlistDefVarParam(vlistID2, varID + numVars * i, cdiEncodeParam(-(varID + numVars * i + 1), 255, 255));
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto numVars = varList1.numVars();

    Varray<double> array;
    Varray2D<double> vardata;
    std::vector<std::vector<size_t>> varnumMissVals;

    auto gridsizeMax = varList1.gridsizeMax();
    array.resize(gridsizeMax);
    vardata.resize(numVars);
    varnumMissVals.resize(numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = varList1.vars[varID].gridsize;
      auto nlevels = varList1.vars[varID].nlevels;
      vardata[varID].resize(gridsize * nlevels);
      varnumMissVals[varID].resize(nlevels);
    }

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

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

        fieldInfoList[fieldID].set(varID, levelID);

        auto gridsize = varList1.vars[varID].gridsize;
        auto offset = gridsize * levelID;
        auto single = &vardata[varID][offset];

        size_t numMissVals;
        cdo_read_field(streamID1, single, &numMissVals);
        varnumMissVals[varID][levelID] = numMissVals;
      }

      for (int i = 0; i < numDuplicates; ++i)
        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();

          auto varID2 = varID + i * numVars;

          auto gridsize = varList1.vars[varID].gridsize;
          auto offset = gridsize * levelID;
          auto single = &vardata[varID][offset];
          auto numMissVals = varnumMissVals[varID][levelID];

          array_copy(gridsize, single, array.data());
          cdo_def_field(streamID2, varID2, levelID);
          cdo_write_field(streamID2, array.data(), numMissVals);
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
