/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selyearidx    selyearidx         Select index of year
*/

#include <cdi.h>

#include "process_int.h"
#include "field_functions.h"

class Selyearidx : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selyearidx",
    // clang-format off
    .operators = { { "selyearidx", 0, 0, SelyearidxHelp },
                   { "seltimeidx", 1, 0, SeltimeidxHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Selyearidx> registration = RegisterEntry<Selyearidx>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;

  bool allSteps;

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    allSteps = cdo_operator_f1(cdo_operator_id());

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    taxisID1 = vlistInqTaxis(vlistID1);

    streamID2 = cdo_open_read(1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    taxisID2 = vlistInqTaxis(vlistID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);
    for (auto &var : varList2.vars) var.memType = MemType::Double;

    auto vlistID3 = vlistDuplicate(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Field field;

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    FieldVector2D varsData1, varsData2;
    field2D_init(varsData1, varList1, FIELD_VEC);
    field2D_init(varsData2, varList2, FIELD_VEC);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto missval = varList2.vars[varID].missval;
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        for (size_t i = 0; i < var.gridsize; ++i) varsData2[varID][levelID].vec_d[i] = missval;
      }
    }

    int tsID = 0;
    int tsID2 = 0;
    int tsID3 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto year1 = taxisInqVdatetime(taxisID1).date.year;

      auto useExactDate = (numFields == 1 && varList1.vars[0].gridsize == 1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field1 = varsData1[varID][levelID];
        cdo_read_field(streamID1, field1);

        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
      }

      int numFields2;
      int numSets = 0;
      while ((numFields2 = cdo_stream_inq_timestep(streamID2, tsID2)))
      {
        auto year = taxisInqVdatetime(taxisID2).date.year;

        if (allSteps == false && year1 != year) break;

        for (int fieldID = 0; fieldID < numFields2; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);
          auto const &var = varList2.vars[varID];
          field.init(var);
          cdo_read_field(streamID2, field);

          auto const &field1 = varsData1[varID][levelID];
          auto &field2 = varsData2[varID][levelID];
          for (size_t i = 0; i < var.gridsize; ++i)
            if (numSets == (int) std::lround(field1.vec_d[i]))
            {
              if (useExactDate) cdo_taxis_copy_timestep(taxisID3, taxisID2);
              field2.vec_d[i] = field.vec_d[i];
            }
        }

        numSets++;
        tsID2++;
      }

      if (numSets)
      {
        if (!useExactDate) cdo_taxis_copy_timestep(taxisID3, taxisID1);
        cdo_def_timestep(streamID3, tsID3);

        for (int fieldID = 0; fieldID < maxFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();
          if (tsID && varList1.vars[varID].isConstant) continue;

          auto &field2 = varsData2[varID][levelID];
          field_num_mv(field2);
          cdo_def_field(streamID3, varID, levelID);
          cdo_write_field(streamID3, field2);
        }

        tsID3++;
      }

      if (numSets == 0)
      {
        cdo_warning("First input stream has more timesteps than the second input stream!");
        break;
      }

      if (numFields2 == 0) break;

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
