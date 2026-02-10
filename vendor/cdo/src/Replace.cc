/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Replace    replace         Replace variables
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"

class Replace : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Replace",
    .operators = { { "replace", ReplaceHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Replace> registration = RegisterEntry<Replace>(module);

private:
  static const int MaxVars = 1024;
  int nchvars = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID3{};
  int taxisID1{ CDI_UNDEFID };

  int vars1[MaxVars], vars2[MaxVars];
  std::vector<std::vector<int>> varlevel;
  std::vector<std::vector<size_t>> varnumMissVals2;

  VarList varList1{};
  VarList varList2{};
  Varray2D<double> vardata2;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);

    streamID2 = cdo_open_read(1);

    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    // compare all variables in vlistID2

    auto numVars1 = varList1.numVars();
    auto numVars2 = varList2.numVars();

    for (int varID2 = 0; varID2 < numVars2; ++varID2)
    {
      auto const &var2 = varList2.vars[varID2];
      int varID1 = 0;
      for (; varID1 < numVars1; ++varID1)
      {
        if (varList1.vars[varID1].name == var2.name) break;
      }

      if (varID1 < numVars1)
      {
        auto const &var1 = varList1.vars[varID1];

        if (var1.gridsize != var2.gridsize) cdo_abort("Variables have different gridsize!");
        if (var1.nlevels < var2.nlevels) cdo_abort("Variables have different number of levels!");

        if (Options::cdoVerbose) cdo_print("Variable %s replaced.", var1.name);

        vars1[nchvars] = varID1;
        vars2[nchvars] = varID2;
        nchvars++;
        if (nchvars > MaxVars) cdo_abort("Internal problem - too many variables!");
      }
      else { cdo_warning("Variable %s not found!", var2.name); }
    }

    if (nchvars)
    {
      vardata2.resize(nchvars);
      varnumMissVals2.resize(nchvars);
      varlevel.resize(nchvars);
      for (int idx = 0; idx < nchvars; idx++)
      {
        auto const &var1 = varList1.vars[vars1[idx]];
        auto const &var2 = varList2.vars[vars2[idx]];

        vardata2[idx].resize(var2.nlevels * var2.gridsize);
        varnumMissVals2[idx].resize(var2.nlevels);
        varlevel[idx].resize(var1.nlevels);
        // for (int levelID = 0; levelID < var1.nlevels; levelID++) varlevel[idx][levelID] = levelID;
        if (var2.nlevels <= var1.nlevels)
        {
          Varray<double> level1(var1.nlevels), level2(var2.nlevels);
          cdo_zaxis_inq_levels(var1.zaxisID, level1.data());
          cdo_zaxis_inq_levels(var2.zaxisID, level2.data());

          for (int levelID = 0; levelID < var1.nlevels; ++levelID) { varlevel[idx][levelID] = -1; }

          for (int l2 = 0; l2 < var2.nlevels; ++l2)
          {
            int l1 = 0;
            for (; l1 < var1.nlevels; ++l1)
              if (is_equal(level2[l2], level1[l1]))
              {
                varlevel[idx][l1] = l2;
                break;
              }

            if (l1 == var1.nlevels) cdo_warning("Variable %s on level %g not found!", var2.name, level2[l2]);
          }
        }
      }
    }

    auto vlistID3 = vlistDuplicate(vlistID1);

    streamID3 = cdo_open_write(2);

    vlistDefTaxis(vlistID3, taxisID3);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Varray<double> array(varList1.gridsizeMax());
    auto numSteps2 = varList2.numSteps();

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      if (tsID == 0 || (numSteps2 != 0 && numSteps2 != 1))
      {
        auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
        if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");

        for (int fieldID = 0; fieldID < numFields2; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);

          for (int idx = 0; idx < nchvars; idx++)
            if (vars2[idx] == varID)
            {
              auto offset = varList2.vars[varID].gridsize * levelID;
              size_t numMissVals = 0;
              cdo_read_field(streamID2, &vardata2[idx][offset], &numMissVals);
              varnumMissVals2[idx][levelID] = numMissVals;
              break;
            }
        }
      }

      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        auto parray = array.data();

        size_t numMissVals = 0;
        int idx = 0;
        for (; idx < nchvars; idx++)
          if (vars1[idx] == varID)
          {
            auto levelID2 = varlevel[idx][levelID];
            if (levelID2 != -1)
            {
              auto offset = varList1.vars[varID].gridsize * levelID2;
              parray = &vardata2[idx][offset];
              numMissVals = varnumMissVals2[idx][levelID2];
              break;
            }
          }

        if (idx == nchvars) cdo_read_field(streamID1, parray, &numMissVals);

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, parray, numMissVals);
      }

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
