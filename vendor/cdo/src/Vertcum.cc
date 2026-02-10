/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertcum    vertcum         Vertical cumulative
      Vertcum    vertcumhl       Vertical cumulative on hybrid sigma half level
*/

#include <cdi.h>

#include "process_int.h"

#define IS_SURFACE_LEVEL(zaxisID) (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

static void
add_vars_mv(size_t gridsize, double missval, Varray<double> const &var1, Varray<double> const &var2, Varray<double> &var3)
{
  auto missval1 = missval;
  auto missval2 = missval;
  for (size_t i = 0; i < gridsize; ++i)
  {
    var3[i] = var2[i];
    if (!fp_is_equal(var1[i], missval1))
    {
      if (!fp_is_equal(var2[i], missval2))
        var3[i] += var1[i];
      else
        var3[i] = var1[i];
    }
  }
}

class Vertcum : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertcum",
    .operators = { { "vertcum" }, { "vertcumhl" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertcum> registration = RegisterEntry<Vertcum>(module);

  int VERTCUMHL;
  int nlevshl = 0;

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operatorID;

  int numVars;

  VarList varList1;
  VarList varList2;

  std::vector<std::vector<size_t>> varnumMissVals;
  Varray3D<double> vardata1, vardata2;

public:
  void
  init() override
  {
    VERTCUMHL = module.get_id("vertcumhl");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    int zaxisIDhl = -1;

    if (operatorID == VERTCUMHL)
    {
      std::vector<double> vct;
      bool lhybrid = false;
      auto numZaxes = vlistNumZaxis(vlistID1);
      for (int i = 0; i < numZaxes; ++i)
      {
        auto zaxisID = vlistZaxis(vlistID1, i);
        auto nlevs = zaxisInqSize(zaxisID);

        if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevs > 1)
        {
          int nvct = zaxisInqVctSize(zaxisID);
          if (nlevs == (nvct / 2 - 1))
          {
            if (lhybrid == false)
            {
              lhybrid = true;
              nlevshl = nlevs + 1;

              vct.resize(nvct);
              zaxisInqVct(zaxisID, vct.data());

              zaxisIDhl = zaxisCreate(ZAXIS_HYBRID_HALF, nlevshl);
              std::vector<double> levels(nlevshl);
              for (int levelID = 0; levelID < nlevshl; ++levelID) levels[levelID] = levelID + 1;
              zaxisDefLevels(zaxisIDhl, levels.data());
              zaxisDefVct(zaxisIDhl, nvct, vct.data());
              vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
            }
            else if (vct.size())
            {
              if (std::memcmp(vct.data(), zaxisInqVctPtr(zaxisID), nvct * sizeof(double)) == 0)
                vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
            }
          }
        }
      }
    }

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    numVars = varList1.numVars();
    varnumMissVals = std::vector<std::vector<size_t>>(numVars);
    vardata1 = Varray3D<double>(numVars);
    vardata2 = Varray3D<double>(numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = varList1.vars[varID].gridsize;
      auto nlevs = varList1.vars[varID].nlevels;
      auto nlevs2 = varList2.vars[varID].nlevels;

      varnumMissVals[varID].resize(nlevs);
      vardata1[varID].resize(nlevs);
      vardata2[varID].resize(nlevs2);
      for (int levelID = 0; levelID < nlevs; ++levelID) vardata1[varID][levelID].resize(gridsize);
      for (int levelID = 0; levelID < nlevs2; ++levelID) vardata2[varID][levelID].resize(gridsize);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

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

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, &vardata1[varID][levelID][0], &varnumMissVals[varID][levelID]);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        auto missval = var2.missval;
        auto gridsize = var2.gridsize;
        auto nlevs2 = var2.nlevels;

        if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
        {
          for (size_t i = 0; i < gridsize; ++i) vardata2[varID][0][i] = 0;
        }
        else
        {
          for (size_t i = 0; i < gridsize; ++i) vardata2[varID][0][i] = vardata1[varID][0][i];
        }

        for (int levelID = 1; levelID < nlevs2; ++levelID)
        {
          if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
            add_vars_mv(gridsize, missval, vardata1[varID][levelID - 1], vardata2[varID][levelID - 1], vardata2[varID][levelID]);
          else
            add_vars_mv(gridsize, missval, vardata1[varID][levelID], vardata2[varID][levelID - 1], vardata2[varID][levelID]);
        }

        if (operatorID == VERTCUMHL && nlevs2 == nlevshl)
        {
          auto const &var1data = vardata2[varID][nlevs2 - 1];
          for (int levelID = 0; levelID < nlevs2; ++levelID)
          {
            auto &var2data = vardata2[varID][levelID];
            for (size_t i = 0; i < gridsize; ++i)
            {
              if (is_not_equal(var1data[i], 0))
                var2data[i] /= var1data[i];
              else
                var2data[i] = 0;
            }
          }
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        auto missval = var2.missval;
        auto gridsize = var2.gridsize;
        auto nlevs2 = var2.nlevels;
        for (int levelID = 0; levelID < nlevs2; ++levelID)
        {
          auto &single = vardata2[varID][levelID];
          auto numMissVals = varray_num_mv(gridsize, single, missval);
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, single.data(), numMissVals);
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

    vlistDestroy(vlistID2);
  }
};
