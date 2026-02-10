/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>

#include "process_int.h"

static void
invertLevDes(int vlistID)
{
  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID1 = vlistZaxis(vlistID, index);
    auto zaxisID2 = zaxisDuplicate(zaxisID1);

    auto nlev = zaxisInqSize(zaxisID1);
    if (nlev <= 1) continue;

    if (zaxisInqLevels(zaxisID1, nullptr))
    {
      Varray<double> yv1(nlev), yv2(nlev);
      zaxisInqLevels(zaxisID1, yv1.data());
      for (int ilev = 0; ilev < nlev; ++ilev) yv2[nlev - ilev - 1] = yv1[ilev];
      zaxisDefLevels(zaxisID2, yv2.data());
    }

    if (zaxisInqLbounds(zaxisID1, nullptr) && zaxisInqUbounds(zaxisID1, nullptr))
    {
      Varray<double> yb1(nlev), yb2(nlev);
      zaxisInqLbounds(zaxisID1, yb1.data());
      for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
      zaxisDefLbounds(zaxisID2, yb2.data());

      zaxisInqUbounds(zaxisID1, yb1.data());
      for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
      zaxisDefUbounds(zaxisID2, yb2.data());
    }

    auto zaxistype = zaxisInqType(zaxisID1);
    if (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF)
    {
      auto vctsize = zaxisInqVctSize(zaxisID1);
      if (vctsize && vctsize % 2 == 0)
      {
        Varray<double> vct1(vctsize), vct2(vctsize);
        zaxisInqVct(zaxisID1, vct1.data());
        for (int i = 0; i < vctsize / 2; ++i)
        {
          vct2[vctsize / 2 - 1 - i] = vct1[i];
          vct2[vctsize - 1 - i] = vct1[vctsize / 2 + i];
        }
        zaxisDefVct(zaxisID2, vctsize, vct2.data());
      }
    }

    vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
  }
}

class Invertlev : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Invertlev",
    .operators = { { "invertlev", InvertlevHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Invertlev> registration = RegisterEntry<Invertlev>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  Varray2D<double> vardata{};

  bool dataIsUnchanged{};
  std::vector<std::vector<size_t>> varnumMissVals{};

  size_t numMissVals{};

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    invertLevDes(vlistID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
    auto numVars = varList1.numVars();

    vardata.resize(numVars);
    varnumMissVals.resize(numVars);

    auto has3dVar = false;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      if (var.nlevels > 1)
      {
        has3dVar = true;
        vardata[varID].resize(var.gridsize * var.nlevels);
        varnumMissVals[varID].resize(var.nlevels);
      }
    }

    if (!has3dVar) cdo_warning("No variables with invertable levels found!");
  }

  void
  run() override
  {
    Varray<double> array(varList1.gridsizeMax());

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

        if (vardata[varID].size())
        {
          auto offset = varList1.vars[varID].gridsize * levelID;
          cdo_read_field(streamID1, &vardata[varID][offset], &numMissVals);
          varnumMissVals[varID][levelID] = numMissVals;
        }
        else
        {
          cdo_def_field(streamID2, varID, levelID);
          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            cdo_read_field(streamID1, array.data(), &numMissVals);
            cdo_write_field(streamID2, array.data(), numMissVals);
          }
        }
      }

      for (int varID = 0, numVars = varList1.numVars(); varID < numVars; ++varID)
      {
        if (vardata[varID].size())
        {
          auto const &var = varList1.vars[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            cdo_def_field(streamID2, varID, levelID);

            auto offset = var.gridsize * (var.nlevels - levelID - 1);
            numMissVals = varnumMissVals[varID][var.nlevels - levelID - 1];

            cdo_write_field(streamID2, &vardata[varID][offset], numMissVals);
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
