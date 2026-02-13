/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>

void genGridIndex(int gridID1, int gridID2, std::vector<long> &index);

class Mergegrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Mergegrid",
    .operators = { { "mergegrid", MergegridHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Mergegrid> registration = RegisterEntry<Mergegrid>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};

  VarList varList1{};
  VarList varList2{};

  size_t gridsize1{};
  size_t gridsize2{};

  std::vector<long> gindex;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);

    streamID2 = cdo_open_read(1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2, CmpVarList::Name | CmpVarList::NumLevels);

    int numDiffGrids = 0;
    for (int index = 1; index < vlistNumGrids(vlistID1); ++index)
      if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) numDiffGrids++;

    if (numDiffGrids > 0) cdo_abort("Too many different grids in %s!", cdo_get_stream_name(0));

    numDiffGrids = 0;
    for (int index = 1; index < vlistNumGrids(vlistID2); ++index)
      if (vlistGrid(vlistID2, 0) != vlistGrid(vlistID2, index)) numDiffGrids++;

    if (numDiffGrids > 0) cdo_abort("Too many different grids in %s!", cdo_get_stream_name(1));

    auto gridID1 = vlistGrid(vlistID1, 0);
    auto gridID2 = vlistGrid(vlistID2, 0);

    gridsize1 = gridInqSize(gridID1);
    gridsize2 = gridInqSize(gridID2);

    gindex.resize(gridsize2);
    genGridIndex(gridID1, gridID2, gindex);

    auto vlistID3 = vlistDuplicate(vlistID1);
    streamID3 = cdo_open_write(2);

    vlistDefTaxis(vlistID3, taxisID3);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Varray<double> array1(gridsize1);
    Varray<double> array2(gridsize2);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");

      if (numFields != numFields2) cdo_abort("Input streams have different number of fields!");

      cdo_def_timestep(streamID3, tsID);

      while (numFields-- > 0)
      {
        (void) cdo_inq_field(streamID2);
        size_t numMissVals2;
        cdo_read_field(streamID2, array2.data(), &numMissVals2);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals1;
        cdo_read_field(streamID1, array1.data(), &numMissVals1);

        auto missval1 = varList1.vars[varID].missval;
        auto missval2 = varList2.vars[varID].missval;

        for (size_t i = 0; i < gridsize2; ++i)
        {
          if (gindex[i] >= 0 && fp_is_not_equal(array2[i], missval2)) { array1[gindex[i]] = array2[i]; }
        }

        if (numMissVals1)
        {
          numMissVals1 = 0;
          for (size_t i = 0; i < gridsize1; ++i)
            if (fp_is_equal(array1[i], missval1)) numMissVals1++;
        }

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, array1.data(), numMissVals1);
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
