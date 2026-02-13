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
#include "field_functions.h"

class Fldrms : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Fldrms",
    .operators = { { "fldrms" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Fldrms> registration = RegisterEntry<Fldrms>();

  int lastgrid = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};

  int gridID3 = -1;

  VarList varList1{};
  VarList varList2{};

  bool needWeights{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    needWeights = true;

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    double slon = 0.0, slat = 0.0;
    gridID3 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID3, 1);
    gridDefYsize(gridID3, 1);
    gridDefXvals(gridID3, &slon);
    gridDefYvals(gridID3, &slat);

    auto numGrids = varList1.numGrids();
    int numDiffGrids = 0;
    for (int index = 1; index < numGrids; ++index)
      if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) numDiffGrids++;

    auto gridID1 = vlistGrid(vlistID1, 0);
    auto gridID2 = vlistGrid(vlistID2, 0);

    if (gridInqSize(gridID1) != gridInqSize(gridID2)) cdo_abort("Fields have different grid size!");

    if (needWeights && gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
      cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID3, index, gridID3);

    if (numDiffGrids > 0) cdo_abort("Too many different grids!");

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Field field1, field2, field3;
    field1.resize(varList1.gridsizeMax());
    if (needWeights) field1.weightv.resize(varList1.gridsizeMax());
    field2.resize(varList1.gridsizeMax());
    field3.resize(1);
    field3.grid = gridID3;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID3, tsID);

      while (numFields--)
      {
        (void) cdo_inq_field(streamID1);
        cdo_read_field(streamID1, field1);
        auto [varID, levelID] = cdo_inq_field(streamID2);
        cdo_read_field(streamID2, field2);

        auto const &var1 = varList1.vars[varID];
        auto const &var2 = varList2.vars[varID];
        field1.grid = var1.gridID;
        field2.grid = var2.gridID;

        if (needWeights && field1.grid != lastgrid)
        {
          lastgrid = field1.grid;
          field1.weightv[0] = 1;
          if (field1.size > 1)
          {
            auto wstatus = gridcell_weights(field1.grid, field1.weightv);
            if (wstatus != 0 && tsID == 0 && levelID == 0)
              cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", var1.name);
          }
        }

        field1.missval = var1.missval;
        field2.missval = var1.missval;
        field3.missval = var1.missval;

        field_rms(field1, field2, field3);

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, field3);
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
