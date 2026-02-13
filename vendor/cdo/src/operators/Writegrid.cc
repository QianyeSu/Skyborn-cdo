/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      writegrid Write grid
*/

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "griddes.h"

class Writegrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Writegrid",
    .operators = { { "writegrid" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Writegrid> registration = RegisterEntry<Writegrid>();

  CdoStreamID streamID;
  int gridID;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID = cdo_open_read(0);
    auto vlistID = cdo_stream_inq_vlist(streamID);

    gridID = vlistGrid(vlistID, 0);
    gridID = generate_full_cell_grid(gridID);

    if (!gridHasCoordinates(gridID)) cdo_abort("Cell corner coordinates missing!");
  }

  void
  run() override
  {
    auto gridsize = gridInqSize(gridID);
    std::vector<int> mask(gridsize);

    if (gridInqMask(gridID, nullptr)) { gridInqMask(gridID, mask.data()); }
    else
    {
      for (size_t i = 0; i < gridsize; ++i) mask[i] = 1;
    }

    write_nc_grid(cdo_get_stream_name(1), gridID, mask.data());
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
