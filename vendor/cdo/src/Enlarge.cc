/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Enlarge    enlarge         Enlarge fields
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "griddes.h"

class Enlarge : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Enlarge",
    .operators = { { "enlarge", EnlargeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Enlarge> registration = RegisterEntry<Enlarge>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  size_t ysize2{};
  size_t xsize2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  size_t gridsize2{};

  VarList varList1{};

public:
  void
  init() override
  {
    operator_check_argc(1);

    auto gridID2 = cdo_define_grid(cdo_operator_argv(0));
    xsize2 = gridInqXsize(gridID2);
    ysize2 = gridInqYsize(gridID2);

    if (Options::cdoVerbose) fprintf(stderr, "gridID2 %d, xsize2 %zu, ysize2 %zu\n", gridID2, xsize2, ysize2);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);

    gridsize2 = gridInqSize(gridID2);
    if (gridsize2 < varList1.gridsizeMax()) cdo_abort("Gridsize of input stream is greater than new gridsize!");

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    bool linfo = true;
    Varray<double> array1(gridsize2);
    Varray<double> array2(gridsize2);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto const &var = varList1.vars[varID];
        auto gridID1 = var.gridID;
        auto gridsize1 = var.gridsize;

        auto xsize1 = gridInqXsize(gridID1);
        auto ysize1 = gridInqYsize(gridID1);
        if (xsize1 == 0) xsize1 = 1;
        if (ysize1 == 0) ysize1 = 1;

        if (xsize1 == 1 && ysize1 == ysize2 && xsize1 * ysize1 == gridsize1)
        {
          if (linfo)
          {
            cdo_print("Enlarge zonal");
            linfo = false;
          }

          for (size_t iy = 0; iy < ysize2; iy++)
            for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[iy];

          if (numMissVals) numMissVals *= xsize2;
        }
        else if (ysize1 == 1 && xsize1 == xsize2 && xsize1 * ysize1 == gridsize1)
        {
          if (linfo)
          {
            cdo_print("Enlarge meridional");
            linfo = false;
          }

          for (size_t iy = 0; iy < ysize2; iy++)
            for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[ix];

          if (numMissVals) numMissVals *= ysize2;
        }
        else
        {
          varray_copy(gridsize1, array1, array2);
          for (size_t i = gridsize1; i < gridsize2; ++i) array2[i] = array1[gridsize1 - 1];

          if (numMissVals && fp_is_equal(array1[gridsize1 - 1], var.missval)) numMissVals += (gridsize2 - gridsize1);
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, array2.data(), numMissVals);
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
