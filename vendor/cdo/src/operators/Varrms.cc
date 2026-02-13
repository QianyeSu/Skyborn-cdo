/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "arithmetic.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "cdi_lockedIO.h"
#include "field_functions.h"

static void
var_rms(Varray<double> const &w, FieldVector const &field1, FieldVector const &field2, Field &field3)
{
  auto is_EQ = fp_is_equal;
  auto grid1 = field1[0].grid;
  auto grid2 = field2[0].grid;
  auto missval1 = field1[0].missval;
  auto missval2 = field2[0].missval;
  double rsum = 0.0, rsumw = 0.0;

  auto nlev = field1.size();
  auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2)) cdo_abort("fields have different size!");

  // if ( numMissVals1 )
  {
    for (size_t k = 0; k < nlev; ++k)
    {
      auto array1 = field1[k].vec_d;
      auto array2 = field2[k].vec_d;
      for (size_t i = 0; i < len; ++i)
      //	  if ( !is_EQ(w[i], missval1) )
      {
        rsum = ADDM(rsum, MULM(w[i], MULM(SUBM(array2[i], array1[i]), SUBM(array2[i], array1[i]))));
        rsumw = ADDM(rsumw, w[i]);
      }
    }
  }
  /*
else
  {
    for ( i = 0; i < len; i++ )
      {
        rsum  += w[i] * array1[i];
        rsumw += w[i];
      }
  }
  */

  auto ravg = SQRTM(DIVM(rsum, rsumw));

  size_t rnumMissVals = 0;
  if (is_EQ(ravg, missval1)) rnumMissVals++;

  field3.vec_d[0] = ravg;
  field3.numMissVals = rnumMissVals;
}

class Varrms : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Varrms",
    .operators = { { "varrms" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Varrms> registration = RegisterEntry<Varrms>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID3{};

  int gridID3{ -1 };

  int numVars{};
  int lastgrid = -1;
  int oldcode = 0;

  bool needWeights = true;

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    double slon = 0.0, slat = 0.0;
    gridID3 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID3, 1);
    gridDefYsize(gridID3, 1);
    gridDefXvals(gridID3, &slon);
    gridDefYvals(gridID3, &slat);

    vlistClearFlag(vlistID1);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) vlistDefFlag(vlistID1, varID, 0, true);

    vlistID3 = vlistCreate();
    cdo_vlist_copy_flag(vlistID3, vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    auto numGrids = vlistNumGrids(vlistID1);
    int index = 0;
    auto gridID1 = vlistGrid(vlistID1, index);

    if (needWeights && gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
      cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

    vlistChangeGridIndex(vlistID3, index, gridID3);
    if (numGrids > 1) cdo_abort("Too many different grids!");

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    Varray<double> weights;
    if (needWeights) weights.resize(varList1.gridsizeMax());
    Field field3;
    field3.resize(1);
    field3.grid = gridID3;

    FieldVector2D varsData1, varsData2;
    field2D_init(varsData1, varList1, FIELD_VEC);
    field2D_init(varsData2, varList2, FIELD_VEC);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          cdo_read_field(streamID1, varsData1[varID][levelID]);
          if (varsData1[varID][levelID].numMissVals) cdo_abort("Missing values unsupported for this operator!");
        }
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);
          cdo_read_field(streamID2, varsData2[varID][levelID]);
          if (varsData2[varID][levelID].numMissVals) cdo_abort("Missing values unsupported for this operator!");
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto wstatus = false;
        auto gridID = varsData1[varID][0].grid;
        if (needWeights && gridID != lastgrid)
        {
          lastgrid = gridID;
          wstatus = gridcell_weights(gridID, weights);
        }
        auto code = vlistInqVarCode(vlistID1, varID);
        if (wstatus != 0 && tsID == 0 && code != oldcode) cdo_warning("Using constant area weights for code %d!", oldcode = code);

        field3.missval = varsData1[varID][0].missval;
        var_rms(weights, varsData1[varID], varsData2[varID], field3);

        cdo_def_field(streamID3, varID, 0);
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

    vlistDestroy(vlistID3);
  }
};
