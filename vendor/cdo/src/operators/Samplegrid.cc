/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Michal Koutek, KMNI
          Uwe Schulzweida

*/

/*
   This module "SampleGrid" contains the following operators:

    samplegrid      Resample current grid with given factor, typically 2 (which will half the resolution);
                    tested on curvilinear and LCC grids;
    subgrid         Similar to selindexbox but this operator works for LCC grids (tested on HARMONIE NWP model).
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "grid_define.h"

static void
sampleData(const double *array1, int gridID1, double *array2, int gridID2, int resampleFactor)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = gridInqXsize(gridID2);
  auto nlat2 = gridInqYsize(gridID2);

  if (cdoDebugExt >= 100)
    cdo_print("%s(): (nlon1: %zu; nlat1: %zu) => (nlon2: %zu; nlat2: %zu); "
              "gridID1: %d; gridID2: %d; resampleFactor: %d)",
              __func__, nlon1, nlat1, nlon2, nlat2, gridID1, gridID2, resampleFactor);

  for (size_t ilat1 = 0; ilat1 < nlat1; ilat1 += resampleFactor)
    for (size_t ilon1 = 0; ilon1 < nlon1; ilon1 += resampleFactor) *array2++ = array1[ilat1 * nlon1 + ilon1];
}

static void
cropData(double *array1, int gridID1, double *array2, int gridID2, int subI0, int subI1, int subJ0, int subJ1)
{
  long nlon1 = gridInqXsize(gridID1);
  long nlon2 = gridInqXsize(gridID2);
  long rowLen = subI1 - subI0 + 1;  // must be same as nlon1

  if (rowLen != nlon2) cdo_abort("cropData() rowLen!= nlon2 [%d != %d]", rowLen, nlon2);

  if (cdoDebugExt >= 10) cdo_print("cropData(%d,%d,%d,%d) ...", subI0, subI1, subJ0, subJ1);

  long array2Idx = 0;
  for (long ilat1 = subJ0; ilat1 <= subJ1; ilat1++)  // copy the last row as well..
  {
    array_copy(rowLen, &array1[ilat1 * nlon1 + subI0], &array2[array2Idx]);
    array2Idx += rowLen;
  }
}

class Samplegrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Samplegrid",
    .operators = { { "samplegrid", 0, 0, "resample factor, typically 2 (which will half the resolution)", SamplegridHelp },
                   { "subgrid", 0, 0, "sub-grid indices: i0,i1,j0,j1", SamplegridHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Samplegrid> registration = RegisterEntry<Samplegrid>();

  int SAMPLEGRID{}, SUBGRID{};
  int resampleFactor = 0;
  int subI0 = 0, subI1 = 0, subJ0 = 0, subJ1 = 0;
  struct sbox_t
  {
    int gridSrcID, gridIDsampled;
  };

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1{};
  VarList varList2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  std::vector<bool> vars;

  int operatorID{};
  int numGrids{};
  std::vector<sbox_t> sbox;

public:
  void
  init() override
  {
    SAMPLEGRID = module.get_id("samplegrid");
    SUBGRID = module.get_id("subgrid");

    operatorID = cdo_operator_id();

    auto nch = cdo_operator_argc();

    if (operatorID == SAMPLEGRID)
    {
      Debug(cdoDebugExt, "samplegrid operator requested..");
      if (nch < 1) cdo_abort("Number of input arguments < 1; At least 1 argument needed: resample-factor (2,3,4, .. etc)");
      resampleFactor = parameter_to_int(cdo_operator_argv(0));

      Debug(cdoDebugExt, "resampleFactor = %d", resampleFactor);
    }
    else if (operatorID == SUBGRID)
    {
      Debug(cdoDebugExt, "subgrid operator requested..");
      if (nch < 4)
        cdo_abort("Number of input arguments < 4; Must specify sub-grid indices: i0,i1,j0,j1; This works only with LCC grid."
                  " For other grids use: selindexbox");
      subI0 = parameter_to_int(cdo_operator_argv(0));
      subI1 = parameter_to_int(cdo_operator_argv(1));
      subJ0 = parameter_to_int(cdo_operator_argv(2));
      subJ1 = parameter_to_int(cdo_operator_argv(3));
    }
    else { cdo_abort("Unknown operator ..."); }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numVars = varList1.numVars();
    vars.resize(numVars, false);

    numGrids = varList1.numGrids();
    Debug(cdoDebugExt, "numGrids=%d", numGrids);

    sbox.resize(numGrids);

    for (int index = 0; index < numGrids; ++index)
    {
      const auto gridSrcID = vlistGrid(vlistID1, index);
      int gridIDsampled = -1;

      if (gridInqSize(gridSrcID) <= 1) continue;

      int gridtype = gridInqType(gridSrcID);
      if (!(gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION || gridtype == GRID_CURVILINEAR
            || gridtype == GRID_GENERIC))
        cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridtype));

      if (operatorID == SAMPLEGRID) { gridIDsampled = cdo_define_sample_grid(gridSrcID, resampleFactor); }
      else if (operatorID == SUBGRID) { gridIDsampled = cdo_define_subgrid_grid(gridSrcID, subI0, subI1, subJ0, subJ1); }

      sbox[index].gridSrcID = gridSrcID;
      sbox[index].gridIDsampled = gridIDsampled;

      // if ( cdoDebugExt>=10 ) cdo_print_griddes(gridSrcID, 1);
      // if ( cdoDebugExt>=10 ) cdo_print_griddes(gridIDsampled, 1);

      vlistChangeGridIndex(vlistID2, index, gridIDsampled);

      for (auto const &var : varList1.vars)
        if (gridSrcID == var.gridID) vars[var.ID] = true;
    }

    Debug(cdoDebugExt,
          [&]()
          {
            if (operatorID == SAMPLEGRID) Debug("Resampled grid has been created.");
            if (operatorID == SUBGRID) Debug("Sub-grid has been created.");
          });

    varList2 = VarList(vlistID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsize1Max = varList1.gridsizeMax();
    if (vlistNumber(vlistID1) != CDI_REAL) gridsize1Max *= 2;
    auto gridsize2Max = varList2.gridsizeMax();
    if (vlistNumber(vlistID2) != CDI_REAL) gridsize2Max *= 2;
    Varray<double> array1(gridsize1Max);
    Varray<double> array2(gridsize2Max);
    Debug(cdoDebugExt, "gridsize = %ld, gridsize2 = %ld", gridsize1Max, gridsize2Max);

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
        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        cdo_def_field(streamID2, varID, levelID);

        Debug(cdoDebugExt >= 20, "Processing field (%d) of %d.", fieldID, numFields);

        if (vars[varID])
        {
          auto gridSrcID = varList1.vars[varID].gridID;

          int index;
          for (index = 0; index < numGrids; ++index)
            if (gridSrcID == sbox[index].gridSrcID) break;

          if (index == numGrids) cdo_abort("Internal problem, grid not found!");

          int gridIDsampled = sbox[index].gridIDsampled;
          gridsize2Max = gridInqSize(gridIDsampled);

          if (operatorID == SAMPLEGRID) { sampleData(array1.data(), gridSrcID, array2.data(), gridIDsampled, resampleFactor); }
          else if (operatorID == SUBGRID)
          {
            cropData(array1.data(), gridSrcID, array2.data(), gridIDsampled, subI0, subI1, subJ0, subJ1);
          }

          if (numMissVals) { numMissVals = varray_num_mv(gridsize2Max, array2, varList2.vars[varID].missval); }

          cdo_write_field(streamID2, array2.data(), numMissVals);
        }
        else { cdo_write_field(streamID2, array1.data(), numMissVals); }
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
