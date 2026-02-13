/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Merstat    merrange        Meridional range
      Merstat    mermin          Meridional minimum
      Merstat    mermax          Meridional maximum
      Merstat    mersum          Meridional sum
      Merstat    mermean         Meridional mean
      Merstat    meravg          Meridional average
      Merstat    merstd          Meridional standard deviation
      Merstat    merstd1         Meridional standard deviation [Normalize by (n-1)]
      Merstat    mervar          Meridional variance
      Merstat    mervar1         Meridional variance [Normalize by (n-1)]
      Merstat    merpctl         Meridional percentiles
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "field_functions.h"

class Merstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Merstat",
    .operators = { { "merrange", FieldFunc_Range, 0, MerstatHelp },
                   { "mermin", FieldFunc_Min, 0, MerstatHelp },
                   { "mermax", FieldFunc_Max, 0, MerstatHelp },
                   { "mersum", FieldFunc_Sum, 0, MerstatHelp },
                   { "mermean", FieldFunc_Meanw, 1, MerstatHelp },
                   { "meravg", FieldFunc_Avgw, 1, MerstatHelp },
                   { "merstd", FieldFunc_Stdw, 1, MerstatHelp },
                   { "merstd1", FieldFunc_Std1w, 1, MerstatHelp },
                   { "mervar", FieldFunc_Varw, 1, MerstatHelp },
                   { "mervar1", FieldFunc_Var1w, 1, MerstatHelp },
                   { "merskew", FieldFunc_Skew, 0, MerstatHelp },
                   { "merkurt", FieldFunc_Kurt, 0, MerstatHelp },
                   { "mermedian", FieldFunc_Median, 0, MerstatHelp },
                   { "merpctl", FieldFunc_Pctl, 0, MerstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Merstat> registration = RegisterEntry<Merstat>();

  int gridID2 = -1, lastgrid = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool needWeights{};
  int operfunc{};

  double pn = 0.0;

  VarList varList1{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);
    needWeights = (cdo_operator_f2(operatorID) != 0);

    if (operfunc == FieldFunc_Pctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
    }
    else { operator_check_argc(0); }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids = vlistNumGrids(vlistID1);
    int numDiffGrids = 0;
    for (int index = 1; index < numGrids; ++index)
      if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) numDiffGrids++;

    if (numDiffGrids > 0) cdo_abort("Too many different grids!");

    auto gridID1 = varList1.vars[0].gridID;
    auto gridType = varList1.vars[0].gridType;
    if (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN || gridType == GRID_GENERIC) { gridID2 = gridToMeridional(gridID1); }
    else { cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridType)); }

    int index = 0;
    vlistChangeGridIndex(vlistID2, index, gridID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridID1 = varList1.vars[0].gridID;
    int nlonmax = gridInqXsize(gridID1);  // max nlon?

    Field field1;
    if (needWeights) field1.weightv.resize(varList1.gridsizeMax());

    Field field2;
    field2.resize(nlonmax);
    field2.grid = gridID2;
    field2.memType = MemType::Double;

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
        auto const &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID1, field1);

        field2.missval = field1.missval;

        auto wstatus = false;
        if (needWeights && field1.grid != lastgrid)
        {
          lastgrid = field1.grid;
          wstatus = gridcell_weights(field1.grid, field1.weightv);
        }

        if (wstatus != 0 && tsID == 0 && levelID == 0)
          cdo_warning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", var.name);

        (operfunc == FieldFunc_Pctl) ? meridional_pctl(field1, field2, pn) : meridional_function(field1, field2, operfunc);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field2);
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
