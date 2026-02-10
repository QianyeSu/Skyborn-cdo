/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Zonstat    zonmin          Zonal minimum
      Zonstat    zonmax          Zonal maximum
      Zonstat    zonrange        Zonal range
      Zonstat    zonsum          Zonal sum
      Zonstat    zonmean         Zonal mean
      Zonstat    zonavg          Zonal average
      Zonstat    zonstd          Zonal standard deviation
      Zonstat    zonstd1         Zonal standard deviation [Normalize by (n-1)]
      Zonstat    zonvar          Zonal variance
      Zonstat    zonvar1         Zonal variance [Normalize by (n-1)]
      Zonstat    zonpctl         Zonal percentiles
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "field_functions.h"

void remap_weights_zonal_mean(int gridID1, int gridID2, Varray2D<size_t> &remapIndices, Varray2D<double> &remapWeights);
void remap_zonal_mean(const Varray2D<size_t> &remapIndices, Varray2D<double> const &remapWeights, Field const &field1,
                      Field &field2);

template <typename T>
static void
reorder_field(Varray<T> &v, std::vector<int> const &hpRingIndices, size_t numIndices)
{
  Varray<T> vtmp(numIndices);
  for (size_t i = 0; i < numIndices; ++i) vtmp[i] = v[hpRingIndices[i]];
  for (size_t i = 0; i < numIndices; ++i) v[i] = vtmp[i];
}

static void
reorder_field(Field &field, std::vector<int> const &hpRingIndices)
{
  auto numIndices = hpRingIndices.size();
  if (numIndices)
  {
    auto func = [&](auto &v) { reorder_field(v, hpRingIndices, numIndices); };
    field_operation(func, field);
  }
}

static int
define_reduced_grid(int gridID1, size_t ysize, std::vector<int> &hpReducedPoints)
{
  auto gridsize = gridInqSize(gridID1);
  auto gridID = gridCreate(GRID_GAUSSIAN_REDUCED, gridsize);
  gridDefYsize(gridID, ysize);
  gridDefReducedPoints(gridID, ysize, hpReducedPoints.data());
  return gridID;
}

class Zonstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Zonstat",
    .operators = { { "zonrange", FieldFunc_Range, 0, ZonstatHelp },
                   { "zonmin", FieldFunc_Min, 0, ZonstatHelp },
                   { "zonmax", FieldFunc_Max, 0, ZonstatHelp },
                   { "zonsum", FieldFunc_Sum, 0, ZonstatHelp },
                   { "zonmean", FieldFunc_Mean, 0, ZonstatHelp },
                   { "zonavg", FieldFunc_Avg, 0, ZonstatHelp },
                   { "zonstd", FieldFunc_Std, 0, ZonstatHelp },
                   { "zonstd1", FieldFunc_Std1, 0, ZonstatHelp },
                   { "zonvar", FieldFunc_Var, 0, ZonstatHelp },
                   { "zonvar1", FieldFunc_Var1, 0, ZonstatHelp },
                   { "zonskew", FieldFunc_Skew, 0, ZonstatHelp },
                   { "zonkurt", FieldFunc_Kurt, 0, ZonstatHelp },
                   { "zonmedian", FieldFunc_Median, 0, ZonstatHelp },
                   { "zonpctl", FieldFunc_Pctl, 0, ZonstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Zonstat> registration = RegisterEntry<Zonstat>(module);
  int gridIDdestroy = -1, gridID1 = -1, gridID2 = -1;
  int zongridID = -1;
  int sourceGridIsRegular = true;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  Field field1{}, field2{};
  std::vector<int> hpRingIndices{};

  VarList varList1{};

  int operfunc{};
  int nlatmax{};

  Varray2D<size_t> remapIndices{};
  Varray2D<double> remapWeights{};

  double pn = 0.0;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);

    if (operfunc == FieldFunc_Pctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
    }
    else if (cdo_operator_argc() == 1 && operfunc == FieldFunc_Mean)
    {
      sourceGridIsRegular = false;
      gridID2 = cdo_define_grid(cdo_operator_argv(0));
      auto gridtype = gridInqType(gridID2);
      if (gridtype != GRID_GAUSSIAN && gridtype != GRID_LONLAT) cdo_abort("Target grid type must be Gaussian or LonLat!");
      if (!gridInqYbounds(gridID2, NULL)) cdo_abort("Target grid cell bounds missing!");
      if (gridInqXsize(gridID2) > 1) cdo_abort("Target grid must be zonal!");
    }
    else { operator_check_argc(0); }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numGrids = varList1.numGrids();
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID = vlistGrid(vlistID1, index);
      auto xsize = gridInqXsize(gridID);
      auto ysize = gridInqYsize(gridID);
      auto gridtype = gridInqType(gridID);
      if (xsize > 1 || gridtype == GRID_GAUSSIAN_REDUCED || is_healpix_grid(gridID))
      {
        if (gridID1 == -1) gridID1 = gridID;
      }
      else
      {
        if (ysize > 1 && zongridID == -1) zongridID = gridID;
      }
    }

    int numDiffGrids = 0;
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID = vlistGrid(vlistID1, index);
      if (zongridID != -1 && zongridID == gridID) continue;
      if (gridID1 != gridID) numDiffGrids++;
    }
    if (zongridID == -1 && gridID1 == -1) cdo_abort("Unsupported grid type!");
    if (numDiffGrids) cdo_abort("Too many different grids!");

    if (gridID1 != -1)
    {
      auto gridtype = gridInqType(gridID1);
      if (sourceGridIsRegular)
      {
        if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED || gridtype == GRID_GENERIC
            || is_healpix_grid(gridID1))
        {
          gridID2 = (zongridID != -1 && gridInqYsize(zongridID) == gridInqYsize(gridID1)) ? zongridID : gridToZonal(gridID1);
        }
        else
        {
          if (operfunc == FieldFunc_Mean)
          {
            cdo_print("Add zonal grid description to calculate a zonal mean for data on non-rectangular grids.");
            cdo_print("A predefined zonal description is zonal_<DY>. DY is the increment of the latitudes in degrees.");
            cdo_print("Example for 2 degree latitude bins:  cdo zonmean,zonal_2 infile outfile");
          }
          cdo_abort("Unsupported gridtype: %s", gridNamePtr(gridtype));
        }
      }
      else
      {
        auto gridID = generate_full_cell_grid(gridID1);
        if (gridID != gridID1) gridIDdestroy = gridID1 = gridID;
      }
    }
    else
    {
      gridID2 = zongridID;
      cdo_warning("Input stream already contains zonal data!");
    }

    if (gridID2 == -1) cdo_abort("Internal error, target grid undefined!");

    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID2, index, gridID2);

    if (Options::cdoChunkType == CDI_UNDEFID)
    {
      for (auto const &var : varList1.vars) cdiDefKeyInt(vlistID2, var.ID, CDI_KEY_CHUNKTYPE, CDI_CHUNK_AUTO);
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    nlatmax = gridInqYsize(gridID2);

    auto zonVar = varList1.vars[0];
    zonVar.gridID = gridID2;
    zonVar.gridsize = nlatmax;
    zonVar.nlevels = 1;
    zonVar.memType = MemType::Double;
    field2.init(zonVar);

    if (!sourceGridIsRegular) remap_weights_zonal_mean(gridID1, gridID2, remapIndices, remapWeights);
    if (is_healpix_grid(gridID1))
    {
      std::vector<int> hpReducedPoints;
      hp_generate_ring_indices(cdo::get_healpix_params(gridID1), gridInqSize(gridID1), hpRingIndices, hpReducedPoints);
      gridIDdestroy = gridID1 = define_reduced_grid(gridID1, gridInqYsize(gridID2), hpReducedPoints);
    }
  }

  void
  step(int tsID, int numFields)
  {
    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsID);

    while (numFields--)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto const &var1 = varList1.vars[varID];
      field1.init(var1);
      cdo_read_field(streamID1, field1);
      field1.grid = gridID1;

      field2.missval = field1.missval;

      if (zongridID != -1 && zongridID == field1.grid) { field_ncopy(nlatmax, field1, field2); }
      else if (sourceGridIsRegular)
      {
        if (is_healpix_grid(var1.gridID)) reorder_field(field1, hpRingIndices);
        (operfunc == FieldFunc_Pctl) ? zonal_pctl(field1, field2, pn) : zonal_function(field1, field2, operfunc);
      }
      else { remap_zonal_mean(remapIndices, remapWeights, field1, field2); }

      cdo_def_field(streamID2, varID, levelID);
      cdo_write_field(streamID2, field2);
    }
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) return;
      step(tsID++, numFields);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (gridIDdestroy != -1) gridDestroy(gridIDdestroy);
  }
};
