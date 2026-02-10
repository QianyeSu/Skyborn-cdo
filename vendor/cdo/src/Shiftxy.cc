/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "matrix_view.h"

static void
shiftx(bool fillCyclic, int numberOfShifts, int nx, int ny, Varray<double> &v1, Varray<double> &v2, double missval)
{
  MatrixView<double> mv1(v1.data(), ny, nx);
  MatrixView<double> mv2(v2.data(), ny, nx);

  for (int i = 0; i < nx; ++i)
  {
    auto isCyclic = false;
    auto ins = i + numberOfShifts % nx;
    while (ins >= nx)
    {
      ins -= nx;
      isCyclic = true;
    }
    while (ins < 0)
    {
      ins += nx;
      isCyclic = true;
    }

    if (!fillCyclic && isCyclic)
    {
      for (int j = 0; j < ny; ++j) mv2[j][ins] = missval;
    }
    else
    {
      for (int j = 0; j < ny; ++j) mv2[j][ins] = mv1[j][i];
    }
  }
}

static void
shifty(bool fillCyclic, int numberOfShifts, int nx, int ny, Varray<double> &v1, Varray<double> &v2, double missval)
{
  MatrixView<double> mv1(v1.data(), ny, nx);
  MatrixView<double> mv2(v2.data(), ny, nx);

  for (int j = 0; j < ny; ++j)
  {
    auto isCyclic = false;
    auto jns = j + numberOfShifts % ny;

    while (jns >= ny)
    {
      jns -= ny;
      isCyclic = true;
    }
    while (jns < 0)
    {
      jns += ny;
      isCyclic = true;
    }

    if (!fillCyclic && isCyclic)
    {
      for (int i = 0; i < nx; ++i) mv2[jns][i] = missval;
    }
    else
    {
      for (int i = 0; i < nx; ++i) mv2[jns][i] = mv1[j][i];
    }
  }
}

static int
shiftx_coord(bool fillCyclic, int numberOfShifts, int gridID1)
{
  auto gridID2 = gridDuplicate(gridID1);

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);
  if (gridInqType(gridID1) != GRID_CURVILINEAR) ny = 1;

  Varray<double> v1(nx * ny), v2(nx * ny);
  gridInqXvals(gridID1, v1.data());
  shiftx(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
  gridDefXvals(gridID2, v2.data());

  if (gridInqXbounds(gridID1, nullptr))
  {
    size_t nv = (gridInqType(gridID1) != GRID_CURVILINEAR) ? 2 : 4;

    Varray<double> bounds(nx * ny * nv);
    gridInqXbounds(gridID1, bounds.data());
    for (size_t k = 0; k < nv; ++k)
    {
      for (size_t i = 0; i < nx * ny; ++i) v1[i] = bounds[i * nv + k];
      shiftx(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
      for (size_t i = 0; i < nx * ny; ++i) bounds[i * nv + k] = v2[i];
    }
    gridDefXbounds(gridID2, bounds.data());
  }

  return gridID2;
}

static int
shifty_coord(bool fillCyclic, int numberOfShifts, int gridID1)
{
  auto gridID2 = gridDuplicate(gridID1);

  auto nx = gridInqXsize(gridID1);
  auto ny = gridInqYsize(gridID1);
  if (gridInqType(gridID1) != GRID_CURVILINEAR) nx = 1;

  Varray<double> v1(nx * ny), v2(nx * ny);
  gridInqYvals(gridID1, v1.data());
  shifty(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
  gridDefYvals(gridID2, v2.data());

  if (gridInqYbounds(gridID1, nullptr))
  {
    size_t nv = (gridInqType(gridID1) != GRID_CURVILINEAR) ? 2 : 4;

    Varray<double> bounds(nx * ny * nv);
    gridInqYbounds(gridID1, bounds.data());
    for (size_t k = 0; k < nv; ++k)
    {
      for (size_t i = 0; i < nx * ny; ++i) v1[i] = bounds[i * nv + k];
      shifty(fillCyclic, numberOfShifts, nx, ny, v1, v2, 0);
      for (size_t i = 0; i < nx * ny; ++i) bounds[i * nv + k] = v2[i];
    }
    gridDefYbounds(gridID2, bounds.data());
  }

  return gridID2;
}

class Shiftxy : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Shiftxy",
    .operators = { { "shiftx", ShiftxyHelp }, { "shifty", ShiftxyHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Shiftxy> registration = RegisterEntry<Shiftxy>(module);

private:
  int SHIFTX{}, SHIFTY{};
  int operatorID{};

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;

  std::vector<bool> vars;

  bool fillCyclic = false;
  bool shiftCoords = false;

  int numberOfShifts = 1;

public:
  void
  init() override
  {
    SHIFTX = module.get_id("shiftx");
    SHIFTY = module.get_id("shifty");

    operatorID = cdo_operator_id();

    if (cdo_operator_argc() > 0)
    {
      numberOfShifts = parameter_to_int(cdo_operator_argv(0));
      auto numArgs = cdo_operator_argc();
      auto const &argList = cdo_get_oper_argv();
      for (int ic = 1; ic < numArgs; ++ic)
      {
        if (argList[ic] == "cyclic") { fillCyclic = true; }
        else if (argList[ic] == "coord") { shiftCoords = true; }
      }
    }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numVars = varList1.numVars();
    vars = std::vector<bool>(numVars, false);

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID1);

      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR
          || (gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL)
          || (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0))
      {
        if (shiftCoords)
        {
          int gridID2 = -1;
          if (operatorID == SHIFTX) { gridID2 = shiftx_coord(fillCyclic, numberOfShifts, gridID1); }
          else if (operatorID == SHIFTY) { gridID2 = shifty_coord(fillCyclic, numberOfShifts, gridID1); }

          vlistChangeGridIndex(vlistID2, index, gridID2);
        }

        for (auto const &var : varList1.vars)
          if (gridID1 == var.gridID) vars[var.ID] = true;
      }
      else if (gridtype == GRID_GENERIC && gridInqXsize(gridID1) <= 1 && gridInqYsize(gridID1) <= 1) {}
      else { cdo_abort("Unsupported grid type: %s", gridNamePtr(gridtype)); }
    }

    {
      int varID = 0;
      for (; varID < numVars; ++varID)
        if (vars[varID]) break;
      if (varID >= numVars) cdo_warning("No variables selected!");
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array1(gridsizeMax);
    Varray<double> array2(gridsizeMax);

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

        if (vars[varID])
        {
          auto const &var1 = varList1.vars[varID];
          auto nx = gridInqXsize(var1.gridID);
          auto ny = gridInqYsize(var1.gridID);

          if (operatorID == SHIFTX) { shiftx(fillCyclic, numberOfShifts, nx, ny, array1, array2, var1.missval); }
          else if (operatorID == SHIFTY) { shifty(fillCyclic, numberOfShifts, nx, ny, array1, array2, var1.missval); }

          numMissVals = varray_num_mv(var1.gridsize, array2, var1.missval);
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
