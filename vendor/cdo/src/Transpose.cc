/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Transpose  transxy         Transpose X/Y
*/

#include <cdi.h>

#include "cpp_lib.h"
/*
#ifdef HAVE_LIB_MDSPAN
#include <mdspan>
#else
*/
#include "matrix_view.h"
// #endif
#include "process_int.h"

void
transxy(int gridID, Varray<double> const &v1, Varray<double> &v2)
{
  auto nx = gridInqXsize(gridID);
  auto ny = gridInqYsize(gridID);
  auto gridsize = gridInqSize(gridID);

  if (gridsize == (nx * ny))
  {
    /*
#ifdef HAVE_LIB_MDSPAN
  printf("using mdspan\n");
  auto mV1 = std::mdspan(v1.data(), ny, nx);
  auto mV2 = std::mdspan(v2.data(), nx, ny);
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) mV2[i, j] = mV1[j, i];
#else
*/
    MatrixView<const double> mV1(v1.data(), ny, nx);
    MatrixView<double> mV2(v2.data(), nx, ny);
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) mV2[i][j] = mV1[j][i];
    // #endif
  }
  else
  {
    for (size_t i = 0; i < gridsize; ++i) v2[i] = v1[i];
  }
}

class Transpose : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Transpose",
    .operators = { { "transxy" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Transpose> registration = RegisterEntry<Transpose>(module);

private:
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);
      auto nx = gridInqXsize(gridID1);
      auto ny = gridInqYsize(gridID1);
      auto gridsize = gridInqSize(gridID1);
      if (gridsize == (nx * ny))
      {
        auto gridID2 = gridCreate(GRID_GENERIC, gridsize);
        gridDefXsize(gridID2, ny);
        gridDefYsize(gridID2, nx);
        vlistChangeGridIndex(vlistID2, index, gridID2);
      }
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Varray<double> array1(varList1.gridsizeMax());
    Varray<double> array2(varList1.gridsizeMax());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        size_t numMissVals;
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto gridID = varList1.vars[varID].gridID;
        transxy(gridID, array1, array2);

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
