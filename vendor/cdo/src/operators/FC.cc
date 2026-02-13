/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      FC         fc2sp           Fourier to spectral
      FC         sp2fc           Spectral to fourier
      FC         fc2gp           Fourier to gridpoint
      FC         gp2fc           Gridpoint to fourier
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#endif
#include "cdo_fftw3.h"

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "specspace.h"
#include "griddes.h"
#include "cdo_omp.h"

static int
vlistGetFirstReg2DGrid(int vlistID)
{
  // find first gaussian grid
  auto numGrids = vlistNumGrids(vlistID);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(vlistID, index);
    if (gridInqType(gridID) == GRID_GAUSSIAN || gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_CURVILINEAR)
      return gridID;
  }

  return -1;
}

class FC : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "FC",
    // clang-format off
    .operators= { { "fc2sp" },
                  { "sp2fc" },
                  { "fc2gp" },
                  { "gp2fc" },
                  { "fourier2grid", 1, 0, nullptr },
                  { "grid2fourier", 1, 0, nullptr } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<FC> registration = RegisterEntry<FC>();

  int FC2SP{}, SP2FC{}, FC2GP{}, GP2FC{}, GRID2FOURIER{}, FOURIER2GRID{};
  int gridID2 = -1;

  int operatorID{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool dataIsUnchanged{};

  Varray<double> array1{};
  Varray<double> array2{};

  std::vector<bool> vars{};

  FC_Transformation fcTrans{};
  SP_Transformation spTrans{};

public:
  void
  init_fc2sp(int &gridID1, int &gridIDsp, int gridIDfc)
  {
    gridID1 = gridIDfc;
    if (gridID1 == -1) { cdo_warning("No fourier data found!"); }
    else
    {
      long nfc = gridInqSize(gridID1);
      long ntr = gridInqTrunc(gridID1);
      long nlat = nfc_to_nlat(nfc, ntr);

      if (gridIDsp != -1)
        if (ntr != gridInqTrunc(gridIDsp)) gridIDsp = -1;

      if (gridIDsp == -1)
      {
        long nsp = (ntr + 1) * (ntr + 2);
        gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
        gridDefTrunc(gridIDsp, ntr);
        gridDefComplexPacking(gridIDsp, 1);
      }

      gridID2 = gridIDsp;
      long nlon = 2 * nlat;
      ntr = gridInqTrunc(gridID2);

      fcTrans.init(nlon, nlat, ntr);
    }
  }

  void
  init_sp2fc(int &gridID1, int gridIDsp, int &gridIDfc)
  {
    gridID1 = gridIDsp;
    if (gridID2 == -1) { cdo_warning("No spectral data found!"); }
    else
    {
      long ntr = gridInqTrunc(gridID1);
      long nlat = ntr_to_nlat(ntr);

      if (gridIDfc != -1)
      {
        if (ntr != gridInqTrunc(gridIDfc)) gridIDfc = -1;
      }

      if (gridIDfc == -1)
      {
        long nfc = 2 * nlat * (ntr + 1);
        gridIDfc = gridCreate(GRID_FOURIER, nfc);
        gridDefTrunc(gridIDfc, ntr);
      }

      gridID2 = gridIDfc;
      long nlon = 2 * nlat;

      fcTrans.init(nlon, nlat, ntr);
    }
  }

  void
  init_gp2fc(int &gridID1, int gridIDgp, int &gridIDfc)
  {
    gridID1 = gridIDgp;
    if (gridID1 == -1) { cdo_warning("No Gaussian grid data found!"); }
    else
    {
      long nlon = gridInqXsize(gridID1);
      long nlat = gridInqYsize(gridID1);
      long ntr = nlat_to_ntr(nlat);

      if (gridIDfc != -1)
        if (ntr != gridInqTrunc(gridIDfc)) gridIDfc = -1;

      if (gridIDfc == -1)
      {
        long nfc = 2 * nlat * (ntr + 1);
        gridIDfc = gridCreate(GRID_FOURIER, nfc);
        gridDefTrunc(gridIDfc, ntr);
      }

      gridID2 = gridIDfc;

      spTrans.init(nlon, nlat, ntr, PolFlag::SP2FC);
    }
  }

  void
  init_fc2gp(int &gridID1, int &gridIDgp, int gridIDfc)
  {
    gridID1 = gridIDfc;
    if (gridID1 == -1) { cdo_warning("No fourier data found!"); }
    else
    {
      long nfc = gridInqSize(gridID1);
      long ntr = gridInqTrunc(gridID1);
      long nlat = nfc_to_nlat(nfc, ntr);

      if (gridIDgp != -1)
      {
        if (nlat != (long) gridInqYsize(gridIDgp)) gridIDgp = -1;
      }

      if (gridIDgp == -1)
      {
        char gridname[20];
        std::snprintf(gridname, sizeof(gridname), "t%dgrid", (int) ntr);

        gridIDgp = grid_from_name(gridname);
      }

      gridID2 = gridIDgp;
      long nlon = gridInqXsize(gridID2);
      nlat = gridInqYsize(gridID2);

      spTrans.init(nlon, nlat, ntr, PolFlag::FC2SP);
    }
  }

  void
  init_fourier2grid(int &gridID1)
  {
    gridID1 = vlistGetFirstReg2DGrid(vlistID1);
    if (gridID1 == -1) { cdo_warning("No regular 2D data found!"); }
    else { gridID2 = gridID1; }
  }

  void
  init_grid2fourier(int &gridID1)
  {
    gridID1 = vlistGetFirstReg2DGrid(vlistID1);
    if (gridID1 == -1) { cdo_warning("No regular 2D data found!"); }
    else { gridID2 = gridID1; }
  }

  void
  init() override
  {
    operator_check_argc(0);

    dataIsUnchanged = data_is_unchanged();

    FC2SP = module.get_id("fc2sp");
    SP2FC = module.get_id("sp2fc");
    FC2GP = module.get_id("fc2gp");
    GP2FC = module.get_id("gp2fc");
    GRID2FOURIER = module.get_id("grid2fourier");
    FOURIER2GRID = module.get_id("fourier2grid");

    operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    int gridIDsp = (operfunc == 0) ? vlist_get_first_spectral_grid(vlistID1) : -1;
    int gridIDgp = (operfunc == 0) ? vlist_get_first_gaussian_grid(vlistID1) : -1;
    int gridIDfc = (operfunc == 0) ? vlist_get_first_fourier_grid(vlistID1) : -1;

    int gridID1 = -1;
    // define output grid
    if (operatorID == FC2SP) { init_fc2sp(gridID1, gridIDsp, gridIDfc); }
    else if (operatorID == SP2FC) { init_sp2fc(gridID1, gridIDsp, gridIDfc); }
    else if (operatorID == GP2FC) { init_gp2fc(gridID1, gridIDgp, gridIDfc); }
    else if (operatorID == FC2GP) { init_fc2gp(gridID1, gridIDgp, gridIDfc); }
    else if (operatorID == FOURIER2GRID) { init_fourier2grid(gridID1); }
    else if (operatorID == GRID2FOURIER) { init_grid2fourier(gridID1); }

    auto numVars = varList1.numVars();
    vars.resize(numVars);
    for (auto const &var : varList1.vars) vars[var.ID] = (gridID1 == var.gridID);

    if (gridID1 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);
    if (operatorID == GRID2FOURIER)
    {
      for (int varID = 0; varID < numVars; ++varID)
        if (vars[varID]) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_CPX64);
    }
    else if (operatorID == FOURIER2GRID)
    {
      for (int varID = 0; varID < numVars; ++varID)
        if (vars[varID]) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
    }

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    auto gridsizeMax = varList1.gridsizeMax();
    if (operatorID == FOURIER2GRID) gridsizeMax *= 2;
    array1.resize(gridsizeMax);

    if (gridID2 != -1)
    {
      auto gridsize = gridInqSize(gridID2);
      if (operatorID == GRID2FOURIER) gridsize *= 2;
      array2.resize(gridsize);
    }
  }

  void
  run() override
  {
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

        if (vars[varID])
        {
          size_t numMissVals;
          cdo_read_field(streamID1, array1.data(), &numMissVals);
          if (numMissVals) cdo_abort("Missing values unsupported for spectral/fourier data!");

          auto gridID1 = varList1.vars[varID].gridID;
          if (operatorID == FC2SP) { four2spec(spTrans, gridID1, array1, gridID2, array2); }
          else if (operatorID == SP2FC) { spec2four(spTrans, gridID1, array1, gridID2, array2); }
          else if (operatorID == FC2GP) { four2grid(fcTrans, gridID1, array1, gridID2, array2); }
          else if (operatorID == GP2FC) { grid2four(fcTrans, gridID1, array1, gridID2, array2); }
          else if (operatorID == FOURIER2GRID) { fourier2grid(gridID1, array1, array2); }
          else if (operatorID == GRID2FOURIER) { grid2fourier(gridID1, array1, gridID2, array2); }

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, array2.data(), numMissVals);
        }
        else
        {
          cdo_def_field(streamID2, varID, levelID);
          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            size_t numMissVals;
            cdo_read_field(streamID1, array1.data(), &numMissVals);
            cdo_write_field(streamID2, array1.data(), numMissVals);
          }
        }
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
