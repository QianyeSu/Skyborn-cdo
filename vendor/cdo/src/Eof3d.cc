/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

     EOF3d        eof3d             3D-EOF in spatial or time space
     EOF3d        eof3dspatial      3D-EOF in spatial space
     EOF3d        eof3dtime         3D-EOF in time space
*/
/*
 * TODO:
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "eigen_solution.h"
#include "eof_mode.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

class EOF3d : public Process
{
  enum
  {
    EOF3D_,
    EOF3D_TIME,
    EOF3D_SPATIAL
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "EOF3d",
    .operators = { { "eof3d", EOF3D_, 0, EOFsHelp },
                   { "eof3dspatial", EOF3D_SPATIAL, 0, EOFsHelp },
                   { "eof3dtime", EOF3D_TIME, 0, EOFsHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 2, OnlyFirst },
  };
  inline static RegisterEntry<EOF3d> registration = RegisterEntry<EOF3d>(module);

  size_t temp_size = 0, npack = 0;
  bool missval_warning = false;

  int calendar = CALENDAR_STANDARD;

  double sumWeights{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };

  VarList varList1{};

  T_WEIGHT_MODE weight_mode{};
  T_EIGEN_MODE eigen_mode{};

  size_t gridsizeMax{};
  int numSteps{};
  int numEigenFunctions{};
  int numFields{};
  int numVars{};

  Varray<double> in;
  Varray2D<int> datacounts;
  Varray3D<double> datafields;
  Varray3D<double> eigenvectors;
  Varray3D<double> eigenvalues;
  Varray<double> weights;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    operator_input_arg("Number of eigen functions to write out");
    numEigenFunctions = parameter_to_int(cdo_operator_argv(0));

    eigen_mode = get_eigenmode();
    weight_mode = get_weightmode();

    // eigenvalues

    if (operfunc == EOF3D_SPATIAL) cdo_abort("Operator not Implemented - use eof3d or eof3dtime instead");

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    // COUNT NUMBER OF TIMESTEPS if EOF3D_ or EOF3D_TIME
    numSteps = varList1.numSteps();
    if (numSteps == -1)
    {
      numSteps = 0;
      while (cdo_stream_inq_timestep(streamID1, numSteps)) numSteps++;

      if (Options::cdoVerbose) cdo_print("Counted %d timeSteps", numSteps);

      cdo_stream_close(streamID1);

      streamID1 = cdo_open_read(0);
      vlistID1 = cdo_stream_inq_vlist(streamID1);
    }
    else if (Options::cdoVerbose) { cdo_print("Found %d timeSteps", numSteps); }

    taxisID1 = vlistInqTaxis(vlistID1);

    // reset the requested number of eigen-function to the maximum if neccessary
    if (numEigenFunctions > numSteps)
    {
      cdo_warning("Solving in time-space:");
      cdo_warning("Number of eigen-functions to write out is bigger than number of time-steps.");
      cdo_warning("Setting n_eig to %d.", numSteps);
      numEigenFunctions = numSteps;
    }

    numVars = varList1.numVars();

    auto gridID1 = varList1.vars[0].gridID;
    gridsizeMax = varList1.gridsizeMax();

    // allocation of temporary fields and output structures

    in = Varray<double>(gridsizeMax);
    datacounts = Varray2D<int>(numVars);
    datafields = Varray3D<double>(numVars);
    eigenvectors = Varray3D<double>(numVars);
    eigenvalues = Varray3D<double>(numVars);

    int maxLevels = 0;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      auto gridsize = varList1.gridsizeMax();
      temp_size = gridsize * var1.nlevels;

      maxLevels = std::max(maxLevels, var1.nlevels);

      datafields[varID].resize(numSteps);
      for (int tsID = 0; tsID < numSteps; ++tsID) datafields[varID][tsID].resize(temp_size, 0.0);

      datacounts[varID].resize(temp_size, 0);

      eigenvectors[varID].resize(numEigenFunctions);
      eigenvalues[varID].resize(numSteps);

      for (int i = 0; i < numSteps; ++i)
      {
        if (i < numEigenFunctions) eigenvectors[varID][i].resize(temp_size, var1.missval);
        eigenvalues[varID][i].resize(1, var1.missval);
      }
    }

    if (Options::cdoVerbose)
      cdo_print("Allocate eigenvalue/eigenvector object with numSteps=%d and gridSize=%zu for processing in %s (%zu Bytes)",
                numSteps, gridsizeMax, "time_space", numSteps * gridsizeMax * 8);

    weights = Varray<double>(maxLevels * gridsizeMax, 1.0);

    if (weight_mode == WEIGHT_ON)
    {
      auto wstatus = gridcell_weights(gridID1, weights);
      if (wstatus != 0)
      {
        weight_mode = WEIGHT_OFF;
        cdo_warning("Using constant grid cell area weights!");
      }
      else
      {
        for (int k = 1; k < maxLevels; ++k)
          for (size_t i = 0; i < gridsizeMax; ++i) weights[k * gridsizeMax + i] = weights[i];
      }
    }
  }

  void
  run() override
  {
    int tsID = 0;

    // read the data and create covariance matrices for each var & level
    while (true)
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        size_t numMissVals;
        cdo_read_field(streamID1, in.data(), &numMissVals);

        auto offset = var1.gridsize * levelID;
        for (size_t i = 0; i < var1.gridsize; ++i)
        {
          if (fp_is_not_equal(in[i], var1.missval))
          {
            datafields[varID][tsID][offset + i] = in[i];
            datacounts[varID][offset + i]++;
          }
          else
          {
            if (datacounts[varID][offset + i] != 0) cdo_abort("Missing values unsupported!");
            if (!missval_warning)
            {
              // cdo_warning("Missing Value Support not checked for this Operator!");
              // cdo_warning("Does not work with changing locations of missing values in time.");
              missval_warning = true;
            }
            datafields[varID][tsID][i + offset] = 0;
          }
        }
      }
      tsID++;
    }

    if (Options::cdoVerbose) cdo_print("Read data for %d variables", numVars);

    Varray<size_t> pack(temp_size);  // TODO

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      temp_size = var1.gridsize * var1.nlevels;

      if (Options::cdoVerbose)
      {
        cdo_print("============================================================================");
        cdo_print("Calculating covariance matrix and SVD for var%d (%s)", varID + 1, var1.name);
      }

      npack = 0;  // TODO already set to 0

      for (size_t i = 0; i < temp_size; ++i)
      {
        if (datacounts[varID][i] > 1)
        {
          pack[npack] = i;
          npack++;
        }
      }

      sumWeights = 1;
      if (weight_mode == WEIGHT_ON)
      {
        sumWeights = 0;
        for (size_t i = 0; i < npack; ++i) sumWeights += weights[pack[i]];
      }

      if (npack < 1)
      {
        cdo_warning("Refusing to calculate EOF from a single time step for var%d (%s)", varID + 1, var1.name);
        continue;
      }

      if (Options::cdoVerbose)
        cdo_print("Allocate covar-matrix with %dx%d elements (npack=%zu) for %s (%zu Bytes)", numSteps, numSteps, npack, var1.name,
                  numSteps * numSteps * sizeof(double));

      Varray<double> eigv(numSteps);
      Varray2D<double> covar(numSteps);
      for (int j1 = 0; j1 < numSteps; ++j1) covar[j1].resize(numSteps);

      if (Options::cdoVerbose)
      {
        cdo_print("varID %d allocated eigv and cov with %dx%d elements", varID + 1, numSteps, numSteps);
        cdo_print("   npack=%zu, nts=%d temp_size=%zu", npack, numSteps, temp_size);
      }

      auto const &data = datafields[varID];
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
      for (int j1 = 0; j1 < numSteps; ++j1)
      {
        auto const &df1p = data[j1];
        for (int j2 = j1; j2 < numSteps; ++j2)
        {
          auto const &df2p = data[j2];
          double sum = 0.0;
          for (size_t i = 0; i < npack; ++i) sum += weights[pack[i] % gridsizeMax] * df1p[pack[i]] * df2p[pack[i]];
          covar[j2][j1] = covar[j1][j2] = sum / sumWeights / numSteps;
        }
      }

      if (Options::cdoVerbose) cdo_print("calculated cov-matrix");

      // SOLVE THE EIGEN PROBLEM
      if (Options::cdoVerbose) cdo_print("Processed correlation matrix for var %d | npack: %zu", varID + 1, numSteps);

      if (eigen_mode == JACOBI)
        parallel_eigen_solution_of_symmetric_matrix(covar, eigv, numSteps, __func__);
      else
        eigen_solution_of_symmetric_matrix(covar, eigv, numSteps, __func__);
      // NOW: covar contains the eigenvectors, eigv the eigenvalues

      if (Options::cdoVerbose) cdo_print("Processed SVD decomposition for var %d from %dx%d matrix", varID + 1, numSteps, numSteps);

      for (int eofID = 0; eofID < numSteps; eofID++) eigenvalues[varID][eofID][0] = eigv[eofID];

      for (int eofID = 0; eofID < numEigenFunctions; eofID++)
      {
        double *eigenvec = eigenvectors[varID][eofID].data();

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
        for (size_t i = 0; i < npack; ++i)
        {
          double sum = 0.0;
          for (int j = 0; j < numSteps; ++j) { sum += datafields[varID][j][pack[i]] * covar[eofID][j]; }
          eigenvec[pack[i]] = sum;
        }

        // NORMALIZING
        double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static) reduction(+ : sum)
#endif
        for (size_t i = 0; i < npack; ++i) sum += weights[pack[i] % gridsizeMax] * eigenvec[pack[i]] * eigenvec[pack[i]];

        if (sum > 0)
        {
          sum = std::sqrt(sum);
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < npack; ++i) eigenvec[pack[i]] /= sum;
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < npack; ++i) eigenvec[pack[i]] = var1.missval;
        }
      }  // for ( eofID = 0; eofID < n_eig; eofID++ )
    }  // for ( varID = 0; varID < numVars; varID++ )

    // write files with eigenvalues (ID3) and eigenvectors (ID2)

    // eigenvalues
    streamID2 = cdo_open_write(1);

    auto vDateTime0 = cdiDateTime_set(10101, 0);

    auto vlistID2 = vlistDuplicate(vlistID1);
    auto taxisID2 = taxisDuplicate(taxisID1);
    cdiDefKeyInt(taxisID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
    taxisDefRdatetime(taxisID2, vDateTime0);
    vlistDefTaxis(vlistID2, taxisID2);

    auto gridID2 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
    double xvals = 0.0, yvals = 0.0;
    gridDefXvals(gridID2, &xvals);
    gridDefYvals(gridID2, &yvals);

    auto numGrids = vlistNumGrids(vlistID2);
    for (int i = 0; i < numGrids; ++i) vlistChangeGridIndex(vlistID2, i, gridID2);

    auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, 1);
    double zValue = 0;
    zaxisDefLevels(zaxisID2, &zValue);
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "zaxis_Reduced");
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_LONGNAME, "Reduced zaxis from EOF3D - only one eigen value per 3D eigen vector");

    auto numZaxes = vlistNumZaxis(vlistID2);
    for (int i = 0; i < numZaxes; ++i) vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

    // eigenvectors
    streamID3 = cdo_open_write(2);

    auto vlistID3 = vlistDuplicate(vlistID1);
    auto taxisID3 = taxisDuplicate(taxisID1);
    cdiDefKeyInt(taxisID3, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
    taxisDefRdatetime(taxisID3, vDateTime0);
    vlistDefTaxis(vlistID3, taxisID3);

    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_vlist(streamID3, vlistID3);

    auto julianDate = julianDate_encode(calendar, vDateTime0);
    for (tsID = 0; tsID < numSteps; ++tsID)
    {
      julianDate = julianDate_add_seconds(julianDate, 60);
      auto vDateTime = julianDate_decode(calendar, julianDate);

      taxisDefVdatetime(taxisID2, vDateTime);
      cdo_def_timestep(streamID2, tsID);

      if (tsID < numEigenFunctions)
      {
        taxisDefVdatetime(taxisID3, vDateTime);
        cdo_def_timestep(streamID3, tsID);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var1 = varList1.vars[varID];
        for (int levelID = 0; levelID < var1.nlevels; ++levelID)
        {
          auto offset = levelID * gridsizeMax;
          if (tsID < numEigenFunctions)
          {
            auto numMissVals = array_num_mv(gridsizeMax, &eigenvectors[varID][tsID][offset], var1.missval);
            cdo_def_field(streamID3, varID, levelID);
            cdo_write_field(streamID3, &eigenvectors[varID][tsID][offset], numMissVals);
          }
        }

        auto numMissVals = (fp_is_equal(eigenvalues[varID][tsID][0], var1.missval)) ? 1 : 0;

        cdo_def_field(streamID2, varID, 0);
        cdo_write_field(streamID2, eigenvalues[varID][tsID].data(), numMissVals);
      }  // for ( varID = 0; ... )
    }  // for ( tsID = 0; ... )
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
