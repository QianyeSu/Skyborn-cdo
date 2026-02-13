/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Eof        eof             EOF in spatial or time space
     Eof        eofspatial      EOF in spatial space
     Eof        eoftime         EOF in time space
*/
/*
 * TODO:
 * Role of the weights for eofs. Should not be mixed up with division with number of contributing values during summation.
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

// No missing value support added so far!

static void
scale_eigvec_grid(Varray<double> &out, int tsID, size_t npack, std::vector<size_t> const &pack, Varray<double> const &weight,
                  Varray2D<double> const &covar, double sumWeights)
{
  for (size_t i = 0; i < npack; ++i) out[pack[i]] = covar[tsID][i] / std::sqrt(weight[pack[i]] / sumWeights);
}

static void
scale_eigvec_time(Varray<double> &out, int tsID, int numSteps, size_t npack, std::vector<size_t> const &pack,
                  Varray<double> const &weight, Varray2D<double> const &covar, Varray2D<double> const &data, double missval,
                  double sumWeights)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, numSteps, tsID, pack, data, covar, out)
#endif
  for (size_t i = 0; i < npack; ++i)
  {
    double sum = 0.0;
    for (int j = 0; j < numSteps; ++j) sum += data[j][i] * covar[tsID][j];
    out[pack[i]] = sum;
  }
  /*
  for ( size_t j = 0; j < numSteps; ++j )
    {
      for ( size_t i = 0; i < npack; ++i )
        out[pack[i]] += data[j][i] * covar[tsID][j];
    }
  */

  // Normalizing
  double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+ : sum) shared(out, weight, pack, npack)
#endif
  for (size_t i = 0; i < npack; ++i)
  {
    // do not need to account for weights as eigenvectors are non-weighted
    sum += weight[pack[i]] * out[pack[i]] * out[pack[i]];
  }

  if (sum > 0.0)
  {
    sum = std::sqrt(sum / sumWeights);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, pack, sum, out)
#endif
    for (size_t i = 0; i < npack; ++i) out[pack[i]] /= sum;
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, pack, out, missval)
#endif
    for (size_t i = 0; i < npack; ++i) out[pack[i]] = missval;
  }
}

class Eof : public Process
{
  enum
  {
    EOF_,
    EOF_TIME,
    EOF_SPATIAL
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Eof",
    // clang-format off
    .operators = { { "eof", EOF_, 0, EofHelp },
                   { "eofspatial", EOF_SPATIAL, 0, EofHelp },
                   { "eoftime", EOF_TIME, 0, EofHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 2, OnlyFirst },
  };
  inline static RegisterEntry<Eof> registration = RegisterEntry<Eof>();

  struct eofdata_t
  {
    bool init = false;
    bool first_call = true;
    Varray<double> eigenValues;
    Varray2D<double> covar;
    Varray2D<double> data;
  };

  size_t numMissVals = 0;
  int gridSpace = 0, timeSpace = 0;
  double sum_w = 1.0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int vlistID1 = -1;
  int vlistID2 = -1;
  int vlistID3 = -1;

  int taxisID1{ CDI_UNDEFID };

  int gridID2{};

  int calendar = CALENDAR_STANDARD;

  T_WEIGHT_MODE weightMode{};
  T_EIGEN_MODE eigenMode{};

  size_t numEigenFunctions = 0;

  int numVars{};
  int numEigen{};
  int numGrids{};
  int numSteps = 1;

  size_t npack = SIZE_MAX;
  VarList varList1;
  std::vector<size_t> pack;
  Varray<double> in;
  std::vector<std::vector<eofdata_t>> eofData2D;
  Varray<double> weights;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    operator_input_arg("Number of eigen functions to write out");
    numEigen = parameter_to_int(cdo_operator_argv(0));

    eigenMode = get_eigenmode();
    weightMode = get_weightmode();

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    varList1 = VarList(vlistID1);

    auto gridID1 = varList1.vars[0].gridID;
    auto gridsizeMax = varList1.gridsizeMax();
    numVars = varList1.numVars();

    numGrids = vlistNumGrids(vlistID1);
    for (int index = 1; index < numGrids; ++index)
      if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) cdo_abort("Too many different grids!");

    // eigenvalues

    // Count number of timesteps if EOF_ or EOF_TIME
    if (operfunc == EOF_ || operfunc == EOF_TIME)
    {
      numSteps = varList1.numSteps();
      if (numSteps == 0) numSteps = 1;
      if (numSteps == -1)
      {
        numSteps = 0;
        while (cdo_stream_inq_timestep(streamID1, numSteps)) numSteps++;

        if (Options::cdoVerbose) cdo_print("Counted %d timeSteps", numSteps);

        cdo_stream_close(streamID1);

        streamID1 = cdo_open_read(0);
        vlistID1 = cdo_stream_inq_vlist(streamID1);
        taxisID1 = vlistInqTaxis(vlistID1);
      }
      else if (Options::cdoVerbose) { cdo_print("Found %d timeSteps", numSteps); }

      if ((size_t) numSteps < gridsizeMax || operfunc == EOF_TIME)
      {
        timeSpace = 1;
        gridSpace = 0;
      }
      else
      {
        timeSpace = 0;
        gridSpace = 1;
      }
    }
    else if (operfunc == EOF_SPATIAL)
    {
      timeSpace = 0;
      gridSpace = 1;
    }

    // reset the requested number of eigen-function to the maximum if neccessary
    if (timeSpace)
    {
      if (numEigen > numSteps)
      {
        cdo_warning("Solving in time-space:");
        cdo_warning("Number of eigen-functions to write out is bigger than number of time-steps.");
        cdo_warning("Setting numEigen to %d.", numSteps);
        cdo_warning("If You want to force a solution in grid-space use operator eofspatial");
        numEigen = numSteps;
      }
      numEigenFunctions = numSteps;
    }
    else if (gridSpace)
    {
      if (((double) gridsizeMax) * gridsizeMax > (double) SIZE_MAX) cdo_abort("Grid space too large!");

      if ((size_t) numEigen > gridsizeMax)
      {
        cdo_warning("Solving in spatial space");
        cdo_warning("Number of eigen-functions to write out is bigger than grid size");
        cdo_warning("Setting numEigen to %zu", gridsizeMax);
        cdo_warning("If You want to force a solution in time-space use operator eoftime");
        numEigen = gridsizeMax;
      }
      numEigenFunctions = gridsizeMax;
    }

    if (Options::cdoVerbose)
      cdo_print("Calculating %d eigenvectors and %zu eigenvalues in %s", numEigen, numEigenFunctions,
                (gridSpace == 1) ? "grid_space" : "time_space");

    weights.resize(gridsizeMax, 1.0);

    if (weightMode == WEIGHT_ON)
    {
      auto wstatus = gridcell_weights(gridID1, weights);
      if (wstatus != 0)
      {
        weightMode = WEIGHT_OFF;
        cdo_warning("Using constant grid cell area weights!");
      }
    }

    // allocation of temporary fields and output structures

    pack.resize(gridsizeMax);
    in.resize(gridsizeMax);
    eofData2D.resize(numVars);

    if (Options::cdoVerbose)
      cdo_print("Allocate eigenvalue/eigenvector object with numSteps=%d and gridSize=%zu (%zu Bytes)", numSteps, gridsizeMax,
                numSteps * gridsizeMax * 8);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto numLevels = varList1.vars[varID].nlevels;
      eofData2D[varID].resize(numLevels);

      if (timeSpace)
        for (int levelID = 0; levelID < numLevels; ++levelID) eofData2D[varID][levelID].data.resize(numSteps);
    }
  }

  void
  run() override
  {
    int tsID = 0;

    // read the data and create covariance matrices for each var & level
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, in.data(), &numMissVals);

        auto const &var1 = varList1.vars[varID];

        if (npack == SIZE_MAX)
        {
          npack = 0;
          for (size_t i = 0; i < var1.gridsize; ++i)
          {
            if (!fp_is_equal(weights[i], 0.0) && !fp_is_equal(weights[i], var1.missval) && !fp_is_equal(in[i], var1.missval))
            {
              pack[npack] = i;
              npack++;
            }
          }

          if (weightMode == WEIGHT_ON)
          {
            sum_w = 0.0;
            for (size_t i = 0; i < npack; ++i) sum_w += weights[pack[i]];
          }
        }

        {
          size_t ipack = 0;
          for (size_t i = 0; i < var1.gridsize; ++i)
          {
            if (!fp_is_equal(weights[i], 0.0) && !fp_is_equal(weights[i], var1.missval) && !fp_is_equal(in[i], var1.missval))
            {
              if (pack[ipack] != i) cdo_abort("Missing values unsupported!");
              ipack++;
            }
          }
          if (ipack != npack) cdo_abort("Missing values unsupported!");
        }

        auto &eofData = eofData2D[varID][levelID];
        if (gridSpace)
        {
          if (!eofData.init)
          {
            if (Options::cdoVerbose)
              cdo_print("Allocate covar-matrix with %zux%zu elements for %s [layer: %d] (%zu Bytes)", npack, npack, var1.name,
                        levelID + 1, npack * npack * 8);
            eofData.covar.resize(npack);
            for (size_t i = 0; i < npack; ++i) eofData.covar[i].resize(npack, 0.0);
          }

          auto &covar = eofData.covar;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(npack, covar, in, pack)
#endif
          for (size_t ipack = 0; ipack < npack; ++ipack)
          {
            auto &covar_i = covar[ipack];
            auto in_i = in[pack[ipack]];
            for (size_t jpack = ipack; jpack < npack; ++jpack) covar_i[jpack] += in_i * in[pack[jpack]];
          }
        }
        else if (timeSpace)
        {
          eofData.data[tsID].resize(npack);
          auto &data = eofData.data[tsID];

          for (size_t ipack = 0; ipack < npack; ipack++) data[ipack] = in[pack[ipack]];
        }

        eofData.init = true;
      }

      tsID++;
    }

    if (gridSpace) numSteps = tsID;

    if (tsID == 1) cdo_abort("File consists of only one timestep!");

    // write files with eigenvalues (ID3) and eigenvectors (ID2)

    // eigenvalues
    streamID2 = cdo_open_write(1);

    auto vDateTime0 = cdiDateTime_set(10101, 0);

    vlistID2 = vlistDuplicate(vlistID1);
    auto taxisID2 = taxisDuplicate(taxisID1);
    cdiDefKeyInt(taxisID2, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
    taxisDefRdatetime(taxisID2, vDateTime0);
    vlistDefTaxis(vlistID2, taxisID2);

    gridID2 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
    double xvals = 0.0, yvals = 0.0;
    gridDefXvals(gridID2, &xvals);
    gridDefYvals(gridID2, &yvals);
    for (int i = 0; i < numGrids; ++i) vlistChangeGridIndex(vlistID2, i, gridID2);

    // eigenvectors
    streamID3 = cdo_open_write(2);

    vlistID3 = vlistDuplicate(vlistID1);
    auto taxisID3 = taxisDuplicate(taxisID1);
    cdiDefKeyInt(taxisID3, CDI_GLOBAL, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
    taxisDefRdatetime(taxisID3, vDateTime0);
    vlistDefTaxis(vlistID3, taxisID3);

    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_vlist(streamID3, vlistID3);

    auto julianDate = julianDate_encode(calendar, vDateTime0);

    Varray<double> &out = in;

    int numStepsOut = (npack < (size_t) numSteps) ? npack : numSteps;

    for (tsID = 0; tsID < numStepsOut; ++tsID)
    {
      julianDate = julianDate_add_seconds(julianDate, 60);
      vDateTime0 = julianDate_decode(calendar, julianDate);

      taxisDefVdatetime(taxisID2, vDateTime0);
      cdo_def_timestep(streamID2, tsID);

      if (tsID < numEigen)
      {
        taxisDefVdatetime(taxisID3, vDateTime0);
        cdo_def_timestep(streamID3, tsID);
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var1 = varList1.vars[varID];

        for (int levelID = 0; levelID < var1.nlevels; ++levelID)
        {
          auto &eofData = eofData2D[varID][levelID];
          auto const &data = eofData.data;
          auto &covar = eofData.covar;

          if (eofData.first_call)
          {
            eofData.first_call = false;

            if (gridSpace)
            {
              if (Options::cdoVerbose) cdo_print("Processing level %d of %s", levelID + 1, var1.name);

              eofData.eigenValues.resize(npack);

              for (size_t ipack = 0; ipack < npack; ++ipack)
              {
                auto i = pack[ipack];
                for (size_t jpack = 0; jpack < npack; ++jpack)
                {
                  if (jpack < ipack) { covar[ipack][jpack] = covar[jpack][ipack]; }
                  else
                  {
                    auto j = pack[jpack];
                    covar[ipack][jpack] = covar[ipack][jpack] *                                    // covariance
                                          std::sqrt(weights[i]) * std::sqrt(weights[j]) / sum_w /  // weights
                                          numSteps;                                                // number of data contributing
                  }
                }
              }
            }
            else if (timeSpace)
            {
              if (Options::cdoVerbose)
                cdo_print("Allocate covar-matrix with %dx%d elements (npack=%zu) for %s [layer: %d] (%zu Bytes)", numSteps,
                          numSteps, npack, var1.name, levelID + 1, numSteps * numSteps * sizeof(double));

              covar.resize(numSteps);
              for (int i = 0; i < numSteps; ++i) covar[i].resize(numSteps);

              eofData.eigenValues.resize(numSteps);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numSteps, data, covar, weights, npack, pack, sum_w) schedule(static)
#endif
              for (int j1 = 0; j1 < numSteps; ++j1)
              {
                auto const &df1p = data[j1];
                for (int j2 = j1; j2 < numSteps; ++j2)
                {
                  auto const &df2p = data[j2];
                  double sum = 0.0;
                  for (size_t i = 0; i < npack; ++i) sum += weights[pack[i]] * df1p[i] * df2p[i];
                  covar[j2][j1] = covar[j1][j2] = sum / sum_w / numSteps;
                }
              }
            }

            // Solve the eigen problem
            auto &eigenValues = eofData.eigenValues;
            if (eigenMode == JACOBI)
              // TODO: use return status (>0 okay, -1 did not converge at all)
              parallel_eigen_solution_of_symmetric_matrix(covar, eigenValues, numEigenFunctions, __func__);
            else
              eigen_solution_of_symmetric_matrix(covar, eigenValues, numEigenFunctions, __func__);

            // NOW: covar contains the eigenvectors, eig_val the eigenvalues

            for (size_t i = 0; i < var1.gridsize; ++i) out[i] = var1.missval;

            // for ( int i = 0; i < n; i++ ) eig_val[i] *= sum_w;
          }  // first_call

          if (tsID < numEigen)
          {
            if (gridSpace)
              scale_eigvec_grid(out, tsID, npack, pack, weights, covar, sum_w);
            else if (timeSpace)
              scale_eigvec_time(out, tsID, numSteps, npack, pack, weights, covar, data, var1.missval, sum_w);

            numMissVals = varray_num_mv(var1.gridsize, out, var1.missval);
            cdo_def_field(streamID3, varID, levelID);
            cdo_write_field(streamID3, out.data(), numMissVals);
          }  // loop numEigen

          auto eigenValues = eofData.eigenValues.data();

          numMissVals = (fp_is_equal(eigenValues[tsID], var1.missval)) ? 1 : 0;

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, &eigenValues[tsID], numMissVals);
        }  // loop nlevs
      }    // loop nvars
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
    vlistDestroy(vlistID3);
  }
};
