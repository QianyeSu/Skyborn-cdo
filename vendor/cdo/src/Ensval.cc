/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
   This module contains the following operators:
   Ensval       enscrps          Ensemble cumulative ranked probability score & decomposition
   Ensval       ensbrs           Ensemble Brier score & decomposition

   The implementation of the decomposition and score calculation as carried out in this routine follows the paper
     Hans Hersbach (2000): Decomposition of the Continuous Ranked Probability
     Score for Ensemble Prediction Systems, in: Weather and Forecasting (15) pp. 559-570
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_options.h"
#include "param_conversion.h"
#include "util_files.h"
#include "util_string.h"

enum OPERTYPE
{
  CRPS,
  BRS
};

enum RESTYPE_BRS
{
  BRS_RES,
  BRS_RELI,
  BRS_RESOL,
  BRS_UNCTY
};

enum RESTYPE_CRPS
{
  CRPS_RES,
  CRPS_RELI,
  CRPS_POT
};

class Ensval : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ensval",
    .operators = { { "enscrps", CRPS, 0, EnsvalHelp }, { "ensbrs", BRS, 0, EnsvalHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, OBASE, NoRestriction },
  };
  inline static RegisterEntry<Ensval> registration = RegisterEntry<Ensval>(module);

  int nostreams = 0, numGrids;
  size_t numMissVals{};
  size_t gridsize = 0;
  int vlistID{};
  int have_miss = 0;
  int stream = 0;
  CdoStreamID streamID = 0;
  double missval = 0;
  double sum_weights = 0;
  double crps_reli = 0, crps_pot = 0, crps = 0;
  double heavyside0{}, heavysideN{};
  double brs_reli{}, brs_resol{}, brs_uncty{}, brs_thresh = 0;

  struct EnsFile
  {
    CdoStreamID streamID;
    int vlistID;
    VarList varList;
    Varray<double> array;
  };

  int numFiles{};

  int operfunc{};
  int taxisID1{ CDI_UNDEFID };
  int gridID2{};

  int nens{};

  Varray<double> alpha, beta, alpha_weights, beta_weights;
  Varray<double> brs_g, brs_o;
  Varray<double> weights;
  // int vlistCheck, gridsizeCheck;

  VarList varList1;

  std::vector<EnsFile> ensFileList;
  std::vector<int> vlistID2;
  std::vector<int> taxisID2;
  std::vector<CdoStreamID> streamID2;

  Varray<double> results;
  Varray<double> val;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    numFiles = cdo_stream_cnt();
    nens = numFiles - 1;

    if (operfunc == CRPS) { nostreams = 3; }
    else if (operfunc == BRS)
    {
      operator_input_arg("Threshold for Brier score?");
      operator_check_argc(1);
      brs_thresh = parameter_to_double(cdo_operator_argv(0));
      nostreams = 4;

      fprintf(stderr, "brs_thres %10.6f\n", brs_thresh);
    }

    // allocate array to hold results
    results = Varray<double>(nostreams);

    // one stream for each value of the decomposition
    streamID2 = std::vector<CdoStreamID>(nostreams);

    vlistID2 = std::vector<int>(nostreams);
    taxisID2 = std::vector<int>(nostreams);
    val = Varray<double>(nens, 0);

    if (operfunc == CRPS)
    {
      alpha.resize(nens + 1, 0);
      beta.resize(nens + 1, 0);
      alpha_weights.resize(nens + 1, 0);
      beta_weights.resize(nens + 1, 0);
    }
    else if (operfunc == BRS)
    {
      brs_g.resize(nens + 1, 0);
      brs_o.resize(nens + 1, 0);
    }

    if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", nens);

    ensFileList = std::vector<EnsFile>(numFiles);

    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      streamID = cdo_open_read(fileIdx);
      vlistID = cdo_stream_inq_vlist(streamID);

      auto &ensFile = ensFileList[fileIdx];
      ensFile.streamID = streamID;
      ensFile.vlistID = vlistID;
      ensFile.varList = VarList(ensFile.vlistID);
    }

    auto streamID1 = ensFileList[0].streamID;

    if (Options::cdoVerbose) cdo_print("Opened %i Input Files for Ensemble Operator", numFiles);

    // check for identical contents of all ensemble members
    auto nvars = vlistNvars(ensFileList[0].vlistID);
    if (Options::cdoVerbose) cdo_print("nvars %i", nvars);

    for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx) { varList_compare(ensFileList[0].varList, ensFileList[fileIdx].varList); }

    auto vlistID1 = ensFileList[0].vlistID;
    taxisID1 = vlistInqTaxis(vlistID1);

    gridID2 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
    double xval = 0, yval = 0;
    gridDefXvals(gridID2, &xval);
    gridDefYvals(gridID2, &yval);

    auto ofilebase = cdo_get_obase();
    auto fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    for (stream = 0; stream < nostreams; stream++)
    {
      const char *typeSuffix = nullptr;
      switch (operfunc)
      {
        case CRPS:
          switch (stream)
          {
            case 0: typeSuffix = "crps"; break;
            case 1: typeSuffix = "crps_reli"; break;
            case 2: typeSuffix = "crps_pot"; break;
          }
          break;
        case BRS:
          switch (stream)
          {
            case 0: typeSuffix = "brs"; break;
            case 1: typeSuffix = "brs_reli"; break;
            case 2: typeSuffix = "brs_reso"; break;
            case 3: typeSuffix = "brs_unct"; break;
          }
          break;
      }

      auto fileName = ofilebase + string_format(".%s", typeSuffix);
      if (fileName.size() > 0) fileName += fileSuffix;

      if (!Options::cdoOverwriteMode && FileUtils::file_exists(fileName) && !FileUtils::user_file_overwrite(fileName))
        cdo_abort("Outputfile %s already exists!", fileName);

      streamID2[stream] = open_write(fileName);

      taxisID2[stream] = taxisDuplicate(taxisID1);
      vlistID2[stream] = vlistDuplicate(vlistID1);

      numGrids = vlistNumGrids(vlistID2[stream]);
      for (int i = 0; i < numGrids; ++i) vlistChangeGridIndex(vlistID2[stream], i, gridID2);

      vlistDefTaxis(vlistID2[stream], taxisID2[stream]);
      cdo_def_vlist(streamID2[stream], vlistID2[stream]);
    }

    if (Options::cdoVerbose) cdo_print(" sum_weights %10.6f", sum_weights);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    int numFields = 0;
    int tsID = 0;
    do {
      auto numFields0 = cdo_stream_inq_timestep(ensFileList[0].streamID, tsID);
      for (int fileIdx = 1; fileIdx < numFiles; ++fileIdx)
      {
        streamID = ensFileList[fileIdx].streamID;
        numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields != numFields0)
          cdo_abort("Number of fields at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(0),
                    cdo_get_stream_name(fileIdx));
      }

      for (stream = 0; stream < nostreams; stream++)
      {
        cdo_taxis_copy_timestep(taxisID2[stream], taxisID1);
        if (numFields0 > 0) cdo_def_timestep(streamID2[stream], tsID);
      }

      while (numFields0-- > 0)
      {
        int varID = 0, levelID = 0;
        for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
        {
          auto &ensFile = ensFileList[fileIdx];
          std::tie(varID, levelID) = cdo_inq_field(ensFile.streamID);

          if (fileIdx == 0)
          {
            gridsize = varList1.vars[varID].gridsize;
            missval = varList1.vars[varID].missval;
            weights.resize(gridsize);
          }

          ensFile.array.resize(gridsize);
          cdo_read_field(ensFile.streamID, ensFile.array.data(), &numMissVals);
        }

        // xsize = gridInqXsize(gridID);
        // ysize = gridInqYsize(gridID);

        /*
        if (xsize > 1 && ysize > 1)
        {
          gridcell_weights(gridID, weights);
          sum_weights = varray_sum(gridsize, weights);
        }
        else
        */
        {
          std::ranges::fill(weights, 1.0 / gridsize);
          sum_weights = 1.0;
        }

        numMissVals = 0;
        heavyside0 = 0;
        heavysideN = 0;

        for (size_t i = 0; i < gridsize; ++i)
        {
          double xa = 0.0;
          have_miss = 0;
          for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
          {
            auto value = ensFileList[fileIdx].array[i];
            if (fileIdx == 0)
              xa = value;  // 1st file contains reference
            else
              val[fileIdx - 1] = value;  // Ensembles start at 2nd file

            if (fp_is_equal(value, missval))
            {
              have_miss = 1;
              break;
            }
          }

          auto &x = val;
          std::ranges::sort(x);  // Sort The Ensemble Array to ascending order

          // only process if no missing value in ensemble
          if (!have_miss && operfunc == CRPS)
          {
            if (xa < x[0])
            { /* Consider outliers            */
              beta[0] += (x[0] - xa) * weights[i];
              heavyside0 += 1.;
            }
            if (xa > x[nens - 1])
            {
              alpha[nens] += (xa - x[nens - 1]) * weights[i];
              alpha_weights[nens] += weights[i];
              heavysideN += 1.;
            }

            // Loop start at zero ==> 1st ensemble (c-indexing)
            for (int k = 0; k < nens - 1; ++k)
            {                     // Cumulate alpha and beta
              if (xa > x[k + 1])  // left of heavyside
                alpha[k + 1] += (x[k + 1] - x[k]) * weights[i];
              else if (xa < x[k])  // right of heavyside
                beta[k + 1] += (x[k + 1] - x[k]) * weights[i];
              else if (x[k + 1] >= xa && xa >= x[k])  // hitting jump pf heavyside (occurs exactly once!)
                beta[k + 1] += (x[k + 1] - xa) * weights[i];
            }
          }
          else if (operfunc == BRS)
          {
            // int occ = xa > brs_thresh? 1 : 0;

            // brs_g[i] - number of enemble members with rank i that
            // forecast event
            //          - event: value > brs_thresh
            //
            if (x[0] > brs_thresh)
              brs_g[0] += weights[i];
            else if (x[nens - 1] < brs_thresh)
              brs_g[nens] += weights[i];
            else
              for (int k = 0; k < nens - 1; ++k)
              {
                if (x[k + 1] >= brs_thresh && brs_thresh >= x[k])
                {
                  brs_g[k + 1] += weights[i];
                  break;
                }
              }

            // brs_o[i] - number of times that the obs is between Ensemble i-1 and i
            if (1)
            {
              if (x[0] > xa)
                brs_o[0] += weights[i];
              else if (x[nens - 1] < xa)
                brs_o[nens] += weights[i];
              else
                for (int k = 0; k < nens - 1; ++k)
                {
                  if (x[k + 1] >= xa && xa >= x[k])
                  {
                    brs_o[k + 1] += weights[i];
                    break;
                  }
                }
            }
          }
        }  // for ( i=0; i<gridsize; i++ )

        if (operfunc == CRPS)
        {
          // First Bin
          double p = 0.0, g = 0.0;
          double o = heavyside0 / gridsize;
          if (o > 0.) g = beta[0] / o;
          crps_reli = g * (o - p) * (o - p);
          crps_pot = g * o * (1. - o);
          crps = g * ((1. - o) * p * p + o * (1. - p) * (1. - p));

          // Middle Bins
          for (int k = 1; k < nens; ++k)
          {
            p = (double) k / (double) nens;

            if (fp_is_not_equal(sum_weights, 1.0))
            {
              alpha[k] /= sum_weights;
              beta[k] /= sum_weights;
            }

            g = alpha[k] + beta[k];
            o = beta[k] / (alpha[k] + beta[k]);

            crps_reli += g * (o - p) * (o - p);
            crps_pot += g * o * (1. - o);
            crps += g * ((1. - o) * p * p + o * (1. - p) * (1. - p));
          }

          // Last Bin
          p = 1.0;
          g = 0.0;
          o = 1.0 - heavysideN / gridsize;
          if (is_not_equal(o, 1.))
          {
            g = alpha[nens] / (1 - o);

            crps_reli += g * (o - p) * (o - p);
            crps_pot += g * o * (1 - o);
            crps += g * ((1 - o) * p * p + o * (1 - p) * (1 - p));
          }
          results[CRPS_RES] = crps;
          results[CRPS_RELI] = crps_reli;
          results[CRPS_POT] = crps_pot;
        }
        else if (operfunc == BRS)
        {
          brs_reli = 0;
          brs_resol = 0;
          brs_uncty = 0;

          double gsum = 0, obar = 0, osum = 0;
          for (int k = 0; k <= nens; ++k)
          {
            obar += brs_g[k] * brs_o[k];
            gsum += brs_g[k];
            osum += brs_o[k];
          }

          if (std::fabs(osum - 1) > 1.e-06 || std::fabs(gsum - 1) > 1.e-06)
          {
            cdo_abort("Internal error - normalization constraint of problem not fulfilled");
            cdo_abort("This is likely due to missing values");
          }

          brs_uncty = obar * (1 - obar);

          for (int k = 0; k <= nens; ++k)
          {
            auto g = brs_g[k];
            auto o = brs_o[k];
            auto p = 1.0 - k / (float) nens;
            // need p = 1 - k/nens here as k=0 if all members forecast
            // event and k=nens if none does so.

            brs_reli += g * (o - p) * (o - p);
            brs_resol += g * (o - obar) * (o - obar);
          }

          results[BRS_RES] = brs_reli - brs_resol + brs_uncty;
          results[BRS_RELI] = brs_reli;
          results[BRS_RESOL] = brs_resol;
          results[BRS_UNCTY] = brs_uncty;

          if (Options::cdoVerbose)
          {
            cdo_print("BRS: obar %12.6g brs  %12.6g reli %12.6g resol %12.6g u %12.6g", obar, brs_reli - brs_resol + brs_uncty,
                      brs_reli, brs_resol, brs_uncty);
          }
        }

        if (Options::cdoVerbose && operfunc == CRPS)
          cdo_print("CRPS:%12.6g reli:%12.6g crps_pot:%12.6g crps:%12.6g", crps, crps_reli, crps_pot, crps_reli + crps_pot);

        for (stream = 0; stream < nostreams; stream++)
        {
          cdo_def_field(streamID2[stream], varID, levelID);
          if (std::isnan(results[stream]))
          {
            results[stream] = missval;
            have_miss = 1;
          }
          cdo_write_field(streamID2[stream], &results[stream], have_miss);
        }

        switch (operfunc)
        {
          case (CRPS):
            std::ranges::fill(alpha, 0.0);
            std::ranges::fill(beta, 0.0);
            heavyside0 = 0;
            heavysideN = 0;
            break;
          case (BRS):
            std::ranges::fill(brs_o, 0.0);
            std::ranges::fill(brs_g, 0.0);
            break;
        }
      }  // while (numFields0-- > 0)
      tsID++;
    } while (numFields);
  }

  void
  close() override
  {
    for (auto &ensFile : ensFileList) { cdo_stream_close(ensFile.streamID); }

    for (stream = 0; stream < nostreams; stream++) { cdo_stream_close(streamID2[stream]); }

    for (stream = 0; stream < nostreams; stream++)
    {
      vlistDestroy(vlistID2[stream]);
      taxisDestroy(taxisID2[stream]);
    }
  }
};
