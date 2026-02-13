/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
   This module contains the following operators:
   Ensstat3       ensrkhistspace   Ensemble ranked histogram averaged over time
   Ensstat3       ensrkhisttime    Ensemble ranked histogram averaged over space
   Ensstat3       ensroccurve      Ensamble Receiver Operating Characteristics
*/

#include <cdi.h>

#include <cstdlib>
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"
#include "util_files.h"
#include "cdo_omp.h"
#include "field_functions.h"

// Defines for rank histogram
enum TDATA_TYPE
{
  TIME,
  SPACE
};

#define time_data TIME
#define space_data SPACE

// Defines for Receiver Operating Characteristics (ROC)
#define DEBUG_ROC 0
enum CONTINGENCY_TYPE
{
  TP,  // TP - True positive  ( event     forecast and     occured)  HIT
  FP,  // FP - False positive ( event     forecast and not occured)  false ALARM
  FN,  // FN - False negative ( event not forecast and     ocurred)  MISSED
  TN   // TN - True negative  ( event not forecast and not ocurred)  CORRECT REJECTION
};

enum ROC_ENUM_TYPE
{
  TPR,  // TPR = True Positive Rate = TP / ( TP + FN )
  FPR   // FNR = False Negtive Rate = FN / ( FP + TN )
};

static double
roc_curve_integrate(Varray2D<double> const &roc, int n)
{
  double area = 0.0;

  for (int i = 1; i <= n; ++i)
  {
    auto x1 = roc[i][FPR];
    auto x0 = roc[i - 1][FPR];
    auto y1 = roc[i][TPR];
    auto y0 = roc[i - 1][TPR];
    auto dx = x1 - x0;
    auto dy = y1 - y0;

    auto step_area = -0.5 * dx * dy - dx * y0;
    area += step_area;
  }

  return area - 0.5;
}

class Ensstat3 : public Process
{
  enum
  {
    func_roc,
    func_rank
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ensstat3",
    .operators = { { "ensrkhistspace", func_rank, space_data, Ensstat2Help },
                   { "ensrkhisttime", func_rank, time_data, Ensstat2Help },
                   { "ensroc", func_roc, 0, Ensstat2Help } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Ensstat3> registration = RegisterEntry<Ensstat3>();

private:
  size_t numMissVals = 0;
  int chksum = 0;  // for check of histogram population
  int have_miss = 0;
  CdoStreamID streamID2 = 0;

  struct EnsFile
  {
    Varray<double> array;
    CdoStreamID streamID;
    int vlistID;
    VarList varList;
  };

  std::vector<EnsFile> ensFileList;
  int numFiles{};
  int operfunc{};
  int datafunc{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numBins = 0;
  int numEns{};

  std::vector<int> varIDs2;

  std::vector<int> hist;
  Varray2D<double> roc;             // receiver operating characteristics table
  Varray2D<int> ctg_tab;            // contingency table and histogram
  Varray<double> uThresh, lThresh;  // thresholds for histograms

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    datafunc = cdo_operator_f2(operatorID);

    if (operfunc == func_roc)
    {
      operator_input_arg("Number of eigen functions to write out");
      numBins = parameter_to_int(cdo_operator_argv(0));
    }

    numFiles = cdo_stream_cnt() - 1;
    numEns = numFiles - 1;

    if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", numFiles);

    std::string ofilename = cdo_get_stream_name(numFiles);

    if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
      cdo_abort("Outputfile %s already exists!", ofilename);

    ensFileList.resize(numFiles);

    for (int k = 0; k < numFiles; ++k)
    {
      auto streamID = cdo_open_read(k);
      auto &ensFile = ensFileList[k];
      ensFile.streamID = streamID;
      ensFile.vlistID = cdo_stream_inq_vlist(streamID);
      ensFile.varList = VarList(ensFile.vlistID);
    }

    // check for identical contents of all ensemble members
    for (int k = 1; k < numFiles; ++k) varList_compare(ensFileList[0].varList, ensFileList[k].varList);

    auto vlistID1 = ensFileList[0].vlistID;
    auto const &varList1 = ensFileList[0].varList;
    auto vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    auto numVars = varList1.numVars();
    varIDs2.resize(numVars);

    auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, numFiles);
    {
      Varray<double> levels(numFiles, 0);
      for (int i = 0; i < numFiles; ++i) levels[i] = i;
      zaxisDefLevels(zaxisID2, levels.data());
      cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "histogram_binID");
    }

    int time_mode = (datafunc == TIME) ? TIME_VARYING : TIME_CONSTANT;

    for (int varID = 0; varID < numVars; ++varID)
    {
      int gridID2;
      /* ******************************************************************** */
      /* numFiles includes the observation, so there are numFiles-1 ensembles */
      /* and exactly numFiles bins, in which the observation could fall       */
      /* ******************************************************************** */
      if (datafunc == TIME)
      {
        double val = 0.0;
        gridID2 = gridCreate(GRID_LONLAT, 1);
        gridDefXsize(gridID2, 1);
        gridDefYsize(gridID2, 1);
        gridDefXvals(gridID2, &val);
        gridDefYvals(gridID2, &val);
      }
      else  // datafunc == SPACE
        gridID2 = varList1.vars[varID].gridID;

      varIDs2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, time_mode);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (varList1.vars[varID].nlevels > 1)
      {
        cdo_warning("More than one level not supported when processing ranked histograms.");
        cdo_warning("Try to use `cdo splitlevel` to split the dataset into levels and apply");
        cdo_warning("the operator seperately to each level.");
        cdo_abort("Exit due to unsupported file structure");
      }
    }

    if (operfunc != func_roc)
    {
      streamID2 = cdo_open_write(numFiles);
      cdo_def_vlist(streamID2, vlistID2);
    }

    auto gridsizeMax = ensFileList[0].varList.gridsizeMax();

    for (int k = 0; k < numFiles; ++k) { ensFileList[k].array.resize(gridsizeMax); }

    if (operfunc == func_roc)
    {
      hist.resize(numBins);

      roc.resize(numBins + 1);
      ctg_tab.resize(numBins + 1);
      for (int i = 0; i <= numBins; ++i)
      {
        roc[i].resize(2, 0);
        ctg_tab[i].resize(4, 0);
      }

      uThresh.resize(numBins);
      lThresh.resize(numBins);
      for (int i = 0; i < numBins; ++i)
      {
        uThresh[i] = ((double) i + 1) / numBins;
        lThresh[i] = (double) i / numBins;
      }
    }
  }

  void
  run() override
  {
    /* *************************************************** */
    /* should each thread be allocating memory locally???? */
    /* ("first touch strategy")                            */
    /* --> #pragma omp parallel for ...                    */
    /* *************************************************** */
    FieldVector field;
    field.resize(Threading::ompNumMaxThreads);
    for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
    {
      field[i].resize(numFiles);
      field[i].weightv.resize(numFiles);
      for (int k = 0; k < numFiles; ++k) field[i].weightv[k] = 1;
    }

    std::vector<std::vector<int>> array2(numFiles + 1);
    if (operfunc == func_rank && datafunc == SPACE)
    {  // need to memorize data for entire grid before writing
      auto gridsizeMax = ensFileList[0].varList.gridsizeMax();
      for (int binID = 0; binID < numFiles; binID++) array2[binID].resize(gridsizeMax, 0);
    }
    else if (operfunc == func_rank)
    {  // can process data separately for each timestep and only need to cumulate values over the grid
      for (int binID = 0; binID < numFiles; binID++) array2[binID].resize(1, 0);
    }

    auto const &varList1 = ensFileList[0].varList;

    int varID = 0, levelID;
    int numFields0;
    int tsID = 0;
    do {
      numFields0 = cdo_stream_inq_timestep(ensFileList[0].streamID, tsID);
      for (int k = 1; k < numFiles; ++k)
      {
        auto streamID = ensFileList[k].streamID;
        int numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields != numFields0)
        {
          if (numFields == 0)
            cdo_abort("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(k));
          else
            cdo_abort("Inconsistent ensemble file, number of fields at time step %d of %s and %s differ!", tsID + 1,
                      cdo_get_stream_name(0), cdo_get_stream_name(k));
        }
      }

      if (operfunc == func_rank && (datafunc == TIME || tsID == 0))
      {
        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        if (numFields0 > 0) cdo_def_timestep(streamID2, tsID);
      }

      for (int fieldID = 0; fieldID < numFields0; ++fieldID)
      {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ensFileList, numFiles) private(numMissVals) lastprivate(varID, levelID)
#endif
        for (int k = 0; k < numFiles; ++k)
        {
          auto &ensFile = ensFileList[k];
          std::tie(varID, levelID) = cdo_inq_field(ensFile.streamID);
          cdo_read_field(ensFile.streamID, ensFile.array.data(), &numMissVals);
        }

        auto gridsizeMax = varList1.vars[varID].gridsize;
        auto missval = varList1.vars[varID].missval;

        numMissVals = 0;
        if (datafunc == TIME && operfunc == func_rank)
          for (int binID = 0; binID < numFiles; binID++) array2[binID][0] = 0;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for (size_t i = 0; i < gridsizeMax; ++i)
        {
          auto ompthID = cdo_omp_get_thread_num();

          field[ompthID].missval = missval;
          field[ompthID].numMissVals = 0;
          have_miss = 0;
          for (int k = 0; k < numFiles; ++k)
          {
            field[ompthID].vec_d[k] = ensFileList[k].array[i];
            if (fp_is_equal(field[ompthID].vec_d[k], missval))
            {
              have_miss = 1;
              break;
            }
          }

          // need to ignore all data for a gridpoint if a single ensemble
          // has a missing value at that gridpoint.
          if (!have_miss)  // only process if no missing value in ensemble
          {
            switch (operfunc)
            {
              case (func_rank):
              { /* ****************/
                /* RANK HISTOGRAM */
                /* ************** */
                // for ( j=0; j<numFiles; j++ )
                //   std::fprintf(stderr,"%5.2g ",field[ompthID].vec_d[j]);
                auto binID = (int) field_rank(field[ompthID]);
                // std::fprintf(stderr,"-->%i\n",binID);

                if (datafunc == SPACE && !have_miss)
                  array2[binID][i]++;
                else if (!have_miss)
                  array2[binID][0]++;

                break;
              }
              case (func_roc):
              { /* ********************************** */
                /* RECEIVER OPERATING CHARACTERISTICS */
                /* ********************************** */
                auto dat = &field[ompthID].vec_d[1];
                auto ival = (field[ompthID].vec_d[0] > 0.5) ? 1 : 0;

                for (int binID = 0; binID < numBins; binID++) hist[binID] = 0;

                for (int j = 0; j < numEns; ++j)
                  for (int binID = 0; binID < numBins; binID++)
                    if (dat[j] >= lThresh[binID] && dat[j] < uThresh[binID]) hist[binID]++;

                chksum = 0;
                for (int binID = 0; binID < numBins; binID++) chksum += hist[binID];

                if (chksum != numEns) std::exit(1);

                int cum = 0;
                if (ival == 1)
                {
                  // all true positives in first bin
                  ctg_tab[0][TP] += numEns;

                  cum += hist[0];
                  int binID = 1;
                  for (; binID < numBins; binID++)
                  {
                    ctg_tab[binID][TP] += numEns - cum;
                    ctg_tab[binID][FN] += cum;
                    cum += hist[binID];
                  }
                  ctg_tab[binID][TP] += numEns - cum;
                  ctg_tab[binID][FN] += cum;
                }
                else if (ival == 0)
                {
                  // all false positives in first bin
                  ctg_tab[0][FP] += numEns;
                  cum += hist[0];
                  int binID = 1;
                  for (; binID < numBins; binID++)
                  {
                    ctg_tab[binID][FP] += numEns - cum;
                    ctg_tab[binID][TN] += cum;
                    cum += hist[binID];
                  }
                  ctg_tab[binID][FP] += numEns - cum;
                  ctg_tab[binID][TN] += cum;
                }
                break;
              }
            }  // switch ( operfunc )
          }    // if ( ! have_miss )
        }      // for ( i=0; i<gridsize; i++ )

        if (datafunc == TIME && operfunc == func_rank)
        {
          for (int binID = 0; binID < numFiles; binID++)
          {
            double val = (double) array2[binID][0];
            //		fprintf(stderr,"%i ",(int)val);
            cdo_def_field(streamID2, varIDs2[varID], binID);
            cdo_write_field(streamID2, &val, numMissVals);
          }
          // std::fprintf(stderr,"\n");
        }
        else if (operfunc == func_roc)
        {
          if (DEBUG_ROC)
          {
            std::fprintf(stderr, "#             :     TP     FP     FN     TN         TPR        FPR\n");

            for (int binID = 0; binID <= numBins; binID++)
            {
              int p = ctg_tab[binID][TP] + ctg_tab[binID][FN];
              int n = ctg_tab[binID][FP] + ctg_tab[binID][TN];
              double tpr = ctg_tab[binID][TP] / (double) p;
              double fpr = ctg_tab[binID][FP] / (double) n;
              chksum += ctg_tab[binID][0] + ctg_tab[binID][1] + ctg_tab[binID][2] + ctg_tab[binID][3];

              roc[binID][TPR] = tpr;
              roc[binID][FPR] = fpr;

              std::fprintf(stderr, "%3i %10.4g: %6i %6i %6i %6i: %10.4g %10.4g\n", binID, (binID < numBins) ? lThresh[binID] : 1,
                           ctg_tab[binID][0], ctg_tab[binID][1], ctg_tab[binID][2], ctg_tab[binID][3], tpr, fpr);
            }
            std::fprintf(stderr, "nbins %10i\n", numBins);
            std::fprintf(stderr, "#ROC CurveArea: %10.6f\n", roc_curve_integrate(roc, numBins));
          }  // if ( DEBUG_ROC )
        }    // else if (operfunc == func_roc )
      }      // for ( fieldID=0; fieldID<numFields0; fieldID++ )
      tsID++;
    }  // do [...]
    while (numFields0 > 0);

    if (operfunc == func_rank)
    {
      int osize = (datafunc == TIME) ? 1 : array2[0].size();
      Varray<double> tmpdoub(osize);

      for (int binID = 0; binID < numFiles; binID++)
      {
        for (int i = 0; i < osize; ++i) tmpdoub[i] = (double) array2[binID][i];

        cdo_def_field(streamID2, varIDs2[varID], binID);
        cdo_write_field(streamID2, tmpdoub.data(), numMissVals);
      }
    }
    else if (operfunc == func_roc)
    {
      std::fprintf(stdout, "#             :     TP     FP     FN     TN         TPR        FPR\n");

      for (int i = 0; i <= numBins; ++i)
      {
        int p = ctg_tab[i][TP] + ctg_tab[i][FN];
        int n = ctg_tab[i][FP] + ctg_tab[i][TN];
        double tpr = ctg_tab[i][TP] / (double) p;
        double fpr = ctg_tab[i][FP] / (double) n;
        // chksum += ctg_tab[i][0] + ctg_tab[i][1] + ctg_tab[i][2] + ctg_tab[i][3];

        roc[i][TPR] = tpr;
        roc[i][FPR] = fpr;

        int sum = ctg_tab[i][TP] + ctg_tab[i][TN] + ctg_tab[i][FP] + ctg_tab[i][FN];

        std::fprintf(stdout, "%3i %10.4g: %6i %6i %6i %6i (%6i): %10.4g %10.4g\n", i, (i < numBins) ? lThresh[i] : 1, ctg_tab[i][0],
                     ctg_tab[i][1], ctg_tab[i][2], ctg_tab[i][3], sum, tpr, fpr);
      }

      std::fprintf(stdout, "#ROC CurveArea: %10.6f\n", roc_curve_integrate(roc, numBins));
    }
  }
  void
  close() override
  {
    for (auto const &ensFile : ensFileList) { cdo_stream_close(ensFile.streamID); }

    if (operfunc != func_roc) cdo_stream_close(streamID2);
  }
};
