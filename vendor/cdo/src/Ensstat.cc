/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ensstat    ensrange        Ensemble range
      Ensstat    ensmin          Ensemble minimum
      Ensstat    ensmax          Ensemble maximum
      Ensstat    enssum          Ensemble sum
      Ensstat    ensmean         Ensemble mean
      Ensstat    ensavg          Ensemble average
      Ensstat    ensstd          Ensemble standard deviation
      Ensstat    ensstd1         Ensemble standard deviation
      Ensstat    ensvar          Ensemble variance
      Ensstat    ensvar1         Ensemble variance
      Ensstat    enspctl         Ensemble percentiles
*/

#include <atomic>

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "workerthread.h"
#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"
#include "util_files.h"
#include "cdo_omp.h"
#include "field_functions.h"

namespace
{
struct EnsFile
{
  VarList varList;
  CdoStreamID streamID;
};
}  // namespace

static void
ensstat(std::vector<EnsFile> const &ensFileList, FieldVector &fieldVector, CdoStreamID streamID2, int varID, int levelID,
        FieldVector &workFields, Varray<double> &array2, Varray<double> &count2, int operfunc, double pn)
{
  int numFiles = ensFileList.size();
  auto withCountData = (count2.size() > 0);

  auto hasMissvals = false;
  for (int k = 0; k < numFiles; ++k)
    if (fieldVector[k].numMissVals > 0) hasMissvals = true;

  auto numVars = ensFileList[0].varList.numVars();
  auto gridsize = ensFileList[0].varList.vars[varID].gridsize;
  auto missval = ensFileList[0].varList.vars[varID].missval;
  auto memType = fieldVector[0].memType;

  std::atomic<size_t> atomicNumMiss{ 0 };
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto ompthID = cdo_omp_get_thread_num();

    auto &work = workFields[ompthID];
    work.missval = missval;
    work.numMissVals = 0;
    if (memType == MemType::Float)
      for (int k = 0; k < numFiles; ++k) work.vec_d[k] = fieldVector[k].vec_f[i];
    else
      for (int k = 0; k < numFiles; ++k) work.vec_d[k] = fieldVector[k].vec_d[i];

    if (hasMissvals)
      for (int k = 0; k < numFiles; ++k)
      {
        if (fp_is_equal(work.vec_d[k], ensFileList[k].varList.vars[varID].missval))
        {
          work.vec_d[k] = missval;
          work.numMissVals++;
        }
      }

    auto lpctl = (operfunc == FieldFunc_Pctl);
    array2[i] = lpctl ? field_pctl(work, pn) : field_function(work, operfunc);

    if (fp_is_equal(array2[i], work.missval)) atomicNumMiss++;

    if (withCountData) count2[i] = numFiles - work.numMissVals;
  }

  size_t numMissVals = atomicNumMiss;

  cdo_def_field(streamID2, varID, levelID);
  cdo_write_field(streamID2, array2.data(), numMissVals);

  if (withCountData)
  {
    cdo_def_field(streamID2, varID + numVars, levelID);
    cdo_write_field(streamID2, count2.data(), 0);
  }
}

class Ensstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ensstat",
    .operators = { { "ensrange", FieldFunc_Range, 0, EnsstatHelp },
                   { "ensmin", FieldFunc_Min, 0, EnsstatHelp },
                   { "ensmax", FieldFunc_Max, 0, EnsstatHelp },
                   { "enssum", FieldFunc_Sum, 0, EnsstatHelp },
                   { "ensmean", FieldFunc_Mean, 0, EnsstatHelp },
                   { "ensavg", FieldFunc_Avg, 0, EnsstatHelp },
                   { "ensvar", FieldFunc_Var, 0, EnsstatHelp },
                   { "ensvar1", FieldFunc_Var1, 0, EnsstatHelp },
                   { "ensstd", FieldFunc_Std, 0, EnsstatHelp },
                   { "ensstd1", FieldFunc_Std1, 0, EnsstatHelp },
                   { "ensskew", FieldFunc_Skew, 0, EnsstatHelp },
                   { "enskurt", FieldFunc_Kurt, 0, EnsstatHelp },
                   { "ensmedian", FieldFunc_Median, 0, EnsstatHelp },
                   { "enspctl", FieldFunc_Pctl, 0, EnsstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Ensstat> registration = RegisterEntry<Ensstat>(module);

private:
  int operfunc{ 0 };
  double pn{ 0 };

  std::vector<EnsFile> ensFileList;

  bool printWarning = false;
  bool printError = false;

  int tsID{ 0 };
  int numFiles;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  CdoStreamID streamID2;

  bool withCountData{ false };

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    auto lpctl = (operfunc == FieldFunc_Pctl);

    auto argc = cdo_operator_argc();
    auto nargc = argc;

    if (lpctl)
    {
      operator_input_arg("percentile number");
      pn = parameter_to_double(cdo_operator_argv(0));
      argc--;
    }

    withCountData = false;
    if (argc == 1)
    {
      if (cdo_operator_argv(nargc - 1) == "count")
        withCountData = true;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(nargc - 1));
    }

    numFiles = cdo_stream_cnt() - 1;

    if (Options::cdoVerbose) cdo_print("Ensemble over %d files.", numFiles);

    cdo::set_numfiles(numFiles + 8);

    std::string ofilename = cdo_get_stream_name(numFiles);

    if (!Options::cdoOverwriteMode && FileUtils::file_exists(ofilename) && !FileUtils::user_file_overwrite(ofilename))
      cdo_abort("Outputfile %s already exists!", ofilename);

    ensFileList.resize(numFiles);

    int vlistID1 = -1;
    for (int k = 0; k < numFiles; ++k)
    {
      auto streamID = cdo_open_read(k);
      auto vlistID = cdo_stream_inq_vlist(streamID);
      ensFileList[k].streamID = streamID;
      ensFileList[k].varList = VarList(vlistID);
      if (k == 0) vlistID1 = vlistID;
    }

    // check that the contents is always the same
    for (int k = 1; k < numFiles; ++k) { varList_compare(ensFileList[0].varList, ensFileList[k].varList); }

    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if (withCountData)
    {
      auto numVars = ensFileList[0].varList.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = ensFileList[0].varList.vars[varID];
        auto name = var.name + "_count";
        auto cvarID = vlistDefVar(vlistID2, var.gridID, var.zaxisID, var.timeType);
        cdiDefKeyString(vlistID2, cvarID, CDI_KEY_NAME, name.c_str());
        vlistDefVarDatatype(vlistID2, cvarID, CDI_DATATYPE_INT16);
        if (cvarID != (varID + numVars)) cdo_abort("Internal error, varIDs do not match!");
      }
    }

    streamID2 = cdo_open_write(numFiles);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsizeMax = ensFileList[0].varList.gridsizeMax();
    Varray<double> array2(gridsizeMax);

    Varray<double> count2;
    if (withCountData) { count2.resize(gridsizeMax); }

    auto workerThread = Options::CDO_task ? std::make_unique<WorkerThread>() : nullptr;
    auto numTasks = Options::CDO_task ? 2 : 1;

    FieldVector workFields(Threading::ompNumMaxThreads);
    for (auto &work : workFields) work.resize(numFiles);

    FieldVector fieldVector[2];
    fieldVector[0].resize(numFiles);
    if (Options::CDO_task) fieldVector[1].resize(numFiles);

    int fieldNum = 0;
    int numFields0;
    do {
      numFields0 = cdo_stream_inq_timestep(ensFileList[0].streamID, tsID);

      for (int k = 1; k < numFiles; ++k)
      {
        auto streamID = ensFileList[k].streamID;
        auto numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields != numFields0)
        {
          if (numFields == 0)
          {
            printWarning = true;
            cdo_warning("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(k));
          }
          else if (numFields0 == 0)
          {
            printWarning = true;
            cdo_warning("Inconsistent ensemble file, too few time steps in %s!", cdo_get_stream_name(0));
          }
          else
          {
            printError = true;
            cdo_warning("Inconsistent ensemble file, number of fields at time step %d of %s and %s differ!", tsID + 1,
                        cdo_get_stream_name(0), cdo_get_stream_name(k));
          }
          return;
        }
      }

      if (numFields0 > 0)
      {
        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID);
      }

      for (int fieldID = 0; fieldID < numFields0; ++fieldID)
      {
        auto taskNum = fieldNum % numTasks;
        auto &fields = fieldVector[taskNum];

        int varID = -1, levelID = -1;
        for (int k = 0; k < numFiles; ++k)
        {
          std::tie(varID, levelID) = cdo_inq_field(ensFileList[k].streamID);
          fields[k].init(ensFileList[k].varList.vars[varID]);
        }
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for (int k = 0; k < numFiles; ++k) { cdo_read_field(ensFileList[k].streamID, fields[k]); }

        if (Options::CDO_task) workerThread->wait();

        std::function<void()> ensstat_task = std::bind(ensstat, std::cref(ensFileList), std::ref(fields), streamID2, varID, levelID,
                                                       std::ref(workFields), std::ref(array2), std::ref(count2), operfunc, pn);

        Options::CDO_task ? workerThread->doAsync(ensstat_task) : ensstat_task();

        fieldNum++;
      }

      if (Options::CDO_task) workerThread->wait();

      tsID++;
    } while (numFields0 > 0);
  }

  void
  close() override
  {
    if (printWarning) cdo_warning("Inconsistent ensemble, processed only the first %d timesteps!", tsID);
    if (printError) cdo_abort("Inconsistent ensemble, processed only the first %d timesteps!", tsID);

    for (auto &ensFile : ensFileList) { cdo_stream_close(ensFile.streamID); }

    cdo_stream_close(streamID2);
  }
};
