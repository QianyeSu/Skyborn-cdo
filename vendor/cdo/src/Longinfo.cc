/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Longinfo       linfo            Long dataset information
*/

#include <cfloat>

#include <cdi.h>

#include "workerthread.h"
#include "cdo_options.h"
#include "cdo_math.h"
#include "process_int.h"
#include "varray.h"
#include "printinfo.h"
#include "cdo_zaxis.h"
#include "field_functions.h"

struct LonginfoStat
{
  double min{ DBL_MAX };
  double max{ -DBL_MAX };
  double sum{ 0.0 };
  size_t nvals{ 0 };
};

static void
field_min_max_sum(Field const &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  auto func = [&](auto const &v, auto n) { return varray_min_max_sum(v, n, mms); };
  mms = field_operation(func, field, field.size);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
}

static size_t
field_min_max_sum_mv(Field const &field, double &min, double &max, double &sum)
{
  auto mms = MinMaxSum(min, max, sum);
  auto func = [&](auto const &v, auto n, double mv) { return varray_min_max_sum_mv(v, n, mms, mv); };
  mms = field_operation(func, field, field.size, field.missval);

  min = mms.min;
  max = mms.max;
  sum = mms.sum;
  return mms.n;
}

static size_t
compute_stat_real(Field const &field, LonginfoStat &infostat, size_t gridsize)
{
  size_t imiss = 0;

  if (field.numMissVals)
  {
    auto nvals = field_min_max_sum_mv(field, infostat.min, infostat.max, infostat.sum);
    imiss = gridsize - nvals;
    infostat.nvals += nvals;
  }
  else if (gridsize == 1)
  {
    infostat.sum = (infostat.nvals == 0) ? field[0] : infostat.sum + field[0];
    infostat.nvals += 1;
  }
  else
  {
    field_min_max_sum(field, infostat.min, infostat.max, infostat.sum);
    infostat.nvals += gridsize;
  }

  return imiss;
}

static void
long_info(Field &field, int tsID, CdiDateTime vDateTime, int numFields, int fieldID, int varID, int levelID, int vlistID,
          CdoVar &var)
{
  char paramstr[32];

  if (fieldID == 0)
  {
    fprintf(stdout, "timestep: %d\n", tsID + 1);
    fprintf(stdout, "\tdateTime: %s\n\n", datetime_to_string(vDateTime).c_str());
  }

  fprintf(stdout, "\tfield: %d of %d\n", fieldID + 1, numFields);

  auto dig = (var.dataType == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

  fprintf(stdout, "\t\tvarIndex: %d\n", varID + 1);
  fprintf(stdout, "\t\tlevelIndex: %d\n", levelID + 1);
  fprintf(stdout, "\t\tlevel: %.*g\n", dig, cdo_zaxis_inq_level(var.zaxisID, levelID));
  fprintf(stdout, "\t\tname: %s\n", var.name.c_str());
  if (var.longname.size()) fprintf(stdout, "\t\tlongname: \"%s\"\n", var.longname.c_str());
  if (var.units.size()) fprintf(stdout, "\t\tunits: \"%s\"\n", var.units.c_str());
  cdiParamToString(var.param, paramstr, sizeof(paramstr));
  if (paramstr[0] && paramstr[0] != '-') fprintf(stdout, "\t\tparam: %s\n", paramstr);

  size_t numNANs = 0;
  if (not Options::fast and not std::isnan(field.missval) and field.numMissVals == 0) { numNANs = field_num_NANs(field); }
  var.counter += numNANs;

  fprintf(stdout, "\t\tdataType: %s\n", cdo::datatype_to_cstr(var.dataType));
  fprintf(stdout, "\t\tmemoryType: %s\n", (var.memType == MemType::Float) ? "float" : "double");
  fprintf(stdout, "\t\tgridsize: %zu\n", var.gridsize);
  fprintf(stdout, "\t\tnumMiss: %zu\n", field.numMissVals);
  fprintf(stdout, "\t\tmissval: %.*g\n", dig, var.missval);

  double addoffset = 0.0, scalefactor = 1.0;
  auto haveAddoffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
  auto haveScalefactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);
  if (haveAddoffset) fprintf(stdout, "\t\taddoffset: %.*g\n", dig, addoffset);
  if (haveScalefactor) fprintf(stdout, "\t\tscalefactor: %.*g\n", dig, scalefactor);

  if (numNANs)
  {
    field.missval = cdo::NaN();
    field.numMissVals = numNANs;
  }

  LonginfoStat infostat;
  auto imiss = compute_stat_real(field, infostat, var.gridsize);
  (void) imiss;

  if (infostat.nvals > 1)
  {
    fprintf(stdout, "\t\trange: %.*g\n", dig, infostat.max - infostat.min);
    fprintf(stdout, "\t\tminimum: %.*g\n", dig, infostat.min);
    fprintf(stdout, "\t\tmaximum: %.*g\n", dig, infostat.max);
    fprintf(stdout, "\t\taverage: %.*g\n", dig, infostat.sum / infostat.nvals);
    // fprintf(stdout, "\t\tmedian: %.*g\n", dig, field_median(field));
    fprintf(stdout, "\t\tstandardDev: %.*g\n", dig, field_std1(field));
    fprintf(stdout, "\t\tskewness: %.*g\n", dig, field_skew(field));
    fprintf(stdout, "\t\tkurtosis: %.*g\n", dig, field_kurt(field));
  }
  else if (infostat.nvals == 1) { fprintf(stdout, "\t\tvalue: %g\n", infostat.sum); }

  cdo_print_attributes(stdout, vlistID, varID, 16);
}

class Longinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Longinfo",
    .operators = { { "linfo" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Longinfo> registration = RegisterEntry<Longinfo>(module);

  CdoStreamID streamID{};

  int taxisID{};
  int vlistID{};

  VarList varList{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    operator_check_argc(0);

    streamID = cdo_open_read(0);
    vlistID = cdo_stream_inq_vlist(streamID);
    taxisID = vlistInqTaxis(vlistID);

    varList = VarList(vlistID);
  }

  void
  run() override
  {
    auto runAsync = (Options::CDO_Async_Read > 0);
    auto workerThread = runAsync ? std::make_unique<WorkerThread>() : nullptr;
    auto numTasks = runAsync ? 2 : 1;

    FieldVector fieldVector(numTasks);

    int numSets = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;
      auto vDateTime = taxisInqVdatetime(taxisID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID);
        auto &var = varList.vars[varID];
        auto taskNum = numSets % numTasks;
        auto &field = fieldVector[taskNum];
        field.init(var);
        cdo_read_field(streamID, field);

        if (runAsync && numSets > 0) { workerThread->wait(); }

        std::function<void()> long_info_task
            = std::bind(long_info, std::ref(field), tsID, vDateTime, numFields, fieldID, varID, levelID, vlistID, std::ref(var));

        runAsync ? workerThread->doAsync(long_info_task) : long_info_task();

        numSets++;
      }

      // if (imiss != numMissVals && numMissVals) cdo_warning("Found %zu of %zu missing values!", imiss, numMissVals);

      tsID++;
    }

    if (runAsync) workerThread->wait();

    for (auto const &var : varList.vars)
    {
      if (var.counter > 0)
      {
        cdo_warning("%s contains %zu NaNs which are not treated as missing values. "
                    "This can lead to incorrect CDO results in all other arithmetic functions!",
                    var.name, var.counter);
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
