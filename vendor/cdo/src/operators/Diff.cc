/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Diff       diff            Compare two datasets
*/

#include <map>
#include <algorithm>
#include <climits>

#include <cdi.h>

#include "workerthread.h"
#include "process_int.h"
#include "mpmo_color.h"
#include "cdo_math.h"
#include "cdo_options.h"
#include "printinfo.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "field_functions.h"
#include "pmlist.h"
#include "progress.h"

namespace
{
struct DiffResult
{
  size_t nvals{};
  size_t ndiff{};
  double absm{};
  double relm{};
  bool dsgn{};
  bool zero{};
};
}  // namespace

static inline void
diff_kernel(double v1, double v2, DiffResult &result)
{
  auto absdiff = std::fabs(v1 - v2);
  if (absdiff > 0.0) result.ndiff++;

  result.absm = std::max(result.absm, absdiff);

  auto vv = v1 * v2;
  if (vv < 0.0)
    result.dsgn = true;
  else if (is_equal(vv, 0.0))
    result.zero = true;
  else
    result.relm = std::max(result.relm, absdiff / std::max(std::fabs(v1), std::fabs(v2)));
}

static void
diff_kernel_mv(double v1, double v2, double missval1, double missval2, DiffResult &result)
{
  auto v1IsNan = std::isnan(v1);
  auto v2IsNan = std::isnan(v2);
  auto v1IsMissval = fp_is_equal(v1, missval1);
  auto v2IsMissval = fp_is_equal(v2, missval2);
  if (v1IsNan != v2IsNan)
  {
    result.ndiff++;
    result.relm = 1.0;
  }
  else if (!v1IsMissval && !v2IsMissval) { diff_kernel(v1, v2, result); }
  else if (v1IsMissval != v2IsMissval)
  {
    result.ndiff++;
    result.relm = 1.0;
  }
}

static DiffResult
diff(size_t n, Field const &field1, Field const &field2)
{
  DiffResult diffParam;
  auto hasMissvals = (field1.numMissVals || field2.numMissVals);
  if (hasMissvals)
  {
    auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2)
    {
      for (size_t i = 0; i < n; ++i) { diff_kernel_mv(v1[i], v2[i], mv1, mv2, diffParam); }
    };
    field_operation2(func, field1, field2, field1.missval, field2.missval);
  }
  else
  {
    auto func = [&](auto const &v1, auto const &v2)
    {
      for (size_t i = 0; i < n; ++i) { diff_kernel(v1[i], v2[i], diffParam); }
    };
    field_operation2(func, field1, field2);
  }

  return diffParam;
}

static inline void
diff_kernel2(double v1, double v2, DiffResult &result)
{
  auto absdiff = std::fabs(v1 - v2);
  if (absdiff > 0.0) result.ndiff++;

  auto vv = v1 * v2;
  if (vv < 0.0) { result.dsgn = true; }
  else if (is_equal(vv, 0.0)) { result.zero = true; }
  else
  {
    auto error = std::fabs(absdiff / v1);
    result.absm = std::max(result.absm, error);
    result.relm += error;
    result.nvals++;
  }
}

static DiffResult
diff2(size_t n, Field const &field1, Field const &field2)
{
  DiffResult diffParam;
  auto hasMissvals = (field1.numMissVals || field2.numMissVals);
  if (hasMissvals)
  {
    /*
    auto func = [&](auto &v1, auto &v2, double mv1, double mv2) {
      for (size_t i = 0; i < n; ++i) { diff_kernel_mv(v1[i], v2[i], mv1, mv2, diffParam); }
    };
    field_operation2(func, field1, field2, field1.missval, field2.missval);
    */
  }
  else
  {
    auto func = [&](auto const &v1, auto const &v2)
    {
      for (size_t i = 0; i < n; ++i) { diff_kernel2(v1[i], v2[i], diffParam); }
    };
    field_operation2(func, field1, field2);
    if (diffParam.nvals > 0) diffParam.relm /= diffParam.nvals;
  }

  return diffParam;
}

static void
use_real_part(Field &field)
{
  auto func = [](auto &v)
  {
    auto n = v.size() / 2;
    for (size_t i = 0; i < n; ++i) v[i] = v[i * 2];
  };

  return field_operation(func, field);
}

namespace
{
struct Parameter
{
  double absLimit{ 0.0 };
  double absLimit2{ 1.e-3 };
  double relLimit{ 1.0 };
  MapFlag mapFlag{ MapFlag::Undefined };
  int maxDiffFields{ INT_MAX };
  int numDiffFields{ 0 };
  int numDiffFields2{ 0 };
  bool printHeader{ true };
};
}  // namespace

static Parameter
get_parameter(void)
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "abslim")   params.absLimit = parameter_to_double(value);
      else if (key == "abslim2")  params.absLimit2 = parameter_to_double(value);
      else if (key == "rellim")   params.relLimit = parameter_to_double(value);
      else if (key == "maxcount") params.maxDiffFields = parameter_to_int(value);
      else if (key == "names")
      {
        if      (value == "left")      params.mapFlag = MapFlag::Left;
        else if (value == "right")     params.mapFlag = MapFlag::Right;
        else if (value == "intersect") params.mapFlag = MapFlag::Intersect;
        else cdo_abort("Invalid value for key >%s< (names=<left/right/intersect>)", key, value);
      }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
print_header(int operfunc, Parameter &params)
{
  params.printHeader = false;

  std::fprintf(stdout, "               Date     Time   Level Gridsize    Miss ");
  std::fprintf(stdout, "   Diff ");
  std::fprintf(stdout, ": S Z  Max_Absdiff Max_Reldiff : ");
  // clang-format off
  if      (operfunc == Func_Name)  std::fprintf(stdout, "Parameter name");
  else if (operfunc == Func_Param) std::fprintf(stdout, "Parameter ID");
  else if (operfunc == Func_Code)  std::fprintf(stdout, "Code number");
  // clang-format on
  std::fprintf(stdout, "\n");
}

static void
print_header2(int operfunc, Parameter &params)
{
  params.printHeader = false;

  std::fprintf(stdout, "               Date     Time   Level Gridsize    Miss ");
  std::fprintf(stdout, "   Diff ");
  std::fprintf(stdout, ": S Z  Max_Error   Mean_Error  : ");
  // clang-format off
  if      (operfunc == Func_Name)  std::fprintf(stdout, "Parameter name");
  else if (operfunc == Func_Param) std::fprintf(stdout, "Parameter ID");
  else if (operfunc == Func_Code)  std::fprintf(stdout, "Code number");
  // clang-format on
  std::fprintf(stdout, "\n");
}

void
print_diff(Field const &field1, Field const &field2, int fieldNumber, CdoVar const &var1, int levelID, int operfunc,
           CdiDateTime const &vDateTime, Parameter &params, int operfunc2, const DiffResult &dr)
{
  if (params.printHeader) (operfunc2 == 0) ? print_header(operfunc, params) : print_header2(operfunc, params);

  std::fprintf(stdout, "%6d ", fieldNumber);
  std::fprintf(stdout, ":");

  auto vdateString = date_to_string(vDateTime.date);
  auto vtimeString = time_to_string(vDateTime.time);

  set_text_color(stdout, MAGENTA);
  std::fprintf(stdout, "%s %s ", vdateString.c_str(), vtimeString.c_str());
  reset_text_color(stdout);
  set_text_color(stdout, GREEN);
  std::fprintf(stdout, "%7g ", cdo_zaxis_inq_level(var1.zaxisID, levelID));
  std::fprintf(stdout, "%8zu %7zu ", var1.gridsize, std::max(field1.numMissVals, field2.numMissVals));
  std::fprintf(stdout, "%7zu ", dr.ndiff);
  reset_text_color(stdout);

  std::fprintf(stdout, ":");
  std::fprintf(stdout, " %c %c ", dr.dsgn ? 'T' : 'F', dr.zero ? 'T' : 'F');
  set_text_color(stdout, BLUE);
  std::fprintf(stdout, "%#12.5g%#12.5g", dr.absm, dr.relm);
  reset_text_color(stdout);
  std::fprintf(stdout, " : ");

  char paramstr[32];
  if (operfunc == Func_Param) cdiParamToString(var1.param, paramstr, sizeof(paramstr));

  set_text_color(stdout, BRIGHT, GREEN);
  // clang-format off
  if      (operfunc == Func_Name)  std::fprintf(stdout, "%-11s", var1.name.c_str());
  else if (operfunc == Func_Param) std::fprintf(stdout, "%-11s", paramstr);
  else if (operfunc == Func_Code)  std::fprintf(stdout, "%4d", var1.code);
  // clang-format on
  reset_text_color(stdout);

  std::fprintf(stdout, "\n");
}

static void
compare_fields(Field &field1, Field &field2, int fieldNumber, CdoVar &var1, int levelID, int operfunc, CdiDateTime vDateTime,
               Parameter &params, int operfunc2, cdo::Progress &progress)
{
  size_t numNANs1 = (Options::fast || std::isnan(field1.missval)) ? 0 : field_num_NANs(field1);
  var1.counter += numNANs1;
  if (numNANs1 && field1.numMissVals == 0)
  {
    field1.missval = cdo::NaN();
    field1.numMissVals = numNANs1;
  }
  size_t numNANs2 = (Options::fast || std::isnan(field2.missval)) ? 0 : field_num_NANs(field2);
  var1.counter += numNANs2;
  if (numNANs2 && field2.numMissVals == 0)
  {
    field2.missval = cdo::NaN();
    field2.numMissVals = numNANs2;
  }

  auto dr = (operfunc2 == 0) ? diff(var1.gridsize, field1, field2) : diff2(var1.gridsize, field1, field2);

  auto checkRelativeLimit = true;
  if (!Options::silentMode || Options::cdoVerbose)
  {
    if (dr.absm > params.absLimit || (checkRelativeLimit && dr.relm >= params.relLimit) || Options::cdoVerbose)
    {
      progress.update(1);
      print_diff(field1, field2, fieldNumber, var1, levelID, operfunc, vDateTime, params, operfunc2, dr);
    }
  }

  if (dr.absm > params.absLimit || (checkRelativeLimit && dr.relm >= params.relLimit)) params.numDiffFields++;
  if (dr.absm > params.absLimit2 || (checkRelativeLimit && dr.relm >= params.relLimit)) params.numDiffFields2++;
}

class Diff : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Diff",
    // clang-format off
    .operators = { { "diff",     Func_Param, 0, DiffHelp },
                   { "diffp",    Func_Param, 0, DiffHelp },
                   { "diffn",    Func_Name,  0, DiffHelp },
                   { "diffc",    Func_Code,  0, DiffHelp },
                   { "difftest", Func_Name,  1, DiffHelp } },
    // clang-format on
    .aliases = { { "diffv", "diffn" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 2, 0, NoRestriction },
  };
  inline static RegisterEntry<Diff> registration = RegisterEntry<Diff>();

private:
  int operfunc{};
  int operfunc2{};

  Parameter params{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID{};

  std::map<int, int> mapOfVarIDs{};

  VarList varList1{};
  VarList varList2{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    operfunc2 = cdo_operator_f2(operatorID);

    params = get_parameter();

    constexpr double rangeMin = -1.e33;
    constexpr double rangeMax = 1.e33;
    if (params.relLimit < rangeMin || params.relLimit > rangeMax) cdo_abort("Rel. limit out of range!");
    if (params.absLimit < rangeMin || params.absLimit > rangeMax) cdo_abort("Abs. limit out of range!");
    if (params.absLimit2 < rangeMin || params.absLimit2 > rangeMax) cdo_abort("Abs2. limit out of range!");

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    if (params.mapFlag == MapFlag::Undefined)
    {
      varList_compare(varList1, varList2);
      for (auto const &var : varList1.vars) mapOfVarIDs[var.ID] = var.ID;
    }
    else { varList_map(varList1, varList2, params.mapFlag, mapOfVarIDs); }

    taxisID = vlistInqTaxis(vlistID1);
  }

  void
  run() override
  {
    auto runAsync = (Options::CDO_Async_Read > 0);
    auto workerThread = runAsync ? std::make_unique<WorkerThread>() : nullptr;
    auto numTasks = runAsync ? 2 : 1;

    FieldVector fieldVector1(numTasks);
    FieldVector fieldVector2(numTasks);

    auto numSteps1 = varList1.numSteps();
    cdo::Progress progress(get_id());

    int numSets = 0;
    int numFields, numFields2;
    int tsID = 0;
    while (true)
    {
      auto stopRead = false;

      numFields = cdo_stream_inq_timestep(streamID1, tsID);
      auto vDateTime = taxisInqVdatetime(taxisID);

      numFields2 = cdo_stream_inq_timestep(streamID2, tsID);

      if (numFields == 0 || numFields2 == 0) break;

      if (numSteps1 > 1 && params.numDiffFields == 0) progress.update((tsID + 1.0) / numSteps1);

      int fieldID2next = 0;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID1, levelID] = cdo_inq_field(streamID1);

        auto it = mapOfVarIDs.find(varID1);
        if (it == mapOfVarIDs.end())
        {
          if (params.mapFlag == MapFlag::Right || params.mapFlag == MapFlag::Intersect) continue;
          cdo_abort("Internal problem (tsID=%d fieldID=%d): varID1=%d not found!", tsID + 1, fieldID + 1, varID1);
        }

        int varID2 = 0;
        for (; fieldID2next < numFields2; ++fieldID2next)
        {
          auto [varID2x, levelID2x] = cdo_inq_field(streamID2);
          varID2 = varID2x;
          if (it->second == varID2)
          {
            ++fieldID2next;
            break;
          }
        }

        if (it->second != varID2 && fieldID2next == numFields2)
          cdo_abort("Internal problem (tsID=%d fieldID=%d): varID2=%d not found in second stream!", tsID + 1, fieldID + 1,
                    it->second);

        auto &var1 = varList1.vars[varID1];
        auto &var2 = varList2.vars[varID2];

        auto taskNum = numSets % numTasks;
        auto &field1 = fieldVector1[taskNum];
        auto &field2 = fieldVector2[taskNum];

        field1.init(var1);
        cdo_read_field(streamID1, field1);
        if (var1.nwpv == CDI_COMP) use_real_part(field1);

        field2.init(var2);
        cdo_read_field(streamID2, field2);
        if (var2.nwpv == CDI_COMP) use_real_part(field2);

        if (runAsync && numSets > 0)
        {
          workerThread->wait();
          // clang-format off
          if (params.numDiffFields >= params.maxDiffFields) { stopRead = true; break; }
          // clang-format on
        }

        std::function<void()> compare_fields_task
            = std::bind(compare_fields, std::ref(field1), std::ref(field2), numSets + 1, std::ref(var1), levelID, operfunc,
                        vDateTime, std::ref(params), operfunc2, std::ref(progress));
        runAsync ? workerThread->doAsync(compare_fields_task) : compare_fields_task();

        if (not runAsync)
        {
          // clang-format off
          if (params.numDiffFields >= params.maxDiffFields) { stopRead = true; break; }
          // clang-format on
        }

        numSets++;
      }

      if (stopRead) break;

      tsID++;
    }

    if (runAsync) workerThread->wait();

    if (params.numDiffFields > 0)
    {
      Options::cdoExitStatus = 1;

      set_text_color(stdout, BRIGHT, RED);
      std::fprintf(stdout, "  %d of %d fields differ", params.numDiffFields, numSets);
      reset_text_color(stdout);
      std::fprintf(stdout, "\n");

      if (params.numDiffFields != params.numDiffFields2 && params.absLimit < params.absLimit2)
        std::fprintf(stdout, "  %d of %d fields differ more than %g\n", params.numDiffFields2, numSets, params.absLimit2);
      //  std::fprintf(stdout, "  %d of %d fields differ more then one thousandth\n", nprec, ngrec);
    }

    if (numFields == 0 && numFields2 > 0) cdo_warning("stream2 has more time steps than stream1!");
    if (numFields > 0 && numFields2 == 0) cdo_warning("stream1 has more time steps than stream2!");

    for (auto const &var : varList1.vars)
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
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
