/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Timstat    timrange        Time range
      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Normalize by (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Normalize by (n-1)]
      Hourstat   hourrange       Hourly range
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Normalize by (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Normalize by (n-1)]
      Daystat    dayrange        Daily range
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Normalize by (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Normalize by (n-1)]
      Monstat    monrange        Monthly range
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Normalize by (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Normalize by (n-1)]
      Yearstat   yearrange       Yearly range
      Yearstat   yearmin         Yearly minimum
      Yearstat   yearmax         Yearly maximum
      Yearstat   yearsum         Yearly sum
      Yearstat   yearmean        Yearly mean
      Yearstat   yearavg         Yearly average
      Yearstat   yearvar         Yearly variance
      Yearstat   yearvar1        Yearly variance [Normalize by (n-1)]
      Yearstat   yearstd         Yearly standard deviation
      Yearstat   yearstd1        Yearly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include <unordered_map>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "workerthread.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "util_date.h"
#include "progress.h"
#include "field_functions.h"
#include "param_conversion.h"
#include "pmlist.h"

static void
vlist_set_frequency(int vlistID, std::string const &frequency)
{
  if (frequency.size()) cdiDefAttTxt(vlistID, CDI_GLOBAL, "frequency", frequency.size(), frequency.c_str());
}

static void
vlist_set_frequency(int vlistID, int compareDate)
{
  // clang-format off
  if      (compareDate == CMP_DAY)   vlist_set_frequency(vlistID, "day");
  else if (compareDate == CMP_MONTH) vlist_set_frequency(vlistID, "mon");
  else if (compareDate == CMP_YEAR)  vlist_set_frequency(vlistID, "year");
  // clang-format on
}

static void
get_parameter(double &vfraction, bool &completeOnly)
{
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
      if      (key == "vfraction")       vfraction = parameter_to_double(value);
      else if (key == "complete_only")   completeOnly = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Timstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timstat",
    // clang-format off
    .operators = { { "timrange", FieldFunc_Range, CMP_DATE, TimstatHelp },
                   { "timmin", FieldFunc_Min, CMP_DATE, TimstatHelp },
                   { "timmax", FieldFunc_Max, CMP_DATE, TimstatHelp },
                   { "timsum", FieldFunc_Sum, CMP_DATE, TimstatHelp },
                   { "timmean", FieldFunc_Mean, CMP_DATE, TimstatHelp },
                   { "timavg", FieldFunc_Avg, CMP_DATE, TimstatHelp },
                   { "timvar", FieldFunc_Var, CMP_DATE, TimstatHelp },
                   { "timvar1", FieldFunc_Var1, CMP_DATE, TimstatHelp },
                   { "timstd", FieldFunc_Std, CMP_DATE, TimstatHelp },
                   { "timstd1", FieldFunc_Std1, CMP_DATE, TimstatHelp },
                   { "timminidx", FieldFunc_Minidx, CMP_DATE, TimstatHelp },
                   { "timmaxidx", FieldFunc_Maxidx, CMP_DATE, TimstatHelp },
                   { "yearrange", FieldFunc_Range, CMP_YEAR, YearstatHelp },
                   { "yearmin", FieldFunc_Min, CMP_YEAR, YearstatHelp },
                   { "yearmax", FieldFunc_Max, CMP_YEAR, YearstatHelp },
                   { "yearsum", FieldFunc_Sum, CMP_YEAR, YearstatHelp },
                   { "yearmean", FieldFunc_Mean, CMP_YEAR, YearstatHelp },
                   { "yearavg", FieldFunc_Avg, CMP_YEAR, YearstatHelp },
                   { "yearvar", FieldFunc_Var, CMP_YEAR, YearstatHelp },
                   { "yearvar1", FieldFunc_Var1, CMP_YEAR, YearstatHelp },
                   { "yearstd", FieldFunc_Std, CMP_YEAR, YearstatHelp },
                   { "yearstd1", FieldFunc_Std1, CMP_YEAR, YearstatHelp },
                   { "yearminidx", FieldFunc_Minidx, CMP_YEAR, YearstatHelp },
                   { "yearmaxidx", FieldFunc_Maxidx, CMP_YEAR, YearstatHelp },
                   { "monrange", FieldFunc_Range, CMP_MONTH, MonstatHelp },
                   { "monmin", FieldFunc_Min, CMP_MONTH, MonstatHelp },
                   { "monmax", FieldFunc_Max, CMP_MONTH, MonstatHelp },
                   { "monsum", FieldFunc_Sum, CMP_MONTH, MonstatHelp },
                   { "monmean", FieldFunc_Mean, CMP_MONTH, MonstatHelp },
                   { "monavg", FieldFunc_Avg, CMP_MONTH, MonstatHelp },
                   { "monvar", FieldFunc_Var, CMP_MONTH, MonstatHelp },
                   { "monvar1", FieldFunc_Var1, CMP_MONTH, MonstatHelp },
                   { "monstd", FieldFunc_Std, CMP_MONTH, MonstatHelp },
                   { "monstd1", FieldFunc_Std1, CMP_MONTH, MonstatHelp },
                   { "dayrange", FieldFunc_Range, CMP_DAY, DaystatHelp },
                   { "daymin", FieldFunc_Min, CMP_DAY, DaystatHelp },
                   { "daymax", FieldFunc_Max, CMP_DAY, DaystatHelp },
                   { "daysum", FieldFunc_Sum, CMP_DAY, DaystatHelp },
                   { "daymean", FieldFunc_Mean, CMP_DAY, DaystatHelp },
                   { "dayavg", FieldFunc_Avg, CMP_DAY, DaystatHelp },
                   { "dayvar", FieldFunc_Var, CMP_DAY, DaystatHelp },
                   { "dayvar1", FieldFunc_Var1, CMP_DAY, DaystatHelp },
                   { "daystd", FieldFunc_Std, CMP_DAY, DaystatHelp },
                   { "daystd1", FieldFunc_Std1, CMP_DAY, DaystatHelp },
                   { "hourrange", FieldFunc_Range, CMP_HOUR, HourstatHelp },
                   { "hourmin", FieldFunc_Min, CMP_HOUR, HourstatHelp },
                   { "hourmax", FieldFunc_Max, CMP_HOUR, HourstatHelp },
                   { "hoursum", FieldFunc_Sum, CMP_HOUR, HourstatHelp },
                   { "hourmean", FieldFunc_Mean, CMP_HOUR, HourstatHelp },
                   { "houravg", FieldFunc_Avg, CMP_HOUR, HourstatHelp },
                   { "hourvar", FieldFunc_Var, CMP_HOUR, HourstatHelp },
                   { "hourvar1", FieldFunc_Var1, CMP_HOUR, HourstatHelp },
                   { "hourstd", FieldFunc_Std, CMP_HOUR, HourstatHelp },
                   { "hourstd1", FieldFunc_Std1, CMP_HOUR, HourstatHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timstat> registration = RegisterEntry<Timstat>();

private:
  static const TimeStat timestatDate{ TimeStat::MEAN };

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID3{};
  int taxisID3{ -1 };

  int compareDate{};

  bool completeOnly{ false };
  bool handleVfraction{ false };
  double vfraction{ -1.0 };

  cdo::StepStat2D stepStat{};

  std::vector<FieldInfo> fieldInfoList{};
  DateTimeList dtlist{};
  VarList varList1{};

  void
  create_diag_stream(int operatorID, int vlistID1, int numVars)
  {
    char fileName[8192];
    std::strcpy(fileName, cdo_operator_name(operatorID));
    std::strcat(fileName, "_");
    std::strcat(fileName, cdo_get_stream_name(1));
    streamID3 = open_write(fileName);

    vlistID3 = vlistDuplicate(vlistID1);

    for (int varID = 0; varID < numVars; ++varID)
    {
      vlistDefVarDatatype(vlistID3, varID, CDI_DATATYPE_INT32);
      vlistDefVarMissval(vlistID3, varID, -1);
      cdiDefKeyString(vlistID3, varID, CDI_KEY_UNITS, "");
      cdiDeleteKey(vlistID3, varID, CDI_KEY_ADDOFFSET);
      cdiDeleteKey(vlistID3, varID, CDI_KEY_SCALEFACTOR);
    }

    taxisID3 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID3);
    vlistDefTaxis(vlistID3, taxisID3);

    cdo_def_vlist(streamID3, vlistID3);
  }

  bool
  check_numSets(std::vector<int> &numSetsList)
  {
    std::unordered_map<int, std::pair<double, std::string>> periodNameMap
        = { { CMP_DAY, { 23.0 / 24.0, "day" } }, { CMP_MONTH, { 28 / 31.0, "month" } }, { CMP_YEAR, { 365 / 366.0, "year" } } };

    if (numSetsList.size() > 1)
    {
      auto start = (numSetsList.size() > 2) ? 1 : 0;
      auto numSetsMin = std::min_element(numSetsList.begin() + start, numSetsList.end() - 1);
      auto [periodFactor, periodName] = periodNameMap[compareDate];
      if ((double) numSetsList.back() < *numSetsMin * periodFactor)
      {
        if (Options::cdoVerbose)
          cdo_warning("Last %s has less steps (%d) than all previous %ss (%d)!", periodName, numSetsList.back(), periodName,
                      *numSetsMin);
        return true;
      }
    }

    return false;
  }

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);
    compareDate = cdo_operator_f2(operatorID);

    stepStat.init(operfunc);

    get_parameter(vfraction, completeOnly);
    handleVfraction = (vfraction >= 0.0 && vfraction <= 1.0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!stepStat.lminmax) vlist_unpack(vlistID2);

    vlist_define_timestep_type(vlistID2, operfunc);

    vlistDefNtsteps(vlistID2, (compareDate == CMP_DATE) ? 1 : -1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    vlist_set_frequency(vlistID2, compareDate);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    if (Options::CDO_diagnostic) create_diag_stream(operatorID, vlistID1, varList1.numVars());

    fieldInfoList.resize(varList1.maxFields());

    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    // if ((Options::CDO_???(--single) == MemType::Float) && stepStat.lmean) VARS_MEMTYPE = FIELD_NAT;
    // if (Options::CDO_Memtype == MemType::Float) VARS_MEMTYPE = FIELD_FLT;
    if (Options::CDO_diagnostic || (handleVfraction && stepStat.lmean)) VARS_MEMTYPE = FIELD_DBL;

    stepStat.alloc(varList1, VARS_MEMTYPE);

    // for (auto &var1 : varList1.vars) var1.memType = stepStat.var1(var1.ID, 0).memType;
  }

  void
  run_sync()
  {
    std::vector<int> numSetsList;
    CdiDateTime vDateTime0{};
    CdiDateTime vDateTimeN{};
    Field field;

    auto numSteps1 = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numSets = 0;
      int numFields = 0;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        if (numSteps1 > 1) progress.update((tsID + 1.0) / numSteps1);

        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime = dtlist.vDateTime(numSets);

        if (numSets == 0) vDateTime0 = vDateTime;

        if (date_is_neq(vDateTime, vDateTime0, compareDate))
        {
          cdo_add_steps(-1);
          break;
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          stepStat.add_field(field, varID, levelID, numSets);
        }

        vDateTimeN = vDateTime;
        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      if (compareDate == CMP_DAY || compareDate == CMP_MONTH || compareDate == CMP_YEAR)
      {
        numSetsList.push_back(numSets);
        if (numFields == 0 && check_numSets(numSetsList) && completeOnly) break;
      }

      cdo::fields_process(fieldInfoList, varList1, stepStat, numSets);

      if (Options::cdoVerbose) cdo_print("%s  numSteps = %d", datetime_to_string(vDateTimeN), numSets);

      if (handleVfraction && stepStat.lmean) cdo::fields_set_missval(fieldInfoList, varList1, stepStat, numSets, vfraction);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo::write_out_stream(streamID2, fieldInfoList, varList1, stepStat, otsID);

      if (Options::CDO_diagnostic)
      {
        dtlist.stat_taxis_def_timestep(taxisID3, numSets);
        cdo::write_diag_stream(streamID3, fieldInfoList, varList1, stepStat, otsID, numSets);
      }

      if (numFields == 0) break;
      otsID++;
    }
  }

  static void
  add_fields2D(const FieldVector2D &fields2D, std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1,
               cdo::StepStat2D &stepStat, int numSets) noexcept
  {
    for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      if (varList1.vars[varID].isConstant) continue;

      stepStat.add_field(fields2D[varID][levelID], varID, levelID, numSets);
    }
  }

  void
  run_async()
  {
    std::vector<int> numSetsList;
    CdiDateTime vDateTime0{};
    CdiDateTime vDateTimeN{};
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    FieldVector3D fields3D(2);
    field2D_init(fields3D[0], varList1, FIELD_VEC | FIELD_NAT);
    field2D_init(fields3D[1], varList1, FIELD_VEC | FIELD_NAT);

    auto useTask = true;
    auto workerThread = useTask ? std::make_unique<WorkerThread>() : nullptr;

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      int numSets = 0;
      int numFields = 0;
      while (true)
      {
        numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

        dtlist.taxis_inq_timestep(taxisID1, numSets);
        auto vDateTime = dtlist.vDateTime(numSets);

        if (numSets == 0) vDateTime0 = vDateTime;

        if (date_is_neq(vDateTime, vDateTime0, compareDate))
        {
          cdo_add_steps(-1);
          break;
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
          cdo_read_field(streamID1, fields3D[numSets % 2][varID][levelID]);
        }

        if (useTask && numSets > 0) workerThread->wait();

        std::function<void()> add_fields2D_task
            = std::bind(add_fields2D, std::cref(fields3D[numSets % 2]), std::cref(fieldInfoList), std::cref(varList1),
                        std::ref(stepStat), numSets);

        useTask ? workerThread->doAsync(add_fields2D_task) : add_fields2D_task();

        vDateTimeN = vDateTime;
        numSets++;
        tsID++;
      }

      if (numFields == 0 && numSets == 0) break;

      if (useTask) workerThread->wait();

      if (compareDate == CMP_DAY || compareDate == CMP_MONTH || compareDate == CMP_YEAR)
      {
        numSetsList.push_back(numSets);
        if (numFields == 0 && check_numSets(numSetsList) && completeOnly) break;
      }

      cdo::fields_process(fieldInfoList, varList1, stepStat, numSets);

      if (Options::cdoVerbose) cdo_print("%s  numSteps = %d", datetime_to_string(vDateTimeN), numSets);

      if (handleVfraction && stepStat.lmean) cdo::fields_set_missval(fieldInfoList, varList1, stepStat, numSets, vfraction);

      dtlist.stat_taxis_def_timestep(taxisID2, numSets);
      cdo::write_out_stream(streamID2, fieldInfoList, varList1, stepStat, otsID);

      if (Options::CDO_diagnostic)
      {
        dtlist.stat_taxis_def_timestep(taxisID3, numSets);
        cdo::write_diag_stream(streamID3, fieldInfoList, varList1, stepStat, otsID, numSets);
      }

      if (numFields == 0) break;
      otsID++;
    }
  }

  void
  run() override
  {
    Options::CDO_Async_Read ? run_async() : run_sync();
  }

  void
  close() override
  {
    if (Options::CDO_diagnostic) cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
