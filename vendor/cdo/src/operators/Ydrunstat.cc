/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida
          Fabian Wachsmann

*/

/*
   This module contains the following operators:

      Ydrunstat    ydrunmin          Multi-year daily running minimum
      Ydrunstat    ydrunmax          Multi-year daily running maximum
      Ydrunstat    ydrunsum          Multi-year daily running sum
      Ydrunstat    ydrunmean         Multi-year daily running mean
      Ydrunstat    ydrunavg          Multi-year daily running average
      Ydrunstat    ydrunvar          Multi-year daily running variance
      Ydrunstat    ydrunvar1         Multi-year daily running variance [Normalize by (n-1)]
      Ydrunstat    ydrunstd          Multi-year daily running standard deviation
      Ydrunstat    ydrunstd1         Multi-year daily running standard deviation [Normalize by (n-1)]
*/

#include "cdi.h"
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "util_string.h"
#include "param_conversion.h"
#include "datetime.h"
#include "field_functions.h"
#include "pmlist.h"

constexpr int MaxDays = 373;

namespace
{
struct YdayStats
{
  int numSets[MaxDays]{};
  CdiDateTime vDateTime[MaxDays]{};
  FieldVector2D varsData1[MaxDays];
  FieldVector2D varsData2[MaxDays];
  int vlistID{};
  VarList varList{};

  explicit YdayStats(int _vlistID) : vlistID{ _vlistID }, varList{ VarList(_vlistID) } {}
};
}  // namespace

static void
ydstat_update(YdayStats &stats, CdiDateTime vDateTime, FieldVector2D const &varsData1, FieldVector2D const &varsData2, int numSets,
              int operfunc)
{
  auto lvarstd = (varsData2.size() > 0);

  auto dayOfYear = decode_day_of_year(vDateTime.date);
  if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

  stats.vDateTime[dayOfYear] = vDateTime;

  if (!stats.varsData1[dayOfYear].size())
  {
    field2D_init(stats.varsData1[dayOfYear], stats.varList, FIELD_VEC);
    if (lvarstd) field2D_init(stats.varsData2[dayOfYear], stats.varList, FIELD_VEC);
  }

  auto numVars = stats.varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var = stats.varList.vars[varID];
    if (var.timeType == TIME_CONSTANT) continue;

    for (int levelID = 0; levelID < var.nlevels; ++levelID)
    {
      auto const &varData1 = varsData1[varID][levelID];
      if (stats.numSets[dayOfYear] == 0)
      {
        field_copy(varData1, stats.varsData1[dayOfYear][varID][levelID]);
        if (lvarstd) field_copy(varsData2[varID][levelID], stats.varsData2[dayOfYear][varID][levelID]);
      }
      else if (lvarstd)
      {
        field2_sum(stats.varsData1[dayOfYear][varID][levelID], varData1);
        field2_sum(stats.varsData2[dayOfYear][varID][levelID], varsData2[varID][levelID]);
      }
      else { field2_function(stats.varsData1[dayOfYear][varID][levelID], varData1, operfunc); }
    }
  }

  stats.numSets[dayOfYear] += numSets;
}

static void
ydstat_finalize(YdayStats &stats, int operfunc)
{
  auto lmean = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
  auto lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
  auto lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
  int divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);

  auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

  for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
    if (stats.numSets[dayOfYear])
    {
      auto numVars = stats.varList.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = stats.varList.vars[varID];
        if (var.timeType == TIME_CONSTANT) continue;

        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto numSets = stats.numSets[dayOfYear];
          auto &rvars1 = stats.varsData1[dayOfYear][varID][levelID];

          if (lmean) { fieldc_div(rvars1, (double) numSets); }
          else if (lvarstd)
          {
            auto const &rvars2 = stats.varsData2[dayOfYear][varID][levelID];
            fieldc_stdvar_func(rvars1, rvars2, numSets, divisor);
          }
        }
      }
    }
}

namespace
{
struct Parameter
{
  int nts{ -1 };  // number of timesteps
  char rm{ 0 };   // Read method (circular)
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs < 1) cdo_abort("Too few arguments!");

  auto argList = cdo_get_oper_argv();
  auto param1 = argList[0];
  if (std::isdigit(param1[0]) && !string_contains(param1, '='))
  {
    params.nts = parameter_to_int(param1);
    argList.erase(argList.begin());
    numArgs--;
  }

  if (numArgs)
  {
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
      if      (key == "nts") params.nts = parameter_to_int(value);
      else if (key == "rm")  params.rm = value[0];
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
check_parameter(const Parameter &parameter)
{
  if (parameter.nts == -1) cdo_abort("Too few parameter!");
  if (parameter.nts <= 0) cdo_abort("Parameter nts must be greater than 0!");
  if (parameter.rm != 0 && parameter.rm != 'c') cdo_abort("Parameter rm must only contain 'c'!");
}

class Ydrunstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ydrunstat",
    .operators = { { "ydrunmin", FieldFunc_Min, 0, YdrunstatHelp },
                   { "ydrunmax", FieldFunc_Max, 0, YdrunstatHelp },
                   { "ydrunsum", FieldFunc_Sum, 0, YdrunstatHelp },
                   { "ydrunmean", FieldFunc_Mean, 0, YdrunstatHelp },
                   { "ydrunavg", FieldFunc_Avg, 0, YdrunstatHelp },
                   { "ydrunstd", FieldFunc_Std, 0, YdrunstatHelp },
                   { "ydrunstd1", FieldFunc_Std1, 0, YdrunstatHelp },
                   { "ydrunvar", FieldFunc_Var, 0, YdrunstatHelp },
                   { "ydrunvar1", FieldFunc_Var1, 0, YdrunstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Ydrunstat> registration = RegisterEntry<Ydrunstat>();

private:
  int operfunc{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };

  char readMethod{ 0 };
  bool lvarstd{};
  int numDates{ -1 };
  int dpy{};

  VarList varList1{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_input_arg("number of timesteps");

    auto params = get_parameter();
    check_parameter(params);
    numDates = params.nts;
    readMethod = params.rm;
    if (Options::cdoVerbose) cdo_print("numDates=%d readMethod=%c", numDates, readMethod);

    auto lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);
    lvarstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Var || operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);
    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    if (!lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    dpy = calendar_dpy(taxisInqCalendar(taxisID1));

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector3D varsData1(numDates + 1);
    FieldVector3D varsData2(numDates + 1);
    for (int its = 0; its < numDates; its++)
    {
      field2D_init(varsData1[its], varList1, FIELD_VEC);
      if (lvarstd) field2D_init(varsData2[its], varList1, FIELD_VEC);
    }

    YdayStats stats = YdayStats(vlistID1);
    std::vector<CdiDateTime> cdiDateTimes(numDates + 1);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    int startYear = 0;
    int tsID = 0;

    for (tsID = 0; tsID < numDates; ++tsID)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) cdo_abort("File has less then %d timesteps!", numDates);

      cdiDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      if (tsID == 0 && readMethod == 'c') startYear = cdiDateTimes[tsID].date.year;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
        auto &rvars1 = varsData1[tsID][varID][levelID];
        cdo_read_field(streamID1, rvars1);

        if (lvarstd)
        {
          field2_moq(varsData2[tsID][varID][levelID], rvars1);
          for (int inp = 0; inp < tsID; ++inp)
            field2_sumsumq(varsData1[inp][varID][levelID], varsData2[inp][varID][levelID], rvars1);
        }
        else
        {
          for (int inp = 0; inp < tsID; ++inp) field2_function(varsData1[inp][varID][levelID], rvars1, operfunc);
        }
      }
    }

    while (true)
    {
      cdiDateTimes[numDates] = datetime_avg(dpy, numDates, cdiDateTimes);

      ydstat_update(stats, cdiDateTimes[numDates], varsData1[0], varsData2[0], numDates, operfunc);

      cdiDateTimes[numDates] = cdiDateTimes[0];
      varsData1[numDates] = varsData1[0];
      if (lvarstd) varsData2[numDates] = varsData2[0];

      for (int inp = 0; inp < numDates; ++inp)
      {
        cdiDateTimes[inp] = cdiDateTimes[inp + 1];
        varsData1[inp] = varsData1[inp + 1];
        if (lvarstd) varsData2[inp] = varsData2[inp + 1];
      }

      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdiDateTimes[numDates - 1] = taxisInqVdatetime(taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &rvars1 = varsData1[numDates - 1][varID][levelID];
        cdo_read_field(streamID1, rvars1);

        if (lvarstd)
        {
          field2_moq(varsData2[numDates - 1][varID][levelID], rvars1);
          for (int inp = 0; inp < numDates - 1; ++inp)
            field2_sumsumq(varsData1[inp][varID][levelID], varsData2[inp][varID][levelID], rvars1);
        }
        else
        {
          for (int inp = 0; inp < numDates - 1; ++inp) field2_function(varsData1[inp][varID][levelID], rvars1, operfunc);
        }
      }

      tsID++;
    }

    cdo_stream_close(streamID1);

    if (readMethod == 'c')
    {
      if (cdo_assert_files_only() == false) cdo_warning("Operators cannot be piped in circular mode");

      auto endYear = cdiDateTimes[numDates - 1].date.year;
      auto cdiStream = streamOpenRead(cdo_get_stream_name(0));
      auto cdiVlistID = streamInqVlist(cdiStream);
      auto cdiTaxisID = vlistInqTaxis(cdiVlistID);
      int missTimes = 0;
      for (missTimes = 0; missTimes < numDates - 1; missTimes++)
      {
        auto numFields = streamInqTimestep(cdiStream, missTimes);
        if (numFields == 0) break;

        cdiDateTimes[numDates - 1] = taxisInqVdatetime(cdiTaxisID);
        cdiDateTimes[numDates - 1].date.year = endYear + 1;

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          int varID, levelID;
          streamInqField(cdiStream, &varID, &levelID);
          auto &rvars1 = varsData1[numDates - 1][varID][levelID];
          streamReadField(cdiStream, rvars1.vec_d.data(), &rvars1.numMissVals);

          if (lvarstd)
          {
            field2_moq(varsData2[numDates - 1][varID][levelID], rvars1);
            for (int inp = 0; inp < numDates - 1; ++inp)
              field2_sumsumq(varsData1[inp][varID][levelID], varsData2[inp][varID][levelID], rvars1);
          }
          else
          {
            for (int inp = 0; inp < numDates - 1; ++inp) field2_function(varsData1[inp][varID][levelID], rvars1, operfunc);
          }
        }

        cdiDateTimes[numDates] = datetime_avg(dpy, numDates, cdiDateTimes);
        auto vDateTime = cdiDateTimes[numDates];
        if (vDateTime.date.year > endYear) vDateTime.date.year = endYear;

        ydstat_update(stats, vDateTime, varsData1[0], varsData2[0], numDates, operfunc);

        cdiDateTimes[numDates] = cdiDateTimes[0];
        varsData1[numDates] = varsData1[0];
        if (lvarstd) varsData2[numDates] = varsData2[0];

        for (int inp = 0; inp < numDates; ++inp)
        {
          cdiDateTimes[inp] = cdiDateTimes[inp + 1];
          varsData1[inp] = varsData1[inp + 1];
          if (lvarstd) varsData2[inp] = varsData2[inp + 1];
        }
      }

      if (missTimes != numDates - 1) cdo_abort("Addding the missing values when using the 'readMethod' method was not possible");

      streamClose(cdiStream);
    }

    ydstat_finalize(stats, operfunc);

    int otsID = 0;

    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (stats.numSets[dayOfYear])
      {
        taxisDefVdatetime(taxisID2, stats.vDateTime[dayOfYear]);
        cdo_def_timestep(streamID2, otsID);

        for (int fieldID = 0; fieldID < maxFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();
          if (otsID && varList1.vars[varID].isConstant) continue;

          auto &rvars1 = stats.varsData1[dayOfYear][varID][levelID];

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, rvars1);
        }

        otsID++;
      }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
  }
};
