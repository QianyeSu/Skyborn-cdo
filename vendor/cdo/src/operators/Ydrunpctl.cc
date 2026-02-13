/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida
          Fabian Wachsmann

*/

/*
   This module contains the following operators:

      Ydrunpctl    ydrunpctl         Multi-year daily running percentiles
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "util_string.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "pmlist.h"
#include "field_functions.h"

namespace
{
struct Parameter
{
  double pn{ -1.0 };  // Percentile number
  int nts{ -1 };      // Number of timesteps
  char rm{ 0 };       // Read method (circular)
  std::string pm;     // Percentile method
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs < 2) cdo_abort("Too few arguments!");

  auto argsList = cdo_get_oper_argv();
  auto param1 = argsList[0];
  auto param2 = argsList[1];
  if (std::isdigit(param1[0]) && !string_contains(param1, '='))
  {
    params.pn = parameter_to_double(param1);
    params.nts = parameter_to_int(param2);
    argsList.erase(argsList.begin(), argsList.begin() + 2);
    numArgs -= 2;
  }

  if (numArgs)
  {
    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argsList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "nts") params.nts = parameter_to_int(value);
      else if (key == "p")   params.pn = parameter_to_double(value);
      else if (key == "rm")  params.rm = value[0];
      else if (key == "pm")  { auto pm = parameter_to_word(value); params.pm = (pm == "r8") ? "rtype8" : pm; }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
check_parameter(const Parameter &parameter)
{
  if (is_equal(parameter.pn, -1)) cdo_abort("Too few parameter!");
  if (parameter.pn < 0 || parameter.pn > 100) cdo_abort("Parameter p=%g out of range (0-100)!", parameter.pn);
  if (parameter.nts == -1) cdo_abort("Too few parameter!");
  if (parameter.nts <= 0) cdo_abort("Parameter nts must be greater than 0!");
  if (parameter.rm != 0 && parameter.rm != 'c') cdo_abort("Parameter rm must only contain 'c'!");
}

constexpr int MaxDays = 373;
class Ydrunpctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ydrunpctl",
    .operators = { { "ydrunpctl", FieldFunc_Pctl, 0, YdrunpctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Ydrunpctl> registration = RegisterEntry<Ydrunpctl>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};
  int taxisID4{};

  int numDates{ -1 };
  int dpy{};
  double pn{ 0.0 };

  char readMethod{ 0 };

  VarList varList1{};

public:
  void
  init() override
  {
    auto params = get_parameter();
    check_parameter(params);
    pn = params.pn;
    numDates = params.nts;
    readMethod = params.rm;
    if (params.pm.size()) percentile_set_method(params.pm);
    if (Options::cdoVerbose) cdo_print("pn=%g numDates=%d readMethod=%c pm=%s", pn, numDates, readMethod, params.pm);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = cdo_stream_inq_vlist(streamID3);
    auto vlistID4 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID4);

    varList1 = VarList(vlistID1);
    VarList varList2(vlistID2);
    VarList varList3(vlistID3);

    varList_compare(varList1, varList2);
    varList_compare(varList1, varList3);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = vlistInqTaxis(vlistID3);
    // TODO - check that time axes 2 and 3 are equal

    taxisID4 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
    vlistDefTaxis(vlistID4, taxisID4);

    dpy = calendar_dpy(taxisInqCalendar(taxisID1));

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);
  }

  void
  run() override
  {
    Field field1, field2;
    FieldVector3D varsData1(numDates + 1);
    for (int its = 0; its < numDates; its++) field2D_init(varsData1[its], varList1, FIELD_VEC | FIELD_NAT);

    std::vector<bool> vars2(MaxDays, false);
    CdiDateTime vDateTimes1[MaxDays]{};
    CdiDateTime vDateTimes2[MaxDays]{};
    HistogramSet hsets[MaxDays];
    int numSets[MaxDays] = { 0 };

    auto numVars = varList1.numVars();
    auto numSteps = varList1.numSteps();
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    FieldVector constFields(maxFields);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      if (numFields != cdo_stream_inq_timestep(streamID3, tsID))
        cdo_abort("Number of fields at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto vDateTime = taxisInqVdatetime(taxisID2);

      if (decode_month_and_day(vDateTime.date) != decode_month_and_day(taxisInqVdatetime(taxisID3).date))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      // if (Options::cdoVerbose) cdo_print("process timestep: %d %d %d", tsID + 1, vdate, vtime);

      auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes2[dayOfYear] = vDateTime;

      if (!vars2[dayOfYear])
      {
        vars2[dayOfYear] = true;
        hsets[dayOfYear].create(numVars, numSteps);

        for (auto const &var : varList1.vars) hsets[dayOfYear].createVarLevels(var.ID, var.nlevels, var.gridsize);
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        auto const &var = varList1.vars[varID];
        field1.init(var);
        cdo_read_field(streamID2, field1);

        (void) cdo_inq_field(streamID3);
        field2.init(var);
        cdo_read_field(streamID3, field2);

        hsets[dayOfYear].defVarLevelBounds(varID, levelID, field1, field2);
      }

      tsID++;
    }

    std::vector<CdiDateTime> cdiDateTimes(numDates + 1);
    for (tsID = 0; tsID < numDates; ++tsID)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) cdo_abort("File has less then %d timesteps!", numDates);

      cdiDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];

        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

        if (tsID == 0 && var.isConstant)
        {
          constFields[fieldID].init(var);
          cdo_read_field(streamID1, constFields[fieldID]);
        }
        else { cdo_read_field(streamID1, varsData1[tsID][varID][levelID]); }
      }
    }

    while (true)
    {
      cdiDateTimes[numDates] = datetime_avg(dpy, numDates, cdiDateTimes);

      auto vDateTime = cdiDateTimes[numDates];

      auto dayOfYear = decode_day_of_year(vDateTime.date);
      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      vDateTimes1[dayOfYear] = vDateTime;

      if (!vars2[dayOfYear])
        cdo_abort("No data for day %d in %s and %s", dayOfYear, cdo_get_stream_name(1), cdo_get_stream_name(2));

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (var.isConstant) continue;

        for (int levelID = 0; levelID < var.nlevels; ++levelID)
          for (int inp = 0; inp < numDates; ++inp)
            hsets[dayOfYear].addVarLevelValues(varID, levelID, varsData1[inp][varID][levelID]);
      }

      cdiDateTimes[numDates] = cdiDateTimes[0];
      varsData1[numDates] = varsData1[0];

      for (int inp = 0; inp < numDates; ++inp)
      {
        cdiDateTimes[inp] = cdiDateTimes[inp + 1];
        varsData1[inp] = varsData1[inp + 1];
      }

      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdiDateTimes[numDates - 1] = taxisInqVdatetime(taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, varsData1[numDates - 1][varID][levelID]);
      }

      numSets[dayOfYear] += numDates;
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
          auto &pvars1 = varsData1[numDates - 1][varID][levelID];
          if (pvars1.memType == MemType::Float)
            streamReadFieldF(cdiStream, pvars1.vec_f.data(), &pvars1.numMissVals);
          else
            streamReadField(cdiStream, pvars1.vec_d.data(), &pvars1.numMissVals);
        }

        cdiDateTimes[numDates] = datetime_avg(dpy, numDates, cdiDateTimes);
        auto vDateTime = cdiDateTimes[numDates];
        if (vDateTime.date.year > endYear) vDateTime.date.year = endYear;

        auto dayOfYear = decode_day_of_year(vDateTime.date);
        if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);
        vDateTimes1[dayOfYear] = vDateTime;

        numSets[dayOfYear] += numDates;

        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList1.vars[varID];
          if (var.isConstant) continue;

          for (int levelID = 0; levelID < var.nlevels; ++levelID)
            for (int inp = 0; inp < numDates; ++inp)
              hsets[dayOfYear].addVarLevelValues(varID, levelID, varsData1[inp][varID][levelID]);
        }

        cdiDateTimes[numDates] = cdiDateTimes[0];
        varsData1[numDates] = varsData1[0];

        for (int inp = 0; inp < numDates; ++inp)
        {
          cdiDateTimes[inp] = cdiDateTimes[inp + 1];
          varsData1[inp] = varsData1[inp + 1];
        }
      }

      if (missTimes != numDates - 1) cdo_abort("Addding the missing values when using the 'readMethod' method was not possible");

      streamClose(cdiStream);
    }

    /*
    int outyear = 1e9;
    for (dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
        {
          int year, month, day;
          cdiDate_decode(vDateTimes1[dayOfYear].date, &year, &month, &day);
          if (year < outyear) outyear = year;
        }

    for (dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
        {
          int year, month, day;
          cdiDate_decode(vDateTimes1[dayOfYear].date, &year, &month, &day);
          vDateTimes1[dayOfYear].date = cdiDate_encode(outyear, month, day);
        }
    */
    int otsID = 0;
    for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
      if (numSets[dayOfYear])
      {
        if (decode_month_and_day(vDateTimes1[dayOfYear].date) != decode_month_and_day(vDateTimes2[dayOfYear].date))
          cdo_abort("Verification dates for day %d of %s and %s differ!", dayOfYear, cdo_get_stream_name(0),
                    cdo_get_stream_name(1));

        taxisDefVdatetime(taxisID4, vDateTimes1[dayOfYear]);
        cdo_def_timestep(streamID4, otsID);

        for (int fieldID = 0; fieldID < maxFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();
          auto const &var = varList1.vars[varID];
          if (otsID && var.isConstant) continue;

          cdo_def_field(streamID4, varID, levelID);

          if (var.isConstant) { cdo_write_field(streamID4, constFields[fieldID]); }
          else
          {
            field1.init(var);
            hsets[dayOfYear].getVarLevelPercentiles(field1, varID, levelID, pn);
            cdo_write_field(streamID4, field1);
          }
        }

        otsID++;
      }
  }

  void
  close() override
  {
    cdo_stream_close(streamID4);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
  }
};
