/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Seltime    seltimestep     Select timesteps
      Seltime    seltime         Select times
      Seltime    selhour         Select hours
      Seltime    selday          Select days
      Seltime    selmonth        Select months
      Seltime    selyear         Select years
      Seltime    selseason       Select seasons
      Seltime    seldate         Select dates
      Seltime    selsmon         Select single month
*/
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_string.h"
#include "param_conversion.h"
#include "field_functions.h"

static std::vector<int>
get_season_list(std::vector<std::string> const &seasonString)
{
  std::vector<int> listArrayInt;
  int imon[13]{ 0 };  // 1-12 !

  int numSel = seasonString.size();
  if (std::isdigit(seasonString[0][0]))
  {
    for (int i = 0; i < numSel; ++i)
    {
      auto ival = parameter_to_int(seasonString[i]);
      // clang-format off
      if      (ival == 1) { imon[12]++; imon[ 1]++; imon[ 2]++; }
      else if (ival == 2) { imon[ 3]++; imon[ 4]++; imon[ 5]++; }
      else if (ival == 3) { imon[ 6]++; imon[ 7]++; imon[ 8]++; }
      else if (ival == 4) { imon[ 9]++; imon[10]++; imon[11]++; }
      else cdo_abort("Season %d not available!", ival);
      // clang-format on
    }
  }
  else
  {
    for (int i = 0; i < numSel; ++i) season_to_months(seasonString[i], imon);
  }

  for (int i = 1; i < 13; ++i)
    if (imon[i]) listArrayInt.push_back(i);

  return listArrayInt;
}

// last update: 2020-02-17 (Oliver Heidmann)
/*
 * Input:
 *      std::vector<std::string> => vector of length 1 or 2 containing string representing a range of dates or a single date to be
 * selected Output: std::vector<double>       => vector of length 1 or 2 containing the to double converted date or date range
 *
 * This function turns a date range string into a list of double values which represent the start and the end of the range.
 * When only one argument is given and it contains no time string (marked with T e.g 2006-01-01T23:00:30) we set the second value
 * of the retrun vector to the end of the day from the first argument.
 * If only one argument is given AND it contains a time string we return a vector of length one which contains only
 * the converted value of the first argument.
 * This function accepts the string "-".
 * This represents the possible maximum start or end date depending on where the string is in the arguments.
 * If the first argument is "-" we set the start date to -999999999 and to 999999999 if the string is in the second.
 */

static std::vector<double>
get_date_list(std::vector<std::string> const &argVec)
{
  if (argVec.size() < 1) cdo_abort("Too few arguments!");
  if (argVec.size() > 2) cdo_abort("Too many arguments!");

  std::vector<double> dateList;
  short containsTime = -1;
  for (size_t i = 0, n = argVec.size(); i < n; ++i)
  {
    double fval = 0.0;
    // if "-" set start date to maximum start (i == 0) or end (i == 1)
    if (argVec[i] == "-") { fval = (i == 0) ? -99999999999. : 99999999999.; }
    // get the double value representing the date and check for time string
    else
    {
      fval = datestr_to_double(argVec[i], i);
      if (i == 0 && string_contains(argVec[i], '-')) containsTime = string_contains(argVec[i], 'T');
    }
    dateList.push_back(fval);
  }
  // if date and only one argument: set second date to first and set time to end of day
  if (argVec.size() == 1) dateList.push_back(dateList[0] + ((containsTime != -1 && containsTime == 0) ? 0.235959 : 0.0));

  return dateList;
}

static std::string
get_argv_string(std::vector<std::string> const &argVec)
{
  std::string s = argVec[0];
  for (size_t i = 1, n = argVec.size(); i < n; ++i) s += "," + argVec[i];
  return s;
}

static bool
argv_has_negativ_values(std::vector<std::string> const &argVec)
{
  for (auto const &argument : argVec)
  {
    int first, last, inc;
    split_intstring(argument, first, last, inc);
    if (first < 0 || last < 0) return true;
  }

  return false;
}

static std::vector<int>
cdo_argv_to_int_timestep(std::vector<std::string> const &argVec, int ntimesteps)
{
  std::vector<int> v;

  for (auto const &argument : argVec)
  {
    int first, last, inc;
    split_intstring(argument, first, last, inc);

    if (inc == -1 && (first < 0 || last < 0)) inc = 1;
    if (first < 0) first += ntimesteps + 1;
    if (last < 0) last += ntimesteps + 1;

    if (inc >= 0)
    {
      for (auto ival = first; ival <= last; ival += inc) v.push_back(ival);
    }
    else
    {
      for (auto ival = first; ival >= last; ival += inc) v.push_back(ival);
    }
  }

  return v;
}

class Seltime : public Process
{
  enum
  {
    func_time,
    func_date,
    func_step,
    func_datetime
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Seltime",
    // clang-format off
    .operators = { { "seltimestep", func_step, 0, "timesteps", SeltimeHelp },
                   { "selyear",     func_date, 0, "years", SeltimeHelp },
                   { "selseason",   func_date, 0, "seasons", SeltimeHelp },
                   { "selmonth",    func_date, 0, "months", SeltimeHelp },
                   { "selday",      func_date, 0, "days", SeltimeHelp },
                   { "selhour",     func_time, 0, "hours", SeltimeHelp },
                   { "seldate",     func_datetime, 0, "startdate and enddate (format: YYYY-MM-DDThh:mm:ss)", SeltimeHelp },
                   { "seltime",     func_time, 0, "times (format: hh:mm:ss)", SeltimeHelp },
                   { "selsmon",     func_date, 0, "month[,nts1[,nts2]]", SeltimeHelp } },
    // clang-format on
    .aliases = { { "selseas", "selseason" }, { "selmon", "selmonth" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Seltime> registration = RegisterEntry<Seltime>();

  int SELTIMESTEP, SELDATE, SELTIME, SELHOUR, SELDAY, SELMONTH, SELYEAR, SELSEASON, SELSMON;
  CdoStreamID streamID1;
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int operfunc;
  int operatorID;

  int numSel = 0;
  int numVars;
  int ncts = 0, nts;
  int nts1 = 0, nts2 = 0;
  int its1 = 0, its2 = 0;
  int iselmax = -1;
  int tsmax = -1;

  double selfval = 0;

  bool isConstOut = false;
  bool process_nts1 = false, process_nts2 = false;
  bool dataIsUnchanged = false;
  bool lnts1;
  bool haveConstVars;

  std::vector<int> intarr;
  std::vector<double> fltarr;

  VarList varList1;
  FieldVector3D varsData;
  std::vector<CdiDateTime> vDateTimes;

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    SELTIMESTEP = module.get_id("seltimestep");
    SELDATE = module.get_id("seldate");
    SELTIME = module.get_id("seltime");
    SELHOUR = module.get_id("selhour");
    SELDAY = module.get_id("selday");
    SELMONTH = module.get_id("selmonth");
    SELYEAR = module.get_id("selyear");
    SELSEASON = module.get_id("selseason");
    SELSMON = module.get_id("selsmon");

    operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_input_arg(cdo_operator_enter(operatorID));
    auto argv = get_argv_string(cdo_get_oper_argv());
    std::string newCommand{};

    auto redirectEnabled = false;  // only for seperating prototype from code (redirects execution to Select.cc)
    auto redirectFound = true;     // for opers that are not redirectable for now
    auto ID = operatorID;
    if (redirectEnabled)
    {
      // clang-format off
      if      (ID == SELTIMESTEP ) { newCommand += "timestep=" + argv; }
      else if (ID == SELDATE)      { newCommand += "date="     + argv; }
      //else if(ID == SELTIME)    { newCommand += "time="      + argv; }  // unimplemented
      else if (ID == SELHOUR)      { newCommand += "hour="     + argv; }
      else if (ID == SELDAY)       { newCommand += "day="      + argv; }
      else if (ID == SELMONTH)     { newCommand += "month="    + argv; }
      else if (ID == SELYEAR)      { newCommand += "year="     + argv; }
      else if (ID == SELSEASON)    { newCommand += "season="   + argv; }
      else if (ID == SELSMON)      { newCommand += "month="    + argv; }
      else                         { redirectFound = false; }
      // clang-format on
    }

    if (redirectFound && redirectEnabled)
    {
      cdo_abort("Redirecting was disabled and does no longer work");
      /*
       * This was a temporary feature and does no longer work with the new module system.
       * Redirecting has to be implemented in antother way or Select needs to be able to handle this kind of
       * input.
       *
       * Debug(Yellow("Redirecting to %s"), newCommand);
       * ((Process *) process)->init_process("select", { newCommand });
       * Select(process);
       * return;
       * // If a redirect was found the entire process is ended through this return! */
    }

    if (operatorID == SELSEASON)
    {
      intarr = get_season_list(cdo_get_oper_argv());
      numSel = intarr.size();
    }
    else if (operatorID == SELDATE)
    {
      fltarr = get_date_list(cdo_get_oper_argv());
      numSel = fltarr.size();
    }
    else if (operatorID == SELTIME)
    {
      numSel = cdo_operator_argc();
      if (numSel < 1) cdo_abort("Too few arguments!");
      intarr.reserve(numSel);
      for (int i = 0; i < numSel; ++i)
      {
        auto currentArgument = cdo_operator_argv(i).c_str();
        if (std::strchr(currentArgument, ':'))
        {
          int hour = 0, minute = 0, second = 0;
          std::sscanf(currentArgument, "%d:%d:%d", &hour, &minute, &second);
          intarr.push_back(cdiEncodeTime(hour, minute, second));
        }
        else { intarr.push_back(parameter_to_int(currentArgument)); }
      }
    }
    else if (operatorID != SELTIMESTEP)
    {
      intarr = cdo_argv_to_intarr(cdo_get_oper_argv());
      numSel = intarr.size();
    }

    if (operatorID != SELTIMESTEP && numSel < 1) cdo_abort("No timestep selected!");

    if (operatorID == SELSMON)
    {
      if (numSel > 1) nts1 = intarr[1];
      nts2 = (numSel > 2) ? intarr[2] : nts1;

      if (numSel > 3) cdo_abort("Too many parameters");

      if (Options::cdoVerbose) cdo_print("mon=%d  nts1=%d  nts2=%d", intarr[0], nts1, nts2);

      numSel = 1;
    }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    // cdo::delete_chunks_dimT(vlistID2);

    varList1 = VarList(vlistID1);

    // add support for negative timestep values
    if (operatorID == SELTIMESTEP)
    {
      if (argv_has_negativ_values(cdo_get_oper_argv()))
      {
        auto numSteps = varList1.numSteps();
        if (numSteps < 0)
        {
          if (cdo_assert_files_only())
          {
            int tsID = 0;
            while (cdo_stream_inq_timestep(streamID1, tsID)) tsID++;
            numSteps = tsID;
            if (Options::cdoVerbose) cdo_print("Found %d timesteps", numSteps);
          }

          cdo_abort("Negative timesteps not supported in CDO pipes!");
        }

        intarr = cdo_argv_to_int_timestep(cdo_get_oper_argv(), numSteps);
      }
      else { intarr = cdo_argv_to_intarr(cdo_get_oper_argv()); }

      numSel = intarr.size();

      if (numSel < 1) cdo_abort("No timestep selected!");
    }

    if (operatorID == SELTIMESTEP)
    {
      std::ranges::sort(intarr);
      auto ip = std::unique(intarr.begin(), intarr.end());
      intarr.resize(std::distance(intarr.begin(), ip));
      numSel = intarr.size();
    }

    if (Options::cdoVerbose)
    {
      for (int i = 0; i < numSel; ++i)
      {
        if (operatorID == SELDATE)
          cdo_print("fltarr entry: %d %14.4f", i + 1, fltarr[i]);
        else
          cdo_print("intarr entry: %d %d", i + 1, intarr[i]);
      }
    }

    auto numStepsOut = (operfunc == func_step) ? numSel : -1;
    vlistDefNtsteps(vlistID2, numStepsOut);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if (operatorID != SELDATE)
      for (int i = 0; i < numSel; ++i) iselmax = std::max(iselmax, intarr[i]);

    numVars = varList1.numVars();
    haveConstVars = (varList1.numConstVars() > 0);

    lnts1 = (operatorID == SELSMON) && (nts1 > 0);

    if (lnts1 || haveConstVars)
    {
      if (lnts1) { vDateTimes.resize(nts1); }
      else { nts1 = 1; }

      varsData.resize(nts1);

      for (int tsID = 0; tsID < nts1; ++tsID)
      {
        field2D_init(varsData[tsID], varList1);

        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList1.vars[varID];
          if (lnts1 || var.isConstant)
          {
            for (int levelID = 0; levelID < var.nlevels; ++levelID) varsData[tsID][varID][levelID].resize(var.gridsize);
          }
        }
      }
    }

    if (operatorID == SELTIMESTEP)
      for (int i = 0; i < numSel; ++i) tsmax = std::max(tsmax, intarr[i]);
  }

  void
  run() override
  {
    std::vector<bool> selfound(numSel, false);
    Field field;

    int indexNext = 0;
    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto vdate = cdiDate_get(vDateTime.date);
      auto vtime = cdiTime_get(vDateTime.time);

      auto copytimestep = false;
      int selival = -1;

      if (operfunc == func_step)
      {
        selival = tsID + 1;
        if (selival > iselmax) break;
      }
      else if (operfunc == func_date)
      {
        int year, month, day;
        cdiDate_decode(vDateTime.date, &year, &month, &day);
        selival = (operatorID == SELYEAR) ? year : (operatorID == SELDAY) ? day : month;
      }
      else if (operfunc == func_time)
      {
        int hour, minute, second, ms;
        cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
        selival = (operatorID == SELHOUR) ? hour : vtime;
      }
      else if (operfunc == func_datetime) { selfval = vdate + vtime / 1000000.0; }

      if (operatorID == SELDATE)
      {
        if (selfval >= fltarr[0] && selfval <= fltarr[numSel - 1])
        {
          copytimestep = true;
          selfound[0] = true;
          selfound[numSel - 1] = true;
        }
        else if (selfval > fltarr[numSel - 1]) { break; }
      }
      else if (operatorID == SELTIMESTEP)
      {
        if (tsID >= tsmax) break;

        int index;
        for (index = indexNext; index < numSel; ++index)
          if (selival == intarr[index]) break;
        if (index < numSel)
        {
          copytimestep = true;
          selfound[index] = true;
          indexNext = index + 1;
        }
      }
      else
      {
        for (int i = 0; i < numSel; ++i)
          if (selival == intarr[i])
          {
            copytimestep = true;
            selfound[i] = true;
            break;
          }
      }

      auto copy_nts2 = false;
      if (operatorID == SELSMON && !copytimestep)
      {
        copy_nts2 = false;

        if (process_nts1)
        {
          process_nts2 = true;
          its2 = 0;
          process_nts1 = false;
        }

        if (process_nts2)
        {
          if (its2++ < nts2)
            copy_nts2 = true;
          else
            process_nts2 = false;
        }
      }

      if (copytimestep || copy_nts2)
      {
        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        if (tsID2 == 0)
        {
          streamID2 = cdo_open_write(1);
          cdo_def_vlist(streamID2, vlistID2);
        }

        if (lnts1 && ncts == 0)
        {
          nts = nts1;
          if (its1 < nts1)
          {
            nts = its1;
            cdo_warning("%d timesteps missing before month %d!", nts1 - its1, intarr[0]);
          }

          for (int it = 0; it < nts; it++)
          {
            taxisDefVdatetime(taxisID2, vDateTimes[it]);
            cdo_def_timestep(streamID2, tsID2++);

            for (int varID = 0; varID < numVars; ++varID)
            {
              auto const &var = varList1.vars[varID];
              if (var.isConstant && tsID2 > 1) continue;
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
              {
                cdo_def_field(streamID2, varID, levelID);
                auto single = varsData[it][varID][levelID].vec_d.data();
                auto numMissVals = varsData[it][varID][levelID].numMissVals;
                cdo_write_field(streamID2, single, numMissVals);
              }
            }
          }

          its1 = 0;
        }

        ncts++;
        if (!process_nts2)
        {
          its2 = 0;
          process_nts1 = true;
        }

        cdo_def_timestep(streamID2, tsID2++);

        if (tsID > 0 && isConstOut)
        {
          isConstOut = false;
          nts = nts1 - 1;
          for (int varID = 0; varID < numVars; ++varID)
          {
            auto const &var = varList1.vars[varID];
            if (var.isConstant)
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
              {
                cdo_def_field(streamID2, varID, levelID);
                auto single = varsData[nts][varID][levelID].vec_d.data();
                auto numMissVals = varsData[nts][varID][levelID].numMissVals;
                cdo_write_field(streamID2, single, numMissVals);
              }
            }
          }
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          cdo_def_field(streamID2, varID, levelID);
          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            auto const &var = varList1.vars[varID];
            field.init(var);
            cdo_read_field(streamID1, field);
            cdo_write_field(streamID2, field);
          }
        }
      }
      else
      {
        ncts = 0;

        if (lnts1 || tsID == 0)
        {
          if (tsID == 0 && haveConstVars && (!lnts1)) isConstOut = true;

          nts = nts1 - 1;
          if (lnts1)
          {
            if (its1 <= nts)
              nts = its1;
            else
              for (int it = 0; it < nts; it++)
              {
                vDateTimes[it] = vDateTimes[it + 1];
                for (int varID = 0; varID < numVars; ++varID)
                {
                  auto const &var = varList1.vars[varID];
                  if (var.isConstant) continue;
                  for (int levelID = 0; levelID < var.nlevels; ++levelID)
                  {
                    varsData[it][varID][levelID].vec_d = varsData[it + 1][varID][levelID].vec_d;
                    varsData[it][varID][levelID].numMissVals = varsData[it + 1][varID][levelID].numMissVals;
                  }
                }
              }

            vDateTimes[nts] = taxisInqVdatetime(taxisID1);

            its1++;
          }

          for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = cdo_inq_field(streamID1);
            auto const &var = varList1.vars[varID];
            if (lnts1 || var.isConstant)
            {
              auto single = varsData[nts][varID][levelID].vec_d.data();
              cdo_read_field(streamID1, single, &varsData[nts][varID][levelID].numMissVals);
            }
          }
        }
      }

      tsID++;
    }

    if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (operatorID == SELSMON)
      if (its2 < nts2) cdo_warning("%d timesteps missing after the last month!", nts2 - its2);

    for (int isel = 0; isel < numSel; isel++)
    {
      if (selfound[isel] == false)
      {
        if (operatorID == SELTIMESTEP)
        {
          int isel2;
          auto lcont = false;
          for (isel2 = isel + 1; isel2 < numSel; isel2++)
            if (selfound[isel2]) break;
          if (isel2 == numSel && (numSel - isel) > 1) lcont = true;

          auto lcont2 = false;
          if (lcont)
          {
            for (isel2 = isel + 1; isel2 < numSel; isel2++)
              if (intarr[isel2 - 1] != intarr[isel2] - 1) break;
            if (isel2 == numSel) lcont2 = true;
          }

          if (lcont2)
          {
            cdo_abort("Timesteps %d-%d not found!", intarr[isel], intarr[numSel - 1]);
            break;
          }
          else
            cdo_abort("Timestep %d not found!", intarr[isel]);
        }
        else if (operatorID == SELDATE)
        {
          if (isel == 0) cdo_warning("Date between %14.4f and %14.4f not found!", fltarr[0], fltarr[numSel - 1]);
        }
        else if (operatorID == SELTIME) { cdo_warning("Time %d not found!", intarr[isel]); }
        else if (operatorID == SELHOUR) { cdo_warning("Hour %d not found!", intarr[isel]); }
        else if (operatorID == SELDAY) { cdo_warning("Day %d not found!", intarr[isel]); }
        else if (operatorID == SELMONTH) { cdo_warning("Month %d not found!", intarr[isel]); }
        else if (operatorID == SELYEAR) { cdo_warning("Year %d not found!", intarr[isel]); }
        else if (operatorID == SELSEASON)
        {
          if (isel < 3) cdo_warning("Month %d not found!", intarr[isel]);
        }
      }
    }
    if (tsID2 == 0) cdo_abort("No timesteps selected!");
  }

  void
  close() override
  {
    vlistDestroy(vlistID2);
  }
};
