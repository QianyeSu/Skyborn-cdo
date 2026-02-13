/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Select      select         Select fields
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "datetime.h"
#include "sellist.h"
#include "printinfo.h"
#include "param_conversion.h"
#include "progress.h"
#include "cdi_lockedIO.h"
#include "cdo_cdi_wrapper.h"

std::vector<bool> cdo_read_timestepmask(std::string const &maskfile);

std::string
dom_to_string(CdiDate date)
{
  int year, month, day;
  cdiDate_decode(date, &year, &month, &day);

  constexpr char const *cmons[] = { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
  if (month < 0 || month > 12) month = 0;
  char cstr[32];
  std::snprintf(cstr, sizeof(cstr), "%d%s", day, cmons[month]);

  return std::string(cstr);
}

static void
write_const_vars(CdoStreamID streamID2, VarList const &varList, Varray2D<double> &varsData2)
{
  for (auto const &var : varList.vars)
  {
    if (varsData2[var.ID].size())
    {
      for (int levelID2 = 0; levelID2 < var.nlevels; ++levelID2)
      {
        auto pdata = &varsData2[var.ID][var.gridsize * levelID2];
        auto numMissVals = array_num_mv(var.gridsize, pdata, var.missval);
        cdo_def_field(streamID2, var.ID, levelID2);
        cdo_write_field(streamID2, pdata, numMissVals);
      }
      varsData2[var.ID].clear();
      varsData2[var.ID].shrink_to_fit();
    }
  }
}

static void
eval_timestepmask(std::string const &maskfile, KVList &kvlist)
{
  auto imask = cdo_read_timestepmask(maskfile);
  int n = imask.size();
  if (n == 0) cdo_abort("Read of timestep mask failed!");

  int nvalues = 0;
  for (int i = 0; i < n; ++i)
    if (imask[i]) nvalues++;
  if (nvalues == 0)
  {
    cdo_print("timestepmask has no values!");
    return;
  }

  KeyValues kv;
  kv.key = "timestep";
  kv.nvalues = nvalues;
  kv.values.resize(nvalues);

  std::vector<char> value(32);
  for (int i = 0, j = 0; i < n; ++i)
    if (imask[i])
    {
      std::snprintf(value.data(), value.size(), "%d", i + 1);
      kv.values[j++] = value.data();
    }

  kvlist.push_back(kv);
}

static bool
has_selected_params(VarList const &varList, int vlistID)
{
  for (auto const &var : varList.vars)
  {
    for (int levelID = 0; levelID < var.nlevels; ++levelID)
      if (vlistInqFlag(vlistID, var.ID, levelID) == true) return true;
  }

  return false;
}

static bool
has_const_vars(VarList const &varList, std::vector<bool> const &processVars)
{
  for (auto const &var : varList.vars)
  {
    if (processVars[var.ID] && var.isConstant) return true;
  }

  return false;
}

static double
datetime_to_double(CdiDateTime dateTime)
{
  auto vdate = cdiDate_get(dateTime.date);
  auto vtime = cdiTime_get(dateTime.time);
  auto fdatetime = ((double) vdate) + ((double) vtime) / 1000000.0;
  return fdatetime;
}

class Select : public Process
{
public:
  using Process::Process;
#ifndef UWES_SELECT_TEST
  inline static CdoModule module = {
    .name = "Select",
    .operators = { { "select", 0, 0, "parameter list", SelectHelp },
                   { "delete",
                     0,
                     0,
                     "parameter list",
                     SelectHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
#else
  inline static CdoModule module = {
    .name = "Select1",
    .operators = { { "select1", 0, 0, "parameter list", SelectHelp }, { "delete1", 0, 0, "parameter list", SelectHelp } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
#endif
  inline static RegisterEntry<Select> registration = RegisterEntry<Select>();

private:
  int SELECT{}, DELETE{};
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;
  int tsID2 = 0;
  int vlistID0 = -1, vlistID2 = -1;
  bool dataIsUnchanged{};

  KVList kvlist{};
  int operatorID{};

public:
  void
  init() override
  {
    SELECT = module.get_id("select");
    DELETE = module.get_id("delete");

    dataIsUnchanged = data_is_unchanged();

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto numArgs = cdo_operator_argc();
    auto const &argList = cdo_get_oper_argv();

    if (numArgs == 0) cdo_abort("Parameter missing!");

    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    auto kv = kvlist.search("timestepmask");
    if (kv && kv->nvalues > 0)
    {
      if (kvlist.search("timestep")) cdo_abort("Parameter timestep and timestepmask can't be combined!");
      eval_timestepmask(kv->values[0], kvlist);
    }
  }

  void
  run() override
  {
    int numVars2 = 0;
    bool hasConstVars = true;
    char paramstr[32];
    char gname[CDI_MAX_NAME];
    char zname[CDI_MAX_NAME];
    int last_year = -999999999;
    Field field;
    double fstartdate = -99999999999.0;
    double fenddate = -99999999999.0;
    int taxisID2 = CDI_UNDEFID;
    int numStepsOut = 0;
    bool doTimeSel = false;
    std::vector<bool> processVars;
    Varray2D<double> varsData2;

    SelectInfo selInfo(kvlist);

    // clang-format off
    SELINFO_ADD_INT(timestep_of_year, "Timestep of year");
    SELINFO_ADD_INT(timestep,         "Timestep");
    SELINFO_ADD_INT(year,             "Year");
    SELINFO_ADD_INT(month,            "Month");
    SELINFO_ADD_INT(day,              "Day");
    SELINFO_ADD_INT(hour,             "Hour");
    SELINFO_ADD_INT(minute,           "Minute");
    SELINFO_ADD_INT(code,             "Code number");
    SELINFO_ADD_INT(levidx,           "Level index");
    SELINFO_ADD_INT(ltype,            "Level type");
    SELINFO_ADD_INT(zaxisnum,         "Zaxis number");
    SELINFO_ADD_INT(gridnum,          "Grid number");
    SELINFO_ADD_FLT(level,            "Level");
    SELINFO_ADD_FLT(levrange,         "Level range");
    SELINFO_ADD_WORD(name,            "Variable name");
    SELINFO_ADD_WORD(param,           "Parameter");
    SELINFO_ADD_WORD(zaxisname,       "Zaxis name");
    SELINFO_ADD_WORD(gridname,        "Grid name");
    SELINFO_ADD_WORD(steptype,        "Time step type");
    SELINFO_ADD_WORD(startdate,       "Start date");
    SELINFO_ADD_WORD(enddate,         "End date");
    SELINFO_ADD_WORD(season,          "Season");
    SELINFO_ADD_WORD(date,            "Date");
    SELINFO_ADD_WORD(timestepmask,    "Timestep mask");
    SELINFO_ADD_WORD(dom,             "Day of month");
    // clang-format on

    if (Options::cdoVerbose) selInfo.print();

    selInfo.verify();

    if (SELINFO_NVAL(levrange) > 0 && SELINFO_NVAL(levrange) != 2) cdo_abort("Key levrange needs two values!");
    if (SELINFO_NVAL(timestepmask) > 1) cdo_abort("Key timestepmask has too many values!");
    (void) (levrange);      // unused
    (void) (timestepmask);  // unused

    auto numFiles = cdo_stream_cnt() - 1;

    DateTimeList dtlist;

    cdo::Progress progress(get_id());

    int tsmax = -1;

    timestep = 0;
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      if (Options::cdoVerbose) cdo_print("Process file: %s", cdo_get_stream_name(fileIdx));

      auto streamID1 = cdo_open_read(fileIdx);

      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1(vlistID1);
      auto numSteps = varList1.numSteps();

      VarList varList2;

      auto copyConstVars = false;

      if (fileIdx == 0 || (fileIdx > 0 && tsID2 == 0))
      {
        auto xresult = true;

        // vlistID0 = vlistDuplicate(vlistID1);

        vlistClearFlag(vlistID1);
        int numVars = varList1.numVars();
        processVars.resize(numVars);

        if (operatorID == DELETE)
        {
          xresult = false;
          for (int varID = 0; varID < numVars; ++varID)
          {
            auto const &var = varList1.vars[varID];
            for (int levelID = 0; levelID < var.nlevels; ++levelID) vlistDefFlag(vlistID1, varID, levelID, true);
          }
        }

        auto findVariable = SELINFO_NVAL(code) || SELINFO_NVAL(name) || SELINFO_NVAL(param);

        auto doVarSel = findVariable || SELINFO_NVAL(ltype) || SELINFO_NVAL(zaxisnum) || SELINFO_NVAL(gridnum)
                        || SELINFO_NVAL(zaxisname) || SELINFO_NVAL(gridname) || SELINFO_NVAL(steptype);

        auto doLevSel = SELINFO_NVAL(level) || SELINFO_NVAL(levrange) || SELINFO_NVAL(levidx);

        doTimeSel = SELINFO_NVAL(date) || SELINFO_NVAL(startdate) || SELINFO_NVAL(enddate) || SELINFO_NVAL(season)
                    || SELINFO_NVAL(timestep_of_year) || SELINFO_NVAL(timestep) || SELINFO_NVAL(year) || SELINFO_NVAL(month)
                    || SELINFO_NVAL(day) || SELINFO_NVAL(hour) || SELINFO_NVAL(minute) || SELINFO_NVAL(dom);

        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList1.vars[varID];

          code = var.code;
          name = var.name.c_str();
          // stdname = var.stdname;
          cdiParamToString(var.param, paramstr, sizeof(paramstr));
          param = paramstr;

          auto zaxisID = var.zaxisID;
          ltype = zaxis_to_ltype(zaxisID);

          zaxisnum = vlistZaxisIndex(vlistID1, zaxisID) + 1;
          zaxisName(zaxisInqType(zaxisID), zname);
          zaxisname = zname;

          gridnum = vlistGridIndex(vlistID1, var.gridID) + 1;
          gridName(gridInqType(var.gridID), gname);
          gridname = gname;

          steptype = var.isConstant ? "constant" : "varying";

          auto found_code = SELINFO_CHECK(code);
          auto found_name = SELINFO_CHECK(name);
          auto found_param = SELINFO_CHECK(param);
          auto found_grid = SELINFO_CHECK(gridnum);
          auto found_gname = SELINFO_CHECK(gridname);
          auto found_ltype = SELINFO_CHECK(ltype);
          auto found_zaxis = SELINFO_CHECK(zaxisnum);
          auto found_zname = SELINFO_CHECK(zaxisname);
          auto found_stype = SELINFO_CHECK(steptype);

          if (SELINFO_NVAL(steptype) && !found_stype)
          {
            steptype = cdo::get_steptype_name(var.stepType);
            found_stype = SELINFO_CHECK(steptype);
          }

          bool lstep = SELINFO_NVAL(steptype) ? found_stype : true;
          bool lvar = (found_code || found_name || found_param);
          bool lgrid = (SELINFO_NVAL(gridnum) || SELINFO_NVAL(gridname)) ? (found_grid || found_gname) : true;
          bool lvert = (SELINFO_NVAL(ltype) || SELINFO_NVAL(zaxisnum) || SELINFO_NVAL(zaxisname))
                           ? (found_ltype || found_zaxis || found_zname)
                           : true;

          processVars[varID] = (lvar && lgrid && lvert && lstep);

          if (!processVars[varID] && !lvar && !findVariable)
          {
            if (found_grid || found_gname) { processVars[varID] = true; }
            else if (found_stype) { processVars[varID] = true; }
            else if (found_ltype || found_zaxis || found_zname) { processVars[varID] = true; }
            else if (!doVarSel && (SELINFO_NVAL(levidx) || SELINFO_NVAL(level) || SELINFO_NVAL(levrange)))
            {
              for (int levelID = 0; levelID < var.nlevels; ++levelID)
              {
                levidx = levelID + 1;
                level = cdo_zaxis_inq_level(zaxisID, levelID);
                if (!processVars[varID] && SELINFO_CHECK_INDEX(levidx, var.nlevels)) processVars[varID] = true;
                if (!processVars[varID] && SELINFO_CHECK(level)) processVars[varID] = true;
                if (!processVars[varID] && SELINFO_CHECK_RANGE(levrange, level)) processVars[varID] = true;
              }
            }
          }
        }

        for (int varID = 0; varID < numVars; ++varID)
        {
          if (processVars[varID])
          {
            auto const &var = varList1.vars[varID];
            if (zaxisInqType(var.zaxisID) == ZAXIS_HYBRID)
            {
              auto psvarid = varList_get_psvarid(varList1, var.zaxisID);
              if (psvarid != -1 && !processVars[psvarid]) processVars[psvarid] = true;
            }
          }
        }

        for (int varID = 0; varID < numVars; ++varID)
        {
          if (processVars[varID])
          {
            auto const &var = varList1.vars[varID];
            for (int levelID = 0; levelID < var.nlevels; ++levelID)
            {
              levidx = levelID + 1;
              level = cdo_zaxis_inq_level(var.zaxisID, levelID);

              if (var.nlevels == 1 && is_equal(level, 0.0))
              {
                SELINFO_CHECK(level);
                vlistDefFlag(vlistID1, varID, levelID, xresult);
              }
              else
              {
                if (SELINFO_NVAL(levidx))
                {
                  if (SELINFO_CHECK_INDEX(levidx, var.nlevels)) vlistDefFlag(vlistID1, varID, levelID, xresult);
                }
                else if (SELINFO_NVAL(level))
                {
                  if (SELINFO_CHECK(level)) vlistDefFlag(vlistID1, varID, levelID, xresult);
                }
                else if (SELINFO_NVAL(levrange))
                {
                  if (SELINFO_CHECK_RANGE(levrange, level)) vlistDefFlag(vlistID1, varID, levelID, xresult);
                }
                else { vlistDefFlag(vlistID1, varID, levelID, xresult); }
              }
            }
          }
        }

        SELINFO_CHECK_FLAG(code);
        SELINFO_CHECK_FLAG(levidx);
        SELINFO_CHECK_FLAG(ltype);
        SELINFO_CHECK_FLAG(zaxisnum);
        SELINFO_CHECK_FLAG(gridnum);
        SELINFO_CHECK_FLAG(level);
        SELINFO_CHECK_FLAG(name);
        SELINFO_CHECK_FLAG(param);
        SELINFO_CHECK_FLAG(zaxisname);
        SELINFO_CHECK_FLAG(gridname);
        SELINFO_CHECK_FLAG(steptype);
        SELINFO_CHECK_RANGE_FLAG(levrange);

        if (has_selected_params(varList1, vlistID1))
        {
          if (doVarSel && doTimeSel) copyConstVars = has_const_vars(varList1, processVars);
        }
        else
        {
          if ((!doVarSel) && (!doLevSel) && doTimeSel)
          {
            copyConstVars = true;

            for (int varID = 0; varID < numVars; ++varID)
            {
              processVars[varID] = true;
              auto const &var = varList1.vars[varID];
              for (int levelID = 0; levelID < var.nlevels; ++levelID) vlistDefFlag(vlistID1, varID, levelID, true);
            }
          }
          else { cdo_abort("No variable selected!"); }
        }

        // if (Options::cdoVerbose) vlistPrint(vlistID1);

        vlistID0 = vlistDuplicate(vlistID1);
        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList1.vars[varID];
          for (int levelID = 0; levelID < var.nlevels; ++levelID)
            vlistDefFlag(vlistID0, varID, levelID, vlistInqFlag(vlistID1, varID, levelID));
        }

        // if (Options::cdoVerbose) vlistPrint(vlistID0);

        vlistID2 = vlistCreate();
        cdo_vlist_copy_flag(vlistID2, vlistID0);

        varList2 = VarList(vlistID2);

        // if (Options::cdoVerbose) vlistPrint(vlistID2);

        taxisID2 = taxisDuplicate(taxisID1);
        numSteps = varList1.numSteps();

        numVars2 = varList2.numVars();

        if (numSteps == 1 && numFiles == 1 && varList2.numVaryingVars() == 0) numSteps = 0;

        numStepsOut = (numFiles == 1 && doTimeSel == false) ? numSteps : -1;
        if (operatorID == SELECT && SELINFO_NVAL(timestep) > 0) numStepsOut = SELINFO_NVAL(timestep);

        if (numStepsOut == 0 && numFiles > 1)
        {
          hasConstVars = false;
          for (int varID = 0; varID < numVars2; ++varID) vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
        }

        auto numStepValues = SELINFO_NVAL(timestep);
        // support for negative timestep values
        if (numStepValues > 0)
        {
          for (int i = 0; i < numStepValues; ++i)
          {
            int ptimestep;
            SELINFO_GET_VAL(timestep, i, &ptimestep);
            if (ptimestep < 0)
            {
              if (numFiles != 1) cdo_abort("Negative timesteps only supported with one input stream!");
              if (numSteps < 0 && cdo_assert_files_only())
              {
                int tsID = 0;
                while (cdo_stream_inq_timestep(streamID1, tsID)) tsID++;
                numSteps = tsID;
                if (Options::cdoVerbose) cdo_print("Found %d timesteps", numSteps);
              }
              if (numSteps > 0)
              {
                if (Options::cdoVerbose) cdo_print("timestep %d changed to %d", ptimestep, ptimestep + numSteps + 1);
                ptimestep += numSteps + 1;
                SELINFO_DEF_VAL(timestep, i, &ptimestep);
              }
            }
          }

          for (int i = 0; i < numStepValues; ++i)
          {
            int ptimestep;
            SELINFO_GET_VAL(timestep, i, &ptimestep);
            tsmax = std::max(tsmax, ptimestep);
          }
        }

        SELINFO_GET_VAL(startdate, 0, &startdate);
        SELINFO_GET_VAL(enddate, 0, &enddate);
        if (SELINFO_NVAL(startdate)) fstartdate = datestr_to_double(startdate, 0);
        if (SELINFO_NVAL(enddate)) fenddate = datestr_to_double(enddate, 1);
      }
      else { varList_compare(VarList(vlistID0), varList1); }

      if (numVars2 == 0)
      {
        cdo_warning("No variable selected!");
        return;
      }

      if (copyConstVars) varsData2.resize(numVars2);

      auto stopReading = false;
      int tsID1 = 0;
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID1, tsID1);
        if (numFields == 0) break;

        auto fstatus = (numSteps > 1) ? fileIdx + (tsID1 + 1.0) / numSteps : fileIdx + 1.0;
        progress.update(fstatus / numFiles);

        timestep++;
        auto copyTimestep = true;

        if (doTimeSel)
        {
          copyTimestep = false;

          if (operatorID == SELECT && SELINFO_NVAL(timestep) > 0)
          {
            if (timestep > tsmax)
            {
              stopReading = true;
              break;
            }
          }

          dtlist.taxis_inq_timestep(taxisID1, 0);
          auto vDateTime = dtlist.vDateTime(0);
          int second, ms;
          cdiDate_decode(vDateTime.date, &year, &month, &day);
          cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);
          (void) (season);  // unused

          if (year != last_year)
          {
            timestep_of_year = 0;
            last_year = year;
          }

          timestep_of_year++;

          if (SELINFO_CHECK(timestep)) copyTimestep = true;
          if (SELINFO_CHECK(timestep_of_year)) copyTimestep = true;

          if (!copyTimestep && SELINFO_NVAL(date) == 0 && SELINFO_NVAL(timestep) == 0 && SELINFO_NVAL(timestep_of_year) == 0
              && SELINFO_NVAL(dom) == 0)
          {
            auto lseason = (SELINFO_NVAL(season) == 0 || SELINFO_CHECK_SEASON(season, month));
            auto lyear = (SELINFO_NVAL(year) == 0 || SELINFO_CHECK(year));
            auto lmonth = (SELINFO_NVAL(month) == 0 || SELINFO_CHECK(month));
            auto lday = (SELINFO_NVAL(day) == 0 || SELINFO_CHECK(day));
            auto lhour = (SELINFO_NVAL(hour) == 0 || SELINFO_CHECK(hour));
            auto lminute = (SELINFO_NVAL(minute) == 0 || SELINFO_CHECK(minute));

            if (lseason && lyear && lmonth && lday && lhour && lminute) copyTimestep = true;
          }

          auto fdatetime = datetime_to_double(vDateTime);

          if (SELINFO_NVAL(enddate))
          {
            copyTimestep = (fdatetime <= fenddate);
            if (fdatetime > fenddate)
            {
              SELINFO_DEF_FLAG(enddate, 0, true);
              if (operatorID == SELECT)
              {
                stopReading = true;
                break;
              }
            }
          }

          if (SELINFO_NVAL(startdate))
          {
            copyTimestep = (fdatetime >= fstartdate);
            if (fdatetime >= fstartdate) SELINFO_DEF_FLAG(startdate, 0, true);
          }

          if (SELINFO_NVAL(date))
          {
            auto datetimeString = datetime_to_string(vDateTime);
            date = datetimeString.c_str();
            if (SELINFO_CHECK_DATE(date)) copyTimestep = true;
          }

          if (SELINFO_NVAL(dom))
          {
            auto domString = dom_to_string(vDateTime.date);
            dom = domString.c_str();
            if (SELINFO_CHECK_DATE(dom)) copyTimestep = true;
          }

          if (operatorID == DELETE) copyTimestep = !copyTimestep;

          if (copyTimestep && fileIdx == 0 && tsID1 == 0) copyConstVars = false;
        }

        if (copyTimestep)
        {
          cdo_taxis_copy_timestep(taxisID2, taxisID1);
          if (streamID2 == CDO_STREAM_UNDEF)
          {
            auto isLastTimestep = ((numFiles == 1) && (numStepsOut > 1) && (numStepsOut == (tsID1 + 1)));
            if (isLastTimestep && tsID2 == 0) numStepsOut = 1;
            vlistDefNtsteps(vlistID2, numStepsOut);
            streamID2 = cdo_open_write(numFiles);
            vlistDefTaxis(vlistID2, taxisID2);
            cdo_def_vlist(streamID2, vlistID2);
          }

          cdo_def_timestep(streamID2, tsID2);

          for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = cdo_inq_field(streamID1);
            if (vlistInqFlag(vlistID0, varID, levelID) == true)
            {
              auto const &var = varList1.vars[varID];

              if (hasConstVars && tsID2 > 0 && tsID1 == 0)
                if (var.isConstant) continue;

              auto varID2 = vlistFindVar(vlistID2, varID);
              auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
              if (copyConstVars && tsID2 == 0) write_const_vars(streamID2, varList2, varsData2);

              cdo_def_field(streamID2, varID2, levelID2);
              if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
              else
              {
                field.init(var);
                cdo_read_field(streamID1, field);
                cdo_write_field(streamID2, field);
              }
            }
          }

          if (copyConstVars && tsID2 == 0) write_const_vars(streamID2, varList2, varsData2);

          tsID2++;
        }
        else if (copyConstVars && fileIdx == 0 && tsID1 == 0)
        {
          for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = cdo_inq_field(streamID1);
            if (vlistInqFlag(vlistID0, varID, levelID) == true)
            {
              auto varID2 = vlistFindVar(vlistID2, varID);
              auto const &var = varList2.vars[varID2];
              if (var.isConstant)
              {
                auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                if (levelID == 0) varsData2[varID2].resize(var.gridsize * var.nlevels);
                size_t numMissVals;
                cdo_read_field(streamID1, &varsData2[varID2][var.gridsize * levelID2], &numMissVals);
              }
            }
          }
        }

        tsID1++;
      }

      cdo_stream_close(streamID1);

      if (stopReading) break;
    }

    SELINFO_CHECK_FLAG(timestep_of_year);
    SELINFO_CHECK_FLAG(timestep);
    SELINFO_CHECK_FLAG(year);
    SELINFO_CHECK_FLAG(month);
    SELINFO_CHECK_FLAG(day);
    SELINFO_CHECK_FLAG(hour);
    SELINFO_CHECK_FLAG(minute);
    SELINFO_CHECK_FLAG(startdate);
    //  SELINFO_CHECK_FLAG(enddate);
    SELINFO_CHECK_FLAG(season);
    SELINFO_CHECK_FLAG(date);
    SELINFO_CHECK_FLAG(dom);
  }

  void
  close() override
  {
    if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);

    vlistDestroy(vlistID0);
    vlistDestroy(vlistID2);

    if (tsID2 == 0) cdo_abort("No timesteps selected!");
  }
};
