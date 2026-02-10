/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Change     chcode          Change code number
      Change     chtabnum        Change GRIB1 parameter table number
      Change     chparam         Change parameter identifier
      Change     chname          Change variable or coordinate name
      Change     chlevel         Change level
      Change     chlevelc        Change level of one code
      Change     chlevelv        Change level of one variable
      Change     chltype         Change GRIB level type
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

using IntPairList = std::vector<std::pair<int, int>>;
using FltPairList = std::vector<std::pair<double, double>>;
using StringPairList = std::vector<std::pair<std::string, std::string>>;

static void
change_code(VarList const &varList1, int vlistID2, IntPairList const &intPairList)
{
  for (auto const &var1 : varList1.vars)
  {
    for (auto const &intPair : intPairList)
      if (var1.code == intPair.first) vlistDefVarCode(vlistID2, var1.ID, intPair.second);
  }
}

static void
change_tabnum(VarList const &varList1, int vlistID2, IntPairList const &intPairList)
{
  for (auto const &var1 : varList1.vars)
  {
    auto tabnum = tableInqNum(vlistInqVarTable(vlistID2, var1.ID));
    for (auto const &intPair : intPairList)
      if (tabnum == intPair.first)
      {
        auto tableID = tableDef(-1, intPair.second, nullptr);
        vlistDefVarTable(vlistID2, var1.ID, tableID);
      }
  }
}

static void
change_param(VarList const &varList1, int vlistID2, StringPairList const &stringPairList)
{
  for (auto const &var1 : varList1.vars)
  {
    if (Options::cdoVerbose)
    {
      int pnum, pcat, pdis;
      cdiDecodeParam(var1.param, &pnum, &pcat, &pdis);
      cdo_print("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);
    }
    auto paramStr = param_to_string(var1.param);
    for (auto const &stringPair : stringPairList)
      if (paramStr == stringPair.first) vlistDefVarParam(vlistID2, var1.ID, string_to_param(stringPair.second));
  }
}

static void
change_name(VarList const &varList1, int vlistID2, StringPairList const &stringPairList)
{
  int numPairs = stringPairList.size();

  std::vector<bool> namefound(numPairs, false);
  for (auto const &var1 : varList1.vars)
  {
    for (int i = 0; i < numPairs; ++i)
      if (var1.name == stringPairList[i].first)
      {
        namefound[i] = true;
        cdiDefKeyString(vlistID2, var1.ID, CDI_KEY_NAME, stringPairList[i].second.c_str());
        break;
      }
  }

  auto searchForGridName = false;
  for (int i = 0; i < numPairs; ++i)
    if (!namefound[i])
    {
      searchForGridName = true;
      break;
    }

  if (searchForGridName)
  {
    auto numGrids = vlistNumGrids(vlistID2);
    for (int index = 0; index < numGrids; ++index)
    {
      int gridID2 = -1;
      auto gridID1 = vlistGrid(vlistID2, index);
      auto xname = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_NAME);
      auto yname = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_NAME);
      auto xfound = false, yfound = false;
      for (int i = 0; i < numPairs; ++i)
      {
        if (!namefound[i])
        {
          if (xname == stringPairList[i].first)
          {
            xfound = true;
            namefound[i] = true;
            if (gridID2 == -1) gridID2 = gridDuplicate(gridID1);
            cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_NAME, stringPairList[i].second.c_str());
          }
        }
        if (!namefound[i])
        {
          if (yname == stringPairList[i].first)
          {
            yfound = true;
            namefound[i] = true;
            if (gridID2 == -1) gridID2 = gridDuplicate(gridID1);
            cdiDefKeyString(gridID2, CDI_YAXIS, CDI_KEY_NAME, stringPairList[i].second.c_str());
          }
        }

        if (xfound && yfound) break;
      }

      if (gridID2 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);
    }
  }

  auto searchForZaxisName = false;
  for (int i = 0; i < numPairs; ++i)
    if (!namefound[i])
    {
      searchForZaxisName = true;
      break;
    }

  if (searchForZaxisName)
  {
    auto numZaxes = vlistNumZaxis(vlistID2);
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID1 = vlistZaxis(vlistID2, index);
      auto varname = cdo::inq_key_string(zaxisID1, CDI_GLOBAL, CDI_KEY_NAME);
      for (int i = 0; i < numPairs; ++i)
      {
        if (!namefound[i])
        {
          if (varname == stringPairList[i].first)
          {
            namefound[i] = true;
            auto zaxisID2 = zaxisDuplicate(zaxisID1);
            cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, stringPairList[i].second.c_str());
            vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
            break;
          }
        }
      }
    }
  }

  for (int i = 0; i < numPairs; ++i)
    if (!namefound[i]) cdo_warning("Variable name %s not found!", stringPairList[i].first);
}

static void
change_unit(VarList const &varList1, int vlistID2, StringPairList const &stringPairList)
{
  for (auto const &var1 : varList1.vars)
  {
    for (auto const &stringPair : stringPairList)
      if (var1.units == stringPair.first) cdiDefKeyString(vlistID2, var1.ID, CDI_KEY_UNITS, stringPair.second.c_str());
  }
}

static void
change_level(int vlistID2, FltPairList const &fltPairList)
{
  auto numZaxes = vlistNumZaxis(vlistID2);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID1 = vlistZaxis(vlistID2, index);
    if (zaxisInqLevels(zaxisID1, nullptr))
    {
      auto numLevels = zaxisInqSize(zaxisID1);
      Varray<double> levels(numLevels);
      zaxisInqLevels(zaxisID1, &levels[0]);

      int nfound = 0;
      for (auto const &fltPair : fltPairList)
      {
        for (int k = 0; k < numLevels; ++k)
          if (std::fabs(levels[k] - fltPair.first) < 0.0001) nfound++;
      }

      if (nfound)
      {
        Varray<double> newlevels = levels;
        auto zaxisID2 = zaxisDuplicate(zaxisID1);
        for (auto const &fltPair : fltPairList)
        {
          for (int k = 0; k < numLevels; ++k)
            if (std::fabs(levels[k] - fltPair.first) < 0.001) newlevels[k] = fltPair.second;
        }

        zaxisDefLevels(zaxisID2, &newlevels[0]);
        vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
      }
    }
  }
}

static void
change_varLevel(int varID, int vlistID2, std::pair<double, double> const &levelPair)
{
  auto zaxisID1 = vlistInqVarZaxis(vlistID2, varID);
  if (zaxisInqLevels(zaxisID1, nullptr))
  {
    auto numLevels = zaxisInqSize(zaxisID1);
    Varray<double> levels(numLevels);
    zaxisInqLevels(zaxisID1, levels.data());

    int nfound = 0;
    for (int k = 0; k < numLevels; ++k)
      if (std::fabs(levels[k] - levelPair.first) < 0.0001) nfound++;

    if (nfound == 0) { cdo_abort("Level %g not found!", levelPair.first); }

    auto zaxisID2 = zaxisDuplicate(zaxisID1);
    for (int k = 0; k < numLevels; ++k)
      if (std::fabs(levels[k] - levelPair.first) < 0.001) levels[k] = levelPair.second;

    zaxisDefLevels(zaxisID2, &levels[0]);
    vlistChangeVarZaxis(vlistID2, varID, zaxisID2);
  }
}

static void
change_levelByCode(int varCode, VarList const &varList1, int vlistID2, std::pair<double, double> const &levelPair)
{
  for (auto const &var1 : varList1.vars)
  {
    if (var1.code == varCode)
    {
      change_varLevel(var1.ID, vlistID2, levelPair);
      return;
    }
  }

  cdo_abort("Code %d not found!", varCode);
}

static void
change_levelByName(std::string const &varName, VarList const &varList1, int vlistID2, std::pair<double, double> const &levelPair)
{
  for (auto const &var1 : varList1.vars)
  {
    if (var1.name == varName)
    {
      change_varLevel(var1.ID, vlistID2, levelPair);
      return;
    }
  }

  cdo_abort("Variable name %s not found!", varName);
}

static void
change_ltype(int vlistID2, IntPairList const &intPairList)
{
  auto numZaxes = vlistNumZaxis(vlistID2);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID1 = vlistZaxis(vlistID2, index);
    auto zaxisID2 = zaxisDuplicate(zaxisID1);
    int ltype = 0;
    cdiInqKeyInt(zaxisID1, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);

    for (auto const &intPair : intPairList)
    {
      if (intPair.first == ltype)
      {
        zaxisChangeType(zaxisID2, ZAXIS_GENERIC);
        cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, intPair.second);
        vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
      }
    }
  }
}

class Change : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Change",
    .operators = { { "chcode", 0, 0, "pairs of old and new code numbers", ChangeHelp },
                   { "chtabnum", 0, 0, "pairs of old and new GRIB1 table numbers", ChangeHelp },
                   { "chparam", 0, 0, "pairs of old and new parameter identifiers", ChangeHelp },
                   { "chname", 0, 0, "pairs of old and new variable names", ChangeHelp },
                   { "chunit", 0, 0, "pairs of old and new variable units", ChangeHelp },
                   { "chlevel", 0, 0, "pairs of old and new levels", ChangeHelp },
                   { "chlevelc", 0, 0, "code number, old and new level", ChangeHelp },
                   { "chlevelv", 0, 0, "variable name, old and new level", ChangeHelp },
                   { "chltype", 0, 0, "pairs of old and new level type", ChangeHelp } },
    .aliases = { { "chvar", "chname" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Change> registration = RegisterEntry<Change>(module);

  int CHCODE{}, CHTABNUM{}, CHPARAM{}, CHNAME{}, CHUNIT{}, CHLEVEL{}, CHLEVELC{}, CHLEVELV{}, CHLTYPE{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};
  int operatorID{};

public:
  void
  process_int_args(std::vector<std::string> const &argList, int vlistID2)
  {
    int numArgs = argList.size();
    if (numArgs % 2) cdo_abort("Odd number of input arguments!");
    int numPairs = numArgs / 2;
    IntPairList intPairList(numPairs);
    for (int i = 0; i < numPairs; ++i)
    {
      intPairList[i].first = parameter_to_int(argList[i * 2]);
      intPairList[i].second = parameter_to_int(argList[i * 2 + 1]);
    }

    if (operatorID == CHCODE) { change_code(varList1, vlistID2, intPairList); }
    else if (operatorID == CHTABNUM) { change_tabnum(varList1, vlistID2, intPairList); }
    else if (operatorID == CHLTYPE) { change_ltype(vlistID2, intPairList); }
  }

  void
  process_string_args(std::vector<std::string> const &argList, int vlistID2)
  {
    int numArgs = argList.size();
    if (numArgs % 2) cdo_abort("Odd number of input arguments!");
    int numPairs = numArgs / 2;
    StringPairList stringPairList(numPairs);
    for (int i = 0; i < numPairs; ++i)
    {
      stringPairList[i].first = parameter_to_word(argList[i * 2]);
      stringPairList[i].second = parameter_to_word(argList[i * 2 + 1]);
    }

    if (operatorID == CHPARAM) { change_param(varList1, vlistID2, stringPairList); }
    else if (operatorID == CHNAME) { change_name(varList1, vlistID2, stringPairList); }
    else if (operatorID == CHUNIT) { change_unit(varList1, vlistID2, stringPairList); }
  }

  void
  process_flt_args(std::vector<std::string> const &argList, int vlistID2)
  {
    int numArgs = argList.size();
    if (numArgs % 2) cdo_abort("Odd number of input arguments!");
    int numPairs = numArgs / 2;
    FltPairList fltPairList(numPairs);
    for (int i = 0; i < numPairs; ++i)
    {
      fltPairList[i].first = parameter_to_double(argList[i * 2]);
      fltPairList[i].second = parameter_to_double(argList[i * 2 + 1]);
    }

    change_level(vlistID2, fltPairList);
  }

  void
  process_flt_pair(std::vector<std::string> const &argList, int vlistID2)
  {
    operator_check_argc(3);

    std::pair<double, double> levelPair;
    levelPair.first = parameter_to_double(argList[1]);
    levelPair.second = parameter_to_double(argList[2]);

    if (operatorID == CHLEVELC) { change_levelByCode(parameter_to_int(argList[0]), varList1, vlistID2, levelPair); }
    else if (operatorID == CHLEVELV) { change_levelByName(parameter_to_word(argList[0]), varList1, vlistID2, levelPair); }
  }

  void
  init() override
  {
    CHCODE = module.get_id("chcode");
    CHTABNUM = module.get_id("chtabnum");
    CHPARAM = module.get_id("chparam");
    CHNAME = module.get_id("chname");
    CHUNIT = module.get_id("chunit");
    CHLEVEL = module.get_id("chlevel");
    CHLEVELC = module.get_id("chlevelc");
    CHLEVELV = module.get_id("chlevelv");
    CHLTYPE = module.get_id("chltype");

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto numArgs = cdo_operator_argc();
    if (numArgs < 2) cdo_abort("Too few arguments!");
    auto const &argList = cdo_get_oper_argv();

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    if (operatorID == CHCODE || operatorID == CHTABNUM || operatorID == CHLTYPE) { process_int_args(argList, vlistID2); }
    else if (operatorID == CHPARAM || operatorID == CHNAME || operatorID == CHUNIT) { process_string_args(argList, vlistID2); }
    else if (operatorID == CHLEVEL) { process_flt_args(argList, vlistID2); }
    else if (operatorID == CHLEVELC || operatorID == CHLEVELV) { process_flt_pair(argList, vlistID2); }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
