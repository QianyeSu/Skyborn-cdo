/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_varlist.h"
#include "cdo_cdi_wrapper.h"
#include "cdo_output.h"
#include "util_string.h"
#include "stdnametable.h"
#include "cdo_vlist.h"

static bool
is_int16_type(int dataType)
{
  return (dataType == CDI_DATATYPE_UINT16 || dataType == CDI_DATATYPE_INT16);
}

static bool
is_int8_type(int dataType)
{
  return (dataType == CDI_DATATYPE_UINT8 || dataType == CDI_DATATYPE_INT8);
}

static bool
is_float_type(int dataType)
{
  return (dataType == CDI_DATATYPE_FLT32 || dataType == CDI_DATATYPE_CPX32);
}

void
cdoVars_init(CdoVars &cdoVars, int vlistID)
{
  auto numVars = vlistNvars(vlistID);
  cdoVars.resize(numVars);

  for (int varID = 0; varID < numVars; ++varID)
  {
    auto &var = cdoVars[varID];
    var.ID = varID;
    var.name = cdo::inq_var_name(vlistID, varID);
    var.longname = cdo::inq_var_longname(vlistID, varID);
    var.units = cdo::inq_var_units(vlistID, varID);
    var.stdname = cdo::inq_key_string(vlistID, varID, CDI_KEY_STDNAME);
    var.gridID = vlistInqVarGrid(vlistID, varID);
    var.zaxisID = vlistInqVarZaxis(vlistID, varID);
    var.timeType = vlistInqVarTimetype(vlistID, varID);
    var.stepType = vlistInqVarTsteptype(vlistID, varID);
    var.gridType = gridInqType(var.gridID);
    var.gridsize = gridInqSize(var.gridID);
    var.zaxisType = zaxisInqType(var.zaxisID);
    var.nlevels = zaxisInqSize(var.zaxisID);
    var.dataType = vlistInqVarDatatype(vlistID, varID);
    var.missval = vlistInqVarMissval(vlistID, varID);
    var.code = vlistInqVarCode(vlistID, varID);
    var.param = vlistInqVarParam(vlistID, varID);
    var.nwpv = vlistInqVarNumber(vlistID, varID);
    var.isConstant = (var.timeType == TIME_CONSTANT);
    double addOffset = 0.0, scaleFactor = 1.0;
    auto haveAddOffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addOffset) == CDI_NOERR);
    auto haveScaleFactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scaleFactor) == CDI_NOERR);
    var.isPacked = (haveAddOffset || haveScaleFactor);
    if (haveAddOffset) var.addOffset = addOffset;
    if (haveScaleFactor) var.scaleFactor = scaleFactor;

    if (Options::CDO_Memtype == MemType::Native)
    {
      auto useFloatType = (var.dataType == CDI_UNDEFID) || is_float_type(var.dataType) || is_int8_type(var.dataType)
                          || (is_int16_type(var.dataType) && !var.isPacked);
      var.memType = useFloatType ? MemType::Float : MemType::Double;
    }
    else { var.memType = Options::CDO_Memtype; }
  }
}

void
varList_set_memtype(VarList &varList, MemType memType)
{
  for (auto &var : varList.vars) var.memType = memType;
}

void
varList_set_unique_memtype(VarList &varList)
{
  auto numVars = varList.numVars();
  if (numVars)
  {
    auto memtype = varList.vars[0].memType;
    int varID;
    for (varID = 1; varID < numVars; ++varID)
    {
      if (varList.vars[varID].memType != memtype) break;
    }
    if (varID < numVars) varList_set_memtype(varList, MemType::Double);
  }
}

void
VarList::set_num_const_vars(CdoVars const &cdoVars)
{
  m_numConstVars = std::ranges::count_if(cdoVars, [](auto const &var) { return (var.timeType == TIME_CONSTANT); });
}

VarList::VarList(int _vlistID) : vlistID(_vlistID)
{
  cdoVars_init(vars, _vlistID);
  m_maxFields = vlistNumFields(_vlistID);
  m_numSteps = vlistNtsteps(_vlistID);
  m_numZaxes = vlistNumZaxis(_vlistID);
  m_numGrids = vlistNumGrids(_vlistID);
  set_num_const_vars(vars);
  set_num_varying_vars(vars);
  m_gridsizeMax = vlistGridsizeMax(_vlistID);
}

void
VarList::set_num_varying_vars(CdoVars const &cdoVars)
{
  m_numVaryingVars = std::ranges::count_if(cdoVars, [](auto const &var) { return (var.timeType == TIME_VARYING); });
}

static bool
table_is_used(VarList const &varList)
{
  auto numVars = varList.numVars();
  auto useTable = false;
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto tableNum = tableInqNum(vlistInqVarTable(varList.vlistID, varID));
    if (tableNum > 0 && tableNum < 255)
    {
      useTable = true;
      break;
    }
  }

  return useTable;
}

static int
get_code(VarIDs const &varIDs, GribCodes const &gribCodes, std::string const &varname, std::string const &stdname)
{
  int code = -1;
  //                                  ECHAM                 ECMWF
  // clang-format off
  if      (-1 == varIDs.sgeopotID && (varname == "geosp" || varname == "z")) code = gribCodes.geopot;
  else if (-1 == varIDs.taID      && (varname == "st"    || varname == "t")) code = gribCodes.ta;
  else if (-1 == varIDs.psID      && (varname == "aps"   || varname == "sp")) code = gribCodes.ps;
  else if (-1 == varIDs.psID      &&  varname == "ps") code = gribCodes.ps;
  else if (-1 == varIDs.lnpsID    && (varname == "lsp"   || varname == "lnsp")) code = gribCodes.lsp;
  else if (-1 == varIDs.lnpsID2   &&  varname == "lnps") code = 777;
  else if (-1 == varIDs.geopotID  &&  stdname == "geopotential_full") code = gribCodes.geopot;
  else if (-1 == varIDs.taID      &&  varname == "t") code = gribCodes.ta;
  else if (-1 == varIDs.husID     &&  varname == "q") code = gribCodes.hus;
  // else if (varname == "clwc") code = 246;
  // else if (varname == "ciwc") code = 247;
  // clang-format on
  return code;
}

VarIDs
varList_search_varIDs(VarList const &varList, int numFullLevels)
{
  VarIDs varIDs;

  auto useTable = table_is_used(varList);
  if (Options::cdoVerbose && useTable) cdo_print("Using code tables!");

  char paramstr[32];
  GribCodes gribCodes;

  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto &var = varList.vars[varID];
    auto numLevels = var.nlevels;
    auto instNum = institutInqCenter(vlistInqVarInstitut(varList.vlistID, varID));
    auto tableNum = tableInqNum(vlistInqVarTable(varList.vlistID, varID));

    auto code = var.code;

    cdiParamToString(var.param, paramstr, sizeof(paramstr));
    int pnum, pcat, pdis;
    cdiDecodeParam(var.param, &pnum, &pcat, &pdis);
    if (pdis >= 0 && pdis < 255) code = -1;

    if (useTable)
    {
      if (tableNum == 2) { gribCodes = wmo_gribcodes(); }
      else if (tableNum == 128 || tableNum == 0 || tableNum == 255) { gribCodes = echam_gribcodes(); }
      //  KNMI: HIRLAM model version 7.2 uses tableNum=1    (LAMH_D11*)
      //  KNMI: HARMONIE model version 36 uses tableNum=1   (grib*) (opreational NWP version)
      //  KNMI: HARMONIE model version 38 uses tableNum=253 (grib,grib_md) and tableNum=1 (grib_sfx) (research version)
      else if (tableNum == 1 || tableNum == 253) { gribCodes = hirlam_harmonie_gribcodes(); }
    }
    else { gribCodes = echam_gribcodes(); }

    if (Options::cdoVerbose)
      cdo_print("Center=%d  TableNum=%d  Code=%d  Param=%s  Varname=%s  varID=%d", instNum, tableNum, code, paramstr, var.name,
                varID);

    if (code <= 0 || code == 255)
    {
      auto varname = string_to_lower(var.name);
      auto stdname = string_to_lower(var.stdname);
      code = stdname_to_echamcode(stdname);
      if (code == -1) code = get_code(varIDs, gribCodes, varname, stdname);
    }

    // clang-format off
    if      (code == gribCodes.geopot  && numLevels == 1)                 varIDs.sgeopotID = varID;
    else if (code == gribCodes.geopot  && numLevels == numFullLevels)     varIDs.geopotID = varID;
    else if (code == gribCodes.ta      && numLevels == numFullLevels)     varIDs.taID = varID;
    else if (code == gribCodes.ps      && numLevels == 1)                 varIDs.psID = varID;
    else if (code == gribCodes.lsp     && numLevels == 1)                 varIDs.lnpsID = varID;
    else if (code == 777               && numLevels == 1)                 varIDs.lnpsID2 = varID;
    else if (code == gribCodes.gheight && numLevels == numFullLevels)     varIDs.gheightID = varID;
    else if (code == gribCodes.gheight && numLevels == numFullLevels + 1) varIDs.gheightID = varID;
    else if (code == gribCodes.hus     && numLevels == numFullLevels)     varIDs.husID = varID;
    // else if (code == 246 && nlevels == numFullLevels) varIDs.clwcID = varID;
    // else if (code == 247 && nlevels == numFullLevels) varIDs.ciwcID = varID;
    // clang-format on
  }

  return varIDs;
}

void
varList_map(VarList const &varList1, VarList const &varList2, MapFlag mapFlag, std::map<int, int> &mapOfVarIDs)
{
  auto cmpFlag{ CmpVarList::All };
  auto numVars1 = varList1.numVars();

  if (mapFlag == MapFlag::Right)
  {
    for (auto const &var2 : varList2.vars)
    {
      int varID1 = 0;
      for (; varID1 < numVars1; ++varID1)
      {
        if (varList1.vars[varID1].name == var2.name) break;
      }
      if (varID1 == numVars1) { cdo_abort("Variable %s not found in first input stream!", var2.name); }
      else { mapOfVarIDs[varID1] = var2.ID; }
    }
  }
  else
  {
    for (auto const &var1 : varList1.vars)
    {
      auto numVars2 = varList2.numVars();
      int varID2 = 0;
      for (; varID2 < numVars2; ++varID2)
      {
        if (var1.name == varList2.vars[varID2].name) break;
      }
      if (varID2 == numVars2)
      {
        if (mapFlag == MapFlag::Intersect) continue;
        cdo_abort("Variable %s not found in second input stream!", var1.name);
      }
      else { mapOfVarIDs[var1.ID] = varID2; }
    }
  }

  if (mapOfVarIDs.empty()) cdo_abort("No variable found that occurs in both streams!");

  if (Options::cdoVerbose)
    for (int varID1 = 0; varID1 < numVars1; ++varID1)
    {
      auto const &var1 = varList1.vars[varID1];
      auto const &it = mapOfVarIDs.find(varID1);
      if (it != mapOfVarIDs.end())
        cdo_print("Variable %d:%s mapped to %d:%s", varID1, var1.name, it->second, varList2.vars[it->second].name);
    }

  if (mapOfVarIDs.size() > 1)
  {
    auto varID2 = mapOfVarIDs.begin()->second;
    for (auto it = ++mapOfVarIDs.begin(); it != mapOfVarIDs.end(); ++it)
    {
      if (it->second < varID2)
        cdo_abort("Variable names must be sorted, use CDO option --sortname to sort the parameter by name (NetCDF only)!");

      varID2 = it->second;
    }
  }

  for (auto it = mapOfVarIDs.begin(); it != mapOfVarIDs.end(); ++it)
  {
    auto varID1 = it->first;
    auto varID2 = it->second;
    auto const &var1 = varList1.vars[varID1];
    auto const &var2 = varList2.vars[varID2];

    if (cmpFlag & CmpVarList::GridSize)
    {
      if (var1.gridsize != var2.gridsize) cdo_abort("Grid size of the input fields do not match!");
    }

    if (cmpFlag & CmpVarList::NumLevels)
    {
      if (zaxis_check_levels(var1.zaxisID, var2.zaxisID) != 0) break;
    }

    if ((cmpFlag & CmpVarList::Grid) && (varID1 == mapOfVarIDs.begin()->first)) { cdo_compare_grids(var1.gridID, var2.gridID); }
  }
}

int
varList_get_psvarid(VarList const &varList, int zaxisID)
{
  auto psname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_PSNAME);
  if (psname.size())
  {
    for (auto const &var : varList.vars)
    {
      if (var.name == psname) return var.ID;
    }
    if (Options::cdoVerbose) cdo_warning("Surface pressure variable not found - %s", psname);
  }

  return -1;
}

static void
cdoVars_check_names(CdoVars const &cdoVars1, CdoVars const &cdoVars2)
{
  int numVars = cdoVars1.size();

  std::vector<std::string> names1(numVars);
  std::vector<std::string> names2(numVars);
  for (int varID = 0; varID < numVars; ++varID) names1[varID] = cdoVars1[varID].name;
  for (int varID = 0; varID < numVars; ++varID) names2[varID] = cdoVars2[varID].name;

  std::ranges::sort(names1);
  std::ranges::sort(names2);

  if (names1 == names2) cdo_print("Use CDO option --sortname to sort the parameter by name (NetCDF only)!");
}

static void
cdoVars_print_missing_vars(CdoVars const &cdoVars1, CdoVars const &cdoVars2)
{
  int numVars1 = cdoVars1.size();
  int numVars2 = cdoVars2.size();

  if (numVars1 > numVars2)
  {
    for (int varID1 = 0; varID1 < numVars1; ++varID1)
    {
      int varID2 = 0;
      for (; varID2 < numVars2; ++varID2)
      {
        if (cdoVars1[varID1].name == cdoVars2[varID2].name) break;
      }
      if (varID2 == numVars2) cdo_print("Variable %s not found in second input stream!", cdoVars1[varID1].name);
    }
  }
  else
  {
    for (int varID2 = 0; varID2 < numVars2; ++varID2)
    {
      int varID1 = 0;
      for (; varID1 < numVars1; ++varID1)
      {
        if (cdoVars1[varID1].name == cdoVars2[varID2].name) break;
      }
      if (varID1 == numVars1) cdo_print("Variable %s not found in first input stream!", cdoVars2[varID2].name);
    }
  }
}

static int
cdoVars_numFields(CdoVars const &cdoVars)
{
  int numFields = 0;
  for (int varID = 0, numVars = cdoVars.size(); varID < numVars; ++varID) numFields += cdoVars[varID].nlevels;
  return numFields;
}

void
varList_compare(VarList const &varList1, VarList const &varList2, int cmpFlag)
{
  auto doCheckNames = false;

  auto numVars = varList1.numVars();
  if (numVars != varList2.numVars())
  {
    cdoVars_print_missing_vars(varList1.vars, varList2.vars);
    cdo_abort("Input streams have different number of variables per timestep!");
  }

  if (cdoVars_numFields(varList1.vars) != cdoVars_numFields(varList2.vars))
    cdo_abort("Input streams have different number of %s per timestep!", (numVars == 1) ? "layers" : "fields");

  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var1 = varList1.vars[varID];
    auto const &var2 = varList2.vars[varID];
    if (numVars > 1)
    {
      if (cmpFlag & CmpVarList::Name)
      {
        if (string_to_lower(var1.name) != string_to_lower(var2.name))
        {
          cdo_warning("Input streams have different parameter names!");
          doCheckNames = true;
          cmpFlag = cmpFlag ^ CmpVarList::Name;
        }
      }
    }

    if (cmpFlag & CmpVarList::GridSize)
    {
      if (var1.gridsize != var2.gridsize) { cdo_abort("Grid size of the input field '%s' do not match!", var1.name); }
    }

    if (cmpFlag & CmpVarList::NumLevels)
    {
      if (zaxis_check_levels(var1.zaxisID, var2.zaxisID) != 0) break;
    }
  }

  if (cmpFlag & CmpVarList::Grid) { cdo_compare_grids(varList1.vars[0].gridID, varList2.vars[0].gridID); }

  if (doCheckNames) cdoVars_check_names(varList1.vars, varList2.vars);
}

void
vlist_compare(int vlistID1, int vlistID2, int cmpFlag)
{
  VarList varList1(vlistID1);
  VarList varList2(vlistID2);
  varList_compare(varList1, varList2, cmpFlag);
}
