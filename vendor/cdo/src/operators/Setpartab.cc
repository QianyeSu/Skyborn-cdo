/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setpartab  setpartab       Set parameter table
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "util_string.h"
#include "table.h"
#include "param_conversion.h"
#include "cdo_cmor.h"
#include "pmlist.h"
#include "convert_units.h"
#include "util_files.h"
#include "cdi_lockedIO.h"
#include "parse_literals.h"

enum pt_mode_t
{
  CODE_NUMBER,
  PARAMETER_ID,
  VARIABLE_NAME,
  STANDARD_NAME
};

void cdo_define_var_units(CmorVar &cmorVar, int vlistID2, int varID, std::string const &units);

void
mapvar(int vlistID, int varID, const KeyValues &kv, CmorVar &cmorVar, bool &hasValidMin, bool &hasValidMax, int ptab,
       bool isnPtmodeName)
{
  const auto key = string_to_lower(kv.key);
  auto const &value = kv.values[0];
  auto lv1 = (kv.nvalues == 1);

  // printf("key=%s  value=>%s<\n", key.c_str(), value.c_str());

  // clang-format off
  if (cmorVar.name == "cdocmor" )
  {
    if ((key == "cn") || (key == "cmor_name"))
    {
      auto name = cdo::inq_var_name(vlistID, varID);
      if (name[0] != 0) cdiDefAttTxt(vlistID, varID, "original_name", (int) name.size(), name.c_str());
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(value.c_str()));
    }
    else if (key == "character_axis")
      cdiDefAttTxt(vlistID, varID, "character_axis", (int) value.size(), value.c_str());
    else if (key == "z_axis")
      cdiDefAttTxt(vlistID, varID, "z_axis", (int) value.size(), value.c_str());
    else if (key == "variable_comment")
      cdiDefAttTxt(vlistID, varID, "variable_comment", (int) value.size(), value.c_str());
    else if (key == "positive")
    {
      if (!isspace(value[0])) cdiDefAttTxt(vlistID, varID, "positive", (int) value.size(), value.c_str());
    }
    else
    {
      if (Options::cdoVerbose) cdo_print("In applying the mapping table:\n          Key: '%s' is ignored.", key);
    }
  }
  if (key == "standard_name")
  {
    if (not lv1) cdo_abort("%s can only have one string value!", key);
    cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, value.c_str());
  }
  else if (key == "long_name")
  {
    if (not lv1) cdo_abort("%s can only have one string value!", key);
    cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, value.c_str());
  }
  else if (key == "units")
  {
    if (not lv1) cdo_abort("%s can only have one string value!", key);
    cdo_define_var_units(cmorVar, vlistID, varID, value);
  }
  else if (key == "filterspec")
  {
    if (not lv1) cdo_abort("%s can only have one string value!", key);
    cdiDefKeyString(vlistID, varID, CDI_KEY_FILTERSPEC, value.c_str());
  }
  else if (key == "name")
  {
    if (isnPtmodeName)
    {
      if (not lv1) cdo_abort("%s can only have one string value!", key);
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(value.c_str()));
    }
  }
  else if (key == "out_name")
  {
    if (not lv1) cdo_abort("%s can only have one string value!", key);
    auto outname = parameter_to_word(value);
    if (cmorVar.name != outname)
    {
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, outname.c_str());
      cdiDefAttTxt(vlistID, varID, "original_name", (int) cmorVar.name.size(), cmorVar.name.c_str());
    }
  }
  else if (lv1 && key == "param")
    vlistDefVarParam(vlistID, varID, string_to_param(parameter_to_word(value)));
  else if (lv1 && key == "out_param")
    vlistDefVarParam(vlistID, varID, string_to_param(parameter_to_word(value)));
  else if (lv1 && key == "code")
    vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter_to_int(value), ptab, 255));
  else if (lv1 && key == "out_code")
    vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter_to_int(value), ptab, 255));
  else if (lv1 && key == "uvRelativeToGrid")
    cdiDefKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, parameter_to_bool(value));
  else if (lv1 && key == "comment")
    cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
  else if (lv1 && key == "chunktype")
    ;
  else if (lv1 && key == "cell_methods")
    cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
  else if (lv1 && key == "cell_measures")
    cdiDefAttTxt(vlistID, varID, key.c_str(), (int) value.size(), value.c_str());
  else if (lv1 && key == "delete")
    cmorVar.remove = parameter_to_bool(value);
  else if (lv1 && key == "convert")
    cmorVar.convert = parameter_to_bool(value);
  else if (lv1 && key == "factor")
  {
    cmorVar.applyFactor = true;
    cmorVar.factor = parameter_to_double(value);
    if (Options::cdoVerbose) cdo_print("%s - scale factor %g", cmorVar.name, cmorVar.factor);
  }
  else if (lv1 && (key == "missval" || key == "missing_value"))
  {
    auto missval = parameter_to_double(value);
    auto missvalOld = vlistInqVarMissval(vlistID, varID);
    if (!fp_is_equal(missval, missvalOld))
    {
      if (Options::cdoVerbose) cdo_print("%s - change missval from %g to %g", cmorVar.name, missvalOld, missval);
      cmorVar.changeMissval = true;
      cmorVar.missvalOld = missvalOld;
      vlistDefVarMissval(vlistID, varID, missval);
    }
  }
  else if (lv1 && key == "valid_min")
  {
    hasValidMin = true;
    cmorVar.valid_min = parameter_to_double(value);
  }
  else if (lv1 && key == "valid_max")
  {
    hasValidMax = true;
    cmorVar.valid_max = parameter_to_double(value);
  }
  else if (lv1 && key == "ok_min_mean_abs")
  {
    cmorVar.check_min_mean_abs = true;
    cmorVar.ok_min_mean_abs = parameter_to_double(value);
  }
  else if (lv1 && key == "ok_max_mean_abs")
  {
    cmorVar.check_max_mean_abs = true;
    cmorVar.ok_max_mean_abs = parameter_to_double(value);
  }
  else if (lv1 && (key == "datatype" || key == "type"))
  {
    auto datatype = cdo::str_to_datatype(parameter_to_word(value));
    if (datatype != -1) vlistDefVarDatatype(vlistID, varID, datatype);
  }
  else if (lv1 && key == "dimensions") {}
  else
  {
    auto const &values = kv.values;
    auto const &rvalue = kv.values[0];
    auto nvalues = kv.nvalues;
    if (nvalues == 1 && rvalue.empty()) nvalues = 0;

    auto dtype = literals_find_datatype(nvalues, values);

    if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
      {
        std::vector<int> ivals(nvalues);
        for (int i = 0; i < nvalues; ++i) ivals[i] = literal_to_int(values[i]);
        cdiDefAttInt(vlistID, varID, key.c_str(), dtype, nvalues, ivals.data());
      }
    else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)
      {
        std::vector<double> dvals(nvalues);
        for (int i = 0; i < nvalues; ++i) dvals[i] = literal_to_double(values[i]);
        cdiDefAttFlt(vlistID, varID, key.c_str(), dtype, nvalues, dvals.data());
      }
    else { cdiDefAttTxt(vlistID, varID, key.c_str(), (int) rvalue.size(), rvalue.c_str()); }
  }
  // clang-format on
}

static void
search_global_missval(PMList &pmlist, bool &hasMissvals, double &missval)
{
  const std::vector<std::string> hentry = { "Header" };
  auto kvlist = pmlist.getKVListVentry(hentry);
  if (kvlist)
  {
    auto kv = kvlist->search("missing_value");
    if (kv && kv->nvalues > 0)
    {
      hasMissvals = true;
      missval = parameter_to_double(kv->values[0]);
    }
  }
}

static void
apply_parameterList(pt_mode_t ptmode, PMList &pmlist, int vlistID2, std::vector<CmorVar> &cmorVars)
{
  // search for global missing value
  auto hasMissvals = false;
  double missval = 0.0;
  search_global_missval(pmlist, hasMissvals, missval);

  const std::vector<std::string> ventry = { "variable_entry", "parameter" };
  char valstr[CDI_MAX_NAME];
  char paramstr[32];
  int codenum = 0;

  VarList varList2(vlistID2);

  int numVarsFound = 0;
  int numVars = cmorVars.size();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto &cmorVar = cmorVars[varID];
    cmorVar.name = cdo::inq_var_name(vlistID2, varID);

    auto const &var2 = varList2.vars[varID];
    if (hasMissvals)
    {
      if (fp_is_not_equal(missval, var2.missval))
      {
        cmorVar.changeMissval = true;
        cmorVar.missvalOld = var2.missval;
        vlistDefVarMissval(vlistID2, varID, missval);
      }
    }

    const KVList *kvlist = nullptr;
    if (ptmode == CODE_NUMBER)
    {
      codenum = var2.code;
      std::snprintf(valstr, sizeof(valstr), "%d", codenum);
      kvlist = pmlist.searchKVListVentry("code", valstr, ventry);
      if (kvlist)
      {
        auto tableID = vlistInqVarTable(vlistID2, varID);
        auto tabnum = tableInqNum(tableID);
        int levtype = 0;
        cdiInqKeyInt(var2.zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &levtype);
        auto table = tabnum;
        auto ltype = levtype;
        {
          auto kv = kvlist->search("table");
          if (kv && kv->nvalues == 1) table = parameter_to_int(kv->values[0]);
        }
        {
          auto kv = kvlist->search("ltype");
          if (kv && kv->nvalues == 1) ltype = parameter_to_int(kv->values[0]);
        }
        if (!(tabnum == table && levtype == ltype)) kvlist = nullptr;
      }
    }
    else if (ptmode == PARAMETER_ID)
    {
      auto param = var2.param;
      cdiParamToString(param, paramstr, sizeof(paramstr));
      std::snprintf(valstr, sizeof(valstr), "%s", paramstr);
      kvlist = pmlist.searchKVListVentry("param", valstr, ventry);
      if (kvlist)
      {
        int levtype = 0;
        cdiInqKeyInt(var2.zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &levtype);
        auto kv = kvlist->search("ltype");
        auto ltype = (kv && kv->nvalues == 1) ? parameter_to_int(kv->values[0]) : levtype;
        if (levtype != ltype) kvlist = nullptr;
      }
    }
    else if (ptmode == VARIABLE_NAME) { kvlist = pmlist.searchKVListVentry("name", cmorVar.name, ventry); }

    if (kvlist)
    {
      numVarsFound++;
      int pnum, ptab, pdum;
      cdiDecodeParam(var2.param, &pnum, &ptab, &pdum);

      auto hasValidMin = false;
      auto hasValidMax = false;
      for (auto const &kv : *kvlist)
      {
        mapvar(vlistID2, varID, kv, cmorVar, hasValidMin, hasValidMax, ptab, (ptmode != VARIABLE_NAME));
      }
      if (hasValidMin && hasValidMax) cmorVar.checkValid = true;
    }
    else if (Options::cdoVerbose)
    {
      // clang-format off
      if      (ptmode == CODE_NUMBER)   cdo_print("Code number %d not found in parameter table!", codenum);
      else if (ptmode == PARAMETER_ID)  cdo_print("Parameter ID %s not found in parameter table!", paramstr);
      else if (ptmode == VARIABLE_NAME) cdo_print("Variable %s not found in parameter table!", cmorVar.name);
      // clang-format on
    }
  }

  if (numVarsFound == 0)
  {
    // clang-format off
    if      (ptmode == CODE_NUMBER)   cdo_warning("None of the input variables has a code number that matches the entries in the parameter table!");
    else if (ptmode == PARAMETER_ID)  cdo_warning("None of the input variables has a parameter ID that matches the entries in the parameter table!");
    else if (ptmode == VARIABLE_NAME) cdo_warning("None of the input variables has a name that matches the entries in the parameter table!");
    // clang-format on
  }
}

static void
apply_codetable(int vlistID2, int tableID)
{
  char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  auto numVars = vlistNvars(vlistID2);
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto param = vlistInqVarParam(vlistID2, varID);
    int pdis, pcat, pnum;
    cdiDecodeParam(param, &pnum, &pcat, &pdis);
    if (pdis == 255)
    {
      auto code = pnum;
      int ltype = 0;
      cdiInqKeyInt(vlistInqVarZaxis(vlistID2, varID), CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
      name[0] = 0;
      longname[0] = 0;
      units[0] = 0;
      tableInqEntry(tableID, code, ltype, name, longname, units);
      if (name[0])
      {
        cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, name);
        if (longname[0]) cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, longname);
        if (units[0]) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, units);
      }
    }
    vlistDefVarTable(vlistID2, varID, tableID);
  }
}

static void
get_tableformat(std::string const &partab, int &tableFormat)
{
  if (FileUtils::file_exists(partab))
  {
    auto fobj = c_fopen(partab, "r");
    if (fobj.get() != nullptr)
    {
      std::fseek(fobj.get(), 0L, SEEK_END);
      auto fsize = (size_t) std::ftell(fobj.get());
      std::vector<char> parbuf(fsize + 1);
      std::fseek(fobj.get(), 0L, SEEK_SET);
      std::fread(parbuf.data(), fsize, 1, fobj.get());
      parbuf[fsize] = 0;
      std::fseek(fobj.get(), 0L, SEEK_SET);

      if (atoi(parbuf.data()) == 0) { tableFormat = 1; }
    }
  }
}

class Setpartab : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setpartab",
    .operators = { { "setcodetab", 0, 0, "parameter code table name", SetHelp },
                   { "setpartabc", 0, 0, "parameter table name", SetpartabHelp },
                   { "setpartabp", 0, 0, "parameter table name", SetpartabHelp },
                   { "setpartabn", 0, 0, "parameter table name", SetpartabHelp } },
    .aliases = { { "setpartab", "setcodetab" }, { "setpartabv", "setpartabn" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setpartab> registration = RegisterEntry<Setpartab>();

private:
  bool deleteVars = false;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  std::vector<CmorVar> cmorVars{};

  VarList varList2{};

public:
  void
  init() override
  {
    auto SETCODETAB = module.get_id("setcodetab");
    auto SETPARTABC = module.get_id("setpartabc");
    auto SETPARTABP = module.get_id("setpartabp");
    auto SETPARTABN = module.get_id("setpartabn");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

    auto convertData = false;
    if (cdo_operator_argc() == 2)
    {
      if (cdo_operator_argv(1) == "convert")
        convertData = true;
      else
        cdo_abort("Unknown parameter: >%s<", cdo_operator_argv(1));
    }

    if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");

    pt_mode_t ptmode = CODE_NUMBER;
    // clang-format off
    if      (operatorID == SETCODETAB) ptmode = CODE_NUMBER;
    else if (operatorID == SETPARTABC) ptmode = CODE_NUMBER;
    else if (operatorID == SETPARTABP) ptmode = PARAMETER_ID;
    else if (operatorID == SETPARTABN) ptmode = VARIABLE_NAME;
    // clang-format on

    int tableID = -1;
    int tableFormat = 0;
    if (ptmode == CODE_NUMBER)
    {
      auto const &partab = cdo_operator_argv(0);
      get_tableformat(partab, tableFormat);
      if (tableFormat == 0) tableID = cdo::define_table(partab);
    }
    else if (ptmode == PARAMETER_ID) { tableFormat = 1; }
    else if (ptmode == VARIABLE_NAME) { tableFormat = 1; }

    if (Options::cdoVerbose) cdo_print("Table format version %d", tableFormat);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);
    // vlistPrint(vlistID2);

    auto numVars = vlistNvars(vlistID2);
    cmorVars.resize(numVars);

    if (convertData)
      for (auto &var : cmorVars) var.convert = true;

    if (tableFormat == 0)
    {
      // for (int varID = 0; varID < numVars; ++varID) vlistDefVarTable(vlistID2, varID, tableID);
      apply_codetable(vlistID2, tableID);
    }
    else
    {
      {
        auto filename = cdo_operator_argv(0);
        auto fobj = c_fopen(filename, "r");
        if (fobj.get() == nullptr) cdo_abort("Open failed on: %s\n", filename);
        PMList pmlist;
        pmlist.read_namelist(fobj.get(), filename);

        apply_parameterList(ptmode, pmlist, vlistID2, cmorVars);
      }

      for (auto const &var : cmorVars)
        if (var.remove)
        {
          deleteVars = true;
          break;
        }

      if (deleteVars)
      {
        vlistClearFlag(vlistID1);
        vlistClearFlag(vlistID2);

        for (int varID = 0; varID < numVars; ++varID)
        {
          auto zaxisID = vlistInqVarZaxis(vlistID2, varID);
          auto numLevels = zaxisInqSize(zaxisID);
          for (int levelID = 0; levelID < numLevels; levelID++)
          {
            vlistDefFlag(vlistID1, varID, levelID, (not cmorVars[varID].remove));
            vlistDefFlag(vlistID2, varID, levelID, (not cmorVars[varID].remove));
          }
        }

        auto vlistIDx = vlistCreate();
        cdo_vlist_copy_flag(vlistIDx, vlistID2);

        vlistDestroy(vlistID2);

        vlistID2 = vlistIDx;
        if (vlistNvars(vlistID2) == 0) cdo_abort("No variable selected!");
      }

      for (auto &var : cmorVars)
      {
        if (!var.convert) var.changeUnits = false;
        if (var.changeUnits)
          cdo::convert_units(&var.ut_converter, &var.changeUnits, (char *) &var.units, (char *) &var.unitsOld, var.name);
      }
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    // vlistPrint(vlistID2);
    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    auto gridsizeMax = varList2.gridsizeMax();
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizeMax *= 2;
    Varray<double> array(gridsizeMax);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      cmor_check_init(cmorVars);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        auto &cmorVar = cmorVars[varID];
        auto varID2 = varID;
        auto levelID2 = levelID;

        if (deleteVars)
        {
          if (cmorVar.remove) continue;

          if (vlistInqFlag(vlistID1, varID, levelID) == true)
          {
            varID2 = vlistFindVar(vlistID2, varID);
            levelID2 = vlistFindLevel(vlistID2, varID, levelID);
          }
        }

        cdo_def_field(streamID2, varID2, levelID2);

        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);

        auto const &var2 = varList2.vars[varID2];
        auto missval = var2.missval;
        auto gridsize = var2.nwpv * var2.gridsize;

        if (numMissVals && cmorVar.changeMissval)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(array[i], cmorVar.missvalOld)) { array[i] = missval; }
          }
        }

        if (cmorVar.applyFactor)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_not_equal(array[i], missval)) { array[i] *= cmorVar.factor; }
          }
        }

#ifdef HAVE_UDUNITS2
        if (cmorVar.changeUnits)
        {
          int nerr = 0;
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_not_equal(array[i], missval))
            {
              array[i] = cv_convert_double((const cv_converter *) cmorVar.ut_converter, array[i]);
              if (ut_get_status() != UT_SUCCESS) nerr++;
            }
          }
          if (nerr)
          {
            cdo_warning("Udunits: Error converting units from [%s] to [%s], parameter: %s", cmorVar.unitsOld, cmorVar.units,
                        cmorVar.name);
            cmorVar.changeUnits = false;
          }
        }
#endif

        cdo_write_field(streamID2, array.data(), numMissVals);
        cmor_check_prep(cmorVar, gridsize, missval, array.data());
      }

      cmor_check_eval(vlistID2, cmorVars);

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

#ifdef HAVE_UDUNITS2
    for (auto &var : cmorVars)
      if (var.changeUnits) cdo::convert_free(var.ut_converter);

    cdo::convert_destroy();
#endif
  }
};
