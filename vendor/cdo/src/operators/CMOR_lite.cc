/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      CMORlite      cmorlite        CMOR lite
*/

#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_cmor.h"
#include "pmlist.h"
#include "convert_units.h"
#include "cdi_lockedIO.h"

void
cdo_define_var_units(CmorVar &cmorVar, int vlistID2, int varID, std::string const &units)
{
  auto unitsOld = cdo::inq_var_units(vlistID2, varID);
  if (units != unitsOld)
  {
    if (unitsOld.size() > 0 && units.size() > 0)
    {
      cmorVar.changeUnits = true;
      std::strcpy(cmorVar.unitsOld, unitsOld.c_str());
      std::strcpy(cmorVar.units, units.c_str());
    }

    cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, units.c_str());
    cdiDefAttTxt(vlistID2, varID, "original_units", (int) unitsOld.size(), unitsOld.c_str());
  }
}

void
cmor_check_init(std::vector<CmorVar> &cmorVars)
{
  for (auto &var : cmorVars)
  {
    if (var.checkValid || var.check_min_mean_abs || var.check_max_mean_abs)
    {
      var.amean = 0;
      var.nvals = 0;
      var.n_lower_min = 0;
      var.n_greater_max = 0;
    }
  }
}

void
cmor_check_eval(int vlistID, std::vector<CmorVar> const &cmorVars)
{
  int numVars = cmorVars.size();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &cmorVar = cmorVars[varID];
    if (cmorVar.checkValid || cmorVar.check_min_mean_abs || cmorVar.check_max_mean_abs)
    {
      auto amean = cmorVar.amean;
      auto nvals = cmorVar.nvals;

      if (nvals > 0) amean /= nvals;

      auto n_lower_min = cmorVar.n_lower_min;
      auto n_greater_max = cmorVar.n_greater_max;

      auto varname = cdo::inq_var_name(vlistID, varID);

      if (n_lower_min > 0)
        cdo_warning("Invalid value(s) detected for variable '%s': %ld values were lower than minimum valid value (%.4g).", varname,
                    n_lower_min, cmorVar.valid_min);
      if (n_greater_max > 0)
        cdo_warning("Invalid value(s) detected for variable '%s': %ld values were greater than maximum valid value (%.4g).",
                    varname, n_greater_max, cmorVar.valid_max);

      if (cmorVar.check_min_mean_abs)
      {
        if (amean < .1 * cmorVar.ok_min_mean_abs)
          cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is lower by more than an order of magnitude than "
                      "minimum allowed: %.4g",
                      varname, amean, cmorVar.ok_min_mean_abs);

        if (amean < cmorVar.ok_min_mean_abs)
          cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is lower than minimum allowed: %.4g", varname, amean,
                      cmorVar.ok_min_mean_abs);
      }

      if (cmorVar.check_max_mean_abs)
      {
        if (amean > 10. * cmorVar.ok_max_mean_abs)
          cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is greater by more than an order of magnitude than "
                      "maximum allowed: %.4g",
                      varname, amean, cmorVar.ok_max_mean_abs);

        if (amean > cmorVar.ok_max_mean_abs)
          cdo_warning("Invalid Absolute Mean for variable '%s' (%.5g) is greater than maximum allowed: %.4g", varname, amean,
                      cmorVar.ok_max_mean_abs);
      }
    }
  }
}

void
cmor_check_prep(CmorVar &var, long gridsize, double missval, const double *const array)
{
  if (var.checkValid || var.check_min_mean_abs || var.check_max_mean_abs)
  {
    double amean = 0;
    long nvals = 0;

    for (long i = 0; i < gridsize; ++i)
    {
      auto aval = array[i];
      if (fp_is_not_equal(aval, missval))
      {
        amean += std::fabs(aval);
        nvals++;
      }
    }

    var.amean += amean;
    var.nvals += nvals;

    long n_lower_min = 0, n_greater_max = 0;

    for (long i = 0; i < gridsize; ++i)
    {
      auto aval = array[i];
      if (fp_is_not_equal(aval, missval))
      {
        if (aval < var.valid_min) n_lower_min++;
        if (aval > var.valid_max) n_greater_max++;
      }
    }

    var.n_lower_min += n_lower_min;
    var.n_greater_max += n_greater_max;
  }
}

static void
search_global_missval(PMList &pmlist, int vlistID2, bool &hasMissvals, double &missval)
{
  static const std::vector<std::string> hentry = { "Header" };
  auto kvlist = pmlist.getKVListVentry(hentry);
  if (kvlist)
  {
    for (auto const &kv : *kvlist)
    {
      auto const &key = kv.key;
      auto const &value = kv.values[0];
      if (kv.nvalues != 1 || value.empty()) continue;

      if (key == "missing_value")
      {
        hasMissvals = true;
        missval = parameter_to_double(value);
      }
      else if (key == "table_id" || key == "modeling_realm" || key == "realm" || key == "project_id" || key == "frequency")
      {
        cdiDefAttTxt(vlistID2, CDI_GLOBAL, key.c_str(), (int) value.size(), value.c_str());
      }
    }
  }
}

static void
apply_cmor_list(PMList &pmlist, int vlistID2, std::vector<CmorVar> &cmorVars)
{
  // search for global missing value
  auto hasMissvals = false;
  double missval = 0.0;
  search_global_missval(pmlist, vlistID2, hasMissvals, missval);

  static const std::vector<std::string> ventry = { "variable_entry", "parameter" };

  int numVars = cmorVars.size();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto &cmorVar = cmorVars[varID];
    cmorVar.name = cdo::inq_var_name(vlistID2, varID);

    if (hasMissvals)
    {
      auto missvalOld = vlistInqVarMissval(vlistID2, varID);
      if (fp_is_not_equal(missval, missvalOld))
      {
        cmorVar.changeMissval = true;
        cmorVar.missvalOld = missvalOld;
        vlistDefVarMissval(vlistID2, varID, missval);
      }
    }

    auto kvlist = pmlist.searchKVListVentry("name", cmorVar.name, ventry);
    if (kvlist)
    {
      auto hasValidMin = false, hasValidMax = false;

      for (auto const &kv : *kvlist)
      {
        auto const &key = kv.key;
        auto const &value = kv.values[0];
        if (kv.nvalues != 1 || value.empty()) continue;
        auto value_cstr = value.c_str();

        // printf("key=%s  value=>%s<\n", key.c_str(), value.c_str());

        // clang-format off
        if      (key == "standard_name") cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, value_cstr);
        else if (key == "long_name")     cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, value_cstr);
        else if (key == "units")         cdo_define_var_units(cmorVar, vlistID2, varID, value);
        else if (key == "name") ;     // cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, parameter_to_word(value));
        else if (key == "out_name")
        {
          auto outname = parameter_to_word(value);
          if (cmorVar.name != outname)
          {
            cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, outname.c_str());
            cdiDefAttTxt(vlistID2, varID, "original_name", cmorVar.name.size(), cmorVar.name.c_str());
          }
        }
        else if (key == "param")         vlistDefVarParam(vlistID2, varID, string_to_param(parameter_to_word(value)));
        else if (key == "out_param")     vlistDefVarParam(vlistID2, varID, string_to_param(parameter_to_word(value)));
        else if (key == "comment")       cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
        else if (key == "cell_methods")  cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
        else if (key == "cell_measures") cdiDefAttTxt(vlistID2, varID, key.c_str(), (int) value.size(), value_cstr);
        else if (key == "delete")        cmorVar.remove = parameter_to_bool(value);
        else if (key == "convert")       cmorVar.convert = parameter_to_bool(value);
        else if (key == "factor")
        {
          cmorVar.applyFactor = true;
          cmorVar.factor = parameter_to_double(value);
          if (Options::cdoVerbose) cdo_print("%s - scale factor %g", cmorVar.name, cmorVar.factor);
        }
        else if (key == "missval" || key == "missing_value")
        {
          missval = parameter_to_double(value);
          auto missvalOld = vlistInqVarMissval(vlistID2, varID);
          if (fp_is_not_equal(missval, missvalOld))
          {
            if (Options::cdoVerbose) cdo_print("%s - change missval from %g to %g", cmorVar.name, missvalOld, missval);
            cmorVar.changeMissval = true;
            cmorVar.missvalOld = missvalOld;
            vlistDefVarMissval(vlistID2, varID, missval);
          }
        }
        else if (key == "valid_min")
        {
          hasValidMin = true;
          cmorVar.valid_min = parameter_to_double(value);
        }
        else if (key == "valid_max")
        {
          hasValidMax = true;
          cmorVar.valid_max = parameter_to_double(value);
        }
        else if (key == "ok_min_mean_abs")
        {
          cmorVar.check_min_mean_abs = true;
          cmorVar.ok_min_mean_abs = parameter_to_double(value);
        }
        else if (key == "ok_max_mean_abs")
        {
          cmorVar.check_max_mean_abs = true;
          cmorVar.ok_max_mean_abs = parameter_to_double(value);
        }
        else if (key == "datatype" || key == "type")
        {
          auto datatype = cdo::str_to_datatype(parameter_to_word(value));
          if (datatype != -1) vlistDefVarDatatype(vlistID2, varID, datatype);
        }
        else
        {
          if (Options::cdoVerbose) cdo_print("Attribute %s:%s not supported!", cmorVar.name, key);
        }
        // clang-format on
      }

      if (hasValidMin && hasValidMax) cmorVar.checkValid = true;
    }
    else { cdo_print("Variable %s not found in CMOR table!", cmorVar.name); }
  }
}

class CMOR_lite : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CMOR_lite",
    .operators = { { "cmorlite", 0, 0, "parameter table name", CmorliteHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<CMOR_lite>();

  bool deleteVars = false;

  CdoStreamID streamID1;
  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int vlistID2{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  std::vector<CmorVar> cmorVars;
  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    Options::CMOR_Mode = 1;
    if (Options::CMOR_Mode) cdiDefGlobal("CMOR_MODE", Options::CMOR_Mode);

    operator_input_arg(cdo_operator_enter(cdo_operator_id()));

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

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    auto numVars = varList1.numVars();
    cmorVars.resize(numVars);

    if (convertData)
      for (auto &var : cmorVars) var.convert = true;

    PMList pmlist;
    {
      auto filename = cdo_operator_argv(0);
      auto fobj = c_fopen(filename, "r");
      if (fobj.get() == nullptr) cdo_abort("Open failed on: %s\n", filename);
      pmlist.read_cmor_table(fobj.get(), filename);
    }

    apply_cmor_list(pmlist, vlistID2, cmorVars);

    varList2 = VarList(vlistID2);

    for (auto &var : cmorVars)
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
        for (int levelID = 0; levelID < varList2.vars[varID].nlevels; levelID++)
        {
          vlistDefFlag(vlistID1, varID, levelID, true);
          vlistDefFlag(vlistID2, varID, levelID, true);
          if (cmorVars[varID].remove)
          {
            vlistDefFlag(vlistID1, varID, levelID, false);
            vlistDefFlag(vlistID2, varID, levelID, false);
          }
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

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    // vlistPrint(vlistID2);
    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
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
        auto gridsize = var2.gridsize;
        if (var2.nwpv != CDI_REAL) gridsize *= 2;

        if (numMissVals && cmorVar.changeMissval)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(array[i], cmorVar.missvalOld)) array[i] = missval;
          }
        }

        if (cmorVar.applyFactor)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_not_equal(array[i], missval)) array[i] *= cmorVar.factor;
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
