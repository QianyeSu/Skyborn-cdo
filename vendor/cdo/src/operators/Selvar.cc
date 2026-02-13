/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selvar     selparam        Select parameters by identifier (format: code.tabnum or pnum.cat.dis)
      Selvar     delparam        Delete parameters by identifier (format: code.tabnum or pnum.cat.dis)
      Selvar     selcode         Select parameters by code number
      Selvar     delcode         Delete parameters by code number
      Selvar     selname         Select parameters by name
      Selvar     delname         Delete parameters by name
      Selvar     selstdname      Select parameters by CF standard name
      Selvar     sellevel        Select levels
      Selvar     sellevidx       Select levels by index
      Selvar     selgrid         Select grids
      Selvar     selzaxis        Select zaxis
      Selvar     seltabnum       Select parameter table number
      Selvar     selltype        Select GRIB level type
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "util_wildcards.h"
#include "cdi_lockedIO.h"
#include "param_conversion.h"

class Selvar : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selvar",
    // clang-format off
    .operators = { { "selparam",     0, 2, "parameters", SelvarHelp },
                   { "selcode",      0, 4, "code numbers", SelvarHelp },
                   { "selname",      0, 2, "variable names", SelvarHelp },
                   { "selstdname",   0, 2, "standard names", SelvarHelp },
                   { "sellevel",     0, 8, "levels", SelvarHelp },
                   { "sellevidx",    0, 4, "index of levels", SelvarHelp },
                   { "selgrid",      0, 4 | 2, "list of grid names or numbers", SelvarHelp },
                   { "selzaxis",     0, 4 | 2, "list of zaxis types or numbers", SelvarHelp },
                   { "selzaxisname", 0, 2, "list of zaxis names", SelvarHelp },
                   { "seltabnum",    0, 4, "table numbers", SelvarHelp },
                   { "delparam",     1, 2 | 1, "parameter", SelvarHelp },
                   { "delcode",      1, 1, "code numbers", SelvarHelp },
                   { "delname",      1, 2 | 1, "variable names", SelvarHelp },
                   { "selltype",     0, 4, "GRIB level types", SelvarHelp } },
    // clang-format on
    .aliases = { { "selvar", "selname" }, { "delvar", "delname" }, { "selgridname", "selgrid" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Selvar> registration = RegisterEntry<Selvar>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

  bool dataIsUnchanged{};

public:
  void
  init() override
  {
    int nsel = 0;
    char paramstr[32];
    char gridname[CDI_MAX_NAME];
    char zaxistypename[CDI_MAX_NAME];
    std::vector<int> intarr;
    std::vector<double> fltarr;

    dataIsUnchanged = data_is_unchanged();

#define INVERTS_SELECTION(id) (cdo_operator_f2(id) & 1)
#define TAKES_STRINGS(id) (cdo_operator_f2(id) & 2)
#define TAKES_INTEGERS(id) (cdo_operator_f2(id) & 4)
#define TAKES_FLOATS(id) (cdo_operator_f2(id) & 8)

    auto SELPARAM = module.get_id("selparam");
    auto SELCODE = module.get_id("selcode");
    auto SELNAME = module.get_id("selname");
    auto SELSTDNAME = module.get_id("selstdname");
    auto SELLEVEL = module.get_id("sellevel");
    auto SELLEVIDX = module.get_id("sellevidx");
    auto SELGRID = module.get_id("selgrid");
    auto SELZAXIS = module.get_id("selzaxis");
    auto SELZAXISNAME = module.get_id("selzaxisname");
    auto SELTABNUM = module.get_id("seltabnum");
    auto DELPARAM = module.get_id("delparam");
    auto DELCODE = module.get_id("delcode");
    auto DELNAME = module.get_id("delname");
    auto SELLTYPE = module.get_id("selltype");

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto doDelete = (cdo_operator_f1(operatorID) == 1);

    auto argsAreNumeric = (cdo_operator_argc() > 0 && std::isdigit(cdo_operator_argv(0)[0]));

    auto const &argList = cdo_get_oper_argv();
    if (TAKES_STRINGS(operatorID) && !(TAKES_INTEGERS(operatorID) && argsAreNumeric))
    {
      nsel = cdo_operator_argc();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("name %d = %s", i + 1, argList[i]);
    }
    else if (TAKES_FLOATS(operatorID))
    {
      fltarr = cdo_argv_to_fltarr(argList);
      nsel = fltarr.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("flt %d = %g", i + 1, fltarr[i]);
    }
    else
    {
      intarr = cdo_argv_to_intarr(argList);
      nsel = intarr.size();

      if (Options::cdoVerbose)
        for (int i = 0; i < nsel; ++i) cdo_print("int %d = %d", i + 1, intarr[i]);
    }

    std::vector<bool> selfound;
    if (nsel)
    {
      selfound.resize(nsel);
      for (int i = 0; i < nsel; ++i) selfound[i] = false;
    }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);
    auto numVars = varList1.numVars();
    std::vector<bool> vars(numVars);

    if (operatorID == SELGRID && !argsAreNumeric && nsel == 1 && argList[0].starts_with("var="))
    {
      int gridnum = 0;
      const char *gridvarname = &argList[0][4];
      if (*gridvarname == 0) cdo_abort("Variable name missing!");

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (var.name == gridvarname)
        {
          gridnum = 1 + vlistGridIndex(vlistID1, var.gridID);
          argsAreNumeric = true;
          intarr.push_back(gridnum);
          break;
        }
      }

      if (!gridnum) cdo_abort("Variable %s not found!", gridvarname);
    }

    vlistClearFlag(vlistID1);
    for (int varID = 0; varID < numVars; ++varID)
    {
      vars[varID] = doDelete;

      auto const &var = varList1.vars[varID];
      auto tabnum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      auto grididx = vlistGridIndex(vlistID1, var.gridID);
      auto zaxisidx = vlistZaxisIndex(vlistID1, var.zaxisID);
      auto numLevels = zaxisInqSize(var.zaxisID);
      gridName(gridInqType(var.gridID), gridname);
      auto zaxisname = cdo::inq_key_string(var.zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
      zaxisName(zaxisInqType(var.zaxisID), zaxistypename);

      cdiParamToString(var.param, paramstr, sizeof(paramstr));

      std::string gridnameStr{ gridname };
      std::string zaxisnameStr{ zaxistypename };
      for (int levelID = 0; levelID < numLevels; levelID++)
      {
        auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);

        if (doDelete) vlistDefFlag(vlistID1, varID, levelID, true);

        for (int isel = 0; isel < nsel; isel++)
        {
          auto found = false;
          // clang-format off
          if      (operatorID == SELCODE) found = (intarr[isel] == var.code);
          else if (operatorID == SELPARAM) found = wildcard_match(paramstr, argList[isel]);
          else if (operatorID == SELNAME) found = wildcard_match(var.name, argList[isel]);
          else if (operatorID == SELSTDNAME) found = wildcard_match(var.stdname, argList[isel]);
          else if (operatorID == SELLEVEL) found = (std::fabs(fltarr[isel] - level) < 0.0001);
          else if (operatorID == SELLEVIDX) found = (intarr[isel] == (levelID + 1));
          else if (operatorID == SELGRID && argsAreNumeric) found = (intarr[isel] == (grididx + 1));
          else if (operatorID == SELGRID && !argsAreNumeric) found = gridnameStr.starts_with(argList[isel]);
          else if (operatorID == SELZAXIS && argsAreNumeric) found = (intarr[isel] == (zaxisidx + 1));
          else if (operatorID == SELZAXIS && !argsAreNumeric) found = zaxisnameStr.starts_with(argList[isel]);
          else if (operatorID == SELZAXISNAME) found = wildcard_match(zaxisname, argList[isel]);
          else if (operatorID == SELTABNUM) found = (intarr[isel] == tabnum);
          else if (operatorID == DELCODE) found = (intarr[isel] == var.code);
          else if (operatorID == DELNAME) found = wildcard_match(var.name, argList[isel]);
          else if (operatorID == DELPARAM) found = (argList[isel] == paramstr);
          else if (operatorID == SELLTYPE) found = (intarr[isel] == zaxis_to_ltype(var.zaxisID));
          // clang-format on

          if (found)
          {
            vlistDefFlag(vlistID1, varID, levelID, !INVERTS_SELECTION(operatorID));
            selfound[isel] = true;
            vars[varID] = !doDelete;
          }
        }
      }
    }

    int npar = 0;
    for (int varID = 0; varID < numVars; ++varID)
      if (vars[varID]) npar++;

    for (int varID = 0; varID < numVars; ++varID)
    {
      if (vars[varID])
      {
        auto const &var = varList1.vars[varID];
        if (zaxisInqType(var.zaxisID) == ZAXIS_HYBRID)
        {
          auto psvarid = varList_get_psvarid(varList1, var.zaxisID);
          if (psvarid != -1 && !vars[psvarid])
          {
            vars[psvarid] = true;
            vlistDefFlag(vlistID1, psvarid, 0, !INVERTS_SELECTION(operatorID));
          }
        }
      }
    }

    for (int isel = 0; isel < nsel; isel++)
    {
      if (selfound[isel] == false)
      {
        // clang-format off
        if      (operatorID == SELCODE || operatorID == DELCODE) cdo_warning("Code number %d not found!", intarr[isel]);
        else if (operatorID == SELPARAM || operatorID == DELPARAM) cdo_warning("Parameter %s not found!", argList[isel]);
        else if (operatorID == SELNAME || operatorID == DELNAME) cdo_warning("Variable name %s not found!", argList[isel]);
        else if (operatorID == SELSTDNAME) cdo_warning("Variable with standard name %s not found!", argList[isel]);
        else if (operatorID == SELLEVEL) cdo_warning("Level %g not found!", fltarr[isel]);
        else if (operatorID == SELLEVIDX) cdo_warning("Level index %d not found!", intarr[isel]);
        else if (operatorID == SELGRID && argsAreNumeric) cdo_warning("Grid %d not found!", intarr[isel]);
        else if (operatorID == SELGRID && !argsAreNumeric) cdo_warning("Grid name %s not found!", argList[isel]);
        else if (operatorID == SELZAXIS && argsAreNumeric) cdo_warning("Zaxis %d not found!", intarr[isel]);
        else if (operatorID == SELZAXIS && !argsAreNumeric) cdo_warning("Zaxis type %s not found!", argList[isel]);
        else if (operatorID == SELZAXISNAME) cdo_warning("Zaxis name %s not found!", argList[isel]);
        else if (operatorID == SELTABNUM) cdo_warning("Table number %d not found!", intarr[isel]);
        else if (operatorID == SELLTYPE) cdo_warning("GRIB level type %d not found!", intarr[isel]);
        // clang-format on
      }
    }

    if (npar == 0) cdo_abort("No variables selected!");

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);

    // if (Options::cdoVerbose) vlistPrint(vlistID2);

    numVars = vlistNvars(vlistID2);
    {
      int varID;
      for (varID = 0; varID < numVars; ++varID)
        if (vlistInqVarTimetype(vlistID2, varID) != TIME_CONSTANT) break;
      if (varID == numVars) vlistDefNtsteps(vlistID2, 0);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

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
        if (vlistInqFlag(vlistID1, varID, levelID) == true)
        {
          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);
          cdo_def_field(streamID2, varID2, levelID2);

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

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
