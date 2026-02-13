/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Expr      expr            Evaluate expressions
      Expr      exprf           Evaluate expressions from script file
      Expr      aexpr           Append evaluated expressions
      Expr      aexprf          Append evaluated expressions from script file
*/
/*
  Operatoren: +, -, *, \, ^, ==, !=, >, <, >=, <=, <=>, &&, ||, ?:
  Functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
  Functions: min, max, avg, std, var
  Constansts: M_PI, M_E
*/

#include <fstream>
#include <algorithm>
#include <cassert>

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "expr.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"
#include "util_string.h"
#include "datetime.h"

void gridcell_areas(int gridID, Varray<double> &array);
int get_surface_ID(int vlistID);  // from Vertstat.cc
struct yy_buffer_state *yy_scan_string(const char *str, void *scanner);

static std::string
exprs_from_argument(std::vector<std::string> const &exprArgv)
{
  std::string exprString{};

  if (exprArgv.size() > 0)
  {
    for (size_t i = 0, n = exprArgv.size(); i < n; ++i)
    {
      if (i > 0) exprString += ",";
      exprString += exprArgv[i];
    }
    if (exprString[exprString.size() - 1] != ';') exprString += ";";
  }
  else { operator_check_argc(1); }

  return exprString;
}

static std::string
exprs_from_file(std::vector<std::string> const &exprArgv)
{
  if (exprArgv.size() != 1) operator_check_argc(1);
  auto exprFile = exprArgv[0];
  std::ifstream stream(exprFile);
  if (!stream.is_open()) cdo_abort("Open failed on %s", exprFile);
  std::stringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
}

constexpr int MaxParams = 4096;

static std::size_t
replace_all(std::string &inout, std::string const &what, std::string const &with)
{
  std::size_t count{};
  for (std::string::size_type pos{}; inout.npos != (pos = inout.find(what.data(), pos, what.length()));
       pos += with.length(), ++count)
  {
    inout.replace(pos, what.length(), with.data(), with.length());
  }
  return count;
}

static std::string
exprs_expand(std::string const &exprString, VarList const &varList)
{
  auto replaceTemplate = false;
  std::string templateName = "_ALL_";

  if (exprString.find(templateName) != std::string::npos)
  {
    replaceTemplate = true;
    for (auto const &var : varList.vars)
    {
      if (templateName == var.name)
      {
        replaceTemplate = false;
        break;
      }
    }
  }

  if (replaceTemplate)
  {
    std::string exprStringNew{};
    auto exprStringArgv = split_string(exprString, ";");
    for (auto const &string : exprStringArgv)
    {
      if (string.find(templateName) == std::string::npos) { exprStringNew += string + ";"; }
      else
      {
        for (auto const &var : varList.vars)
        {
          auto tmpString = string;
          replace_all(tmpString, templateName, var.name);
          exprStringNew += tmpString + ";";
        }
      }
    }

    return exprStringNew;
  }

  return exprString;
}

static void
params_init(std::vector<ParamEntry> &params, VarList const &varList)
{
  for (auto const &var : varList.vars)
  {
    auto &param = params[var.ID];

    param.type = ParamType::VAR;
    param.isValid = true;
    param.hasMV = true;
    param.gridID = var.gridID;
    param.zaxisID = var.zaxisID;
    param.datatype = var.dataType;
    param.steptype = var.timeType;
    param.nlat = gridInqYsize(var.gridID);
    param.ngp = var.gridsize;
    param.nlev = var.nlevels;
    param.missval = var.missval;
    param.name = var.name;
    if (var.longname.size()) param.longname = var.longname;
    if (var.units.size()) param.units = var.units;
    if (var.stdname.size()) param.stdname = var.stdname;
  }
}

static void
params_delete(std::vector<ParamEntry> const &params)
{
  for (auto const &param : params)
  {
    if (param.data) delete[] param.data;
  }
}

static void
params_add_coord(ParseParamType &parseArg, int coord, int cdiID, size_t size, std::string const &units, std::string const &longname)
{
  auto ncoords = parseArg.numCoords;
  if (ncoords >= parseArg.maxCoords) cdo_abort("Too many coordinates (limit=%d)", parseArg.maxCoords);

  auto &coords = parseArg.coords[ncoords];
  coords.needed = false;
  coords.coord = coord;
  coords.cdiID = cdiID;
  coords.size = size;
  if (units.size()) coords.units = units;
  if (longname.size()) coords.longname = longname;

  parseArg.numCoords++;
}

static void
params_add_coordinates(int vlistID, ParseParamType &parseArg)
{
  auto numGrids = vlistNumGrids(vlistID);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(vlistID, index);
    auto size = gridInqSize(gridID);
    auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
    auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);

    params_add_coord(parseArg, 'x', gridID, size, xunits, "longitude");
    params_add_coord(parseArg, 'y', gridID, size, yunits, "latitude");
    params_add_coord(parseArg, 'a', gridID, size, "m^2", "grid cell area");
    params_add_coord(parseArg, 'w', gridID, size, "", "grid cell area weights");
    params_add_coord(parseArg, 'g', gridID, size, "", "grid cell indices");
  }

  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    auto size = zaxisInqSize(zaxisID);
    auto zunits = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
    auto longname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME);
    params_add_coord(parseArg, 'z', zaxisID, size, zunits, longname);
    params_add_coord(parseArg, 'i', zaxisID, size, zunits, "level index");
    params_add_coord(parseArg, 'd', zaxisID, size, zunits, "delta z");
  }
}

static int
params_add_ts(ParseParamType &parseArg)
{
  auto &params = parseArg.params;

  auto varID = parseArg.numParams;
  if (varID >= parseArg.maxParams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxParams);

  auto &param = params[varID];
  param.name = "_timestep_info";
  param.gridID = parseArg.pointID;
  param.zaxisID = parseArg.surfaceID;
  param.steptype = TIME_VARYING;
  param.ngp = CoordIndex::LEN;
  param.nlev = 1;

  parseArg.numParams++;

  return varID;
}

static void
parse_param_init(ParseParamType &parseArg, int vlistID, int pointID, int zonalID, int surfaceID)
{
  auto numVars = vlistNvars(vlistID);
  auto numGrids = vlistNumGrids(vlistID);
  auto numZaxes = vlistNumZaxis(vlistID);
  auto maxCoords = numGrids * 5 + numZaxes * 3;

  parseArg.maxParams = MaxParams;
  parseArg.params.resize(MaxParams);
  parseArg.numParams = numVars;
  parseArg.numVars1 = numVars;
  parseArg.init = true;
  parseArg.debug = (Options::cdoVerbose != 0);
  parseArg.pointID = pointID;
  parseArg.zonalID = zonalID;
  parseArg.surfaceID = surfaceID;
  parseArg.needed.resize(numVars);
  parseArg.coords.resize(maxCoords);
  parseArg.maxCoords = maxCoords;
  parseArg.numCoords = 0;
}

static int
genZonalID(int vlistID)
{
  int zonalID = -1;

  auto numGrids = vlistNumGrids(vlistID);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(vlistID, index);
    auto gridtype = gridInqType(gridID);
    if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC)
      if (gridInqXsize(gridID) > 1 && gridInqYsize(gridID) >= 1)
      {
        zonalID = gridToZonal(gridID);
        break;
      }
  }

  return zonalID;
}

static void
set_date_and_time(ParamEntry &varts, int calendar, int tsID, CdiDateTime const &vDateTime0, CdiDateTime const &vDateTime)
{
  double jdelta = 0.0;

  if (tsID)
  {
    auto julianDate0 = julianDate_encode(calendar, vDateTime0);
    auto julianDate = julianDate_encode(calendar, vDateTime);
    jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
  }

  varts.data[CoordIndex::TIMESTEP] = tsID + 1;
  varts.data[CoordIndex::DATE] = cdiDate_get(vDateTime.date);
  varts.data[CoordIndex::TIME] = cdiTime_get(vDateTime.time);
  varts.data[CoordIndex::DELTAT] = jdelta;

  int year, mon, day;
  int hour, minute, second, ms;
  cdiDate_decode(vDateTime.date, &year, &mon, &day);
  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);

  varts.data[CoordIndex::DAY] = day;
  varts.data[CoordIndex::MONTH] = mon;
  varts.data[CoordIndex::YEAR] = year;
  varts.data[CoordIndex::SECOND] = second;
  varts.data[CoordIndex::MINUTE] = minute;
  varts.data[CoordIndex::HOUR] = hour;

  varts.data[CoordIndex::CALENDAR] = calendar;
  varts.data[CoordIndex::DOY] = day_of_year(calendar, cdiDate_get(vDateTime.date));
  varts.data[CoordIndex::DPY] = days_per_year(calendar, year);
}

class Expr : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Expr",
    .operators = { { "expr", 1, 1, "expressions", ExprHelp },
                   { "exprf", 1, 0, "exprscriptfilename", ExprHelp },
                   { "aexpr", 0, 1, "expressions", ExprHelp },
                   { "aexprf", 0, 0, "exprscriptfilename", ExprHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Expr> registration = RegisterEntry<Expr>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vartsID{};

  int numVars1{};
  int numVars2{};

  ParseParamType parseArg{};
  CdiDateTime vDateTime0{};

  std::vector<int> varIDmap{};
  std::string exprString{};

  int pointID{};
  int zonalID{};
  int surfaceID{};
  int calendar{};

  void
  parse_expressions()
  {
    void *scanner = nullptr;
    yylex_init(&scanner);
    yyset_extra(&parseArg, scanner);
    yy_scan_string(exprString.c_str(), scanner);
    yyparse(parseArg, scanner);
    yylex_destroy(scanner);
  }

  void
  allocate_params()
  {
    auto &params = parseArg.params;
    for (int varID = 0; varID < numVars1; ++varID)
    {
      if (parseArg.needed[varID])
      {
        auto nItems = std::max((size_t) 4, params[varID].ngp * params[varID].nlev);
        params[varID].data = new double[nItems];
      }
    }

    for (int varID = parseArg.numVars1, n = parseArg.numParams; varID < n; ++varID)
    {
      auto nItems = std::max((size_t) 4, params[varID].ngp * params[varID].nlev);
      params[varID].data = new double[nItems];
    }
  }

  void
  read_coordinates_hgrid(Varray<double> &cdata, size_t csize, int cdiID, int coord)
  {
    {
      auto gridID = cdiID;
      auto ngp = csize;
      cdata.resize(ngp);
      if (coord == 'x' || coord == 'y')
      {
        gridID = generate_full_point_grid(gridID);
        if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

        if (coord == 'x') gridInqXvals(gridID, cdata.data());
        if (coord == 'y') gridInqYvals(gridID, cdata.data());

        if (gridID != cdiID) gridDestroy(gridID);
      }
      else if (coord == 'a') { gridcell_areas(gridID, cdata); }
      else if (coord == 'w')
      {
        cdata[0] = 1;
        if (ngp > 1)
        {
          auto wstatus = gridcell_weights(gridID, cdata);
          if (wstatus) cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
        }
      }
      else if (coord == 'g')
      {
        for (size_t k = 0; k < ngp; ++k) cdata[k] = k + 1;
      }
    }
  }

  void
  read_coordinates_vgrid(Varray<double> &cdata, size_t csize, int cdiID, int coord)
  {
    {
      auto zaxisID = cdiID;
      auto nlev = csize;
      cdata.resize(nlev);
      if (coord == 'z') { cdo_zaxis_inq_levels(zaxisID, cdata.data()); }
      else if (coord == 'i')
      {
        for (size_t k = 0; k < nlev; ++k) cdata[k] = k + 1;
        cdo_zaxis_inq_levels(zaxisID, cdata.data());
      }
      else if (coord == 'd')
      {
        std::ranges::fill(cdata, 1.0);
        if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
        {
          std::vector<double> lbounds(nlev), ubounds(nlev);
          zaxisInqLbounds(zaxisID, lbounds.data());
          zaxisInqUbounds(zaxisID, ubounds.data());
          for (size_t k = 0; k < nlev; ++k) cdata[k] = ubounds[k] - lbounds[k];
        }
      }
    }
  }

  void
  read_coordinates()
  {
    for (int i = 0, n = parseArg.numCoords; i < n; ++i)
    {
      if (parseArg.coords[i].needed)
      {
        auto &cdata = parseArg.coords[i].data;
        auto csize = parseArg.coords[i].size;
        auto cdiID = parseArg.coords[i].cdiID;
        auto coord = parseArg.coords[i].coord;
        if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w' || coord == 'g')
        {
          read_coordinates_hgrid(cdata, csize, cdiID, coord);
        }
        else if (coord == 'z' || coord == 'i' || coord == 'd') { read_coordinates_vgrid(cdata, csize, cdiID, coord); }
        else { cdo_abort("Computation of coordinate %c not implemented!", coord); }
      }
    }
  }

  void
  copy_coordinates()
  {
    auto &params = parseArg.params;
    for (int varID = parseArg.numVars1, n = parseArg.numParams; varID < n; ++varID)
    {
      auto coord = params[varID].coord;
      if (coord)
      {
        if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w' || coord == 'g')
        {
          auto coordID = params_get_coord_ID(parseArg, coord, params[varID].gridID);
          auto gridID = parseArg.coords[coordID].cdiID;
          auto ngp = parseArg.coords[coordID].size;
          auto const &cdata = parseArg.coords[coordID].data;
          assert(gridID == params[varID].gridID);
          assert(!cdata.empty());

          array_copy(ngp, cdata.data(), params[varID].data);
        }
        else if (coord == 'z' || coord == 'i' || coord == 'd')
        {
          auto coordID = params_get_coord_ID(parseArg, coord, params[varID].zaxisID);
          auto zaxisID = parseArg.coords[coordID].cdiID;
          auto nlev = parseArg.coords[coordID].size;
          auto const &cdata = parseArg.coords[coordID].data;
          assert(zaxisID == params[varID].zaxisID);
          assert(!cdata.empty());

          array_copy(nlev, cdata.data(), params[varID].data);
        }
        else
          cdo_abort("Computation of coordinate %c not implemented!", coord);
      }
    }
  }

  void
  init_output_params(std::vector<ParamEntry> &params)
  {
    for (int varID = 0; varID < numVars2; ++varID)
    {
      auto pidx = varIDmap[varID];
      if (pidx < numVars1) continue;

      auto &param = params[pidx];
      param.numMissVals = 0;
      std::ranges::fill_n(param.data, param.ngp * param.nlev, 0.0);
    }
  }

  void
  add_new_output_params()
  {
    auto const &params = parseArg.params;
    //  printf("parseArg.nparams %d\n", parseArg.nparams);
    for (int pidx = 0, n = parseArg.numParams; pidx < n; pidx++)
    {
      auto const &param = params[pidx];
      if (pidx < numVars1 && !param.select) continue;
      if (pidx >= numVars1)
      {
        if (param.type == ParamType::CONST) continue;
        if (param.name[0] == '_') continue;
        if (param.remove) continue;
        if (param.coord) continue;
      }

      // printf("gridID %d zaxisID %d\n",  param.gridID, param.zaxisID);
      auto varID = vlistDefVar(vlistID2, param.gridID, param.zaxisID, param.steptype);
      cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, param.name.c_str());
      // printf("add: %d %s %d levs %d\n", pidx,  param.name.c_str(), varID, zaxisInqSize(param.zaxisID));
      if (param.hasMV) vlistDefVarMissval(vlistID2, varID, param.missval);
      if (param.units.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, param.units.c_str());
      if (param.longname.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, param.longname.c_str());
      if (param.stdname.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, param.stdname.c_str());
      if (param.name.size() > 3 && param.name.rfind("var", 0) == 0)
      {
        if (std::isdigit(param.name[3]))
        {
          auto code = std::atoi(param.name.c_str() + 3);
          vlistDefVarCode(vlistID2, varID, code);
        }
      }
      varIDmap[varID] = pidx;
    }
  }

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    bool replacesVariables = cdo_operator_f1(operatorID);
    bool readsCommandLine = cdo_operator_f2(operatorID);

    operator_input_arg(cdo_operator_enter(operatorID));

    auto const &argList = cdo_get_oper_argv();

    exprString = readsCommandLine ? exprs_from_argument(argList) : exprs_from_file(argList);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1(vlistID1);

    exprString = exprs_expand(exprString, varList1);
    if (Options::cdoVerbose) cdo_print(exprString);

    numVars1 = varList1.numVars();

    pointID = gridCreate(GRID_GENERIC, 1);
    zonalID = genZonalID(vlistID1);
    surfaceID = get_surface_ID(vlistID1);

    parse_param_init(parseArg, vlistID1, pointID, zonalID, surfaceID);

    auto &params = parseArg.params;
    params_init(parseArg.params, varList1);

    // Set all input variables to 'needed' if replacing is switched off
    for (int varID = 0; varID < numVars1; ++varID) parseArg.needed[varID] = !replacesVariables;

    // init function rand()
    std::srand(Options::Random_Seed);

    vartsID = params_add_ts(parseArg);
    parseArg.tsID = vartsID;
    params_add_coordinates(vlistID1, parseArg);

    CDO_parser_errorno = 0;
    parse_expressions();
    if (CDO_parser_errorno != 0) cdo_abort("Syntax error!");

    parseArg.init = false;

    if (Options::cdoVerbose)
      for (int varID = 0; varID < numVars1; ++varID)
        if (parseArg.needed[varID]) cdo_print("Needed var: %d %s", varID, params[varID].name);

    if (Options::cdoVerbose)
      for (int varID = 0, n = parseArg.numParams; varID < n; ++varID)
        cdo_print("var: %d %s ngp=%zu nlev=%zu coord=%c", varID, params[varID].name, params[varID].ngp, params[varID].nlev,
                  (params[varID].coord == 0) ? ' ' : params[varID].coord);

    varIDmap.resize(parseArg.numParams);

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));
    vlistClearFlag(vlistID1);
    if (!replacesVariables)
    {
      int pidx = 0;
      for (int varID = 0; varID < numVars1; ++varID)
      {
        params[varID].select = false;
        if (!params[varID].remove)
        {
          varIDmap[pidx++] = varID;
          auto nlevels = varList1.vars[varID].nlevels;
          // printf("Replace %d nlevs %d\n", varID, nlevels);
          for (int levID = 0; levID < nlevels; levID++) vlistDefFlag(vlistID1, varID, levID, true);
        }
      }
    }
    cdo_vlist_copy_flag(vlistID2, vlistID1);  // Copy global attributes

    add_new_output_params();

    if (Options::cdoVerbose)
    {
      for (int varID = 0; varID < numVars1; ++varID)
        if (parseArg.needed[varID]) cdo_print("needed: %d %s", varID, parseArg.params[varID].name);
      cdo_print("numVars1=%d, numVars2=%d", numVars1, numVars2);
    }

    numVars2 = vlistNvars(vlistID2);
    if (numVars2 == 0) cdo_abort("No output variable found!");

    allocate_params();

    read_coordinates();
    copy_coordinates();

    if (Options::cdoVerbose) vlistPrint(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    calendar = taxisInqCalendar(taxisID1);
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto &params = parseArg.params;
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      set_date_and_time(params[vartsID], calendar, tsID, vDateTime0, vDateTime);

      vDateTime0 = vDateTime;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      cdo_def_timestep(streamID2, tsID);

      // for (int varID = 0; varID < numVars1; ++varID) printf(">>> %s %d\n", params[varID].name.c_str(), params[varID].isValid);
      for (int varID = 0; varID < numVars1; ++varID) params[varID].isValid = true;
      for (int varID = 0; varID < numVars1; ++varID)
        if (tsID == 0 || params[varID].steptype != TIME_CONSTANT) params[varID].numMissVals = 0;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (parseArg.needed[varID])
        {
          auto offset = params[varID].ngp * levelID;
          auto vardata = params[varID].data + offset;
          size_t numMissVals;
          cdo_read_field(streamID1, vardata, &numMissVals);
          params[varID].numMissVals += numMissVals;

          if (numMissVals > 0) cdo_check_missval(params[varID].missval, params[varID].name);
        }
      }

      init_output_params(params);

      parse_expressions();

      for (int varID = 0; varID < numVars2; ++varID)
      {
        auto pidx = varIDmap[varID];
        if (tsID > 0 && params[pidx].steptype == TIME_CONSTANT) continue;

        auto missval = vlistInqVarMissval(vlistID2, varID);
        auto ngp = params[pidx].ngp;
        auto nlev = (int) params[pidx].nlev;
        for (int levelID = 0; levelID < nlev; ++levelID)
        {
          auto offset = ngp * levelID;
          double *vardata = params[pidx].data + offset;
          auto numMissVals = array_num_mv(ngp, vardata, missval);
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, vardata, numMissVals);
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

    params_delete(parseArg.params);
  }
};
