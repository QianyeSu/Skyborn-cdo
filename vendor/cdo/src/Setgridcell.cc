/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setgridcell    setgridcell      Set grid cells to value
*/

#include <limits>
#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"

template <typename T>
static size_t
set_value(size_t gridsize, Varray<T> &array, double mv, double val)
{
  T missval = mv;
  T value = val;
  for (size_t i = 0; i < gridsize; ++i) { array[i] = value; }

  return varray_num_mv(gridsize, array, missval);
}

static void
set_value(Field &field, double value)
{
  auto func = [&](auto &v) { field.numMissVals = set_value(field.size, v, field.missval, value); };
  field_operation(func, field);
}

template <typename T>
static size_t
set_value(size_t gridsize, Varray<T> &array, double mv, double val, Varray<size_t> const &cells)
{
  T missval = mv;
  T value = val;
  auto numCells = cells.size();
  for (size_t i = 0; i < numCells; ++i) { array[cells[i] - 1] = value; }

  return varray_num_mv(gridsize, array, missval);
}

static void
set_value(Field &field, double value, Varray<size_t> const &cells)
{
  auto func = [&](auto &v) { field.numMissVals = set_value(field.size, v, field.missval, value, cells); };
  field_operation(func, field);
}

static void
read_index_from_maskfile(std::string const &maskFileName, Varray<size_t> &cells)
{
  std::vector<bool> cdo_read_mask(const char *maskFileName);
  auto mask = cdo_read_mask(maskFileName.c_str());
  size_t nind = 0;
  for (size_t i = 0; i < mask.size(); ++i)
    if (mask[i]) nind++;
  if (nind == 0) cdo_abort("Mask is empty!");

  cells.resize(nind);
  nind = 0;
  for (size_t i = 0; i < mask.size(); ++i)
    if (mask[i]) cells[nind++] = i + 1;

  if (nind == 0) cdo_abort("Mask file %s generates no input!", cdo_operator_argv(0));
}

static void
setgridcell_get_parameter(double &constant, Varray<size_t> &cells, std::vector<std::string> &varnames, std::string &maskfile)
{
  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &values = kv.values;
      auto const &value = kv.values[0];
      int nvalues = kv.nvalues;

      if (key == "value")
      {
        if (nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
        constant = parameter_to_double(value);
      }
      else if (key == "cell")
      {
        cells.resize(nvalues);
        for (int i = 0; i < nvalues; ++i) cells[i] = parameter_to_size_t(values[i]);
      }
      else if (key == "name")
      {
        varnames.resize(nvalues);
        for (int i = 0; i < nvalues; ++i) varnames[i] = values[i];
      }
      else if (key == "mask")
      {
        if (nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
        maskfile = value;
      }
      else { cdo_abort("Invalid parameter key >%s<!", key); }
    }
  }
}

class Setgridcell : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setgridcell",
    .operators = { { "setgridcell", 0, 0, "value=constant[,cell=gridcellindices(1-N)]", SetgridcellHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setgridcell> registration = RegisterEntry<Setgridcell>(module);

  Field field;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numVars{};
  VarList varList1;
  Varray<bool> selectVars;

  std::vector<std::string> varnames;
  Varray<size_t> cells;
  double value = DBL_MAX;
  size_t maxIndex = 0;
  size_t numCells{};

public:
  void
  init() override
  {
    operator_input_arg(cdo_operator_enter(0));

    auto nparam = cdo_operator_argc();
    if (nparam == 0) cdo_abort("Parameter missing!");

    std::string maskFileName;
    setgridcell_get_parameter(value, cells, varnames, maskFileName);
    if (fp_is_equal(value, DBL_MAX)) cdo_abort("Parameter <values> not set!");

    if (cells.size() && maskFileName.size()) cdo_abort("Too many parameter, choose either cell or mask!");

    if (maskFileName.size()) read_index_from_maskfile(maskFileName, cells);
    numCells = cells.size();
    size_t minIndex = std::numeric_limits<size_t>::max();

    for (size_t i = 0; i < numCells; ++i)
    {
      minIndex = std::min(minIndex, cells[i]);
      maxIndex = std::max(maxIndex, cells[i]);
    }

    if (numCells && minIndex == 0) cdo_abort("Min cell index is 0, muss be >= 1!");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    selectVars = Varray<bool>(numVars, !varnames.size());
    if (varnames.size())
    {
      for (auto const &varname : varnames)
      {
        auto varFound = false;
        for (int varID = 0; varID < numVars; ++varID)
        {
          if (varname == varList1.vars[varID].name)
          {
            selectVars[varID] = true;
            varFound = true;
            break;
          }
        }
        if (!varFound) cdo_abort("Variable %s not found!", varname);
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
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
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        if (selectVars[varID])
        {
          if (numCells)
          {
            if (maxIndex > field.size) cdo_abort("Max cell index (%zu) > gridsize (%zu)!", maxIndex, field.size);
            set_value(field, value, cells);
          }
          else { set_value(field, value); }
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
