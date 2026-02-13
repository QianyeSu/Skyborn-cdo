/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setmiss    setmissval      Set a new missing value
      Setmiss    setctomiss      Set constant to missing value
      Setmiss    setmisstoc      Set missing value to constant
      Setmiss    setrtomiss      Set range to missing value
      Setmiss    setvrange       Set range of valid value
*/

#include <cmath>
#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "progress.h"

template <typename T>
static size_t
set_missval(size_t gridsize, Varray<T> &array, double mv, double mvNew)
{
  T missval = mv;
  T missvalNew = mvNew;
  size_t numMissVals = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (fp_is_equal(array[i], missval) || fp_is_equal(array[i], (float) missval) || fp_is_equal(array[i], missvalNew)
        || fp_is_equal(array[i], (float) missvalNew))
    {
      array[i] = missvalNew;
      numMissVals++;
    }

  return numMissVals;
}

static void
set_missval(Field &field, double missvalNew)
{
  auto func = [&](auto &v, auto size, double mv) { return set_missval(size, v, mv, missvalNew); };
  field.numMissVals = field_operation(func, field, field.size, field.missval);
}

template <typename T>
static size_t
set_const_to_miss(size_t gridsize, Varray<T> &array, double _missval, double _rconst)
{
  T missval = _missval;
  T rconst = _rconst;
  size_t numMissVals = 0;
  if (std::isnan(_rconst))
  {
    for (size_t i = 0; i < gridsize; ++i)
      if (std::isnan(array[i]))
      {
        array[i] = missval;
        numMissVals++;
      }
  }
  else if (std::isinf(_rconst))
  {
    for (size_t i = 0; i < gridsize; ++i)
      if (std::isinf(array[i]))
      {
        array[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < gridsize; ++i)
      if (fp_is_not_equal(array[i], missval) && (fp_is_equal(array[i], rconst) || fp_is_equal(array[i], (float) rconst)))
      {
        array[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}

static void
set_const_to_miss(Field &field, double rconst)
{
  auto func = [&](auto &v, auto size, double mv) { return set_const_to_miss(size, v, mv, rconst); };
  field.numMissVals += field_operation(func, field, field.size, field.missval);
}

template <typename T>
static size_t
set_miss_to_const(size_t gridsize, Varray<T> &array, double _missval, double _rconst)
{
  T missval = _missval;
  T rconst = _rconst;
  for (size_t i = 0; i < gridsize; ++i)
    if (fp_is_equal(array[i], missval) || fp_is_equal(array[i], (float) missval)) { array[i] = rconst; }

  return 0;
}

void
set_miss_to_const(Field &field, double rconst)
{
  auto func = [&](auto &v, auto size, double mv) { return set_miss_to_const(size, v, mv, rconst); };
  field.numMissVals = field_operation(func, field, field.size, field.missval);
}

template <typename T>
static size_t
set_range_to_miss(size_t gridsize, Varray<T> &array, double _missval, double _rmin, double _rmax)
{
  T missval = _missval;
  T rmin = _rmin;
  T rmax = _rmax;
  size_t numMissVals = 0;
  for (size_t i = 0; i < gridsize; ++i)
    if (array[i] >= rmin && array[i] <= rmax)
    {
      array[i] = missval;
      numMissVals++;
    }

  return numMissVals;
}

static void
set_range_to_miss(Field &field, double rmin, double rmax)
{
  auto func = [&](auto &v, auto size, double mv) { return set_range_to_miss(size, v, mv, rmin, rmax); };
  field.numMissVals += field_operation(func, field, field.size, field.missval);
}

template <typename T>
static size_t
set_valid_range(size_t gridsize, Varray<T> &array, double _missval, double _rmin, double _rmax)
{
  T missval = _missval;
  T rmin = _rmin;
  T rmax = _rmax;
  for (size_t i = 0; i < gridsize; ++i)
    if (array[i] < rmin || array[i] > rmax) array[i] = missval;

  auto numMissVals = varray_num_mv(gridsize, array, missval);

  return numMissVals;
}

static void
set_valid_range(Field &field, double rmin, double rmax)
{
  auto func = [&](auto &v, auto size, double mv) { return set_valid_range(size, v, mv, rmin, rmax); };
  field.numMissVals = field_operation(func, field, field.size, field.missval);
}

class Setmiss : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setmiss",
    .operators = { { "setmissval", 0, 0, "missing value", SetmissHelp },
                   { "setctomiss", 0, 0, "constant", SetmissHelp },
                   { "setmisstoc", 0, 0, "constant", SetmissHelp },
                   { "setrtomiss", 0, 0, "range(min,max)", SetmissHelp },
                   { "setvrange", 0, 0, "range(min,max)", SetmissHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setmiss> registration = RegisterEntry<Setmiss>();

private:
  int SETMISSVAL{}, SETCTOMISS{}, SETMISSTOC{}, SETRTOMISS{}, SETVRANGE{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  int operatorID{};

  double rconst = 0.0;
  double rmin = 0.0;
  double rmax = 0.0;
  double newMissval = 0.0;

public:
  void
  init() override
  {
    SETMISSVAL = module.get_id("setmissval");
    SETCTOMISS = module.get_id("setctomiss");
    SETMISSTOC = module.get_id("setmisstoc");
    SETRTOMISS = module.get_id("setrtomiss");
    SETVRANGE = module.get_id("setvrange");

    operatorID = cdo_operator_id();

    if (operatorID == SETMISSVAL)
    {
      operator_check_argc(1);
      newMissval = parameter_to_double(cdo_operator_argv(0));
    }
    else if (operatorID == SETCTOMISS || operatorID == SETMISSTOC)
    {
      operator_check_argc(1);
      rconst = parameter_to_double(cdo_operator_argv(0));
    }
    else
    {
      operator_check_argc(2);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    if (operatorID == SETMISSVAL)
    {
      auto numVars = vlistNvars(vlistID2);
      for (int varID = 0; varID < numVars; ++varID) vlistDefVarMissval(vlistID2, varID, newMissval);
    }
    else if (operatorID == SETMISSTOC)
    {
      auto numVars = vlistNvars(vlistID2);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (fp_is_equal(rconst, var.missval))
        {
          cdo_warning("Missing value and constant have the same value!");
          break;
        }
      }
    }

    /*
    if (operatorID == SETVRANGE)
      {
        double range[2] = {rmin, rmax};
        nvars = vlistNvars(vlistID2);
        for (varID = 0; varID < nvars; ++varID)
          cdiDefAttFlt(vlistID2, varID, "valid_range", CDI_DATATYPE_FLT64, 2, range);
      }
    */
    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        // clang-format off
        if      (operatorID == SETMISSVAL) set_missval(field, newMissval);
        else if (operatorID == SETCTOMISS) set_const_to_miss(field, rconst);
        else if (operatorID == SETMISSTOC) set_miss_to_const(field, rconst);
        else if (operatorID == SETRTOMISS) set_range_to_miss(field, rmin, rmax);
        else if (operatorID == SETVRANGE)  set_valid_range(field, rmin, rmax);
        // clang-format on

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
