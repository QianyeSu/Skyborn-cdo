/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Compc      eqc             Equal constant
      Compc      nec             Not equal constant
      Compc      lec             Less equal constant
      Compc      ltc             Less then constant
      Compc      gec             Greater equal constant
      Compc      gtc             Greater then constant
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static auto func_compc = [](auto hasMissvals, auto n, auto mv, auto &v, const auto cVal, auto binary_operator)
{
  if (hasMissvals)
  {
    auto constantIsMissval = fp_is_equal(cVal, mv);
    if (std::isnan(mv))
      for (size_t i = 0; i < n; ++i) v[i] = (fp_is_equal(v[i], mv) || constantIsMissval) ? mv : binary_operator(v[i], cVal);
    else
      for (size_t i = 0; i < n; ++i) v[i] = (is_equal(v[i], mv) || constantIsMissval) ? mv : binary_operator(v[i], cVal);
  }
  else
  {
    for (size_t i = 0; i < n; ++i) v[i] = binary_operator(v[i], cVal);
  }
};

template <typename T>
static void
comp_function(int operFunc, bool hasMissvals, size_t ngp, double missval, Varray<T> &v, double value)
{
  T mv = missval;
  T rconst = value;
  // clang-format off
  if      (operFunc == FieldFunc_EQ) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_EQ);
  else if (operFunc == FieldFunc_NE) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_NE);
  else if (operFunc == FieldFunc_LE) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_LE);
  else if (operFunc == FieldFunc_LT) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_LT);
  else if (operFunc == FieldFunc_GE) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_GE);
  else if (operFunc == FieldFunc_GT) func_compc(hasMissvals, ngp, mv, v, rconst, binary_op_GT);
  else cdo_abort("Operator not implemented!");
  // clang-format on
}

static void
comp_function(Field &field, int operFunc, bool hasMissvals, double rconst)
{
  auto func = [&](auto &v, auto n, double mv) { comp_function(operFunc, hasMissvals, n, mv, v, rconst); };
  field_operation(func, field, field.size, field.missval);
}

class Compc : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Compc",
    .operators = { { "eqc", FieldFunc_EQ, 0, CompcHelp },
                   { "nec", FieldFunc_NE, 0, CompcHelp },
                   { "lec", FieldFunc_LE, 0, CompcHelp },
                   { "ltc", FieldFunc_LT, 0, CompcHelp },
                   { "gec", FieldFunc_GE, 0, CompcHelp },
                   { "gtc", FieldFunc_GT, 0, CompcHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Compc> registration = RegisterEntry<Compc>();

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operFunc{};
  double rconst{};

  VarList varList1;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operFunc = cdo_operator_f1(operatorID);

    operator_input_arg("constant value");
    rconst = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    vlistDestroy(vlistID2);
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
        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        auto missval = field.missval;
        auto numMissVals = field.numMissVals;
        auto hasMissvals = (numMissVals > 0 || fp_is_equal(rconst, missval));

        if (numMissVals > 0) cdo_check_missval(missval, var.name);

        comp_function(field, operFunc, hasMissvals, rconst);

        if (hasMissvals > 0) field_num_mv(field);

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
