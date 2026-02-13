/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Condc      ifthenc         If then constant
      Condc      ifnotthenc      If not then constant
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"

template <typename T>
static void
operator_IFTHENC(size_t n, double missval, Varray<T> &v, double value)
{
  T mv = missval;
  T cVal = value;
  constexpr T nullVal{ 0 };
  if (std::isnan(mv))
    for (size_t i = 0; i < n; ++i) v[i] = (fp_is_not_equal(v[i], mv) && fp_is_not_equal(v[i], nullVal)) ? cVal : mv;
  else
    for (size_t i = 0; i < n; ++i) v[i] = (!is_equal(v[i], mv) && !is_equal(v[i], nullVal)) ? cVal : mv;
}

template <typename T>
static void
operator_IFNOTTHENC(size_t n, double missval, Varray<T> &v, double value)
{
  T mv = missval;
  T cVal = value;
  constexpr T nullVal{ 0 };
  if (std::isnan(mv))
    for (size_t i = 0; i < n; ++i) v[i] = (fp_is_not_equal(v[i], mv) && fp_is_equal(v[i], nullVal)) ? cVal : mv;
  else
    for (size_t i = 0; i < n; ++i) v[i] = (!is_equal(v[i], mv) && is_equal(v[i], nullVal)) ? cVal : mv;
}

class Condc : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Condc",
    .operators = { { "ifthenc", CondcHelp }, { "ifnotthenc", CondcHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Condc> registration = RegisterEntry<Condc>();

private:
  int IFTHENC{}, IFNOTTHENC{};
  int taxisID1{ CDI_UNDEFID };
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID2{ CDI_UNDEFID };

  int operatorID{};
  double rconst{};

  VarList varList1{};

public:
  void
  init() override
  {
    IFTHENC = module.get_id("ifthenc");
    IFNOTTHENC = module.get_id("ifnotthenc");

    operatorID = cdo_operator_id();

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

        auto numMissVals = field.numMissVals;
        if (numMissVals > 0) cdo_check_missval(field.missval, var.name);

        auto func_ifthenc = [&](auto &v, auto n, double mv) { operator_IFTHENC(n, mv, v, rconst); };
        auto func_ifnotthenc = [&](auto &v, auto n, double mv) { operator_IFNOTTHENC(n, mv, v, rconst); };

        // clang-format off
        if      (operatorID == IFTHENC)    field_operation(func_ifthenc, field, field.size, field.missval);
        else if (operatorID == IFNOTTHENC) field_operation(func_ifnotthenc, field, field.size, field.missval);
        else cdo_abort("Operator not implemented!");
        // clang-format on

        if (numMissVals > 0) field_num_mv(field);

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
