/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithc     addc            Add by constant
      Arithc     subc            Subtract by constant
      Arithc     mulc            Multiply by constant
      Arithc     divc            Divide by constant
      Arithc     minc            Minimum of a field and a constant
      Arithc     maxc            Maximum of a field and a constant
      Arithc     mod             Modulo operator
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "field_functions.h"

static void
fill_vars(VarList const &varList, std::vector<bool> &vars)
{
  auto numVars = varList.numVars();

  if (Options::cdo_num_varnames() > 0)
    {
      auto found = false;
      for (int varID = 0; varID < numVars; ++varID)
        {
          vars[varID] = false;
          for (size_t i = 0; i < Options::cdo_num_varnames(); ++i)
            if (varList.vars[varID].name == Options::cdoVarnames[i])
              {
                vars[varID] = true;
                found = true;
                break;
              }
        }

      if (!found) cdo_abort("Variable %s%s not found!", Options::cdoVarnames[0], (Options::cdo_num_varnames() > 1) ? ",..." : "");
    }
  else
    {
      for (int varID = 0; varID < numVars; ++varID) vars[varID] = true;
    }
}

class Arithc : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Arithc",
    .operators = { { "addc", FieldFunc_Add, 1, "constant value", ArithcHelp },
                   { "subc", FieldFunc_Sub, 1, "constant value", ArithcHelp },
                   { "mulc", FieldFunc_Mul, 1, "constant value", ArithcHelp },
                   { "divc", FieldFunc_Div, 1, "constant value", ArithcHelp },
                   { "minc", FieldFunc_Min, 0, "constant value", ArithcHelp },
                   { "maxc", FieldFunc_Max, 0, "constant value", ArithcHelp },
                   { "mod", FieldFunc_Mod, 0, "divisor", ArithcHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Arithc> registration = RegisterEntry<Arithc>(module);

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };

  int operfunc;
  double rconstcplx[2];
  double rconst;

  VarList varList1;
  std::vector<bool> vars;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    bool opercplx = cdo_operator_f2(operatorID);

    operator_input_arg(cdo_operator_enter(operatorID));
    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
    if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
    rconst = parameter_to_double(cdo_operator_argv(0));

    rconstcplx[0] = rconst;
    rconstcplx[1] = 0.0;

    if (cdo_operator_argc() == 2) rconstcplx[1] = parameter_to_double(cdo_operator_argv(1));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    varList1 = VarList(vlistID1);

    auto numVars = vlistNvars(vlistID1);
    // for (auto &var : varList1.vars) var.memType = MemType::Double;

    vars = std::vector<bool>(numVars);
    fill_vars(varList1, vars);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    vlistDestroy(vlistID2);

    auto nwpv = (vlistNumber(vlistID1) == CDI_COMP) ? 2 : 1;
    if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");
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
            field.init(varList1.vars[varID]);
            cdo_read_field(streamID1, field);

            if (vars[varID])
              {
                if (field.nwpv == 2)
                  fieldc_function_complex(field, rconstcplx, operfunc);
                else
                  fieldc_function(field, rconst, operfunc);

                // recalculate number of missing values
                field_num_mv(field);
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
