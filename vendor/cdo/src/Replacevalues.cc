/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setvals     setvals       Set list of old values to new values
      Setrtoc     setrtoc       Set range to new value
      Setrtoc2    setrtoc2      Set range to new value others to value2
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"

class Replacevalues : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Replacevalues",
    .operators = { { "setvals", 0, 0, "I1,O1,...,In,On", ReplacevaluesHelp },
                   { "setrtoc", 0, 0, "range(min,max),value", ReplacevaluesHelp },
                   { "setrtoc2", 0, 0, "range(min,max),value1,value2", ReplacevaluesHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Replacevalues> registration = RegisterEntry<Replacevalues>(module);

  int SETVALS{}, SETRTOC{}, SETRTOC2{};
  int nvals = 0;
  std::vector<double> fltarr;
  double rmin = 0, rmax = 0;
  double newval = 0, newval2 = 0;

  int operatorID{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

public:
  void
  init() override
  {
    SETVALS = module.get_id("setvals");
    SETRTOC = module.get_id("setrtoc");
    SETRTOC2 = module.get_id("setrtoc2");

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    if (operatorID == SETVALS)
    {
      fltarr = cdo_argv_to_fltarr(cdo_get_oper_argv());
      nvals = fltarr.size();
      if (nvals < 2) cdo_abort("Too few arguments!");
      if (nvals % 2 != 0) cdo_abort("Need pairs of arguments!");
      nvals = nvals / 2;
    }
    else if (operatorID == SETRTOC)
    {
      operator_check_argc(3);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
      newval = parameter_to_double(cdo_operator_argv(2));
    }
    else if (operatorID == SETRTOC2)
    {
      operator_check_argc(4);
      rmin = parameter_to_double(cdo_operator_argv(0));
      rmax = parameter_to_double(cdo_operator_argv(1));
      newval = parameter_to_double(cdo_operator_argv(2));
      newval2 = parameter_to_double(cdo_operator_argv(3));
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Varray<double> array(varList1.gridsizeMax());

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
        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);
        auto const &var = varList1.vars[varID];

        auto gridsize = var.gridsize;
        auto missval = var.missval;

        if (operatorID == SETVALS)
        {
          for (size_t i = 0; i < gridsize; ++i)
            if (fp_is_not_equal(array[i], missval))
            {
              for (int j = 0; j < nvals; ++j)
              {
                if (fp_is_equal(array[i], fltarr[j * 2]))
                {
                  array[i] = fltarr[j * 2 + 1];
                  break;
                }
              }
            }
        }
        else if (operatorID == SETRTOC)
        {
          for (size_t i = 0; i < gridsize; ++i)
            if (fp_is_not_equal(array[i], missval))
            {
              if (array[i] >= rmin && array[i] <= rmax) array[i] = newval;
            }
        }
        else if (operatorID == SETRTOC2)
        {
          for (size_t i = 0; i < gridsize; ++i)
            if (fp_is_not_equal(array[i], missval)) { array[i] = (array[i] >= rmin && array[i] <= rmax) ? newval : newval2; }
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, array.data(), numMissVals);
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
