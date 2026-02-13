/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "statistic.h"

class Tests : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Tests",
    .operators = { { "normal" },
                   { "studentt", 0, 0, "degree of freedom" },
                   { "chisquare", 0, 0, "degree of freedom" },
                   { "beta", 0, 0, "p and q" },
                   { "fisher", 0, 0, "degree of freedom of nominator and of denominator" } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Tests> registration = RegisterEntry<Tests>();

  int NORMAL, STUDENTT, CHISQUARE, BETA, FISHER;
  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operatorID;

  double degreeOfFreedom = 0, p = 0, q = 0, n = 0, d = 0;

  VarList varList1;

public:
  void
  init() override
  {
    NORMAL = module.get_id("normal");
    STUDENTT = module.get_id("studentt");
    CHISQUARE = module.get_id("chisquare");
    BETA = module.get_id("beta");
    FISHER = module.get_id("fisher");

    operatorID = cdo_operator_id();

    if (operatorID == STUDENTT || operatorID == CHISQUARE)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(1);

      degreeOfFreedom = parameter_to_double(cdo_operator_argv(0));
      if (degreeOfFreedom <= 0) cdo_abort("degree of freedom must be positive!");
    }
    else if (operatorID == BETA)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(2);

      p = parameter_to_double(cdo_operator_argv(0));
      q = parameter_to_double(cdo_operator_argv(1));

      if (p <= 0 || q <= 0) cdo_abort("p and q must be positive!");
    }
    else if (operatorID == FISHER)
    {
      operator_input_arg(cdo_operator_enter(operatorID));

      operator_check_argc(2);

      n = parameter_to_double(cdo_operator_argv(0));
      d = parameter_to_double(cdo_operator_argv(1));
      if (n <= 0 || d <= 0) cdo_abort("both degrees must be positive!");
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

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

    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array1(gridsizeMax);
    Varray<double> array2(gridsizeMax);

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
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto gridsize = varList1.vars[varID].gridsize;
        auto missval = varList1.vars[varID].missval;

        if (operatorID == NORMAL)
        {
          for (size_t i = 0; i < gridsize; ++i) array2[i] = fp_is_equal(array1[i], missval) ? missval : cdo::normal(array1[i]);
        }
        else if (operatorID == STUDENTT)
        {
          for (size_t i = 0; i < gridsize; ++i)
            array2[i] = fp_is_equal(array1[i], missval) ? missval : cdo::student_t(degreeOfFreedom, array1[i]);
        }
        else if (operatorID == CHISQUARE)
        {
          for (size_t i = 0; i < gridsize; ++i)
            array2[i] = fp_is_equal(array1[i], missval) ? missval : cdo::chi_square(degreeOfFreedom, array1[i]);
        }
        else if (operatorID == BETA)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (array1[i] < 0 || array1[i] > 1) cdo_abort("Value out of range (0-1)!");

            array2[i] = fp_is_equal(array1[i], missval) ? missval : cdo::beta_distr(p, q, array1[i]);
          }
        }
        else if (operatorID == FISHER)
        {
          for (size_t i = 0; i < gridsize; ++i)
            array2[i] = fp_is_equal(array1[i], missval) ? missval : cdo::fisher(n, d, array1[i]);
        }
        else { cdo_abort("Internal problem, operator not implemented!"); }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, array2.data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
