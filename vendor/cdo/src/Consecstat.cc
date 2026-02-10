/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: 2010 Ralf Mueller

*/

/*
   This module contains the following operators:

      Consectstep  consecsum  For each timestep, the current number of
                              onsecutive timsteps is counted
      Consectstep  consects   For each period of consecutive timesteps, only
                              count its length + last contributing timesteps

   =============================================================================
   Created:  04/08/2010 11:58:01 AM
    Author:  Ralf Mueller (ram), ralf.mueller@mpimet.mpg.de
   Company:  Max-Planck-Institute for Meteorology
   =============================================================================
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "field_functions.h"

#define SWITCHWARN "Hit default case! This should never happen (%s).\n"

static void
selEndOfPeriod(Field &periods, Field const &history, Field const &current, int isLastTimestep)
{
  auto pmissval = periods.missval;
  auto &parray = periods.vec_d;
  auto hmissval = history.missval;
  auto const &harray = history.vec_d;
  auto cmissval = current.missval;
  auto const &carray = current.vec_d;

  auto len = gridInqSize(periods.grid);
  if (len != gridInqSize(current.grid) || (gridInqSize(current.grid) != gridInqSize(history.grid)))
    cdo_abort("Fields have different gridsize (%s)", __func__);

  if (!isLastTimestep)
    {
      if (current.numMissVals || history.numMissVals)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i)
            {
              if (fp_is_not_equal(harray[i], hmissval))
                {
                  if (fp_is_not_equal(carray[i], cmissval))
                    parray[i] = (fp_is_equal(carray[i], 0.0) && is_not_equal(harray[i], 0.0)) ? harray[i] : pmissval;
                  else  // fp_is_equal(carray[i], cmissval)
                    parray[i] = (is_not_equal(harray[i], 0.0)) ? harray[i] : pmissval;
                }
              else /* fp_is_equal(harray[i], hmissval) */ { parray[i] = pmissval; }
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < len; ++i)
            {
              parray[i] = (fp_is_equal(carray[i], 0.0) && is_not_equal(harray[i], 0.0)) ? harray[i] : pmissval;
            }
        }
    }
  else
    {
      if (current.numMissVals)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i)
            {
              if (!fp_is_equal(carray[i], cmissval))
                parray[i] = (fp_is_equal(carray[i], 0.0)) ? pmissval : carray[i];
              else  // fp_is_equal(carray[i], cmissval)
                parray[i] = pmissval;
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (size_t i = 0; i < len; ++i) parray[i] = fp_is_equal(carray[i], 0.0) ? pmissval : carray[i];
        }
    }

  periods.numMissVals = varray_num_mv(len, parray, pmissval);
}

class Consecstat : public Process
{
  enum
  {
    CONSECSUM,
    CONSECTS
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Consecstat",
    // clang-format off
    .operators = { { "consects",  CONSECTS,  0, ConsecstatHelp },
                   { "consecsum", CONSECSUM, 0, "refval", ConsecstatHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Consecstat> registration = RegisterEntry<Consecstat>(module);

  CdiDateTime vDateTime{};
  CdiDateTime histDateTime{};
  double refval = 0.0;

  CdoStreamID istreamID;
  CdoStreamID ostreamID;

  int ivlistID;
  int itaxisID;

  int ovlistID;
  int otaxisID;

  int operfunc;

  VarList varList1;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    if (operfunc == CONSECSUM)
      if (cdo_operator_argc() > 0) refval = parameter_to_double(cdo_operator_argv(0));

    istreamID = cdo_open_read(0);

    ivlistID = cdo_stream_inq_vlist(istreamID);
    itaxisID = vlistInqTaxis(ivlistID);
    ovlistID = vlistDuplicate(ivlistID);
    otaxisID = taxisDuplicate(itaxisID);
    vlistDefTaxis(ovlistID, otaxisID);

    varList1 = VarList(ivlistID);
    for (auto &var : varList1.vars) var.memType = MemType::Double;

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyString(ovlistID, varID, CDI_KEY_UNITS, "steps");  // TODO

    ostreamID = cdo_open_write(1);
    cdo_def_vlist(ostreamID, ovlistID);
  }

  void
  run() override
  {
    Field field;

    FieldVector2D varsData, histData, periodsData;
    field2D_init(varsData, varList1, FIELD_VEC, 0);
    if (operfunc == CONSECTS) field2D_init(histData, varList1, FIELD_VEC);
    if (operfunc == CONSECTS) field2D_init(periodsData, varList1, FIELD_VEC);

    int itsID = 0;
    int otsID = 0;
    while (true)
      {
        auto numFields = cdo_stream_inq_timestep(istreamID, itsID);
        if (numFields == 0) break;

        vDateTime = taxisInqVdatetime(itaxisID);
        switch (operfunc)
          {
          case CONSECSUM:
            taxisDefVdatetime(otaxisID, vDateTime);
            cdo_def_timestep(ostreamID, otsID);
            break;
          case CONSECTS:
            if (itsID != 0)
              {
                taxisDefVdatetime(otaxisID, histDateTime);
                cdo_def_timestep(ostreamID, otsID - 1);
              }
            break;
          default: printf(SWITCHWARN, __func__); break;
          }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = cdo_inq_field(istreamID);
            auto const &var1 = varList1.vars[varID];
            field.init(var1);
            cdo_read_field(istreamID, field);

            auto &varData = varsData[varID][levelID];
            field2_sumtr(varData, field, refval);

            switch (operfunc)
              {
              case CONSECSUM:
                cdo_def_field(ostreamID, varID, levelID);
                cdo_write_field(ostreamID, varData);
                break;
              case CONSECTS:
                if (itsID != 0)
                  {
                    selEndOfPeriod(periodsData[varID][levelID], histData[varID][levelID], varData, false);
                    cdo_def_field(ostreamID, varID, levelID);
                    cdo_write_field(ostreamID, periodsData[varID][levelID]);
                  }
                histData[varID][levelID].vec_d = varData.vec_d;
                break;
              default: printf(SWITCHWARN, __func__); break;
              }
          }

        histDateTime = vDateTime;
        itsID++;
        otsID++;
      }

    if (operfunc == CONSECTS)  // Save the last timestep
      {
        taxisDefVdatetime(otaxisID, vDateTime);
        cdo_def_timestep(ostreamID, otsID - 1);

        auto numVars = varList1.numVars();
        for (int varID = 0; varID < numVars; ++varID)
          {
            auto nlevels = varList1.vars[varID].nlevels;
            for (int levelID = 0; levelID < nlevels; ++levelID)
              {
                selEndOfPeriod(periodsData[varID][levelID], histData[varID][levelID], varsData[varID][levelID], true);
                cdo_def_field(ostreamID, varID, levelID);
                cdo_write_field(ostreamID, periodsData[varID][levelID]);
              }
          }
      }
  }

  void
  close() override
  {
    cdo_stream_close(istreamID);
    cdo_stream_close(ostreamID);
  }
};
