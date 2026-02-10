/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "process_int.h"
#include "printinfo.h"
#include "cdo_zaxis.h"

class Pinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Pinfo",
    .operators = { { "pinfo" }, { "pinfov" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Pinfo> registration = RegisterEntry<Pinfo>(module);

  int PINFO{}, PINFOV{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operatorID{};

  size_t imiss = 0;
  double arrmin{}, arrmax{}, arrmean{};

  VarList varList1;

public:
  void
  init() override
  {
    PINFO = module.get_id("pinfo");
    PINFOV = module.get_id("pinfov");

    (void) (PINFO);  // unused

    operatorID = cdo_operator_id();

    operator_check_argc(0);

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
    Varray<double> array1(varList1.gridsizeMax()), array2(varList1.gridsizeMax());

    int indg = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto vdateString = date_to_string(vDateTime.date);
      auto vtimeString = time_to_string(vDateTime.time);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        if (tsID == 0 && fieldID == 0)
        {
          if (operatorID == PINFOV)
            fprintf(stdout,
                    "   Rec :       Date  Time    Varname     Level    Size    Miss :     Minimum        Mean     Maximum\n");
          else
            fprintf(stdout, "   Rec :       Date  Time    Code  Level    Size    Miss :     Minimum        Mean     Maximum\n");
        }

        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto const &var = varList1.vars[varID];
        indg += 1;
        auto gridsize = var.gridsize;

        fprintf(stdout, "%6d :%s %s ", indg, vdateString.c_str(), vtimeString.c_str());
        if (operatorID == PINFOV)
          fprintf(stdout, "%-8s ", var.name.c_str());
        else
          fprintf(stdout, "%3d", var.code);

        auto level = cdo_zaxis_inq_level(var.zaxisID, levelID);
        fprintf(stdout, " %7g ", level);

        fprintf(stdout, "%7zu %7zu :", gridsize, numMissVals);

        if (gridInqType(var.gridID) == GRID_SPECTRAL || (gridsize == 1 && numMissVals == 0))
        {
          fprintf(stdout, "            %#12.5g\n", array1[0]);
        }
        else
        {
          if (numMissVals)
          {
            auto mmm = varray_min_max_mean_mv(array1, gridsize, var.missval);
            arrmin = mmm.min;
            arrmax = mmm.max;
            arrmean = mmm.mean;
            auto ivals = mmm.n;
            imiss = gridsize - ivals;
            gridsize = ivals;
          }
          else
          {
            auto mmm = varray_min_max_mean(array1, gridsize);
            arrmin = mmm.min;
            arrmax = mmm.max;
            arrmean = mmm.mean;
          }

          if (gridsize) { fprintf(stdout, "%#12.5g%#12.5g%#12.5g\n", arrmin, arrmean, arrmax); }
          else { fprintf(stdout, "                     nan\n"); }

          if (imiss != numMissVals && numMissVals) fprintf(stdout, "Found %zu of %zu missing values!\n", imiss, numMissVals);
        }

        varray_copy(gridsize, array1, array2);

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
  }
};
