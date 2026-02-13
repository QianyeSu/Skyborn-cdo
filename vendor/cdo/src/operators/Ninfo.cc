/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ninfo      npar            Number of parameters
      Ninfo      nlevel          Number of levels
      Ninfo      nyear           Number of years
      Ninfo      nmon            Number of months
      Ninfo      ndate           Number of dates
      Ninfo      ntime           Number of timesteps
      Ninfo      ngridpoints     Number of gridpoints
      Ninfo      ngrids          Number of grids
*/

#include <cdi.h>

#include "process_int.h"

class Ninfo : public Process
{
  enum
  {
    NYEAR,
    NMON,
    NDATE,
    NTIME,
    NPAR,
    NLEVEL,
    NGRIDPOINTS,
    NGRIDS
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ninfo",
    .operators = { { "nyear", NYEAR, 0, NinfoHelp },
                   { "nmon", NMON, 0, NinfoHelp },
                   { "ndate", NDATE, 0, NinfoHelp },
                   { "ntime", NTIME, 0, NinfoHelp },
                   { "ncode", NinfoHelp },
                   { "npar", NPAR, 0, NinfoHelp },
                   { "nlevel", NLEVEL, 0, NinfoHelp },
                   { "ngridpoints", NGRIDPOINTS, 0, NinfoHelp },
                   { "ngrids", NGRIDS, 0, NinfoHelp } },
    .aliases = { { "nvar", "npar" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Ninfo> registration = RegisterEntry<Ninfo>();

private:
  int operfunc{};
  CdoStreamID streamID{};
  int taxisID{};

  VarList varList{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID = cdo_open_read(0);
    auto vlistID = cdo_stream_inq_vlist(streamID);
    varList = VarList(vlistID);

    taxisID = vlistInqTaxis(vlistID);
  }

  void
  run() override
  {
    auto numVars = varList.numVars();
    auto numSteps = varList.numSteps();
    auto numGrids = varList.numGrids();

    switch (operfunc)
    {
      case NYEAR:
      {
        int nyear = 0;
        if (numSteps != 0)
        {
          int year0 = 0;
          int tsID = 0;
          while (cdo_stream_inq_timestep(streamID, tsID))
          {
            int year = taxisInqVdatetime(taxisID).date.year;
            if (tsID == 0 || year0 != year)
            {
              year0 = year;
              nyear++;
            }

            tsID++;
          }
        }
        std::fprintf(stdout, "%d\n", nyear);
        break;
      }
      case NMON:
      {
        int nmonth = 0;
        if (numSteps != 0)
        {
          int month0 = 0;
          int tsID = 0;
          while (cdo_stream_inq_timestep(streamID, tsID))
          {
            int month = taxisInqVdatetime(taxisID).date.month;
            if (tsID == 0 || month0 != month)
            {
              month0 = month;
              nmonth++;
            }

            tsID++;
          }
        }
        std::fprintf(stdout, "%d\n", nmonth);
        break;
      }
      case NDATE:
      {
        CdiDate date0{};
        int ndate = 0;
        if (numSteps != 0)
        {
          int tsID = 0;
          while (cdo_stream_inq_timestep(streamID, tsID))
          {
            auto vDate = taxisInqVdatetime(taxisID).date;
            if (tsID == 0 || !cdiDate_isEQ(date0, vDate))
            {
              date0 = vDate;
              ndate++;
            }

            tsID++;
          }
        }
        std::fprintf(stdout, "%d\n", ndate);
        break;
      }
      case NTIME:
      {
        int tsID = (numSteps > 0) ? numSteps : 0;
        if (tsID == 0)
          while (cdo_stream_inq_timestep(streamID, tsID)) tsID++;
        std::fprintf(stdout, "%d\n", tsID);
        break;
      }
      case NPAR: std::fprintf(stdout, "%d\n", numVars); break;
      case NLEVEL:
        for (auto const &var : varList.vars) { std::fprintf(stdout, "%d\n", var.nlevels); }
        break;
      case NGRIDPOINTS:
        for (auto const &var : varList.vars) { std::fprintf(stdout, "%zu\n", var.gridsize); }
        break;
      case NGRIDS: std::fprintf(stdout, "%d\n", numGrids); break;
      default: cdo_abort("operator not implemented!"); break;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
