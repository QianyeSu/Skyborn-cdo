/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_timer.h"
#include "param_conversion.h"
#include "process_int.h"
#include "util_files.h"

static void
print_stat(const char *sinfo, MemType memtype, int datatype, int filetype, off_t nvalues, double dataSize, double fileSize,
           double tw)
{
  nvalues /= 1000000;
  dataSize /= 1024. * 1024. * 1024.;

  cdo_print("%s Read %.1f GB of %d bit floats from %s %s, %.1f MVal/s", sinfo, dataSize, (memtype == MemType::Float) ? 32 : 64,
            cdo::datatype_to_cstr(datatype), cdo::filetype_to_cstr(filetype), (tw > 0) ? nvalues / tw : -1);

  fileSize /= 1024. * 1024. * 1024.;
  cdo_print("%s Read %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, fileSize, tw, (tw > 0) ? 1024 * fileSize / tw : -1);
}

class CDIread : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CDIread",
    .operators = { { "cdiread" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<CDIread> registration = RegisterEntry<CDIread>(module);

private:
  MemType memtype = Options::CDO_Memtype;
  int filetype = -1, datatype = -1;
  char sinfo[64];
  off_t nvalues = 0;
  double fileSize = 0, dataSize = 0;
  double runTimeSum = 0.0;

  int numRuns{};

public:
  void
  init() override
  {
    sinfo[0] = 0;

    if (Options::cdoVerbose) cdo_print("parameter: <nruns>");

    if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

    numRuns = (cdo_operator_argc() == 1) ? parameter_to_int(cdo_operator_argv(0)) : 1;
    if (numRuns < 0) numRuns = 0;
    if (numRuns > 99) numRuns = 99;

    if (Options::cdoVerbose) cdo_print("nruns      : %d", numRuns);

    // vlistDefNtsteps(vlistID, 1);
  }

  void
  run() override
  {
    for (int irun = 0; irun < numRuns; ++irun)
    {
      cdo::timer runTimer;
      dataSize = 0;
      nvalues = 0;

      auto streamID = cdo_open_read(0);
      auto vlistID = cdo_stream_inq_vlist(streamID);

      VarList varList(vlistID);

      filetype = cdo_inq_filetype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);

      auto gridsizeMax = varList.gridsizeMax();

      Varray<float> farray;
      Varray<double> darray;
      if (memtype == MemType::Float)
        farray.resize(gridsizeMax);
      else
        darray.resize(gridsizeMax);

      int tsID = 0;
      while (true)
      {
        cdo::timer stepTimer;
        auto numFields = cdo_stream_inq_timestep(streamID, tsID);
        if (numFields == 0) break;

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID);

          auto gridsize = varList.vars[varID].gridsize;
          nvalues += gridsize;

          size_t numMissVals;
          if (memtype == MemType::Float)
          {
            cdo_read_field_f(streamID, farray.data(), &numMissVals);
            dataSize += gridsize * 4;
          }
          else
          {
            cdo_read_field(streamID, darray.data(), &numMissVals);
            dataSize += gridsize * 8;
          }
        }

        if (Options::cdoVerbose) { cdo_print("Timestep %d: %.3f seconds", tsID + 1, stepTimer.elapsed()); }

        tsID++;
      }

      cdo_stream_close(streamID);

      auto runTime = runTimer.elapsed();
      runTimeSum += runTime;

      fileSize = (double) FileUtils::size(cdo_get_stream_name(0));

      if (numRuns > 1) std::snprintf(sinfo, sizeof(sinfo), "(run %d)", irun + 1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, dataSize, fileSize, runTime);
    }

    if (numRuns > 1) print_stat("(mean)", memtype, datatype, filetype, nvalues, dataSize, fileSize, runTimeSum / numRuns);
  }

  void
  close() override
  {
  }
};
