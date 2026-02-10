/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

class Histogram : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Histogram",
    .operators = { { "histcount", HistogramHelp },
                   { "histsum", HistogramHelp },
                   { "histmean", HistogramHelp },
                   { "histfreq", HistogramHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Histogram> registration = RegisterEntry<Histogram>(module);

private:
  int HISTCOUNT{}, HISTSUM{}, HISTMEAN{}, HISTFREQ{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;

  int numBins{};

  std::vector<double> fltarr;
  int operatorID{};

public:
  void
  init() override
  {
    HISTCOUNT = module.get_id("histcount");
    HISTSUM = module.get_id("histsum");
    HISTMEAN = module.get_id("histmean");
    HISTFREQ = module.get_id("histfreq");

    (void) (HISTSUM);  // unused

    operatorID = cdo_operator_id();

    operator_input_arg("bins");

    fltarr = cdo_argv_to_fltarr(cdo_get_oper_argv());
    numBins = fltarr.size() - 1;
    if (numBins < 1) cdo_abort("Too few arguments!");

    if (Options::cdoVerbose)
    {
      printf("numBins = %d\n", numBins);
      for (int i = 0; i < numBins; ++i) printf("flt %d = %g\n", i + 1, fltarr[i]);
    }

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);

    auto vlistID2 = vlistDuplicate(vlistID1);

    /* create zaxis for output bins */
    auto zaxisID2 = zaxisCreate(ZAXIS_GENERIC, numBins);

    zaxisDefLevels(zaxisID2, &fltarr[0]);
    zaxisDefLbounds(zaxisID2, &fltarr[0]);
    zaxisDefUbounds(zaxisID2, &fltarr[1]);

    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, "bin");
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_LONGNAME, "histogram bins");
    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_UNITS, "level");

    // check zaxis: only 2D fields allowed
    auto numZaxes = varList1.numZaxes();
    for (int index = 0; index < numZaxes; ++index)
    {
      auto zaxisID = vlistZaxis(vlistID1, index);
      auto numLevels = zaxisInqSize(zaxisID);
      if (numLevels > 1) cdo_abort("Found 3D field with %d levels. Only 2D fields allowed!", numLevels);
      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
    }

    streamID2 = cdo_open_write(1);

    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    cdo_def_vlist(streamID2, vlistID2);

    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    auto numVars = varList2.numVars();
    Varray2D<double> vardata(numVars);
    Varray2D<double> varcount(numVars);
    Varray2D<double> vartcount(numVars);
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = varList1.vars[varID].gridsize;
      vardata[varID].resize(numBins * gridsize, 0);
      varcount[varID].resize(numBins * gridsize, 0);
      vartcount[varID].resize(gridsize, 0);
    }

    Varray<double> array(varList1.gridsizeMax());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID++);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);

        auto gridsize = varList1.vars[varID].gridsize;
        auto missval = varList1.vars[varID].missval;

        numMissVals = 0;
        for (size_t i = 0; i < gridsize; ++i)
        {
          if (!fp_is_equal(array[i], missval))
          {
            vartcount[varID][i] += 1;
            int index = 0;
            while (index < numBins)
            {
              auto offset = gridsize * index;
              if (!fp_is_equal(vardata[varID][offset + i], missval) && array[i] >= fltarr[index] && array[i] < fltarr[index + 1])
              {
                vardata[varID][offset + i] += array[i];
                varcount[varID][offset + i] += 1;
                break;
              }
              index++;
            }
          }
          else
          {
            numMissVals++;  // missing value
          }
        }
      }
    }

    cdo_def_timestep(streamID2, 0);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = varList2.vars[varID].gridsize;
      auto missval = varList2.vars[varID].missval;

      // fix missing values

      for (int index = 0; index < numBins; ++index)
      {
        size_t numMissVals = 0;
        auto offset = gridsize * index;

        for (size_t i = 0; i < gridsize; ++i)
        {
          if (vartcount[varID][i] > 0)
          {
            if (operatorID == HISTMEAN || operatorID == HISTFREQ)
            {
              if (varcount[varID][offset + i] > 0)
              {
                if (operatorID == HISTMEAN)
                  vardata[varID][offset + i] /= varcount[varID][offset + i];
                else
                  vardata[varID][offset + i] = varcount[varID][offset + i] / vartcount[varID][i];
              }
            }
          }
          else
          {
            numMissVals++;
            varcount[varID][offset + i] = missval;
            vardata[varID][offset + i] = missval;
          }
        }

        cdo_def_field(streamID2, varID, index);

        if (operatorID == HISTCOUNT)
          cdo_write_field(streamID2, &varcount[varID][offset], numMissVals);
        else
          cdo_write_field(streamID2, &vardata[varID][offset], numMissVals);
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
