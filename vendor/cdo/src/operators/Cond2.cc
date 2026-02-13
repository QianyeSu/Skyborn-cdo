/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Cond2      ifthenelse      If then else
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"

class Cond2 : public Process
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Cond2",
    .operators = { { "ifthenelse", Cond2Help } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Cond2> registration = RegisterEntry<Cond2>();

private:
  int filltype = FILL_NONE;
  double missval1 = -9.E33;
  size_t numMissVals1 = 0;
  Varray2D<size_t> varnumMissVals1;
  Varray2D<double> vardata1;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  CdoStreamID streamID4;
  int taxisID2{ CDI_UNDEFID };
  int taxisID4{};

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = cdo_stream_inq_vlist(streamID3);
    auto vlistID4 = vlistDuplicate(vlistID2);

    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID4 = taxisDuplicate(taxisID2);
    vlistDefTaxis(vlistID4, taxisID4);

    auto ntsteps1 = vlistNtsteps(vlistID1);
    auto ntsteps2 = vlistNtsteps(vlistID2);
    if (ntsteps1 == 0) ntsteps1 = 1;
    if (ntsteps2 == 0) ntsteps2 = 1;

    if (vlistNumFields(vlistID1) == 1 && vlistNumFields(vlistID2) != 1)
    {
      filltype = FILL_REC;
      cdo_print("Filling up stream1 >%s< by copying the first field.", cdo_get_stream_name(0));
    }

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    if (filltype == FILL_NONE) varList_compare(varList1, varList2, CmpVarList::Dim);

    varList_compare(varList2, VarList(vlistID3), CmpVarList::Dim);

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);

    auto gridsizeMax = varList1.gridsizeMax();
    if (filltype == FILL_REC && gridsizeMax != gridInqSize(vlistGrid(vlistID1, 0)))
      cdo_abort("Stream1 >%s< has wrong gridsize!", cdo_get_stream_name(0));

    if (Options::cdoVerbose)
      cdo_print("Number of timesteps: file1 %d, file2 %d, file3 %d", ntsteps1, ntsteps2, vlistNtsteps(vlistID3));

    if (filltype == FILL_NONE)
    {
      if (ntsteps1 == 1 && ntsteps2 != 1)
      {
        filltype = FILL_TS;
        cdo_print("Filling up stream1 >%s< by copying the first timestep.", cdo_get_stream_name(0));

        cdo_fill_ts(vlistID1, vardata1, varnumMissVals1);
      }
    }
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array1(gridsizeMax);
    Varray<double> array2(gridsizeMax);
    Varray<double> array3(gridsizeMax);
    Varray<double> array4(gridsizeMax);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      auto numFields3 = cdo_stream_inq_timestep(streamID3, tsID);
      if (numFields3 == 0) cdo_abort("Input streams have different number of timesteps!");

      if (tsID == 0 || filltype == FILL_NONE)
      {
        auto numFields2 = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");
      }

      cdo_taxis_copy_timestep(taxisID4, taxisID2);
      cdo_def_timestep(streamID4, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        (void) cdo_inq_field(streamID2);
        size_t numMissVals;
        cdo_read_field(streamID2, &array2[0], &numMissVals);

        auto [varID, levelID] = cdo_inq_field(streamID3);
        cdo_read_field(streamID3, &array3[0], &numMissVals);

        if (tsID == 0 || filltype == FILL_NONE)
        {
          if (fieldID == 0 || filltype != FILL_REC)
          {
            auto [varIDx, levelIDx] = cdo_inq_field(streamID1);
            cdo_read_field(streamID1, &array1[0], &numMissVals1);
            varID = varIDx;
            levelID = levelIDx;
          }

          if (filltype == FILL_TS)
          {
            auto gridsize = varList1.vars[varID].gridsize;
            auto offset = gridsize * levelID;
            array_copy(gridsize, &array1[0], &vardata1[varID][offset]);
            varnumMissVals1[varID][levelID] = numMissVals1;
          }
        }
        else if (filltype == FILL_TS)
        {
          auto gridsize = varList1.vars[varID].gridsize;
          auto offset = gridsize * levelID;
          array_copy(gridsize, &vardata1[varID][offset], &array1[0]);
          numMissVals1 = varnumMissVals1[varID][levelID];
        }

        auto const &var1 = varList1.vars[varID];
        auto const &var2 = varList2.vars[varID];
        auto gridsize = var2.gridsize;
        auto missval2 = var2.missval;
        if (fieldID == 0 || filltype != FILL_REC) missval1 = var1.missval;

        if (numMissVals1 > 0) cdo_check_missval(missval1, var1.name);

        for (size_t i = 0; i < gridsize; ++i)
          array4[i] = fp_is_equal(array1[i], missval1) ? missval2 : !fp_is_equal(array1[i], 0.) ? array2[i] : array3[i];

        numMissVals = varray_num_mv(gridsize, array4, missval2);
        cdo_def_field(streamID4, varID, levelID);
        cdo_write_field(streamID4, array4.data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID4);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
