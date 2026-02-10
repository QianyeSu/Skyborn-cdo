/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Cond       ifthen          If then
      Cond       ifnotthen       If not then
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"

static void
operator_IFTHEN(size_t n, double mv1, double mv2, Varray<double> const &vIn1, Varray<double> const &vIn2, Varray<double> &vOut)
{
  if (std::isnan(mv1))
    for (size_t i = 0; i < n; ++i) vOut[i] = (fp_is_not_equal(vIn1[i], mv1) && fp_is_not_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn1[i], mv1) && !is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
}

static void
operator_IFNOTTHEN(size_t n, double mv1, double mv2, Varray<double> const &vIn1, Varray<double> const &vIn2, Varray<double> &vOut)
{
  if (std::isnan(mv1))
    for (size_t i = 0; i < n; ++i) vOut[i] = (fp_is_not_equal(vIn1[i], mv1) && fp_is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
  else
    for (size_t i = 0; i < n; ++i) vOut[i] = (!is_equal(vIn1[i], mv1) && is_equal(vIn1[i], 0.0)) ? vIn2[i] : mv2;
}

class Cond : public Process
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
    .name = "Cond",
    .operators = { { "ifthen", CondHelp }, { "ifnotthen", CondHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Cond> registration = RegisterEntry<Cond>(module);

private:
  int IFTHEN{}, IFNOTTHEN{};
  int filltype = FILL_NONE;
  double missval1 = -9.E33;
  Varray2D<size_t> varnumMissVals1;
  Varray2D<double> vardata1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};

  int operatorID{};

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    IFTHEN = module.get_id("ifthen");
    IFNOTTHEN = module.get_id("ifnotthen");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID2);

    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID2);
    vlistDefTaxis(vlistID3, taxisID3);

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

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    if (filltype == FILL_REC && varList2.gridsizeMax() != gridInqSize(vlistGrid(vlistID1, 0)))
      cdo_abort("Stream1 >%s< has wrong gridsize!", cdo_get_stream_name(0));

    if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

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
    Varray<double> varrayIn1(varList2.gridsizeMax());
    Varray<double> varrayIn2(varList2.gridsizeMax());
    Varray<double> varrayOut(varList2.gridsizeMax());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      if (tsID == 0 || filltype == FILL_NONE)
      {
        auto numFields2 = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID2);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        size_t numMissVals1 = 0, numMissVals2;
        cdo_read_field(streamID2, &varrayIn2[0], &numMissVals2);

        if (tsID == 0 || filltype == FILL_NONE)
        {
          if (fieldID == 0 || filltype != FILL_REC)
          {
            auto [varIDx, levelIDx] = cdo_inq_field(streamID1);
            cdo_read_field(streamID1, &varrayIn1[0], &numMissVals1);
            varID = varIDx;
            levelID = levelIDx;
          }

          if (filltype == FILL_TS)
          {
            auto gridsize = varList1.vars[varID].gridsize;
            auto offset = gridsize * levelID;
            array_copy(gridsize, &varrayIn1[0], &vardata1[varID][offset]);
            varnumMissVals1[varID][levelID] = numMissVals1;
          }
        }
        else if (filltype == FILL_TS)
        {
          auto gridsize = varList1.vars[varID].gridsize;
          auto offset = gridsize * levelID;
          array_copy(gridsize, &vardata1[varID][offset], &varrayIn1[0]);
          numMissVals1 = varnumMissVals1[varID][levelID];
        }

        auto const &var1 = varList1.vars[varID];
        auto const &var2 = varList2.vars[varID];
        auto ngp = var2.gridsize;
        auto missval2 = var2.missval;
        if (fieldID == 0 || filltype != FILL_REC) missval1 = var1.missval;

        if (numMissVals1 > 0) cdo_check_missval(missval1, var1.name);

        // clang-format off
        if      (operatorID == IFTHEN)    operator_IFTHEN(ngp, missval1, missval2, varrayIn1, varrayIn2, varrayOut);
        else if (operatorID == IFNOTTHEN) operator_IFNOTTHEN(ngp, missval1, missval2, varrayIn1, varrayIn2, varrayOut);
        else cdo_abort("Operator not implemented!");
        // clang-format on

        auto numMissVals3 = varray_num_mv(ngp, varrayOut, missval2);
        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, varrayOut.data(), numMissVals3);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
