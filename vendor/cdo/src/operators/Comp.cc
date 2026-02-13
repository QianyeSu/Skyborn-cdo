/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Comp       eq              Equal
      Comp       ne              Not equal
      Comp       le              Less equal
      Comp       lt              Less than
      Comp       ge              Greater equal
      Comp       gt              Greater than
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_fill.h"
#include "field_functions.h"

static auto func_comp
    = [](auto hasMissvals, auto n, auto mv1, auto mv2, auto const &vIn1, auto const &vIn2, auto &vOut, auto binary_operator)
{
  if (hasMissvals)
  {
    if (std::isnan(mv1) || std::isnan(mv2))
      for (size_t i = 0; i < n; ++i)
        vOut[i] = (fp_is_equal(vIn1[i], mv1) || fp_is_equal(vIn2[i], mv2)) ? mv1 : binary_operator(vIn1[i], vIn2[i]);
    else
      for (size_t i = 0; i < n; ++i)
        vOut[i] = (is_equal(vIn1[i], mv1) || is_equal(vIn2[i], mv2)) ? mv1 : binary_operator(vIn1[i], vIn2[i]);
  }
  else
  {
    for (size_t i = 0; i < n; ++i) vOut[i] = binary_operator(vIn1[i], vIn2[i]);
  }
};

class Comp : public Process
{
  enum struct FillType
  {
    NONE,
    TS,
    REC,
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Comp",
    .operators = { { "eq", FieldFunc_EQ, 0, CompHelp },
                   { "ne", FieldFunc_NE, 0, CompHelp },
                   { "le", FieldFunc_LE, 0, CompHelp },
                   { "lt", FieldFunc_LT, 0, CompHelp },
                   { "ge", FieldFunc_GE, 0, CompHelp },
                   { "gt", FieldFunc_GT, 0, CompHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Comp> registration = RegisterEntry<Comp>();

  FillType fillType{ FillType::NONE };
  Varray2D<double> vardata;

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  int taxisID3{};

  int operFunc{};

  double *arrayx1{};
  double *arrayx2{};

  VarList varList1;
  VarList varList2;

  Varray<double> vaIn1, vaIn2, vaOut;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operFunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    auto taxisID2 = vlistInqTaxis(vlistID2);

    auto ntsteps1 = vlistNtsteps(vlistID1);
    auto ntsteps2 = vlistNtsteps(vlistID2);
    if (ntsteps1 == 0) ntsteps1 = 1;
    if (ntsteps2 == 0) ntsteps2 = 1;

    auto fillstream1 = false;

    if (vlistNumFields(vlistID1) != 1 && vlistNumFields(vlistID2) == 1)
    {
      fillType = FillType::REC;
      cdo_print("Filling up stream2 >%s< by copying the first field.", cdo_get_stream_name(1));
      if (ntsteps2 != 1) cdo_abort("stream2 has more than 1 timestep!");
    }
    else if (vlistNumFields(vlistID1) == 1 && vlistNumFields(vlistID2) != 1)
    {
      fillType = FillType::REC;
      cdo_print("Filling up stream1 >%s< by copying the first field.", cdo_get_stream_name(0));
      if (ntsteps1 != 1) cdo_abort("stream1 has more than 1 timestep!");
      fillstream1 = true;
      std::swap(streamID1, streamID2);
      std::swap(vlistID1, vlistID2);
      std::swap(taxisID1, taxisID2);
    }

    if (fillType == FillType::NONE) vlist_compare(vlistID1, vlistID2, CmpVarList::All);

    if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

    if (fillType == FillType::NONE)
    {
      if (ntsteps1 != 1 && ntsteps2 == 1)
      {
        fillType = FillType::TS;
        cdo_print("Filling up stream2 >%s< by copying the first timestep.", cdo_get_stream_name(1));
      }
      else if (ntsteps1 == 1 && ntsteps2 != 1)
      {
        fillType = FillType::TS;
        cdo_print("Filling up stream1 >%s< by copying the first timestep.", cdo_get_stream_name(0));
        fillstream1 = true;
        std::swap(streamID1, streamID2);
        std::swap(vlistID1, vlistID2);
        std::swap(taxisID1, taxisID2);
      }

      if (fillType == FillType::TS) cdo_fill_ts(vlistID2, vardata);
    }

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    vaIn1.resize(varList1.gridsizeMax());
    vaIn2.resize(varList1.gridsizeMax());
    vaOut.resize(varList1.gridsizeMax());

    arrayx1 = vaIn1.data();
    arrayx2 = vaIn2.data();
    if (fillstream1) { std::swap(arrayx1, arrayx2); }

    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID3);

    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (tsID == 0 || fillType == FillType::NONE)
      {
        auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
        if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        size_t numMissVals1 = 0, numMissVals2 = 0;
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, arrayx1, &numMissVals1);

        if (tsID == 0 || fillType == FillType::NONE)
        {
          if (fieldID == 0 || fillType != FillType::REC)
          {
            auto [varIDx, levelIDx] = cdo_inq_field(streamID2);
            cdo_read_field(streamID2, arrayx2, &numMissVals2);
            varID = varIDx;
            levelID = levelIDx;
          }

          if (fillType == FillType::TS)
          {
            auto gridsize = varList2.vars[varID].gridsize;
            auto offset = gridsize * levelID;
            array_copy(gridsize, arrayx2, &vardata[varID][offset]);
          }
        }
        else if (fillType == FillType::TS)
        {
          auto gridsize = varList2.vars[varID].gridsize;
          auto offset = gridsize * levelID;
          array_copy(gridsize, &vardata[varID][offset], arrayx2);
        }

        auto const &var1 = varList1.vars[varID];
        auto datatype1 = var1.dataType;
        auto gridsize1 = var1.gridsize;
        auto missval1 = var1.missval;

        auto xvarID = (fillType == FillType::REC) ? 0 : varID;
        auto const &var2 = varList2.vars[xvarID];
        auto datatype2 = var2.dataType;
        auto gridsize2 = var2.gridsize;
        auto missval2 = var2.missval;

        if (gridsize1 != gridsize2)
          cdo_abort("Streams have different gridsize (gridsize1 = %zu; gridsize2 = %zu)!", gridsize1, gridsize2);

        auto ngp = gridsize1;

        if (datatype1 != datatype2)
        {
          if (datatype1 == CDI_DATATYPE_FLT32 && datatype2 == CDI_DATATYPE_FLT64)
          {
            missval2 = (float) missval2;
            for (size_t i = 0; i < ngp; ++i) vaIn2[i] = (float) vaIn2[i];
          }
          else if (datatype1 == CDI_DATATYPE_FLT64 && datatype2 == CDI_DATATYPE_FLT32)
          {
            missval1 = (float) missval1;
            for (size_t i = 0; i < ngp; ++i) vaIn1[i] = (float) vaIn1[i];
          }
        }

        if (numMissVals1 > 0) cdo_check_missval(missval1, varList1.vars[varID].name);
        // if (numMissVals2 > 0) cdo_check_missval(missval2, varList2.vars[varID].name);

        auto hasMissvals = (numMissVals1 > 0 || numMissVals2 > 0);
        // clang-format off
        if      (operFunc == FieldFunc_EQ) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_EQ);
        else if (operFunc == FieldFunc_NE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_NE);
        else if (operFunc == FieldFunc_LE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_LE);
        else if (operFunc == FieldFunc_LT) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_LT);
        else if (operFunc == FieldFunc_GE) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_GE);
        else if (operFunc == FieldFunc_GT) func_comp(hasMissvals, ngp, missval1, missval2, vaIn1, vaIn2, vaOut, binary_op_GT);
        else cdo_abort("Operator not implemented!");
        // clang-format on

        auto numMissValsOut = varray_num_mv(ngp, vaOut, missval1);
        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, vaOut.data(), numMissValsOut);
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
