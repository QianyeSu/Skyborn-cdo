/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arith      add             Add two fields
      Arith      sub             Subtract two fields
      Arith      mul             Multiply two fields
      Arith      div             Divide two fields
      Arith      min             Minimum of two fields
      Arith      max             Maximum of two fields
      Arith      atan2           Arc tangent of two fields
*/

#include <cdi.h>

#include <utility>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_options.h"
#include "cdo_fill.h"
#include "cdo_cdi_wrapper.h"
#include "field_functions.h"

class Arith : public Process
{
  enum struct FillType
  {
    NONE,
    TS,
    VAR,
    VARTS,
    FULL_FILE
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Arith",
    .operators = { { "add", FieldFunc_Add, 1, ArithHelp },
                   { "sub", FieldFunc_Sub, 1, ArithHelp },
                   { "mul", FieldFunc_Mul, 1, ArithHelp },
                   { "div", FieldFunc_Div, 1, ArithHelp },
                   { "min", FieldFunc_Min, 0, ArithHelp },
                   { "max", FieldFunc_Max, 0, ArithHelp },
                   { "atan2", FieldFunc_Atan2, 0, ArithHelp },
                   { "setmiss", FieldFunc_Setmiss, 0, ArithHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Arith> registration = RegisterEntry<Arith>();

  FillType fillType{ FillType::NONE };
  int nlevels2 = 1;
  int levelID2 = -1;
  Varray2D<size_t> varnumMissVals;
  Varray2D<double> vardata;
  std::vector<size_t> varnumMissVals2;
  Varray<double> vardata2;

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  CdoStreamID streamID2x;

  CdoStreamID streamID3;
  int vlistID3{};
  int taxisID3{};

  int operfunc{};

  size_t nwpv{};

  bool fillStream1 = false, fillStream2 = false, fillStream1x = false;

  Field field1, field2;
  Field *fieldx1{};
  Field *fieldx2{};

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    bool opercplx = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    streamID2x = streamID2;

    fieldx1 = &field1;
    fieldx2 = &field2;

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID1x = vlistID1;
    auto vlistID2x = vlistID2;

    if (Options::cdoVerbose) vlistPrint(vlistID1);
    if (Options::cdoVerbose) vlistPrint(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    auto taxisID2 = vlistInqTaxis(vlistID2);

    auto ntsteps1 = vlistNtsteps(vlistID1);
    auto ntsteps2 = vlistNtsteps(vlistID2);
    if (ntsteps1 == 0) ntsteps1 = 1;
    if (ntsteps2 == 0) ntsteps2 = 1;

    auto nvars1 = vlistNvars(vlistID1);
    auto nvars2 = vlistNvars(vlistID2);

    if (nvars1 == 1 && nvars2 == 1)
    {
      fillStream2 = (vlistNumFields(vlistID1) != 1 && vlistNumFields(vlistID2) == 1);
      fillStream1 = (vlistNumFields(vlistID1) == 1 && vlistNumFields(vlistID2) != 1);
      if (fillStream1 && ntsteps1 != 1 && ntsteps2 == 1)
      {
        fillStream1 = false;
        fillStream2 = true;
        fillStream1x = true;
      }
    }
    else
    {
      fillStream2 = (nvars1 != 1 && nvars2 == 1);
      fillStream1 = (nvars1 == 1 && nvars2 != 1);
    }

    if (fillStream1x)
    {
      nlevels2 = vlist_compare_x(vlistID2, vlistID1, CmpVarList::Dim);

      fillType = FillType::NONE;
      cdo_print("Filling up stream1 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(0));
    }
    else if (fillStream2)
    {
      nlevels2 = vlist_compare_x(vlistID1, vlistID2, CmpVarList::Dim);

      if (ntsteps1 != 1 && ntsteps2 == 1)
      {
        fillType = FillType::VAR;
        cdo_print("Filling up stream2 >%s< by copying the first variable.", cdo_get_stream_name(1));
      }
      else
      {
        fillType = FillType::VARTS;
        cdo_print("Filling up stream2 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(1));
      }
    }
    else if (fillStream1)
    {
      nlevels2 = vlist_compare_x(vlistID2, vlistID1, CmpVarList::Dim);

      if (ntsteps1 == 1 && ntsteps2 != 1)
      {
        fillType = FillType::VAR;
        cdo_print("Filling up stream1 >%s< by copying the first variable.", cdo_get_stream_name(0));
      }
      else
      {
        fillType = FillType::VARTS;
        cdo_print("Filling up stream1 >%s< by copying the first variable of each timestep.", cdo_get_stream_name(0));
      }

      std::swap(streamID1, streamID2);
      vlistID1x = vlistID2;
      vlistID2x = vlistID1;
      std::swap(taxisID1, taxisID2);
      fieldx1 = &field2;
      fieldx2 = &field1;
    }

    if (fillStream1x == false && fillType == FillType::NONE) vlist_compare(vlistID1, vlistID2, CmpVarList::All);

    nwpv = (vlistNumber(vlistID1) == CDI_COMP && vlistNumber(vlistID2) == CDI_COMP) ? 2 : 1;
    if (nwpv == 2 && !opercplx) cdo_abort("Fields with complex numbers are not supported by this operator!");
    auto gridsizeMax = nwpv * cdo_vlist_gridsizemax(vlistID1x);

    field1.resize(gridsizeMax);
    field2.resize(gridsizeMax);

    if (fillStream1x || fillType == FillType::VAR || fillType == FillType::VARTS)
    {
      vardata2.resize(gridsizeMax * nlevels2);
      varnumMissVals2.resize(nlevels2);
    }

    if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

    bool streamsSwaped = false;
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

        streamsSwaped = true;
        vlistDefNtsteps(vlistID1, ntsteps2);
        std::swap(streamID1, streamID2);
        vlistID1x = vlistID2;
        vlistID2x = vlistID1;
        std::swap(taxisID1, taxisID2);
        fieldx1 = &field2;
        fieldx2 = &field1;
      }

      if (fillType == FillType::TS) cdo_fill_ts(vlistID2x, vardata, varnumMissVals);
    }

    vlistID3 = vlistDuplicate(streamsSwaped ? vlistID1 : vlistID1x);
    vlist_unpack(vlistID3);
    if (streamsSwaped)
    {
      auto nvars = vlistNvars(vlistID1);
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarMissval(vlistID3, varID, vlistInqVarMissval(vlistID1, varID));
      for (int varID = 0; varID < nvars; ++varID) vlistDefVarTimetype(vlistID3, varID, vlistInqVarTimetype(vlistID1x, varID));
    }

    if (fillStream1x)
    {
      auto zaxisID2 = vlistZaxis(vlistID2x, 0);
      vlistChangeZaxisIndex(vlistID3, 0, zaxisID2);
    }

    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    varList1 = VarList(vlistID1x);
    varList2 = VarList(vlistID2x);
  }

  void
  run() override
  {
    int numFields1;
    int numFields2;
    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      numFields1 = cdo_stream_inq_timestep(streamID1, tsID);

      numFields2 = 0;
      if (tsID == 0 || fillType == FillType::NONE || fillType == FillType::FULL_FILE || fillType == FillType::VARTS)
      {
        numFields2 = cdo_stream_inq_timestep(streamID2, tsID2);
        if (numFields2 == 0)
        {
          if (numFields1 == 0) break;

          if (fillType == FillType::NONE && streamID2x == streamID2)
          {
            fillType = FillType::FULL_FILE;
            cdo_print("Filling up stream2 >%s< by copying all timesteps.", cdo_get_stream_name(1));
          }

          if (fillType == FillType::FULL_FILE)
          {
            cdo_stream_close(streamID2);

            if (stream_is_pipe(1)) cdo_abort("infile2 cannot be a pipe in fill mode!");

            streamID2 = cdo_open_read(1);
            streamID2x = streamID2;

            (void) cdo_stream_inq_vlist(streamID2);

            tsID2 = 0;
            numFields2 = cdo_stream_inq_timestep(streamID2, tsID2);
            if (numFields2 == 0) cdo_abort("Empty input stream %s!", cdo_get_stream_name(1));
          }
          else
            cdo_abort("Input streams have different number of timesteps!");
        }
      }

      if (numFields1 == 0 || (numFields2 == 0 && fillType != FillType::TS && fillType != FillType::VAR)) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      auto numFields = fillStream1x ? numFields2 : numFields1;
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto lread1 = true;
        if (fillStream1x && fieldID > 0) lread1 = false;
        int varID = -1, levelID;
        if (lread1)
        {
          std::tie(varID, levelID) = cdo_inq_field(streamID1);
          cdo_read_field(streamID1, *fieldx1);

          if (fillStream1x)
          {
            auto gridsize = nwpv * varList1.vars[varID].gridsize;
            array_copy(gridsize, fieldx1->vec_d.data(), &vardata2[0]);
            varnumMissVals2[0] = fieldx1->numMissVals;
          }
        }

        if (fillStream1x) levelID = fieldID;

        auto varID2 = varID;

        if (tsID == 0 || fillType == FillType::NONE || fillType == FillType::FULL_FILE || fillType == FillType::VARTS)
        {
          auto lstatus = (nlevels2 > 1) ? (varID == 0) : (fieldID == 0);
          if (lstatus || (fillType != FillType::VAR && fillType != FillType::VARTS))
          {
            std::tie(varID2, levelID2) = cdo_inq_field(streamID2);
            cdo_read_field(streamID2, *fieldx2);
            if (varID != varID2) cdo_abort("Internal error, varIDs of input streams differ!");
            if (fillStream1x == false && levelID != levelID2) cdo_abort("Internal error, levelIDs of input streams differ!");
          }

          if (fillType == FillType::TS)
          {
            auto gridsize = nwpv * varList2.vars[varID].gridsize;
            auto offset = gridsize * levelID;
            array_copy(gridsize, fieldx2->vec_d.data(), &vardata[varID][offset]);
            varnumMissVals[varID][levelID] = fieldx2->numMissVals;
          }
          else if (lstatus && (fillType == FillType::VAR || fillType == FillType::VARTS))
          {
            auto gridsize = nwpv * varList2.vars[0].gridsize;
            auto offset = gridsize * levelID2;
            array_copy(gridsize, fieldx2->vec_d.data(), &vardata2[offset]);
            varnumMissVals2[levelID2] = fieldx2->numMissVals;
          }
        }
        else if (fillType == FillType::TS)
        {
          auto gridsize = nwpv * varList2.vars[varID2].gridsize;
          auto offset = gridsize * levelID;
          array_copy(gridsize, &vardata[varID][offset], fieldx2->vec_d.data());
          fieldx2->numMissVals = varnumMissVals[varID][levelID];
        }

        if (fillStream1x)
        {
          auto gridsize = nwpv * varList1.vars[0].gridsize;
          array_copy(gridsize, &vardata2[0], fieldx1->vec_d.data());
          fieldx1->numMissVals = varnumMissVals2[0];
        }

        fieldx1->grid = varList1.vars[varID].gridID;
        fieldx1->missval = varList1.vars[varID].missval;
        fieldx1->nwpv = varList1.vars[varID].nwpv;

        if (fillType == FillType::VAR || fillType == FillType::VARTS)
        {
          levelID2 = (nlevels2 > 1) ? levelID : 0;
          auto gridsize = nwpv * varList2.vars[0].gridsize;
          auto offset = gridsize * levelID2;
          array_copy(gridsize, &vardata2[offset], fieldx2->vec_d.data());
          fieldx2->numMissVals = varnumMissVals2[levelID2];
          fieldx2->grid = varList2.vars[0].gridID;
          fieldx2->missval = varList2.vars[0].missval;
          fieldx2->nwpv = varList2.vars[0].nwpv;
        }
        else
        {
          fieldx2->grid = varList2.vars[varID2].gridID;
          fieldx2->missval = varList2.vars[varID2].missval;
          fieldx2->nwpv = varList2.vars[varID2].nwpv;
        }

        auto field2_func = (nwpv == 2) ? field2_function_complex : field2_function;
        field2_func(field1, field2, operfunc);

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, field1);
      }

      tsID++;
      tsID2++;
    }

    if (numFields1 == 0 && numFields2 > 0) cdo_warning("stream2 has more time steps than stream1!");
    // if (numFields > 0 && numFields2 == 0) cdo_warning("stream1 has more time steps than stream2!");
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID3);
  }
};
