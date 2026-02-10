/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yeararith  yearadd         Add yearly time series
      Yeararith  yearsub         Subtract yearly time series
      Yeararith  yearmul         Multiply yearly time series
      Yeararith  yeardiv         Divide yearly time series
*/

#include <cdi.h>

#include <climits>

#include "cdo_vlist.h"
#include "process_int.h"
#include "field_functions.h"

class Yeararith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Yeararith",
    .operators = { { "yearadd", FieldFunc_Add, 0, YeararithHelp },
                   { "yearsub", FieldFunc_Sub, 0, YeararithHelp },
                   { "yearmul", FieldFunc_Mul, 0, YeararithHelp },
                   { "yeardiv", FieldFunc_Div, 0, YeararithHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Yeararith> registration = RegisterEntry<Yeararith>(module);

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;

  FieldVector2D varsData2;

  int operfunc;
  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID3);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    field2D_init(varsData2, varList2, FIELD_VEC | FIELD_NAT);
  }

  void
  run() override
  {
    Field field;

    int year0 = -INT_MAX + 1;
    int year2last = 0;
    int tsID2 = 0;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto year = taxisInqVdatetime(taxisID1).date.year;
      if (year > year0)
      {
        auto lfound = false;
        while (true)
        {
          auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID2);
          if (numFields2 == 0) break;

          tsID2++;
          auto year2 = taxisInqVdatetime(taxisID2).date.year;
          if (year == year2)
          {
            lfound = true;
            year0 = year;
            while (numFields2--)
            {
              auto [varID, levelID] = cdo_inq_field(streamID2);
              cdo_read_field(streamID2, varsData2[varID][levelID]);
            }
            break;
          }

          if (tsID2 > 1 && year2 <= year2last) cdo_abort("stream2 doesn't contain yearly data!");
          year2last = year2;
        }

        if (!lfound) cdo_abort("Data of year %d not found in stream2!", year);
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        field2_function(field, varsData2[varID][levelID], operfunc);

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, field);
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
