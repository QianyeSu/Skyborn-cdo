/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yhourarith  yhouradd         Add multi-year hourly time series
      Yhourarith  yhoursub         Subtract multi-year hourly time series
      Yhourarith  yhourmul         Multiply multi-year hourly time series
      Yhourarith  yhourdiv         Divide multi-year hourly time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "field_functions.h"

constexpr int MaxHours = 9301;  // 31*12*25 + 1

class Yhourarith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Yhourarith",
    .operators = { { "yhouradd", FieldFunc_Add, 0, YhourarithHelp },
                   { "yhoursub", FieldFunc_Sub, 0, YhourarithHelp },
                   { "yhourmul", FieldFunc_Mul, 0, YhourarithHelp },
                   { "yhourdiv", FieldFunc_Div, 0, YhourarithHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Yhourarith> registration = RegisterEntry<Yhourarith>();

  int operfunc;
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;

  int vlistID2{ CDI_UNDEFID };

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
    vlistID2 = cdo_stream_inq_vlist(streamID2);
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
  }

  void
  run() override
  {
    Field field;
    FieldVector2D varsData2[MaxHours];

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      auto hourOfYear = decode_hour_of_year(taxisInqVdatetime(taxisID2), MaxHours);
      if (varsData2[hourOfYear].size() > 0) cdo_abort("Hour of year index %d already allocated!", hourOfYear);

      field2D_init(varsData2[hourOfYear], varList2, FIELD_VEC | FIELD_NAT);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        cdo_read_field(streamID2, varsData2[hourOfYear][varID][levelID]);
      }

      tsID++;
    }

    cdo_stream_close(streamID2);

    tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto hourOfYear = decode_hour_of_year(taxisInqVdatetime(taxisID1), MaxHours);
      if (varsData2[hourOfYear].size() == 0) cdo_abort("Hour of year index %d not found!", hourOfYear);

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        field2_function(field, varsData2[hourOfYear][varID][levelID], operfunc);

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
    cdo_stream_close(streamID1);
  }
};
