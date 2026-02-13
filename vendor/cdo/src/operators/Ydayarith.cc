/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ydayarith  ydayadd         Add multi-year daily time series
      Ydayarith  ydaysub         Subtract multi-year daily time series
      Ydayarith  ydaymul         Multiply multi-year daily time series
      Ydayarith  ydaydiv         Divide multi-year daily time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

class Ydayarith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ydayarith",
    .operators = { { "ydayadd", FieldFunc_Add, 0, YdayarithHelp },
                   { "ydaysub", FieldFunc_Sub, 0, YdayarithHelp },
                   { "ydaymul", FieldFunc_Mul, 0, YdayarithHelp },
                   { "ydaydiv", FieldFunc_Div, 0, YdayarithHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Ydayarith> registration = RegisterEntry<Ydayarith>();

private:
  VarList varList1;
  VarList varList2;

  int operfunc{};

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{ CDI_UNDEFID };

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
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
    constexpr int MaxDays = 373;  //~31*12
    FieldVector2D varsData2[MaxDays];
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID2);
      auto dayOfYear = decode_day_of_year(vDateTime.date);
      // assert(dayOfYear < 1 || dayOfYear >= MaxDays);
      if (dayOfYear == 0)
      {
        cdo_error("Day of year %d out of range (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
        return;
      }
      if (varsData2[dayOfYear].size() > 0)
      {
        cdo_error("Day of year index %d already allocated (date=%s)! Each day of year must only exist once", dayOfYear,
                  date_to_string(vDateTime.date));
        return;
      }

      field2D_init(varsData2[dayOfYear], varList2, FIELD_VEC | FIELD_NAT);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        cdo_read_field(streamID2, varsData2[dayOfYear][varID][levelID]);
      }

      tsID++;
    }

    tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto dayOfYear = decode_day_of_year(vDateTime.date);
      // assert(dayOfYear < 1 || dayOfYear >= MaxDays);
      if (dayOfYear == 0)
      {
        cdo_error("Day of year %d out of range (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
        return;
      }
      if (varsData2[dayOfYear].size() == 0)
      {
        cdo_error("Day of year index %d not found (date=%s)!", dayOfYear, date_to_string(vDateTime.date));
        return;
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        field2_function(field, varsData2[dayOfYear][varID][levelID], operfunc);

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, field);
      }
      tsID++;
    }
    return;
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);
  }
};
