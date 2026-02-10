/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
      Ymonarith  yseasadd        Add multi-year seasonal time series
      Ymonarith  yseassub        Subtract multi-year seasonal time series
      Ymonarith  yseasmul        Multiply multi-year seasonal time series
      Ymonarith  yseasdiv        Divide multi-year seasonal time series
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "cdo_season.h"
#include "process_int.h"
#include "field_functions.h"

constexpr int MaxMonths = 12;

static const char *monthNames[]
    = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

static int
get_month_index(const CdiDate vDate)
{
  int month = vDate.month;
  if (month < 1 || month > MaxMonths) cdo_abort("Month %d out of range!", month);
  return --month;
}

enum
{
  MONTHLY,
  SEASONAL
};

class Ymonarith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ymonarith",
    .operators = { { "ymonadd", FieldFunc_Add, MONTHLY, YmonarithHelp },
                   { "ymonsub", FieldFunc_Sub, MONTHLY, YmonarithHelp },
                   { "ymonmul", FieldFunc_Mul, MONTHLY, YmonarithHelp },
                   { "ymondiv", FieldFunc_Div, MONTHLY, YmonarithHelp },
                   { "yseasadd", FieldFunc_Add, SEASONAL, YseasarithHelp },
                   { "yseassub", FieldFunc_Sub, SEASONAL, YseasarithHelp },
                   { "yseasmul", FieldFunc_Mul, SEASONAL, YseasarithHelp },
                   { "yseasdiv", FieldFunc_Div, SEASONAL, YseasarithHelp } },
    .aliases = { { "anomaly", "ymonsub" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Ymonarith> registration = RegisterEntry<Ymonarith>(module);

  int operfunc;
  int opertype;

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;

  int vlistID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;
  Field field;
  FieldVector2D varsData2[MaxMonths];

  void
  process_data(int tsID, int numFields, int mon)
  {
    cdo_taxis_copy_timestep(taxisID3, taxisID1);
    cdo_def_timestep(streamID3, tsID);

    while (numFields--)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      field.init(varList1.vars[varID]);
      cdo_read_field(streamID1, field);

      field2_function(field, varsData2[mon][varID][levelID], operfunc);

      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, field);
    }
  }

  void
  load_data_from_stream(int mon, int numFields)
  {
    field2D_init(varsData2[mon], varList2, FIELD_VEC | FIELD_NAT);

    while (numFields--)
    {
      auto [varID, levelID] = cdo_inq_field(streamID2);
      cdo_read_field(streamID2, varsData2[mon][varID][levelID]);
    }
  }

  void
  run_monthly()
  {
    for (int numFields, tsID = 0; (numFields = cdo_stream_inq_timestep(streamID2, tsID)) != 0; tsID++)
    {
      auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
      if (varsData2[mon].size())
        cdo_abort("%s already allocated! The second input file must contain monthly mean values for a maximum of one year.",
                  monthNames[mon]);

      load_data_from_stream(mon, numFields);
    }

    for (int numFields, tsID = 0; (numFields = cdo_stream_inq_timestep(streamID1, tsID)) != 0; tsID++)
    {
      auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
      if (varsData2[mon].size() == 0)
        cdo_abort("%s not found! The second input file must contain monthly mean values for a maximum of one year.",
                  monthNames[mon]);

      process_data(tsID, numFields, mon);
    }
  }

  void
  run_seasonal()
  {
    auto seasonNames = get_season_name();
    for (int numFields, tsID = 0; (numFields = cdo_stream_inq_timestep(streamID2, tsID)) != 0; tsID++)
    {
      auto mon = get_month_index(taxisInqVdatetime(taxisID2).date);
      auto season = month_to_season(mon + 1);
      if (varsData2[season].size()) cdo_abort("Season %s already allocated!", seasonNames[season]);

      load_data_from_stream(season, numFields);
    }

    for (int numFields, tsID = 0; (numFields = cdo_stream_inq_timestep(streamID1, tsID)) != 0; tsID++)
    {
      auto mon = get_month_index(taxisInqVdatetime(taxisID1).date);
      auto season = month_to_season(mon + 1);
      if (varsData2[season].size() == 0) cdo_abort("Season %s not found!", seasonNames[season]);

      process_data(tsID, numFields, season);
    }
  }

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    opertype = cdo_operator_f2(operatorID);
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
    if (opertype == SEASONAL)
      run_seasonal();
    else
      run_monthly();
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);
  }
};
