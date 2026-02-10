/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Dayarith  dayadd         Add daily time series
      Dayarith  daysub         Subtract daily time series
      Dayarith  daymul         Multiply daily time series
      Dayarith  daydiv         Divide daily time series
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "field_functions.h"

class Dayarith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Dayarith",
    .operators = { { "dayadd", FieldFunc_Add, 0, DayarithHelp },
                   { "daysub", FieldFunc_Sub, 0, DayarithHelp },
                   { "daymul", FieldFunc_Mul, 0, DayarithHelp },
                   { "daydiv", FieldFunc_Div, 0, DayarithHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Dayarith> registration = RegisterEntry<Dayarith>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  CdoStreamID streamID3;
  int taxisID3{};

  int operfunc{};

  VarList varList1{};
  Field field;
  FieldVector2D varsData2;

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
    VarList varList2(vlistID2);
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
    CdiDate vDate2{};
    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDate1 = taxisInqVdatetime(taxisID1).date;

      if (!cdiDate_isEQ(vDate1, vDate2))
      {
        int year1, month1, day1;
        cdiDate_decode(vDate1, &year1, &month1, &day1);
        if (Options::cdoVerbose) cdo_print("Process: Year=%4d  Month=%2d  Day=%2d", year1, month1, day1);

        auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID2);
        if (numFields2 == 0) cdo_abort("Missing year=%4d mon=%2d day=%2d in %s!", year1, month1, day1, cdo_get_stream_name(1));

        vDate2 = taxisInqVdatetime(taxisID2).date;
        if (!cdiDate_isEQ(vDate1, vDate2))
        {
          int year2, month2, day2;
          cdiDate_decode(vDate2, &year2, &month2, &day2);
          cdo_abort("Timestep %d in %s has wrong date! Current year=%4d mon=%2d day=%2d, expected year=%4d mon=%2d day=%2d",
                    tsID2 + 1, cdo_get_stream_name(1), year2, month2, day2, year1, month1, day1);
        }

        for (int fieldID = 0; fieldID < numFields2; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);
          cdo_read_field(streamID2, varsData2[varID][levelID]);
        }

        tsID2++;
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
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
