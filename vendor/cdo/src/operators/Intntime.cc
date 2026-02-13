/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intntime   intntime        Time interpolation
*/

#include "cdi.h"
#include "julian_date.h"

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "printinfo.h"
#include "field_functions.h"

void interp_time(double fac1, double fac2, Field const &field1, Field const &field2, Field &field3, bool withMissval);

class Intntime : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intntime",
    .operators = { { "intntime", InttimeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Intntime> registration = RegisterEntry<Intntime>();

  int curFirst = 0, curSecond = 1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numts{};
  VarList varList1{};

public:
  void
  init() override
  {
    operator_input_arg("number of timesteps between 2 timesteps");
    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");

    numts = parameter_to_int(cdo_operator_argv(0));
    if (numts < 2) cdo_abort("parameter must be greater than 1!");

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);

    vlistDefNtsteps(vlistID2, -1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field3;
    FieldVector2D varsData[2];
    field2D_init(varsData[0], varList1, FIELD_VEC | FIELD_NAT);
    field2D_init(varsData[1], varList1, FIELD_VEC | FIELD_NAT);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    auto calendar = taxisInqCalendar(taxisID1);

    int tsID = 0;
    int tsIDo = 0;

    auto numFields = cdo_stream_inq_timestep(streamID1, tsID++);
    auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

    cdo_taxis_copy_timestep(taxisID2, taxisID1);
    cdo_def_timestep(streamID2, tsIDo++);
    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto &field = varsData[curFirst][varID][levelID];
      cdo_read_field(streamID1, field);

      cdo_def_field(streamID2, varID, levelID);
      cdo_write_field(streamID2, field);
    }

    while (true)
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID++);
      if (numFields == 0) break;

      auto vDateTime2 = taxisInqVdatetime(taxisID1);
      auto julianDate2 = julianDate_encode(calendar, vDateTime2);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        fieldInfoList[fieldID].set(varID, levelID);
        auto &field = varsData[curSecond][varID][levelID];
        cdo_read_field(streamID1, field);
      }

      for (int it = 1; it < numts; it++)
      {
        auto seconds = it * julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / numts;
        auto julianDate = julianDate_add_seconds(julianDate1, std::lround(seconds));
        auto dt = julianDate_decode(calendar, julianDate);

        if (Options::cdoVerbose) cdo_print("%s %s", date_to_string(dt.date), time_to_string(dt.time));

        taxisDefVdatetime(taxisID2, dt);
        cdo_def_timestep(streamID2, tsIDo++);

        auto diff = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1));
        auto fac1 = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate)) / diff;
        auto fac2 = julianDate_to_seconds(julianDate_sub(julianDate, julianDate1)) / diff;

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();

          auto const &field1 = varsData[curFirst][varID][levelID];
          auto const &field2 = varsData[curSecond][varID][levelID];

          field3.init(varList1.vars[varID]);

          auto withMissval = (field1.numMissVals || field2.numMissVals);
          interp_time(fac1, fac2, field1, field2, field3, withMissval);

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, field3);
        }
      }

      taxisDefVdatetime(taxisID2, vDateTime2);
      cdo_def_timestep(streamID2, tsIDo++);
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        auto &field = varsData[curSecond][varID][levelID];
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      julianDate1 = julianDate2;
      std::swap(curFirst, curSecond);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
