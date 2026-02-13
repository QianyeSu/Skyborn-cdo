/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/

/*
   This module contains the following operators:

      Wct     wct          Compute the windchill temperature (degree C)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"

static const int FIRST_VAR = 0;

static double
windchillTemperature(double t, double ff, double missval)
{
  constexpr double tmax = 33.0;
  constexpr double vmin = 1.39;  // minimum wind speed (m/s)

  return ff < vmin || t > tmax ? missval : tmax + (t - tmax) * (0.478 + 0.237 * (std::sqrt(ff * 3.6) - 0.0124 * ff * 3.6));
}

static void
farexpr(Field &field1, Field &field2, double (*expression)(double, double, double))
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (fp_is_equal(field1.vec_d[i], missval1) || fp_is_equal(field2.vec_d[i], missval2))
        field1.vec_d[i] = missval1;
      else
        field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], missval1);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], missval1);
  }

  field_num_mv(field1);
}

class Wct : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Wct",
    .operators = { { "wct", WctHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Wct> registration = RegisterEntry<Wct>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1{ CDI_UNDEFID };
  int taxisID3;

  VarList varList1;
  VarList varList2;

  int varID3;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2, CmpVarList::Dim);
    for (auto &var : varList1.vars) var.memType = MemType::Double;
    for (auto &var : varList2.vars) var.memType = MemType::Double;

    taxisID1 = vlistInqTaxis(vlistID1);

    if (Options::cdoVerbose) cdo_print("Number of timesteps: file1 %d, file2 %d", varList1.numSteps(), varList2.numSteps());

    auto vlistID3 = vlistCreate();
    auto gridID = varList1.vars[FIRST_VAR].gridID;
    auto zaxisID = varList1.vars[FIRST_VAR].zaxisID;
    varID3 = vlistDefVar(vlistID3, gridID, zaxisID, TIME_VARYING);

    taxisID3 = cdo_taxis_create(TAXIS_RELATIVE);
    taxisDefTunit(taxisID3, TUNIT_MINUTE);
    taxisDefCalendar(taxisID3, CALENDAR_STANDARD);
    taxisDefRdatetime(taxisID3, cdiDateTime_set(19550101, 0));
    vlistDefTaxis(vlistID3, taxisID3);

    constexpr char WCT_NAME[] = "wind_chill_temperature";
    constexpr char WCT_LONGNAME[] = "Windchill temperature describes the fact that low temperatures are felt "
                                    "to be even lower in case of wind. It is based on the rate of heat loss "
                                    "from exposed skin caused by wind and cold. It is calculated according "
                                    "to the empirical formula: 33 + (T - 33) * (0.478 + 0.237 * ( "
                                    "SQRT(ff*3.6) - 0.0124 * ff * 3.6)) with T  = air temperature in "
                                    "degree Celsius, ff = 10 m wind speed in m/s. Windchill temperature is "
                                    "only defined for temperatures at or below 33 degree Celsius and wind "
                                    "speeds above 1.39 m/s. It is mainly used for freezing temperatures.";
    constexpr char WCT_UNITS[] = "Celsius";
    cdiDefKeyString(vlistID3, varID3, CDI_KEY_NAME, WCT_NAME);
    cdiDefKeyString(vlistID3, varID3, CDI_KEY_LONGNAME, WCT_LONGNAME);
    cdiDefKeyString(vlistID3, varID3, CDI_KEY_UNITS, WCT_UNITS);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Field field1, field2;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0) cdo_abort("Input streams have different number of timesteps!");

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID1, levelID1] = cdo_inq_field(streamID1);
        field1.init(varList1.vars[varID1]);
        cdo_read_field(streamID1, field1);

        auto [varID2, levelID2] = cdo_inq_field(streamID2);
        field2.init(varList2.vars[varID2]);
        cdo_read_field(streamID2, field2);

        if (varID1 != varID2 || levelID1 != levelID2) cdo_abort("Input streams have different structure!");

        if (varID1 != FIRST_VAR) continue;

        farexpr(field1, field2, windchillTemperature);

        cdo_def_field(streamID3, varID3, levelID1);
        cdo_write_field(streamID3, field1);
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
