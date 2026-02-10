/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Hi      hi           Compute the humidity index
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"

static const int FIRST_VAR = 0;

static double
humidityIndex(double t, double e, double r, double missval)
{
  constexpr double tmin = 26.0;
  constexpr double rmin = 40.0;

  if (t < tmin || r < rmin) return missval;

  return t + (5.0 / 9.0) * ((0.01 * r * e * 6.112 * std::pow(10.0, (7.5 * t) / (237.7 + t))) - 10.0);
}

static void
farexpr(Field &field1, Field const &field2, Field const &field3, double (*expression)(double, double, double, double))
{
  auto grid1 = field1.grid;
  auto grid2 = field2.grid;
  auto grid3 = field3.grid;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto missval3 = field3.missval;

  auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2) || len != gridInqSize(grid3)) cdo_abort("Fields have different gridsize (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals || field3.numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (fp_is_equal(field1.vec_d[i], missval1) || fp_is_equal(field2.vec_d[i], missval2)
          || fp_is_equal(field3.vec_d[i], missval3))
        field1.vec_d[i] = missval1;
      else
        field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], field3.vec_d[i], missval1);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) field1.vec_d[i] = expression(field1.vec_d[i], field2.vec_d[i], field3.vec_d[i], missval1);
  }

  field_num_mv(field1);
}

class Hi : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Hi",
    .operators = { { "hi" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Hi> registration = RegisterEntry<Hi>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };
  int vlistID3{};
  int vlistID4{};

  int taxisID1{ CDI_UNDEFID };
  // int taxisID2
  // int taxisID3
  int taxisID4{};
  int varID4{};

  VarList varList1{};
  VarList varList2{};
  VarList varList3{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    vlistID3 = cdo_stream_inq_vlist(streamID3);

    taxisID1 = vlistInqTaxis(vlistID1);
    // taxisID2 = vlistInqTaxis(vlistID2);
    // taxisID3 = vlistInqTaxis(vlistID3);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList3 = VarList(vlistID3);

    varList_compare(varList1, varList2);
    varList_compare(varList1, varList3);

    for (auto &var : varList1.vars) var.memType = MemType::Double;
    for (auto &var : varList2.vars) var.memType = MemType::Double;
    for (auto &var : varList3.vars) var.memType = MemType::Double;

    if (Options::cdoVerbose)
      cdo_print("Number of timesteps: file1 %d, file2 %d, file3 %d", varList1.numSteps(), varList2.numSteps(), varList3.numSteps());

    vlistID4 = vlistCreate();
    auto gridID = varList1.vars[FIRST_VAR].gridID;
    auto zaxisID = varList1.vars[FIRST_VAR].zaxisID;
    varID4 = vlistDefVar(vlistID4, gridID, zaxisID, TIME_VARYING);

    taxisID4 = cdo_taxis_create(TAXIS_RELATIVE);
    taxisDefTunit(taxisID4, TUNIT_MINUTE);
    taxisDefCalendar(taxisID4, CALENDAR_STANDARD);
    taxisDefRdatetime(taxisID4, cdiDateTime_set(19550101, 0));
    vlistDefTaxis(vlistID4, taxisID4);

    constexpr char HI_NAME[] = "hum_index";
    constexpr char HI_LONGNAME[]
        = "Humindex describes empirically in units of temperature how the temperature and humidity influence the wellness of a "
          "human "
          "being. HI = T + 5/9 * (A - 10) with A = e * (6.112 * 10 ** ((7.5 * T)/(237.7 + T)) * R), T  = air temperature in degree "
          "Celsius, R = relative humidity, e = vapour pressure. Humindex is only defined for temperatures of at least 26 degree "
          "Celsius and relative humidity of at least 40 percent.";
    constexpr char HI_UNITS[] = "Celsius";

    cdiDefKeyString(vlistID4, varID4, CDI_KEY_NAME, HI_NAME);
    cdiDefKeyString(vlistID4, varID4, CDI_KEY_LONGNAME, HI_LONGNAME);
    cdiDefKeyString(vlistID4, varID4, CDI_KEY_UNITS, HI_UNITS);

    streamID4 = cdo_open_write(3);

    cdo_def_vlist(streamID4, vlistID4);
  }

  void
  run() override
  {
    Field field1, field2, field3;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      auto numFields3 = cdo_stream_inq_timestep(streamID3, tsID);
      if (numFields2 == 0 || numFields3 == 0) cdo_abort("Input streams have different number of timesteps!");

      cdo_taxis_copy_timestep(taxisID4, taxisID1);
      cdo_def_timestep(streamID4, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID1, levelID1] = cdo_inq_field(streamID1);
        field1.init(varList1.vars[varID1]);
        cdo_read_field(streamID1, field1);

        auto [varID2, levelID2] = cdo_inq_field(streamID2);
        field2.init(varList2.vars[varID2]);
        cdo_read_field(streamID2, field2);

        auto [varID3, levelID3] = cdo_inq_field(streamID3);
        field3.init(varList3.vars[varID3]);
        cdo_read_field(streamID3, field3);

        if (varID1 != varID2 || varID1 != varID3 || levelID1 != levelID2 || levelID1 != levelID3)
          cdo_abort("Input streams have different structure!");

        if (varID1 != FIRST_VAR) continue;

        farexpr(field1, field2, field3, humidityIndex);

        cdo_def_field(streamID4, varID4, levelID1);
        cdo_write_field(streamID4, field1);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID4);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
