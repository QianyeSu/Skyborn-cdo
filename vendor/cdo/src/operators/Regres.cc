/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Regres      regres           Regression
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_options.h"
#include "field_trend.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "field_functions.h"

// Same code as Trend!

static bool
get_parameter()
{
  bool tstepIsEqual = true;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "equal") tstepIsEqual = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return tstepIsEqual;
}

class Regres : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Regres",
    .operators = { { "regres", RegresHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Regres> registration = RegisterEntry<Regres>();

  CdoStreamID streamID1{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  bool tstepIsEqual{ true };

  VarList varList1{};
  VarList varList2{};

public:
  void
  init() override
  {
    tstepIsEqual = get_parameter();

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlistDefNtsteps(vlistID2, 1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID3 = cdo_open_write(1);
    cdo_def_vlist(streamID3, vlistID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
  }

  void
  run() override
  {
    for (auto &var : varList2.vars) var.memType = MemType::Double;

    Field field1, field2, field3;

    constexpr size_t numWork = 5;
    FieldVector3D work(numWork);
    for (auto &w : work) field2D_init(w, varList1, FIELD_VEC, 0);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    auto calendar = taxisInqCalendar(taxisID1);

    CheckTimeIncr checkTimeIncr;
    JulianDate julianDate0;
    double deltat1 = 0;
    CdiDateTime vDateTime{};
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        fieldInfoList[fieldID].set(varID, levelID);

        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

        calc_trend_sum(work, field1, zj, varID, levelID);
      }

      tsID++;
    }

    taxisDefVdatetime(taxisID2, vDateTime);
    cdo_def_timestep(streamID3, 0);

    for (int fieldID = 0; fieldID < maxFields; ++fieldID)
    {
      auto [varID, levelID] = fieldInfoList[fieldID].get();

      auto const &var = varList2.vars[varID];
      field2.init(var);
      field3.init(var);

      calc_trend_param(work, field2, field3, varID, levelID);

      field_num_mv(field3);

      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, field3);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);
  }
};
