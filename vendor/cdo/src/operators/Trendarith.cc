/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

   Trendarith    addtrend        Add trend
   Trendarith    subtrend        Subtract trend
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "cdo_omp.h"
#include "datetime.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "field_functions.h"
#include "arithmetic.h"

template <typename T>
static void
add_trend(double zj, Varray<T> &v1, Varray<double> const &v2, Varray<double> const &v3, size_t n, double mv)
{
  auto missval1 = mv;
  auto missval2 = mv;

  auto add_kernel = [&](auto is_EQ)
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n; ++i)
    {
      auto tmp = (is_EQ(v2[i], mv) || is_EQ(v3[i], mv)) ? mv : (v2[i] + v3[i] * zj);
      v1[i] = ADDM(v1[i], tmp);
    }
  };

  std::isnan(mv) ? add_kernel(fp_is_equal) : add_kernel(is_equal);
}

static void
add_trend(double zj, Field &field1, Field const &field2, Field const &field3)
{
  auto func = [&](auto &v1) { add_trend(zj, v1, field2.vec_d, field3.vec_d, field1.size, field1.missval); };
  field_operation(func, field1);
}

template <typename T>
static void
sub_trend(double zj, Varray<T> &v1, Varray<double> const &v2, Varray<double> const &v3, size_t n, double mv)
{
  auto missval1 = mv;
  auto missval2 = mv;

  auto sub_kernel = [&](auto is_EQ)
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n; ++i)
    {
      auto tmp = (is_EQ(v2[i], mv) || is_EQ(v3[i], mv)) ? mv : (v2[i] + v3[i] * zj);
      v1[i] = SUBM(v1[i], tmp);
    }
  };

  std::isnan(mv) ? sub_kernel(fp_is_equal) : sub_kernel(is_equal);
}

static void
sub_trend(double zj, Field &field1, Field const &field2, Field const &field3)
{
  auto func = [&](auto &v1) { sub_trend(zj, v1, field2.vec_d, field3.vec_d, field1.size, field1.missval); };
  field_operation(func, field1);
}

static void
get_parameter(bool &tstepIsEqual)
{
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
      if (key == "equal") tstepIsEqual = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Trendarith : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Trendarith",
    // clang-format off
    .operators = { { "addtrend", FieldFunc_Add, 0, TrendarithHelp },
                   { "subtrend", FieldFunc_Sub, 0, TrendarithHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, NoRestriction },
  };
  inline static RegisterEntry<Trendarith> registration = RegisterEntry<Trendarith>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  CdoStreamID streamID4{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID4{};

  int operfunc{};

  bool tstepIsEqual{};

  VarList varList1{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    tstepIsEqual = true;
    get_parameter(tstepIsEqual);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);
    streamID3 = cdo_open_read(2);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = cdo_stream_inq_vlist(streamID3);
    auto vlistID4 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID4);

    varList1 = VarList(vlistID1);
    VarList varList2(vlistID2);
    VarList varList3(vlistID3);

    varList_compare(varList1, varList2);
    varList_compare(varList1, varList3);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID4 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID4, taxisID4);

    streamID4 = cdo_open_write(3);
    cdo_def_vlist(streamID4, vlistID4);
  }

  void
  run() override
  {
    FieldVector2D vars2, vars3;
    field2D_init(vars2, varList1, FIELD_VEC);
    field2D_init(vars3, varList1, FIELD_VEC);

    {
      auto numFields = cdo_stream_inq_timestep(streamID2, 0);
      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        cdo_read_field(streamID2, vars2[varID][levelID]);
      }
    }

    {
      auto numFields = cdo_stream_inq_timestep(streamID3, 0);
      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID3);
        cdo_read_field(streamID3, vars3[varID][levelID]);
      }
    }

    auto calendar = taxisInqCalendar(taxisID1);
    CheckTimeIncr checkTimeIncr;
    JulianDate julianDate0;
    double deltat1 = 0.0;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      if (tstepIsEqual) check_time_increment(tsID, calendar, vDateTime, checkTimeIncr);
      auto zj = tstepIsEqual ? (double) tsID : delta_time_step_0(tsID, calendar, vDateTime, julianDate0, deltat1);

      cdo_taxis_copy_timestep(taxisID4, taxisID1);
      cdo_def_timestep(streamID4, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
        if (numSteps > 0) progress.update(fstatus);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        field.init(var1);
        cdo_read_field(streamID1, field);

        if (operfunc == FieldFunc_Add)
          add_trend(zj, field, vars2[varID][levelID], vars3[varID][levelID]);
        else
          sub_trend(zj, field, vars2[varID][levelID], vars3[varID][levelID]);

        cdo_def_field(streamID4, varID, levelID);
        cdo_write_field(streamID4, field);
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
