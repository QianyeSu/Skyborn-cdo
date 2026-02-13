/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Ymoncomp   ymoneq         Equal
      Ymoncomp   ymonne         Not equal
      Ymoncomp   ymonle         Less equal
      Ymoncomp   ymonlt         Less than
      Ymoncomp   ymonge         Greater equal
      Ymoncomp   ymongt         Greater than
      Yseascomp  yseaseq        Equal
      Yseascomp  yseasne        Not equal
      Yseascomp  yseasle        Less equal
      Yseascomp  yseaslt        Less than
      Yseascomp  yseasge        Greater equal
      Yseascomp  yseasgt        Greater than
*/

#include <cdi.h>

#include "cdo_season.h"
#include "cdo_varlist.h"
#include "process_int.h"
#include "field_functions.h"

template <typename T1, typename T2>
static void
comp_function(Varray<T1> &v1, Varray<T2> const &v2, double missval1, double missval2, size_t n, int operFunc, bool hasMissvals)
{
  T1 mv1 = missval1;
  T2 mv2 = missval2;

  auto func_comp = [&](auto binary_operator)
  {
    if (hasMissvals)
    {
      if (std::isnan(mv1) || std::isnan(mv2))
        for (size_t i = 0; i < n; ++i)
          v1[i] = (fp_is_equal(v1[i], mv1) || fp_is_equal(v2[i], mv2)) ? mv1 : binary_operator(v1[i], v2[i]);
      else
        for (size_t i = 0; i < n; ++i) v1[i] = (is_equal(v1[i], mv1) || is_equal(v2[i], mv2)) ? mv1 : binary_operator(v1[i], v2[i]);
    }
    else
    {
      for (size_t i = 0; i < n; ++i) v1[i] = binary_operator(v1[i], v2[i]);
    }
  };

  // clang-format off
  if      (operFunc == FieldFunc_EQ) func_comp(binary_op_EQ);
  else if (operFunc == FieldFunc_NE) func_comp(binary_op_NE);
  else if (operFunc == FieldFunc_LE) func_comp(binary_op_LE);
  else if (operFunc == FieldFunc_LT) func_comp(binary_op_LT);
  else if (operFunc == FieldFunc_GE) func_comp(binary_op_GE);
  else if (operFunc == FieldFunc_GT) func_comp(binary_op_GT);
  else cdo_abort("Operator not implemented!");
  // clang-format on
}

static void
comp_function(Field &f1, Field const &f2, int operFunc)
{
  auto hasMissvals = (f1.numMissVals > 0 || f2.numMissVals > 0);
  auto func = [&](auto &v1, auto const &v2) { comp_function(v1, v2, f1.missval, f2.missval, f1.size, operFunc, hasMissvals); };
  field_operation2(func, f1, f2);
}

constexpr int MaxMonths = 12;

static const char *monthNames[]
    = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

static int
get_month_index(CdiDate vDate)
{
  int month = vDate.month;
  if (month < 1 || month > MaxMonths) cdo_abort("Month %d out of range!", month);
  month--;
  return month;
}

enum
{
  MONTHLY,
  SEASONAL
};

class Ymoncomp : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Ymoncomp",
    .operators = { { "ymoneq", FieldFunc_EQ, MONTHLY, YmoncompHelp },
                   { "ymonne", FieldFunc_NE, MONTHLY, YmoncompHelp },
                   { "ymonle", FieldFunc_LE, MONTHLY, YmoncompHelp },
                   { "ymonlt", FieldFunc_LT, MONTHLY, YmoncompHelp },
                   { "ymonge", FieldFunc_GE, MONTHLY, YmoncompHelp },
                   { "ymongt", FieldFunc_GT, MONTHLY, YmoncompHelp },
                   { "yseaseq", FieldFunc_EQ, SEASONAL, YseascompHelp },
                   { "yseasne", FieldFunc_NE, SEASONAL, YseascompHelp },
                   { "yseasle", FieldFunc_LE, SEASONAL, YseascompHelp },
                   { "yseaslt", FieldFunc_LT, SEASONAL, YseascompHelp },
                   { "yseasge", FieldFunc_GE, SEASONAL, YseascompHelp },
                   { "yseasgt", FieldFunc_GT, SEASONAL, YseascompHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Ymoncomp> registration = RegisterEntry<Ymoncomp>();

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

      comp_function(field, varsData2[mon][varID][levelID], operfunc);

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

      load_data_from_stream(mon, numFields);
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
    operfunc = cdo_operator_f1(operatorID);
    opertype = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    vlistDestroy(vlistID3);
  }

  void
  run() override
  {
    (opertype == SEASONAL) ? run_seasonal() : run_monthly();
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID1);
  }
};
