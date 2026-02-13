/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithdays  muldpm          Multiply with days per month
      Arithdays  divdpm          Divide by days per month
      Arithdays  muldpy          Multiply with days per year
      Arithdays  divdpy          Divide by days per year
      Arithdays  muldoy          Multiply with day of year
*/

#include <cdi.h>
#include "calendar.h"

#include "cdo_options.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

static double
dayofyear(int calendar, CdiDateTime const &vDateTime)
{
  constexpr int month_360[12] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };
  constexpr int month_365[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  constexpr int month_366[12] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  int year, month, day;
  int hour, minute, second, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  cdiTime_decode(vDateTime.time, &hour, &minute, &second, &ms);

  auto dpy = days_per_year(calendar, year);

  double doy = 0.0;
  for (int im = 1; im < month; ++im)
  {
    auto dpm = (dpy == 360) ? month_360 : ((dpy == 365) ? month_365 : month_366);
    if (im <= 12) doy += dpm[im - 1];
  }

  doy += (day - 1);
  doy += (second + minute * 60 + hour * 3600) / 86400.0;

  if (Options::cdoVerbose) cdo_print("vDateTime, dpy, doy: %s %d %g", datetime_to_string(vDateTime), dpy, doy);

  return doy;
}

class Arithdays : public Process
{
  enum
  {
    Func_Month = 1,
    Func_Year,
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Arithdays",
    .operators = { { "muldpm", FieldFunc_Mul, Func_Month, ArithdaysHelp },
                   { "divdpm", FieldFunc_Div, Func_Month, ArithdaysHelp },
                   { "muldpy", FieldFunc_Mul, Func_Year, ArithdaysHelp },
                   { "divdpy", FieldFunc_Div, Func_Year, ArithdaysHelp },
                   { "muldoy", FieldFunc_Mul, 0, ArithdaysHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Arithdays> registration = RegisterEntry<Arithdays>();

private:
  int MULDOY{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operatorID{};
  int operfunc{};
  int operfunc2{};
  int calendar{};

  VarList varList1;

  double
  calc_days(CdiDateTime const &vDateTime)
  {
    int year, month, day;
    cdiDate_decode(vDateTime.date, &year, &month, &day);

    double rconst;
    if (operatorID == MULDOY)
      rconst = dayofyear(calendar, vDateTime);
    else
      rconst = (operfunc2 == Func_Month) ? days_per_month(calendar, year, month) : days_per_year(calendar, year);

    if (Options::cdoVerbose) cdo_print("calendar %d  year %d  month %d  result %g", calendar, year, month, rconst);

    return rconst;
  }

public:
  void
  init() override
  {
    MULDOY = module.get_id("muldoy");

    operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    operfunc2 = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    calendar = taxisInqCalendar(taxisID1);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      auto rconst = calc_days(vDateTime);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        fieldc_function(field, rconst, operfunc);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
