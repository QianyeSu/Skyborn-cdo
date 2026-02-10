/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Inttime    inttime         Time interpolation
*/

#include "cdi.h"
#include "julian_date.h"

#include <utility>

#include "cdo_options.h"
#include "cdo_omp.h"
#include "process_int.h"
#include "datetime.h"
#include "printinfo.h"
#include "field_functions.h"

template <typename T>
size_t
interp_time(double fac1, double fac2, size_t n, Varray<T> const &v1, Varray<T> const &v2, Varray<T> &v3, bool withMissval,
            T missval)
{
  size_t numMissVals3 = 0;

  if (withMissval)
  {
    for (size_t i = 0; i < n; ++i)
    {
      auto v1IsMissval = fp_is_equal(v1[i], missval);
      auto v2IsMissval = fp_is_equal(v2[i], missval);
      if (!v1IsMissval && !v2IsMissval) { v3[i] = v1[i] * fac1 + v2[i] * fac2; }
      else if (fac2 >= 0.5 && v1IsMissval && !v2IsMissval) { v3[i] = v2[i]; }
      else if (fac1 >= 0.5 && v2IsMissval && !v1IsMissval) { v3[i] = v1[i]; }
      else
      {
        v3[i] = missval;
        numMissVals3++;
      }
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
    for (size_t i = 0; i < n; ++i) { v3[i] = v1[i] * fac1 + v2[i] * fac2; }
  }

  return numMissVals3;
}

void
interp_time(double fac1, double fac2, Field const &field1, Field const &field2, Field &field3, bool withMissval)
{
  if (field1.memType != field3.memType) cdo_abort("Interal error, memType of field1 and field3 differ!");
  if (field3.memType == MemType::Float)
    field3.numMissVals
        = interp_time(fac1, fac2, field3.gridsize, field1.vec_f, field2.vec_f, field3.vec_f, withMissval, (float) field3.missval);
  else
    field3.numMissVals
        = interp_time(fac1, fac2, field3.gridsize, field1.vec_d, field2.vec_d, field3.vec_d, withMissval, field3.missval);
}

static void
julianDate_add_increment(JulianDate &julianDate, int64_t ijulinc, int calendar, int timeUnits)
{
  if (timeUnits == TUNIT_MONTH || timeUnits == TUNIT_YEAR)
  {
    auto vDateTime = julianDate_decode(calendar, julianDate);

    int year, month, day;
    cdiDate_decode(vDateTime.date, &year, &month, &day);

    month += (int) ijulinc;
    adjust_month_and_year(month, year);

    vDateTime.date = cdiDate_encode(year, month, day);
    julianDate = julianDate_encode(calendar, vDateTime);
  }
  else { julianDate = julianDate_add_seconds(julianDate, ijulinc); }
}

static void
adjust_time_units(int taxisID, int timeUnitsOut)
{
  auto timeUnitsIn = taxisInqTunit(taxisID);
  if ((timeUnitsIn == TUNIT_MONTH || timeUnitsIn == TUNIT_YEAR) && timeUnitsOut != TUNIT_MONTH && timeUnitsIn != TUNIT_YEAR)
  {
    taxisDefTunit(taxisID, timeUnitsOut);
  }
}

class Inttime : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Inttime",
    .operators = { { "inttime", InttimeHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Inttime> registration = RegisterEntry<Inttime>(module);

  CdiDateTime sDateTime{};
  int incrPeriod = 0, incrUnits = 3600, timeUnits = TUNIT_HOUR;

  CdoStreamID streamID1{};
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int curFirst = 0, curSecond = 1;

  int64_t ijulinc{};

  VarList varList1{};

public:
  void
  init() override
  {
    operator_input_arg("date,time<,increment> (format YYYY-MM-DD,hh:mm:ss)");
    if (cdo_operator_argc() < 2) cdo_abort("Too few arguments!");

    sDateTime.date = decode_datestring(cdo_operator_argv(0));
    sDateTime.time = decode_timestring(cdo_operator_argv(1));
    if (cdo_operator_argc() == 3) decode_timeunits(cdo_operator_argv(2), incrPeriod, incrUnits, timeUnits);

    // increment in seconds
    ijulinc = (int64_t) incrPeriod * incrUnits;

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    if (ijulinc == 0) vlistDefNtsteps(vlistID2, 1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    adjust_time_units(taxisID2, timeUnits);
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);
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
    auto julianDate = julianDate_encode(calendar, sDateTime);

    if (Options::cdoVerbose)
    {
      cdo_print("Start date/time %s", datetime_to_string(sDateTime));
      cdo_print("julianDate = %f", julianDate_to_seconds(julianDate));
      cdo_print("ijulinc = %lld", ijulinc);
    }

    int tsID = 0;
    int tsIDo = 0;

    auto numFields = cdo_stream_inq_timestep(streamID1, tsID++);
    auto vDateTime1 = taxisInqVdatetime(taxisID1);
    auto julianDate1 = julianDate_encode(calendar, vDateTime1);
    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto &field = varsData[curFirst][varID][levelID];
      cdo_read_field(streamID1, field);
    }

    if (Options::cdoVerbose)
    {
      cdo_print("Dataset begins on %s", datetime_to_string(vDateTime1));
      cdo_print("julianDate1 = %f", julianDate_to_seconds(julianDate1));
    }

    if (julianDate_to_seconds(julianDate1) > julianDate_to_seconds(julianDate))
    {
      cdo_print("Dataset begins on %s", datetime_to_string(vDateTime1));
      cdo_warning("The start time %s is before the beginning of the dataset!", datetime_to_string(sDateTime));
    }

    while (julianDate_to_seconds(julianDate1) <= julianDate_to_seconds(julianDate))
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID++);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);
      auto julianDate2 = julianDate_encode(calendar, vDateTime);
      if (Options::cdoVerbose)
      {
        cdo_print("date/time: %s", datetime_to_string(vDateTime));
        cdo_print("julianDate2 = %f", julianDate_to_seconds(julianDate2));
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        fieldInfoList[fieldID].set(varID, levelID);
        auto &field = varsData[curSecond][varID][levelID];
        cdo_read_field(streamID1, field);
      }

      while (julianDate_to_seconds(julianDate) <= julianDate_to_seconds(julianDate2))
      {
        if (julianDate_to_seconds(julianDate) >= julianDate_to_seconds(julianDate1)
            && julianDate_to_seconds(julianDate) <= julianDate_to_seconds(julianDate2))
        {
          auto dt = julianDate_decode(calendar, julianDate);

          if (Options::cdoVerbose)
            cdo_print("%s %s  %f  %d", date_to_string(dt.date), time_to_string(dt.time), julianDate_to_seconds(julianDate),
                      calendar);

          if (streamID2 == CDO_STREAM_UNDEF)
          {
            streamID2 = cdo_open_write(1);
            cdo_def_vlist(streamID2, vlistID2);
          }

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

        if (ijulinc == 0) break;

        julianDate_add_increment(julianDate, ijulinc, calendar, timeUnits);
      }

      julianDate1 = julianDate2;
      std::swap(curFirst, curSecond);
    }

    if (tsIDo == 0) cdo_warning("Start date/time %s out of range, no time steps interpolated!", datetime_to_string(sDateTime));
  }

  void
  close() override
  {
    if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
