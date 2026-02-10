/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Deltat    deltat         Delta t
*/

#include "cdi.h"
#include "julian_date.h"
#include "process_int.h"
#include "field_functions.h"

template <typename T>
static void
varray_deltat(size_t len, Varray<T> const &v0, Varray<T> const &v1, Varray<T> &v2, double factor, double mv)
{
  assert(len > 0);
  assert(v0.size() > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v0.size());
  assert(len <= v1.size());
  assert(len <= v2.size());

  for (size_t i = 0; i < len; ++i) { v2[i] = (fp_is_equal(v0[i], mv) || fp_is_equal(v1[i], mv)) ? mv : (v1[i] - v0[i]) * factor; }
}

template <typename T>
static void
varray_deltat(size_t len, Varray<T> const &v0, Varray<T> const &v1, Varray<T> &v2, double factor)
{
  assert(len > 0);

  for (size_t i = 0; i < len; ++i) v2[i] = (v1[i] - v0[i]) * factor;
}

void
field_deltat(Field const &field0, Field const &field1, Field &field2, double factor)
{
  if (field1.memType != field2.memType) cdo_abort("Interal error, memType of field1 and field2 differ!");

  if (field1.numMissVals || field0.numMissVals)
  {
    if (field1.memType == MemType::Float)
      varray_deltat(field1.size, field0.vec_f, field1.vec_f, field2.vec_f, factor, field1.missval);
    else
      varray_deltat(field1.size, field0.vec_d, field1.vec_d, field2.vec_d, factor, field1.missval);

    field2.numMissVals = field_num_mv(field2);
  }
  else
  {
    if (field1.memType == MemType::Float)
      varray_deltat(field1.size, field0.vec_f, field1.vec_f, field2.vec_f, factor);
    else
      varray_deltat(field1.size, field0.vec_d, field1.vec_d, field2.vec_d, factor);
  }
}

class Deltat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Deltat",
    .operators = { { "deltat", DeltatHelp }, { "timederivative", 0, 1, DeltatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Deltat> registration = RegisterEntry<Deltat>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int calendar{ CDI_UNDEFID };
  VarList varList1{};
  bool ldivdt{ false };

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    ldivdt = cdo_operator_f2(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    calendar = taxisInqCalendar(taxisID1);

    varList1 = VarList(vlistID1);
    auto numSteps = varList1.numSteps();
    if (numSteps > 1) vlistDefNtsteps(vlistID2, numSteps - 1);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field1, field2;
    FieldVector2D varsData;
    field2D_init(varsData, varList1, FIELD_VEC | FIELD_NAT);

    int tsID = 0;
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID);

    auto julianDate0 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      cdo_read_field(streamID1, varsData[varID][levelID]);
    }

    tsID++;
    int tsID2 = 0;
    while (true)
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));
      auto idtInSec = ldivdt ? 1.0 / julianDate_to_seconds(julianDate_sub(julianDate1, julianDate0)) : 1.0;
      julianDate0 = julianDate1;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID2);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        field1.init(var1);
        cdo_read_field(streamID1, field1);

        auto &field0 = varsData[varID][levelID];

        field2.init(var1);
        field_deltat(field0, field1, field2, idtInSec);

        field_copy(field1, field0);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field2);
      }

      tsID++;
      tsID2++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
