/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "field_functions.h"

template <typename T>
static void
field_scale(Varray<T> &v, size_t n, size_t numMissVals, double missval, Varray<double> const &scale)
{
  if (numMissVals)
  {
    for (size_t i = 0; i < n; ++i)
      if (fp_is_not_equal(v[i], missval)) v[i] *= scale[i];
  }
  else
  {
    for (size_t i = 0; i < n; ++i) v[i] *= scale[i];
  }
}

static void
field_scale(Field &field, Varray<double> const &scale)
{
  auto func = [&](auto &v, auto n, auto numMissVals, auto missval) { field_scale(v, n, numMissVals, missval, scale); };
  field_operation(func, field, field.gridsize, field.numMissVals, field.missval);
}

class Arithlat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Arithlat",
    // clang-format off
    .operators = { { "mulcoslat", FieldFunc_Mul, 0, ArithlatHelp },
                   { "divcoslat", FieldFunc_Div, 0, ArithlatHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Arithlat> registration = RegisterEntry<Arithlat>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operfunc{};

  VarList varList1;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Varray<double> scale;
    Field field;

    int gridID0 = -1;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        field.init(var1);
        cdo_read_field(streamID1, field);

        auto gridID = var1.gridID;
        auto gridsize = var1.gridsize;

        if (gridID != gridID0)
        {
          gridID0 = gridID;
          gridID = generate_full_point_grid(gridID);
          if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

          scale.resize(gridsize);
          gridInqYvals(gridID, scale.data());

          // Convert lat/lon units if required
          cdo_grid_to_radian(gridID, CDI_YAXIS, scale, "grid latitudes");

          if (operfunc == FieldFunc_Mul)
            for (size_t i = 0; i < gridsize; ++i) scale[i] = std::cos(scale[i]);
          else
            for (size_t i = 0; i < gridsize; ++i) scale[i] = 1.0 / std::cos(scale[i]);

          if (Options::cdoVerbose)
            for (int i = 0; i < 10; ++i) cdo_print("coslat  %3d  %g", i + 1, scale[i]);
        }

        field_scale(field, scale);

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
