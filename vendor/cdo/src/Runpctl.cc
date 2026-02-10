/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Quast

*/

/*
   This module contains the following operators:

      Runpctl    runpctl         Running percentiles
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "percentiles.h"
#include "datetime.h"
#include "field_functions.h"
#include "cdo_omp.h"

template <typename T>
static size_t
runpctl(double pn, int ndates, size_t gridsize, Varray<T> &v2, double mv, const FieldVector3D &vars1, int varID, int levelID,
        MemType memType)
{
  T missval = mv;
  size_t numMissVals = 0;
  Varray2D<T> array_2D(Threading::ompNumMaxThreads, Varray<T>(ndates));

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto ompthID = cdo_omp_get_thread_num();
    auto &array = array_2D[ompthID];

    int j = 0;

    if (memType == MemType::Float)
    {
      for (int inp = 0; inp < ndates; ++inp)
      {
        auto val = vars1[inp][varID][levelID].vec_f[i];
        if (fp_is_not_equal(val, static_cast<float>(missval))) array[j++] = val;
      }
    }
    else
    {
      for (int inp = 0; inp < ndates; ++inp)
      {
        auto val = vars1[inp][varID][levelID].vec_d[i];
        if (fp_is_not_equal(val, missval)) array[j++] = val;
      }
    }
    /*
    for (int inp = 0; inp < ndates; ++inp)
      {
        auto func = [&](auto &v) {
          auto val = v[i];
          if (fp_is_not_equal(val, missval)) array[j++] = val;
        };
        field_operation(func, vars1[inp][varID][levelID]);
      }
    */
    if (j > 0) { v2[i] = percentile(array.data(), j, pn); }
    else
    {
      v2[i] = missval;
      numMissVals++;
    }
  }

  return numMissVals;
}

static void
runpctl(double pn, int ndates, Field &field1, const FieldVector3D &vars1, int varID, int levelID)
{
  auto func = [&](auto &v)
  { field1.numMissVals = runpctl(pn, ndates, field1.gridsize, v, field1.missval, vars1, varID, levelID, field1.memType); };
  field_operation(func, field1);
}

class Runpctl : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Runpctl",
    .operators = { { "runpctl", RunpctlHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Runpctl> registration = RegisterEntry<Runpctl>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };
  VarList varList1{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  FieldVector3D varsData1;

  DateTimeList dtlist{};
  double pn{};
  int ndates{};
  int maxFields{};
  int tsID{};
  std::vector<FieldInfo> fieldInfoList;

public:
  void
  init() override
  {
    constexpr auto timestatDate{ TimeStat::MEAN };

    operator_input_arg("percentile number, number of timesteps");
    operator_check_argc(2);
    pn = parameter_to_double(cdo_operator_argv(0));
    ndates = parameter_to_int(cdo_operator_argv(1));

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisWithBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);

    maxFields = varList1.maxFields();
    fieldInfoList = std::vector<FieldInfo>(maxFields);

    dtlist.set_stat(timestatDate);
    dtlist.set_calendar(taxisInqCalendar(taxisID1));

    varsData1 = FieldVector3D(ndates + 1);
    for (int its = 0; its < ndates; its++) field2D_init(varsData1[its], varList1);
  }

  void
  write_fields(int otsID)
  {
    dtlist.stat_taxis_def_timestep(taxisID2, ndates);
    cdo_def_timestep(streamID2, otsID);
    for (int fieldID = 0; fieldID < maxFields; ++fieldID)
    {
      auto [varID, levelID] = fieldInfoList[fieldID].get();
      if (otsID && varList1.vars[varID].isConstant) continue;

      cdo_def_field(streamID2, varID, levelID);
      auto &field1 = varsData1[0][varID][levelID];
      cdo_write_field(streamID2, field1);
    }
  }

  void
  run() override
  {
    for (tsID = 0; tsID < ndates; ++tsID)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) cdo_abort("File has less than %d timesteps!", ndates);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

        auto &field = varsData1[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }
    }
    int otsID = 0;
    while (true)
    {
      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        if (varList1.vars[varID].isConstant) continue;

        auto nlevels = varList1.vars[varID].nlevels;
        for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          auto &field1 = varsData1[0][varID][levelID];
          runpctl(pn, ndates, field1, varsData1, varID, levelID);
        }
      }

      write_fields(otsID);
      otsID++;

      dtlist.shift();

      varsData1[ndates] = varsData1[0];
      for (int inp = 0; inp < ndates; ++inp) varsData1[inp] = varsData1[inp + 1];

      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      dtlist.taxis_inq_timestep(taxisID1, ndates - 1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &fieldN = varsData1[ndates - 1][varID][levelID];
        cdo_read_field(streamID1, fieldN);
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
