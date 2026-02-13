/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Tstepcount  tstepcount  Count number of timesteps
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "cdo_omp.h"
#include "field_functions.h"

template <typename T>
static T
tstepcount(long nts, T missval, Varray<T> const &v, T refval)
{
  if (fp_is_equal(refval, missval)) return missval;

  long j;
  long n = 0;
  for (j = 0; j < nts; ++j)
  {
    n++;
    if (fp_is_equal(v[j], refval)) break;
  }

  return (j == nts) ? missval : (T) n;
}

class Tstepcount : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Tstepcount",
    .operators = { { "tstepcount" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Tstepcount> registration = RegisterEntry<Tstepcount>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numVars{};
  double refval{};

  VarList varList1{};

public:
  void
  init() override
  {
    refval = (cdo_operator_argc() == 1) ? parameter_to_double(cdo_operator_argv(0)) : 0.0;

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlistDefNtsteps(vlistID2, 1);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "steps");

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    FieldVector3D varsData;
    CdiDateTime vDateTime{};

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);

      vDateTime = taxisInqVdatetime(taxisID1);

      field2D_init(varsData[tsID], varList1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field = varsData[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }

      tsID++;
    }

    int nts = tsID;

    std::vector<Field> fields(Threading::ompNumMaxThreads);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      auto memType = var1.memType;
      auto missval = var1.missval;
      auto gridsize = var1.gridsize;
      for (int levelID = 0; levelID < var1.nlevels; ++levelID)
      {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
        for (size_t i = 0; i < gridsize; ++i)
        {
          auto ompthID = cdo_omp_get_thread_num();

          if (memType == MemType::Float)
          {
            auto &v = fields[ompthID].vec_f;
            v.resize(nts);
            for (int t = 0; t < nts; ++t) v[t] = varsData[t][varID][levelID].vec_f[i];

            auto count = tstepcount(nts, (float) missval, v, (float) refval);

            varsData[0][varID][levelID].vec_f[i] = count;
          }
          else
          {
            auto &v = fields[ompthID].vec_d;
            v.resize(nts);
            for (int t = 0; t < nts; ++t) v[t] = varsData[t][varID][levelID].vec_d[i];

            auto count = tstepcount(nts, missval, v, refval);

            varsData[0][varID][levelID].vec_d[i] = count;
          }
        }
      }
    }

    taxisDefVdatetime(taxisID2, vDateTime);
    cdo_def_timestep(streamID2, 0);

    for (int varID = 0; varID < numVars; ++varID)
    {
      for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
      {
        cdo_def_field(streamID2, varID, levelID);
        auto &field1 = varsData[0][varID][levelID];
        field_num_mv(field1);
        cdo_write_field(streamID2, field1);
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
