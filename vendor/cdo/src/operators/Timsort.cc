/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Timsort    timsort         Sort over the time
*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_options.h"
#include "cdo_omp.h"
#include "field_functions.h"

class Timsort : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timsort",
    .operators = { { "timsort", TimsortHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Timsort> registration = RegisterEntry<Timsort>();

  int nalloc = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1;

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector3D varsData;
    std::vector<CdiDateTime> vDateTimes;

    auto numVars = varList1.numVars();

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (tsID >= nalloc)
      {
        constexpr int NALLOC_INC = 1024;
        nalloc += NALLOC_INC;
        vDateTimes.resize(nalloc);
        varsData.resize(nalloc);
      }

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

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
      auto const &var = varList1.vars[varID];

      if (var.isConstant) continue;

      auto memType = var.memType;
      auto gridsize = var.gridsize;
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
        for (size_t i = 0; i < gridsize; ++i)
        {
          auto ompthID = cdo_omp_get_thread_num();

          if (memType == MemType::Float)
          {
            auto &v = fields[ompthID].vec_f;
            v.resize(nts);
            for (int t = 0; t < nts; ++t) v[t] = varsData[t][varID][levelID].vec_f[i];

            std::ranges::sort(v);

            for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_f[i] = v[t];
          }
          else
          {
            auto &v = fields[ompthID].vec_d;
            v.resize(nts);
            for (int t = 0; t < nts; ++t) v[t] = varsData[t][varID][levelID].vec_d[i];

            std::ranges::sort(v);

            for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_d[i] = v[t];
          }
        }
      }
    }

    for (tsID = 0; tsID < nts; ++tsID)
    {
      taxisDefVdatetime(taxisID2, vDateTimes[tsID]);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto &field = varsData[tsID][varID][levelID];
          if (field.hasData())
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field);
          }
        }
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
