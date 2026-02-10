/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

     Sorttimestamp    sorttimestamp         Sort all timesteps
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "printinfo.h"
#include "field_functions.h"

bool getenv_skip_same_time();

struct TimeInfo
{
  int index;
  double datetime;
};

class Sorttimestamp : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Sorttimestamp",
    .operators = { { "sorttimestamp" }, { "sorttaxis" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { -1, 1, NoRestriction },
  };
  inline static RegisterEntry<Sorttimestamp> registration = RegisterEntry<Sorttimestamp>(module);

  int lasttsID = -1;
  int nalloc = 0;
  int vlistID2 = -1, taxisID2 = -1;
  int numVars = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  bool unique = false;
  bool skipSameTime{};

public:
  void
  init() override
  {
    skipSameTime = getenv_skip_same_time();

    if (cdo_operator_argc() == 1)
    {
      if (cdo_operator_argv(0) == "unique")
        unique = true;
      else
        cdo_abort("Unexpected parameter %s!", cdo_operator_argv(0));
    }
    if (cdo_operator_argc() > 1) operator_check_argc(1);
  }

  void
  run() override
  {
    FieldVector3D varsData;
    std::vector<CdiDateTime> vDateTimes;
    auto numFiles = cdo_stream_cnt() - 1;

    int xtsID = 0;
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
    {
      streamID1 = cdo_open_read(fileIdx);
      auto vlistID1 = cdo_stream_inq_vlist(streamID1);
      auto taxisID1 = vlistInqTaxis(vlistID1);

      VarList varList1(vlistID1);

      if (fileIdx == 0)
      {
        vlistID2 = vlistDuplicate(vlistID1);
        taxisID2 = taxisDuplicate(taxisID1);
        if (taxisHasBounds(taxisID2))
        {
          cdo_warning("Time bounds unsupported by this operator, removed!");
          taxisDeleteBounds(taxisID2);
        }
      }
      else { vlist_compare(vlistID2, vlistID1, CmpVarList::All); }

      numVars = varList1.numVars();

      int tsID = 0;
      while (true)
      {
        auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
        if (numFields == 0) break;

        if (xtsID >= nalloc)
        {
          constexpr int NALLOC_INC = 1024;
          nalloc += NALLOC_INC;
          vDateTimes.resize(nalloc);
          varsData.resize(nalloc);
        }

        vDateTimes[xtsID] = taxisInqVdatetime(taxisID1);

        field2D_init(varsData[xtsID], varList1);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          auto &field = varsData[xtsID][varID][levelID];
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
        }

        tsID++;
        xtsID++;
      }

      cdo_stream_close(streamID1);
    }

    int nts = xtsID;

    std::vector<TimeInfo> timeinfo(nts);
    auto calendar = taxisInqCalendar(taxisID2);

    for (int tsID = 0; tsID < nts; ++tsID)
    {
      auto julday = date_to_julday(calendar, cdiDate_get(vDateTimes[tsID].date));
      auto secofday = time_to_sec(cdiTime_get(vDateTimes[tsID].time));
      double jdatetime = julday + secofday / 86400.0;
      timeinfo[tsID].index = tsID;
      timeinfo[tsID].datetime = jdatetime;
    }

    std::ranges::stable_sort(timeinfo, {}, &TimeInfo::datetime);

    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(numFiles);
    cdo_def_vlist(streamID2, vlistID2);

    int tsID2 = 0;
    for (int tsID = 0; tsID < nts; ++tsID)
    {
      xtsID = timeinfo[tsID].index;

      if (tsID > 0 && is_equal(timeinfo[tsID].datetime, timeinfo[lasttsID].datetime))
      {
        if (skipSameTime)
        {
          if (Options::cdoVerbose)
            cdo_print("Timestep %4d %s already exists, skipped!", xtsID + 1, datetime_to_string(vDateTimes[xtsID]));
          continue;
        }

        if (unique)
        {
          auto lskip = false;
          auto xtsID2 = timeinfo[lasttsID].index;
          auto const &field1 = varsData[xtsID][0][0];
          auto const &field2 = varsData[xtsID2][0][0];
          if (field1.memType == MemType::Float)
          {
            if (field1.vec_f == field2.vec_f) lskip = true;
          }
          else
          {
            if (field1.vec_d == field2.vec_d) lskip = true;
          }

          if (lskip)
          {
            if (Options::cdoVerbose)
              cdo_print("Timestep %4d %s already exists with the same data, skipped!", xtsID + 1,
                        datetime_to_string(vDateTimes[xtsID]));
            continue;
          }
        }
      }

      lasttsID = tsID;

      taxisDefVdatetime(taxisID2, vDateTimes[xtsID]);
      cdo_def_timestep(streamID2, tsID2++);

      VarList varList2(vlistID2);

      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
        {
          auto &field = varsData[xtsID][varID][levelID];
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
  }
};
