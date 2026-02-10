/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Yseasstat  yseasrange      Multi-year seasonal range
      Yseasstat  yseasmin        Multi-year seasonal minimum
      Yseasstat  yseasmax        Multi-year seasonal maximum
      Yseasstat  yseassum        Multi-year seasonal sum
      Yseasstat  yseasmean       Multi-year seasonal mean
      Yseasstat  yseasavg        Multi-year seasonal average
      Yseasstat  yseasvar        Multi-year seasonal variance
      Yseasstat  yseasvar1       Multi-year seasonal variance [Normalize by (n-1)]
      Yseasstat  yseasstd        Multi-year seasonal standard deviation
      Yseasstat  yseasstd1       Multi-year seasonal standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_stepstat.h"
#include "cdo_season.h"
#include "datetime.h"
#include "process_int.h"
#include "progress.h"
#include "field_functions.h"

class Yseasstat : public Process
{
  enum
  {
    MaxSeasons = 4
  };
  int seas_numSets[MaxSeasons]{};

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Yseasstat",
    .operators = { { "yseasrange", FieldFunc_Range, 0, YseasstatHelp },
                   { "yseasmin", FieldFunc_Min, 0, YseasstatHelp },
                   { "yseasmax", FieldFunc_Max, 0, YseasstatHelp },
                   { "yseassum", FieldFunc_Sum, 0, YseasstatHelp },
                   { "yseasmean", FieldFunc_Mean, 0, YseasstatHelp },
                   { "yseasavg", FieldFunc_Avg, 0, YseasstatHelp },
                   { "yseasstd", FieldFunc_Std, 0, YseasstatHelp },
                   { "yseasstd1", FieldFunc_Std1, 0, YseasstatHelp },
                   { "yseasvar", FieldFunc_Var, 0, YseasstatHelp },
                   { "yseasvar1", FieldFunc_Var1, 0, YseasstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Yseasstat> registration = RegisterEntry<Yseasstat>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  int maxFields{};
  std::vector<FieldInfo> fieldInfoList{};

  cdo::StepStat3D stepStat{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);

    stepStat.init(operfunc);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    if (!stepStat.lminmax) vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    maxFields = varList1.maxFields();
    fieldInfoList = std::vector<FieldInfo>(maxFields);
  }

  void
  run() override
  {
    Field field;
    CdiDateTime vDateTimes[MaxSeasons]{};
    FieldVector2D varsData1[MaxSeasons], varsData2[MaxSeasons], samp1[MaxSeasons];

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    stepStat.set_dimlen0(MaxSeasons);

    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    int otsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (numSteps > 1) progress.update((tsID + 1.0) / numSteps);

      auto vDateTime = taxisInqVdatetime(taxisID1);

      auto season = month_to_season(decode_month(vDateTime.date));

      set_date_time(vDateTimes[season], vDateTime);

      if (!stepStat.var1(season).size()) { stepStat.alloc(season, varList1, VARS_MEMTYPE); }

      auto numSets = seas_numSets[season];
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
        stepStat.add_field(field, season, varID, levelID, numSets);
      }

      seas_numSets[season]++;
      tsID++;
    }

    for (int season = 0; season < MaxSeasons; ++season)
      if (seas_numSets[season])
      {
        auto numSets = seas_numSets[season];

        cdo::fields_process_3D(season, fieldInfoList, varList1, stepStat, numSets);

        taxisDefVdatetime(taxisID2, vDateTimes[season]);
        cdo_def_timestep(streamID2, otsID);

        for (int fieldID = 0; fieldID < maxFields; ++fieldID)
        {
          auto [varID, levelID] = fieldInfoList[fieldID].get();
          if (otsID && varList1.vars[varID].isConstant) continue;

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, stepStat.var1(season, varID, levelID));
        }

        otsID++;
      }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
