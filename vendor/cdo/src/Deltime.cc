/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "util_string.h"
#include "cdo_options.h"

constexpr char const *cmons[] = { "", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec" };

class Deltime : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Deltime",
    .operators = { { "delday" }, { "del29feb" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Deltime> registration = RegisterEntry<Deltime>(module);

private:
  int DELDAY{}, DEL29FEB{};
  int dday{}, dmon{};

  CdoStreamID streamID1{};
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2{};
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int nfound = 0;
  bool dataIsUnchanged{};

  VarList varList1{};

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    DELDAY = module.get_id("delday");
    DEL29FEB = module.get_id("del29feb");

    (void) (DELDAY);  // unused

    auto operatorID = cdo_operator_id();

    if (operatorID == DEL29FEB)
    {
      dday = 29;
      dmon = 2;
      operator_check_argc(0);
    }
    else
    {
      auto nsel = cdo_operator_argc();
      operator_check_argc(1);
      auto sarg = cdo_operator_argv(0).c_str();
      dday = atoi(sarg);
      dmon = 0;
      while (std::isdigit(*sarg)) sarg++;
      if (std::isalpha(*sarg))
      {
        char smon[32];
        std::strncpy(smon, sarg, sizeof(smon) - 1);
        smon[sizeof(smon) - 1] = 0;
        cstr_to_lower(smon);
        int im = 0;
        for (; im < 12; ++im)
          if (std::memcmp(smon, cmons[im + 1], 3) == 0) break;

        if (im < 12) dmon = im + 1;
      }
    }

    if (Options::cdoVerbose) cdo_print("delete day %d%s", dday, cmons[dmon]);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisDefCalendar(taxisID2, CALENDAR_365DAYS);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      int year, month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);

      auto copyTimestep = true;
      if (day == dday && (month == dmon || dmon == 0))
      {
        nfound++;
        copyTimestep = false;
        if (Options::cdoVerbose) cdo_print("Delete %4.4d-%2.2d-%2.2d at timestep %d", year, month, day, tsID + 1);
      }

      if (copyTimestep)
      {
        cdo_taxis_copy_timestep(taxisID2, taxisID1);
        cdo_def_timestep(streamID2, tsID2++);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          cdo_def_field(streamID2, varID, levelID);
          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            field.init(varList1.vars[varID]);
            cdo_read_field(streamID1, field);
            cdo_write_field(streamID2, field);
          }
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (nfound == 0) cdo_warning("Day %d%s not found!", dday, cmons[dmon]);

    vlistDestroy(vlistID2);
  }
};
