/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author:

*/

#include <fstream>

#include <cdi.h>

#include "process_int.h"
#include "cdo_zaxis.h"

class Setrcaname : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setrcaname",
    .operators = { { "setrcaname" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Setrcaname>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int nvars{};
  bool dataIsUnchanged{};

  VarList varList1{};

  void
  read_rca(std::string const &filename, int p_nvars, int p_vlistID2)
  {
    std::ifstream file(filename);
    if (!file.is_open()) cdo_abort("Open failed on: %s\n", filename);

    std::string line;
    while (std::getline(file, line))
    {
      char sname[CDI_MAX_NAME], sdescription[CDI_MAX_NAME], sunits[CDI_MAX_NAME];
      int scode, sltype, slevel;
      std::sscanf(line.c_str(), "%d\t%d\t%d\t%s\t%s\t%s", &scode, &sltype, &slevel, sname, sdescription, sunits);
      /*
      printf("%s\n", line);
      printf("%d:%d:%d:%s:%s:%s\n", scode, sltype, slevel, sname,
      sdescription, sunits);
      */
      for (int varID = 0; varID < p_nvars; ++varID)
      {
        auto code = vlistInqVarCode(p_vlistID2, varID);
        auto zaxisID = vlistInqVarZaxis(p_vlistID2, varID);
        auto nlev = zaxisInqSize(zaxisID);

        auto ltype = zaxis_to_ltype(zaxisID);

        if (code == scode)
        {
          if (ltype == 105)
          {
            if (nlev != 1)
            {
              cdo_warning("Number of levels should be 1 for level type 105!");
              cdo_warning("Maybe environment variable SPLIT_LTYPE_105 is not set.");
              continue;
            }
            auto level = (int) cdo_zaxis_inq_level(zaxisID, 0);
            if (sltype == 105 && slevel == level)
            {
              cdiDefKeyString(p_vlistID2, varID, CDI_KEY_NAME, sname);
              cdiDefKeyString(p_vlistID2, varID, CDI_KEY_LONGNAME, sdescription);
              cdiDefKeyString(p_vlistID2, varID, CDI_KEY_UNITS, sunits);
              break;
            }
          }
          else if (sltype != 105)
          {
            cdiDefKeyString(p_vlistID2, varID, CDI_KEY_NAME, sname);
            cdiDefKeyString(p_vlistID2, varID, CDI_KEY_LONGNAME, sdescription);
            cdiDefKeyString(p_vlistID2, varID, CDI_KEY_UNITS, sunits);
            break;
          }
        }
      }
    }

    file.close();
  }

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    operator_input_arg("file name with RCA names");
    auto rcsnames = cdo_operator_argv(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    nvars = vlistNvars(vlistID2);

    read_rca(rcsnames, nvars, vlistID2);

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
    Field field{};

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
        cdo_def_field(streamID2, varID, levelID);

        if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
        else
        {
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID1, field);
          cdo_write_field(streamID2, field);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
