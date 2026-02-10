/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setzaxis   setzaxis        Set zaxis
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_zaxis.h"

int
getkeyval_dp(std::string const &keyval, std::string const &key, double *val)
{
  int status = 0;

  if (keyval.starts_with(key))
  {
    if (keyval.size() > (key.size() + 2) && keyval[key.size()] == '=')
    {
      *val = parameter_to_double(keyval.substr(key.size() + 1));
      status = 1;
    }
    else { cdo_abort("Syntax error for parameter %s!", keyval); }
  }

  return status;
}

class Setzaxis : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setzaxis",
    .operators = { { "setzaxis", 0, 0, "zaxis description file", SetzaxisHelp }, { "genlevelbounds", SetzaxisHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setzaxis> registration = RegisterEntry<Setzaxis>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

public:
  void
  init() override
  {
    auto SETZAXIS = module.get_id("setzaxis");
    auto GENLEVELBOUNDS = module.get_id("genlevelbounds");

    auto operatorID = cdo_operator_id();

    int zaxisID2 = -1;
    bool hasZtop = false, hasZbot = false;
    double ztop = 0, zbot = 0;

    if (operatorID == SETZAXIS)
    {
      operator_input_arg(cdo_operator_enter(operatorID));
      operator_check_argc(1);
      zaxisID2 = cdo_define_zaxis(cdo_operator_argv(0));
    }
    else if (operatorID == GENLEVELBOUNDS)
    {
      auto numArgs = cdo_operator_argc();
      auto const &argList = cdo_get_oper_argv();
      for (int i = 0; i < numArgs; ++i)
      {
        if (Options::cdoVerbose) cdo_print("keyval[%d]: %s", i + 1, argList[i]);

        if (!hasZbot && getkeyval_dp(argList[i], "zbot", &zbot))
          hasZbot = true;
        else if (!hasZtop && getkeyval_dp(argList[i], "ztop", &ztop))
          hasZtop = true;
        else
          cdo_abort("Parameter >%s< unsupported! Supported parameter are: zbot, ztop", argList[i]);
      }
    }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    if (operatorID == SETZAXIS)
    {
      int found = 0;
      auto numZaxes = vlistNumZaxis(vlistID1);
      for (int index = 0; index < numZaxes; ++index)
      {
        auto zaxisID1 = vlistZaxis(vlistID1, index);
        if (zaxisInqSize(zaxisID1) == zaxisInqSize(zaxisID2))
        {
          vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
          found++;
        }
      }
      if (!found) cdo_warning("No zaxis with %d levels found!", zaxisInqSize(zaxisID2));
    }
    else if (operatorID == GENLEVELBOUNDS)
    {
      auto numZaxes = vlistNumZaxis(vlistID1);
      for (int index = 0; index < numZaxes; ++index)
      {
        auto zaxisID1 = vlistZaxis(vlistID1, index);
        auto nlevels = zaxisInqSize(zaxisID1);
        if (nlevels > 1)
        {
          Varray<double> levels(nlevels), lbounds(nlevels), ubounds(nlevels);

          cdo_zaxis_inq_levels(zaxisID1, levels.data());
          zaxisID2 = zaxisDuplicate(zaxisID1);
          if (!zaxisInqLevels(zaxisID1, nullptr)) zaxisDefLevels(zaxisID2, levels.data());

          gen_layer_bounds(nlevels, levels, lbounds, ubounds);

          auto isPositive = !(levels[0] < 0.0 && levels[nlevels - 1] < 0.0);
          auto isReverse = (levels[0] > levels[nlevels - 1]);
          auto positveIsDown = (isPositive && !isReverse && positive_is_down(zaxisID1));

          if (hasZbot)
          {
            if (positveIsDown)
              ubounds[nlevels - 1] = zbot;
            else
              lbounds[0] = zbot;
          }
          if (hasZtop)
          {
            if (positveIsDown)
              lbounds[0] = ztop;
            else
              ubounds[nlevels - 1] = ztop;
          }
          zaxisDefLbounds(zaxisID2, lbounds.data());
          zaxisDefUbounds(zaxisID2, ubounds.data());
          vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
        }
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;
    VarList varList1(vlistID1);

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

        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
