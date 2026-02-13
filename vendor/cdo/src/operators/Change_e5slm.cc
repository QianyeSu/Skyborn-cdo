/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Change_e5slm      change_e5slm          Change ECHAM5 sea land mask
*/

#include <cdi.h>

#include "process_int.h"
#include "cdi_lockedIO.h"

static void
set_code(std::string const &varname, int &code)
{
  // clang-format off
  if      (varname == "SLM")       code = 172;
  else if (varname == "ALAKE")     code = 99;
  else if (varname == "WS")        code = 140;
  else if (varname == "AZ0")       code = 173;
  else if (varname == "ALB")       code = 174;
  else if (varname == "VGRAT")     code = 198;
  else if (varname == "FOREST")    code = 212;
  else if (varname == "FAO")       code = 226;
  else if (varname == "WSMX")      code = 229;
  else if (varname == "GLAC")      code = 232;
  else if (varname == "VLTCLIM")   code = 71;
  else if (varname == "VGRATCLIM") code = 70;
  // clang-format on
}

class Change_e5slm : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Change_e5slm",
    .operators = { { "change_e5slm" }, { "change_e5lsm" }, { "change_e5mask" } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Change_e5slm> registration = RegisterEntry<Change_e5slm>();

private:
  int numFields{};
  size_t numMissVals{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  size_t gridsize{};

  Varray<double> array;
  Varray<double> cland;
  Varray<bool> lsea;

  std::vector<short> codes;

public:
  void
  init() override
  {
    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    VarList varList1(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);

    auto vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    /* get filename of SLM */
    operator_input_arg("filename of the sea land mask");
    operator_check_argc(1);
    const char *fn_slm = cdo_operator_argv(0).c_str();

    /* read SLM */
    auto streamIDslm = stream_open_read_locked(fn_slm);
    auto vlistIDslm = streamInqVlist(streamIDslm);

    {
      VarList varListSLM(vlistIDslm);
      gridsize = varListSLM.vars[0].gridsize;
    }

    array.resize(gridsize);
    cland.resize(gridsize);
    lsea.resize(gridsize);

    streamInqTimestep(streamIDslm, 0);

    {
      int varID, levelID;
      streamInqField(streamIDslm, &varID, &levelID);
      streamReadField(streamIDslm, cland.data(), &numMissVals);
    }

    if (numMissVals) cdo_abort("SLM with missing values are unsupported!");

    auto mm = varray_min_max(cland);
    if (mm.min < 0 || mm.max > 1) cdo_warning("Values of SLM out of bounds! (minval=%g, maxval=%g)", mm.min, mm.max);

    streamClose(streamIDslm);

    for (size_t i = 0; i < gridsize; ++i) lsea[i] = cland[i] <= 0;

    auto numVars = varList1.numVars();
    codes.resize(numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      if (gridsize != var.gridsize) cdo_abort("gridsize differ!");

      auto code = var.code;
      if (code < 0) { set_code(var.name, code); }
      codes[varID] = code;
    }
  }

  void
  run() override
  {
    int tsID = 0;
    while ((numFields = cdo_stream_inq_timestep(streamID1, tsID)))
    {
      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, array.data(), &numMissVals);

        auto code = codes[varID];
        if (code == 172)
        {
          cdo_print("SLM changed!");
          for (size_t i = 0; i < gridsize; ++i) array[i] = cland[i];
        }
        else if (code == 99)
        {
          cdo_print("ALAKE set all values to zero!");
          for (size_t i = 0; i < gridsize; ++i) array[i] = 0;
        }
        else if (code == 232)
        {
          cdo_print("GLAC set sea points to %g!", array[0]);
          for (size_t i = 0; i < gridsize; ++i)
            if (cland[i] < 0.5) array[i] = array[0];
        }
        else if (code == 70 || code == 71 || code == 140 || code == 173 || code == 174 || code == 198 || code == 200 || code == 212
                 || code == 226 || code == 229)
        {
          cdo_print("Code %d set sea points to %g!", code, array[0]);
          for (size_t i = 0; i < gridsize; ++i)
            if (lsea[i]) array[i] = array[0];
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, array.data(), numMissVals);
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
