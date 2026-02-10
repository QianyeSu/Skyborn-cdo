/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "util_string.h"
#include "cdo_options.h"
#include "chunkspec.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "mpim_grid.h"

namespace
{
struct Parameter
{
  std::string filename{};
};
}  // namespace

static Parameter
get_parameter()
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if (key == "filename")  params.filename = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static std::vector<std::string>
get_vars_chunkspec(VarList const &varList, std::string const &filename)
{
  auto numVars = varList.numVars();
  std::vector<std::string> varsChunkSpec(numVars);

  if (filename.size())
  {
    auto fp = std::fopen(filename.c_str(), "r");
    if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
    PMList pmlist;
    pmlist.read_namelist(fp, filename.c_str());
    auto &kvlist = pmlist.front();
    std::fclose(fp);
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

      for (auto const &var : varList.vars)
      {
        if (key == var.name)
        {
          auto const &value = kv.values[0];
          auto chunkSpec = parameter_to_word(value);
          varsChunkSpec[var.ID] = chunkSpec;
        }
      }
    }
  }

  return varsChunkSpec;
}

static void
set_key_chunkspec(int maxSteps, CdoVar const &var, int vlistID, int varID, std::string const &chunkSpecString)
{
  auto chunkSpecIn = cdo::get_chunkspec(vlistID, varID);
  auto chunkSpec = cdo::parse_chunkspec_parameter(chunkSpecString);

  if (is_unstruct_grid(var.gridID))
  {
    if (chunkSpec.y) cdo_abort("%s: chunkSpec of y=%d not available for unstructured grids", var.name, chunkSpec.z);
    if (chunkSpec.x > 0 && (SizeType) chunkSpec.x > var.gridsize)
      cdo_abort("%s: chunkSpec of x=%d is greater than gridSize=%d", var.name, chunkSpec.x, var.gridsize);
  }
  else
  {
    int xsize = gridInqXsize(var.gridID);
    int ysize = gridInqYsize(var.gridID);
    if (chunkSpec.x > 0 && chunkSpec.x > xsize)
      cdo_error("%s: chunkSpec of x=%d is greater than xsize=%d", var.name, chunkSpec.x, xsize);
    if (chunkSpec.y > 0 && chunkSpec.y > ysize)
      cdo_error("%s: chunkSpec of y=%d is greater than xsize=%d", var.name, chunkSpec.y, ysize);
  }

  if (chunkSpec.z > 0 && chunkSpec.z > var.nlevels)
    cdo_error("%s: chunkSpec of z=%d is greater than numLevels=%d", var.name, chunkSpec.z, var.nlevels);

  if (chunkSpec.t > 0 && chunkSpec.t > maxSteps)
    cdo_warning("%s: chunkSpec of t=%d is greater than numSteps=%d", var.name, chunkSpec.t, maxSteps);

  if (chunkSpec.t) cdiDefKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE_DIMT, chunkSpec.t);
  if (chunkSpec.z) cdiDefKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE_DIMZ, chunkSpec.z);
  if (chunkSpec.y) cdiDefKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE_DIMY, chunkSpec.y);
  if (chunkSpec.x) cdiDefKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE_DIMX, chunkSpec.x);

  auto chunkSpecOut = cdo::get_chunkspec(vlistID, varID);
  if (Options::cdoVerbose)
    cdo_print("%s: chunkSpecIn x=%d y=%d z=%d t=%d   chunkSpecOut x=%d y=%d z=%d t=%d", var.name.c_str(), chunkSpecIn.x,
              chunkSpecIn.y, chunkSpecIn.z, chunkSpecIn.t, chunkSpecOut.x, chunkSpecOut.y, chunkSpecOut.z, chunkSpecOut.t);
}

class Setchunkspec : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setchunkspec",
    .operators = { { "setchunkspec", SetchunkspecHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setchunkspec> registration = RegisterEntry<Setchunkspec>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1;

public:
  void
  init() override
  {
    auto params = get_parameter();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    varList1 = VarList(vlistID1);

    auto varsChunkSpec = get_vars_chunkspec(varList1, params.filename);

    vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto numSteps = varList1.numSteps();
    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      if (varsChunkSpec[varID].size() > 0)
      {
        auto chunkSpecString = string_to_lower(varsChunkSpec[varID]);
        set_key_chunkspec(numSteps, varList1.vars[varID], vlistID2, varID, chunkSpecString);
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        if (numSteps > 0)
        {
          auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
          progress.update(fstatus);
        }

        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        auto const &var = varList1.vars[varID];
        field.init(var);
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

    vlistDestroy(vlistID2);
  }
};
