/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_string.h"
#include "field_functions.h"

template <typename T>
static void
calc_full_depth(size_t gridsize, size_t nlevels, Varray<T> const &thick_c, Varray<T> const &stretch_c, Varray<T> const &zos,
                Varray<T> &fullDepth)
{
  std::ranges::fill_n(fullDepth.begin(), gridsize, 0.0);

  for (size_t k = 1; k < nlevels; ++k)
  {
    auto depth = &fullDepth[k * gridsize];
    auto depthm1 = &fullDepth[(k - 1) * gridsize];
    auto thickm1 = &thick_c[(k - 1) * gridsize];
    for (size_t i = 0; i < gridsize; ++i) depth[i] = depthm1[i] + stretch_c[i] * thickm1[i];
  }

  for (size_t k = 0; k < nlevels; ++k)
  {
    auto depth = &fullDepth[k * gridsize];
    auto thick = &thick_c[k * gridsize];
    for (size_t i = 0; i < gridsize; ++i) depth[i] += 0.5 * stretch_c[i] * thick[i] - zos[i];
  }
}

static void
calc_full_depth(const Field3D &thick_c, const Field3D &stretch_c, const Field3D &zos, Field3D &fullDepth)
{
  auto gridsize = thick_c.gridsize;
  auto nlevels = thick_c.nlevels;
  auto memType = thick_c.memType;
  if (memType == MemType::Float)
    calc_full_depth(gridsize, nlevels, thick_c.vec_f, stretch_c.vec_f, zos.vec_f, fullDepth.vec_f);
  else
    calc_full_depth(gridsize, nlevels, thick_c.vec_d, stretch_c.vec_d, zos.vec_d, fullDepth.vec_d);
}

class Depth : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Depth",
    .operators = { { "zsdepth" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Depth> registration = RegisterEntry<Depth>(module);

private:
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numLevels{};
  int depthID{};
  int thickID = -1, zosID = -1, stretchID = -1, draftaveID = -1;

  Field3DVector vardata1;
  Field3D fullDepth;

public:
  void
  init() override
  {
    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);

    VarList varList1(vlistID1);
    varList_set_unique_memtype(varList1);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto varname = string_to_lower(varList1.vars[varID].name);

      // clang-format off
      if      (varname == "prism_thick_c") thickID = varID;
      else if (varname == "stretch_c")     stretchID = varID;
      else if (varname == "zos")           zosID = varID;
      else if (varname == "draftave")      draftaveID = varID;
      // clang-format on
    }

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
      if (-1 != thickID)    cdo_print("  %s -> %s", "prism thickness at cells", varList1.vars[thickID].name);
      if (-1 != stretchID)  cdo_print("  %s -> %s", "zstar surface stretch at cell center", varList1.vars[stretchID].name);
      if (-1 != zosID)      cdo_print("  %s -> %s", "zstar sfc elevation at cell center", varList1.vars[zosID].name);
      if (-1 != draftaveID) cdo_print("  %s -> %s", "draftave", varList1.vars[draftaveID].name);
      // clang-format on
    }

    if (-1 == thickID) cdo_abort("prism_thick_c not found!");
    if (-1 == stretchID) cdo_abort("stretch_c not found!");
    if (-1 == zosID) cdo_abort("zos not found!");
    if (-1 == draftaveID) cdo_warning("draftave not found, set to zero!");

    auto zaxisID = varList1.vars[thickID].zaxisID;
    numLevels = varList1.vars[thickID].nlevels;

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids > 1) cdo_abort("Too many different grids!");

    int index = 0;
    auto gridID = vlistGrid(vlistID1, index);

    auto vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());
    vlistDefTaxis(vlistID2, taxisID2);

    depthID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);

    cdiDefKeyString(vlistID2, depthID, CDI_KEY_NAME, "depth_c");
    cdiDefKeyString(vlistID2, depthID, CDI_KEY_STDNAME, "depth");
    cdiDefKeyString(vlistID2, depthID, CDI_KEY_LONGNAME, "depth_below_sea");

    vardata1 = Field3DVector(numVars);

    for (int varID = 0; varID < numVars; ++varID) vardata1[varID].init(varList1.vars[varID]);

    fullDepth.init(varList1.vars[thickID]);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
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
        size_t numMissVals;
        cdo_read_field(streamID1, vardata1[varID], levelID, &numMissVals);
        if (numMissVals) cdo_abort("Missing values unsupported!");
      }

      if (-1 != draftaveID) field2_add(vardata1[zosID], vardata1[draftaveID]);
      calc_full_depth(vardata1[thickID], vardata1[stretchID], vardata1[zosID], fullDepth);

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        cdo_def_field(streamID2, depthID, levelID);
        cdo_write_field(streamID2, fullDepth, levelID, 0);
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
