/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Rotuv      rotuvb          Backward rotation
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>

static void
rot_uv_back(int gridID, Varray<double> &us, Varray<double> &vs)
{
  double xpole = 0, ypole = 0, angle = 0;
  if (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL)
    gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);

  Varray<double> xvals(nlon), yvals(nlat);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, 1, &angle, "angle");
  cdo_grid_to_degree(gridID, CDI_XAXIS, 1, &xpole, "xpole");
  cdo_grid_to_degree(gridID, CDI_YAXIS, 1, &ypole, "ypole");
  cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

  if (xpole > 180) xpole -= 360;
  if (angle > 180) angle -= 360;

  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++)
    {
      auto i = ilat * nlon + ilon;
      auto xval = lamrot_to_lam(yvals[ilat], xvals[ilon], ypole, xpole, angle);
      auto yval = phirot_to_phi(yvals[ilat], xvals[ilon], ypole, angle);
      usvs_to_uv(us[i], vs[i], yval, xval, ypole, xpole, &us[i], &vs[i]);
    }
}

#define MAXARG 16384

class Rotuv : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Rotuv",
    .operators = { { "rotuvb", RotuvHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Rotuv> registration = RegisterEntry<Rotuv>();
  int chcodes[MAXARG];
  const char *chvars[MAXARG];

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int nch{};
  int numVars{};

  bool lvar{};

  VarList varList1{};
  Varray3D<double> varsData;
  std::vector<std::vector<size_t>> varnumMissVals;

public:
  void
  init() override
  {
    operator_input_arg("pairs of u and v in the rotated system");

    nch = cdo_operator_argc();
    if (nch % 2) cdo_abort("Odd number of input arguments!");

    lvar = false;  // We have a list of codes
    int len = (int) cdo_operator_argv(0).size();
    int ix = (cdo_operator_argv(0)[0] == '-') ? 1 : 0;
    for (int i = ix; i < len; ++i)
      if (!std::isdigit(cdo_operator_argv(0)[i]))
      {
        lvar = true;  // We have a list of variables
        break;
      }

    if (lvar)
    {
      for (int i = 0; i < nch; ++i) chvars[i] = cdo_operator_argv(i).c_str();
    }
    else
    {
      for (int i = 0; i < nch; ++i) chcodes[i] = parameter_to_int(cdo_operator_argv(i));
    }

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();

    varnumMissVals.resize(numVars);
    varsData.resize(numVars);

    bool lfound[MAXARG];
    for (int i = 0; i < nch; ++i) lfound[i] = false;

    if (lvar)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        for (int i = 0; i < nch; ++i)
          if (varList1.vars[varID].name == chvars[i]) lfound[i] = true;
      }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Variable %s not found!", chvars[i]);
    }
    else
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto code = varList1.vars[varID].code;
        for (int i = 0; i < nch; ++i)
          if (code == chcodes[i]) lfound[i] = true;
      }
      for (int i = 0; i < nch; ++i)
        if (!lfound[i]) cdo_abort("Code %d not found!", chcodes[i]);
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridID = varList1.vars[varID].gridID;
      if (!(gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL))
        cdo_abort("Only rotated lon/lat grids supported!");

      auto gridsize = gridInqSize(gridID);
      auto nlevels = varList1.vars[varID].nlevels;
      varnumMissVals[varID].resize(nlevels);
      varsData[varID].resize(nlevels);
      for (int levelID = 0; levelID < nlevels; ++levelID) varsData[varID][levelID].resize(gridsize);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

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

        fieldInfoList[fieldID].set(varID, levelID);

        cdo_read_field(streamID1, varsData[varID][levelID].data(), &varnumMissVals[varID][levelID]);
        if (varnumMissVals[varID][levelID]) cdo_abort("Missing values unsupported for this operator!");
      }

      for (int i = 0; i < nch; i += 2)
      {
        int varID;
        for (varID = 0; varID < numVars; ++varID)
        {
          if (lvar)
          {
            if (varList1.vars[varID].name == chvars[i]) break;
          }
          else
          {
            if (varList1.vars[varID].code == chcodes[i]) break;
          }
        }

        if (varID == numVars) cdo_abort("u-wind not found!");

        auto usvarID = varID;

        for (varID = 0; varID < numVars; ++varID)
        {
          if (lvar)
          {
            if (varList1.vars[varID].name == chvars[i + 1]) break;
          }
          else
          {
            if (varList1.vars[varID].code == chcodes[i + 1]) break;
          }
        }

        if (varID == numVars) cdo_abort("v-wind not found!");

        auto vsvarID = varID;
        auto const &usVar = varList1.vars[usvarID];
        auto const &vsVar = varList1.vars[vsvarID];

        if (Options::cdoVerbose)
        {
          if (lvar)
            cdo_print("Using var %s [%s](u) and var %s [%s](v)", usVar.name, chvars[i], vsVar.name, chvars[i + 1]);
          else
            cdo_print("Using code %d [%d](u) and code %d [%d](v)", usVar.code, chcodes[i], vsVar.code, chcodes[i + 1]);
        }

        auto gridID = usVar.gridID;
        auto nlevels1 = usVar.nlevels;
        auto nlevels2 = vsVar.nlevels;
        if (nlevels1 != nlevels2) cdo_abort("u-wind and v-wind have different number of levels!");

        for (int levelID = 0; levelID < nlevels1; ++levelID)
          rot_uv_back(gridID, varsData[usvarID][levelID], varsData[vsvarID][levelID]);
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = fieldInfoList[fieldID].get();
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, varsData[varID][levelID].data(), varnumMissVals[varID][levelID]);
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
