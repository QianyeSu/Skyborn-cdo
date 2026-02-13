/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Setbox     setclonlatbox   Set lon/lat box to constant
      Setbox     setcindexbox    Set index box to constant
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "selboxinfo.h"

static void
setcbox(double constant, double *array, int gridID, const SelboxInfo &selboxInfo)
{
  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  for (long ilat = 0; ilat < nlat; ilat++)
    for (long ilon = 0; ilon < nlon; ilon++)
      if ((lat1 <= ilat && ilat <= lat2 && ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))))
      {
        array[nlon * ilat + ilon] = constant;
      }
}

static int
get_gridID(int vlistID1, bool operIndexBox)
{
  std::vector<int> gridsFound;

  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    if (gridInqSize(gridID1) == 1) continue;

    auto gridtype = gridInqType(gridID1);
    auto projtype = gridInqProjType(gridID1);

    auto isReg2dGeoGrid = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR);
    auto projHasGeoCoords = (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_RLL);

    if (isReg2dGeoGrid || projHasGeoCoords || (operIndexBox && (gridtype == GRID_GENERIC || gridtype == GRID_PROJECTION)))
    {
      gridsFound.push_back(gridID1);
    }
    else
    {
      if (gridInqSize(gridID1) > 2) cdo_warning("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
  }

  if (gridsFound.size() == 0) cdo_abort("No processable grid found!");
  if (gridsFound.size() > 1) cdo_abort("Too many different grids!");

  auto gridID = gridsFound[0];
  return gridID;
}

static std::vector<bool>
get_processVars(VarList const &varList, int gridID)
{
  auto numVars = varList.numVars();

  std::vector<bool> processVars(numVars, false);

  int varID;
  for (varID = 0; varID < numVars; ++varID)
    if (gridID == varList.vars[varID].gridID) processVars[varID] = true;

  for (varID = 0; varID < numVars; ++varID)
    if (processVars[varID]) break;

  if (varID >= numVars) cdo_abort("No processable variable found!");

  return processVars;
}

class Setbox : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Setbox",
    .operators
    = { { "setclonlatbox", 0, 0, "constant, western and eastern longitude and southern and northern latitude", SetboxHelp },
        { "setcindexbox", 0, 0, "constant, index of first and last longitude and index of first and last latitude", SetboxHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Setbox> registration = RegisterEntry<Setbox>();

  int SETCLONLATBOX{}, SETCINDEXBOX{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  int gridID{};

  double constant{};

  std::vector<bool> processVars;
  SelboxInfo selboxInfo;

public:
  void
  init() override
  {
    SETCLONLATBOX = module.get_id("setclonlatbox");
    SETCINDEXBOX = module.get_id("setcindexbox");

    (void) SETCLONLATBOX;

    auto operatorID = cdo_operator_id();
    auto operIndexBox = (operatorID == SETCINDEXBOX);

    operator_input_arg(cdo_operator_enter(operatorID));

    constant = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    gridID = get_gridID(vlistID1, operIndexBox);

    processVars = get_processVars(varList1, gridID);

    operator_input_arg(cdo_operator_enter(operatorID));

    selboxInfo = operIndexBox ? gen_index_selbox(1, gridID) : gen_lonlat_selbox(1, gridID);

    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsize = gridInqSize(gridID);
    Varray<double> array(gridsize);

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

        if (processVars[varID])
        {
          size_t numMissVals;
          cdo_read_field(streamID1, array.data(), &numMissVals);

          setcbox(constant, array.data(), gridID, selboxInfo);

          auto missval = varList1.vars[varID].missval;
          numMissVals = varray_num_mv(gridsize, array, missval);
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, array.data(), numMissVals);
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
  }
};
