/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_zaxis.h"

class Gengrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Gengrid",
    .operators = { { "gengrid" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Gengrid> registration = RegisterEntry<Gengrid>(module);

  CdoStreamID streamID3{};
  size_t gridsize{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    auto streamID1 = cdo_open_read(0);
    auto streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    auto gridID1 = vlistGrid(vlistID1, 0);
    auto gridID2 = vlistGrid(vlistID2, 0);

    if (gridInqSize(gridID1) != gridInqSize(gridID2)) cdo_abort("Arrays have different grid size!");

    gridsize = gridInqSize(gridID1);
    auto xsize = gridInqXsize(gridID1);
    auto ysize = gridInqYsize(gridID1);

    Varray<double> array1(gridsize), array2(gridsize);

    cdo_stream_inq_timestep(streamID1, 0);
    cdo_stream_inq_timestep(streamID2, 0);

    size_t numMissVals1, numMissVals2;
    (void) cdo_inq_field(streamID1);
    cdo_read_field(streamID1, array1.data(), &numMissVals1);
    (void) cdo_inq_field(streamID2);
    cdo_read_field(streamID2, array2.data(), &numMissVals2);

    auto datatype = vlistInqVarDatatype(vlistID1, 0);

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    if (numMissVals1 || numMissVals2) cdo_abort("Missing values unsupported!");

    auto gridType = (ysize == 0) ? GRID_UNSTRUCTURED : GRID_CURVILINEAR;
    auto gridID3 = gridCreate(gridType, gridsize);

    if (Options::cdoVerbose) cdo_print("xsize %zu  ysize %zu", xsize, ysize);
    if (gridType == GRID_CURVILINEAR && xsize * ysize != gridsize) cdo_abort("xsize*ysize != gridsize");

    if (gridType == GRID_CURVILINEAR) gridDefXsize(gridID3, xsize);
    if (gridType == GRID_CURVILINEAR) gridDefYsize(gridID3, ysize);
    gridDefXvals(gridID3, array1.data());
    gridDefYvals(gridID3, array2.data());

    cdiDefKeyInt(gridID3, CDI_GLOBAL, CDI_KEY_DATATYPE, (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32);

    auto xmm = varray_min_max(array1);
    auto ymm = varray_min_max(array2);

    if (Options::cdoVerbose)
      cdo_print("xminval = %g, xmaxval = %g, yminval = %g, ymaxval = %g", xmm.min, xmm.max, ymm.min, ymm.max);

    // check units
    if (xmm.min > -4 && xmm.max < 8 && ymm.min > -2 && ymm.max < 2)
    {
      cdiDefKeyString(gridID3, CDI_XAXIS, CDI_KEY_UNITS, "radians");
      cdiDefKeyString(gridID3, CDI_YAXIS, CDI_KEY_UNITS, "radians");
    }
    else if (xmm.min > -181 && xmm.max < 361 && ymm.min > -91 && ymm.max < 91)
    {
      // default is degrees
    }
    else { cdo_abort("Units undefined!"); }

    auto zaxisID3 = zaxis_from_name("surface");

    auto vlistID3 = vlistCreate();
    vlistDefVar(vlistID3, gridID3, zaxisID3, TIME_CONSTANT);
    cdiDefKeyString(vlistID3, 0, CDI_KEY_NAME, "dummy");
    vlistDefVarDatatype(vlistID3, 0, CDI_DATATYPE_INT8);

    auto taxisID3 = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID3 = cdo_open_write(2);

    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    int tsID = 0;
    cdo_def_timestep(streamID3, tsID);

    cdo_def_field(streamID3, 0, 0);
    Varray<double> array3(gridsize, 0);
    cdo_write_field(streamID3, array3.data(), gridsize);
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
  }
};
