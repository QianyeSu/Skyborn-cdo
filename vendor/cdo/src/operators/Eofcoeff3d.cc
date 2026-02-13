/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge

*/

/*
  This module contains the following operators:

  Eofcoeff3d             eofcoeff3d             process eof coefficients
*/

#include "cdi.h"

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "util_files.h"
#include "util_string.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

class Eofcoeff3d : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Eofcoeff3d",
    .operators = { { "eofcoeff3d", EofcoeffHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Eofcoeff3d> registration = RegisterEntry<Eofcoeff3d>();

  double missval1 = -999, missval2 = -999;
  int numFields{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;

  int gridID1{};

  int taxisID2{ CDI_UNDEFID };

  size_t gridsizeMax{};
  int numVars{};
  int numLevels{};
  int taxisID3{};

  std::string fileSuffix;

  std::vector<CdoStreamID> streamIDs;

public:
  void
  init() override
  {
    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = cdo_stream_inq_vlist(streamID2);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    // taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID2);

    gridID1 = varList1.vars[0].gridID;
    int gridID2 = varList2.vars[0].gridID;

    gridsizeMax = varList1.gridsizeMax();
    if (gridsizeMax != varList2.gridsizeMax()) cdo_abort("Gridsize of input files does not match!");

    if (vlistNumGrids(vlistID2) > 1 || vlistNumGrids(vlistID1) > 1) cdo_abort("Too many different grids in input");

    numVars = (varList1.numVars() == varList2.numVars()) ? varList1.numVars() : -1;
    numLevels = varList1.vars[0].nlevels;

    cdo_compare_grids(gridID1, gridID2);

    fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));
  }

  void
  run() override
  {
    FieldVector3D eof(numVars);
    for (int varID = 0; varID < numVars; ++varID) eof[varID].resize(numLevels);

    int neof = 0;

    {
      int eofID = 0;
      while (1)
      {
        numFields = cdo_stream_inq_timestep(streamID1, eofID);
        if (numFields == 0) break;

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID1);
          missval1 = varList1.vars[varID].missval;
          eof[varID][levelID].resize((eofID == 0) ? 1 : eofID + 1);
          eof[varID][levelID][eofID].grid = gridID1;
          eof[varID][levelID][eofID].missval = missval1;
          eof[varID][levelID][eofID].resize(gridsizeMax);
          std::ranges::fill(eof[varID][levelID][eofID].vec_d, missval1);

          if (varID >= numVars) cdo_abort("Internal error - varID >= nvars");
          if (levelID >= numLevels) cdo_abort("Internal error - levelID >= nlevs");

          cdo_read_field(streamID1, eof[varID][levelID][eofID]);
        }
        eofID++;
      }

      neof = eofID;
    }

    if (Options::cdoVerbose) cdo_print("%s contains %i eof's", cdo_get_stream_name(0), neof);
    // Create 1x1 Grid for output
    auto gridID3 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID3, 1);
    gridDefYsize(gridID3, 1);
    double xvals = 0, yvals = 0;
    gridDefXvals(gridID3, &xvals);
    gridDefYvals(gridID3, &yvals);

    double zvals = 0.;
    auto zaxisID3 = zaxisCreate(ZAXIS_GENERIC, 1);
    zaxisDefLevels(zaxisID3, &zvals);
    cdiDefKeyString(zaxisID3, CDI_GLOBAL, CDI_KEY_NAME, "zaxis_Reduced");
    cdiDefKeyString(zaxisID3, CDI_GLOBAL, CDI_KEY_LONGNAME,
                    "Reduced zaxis from EOF3D - only one coefficient per 3D eigenvector and time step");

    // Create var-list and time-axis for output

    auto vlistID3 = vlistDuplicate(vlistID2);

    auto numGrids = vlistNumGrids(vlistID3);
    for (int i = 0; i < numGrids; ++i) vlistChangeGridIndex(vlistID3, i, gridID3);

    auto numZaxes = vlistNumZaxis(vlistID3);
    for (int i = 0; i < numZaxes; ++i) vlistChangeZaxisIndex(vlistID3, i, zaxisID3);

    vlistDefTaxis(vlistID3, taxisID3);
    for (int varID = 0; varID < numVars; ++varID) vlistDefVarTimetype(vlistID3, varID, TIME_VARYING);

    // open streams for eofcoeff output

    streamIDs.resize(neof);
    for (int eofID = 0; eofID < neof; eofID++)
    {
      auto fileName = cdo_get_obase() + string_format("%5.5i", eofID);
      if (fileSuffix.size() > 0) fileName += fileSuffix;

      streamIDs[eofID] = open_write(fileName);
      cdo_def_vlist(streamIDs[eofID], vlistID3);
    }

    // ALLOCATE temporary fields for data read and write
    Field in;
    in.resize(gridsizeMax);
    in.grid = gridID1;
    FieldVector2D out(numVars);
    for (int varID = 0; varID < numVars; ++varID)
    {
      out[varID].resize(neof);
      for (int eofID = 0; eofID < neof; eofID++)
      {
        out[varID][eofID].missval = missval1;
        out[varID][eofID].numMissVals = 0;
        out[varID][eofID].resize(1);
      }
    }

    int tsID = 0;
    while (1)
    {
      numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      for (int varID = 0; varID < numVars; ++varID)
        for (int eofID = 0; eofID < neof; eofID++)
        {
          out[varID][eofID].vec_d[0] = 0;
          out[varID][eofID].grid = gridID3;
          out[varID][eofID].missval = missval2;
        }

      cdo_taxis_copy_timestep(taxisID3, taxisID2);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID2);
        cdo_read_field(streamID2, in);
        missval2 = varList2.vars[varID].missval;

        for (int eofID = 0; eofID < neof; eofID++)
        {
          if (fieldID == 0) cdo_def_timestep(streamIDs[eofID], tsID);

          // size_t numMissVals = 0;
          for (size_t i = 0; i < gridsizeMax; ++i)
          {
            if (!fp_is_equal(in.vec_d[i], missval2) && !fp_is_equal(eof[varID][levelID][eofID].vec_d[i], missval1))
            {
              // out[varID][eofID].vec_d[0] += w[i]*in.vec_d[i]*eof[varID][levelID][eofID].vec_d[i];
              out[varID][eofID].vec_d[0] += in.vec_d[i] * eof[varID][levelID][eofID].vec_d[i];
            }
            // else { numMissVals += 1; }
          }
          /*
          if ( numMissVals ) {
            out[varID][eofID].numMissVals=1;
            out[varID][eofID].vec_d[0] = missval2;
          }
          */
        }

        if (varID >= numVars) cdo_abort("Internal error - varID >= nvars");
        if (levelID >= numLevels) cdo_abort("Internal error - levelID >= nlevs");
      }

      for (int eofID = 0; eofID < neof; eofID++)
      {
        for (int varID = 0; varID < numVars; ++varID)
        {
          cdo_def_field(streamIDs[eofID], varID, 0);
          cdo_write_field(streamIDs[eofID], out[varID][eofID]);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    for (auto const &streamID : streamIDs) cdo_stream_close(streamID);

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
