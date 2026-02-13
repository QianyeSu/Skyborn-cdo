/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Ralf Müller
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
      Intlevel   intlevelx3d     Linear level interpolation on a 3d vertical coordinates variable with extrapolation
*/

#include <cdi.h>

#include "cdo_omp.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdi_lockedIO.h"

void vert_interp_lev3d(size_t gridsize, int nlev1, double missval, const Field3D &field1, Field3D &field2, int nlev2,
                       Varray<int> const &lev_idx, Varray<float> const &lev_wgt);
void vert_gen_weights(int expol, int nlev1, Varray<double> const &lev1, int nlev2, Varray<double> const &lev2, Varray<int> &lev_idx,
                      Varray<float> &lev_wgt);
bool levelDirUp(int nlev, const double *const lev);
bool levelDirDown(int nlev, const double *const lev);

/*
 * Create weights for the 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d data fields.
 *
 * 3d version of vert_gen_weights() (src/Intlevel.cc)
 */
static void
vert_gen_weights3d(bool expol, size_t gridsize, int nlev1, Varray<float> const &xlev1, int nlev2, Varray<float> const &xlev2,
                   Varray<int> &xlev_idx, Varray<float> &xlev_wgt)
{
  auto nthreads = Threading::ompNumMaxThreads;
  Varray2D<double> lev1p2(nthreads, Varray<double>(nlev1 + 2));
  Varray2D<double> lev2(nthreads, Varray<double>(nlev2));
  Varray2D<float> lev_wgt(nthreads, Varray<float>(nlev2));
  Varray2D<int> lev_idx(nthreads, Varray<int>(nlev2));

  // Check monotony of vertical levels
  for (int k = 0; k < nlev1; ++k) lev1p2[0][k] = xlev1[k * gridsize];
  auto lup = levelDirUp(nlev1, lev1p2[0].data());
  auto ldown = levelDirDown(nlev1, lev1p2[0].data());
  if (!lup && !ldown) cdo_abort("Non monotonic zaxis!");
  double level_0 = lup ? -1.e33 : 1.e33;
  double level_N = lup ? 1.e33 : -1.e33;

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto ompthID = cdo_omp_get_thread_num();

    lev1p2[ompthID][0] = level_0;
    lev1p2[ompthID][nlev1 + 1] = level_N;
    for (int k = 0; k < nlev1; ++k) lev1p2[ompthID][k + 1] = xlev1[k * gridsize + i];
    for (int k = 0; k < nlev2; ++k) lev2[ompthID][k] = xlev2[k * gridsize + i];

    vert_gen_weights(expol, nlev1 + 2, lev1p2[ompthID], nlev2, lev2[ompthID], lev_idx[ompthID], lev_wgt[ompthID]);

    for (int k = 0; k < nlev2; ++k) xlev_idx[k * gridsize + i] = lev_idx[ompthID][k];
    for (int k = 0; k < nlev2; ++k) xlev_wgt[k * gridsize + i] = lev_wgt[ompthID][k];
  }
}

static void
vlist_copy_var_attributes(const CdoVar &var0, int vlistID3, int oz3dvarID)
{
  cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_NAME, var0.name.c_str());
  if (var0.longname.size()) cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_LONGNAME, var0.longname.c_str());
  if (var0.units.size()) cdiDefKeyString(vlistID3, oz3dvarID, CDI_KEY_UNITS, var0.units.c_str());
}

static CdoVar
read_source_coordinate(int streamNumber, Varray<float> &zlevelsIn)
{
  auto streamID2 = cdo_open_read(streamNumber);  // 3d vertical source coordinate
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);

  VarList varList(vlistID2);

  if (varList.numVars() != 1) cdo_abort("infile2: Only one single variable is allowed!");

  auto const &var0 = varList.vars[0];

  zlevelsIn.resize(var0.gridsize * var0.nlevels);

  auto numFields = cdo_stream_inq_timestep(streamID2, 0);
  if (Options::cdoVerbose) cdo_print("%d fields input 3d vertical height", numFields);

  while (numFields-- > 0)
  {
    auto [varID, levelID] = cdo_inq_field(streamID2);
    auto offset = var0.gridsize * levelID;
    size_t numMissVals;
    cdo_read_field_f(streamID2, &zlevelsIn[offset], &numMissVals);
    if (0 != numMissVals) cdo_abort("Input vertical coordinate variables are not allowed to contain missing values.");
  }

  cdo_stream_close(streamID2);

  return var0;
}

static CdoVar
read_target_coordinate(std::string const &fileName, Varray<float> &zlevelsOut)
{
  auto streamID0 = stream_open_read_locked(fileName.c_str());  // 3d vertical target coordinate
  auto vlistID0 = streamInqVlist(streamID0);

  VarList varList(vlistID0);

  if (varList.numVars() != 1) cdo_abort("tgtcoordinate: Only one single variable is allowed!");

  auto const &var0 = varList.vars[0];

  zlevelsOut.resize(var0.gridsize * var0.nlevels);

  auto numFields = streamInqTimestep(streamID0, 0);
  if (Options::cdoVerbose) cdo_print("%d fields target 3d vertical height and gridsize %zu", numFields, var0.gridsize);

  while (numFields-- > 0)
  {
    int varID, levelID;
    streamInqField(streamID0, &varID, &levelID);
    auto offset = var0.gridsize * levelID;
    size_t numMissVals;
    streamReadFieldF(streamID0, &zlevelsOut[offset], &numMissVals);
    if (0 != numMissVals) cdo_abort("Output vertical coordinate variables are not allowd to contain missing values.");
  }

  streamClose(streamID0);

  return var0;
}

class Intlevel3d : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intlevel3d",
    .operators = { { "intlevel3d", Intlevel3dHelp }, { "intlevelx3d", Intlevel3dHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Intlevel3d> registration = RegisterEntry<Intlevel3d>();
  int INTLEVEL3D{}, INTLEVELX3D{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};
  int vlistID3{};

  int numVars{};

  int nlevi{};
  int nlevo{};

  int zaxisID1 = -1;

  int oz3dvarID{};

  size_t gridSize;

  MemType memType{};

  VarList varList1{};
  VarList varList3{};

  Varray<int> lev_idx;
  Varray<float> lev_wgt;

  Varray<float> zlevelsOut;

public:
  void
  init() override
  {
    INTLEVEL3D = module.get_id("intlevel3d");
    INTLEVELX3D = module.get_id("intlevelx3d");

    (void) INTLEVEL3D;  // unused

    auto operatorID = cdo_operator_id();
    auto expol = (operatorID == INTLEVELX3D);

    operator_input_arg("tgtcoordinate");

    streamID1 = cdo_open_read(0);  // input data

    // Read filename from Parameter
    operator_input_arg("filename for vertical target coordinates variable");
    operator_check_argc(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    vlistID3 = vlistDuplicate(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID1);

    varList1 = VarList(vlistID1);
    varList_set_unique_memtype(varList1);
    memType = varList1.vars[0].memType;

    // Read 3d source coordinate (streamID2)
    Varray<float> zlevelsIn;
    auto zvari = read_source_coordinate(1, zlevelsIn);
    auto gridsizei = zvari.gridsize;
    nlevi = zvari.nlevels;

    // Read 3d target coordinate (streamID0)
    auto zvaro = read_target_coordinate(cdo_operator_argv(0), zlevelsOut);

    auto gridsizeo = zvaro.gridsize;  // horizontal gridsize of output z coordinate
    nlevo = zvaro.nlevels;            // number of output levels for later use
    auto gridID3 = zvaro.gridID;
    auto zaxisID3 = zvaro.zaxisID;

    // gridsize of input and output vertical coordinate must be equal (later use of gridsizeo ONLY)
    if (gridsizei != gridsizeo) cdo_abort("Input and output vertical coordinate must have the same gridsize!");

    gridSize = gridsizeo;

    /*
     * Check for the correct vertical axis in the input: Variables with the same
     * number of levels as the input vertical levels from operators parameter (streamID0).
     * Variables with a different z-axis should be copied into output.
     */
    auto numZaxes = vlistNumZaxis(vlistID1);
    int i;
    for (i = 0; i < numZaxes; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevel = zaxisInqSize(zaxisID);
      if (nlevel == nlevi)
      {
        zaxisID1 = zaxisID;
        break;
      }
    }
    if (i == numZaxes) cdo_abort("No processable variable found (vertical coordinate differ)!");

    auto numGrids = vlistNumGrids(vlistID1);
    for (i = 0; i < numGrids; ++i)
    {
      auto gridsize = gridInqSize(vlistGrid(vlistID1, i));
      if (gridsize == gridSize) break;
    }
    if (i == numZaxes) cdo_abort("No processable variable found (grid coordinate differ)!");

    // Create weights for later interpolation - assumption: input vertical correct is constant in time
    lev_idx.resize(nlevo * gridSize);
    lev_wgt.resize(nlevo * gridSize);

    vert_gen_weights3d(expol, gridSize, nlevi, zlevelsIn, nlevo, zlevelsOut, lev_idx, lev_wgt);
    varray_free(zlevelsIn);

    for (int index = 0; index < numZaxes; ++index)
      if (zaxisID1 == vlistZaxis(vlistID1, index)) vlistChangeZaxisIndex(vlistID3, index, zaxisID3);

    // add the vertical output field to the output stream
    oz3dvarID = vlistDefVar(vlistID3, gridID3, zaxisID3, TIME_VARYING);
    vlist_copy_var_attributes(zvaro, vlistID3, oz3dvarID);

    streamID3 = cdo_open_write(2);  // output stream
    cdo_def_vlist(streamID3, vlistID3);

    varList3 = VarList(vlistID3);
    varList_set_memtype(varList3, memType);
  }

  void
  run() override
  {
    numVars = varList1.numVars();

    std::vector<bool> processVars(numVars);
    std::vector<bool> interpVars(numVars, false);              // marker for variables to be interpolated
    std::vector<std::vector<size_t>> varnumMissVals(numVars);  // can for missing values of arbitrary variables
    Field3DVector vardata1(numVars);
    Field3DVector vardata2(numVars);

    auto maxLevels = std::max(nlevi, nlevo);
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      auto zaxisID = var1.zaxisID;
      auto gridsize = var1.gridsize;
      auto nlevel = var1.nlevels;

      vardata1[varID].init(var1);

      /*  variabls for interpolation:
       *  * have the required vertical axis, i.e. the correct number of levels (nlevi)
       *  * have the same number of horizontal grid points (i.e. same gridSize) like the two vertical coordinates
       *  * are NOT the output vertical coordinates itself
       */
      interpVars[varID] = (zaxisID == zaxisID1 && varID != oz3dvarID && gridsize == gridSize);
      if (interpVars[varID])
      {
        varnumMissVals[varID].resize(maxLevels, 0);
        vardata2[varID].init(varList3.vars[varID]);
      }
      else
      {
        varnumMissVals[varID].resize(nlevel);
        if (Options::cdoVerbose) cdo_print("Ignore variable %s (levels=%d gridsize=%zu)!", var1.name, nlevel, gridsize);
      }
    }

    {
      int varID;
      for (varID = 0; varID < numVars; ++varID)
        if (interpVars[varID]) break;
      if (varID == numVars) cdo_abort("No processable variable found!");
    }

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (int varID = 0; varID < numVars; ++varID) processVars[varID] = false;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      // Read the whole 3d data field
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_read_field(streamID1, vardata1[varID], levelID, &varnumMissVals[varID][levelID]);
        processVars[varID] = true;
      }

      // Perform the interpolation on all valid data variables
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var1 = varList1.vars[varID];
        if (processVars[varID] && interpVars[varID])
        {
          auto gridsize = var1.gridsize;
          auto missval = var1.missval;

          vert_interp_lev3d(gridsize, nlevi, missval, vardata1[varID], vardata2[varID], nlevo, lev_idx, lev_wgt);

          for (int levelID = 0; levelID < nlevo; ++levelID)
          {
            auto offset = gridsize * levelID;
            auto func = [&](auto const &v) { varnumMissVals[varID][levelID] = array_num_mv(gridsize, &v[offset], missval); };
            field_operation(func, vardata2[varID]);
          }
        }
        else
        {
          if (Options::cdoVerbose && tsID <= 1) cdo_print("Perform no interpolation on variable %s", var1.name);
        }
      }

      // write the output
      for (int varID = 0; varID < numVars; ++varID)
      {
        if (processVars[varID])
        {
          for (int levelID = 0; levelID < varList3.vars[varID].nlevels; ++levelID)
          {
            cdo_def_field(streamID3, varID, levelID);
            cdo_write_field(streamID3, interpVars[varID] ? vardata2[varID] : vardata1[varID], levelID,
                            varnumMissVals[varID][levelID]);
          }
        }
      }

      // copy output z coordinate to output stream
      for (int levelID = 0; levelID < varList3.vars[oz3dvarID].nlevels; ++levelID)
      {
        auto offset = varList3.vars[oz3dvarID].gridsize * levelID;
        cdo_def_field(streamID3, oz3dvarID, levelID);
        cdo_write_field_f(streamID3, &zlevelsOut[offset], 0);
      }

      tsID++;
    }
  }
  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID3);

    vlistDestroy(vlistID3);
  }
};
