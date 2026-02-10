/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Remapeta     remapeta          Model to model level interpolation
*/

#include <fstream>

#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "hetaeta.h"
#include "vertical_interp.h"
#include "stdnametable.h"
#include "const.h"
#include "cdo_options.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"

static void
set_missval(Vmask const &imiss, double missval, auto &sarray)
{
  auto n = imiss.size();
  if (n > 0)
  {
    assert(n <= sarray.size());

    for (size_t i = 0; i < n; ++i)
    {
      if (imiss[i]) sarray[i] = missval;
    }
  }
}

static void
correct_humidity(auto q, double q_min)
{
  auto n = q.size();
  for (size_t i = 0; i < n; ++i)
  {
    if (q[i] < q_min) q[i] = q_min;
  }
}

static long
ncctop(double cptop, long nlev, long nlevp1, const double *vct_a, const double *vct_b)
{
  /*
    Description:
    Defines highest level *ncctop* where condensation is allowed.

    Author:

    E. Roeckner, MPI, October 2001
  */
  long nctop = 0;
  Varray<double> zph(nlevp1), zp(nlev);
  // double    cptop  =  1000.;   /* min. pressure level for cond. */

  // half level pressure values, assuming 101320. Pa surface pressure

  for (long jk = 0; jk < nlevp1; ++jk)
  {
    auto za = vct_a[jk];
    auto zb = vct_b[jk];
    zph[jk] = za + zb * 101320.;
  }

  // full level pressure

  for (long jk = 0; jk < nlev; ++jk) zp[jk] = (zph[jk] + zph[jk + 1]) * 0.5;

  // search for pressure level cptop (Pa)

  for (long jk = 0; jk < nlev; ++jk)
  {
    nctop = jk;
    if (zp[jk] >= cptop) break;
  }

  return nctop;
}

static Varray<double>
vct_from_file(std::string const &filename)
{
  constexpr int maxvct = 8192;
  Varray<double> vct;
  vct.resize(maxvct);

  std::ifstream file(filename);
  if (!file.is_open()) cdo_abort("Open failed on: %s\n", filename);

  int i = 0;
  std::string line;
  while (std::getline(file, line))
  {
    if (line[0] == '#' || line[0] == '\0') continue;

    char *lineCpy = strdup(line.c_str());
    char *pline = lineCpy;
    auto num = (int) std::strtod(pline, &pline);
    if (pline == nullptr) cdo_abort("Format error in VCT file %s!", filename);
    if (num != i) cdo_warning("Inconsistent VCT file, entry %d is %d.", i, num);

    if (i + maxvct / 2 >= maxvct - 1) cdo_abort("Too many values in VCT file!");

    vct[i] = std::strtod(pline, &pline);
    if (pline == nullptr) cdo_abort("Format error in VCT file %s!", filename);

    vct[i + maxvct / 2] = std::strtod(pline, &pline);

    std::free(lineCpy);

    i++;
  }

  file.close();

  auto nvct = 2 * i;
  auto nlevh = i - 1;

  for (i = 0; i < nlevh + 1; ++i) vct[i + nvct / 2] = vct[i + maxvct / 2];

  vct.resize(nvct);

  return vct;
}

template <typename T>
static void
vertSum(Varray<double> &sum, Varray<T> const &var3d, size_t gridsize, size_t nlevels)
{
  for (size_t i = 0; i < gridsize; ++i) sum[i] = 0;

  for (size_t k = 0; k < nlevels; ++k)
    for (size_t i = 0; i < gridsize; ++i) { sum[i] += var3d[k * gridsize + i]; }
}

template <typename T>
static void
vertSumw(Varray<double> &sum, Varray<T> const &var3d, size_t gridsize, size_t nlevels, Varray<double> const &deltap)
{
  for (size_t i = 0; i < gridsize; ++i) sum[i] = 0;

  for (size_t k = 0; k < nlevels; ++k)
    for (size_t i = 0; i < gridsize; ++i) { sum[i] += var3d[k * gridsize + i] * deltap[k * gridsize + i]; }
}

template <typename T>
void
field_copy_array(size_t len, Field const &field, T *array)
{
  auto func = [&](auto const &v)
  {
    for (size_t i = 0; i < len; ++i) array[i] = v[i];
  };
  field_operation(func, field);
}

constexpr int MaxVars3D = 1024;

template <typename T>
void
resize(Varray<T> &arr, size_t new_size)
{
  arr.resize(new_size);
}

class Remapeta : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Remapeta",
    .operators = { { "remapeta", 0, 0, "VCT filename", RemapetaHelp },
                   { "remapeta_s", 0, 0, "VCT filename", RemapetaHelp },
                   { "remapeta_z", 0, 0, "VCT filename", RemapetaHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Remapeta>(module);

  int REMAPETA_S, REMAPETA_Z;
  double cconst = 1.0E-6;
  size_t nfis2gp = 0;
  int nvars3D = 0;
  Varray<double> fis2;
  size_t numMissValsout = 0;
  bool lfis2 = false;
  int varids[MaxVars3D];
  Vmask imiss;
  double missval = 0.0;
  double cptop = 0.0;  // min. pressure level for cond.

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operatorID;
  int presID;

  VarIDs varIDs;

  Varray<double> ps2;

  Varray<float> f_t1, f_t2;
  Varray<float> f_q1, f_q2;

  Varray<double> d_t1, d_t2;
  Varray<double> d_q1, d_q2;

  Varray<double> tscor, pscor, secor;
  Varray<double> fis1;
  Varray<double> ps1;

  int zaxisID_ML = -1;

  int numFullLevels2;
  int numFullLevels1;

  size_t gridsize;

  bool ltq;

  VarList varList1;
  VarList varList2;

  MemType memType;

  Varray<double> sum1, sum2;
  Varray<double> deltap1, deltap2;
  Varray<double> halfPress1, halfPress2;

  Varray<double> vct1;
  Varray<double> vct2;

  const double *a1;
  const double *b1;

  const double *a2;
  const double *b2;

  template <typename T>
  void
  do_work(long nctop, Varray<T> &t2, Varray<T> &q2, Varray2D<T> &vars1, Varray2D<T> &vars2)
  {
    if (ltq)
    {
      int varID = varIDs.taID;
      for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
      {
        std::span sarray(&t2[gridsize * levelID], gridsize);

        auto mm = array_min_max_mask(sarray.data(), gridsize, imiss);
        if (mm.min < MIN_T || mm.max > MAX_T)
          cdo_warning("Output temperature at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);

        set_missval(imiss, missval, sarray);
        cdo_def_field(streamID2, varID, levelID);

        if (memType == MemType::Float)
          cdo_write_field_f(streamID2, (float *) sarray.data(), numMissValsout);
        else
          cdo_write_field(streamID2, (double *) sarray.data(), numMissValsout);
      }

      varID = varIDs.husID;
      for (int levelID = 0; levelID < varList2.vars[varID].nlevels; ++levelID)
      {
        std::span sarray(&q2[gridsize * levelID], gridsize);

        correct_humidity(sarray, MIN_Q);

        if (levelID < nctop)
          for (size_t i = 0; i < gridsize; ++i) sarray[i] = cconst;

        auto mm = array_min_max_mask(sarray.data(), gridsize, imiss);
        if (mm.min < MIN_Q || mm.max > MAX_Q)
          cdo_warning("Output humidity at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);

        set_missval(imiss, missval, sarray);
        cdo_def_field(streamID2, varID, levelID);

        if (memType == MemType::Float)
          cdo_write_field_f(streamID2, (float *) sarray.data(), numMissValsout);
        else
          cdo_write_field(streamID2, (double *) sarray.data(), numMissValsout);
      }
    }

    for (int iv = 0; iv < nvars3D; ++iv)
    {
      int varID = varids[iv];

      auto nlevels = varList2.vars[varID].nlevels;

      if (operatorID == REMAPETA_S)
      {
        vertSum(sum1, vars1[iv], gridsize, numFullLevels1);
        vertSum(sum2, vars2[iv], gridsize, numFullLevels2);
      }
      else if (operatorID == REMAPETA_Z)
      {
        vct_to_hybrid_pressure((double *) nullptr, halfPress1.data(), vct1, ps1.data(), numFullLevels1, gridsize);
        for (int k = 0; k < numFullLevels1; ++k)
          for (size_t i = 0; i < gridsize; ++i)
          {
            deltap1[k * gridsize + i] = halfPress1[(k + 1) * gridsize + i] - halfPress1[k * gridsize + i];
            deltap1[k * gridsize + i] = std::log(deltap1[k * gridsize + i]);
          }
        vertSumw(sum1, vars1[iv], gridsize, numFullLevels1, deltap1);

        vct_to_hybrid_pressure((double *) nullptr, halfPress2.data(), vct2, ps1.data(), numFullLevels2, gridsize);
        for (int k = 0; k < numFullLevels2; ++k)
          for (size_t i = 0; i < gridsize; ++i)
          {
            deltap2[k * gridsize + i] = halfPress2[(k + 1) * gridsize + i] - halfPress2[k * gridsize + i];
            deltap2[k * gridsize + i] = std::log(deltap2[k * gridsize + i]);
          }
        vertSumw(sum2, vars2[iv], gridsize, numFullLevels2, deltap2);
      }

      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        std::span sarray(&vars2[iv][gridsize * levelID], gridsize);

        if (operatorID == REMAPETA_S || operatorID == REMAPETA_Z)
          for (size_t i = 0; i < gridsize; ++i) sarray[i] *= sum1[i] / sum2[i];

        set_missval(imiss, missval, sarray);
        cdo_def_field(streamID2, varID, levelID);

        if (memType == MemType::Float)
          cdo_write_field_f(streamID2, (float *) sarray.data(), numMissValsout);
        else
          cdo_write_field(streamID2, (double *) sarray.data(), numMissValsout);
      }
    }
  }

  void
  read_fis(std::string const &fileName)
  {
    auto streamID = stream_open_read_locked(fileName.c_str());
    auto vlistID = streamInqVlist(streamID);
    VarList varList(vlistID);

    int varID, levelID;
    streamInqField(streamID, &varID, &levelID);

    auto const &var = varList.vars[0];

    nfis2gp = var.gridsize;
    fis2.resize(nfis2gp);

    size_t numMissVals;
    streamReadField(streamID, fis2.data(), &numMissVals);

    if (numMissVals)
    {
      imiss.resize(nfis2gp);
      for (size_t i = 0; i < nfis2gp; ++i) imiss[i] = fp_is_equal(fis2[i], var.missval);

      numMissValsout = numMissVals;
    }

    // check range of surface_geopotential
    auto mm = array_min_max_mask(fis2.data(), nfis2gp, imiss);
    if (mm.min < MIN_FIS || mm.max > MAX_FIS)
      cdo_warning("%s out of range (min=%g max=%g)!", var_stdname(surface_geopotential), mm.min, mm.max);

    if (mm.min < -1.e10 || mm.max > 1.e10) cdo_abort("%s out of range!", var_stdname(surface_geopotential));

    streamClose(streamID);
  }

public:
  void
  init() override
  {
    memType = Options::CDO_Memtype;

    REMAPETA_S = module.get_id("remapeta_s");
    REMAPETA_Z = module.get_id("remapeta_z");

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    const auto envstr = getenv("REMAPETA_PTOP");
    if (envstr)
    {
      auto fval = atof(envstr);
      if (fval > 0.0)
      {
        cptop = fval;
        cdo_print("Set REMAPETA_PTOP to %g", cptop);
      }
    }

    vct2 = vct_from_file(cdo_operator_argv(0));
    int nvct2 = vct2.size();
    numFullLevels2 = nvct2 / 2 - 1;

    a2 = vct2.data();
    b2 = vct2.data() + nvct2 / 2;
    if (Options::cdoVerbose)
      for (int i = 0; i < numFullLevels2 + 1; ++i) cdo_print("vct2: %5d %25.17f %25.17f", i, vct2[i], vct2[nvct2 / 2 + i]);

    streamID1 = cdo_open_read(0);

    if (cdo_operator_argc() == 2)
    {
      lfis2 = true;
      read_fis(cdo_operator_argv(1));
    }

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    vlist_unpack(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    auto gridID0 = vlistGrid(vlistID1, 0);
    if (gridInqType(gridID0) == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

    gridsize = vlist_check_gridsize(vlistID1);

    auto zaxisID2 = zaxisCreate(ZAXIS_HYBRID, numFullLevels2);

    {
      Varray<double> lev2(numFullLevels2);
      for (int i = 0; i < numFullLevels2; ++i) lev2[i] = i + 1;
      zaxisDefLevels(zaxisID2, lev2.data());
    }

    if (nvct2 == 0) cdo_abort("Internal problem, vct2 undefined!");
    zaxisDefVct(zaxisID2, nvct2, vct2.data());

    auto surfaceID = zaxis_from_name("surface");

    int numHybridLevels = 0, numHalfLevels1 = 0;
    vct1 = vlist_read_vct(vlistID1, zaxisID_ML, numHybridLevels, numFullLevels1, numHalfLevels1);
    int nvct1 = vct1.size();

    vlist_change_hybrid_zaxis(vlistID1, vlistID2, zaxisID_ML, zaxisID2);

    auto numZaxes = varList1.numZaxes();
    for (int i = 0; i < numZaxes; ++i)
    {
      auto zaxisID = vlistZaxis(vlistID1, i);
      auto nlevels = zaxisInqSize(zaxisID);
      if (zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevels == 1) vlistChangeZaxisIndex(vlistID2, i, surfaceID);
    }

    a1 = vct1.data();
    b1 = vct1.data() + nvct1 / 2;
    if (Options::cdoVerbose)
      for (int i = 0; i < nvct1 / 2; ++i) cdo_print("vct1: %5d %25.17f %25.17f", i, vct1[i], vct1[nvct1 / 2 + i]);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList2 = VarList(vlistID2);

    if (zaxisID_ML == -1) cdo_warning("No 3D variable with hybrid sigma pressure coordinate found!");

    auto numVars = varList1.numVars();

    varIDs = varList_search_varIDs(varList1, numFullLevels1);

    if (Options::cdoVerbose)
    {
      cdo_print("Found:");
      // clang-format off
        if (-1 != varIDs.taID)    cdo_print("  %s", var_stdname(air_temperature));
        if (-1 != varIDs.psID)      cdo_print("  %s", var_stdname(surface_air_pressure));
        if (-1 != varIDs.lnpsID)    cdo_print("  LOG(%s)", var_stdname(surface_air_pressure));
        if (-1 != varIDs.sgeopotID) cdo_print("  %s", var_stdname(surface_geopotential));
        if (-1 != varIDs.husID)     cdo_print("  %s", var_stdname(specific_humidity));
      // clang-format on
    }

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];

      if (var.gridType == GRID_SPECTRAL && var.zaxisType == ZAXIS_HYBRID) cdo_abort("Spectral data on model level unsupported!");

      if (var.gridType == GRID_SPECTRAL) cdo_abort("Spectral data unsupported!");

      if (var.zaxisType == ZAXIS_HYBRID && zaxisID_ML != -1 && var.nlevels == numFullLevels1)
      {
        if (!(varID == varIDs.taID || varID == varIDs.husID)) varids[nvars3D++] = varID;
      }
      else
      {
        if (varID == varIDs.taID) varIDs.taID = -1;
        if (varID == varIDs.husID) varIDs.husID = -1;
      }
    }

    ltq = (varIDs.taID != -1 && varIDs.husID != -1);

    if (!ltq)
    {
      if (varIDs.taID != -1) cdo_abort("Temperature without humidity unsupported!");
      if (varIDs.husID != -1) cdo_abort("Humidity without temperature unsupported!");
    }

    if (operatorID == REMAPETA_S || operatorID == REMAPETA_Z)
    {
      sum1.resize(gridsize);
      sum2.resize(gridsize);
    }

    if (operatorID == REMAPETA_Z)
    {
      deltap1.resize(gridsize * numFullLevels1);
      deltap2.resize(gridsize * numFullLevels2);
      halfPress1.resize(gridsize * (numFullLevels1 + 1));
      halfPress2.resize(gridsize * (numFullLevels2 + 1));
    }

    ps1 = Varray<double>(gridsize);
    fis1 = Varray<double>(gridsize);

    if (!lfis2) fis2.resize(gridsize);
    if (lfis2 && gridsize != nfis2gp) cdo_abort("Orographies have different grid size!");

    ps2 = Varray<double>(gridsize);

    if (ltq)
    {
      tscor.resize(gridsize);
      pscor.resize(gridsize);
      secor.resize(gridsize);

      auto gs1 = gridsize * numFullLevels1;
      auto gs2 = gridsize * numFullLevels2;
      (memType == MemType::Float) ? resize(f_t1, gs1) : resize(d_t1, gs1);
      (memType == MemType::Float) ? resize(f_q1, gs1) : resize(d_q1, gs1);

      (memType == MemType::Float) ? resize(f_t2, gs2) : resize(d_t2, gs2);
      (memType == MemType::Float) ? resize(f_q2, gs2) : resize(d_q2, gs2);
    }

    if (zaxisID_ML != -1 && varIDs.sgeopotID == -1)
    {
      std::ranges::fill(fis1, 0.0);
      if (ltq) cdo_warning("%s not found - set to zero!", var_stdname(surface_geopotential));
    }

    presID = varIDs.lnpsID;
    if (zaxisID_ML != -1 && varIDs.lnpsID == -1)
    {
      if (varIDs.psID == -1)
        cdo_abort("%s not found!", var_stdname(surface_air_pressure));
      else
        presID = varIDs.psID;
    }

    if (Options::cdoVerbose)
    {
      if (presID == varIDs.lnpsID)
        cdo_print("using LOG(%s)", var_stdname(surface_air_pressure));
      else
        cdo_print("using %s", var_stdname(surface_air_pressure));
    }

    if (Options::cdoVerbose) cdo_print("nvars3D = %d   ltq = %d", nvars3D, (int) ltq);
  }

  void
  run() override
  {
    Field field;
    Varray2D<float> f_vars1, f_vars2;
    Varray2D<double> d_vars1, d_vars2;
    auto numVars = varList1.numVars();
    if (nvars3D)
    {
      (memType == MemType::Float) ? resize(f_vars1, numVars) : resize(d_vars1, numVars);
      (memType == MemType::Float) ? resize(f_vars2, numVars) : resize(d_vars2, numVars);

      if (memType == MemType::Float)
      {
        for (int varID = 0; varID < nvars3D; ++varID) f_vars1[varID].resize(gridsize * numFullLevels1);
        for (int varID = 0; varID < nvars3D; ++varID) f_vars2[varID].resize(gridsize * numFullLevels2);
      }
      else
      {
        for (int varID = 0; varID < nvars3D; ++varID) d_vars1[varID].resize(gridsize * numFullLevels1);
        for (int varID = 0; varID < nvars3D; ++varID) d_vars2[varID].resize(gridsize * numFullLevels2);
      }
    }

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
        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        if (zaxisID_ML != -1)
        {
          auto offset = gridsize * levelID;

          if (varID == varIDs.sgeopotID)
            field_copy_array(gridsize, field, fis1.data());
          else if (varID == presID)
          {
            if (varIDs.lnpsID != -1)
            {
              if (field.memType == MemType::Float)
                for (size_t i = 0; i < gridsize; ++i) ps1[i] = std::exp((double) field.vec_f[i]);
              else
                for (size_t i = 0; i < gridsize; ++i) ps1[i] = std::exp(field.vec_d[i]);
            }
            else if (varIDs.psID != -1)
              field_copy_array(gridsize, field, ps1.data());
          }
          else if (ltq && varID == varIDs.taID)
          {
            if (memType == MemType::Float) { field_copy_array(gridsize, field, &f_t1[offset]); }
            else { field_copy_array(gridsize, field, &d_t1[offset]); }
          }
          else if (ltq && varID == varIDs.husID)
          {
            if (memType == MemType::Float) { field_copy_array(gridsize, field, &f_q1[offset]); }
            else { field_copy_array(gridsize, field, &d_q1[offset]); }
          }
          // else if ( var.zaxisID == zaxisID_ML )
          else if (var.zaxisType == ZAXIS_HYBRID && var.nlevels == numFullLevels1)
          {
            int i;
            for (i = 0; i < nvars3D; ++i)
              if (varID == varids[i]) break;

            if (i == nvars3D) cdo_abort("Internal error, 3D variable not found!");

            if (memType == MemType::Float) { field_copy_array(gridsize, field, &f_vars1[i][offset]); }
            else { field_copy_array(gridsize, field, &d_vars1[i][offset]); }
          }
          else
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field);
          }
        }
        else
        {
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, field);
        }
      }

      if (zaxisID_ML != -1)
      {
        // check range of psProg
        auto mm = array_min_max_mask(ps1.data(), gridsize, imiss);
        if (mm.min < MIN_PS || mm.max > MAX_PS) cdo_warning("Surface pressure out of range (min=%g max=%g)!", mm.min, mm.max);

        // check range of geop
        mm = array_min_max_mask(fis1.data(), gridsize, imiss);
        if (mm.min < MIN_FIS || mm.max > MAX_FIS) cdo_warning("Orography out of range (min=%g max=%g)!", mm.min, mm.max);
      }

      if (!lfis2)
        for (size_t i = 0; i < gridsize; ++i) fis2[i] = fis1[i];

      if (ltq)
      {
        int varID = varIDs.taID;
        auto const &var1 = varList1.vars[varID];
        for (int levelID = 0; levelID < var1.nlevels; ++levelID)
        {
          auto offset = gridsize * levelID;
          auto mm = [&]()
          {
            if (memType == MemType::Float) { return array_min_max_mask(&f_t1[offset], gridsize, imiss); }
            else { return array_min_max_mask(&d_t1[offset], gridsize, imiss); }
          }();

          if (mm.min < MIN_T || mm.max > MAX_T)
            cdo_warning("Input temperature at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);
        }

        varID = varIDs.husID;
        for (int levelID = 0; levelID < var1.nlevels; ++levelID)
        {
          auto offset = gridsize * levelID;

          if (memType == MemType::Float)
            correct_humidity(std::span(&f_q1[offset], gridsize), MIN_Q);
          else
            correct_humidity(std::span(&d_q1[offset], gridsize), MIN_Q);

          auto mm = [&]()
          {
            if (memType == MemType::Float) { return array_min_max_mask(&f_q1[offset], gridsize, imiss); }
            else { return array_min_max_mask(&d_q1[offset], gridsize, imiss); }
          }();
          if (mm.min < MIN_Q || mm.max > MAX_Q)
            cdo_warning("Input humidity at level %d out of range (min=%g max=%g)!", levelID + 1, mm.min, mm.max);
        }
      }

      if (nvars3D || ltq)
      {
        if (memType == MemType::Float)
          hetaeta(ltq, gridsize, imiss, numFullLevels1, a1, b1, fis1, ps1, f_t1, f_q1, numFullLevels2, a2, b2, fis2, ps2, f_t2,
                  f_q2, nvars3D, f_vars1, f_vars2, tscor, pscor, secor);
        else
          hetaeta(ltq, gridsize, imiss, numFullLevels1, a1, b1, fis1, ps1, d_t1, d_q1, numFullLevels2, a2, b2, fis2, ps2, d_t2,
                  d_q2, nvars3D, d_vars1, d_vars2, tscor, pscor, secor);
      }

      long nctop = (cptop > 0) ? ncctop(cptop, (long) numFullLevels2, (long) numFullLevels2 + 1, a2, b2) : 0;

      if (zaxisID_ML != -1 && varIDs.sgeopotID != -1)
      {
        int varID = varIDs.sgeopotID;
        int levelID = 0;
        set_missval(imiss, missval, fis2);
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, fis2.data(), numMissValsout);
      }

      if (zaxisID_ML != -1 && varIDs.lnpsID != -1)
        for (size_t i = 0; i < gridsize; ++i) ps2[i] = std::log(ps2[i]);

      if (zaxisID_ML != -1 && presID != -1)
      {
        int varID = presID;
        int levelID = 0;
        set_missval(imiss, missval, ps2);
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, ps2.data(), numMissValsout);
      }

      if (memType == MemType::Float) { do_work(nctop, f_t2, f_q2, f_vars1, f_vars2); }
      else { do_work(nctop, d_t2, d_q2, d_vars1, d_vars2); }

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
