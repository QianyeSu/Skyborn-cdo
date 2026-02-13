/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>
#include <cdi.h>

#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "util_string.h"

/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

static void
layer_cloud(Varray<double> const &cc, Varray<double> &ll, long maxLevIndex, long minLevIndex, long dimgp)
{
  constexpr double ZEPSEC = 1.0 - 1.0e-12;

  for (long i = 0; i < dimgp; ++i) ll[i] = 1.0 - cc[i + maxLevIndex * dimgp];

  for (long k = maxLevIndex + 1; k <= minLevIndex; ++k)
  {
    for (long i = 0; i < dimgp; ++i)
    {
      double maxval = std::max(cc[i + (k - 1) * dimgp], cc[i + k * dimgp]);
      double minval = std::min(cc[i + (k - 1) * dimgp], ZEPSEC);
      ll[i] *= (1.0 - maxval) / (1.0 - minval);
    }
  }

  for (long i = 0; i < dimgp; ++i) ll[i] = 1.0 - ll[i];
}

static void
vct2plev(Varray<double> const &vct, Varray<double> &plevs, long nlevels)
{
  constexpr double SCALESLP = 101325.0;
  for (long k = 0; k < nlevels; ++k) plevs[k] = vct[k] + vct[k + nlevels] * SCALESLP;
}

static void
hl_index(long &maxLevIndex, long &minLevIndex, double pmax, double pmin, long nlevels, Varray<double> const &levels)
{
  maxLevIndex = -1;
  minLevIndex = -1;

  for (long k = 0; k < nlevels; ++k)
    if (levels[k] > pmax)
    {
      maxLevIndex = k - 1;
      break;
    }

  for (long k = nlevels - 1; k >= 0; --k)
    if (levels[k] < pmin)
    {
      minLevIndex = k;
      break;
    }
}

static void
pl_index(long &maxLevIndex, long &minLevIndex, double pmax, double pmin, long nlevels, Varray<double> const &levels)
{
  maxLevIndex = -1;
  minLevIndex = -1;

  for (long k = 0; k < nlevels; ++k)
    if (levels[k] >= pmax)
    {
      maxLevIndex = k;
      break;
    }

  for (long k = nlevels - 1; k >= 0; --k)
    if (levels[k] < pmin)
    {
      minLevIndex = k;
      break;
    }
}

class Cloudlayer : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Cloudlayer",
    .operators = { { "cloudlayer" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Cloudlayer> registration = RegisterEntry<Cloudlayer>();

private:
  static const int MaxCloudLayers = 3;
  int gridID{}, zaxisID{};
  bool zrev = false;
  int aclcacID = -1;
  int numVars2 = 0;
  int aclcac_code_found = 0;
  long kmin[MaxCloudLayers] = { -1, -1, -1 }, kmax[MaxCloudLayers] = { -1, -1, -1 };
  double sfclevel = 0;
  double pmin = 0, pmax = 0;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int numLevels{};
  size_t gridsize{};
  double missval{};

  void
  define_cld_lay(int surfaceID)
  {
    auto varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
    vlistDefVarParam(vlistID2, varID, cdiEncodeParam(33, 128, 255));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "cld_lay");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "cloud layer");
    vlistDefVarMissval(vlistID2, varID, missval);
  }

  void
  define_cld_vars(int surfaceID)
  {
    auto varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
    vlistDefVarParam(vlistID2, varID, cdiEncodeParam(34, 128, 255));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "low_cld");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "low cloud");
    vlistDefVarMissval(vlistID2, varID, missval);

    varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
    vlistDefVarParam(vlistID2, varID, cdiEncodeParam(35, 128, 255));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "mid_cld");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "mid cloud");
    vlistDefVarMissval(vlistID2, varID, missval);

    varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
    vlistDefVarParam(vlistID2, varID, cdiEncodeParam(36, 128, 255));
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "hih_cld");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "high cloud");
    vlistDefVarMissval(vlistID2, varID, missval);
  }

public:
  void
  init() override
  {
    if (cdo_operator_argc() > 0)
    {
      operator_check_argc(2);
      numVars2 = 1;
      pmin = parameter_to_double(cdo_operator_argv(0));
      pmax = parameter_to_double(cdo_operator_argv(1));
    }
    else { numVars2 = MaxCloudLayers; }

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1(vlistID1);

    gridsize = vlist_check_gridsize(vlistID1);

    auto aclcac_code = 223;

    for (auto const &var1 : varList1.vars)
    {
      zaxisID = var1.zaxisID;
      auto code = var1.code;

      if (code <= 0)
      {
        if (string_to_lower(var1.name) == "aclcac") code = 223;
      }

      if (code == aclcac_code)
      {
        aclcac_code_found = 1;
        if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE || zaxisInqType(zaxisID) == ZAXIS_HYBRID)
        {
          aclcacID = var1.ID;
          break;
        }
      }
    }

    if (aclcacID == -1)
    {
      cdo_abort("Cloud cover (parameter 223) not found%s!", aclcac_code_found ? " on pressure or hybrid levels" : "");
    }

    auto const &aclcacVar = varList1.vars[aclcacID];
    missval = aclcacVar.missval;
    gridID = aclcacVar.gridID;
    zaxisID = aclcacVar.zaxisID;
    numLevels = aclcacVar.nlevels;
    auto nhlev = numLevels + 1;

    if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
    {
      Varray<double> plevs(numLevels);
      zaxisInqLevels(zaxisID, plevs.data());
      if (plevs[0] > plevs[numLevels - 1])
      {
        zrev = true;
        for (int levelID = 0; levelID < numLevels / 2; ++levelID) std::swap(plevs[levelID], plevs[numLevels - 1 - levelID]);
      }
      /*
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          printf("level %d %g\n", levelID, plevs[levelID]);
        }
      */
      if (numVars2 == 1) { pl_index(kmax[0], kmin[0], pmin, pmax, numLevels, plevs); }
      else
      {
        pl_index(kmax[2], kmin[2], 5000., 44000., numLevels, plevs);
        pl_index(kmax[1], kmin[1], 46000., 73000., numLevels, plevs);
        pl_index(kmax[0], kmin[0], 75000., 101300., numLevels, plevs);
      }
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID)
    {
      int nvct = zaxisInqVctSize(zaxisID);
      if (numLevels == (nvct / 2 - 1))
      {
        Varray<double> vct(nvct);
        zaxisInqVct(zaxisID, vct.data());

        auto nlevs = numLevels + 1;
        Varray<double> plevs(nlevs);
        vct2plev(vct, plevs, nlevs);

        if (numVars2 == 1) { hl_index(kmax[0], kmin[0], pmin, pmax, nhlev, plevs); }
        else
        {
          hl_index(kmax[2], kmin[2], 5000., 44000., nhlev, plevs);
          hl_index(kmax[1], kmin[1], 46000., 73000., nhlev, plevs);
          hl_index(kmax[0], kmin[0], 75000., 101300., nhlev, plevs);
        }
      }
      else
        cdo_abort("Unsupported vertical coordinate table format!");
    }
    else
      cdo_abort("Unsupported Z-Axis type!");

    auto surfaceID = zaxisCreate(ZAXIS_SURFACE, 1);
    zaxisDefLevels(surfaceID, &sfclevel);

    vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

    (numVars2 == 1) ? define_cld_lay(surfaceID) : define_cld_vars(surfaceID);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Varray<double> aclcac(gridsize * numLevels);
    Varray<double> cloud[MaxCloudLayers];
    for (int varID = 0; varID < numVars2; ++varID) cloud[varID].resize(gridsize);

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

        size_t offset = zrev ? (numLevels - 1 - levelID) * gridsize : levelID * gridsize;

        if (varID == aclcacID)
        {
          size_t numMissVals;
          cdo_read_field(streamID1, aclcac.data() + offset, &numMissVals);
          if (numMissVals != 0) cdo_abort("Missing values unsupported!");
        }
      }

      for (int varID = 0; varID < numVars2; ++varID)
      {
        for (size_t i = 0; i < gridsize; ++i) cloud[varID][i] = missval;
      }

      for (int varID = 0; varID < numVars2; ++varID)
      {
        if (kmax[varID] != -1 && kmin[varID] != -1) layer_cloud(aclcac, cloud[varID], kmax[varID], kmin[varID], gridsize);
      }

      for (int varID = 0; varID < numVars2; ++varID)
      {
        auto numMissVals = varray_num_mv(gridsize, cloud[varID], missval);

        cdo_def_field(streamID2, varID, 0);
        cdo_write_field(streamID2, cloud[varID].data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
