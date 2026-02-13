/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>
#include <climits>

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "libncl.h"
#include "pmlist.h"

static void
uv2dv_cfd_W(double missval, double *u, double *v, double *lon, double *lat, size_t nlon, size_t nlat, size_t nlev, int boundOpt,
            double *div)
{
  int ierror;

  // Test dimension sizes.
  if ((nlon > INT_MAX) || (nlat > INT_MAX)) cdo_abort("nlat and/or nlon is greater than INT_MAX!");

  int inlon = (int) nlon;
  int inlat = (int) nlat;

  size_t gridsize_uv = nlat * nlon;

  for (size_t k = 0; k < nlev; ++k)
  {
    double *tmp_u = u + k * gridsize_uv;
    double *tmp_v = v + k * gridsize_uv;
    double *tmp_div = div + k * gridsize_uv;
    // Init output array.
    std::ranges::fill_n(tmp_div, gridsize_uv, 0.0);
    // Call the Fortran routine.
#ifdef HAVE_CF_INTERFACE
    DDVFIDF(tmp_u, tmp_v, lat, lon, inlon, inlat, missval, boundOpt, tmp_div, ierror);
#else
    cdo_abort("Fortran support not compiled in!");
#endif
  }
}

static void
uv2vr_cfd_W(double missval, double *u, double *v, double *lon, double *lat, size_t nlon, size_t nlat, size_t nlev, int boundOpt,
            double *vort)
{
  int ierror;

  // Test dimension sizes.
  if ((nlon > INT_MAX) || (nlat > INT_MAX)) cdo_abort("nlat and/or nlon is greater than INT_MAX!");

  int inlon = (int) nlon;
  int inlat = (int) nlat;

  size_t gridsize_uv = nlat * nlon;

  for (size_t k = 0; k < nlev; ++k)
  {
    double *tmp_u = u + k * gridsize_uv;
    double *tmp_v = v + k * gridsize_uv;
    double *tmp_vort = vort + k * gridsize_uv;
    // Init output array.
    std::ranges::fill_n(tmp_vort, gridsize_uv, 0.0);
    // Call the Fortran routine.
#ifdef HAVE_CF_INTERFACE
    DVRFIDF(tmp_u, tmp_v, lat, lon, inlon, inlat, missval, boundOpt, tmp_vort, ierror);
#else
    cdo_abort("Fortran support not compiled in!");
#endif
  }
}

static int
find_name(VarList const &varList, std::string const &name)
{
  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    if (name == varList.vars[varID].name) return varID;
  }

  return CDI_UNDEFID;
}

enum struct OutMode
{
  NEW,
  APPEND,
  REPLACE
};

// Parameter
static OutMode outMode(OutMode::NEW);
static int boundOpt = -1;
static std::string uName;
static std::string vName;

static void
print_parameter(void)
{
  cdo_print("u=%s, v=%s, boundOpt=%d, outMode=%s", uName, vName, boundOpt,
            (outMode == OutMode::NEW)      ? "new"
            : (outMode == OutMode::APPEND) ? "append"
                                           : "replace");
}

static void
set_parameter(void)
{
  uName = "u";
  vName = "v";

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = "PARAMETER";
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "u") uName = value;
      else if (key == "v") vName = value;
      else if (key == "boundOpt") boundOpt = parameter_to_int(value);
      else if (key == "outMode")
      {
        if      (value == "new")     outMode = OutMode::NEW;
        else if (value == "append")  outMode = OutMode::APPEND;
        else if (value == "replace") outMode = OutMode::REPLACE;
        else cdo_abort("Invalid parameter key value: outMode=%s (valid are: new/append/replace)", value);
      }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  if (Options::cdoVerbose) print_parameter();
}

class NCL_wind : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "NCL_wind",
    // clang-format off
    .operators = { { "uv2dv_cfd", 0, 0, "[u,v,boundsOpt,outMode]", Ncl_windHelp },
                   { "uv2vr_cfd", 0, 0, "[u,v,boundsOpt,outMode]", Ncl_windHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<NCL_wind> registration = RegisterEntry<NCL_wind>();

  int UV2DV_CFD{}, UV2VR_CFD{};

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int operatorID{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int varIDu{};
  int varIDv{};
  int varIDo{};
  int nlev{};

  size_t gridsizeuv{};
  size_t nlon{};
  size_t nlat{};

  Varray<double> lon;
  Varray<double> lat;

  VarList varList1;

public:
  void
  init() override
  {
    UV2DV_CFD = module.get_id("uv2dv_cfd");
    UV2VR_CFD = module.get_id("uv2vr_cfd");

    operatorID = cdo_operator_id();

    set_parameter();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    vlistID2 = CDI_UNDEFID;
    if (outMode == OutMode::NEW)
      vlistID2 = vlistCreate();
    else if (outMode == OutMode::APPEND)
      vlistID2 = vlistDuplicate(vlistID1);
    else
      cdo_abort("outMode=%d unsupported!", outMode);

    varList1 = VarList(vlistID1);

    varIDu = find_name(varList1, uName);
    varIDv = find_name(varList1, vName);

    if (varIDu == CDI_UNDEFID) cdo_abort("%s not found!", uName);
    if (varIDv == CDI_UNDEFID) cdo_abort("%s not found!", vName);

    auto const &varU = varList1.vars[varIDu];
    auto const &varV = varList1.vars[varIDv];

    auto gridIDu = varU.gridID;
    auto gridIDv = varV.gridID;
    auto gridtype = gridInqType(gridIDu);
    gridsizeuv = gridInqSize(gridIDu);

    if (!((gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) && gridtype == gridInqType(gridIDv)))
      cdo_abort("u and v must be on a regular lonlat or Gaussian grid!");

    if (gridsizeuv != gridInqSize(gridIDv)) cdo_abort("u and v must have the same grid size!");

    if (boundOpt == -1) boundOpt = gridIsCircular(gridIDu) ? 1 : 0;
    if (Options::cdoVerbose) print_parameter();
    if (boundOpt < 0 || boundOpt > 3) cdo_abort("Parameter boundOpt=%d out of bounds (0-3)!", boundOpt);

    nlon = gridInqXsize(gridIDu);
    nlat = gridInqYsize(gridIDu);

    auto zaxisIDu = varU.zaxisID;
    nlev = zaxisInqSize(zaxisIDu);

    if (nlev != zaxisInqSize(varU.zaxisID)) cdo_abort("u and v must have the same number of level!");

    varIDo = vlistDefVar(vlistID2, gridIDu, zaxisIDu, varU.timeType);
    if (operatorID == UV2DV_CFD)
    {
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_NAME, "d");
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_LONGNAME, "divergence");
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_UNITS, "1/s");
    }
    else if (operatorID == UV2VR_CFD)
    {
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_NAME, "vo");
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_LONGNAME, "vorticity");
      cdiDefKeyString(vlistID2, varIDo, CDI_KEY_UNITS, "1/s");
    }

    vlistDefVarMissval(vlistID2, varIDo, varU.missval);

    lon.resize(nlon);
    lat.resize(nlat);
    gridInqXvals(gridIDu, &lon[0]);
    gridInqYvals(gridIDu, &lat[0]);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Varray<double> array(varList1.gridsizeMax());
    Varray<double> arrayu(nlev * gridsizeuv);
    Varray<double> arrayv(nlev * gridsizeuv);
    Varray<double> arrayo(nlev * gridsizeuv);

    auto const &varU = varList1.vars[varIDu];
    auto const &varV = varList1.vars[varIDv];
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      size_t numMissValsu = 0, numMissValsv = 0;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, &array[0], &numMissVals);

        if (varID == varIDu || varID == varIDv)
        {
          if (varID == varIDu)
          {
            std::copy_n(&array[0], gridsizeuv, &arrayu[levelID * gridsizeuv]);
            numMissValsu += numMissVals;
          }
          if (varID == varIDv)
          {
            std::copy_n(&array[0], gridsizeuv, &arrayv[levelID * gridsizeuv]);
            numMissValsv += numMissVals;
          }
        }

        if (outMode == OutMode::APPEND)
        {
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, &array[0], numMissVals);
        }
      }

      if (numMissValsu != numMissValsv)
      {
        cdo_abort("u and v have different number of missing values!");
        if (numMissValsu && fp_is_not_equal(varU.missval, varV.missval))
        {
          for (int levelID = 0; levelID < nlev; ++levelID)
          {
            auto parray = &arrayv[levelID * gridsizeuv];
            for (size_t i = 0; i < gridsizeuv; ++i)
              if (fp_is_equal(parray[i], varV.missval)) parray[i] = varU.missval;
          }
        }
      }

      if (operatorID == UV2DV_CFD)
        uv2dv_cfd_W(varU.missval, &arrayu[0], &arrayv[0], &lon[0], &lat[0], nlon, nlat, nlev, boundOpt, &arrayo[0]);
      else if (operatorID == UV2VR_CFD)
        uv2vr_cfd_W(varU.missval, &arrayu[0], &arrayv[0], &lon[0], &lat[0], nlon, nlat, nlev, boundOpt, &arrayo[0]);

      for (int levelID = 0; levelID < nlev; ++levelID)
      {
        auto parray = &arrayo[levelID * gridsizeuv];
        auto numMissVals = array_num_mv(gridsizeuv, parray, varU.missval);
        cdo_def_field(streamID2, varIDo, levelID);
        cdo_write_field(streamID2, parray, numMissVals);
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
