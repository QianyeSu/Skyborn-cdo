/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Adisit      adisit          compute insitu from potential temperature
      Adisit      adipot          compute potential from insitu temperature
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_zaxis.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include "util_string.h"

/*
  Transformation from potential to in situ temperature according to Bryden, 1973,
  "New polynomials for thermal expansion, adiabatic temperature gradient and potential temperature of sea water".
  Deep Sea Research and Oceanographic Abstracts. 20, 401-408 (GILL P.602), which gives the inverse
  transformation for an approximate value, all terms linear in t are taken after that one newton step.
  For the check value 8.4678516 the accuracy is 0.2 mikrokelvin.
*/

// compute insitu temperature from potential temperature
static inline double
adisit_kernel(double tpot, double sal, double p)
{
  constexpr double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9, a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
                   a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10, a_d = 4.1057E-9, a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double qc = p * (a_a1 + p * (a_c1 - a_e1 * p));
  double qv = p * (a_b1 - a_d * p);
  double dc = 1.0 + p * (-a_a2 + p * (a_c2 - a_e2 * p));
  double dv = a_b2 * p;
  double qnq = -p * (-a_a3 + p * a_c3);
  double qn3 = -p * a_a4;

  double tpo = tpot;
  double qvs = qv * (sal - 35.0) + qc;
  double dvs = dv * (sal - 35.0) + dc;
  double t = (tpo + qvs) / dvs;
  double fne = -qvs + t * (dvs + t * (qnq + t * qn3)) - tpo;
  double fst = dvs + t * (2.0 * qnq + 3.0 * qn3 * t);
  t = t - fne / fst;

  return t;
}

// compute potential temperature from insitu temperature
// Ref: Gill, p. 602, Section A3.5:Potential Temperature
static inline double
adipot_kernel(double t, double s, double p)
{
  constexpr double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9, a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
                   a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10, a_d = 4.1057E-9, a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double s_rel = s - 35.0;

  double aa = (a_a1 + t * (a_a2 - t * (a_a3 - a_a4 * t)));
  double bb = s_rel * (a_b1 - a_b2 * t);
  double cc = (a_c1 + t * (-a_c2 + a_c3 * t));
  double cc1 = a_d * s_rel;
  double dd = (-a_e1 + a_e2 * t);

  double tpot = t - p * (aa + bb + p * (cc - cc1 + p * dd));

  return tpot;
}

static void
calc_adisit(size_t gridsize, size_t nlevel, Varray<double> const &pressure, FieldVector const &tho, FieldVector const &sao,
            FieldVector &tis)
{
  // pressure units: hPa
  // tho units:      Celsius
  // sao units:      psu

  for (size_t levelID = 0; levelID < nlevel; ++levelID)
  {
    auto const &thovec = tho[levelID].vec_d;
    auto const &saovec = sao[levelID].vec_d;
    auto &tisvec = tis[levelID].vec_d;
    auto thoMissval = tho[levelID].missval;
    auto saoMissval = sao[levelID].missval;
    auto tisMissval = tis[levelID].missval;
    for (size_t i = 0; i < gridsize; ++i)
    {
      auto isMissing = (fp_is_equal(thovec[i], thoMissval) || fp_is_equal(saovec[i], saoMissval));
      tisvec[i] = isMissing ? tisMissval : adisit_kernel(thovec[i], saovec[i], pressure[levelID]);
    }
  }
}

static void
calc_adipot(size_t gridsize, size_t nlevel, Varray<double> const &pressure, FieldVector const &t, FieldVector const &s,
            FieldVector &tpot)
{
  // pressure units: hPa
  // t units:        Celsius
  // s units:        psu

  for (size_t levelID = 0; levelID < nlevel; ++levelID)
  {
    auto const &tvec = t[levelID].vec_d;
    auto const &svec = s[levelID].vec_d;
    auto &tpotvec = tpot[levelID].vec_d;
    auto tMissval = t[levelID].missval;
    auto sMissval = s[levelID].missval;
    auto tpotMissval = tpot[levelID].missval;
    for (size_t i = 0; i < gridsize; ++i)
    {
      auto isMissing = (fp_is_equal(tvec[i], tMissval) || fp_is_equal(svec[i], sMissval));
      tpotvec[i] = isMissing ? tpotMissval : adipot_kernel(tvec[i], svec[i], pressure[levelID]);
    }
  }
}

int
get_code(const CdoVar &var, std::string const &cname)
{
  auto code = var.code;
  if (code <= 0)
  {
    auto varname = string_to_lower(var.name);
    auto stdname = string_to_lower(var.stdname);

    if (varname == "s" || varname == "so" || stdname == "sea_water_salinity") { code = 5; }
    else if (varname == "t" || varname == "to") { code = 2; }

    if (stdname == cname) code = 2;
  }

  return code;
}

namespace
{
struct IOSettings
{
  CdoStreamID streamID2;
  int vlistID2{ CDI_UNDEFID };
  size_t gridsize;
  int numLevels;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int tisID2;
  int saoID2;
};
}  // namespace

IOSettings
configureOutput(std::function<void(int, int)> const &outputSettingFunc, VarList const &varList, int vlistID, int thoID, int saoID,
                Varray<double> &pressure)
{
  double pin = (cdo_operator_argc() == 1) ? parameter_to_double(cdo_operator_argv(0)) : -1.0;

  auto const &vars = varList.vars;

  auto units = vars[thoID].units;
  if (units.empty()) units = "Celcius";

  auto gridID = vlistGrid(vlistID, 0);
  auto gridsize = vlist_check_gridsize(vlistID);

  auto nlevels1 = vars[saoID].nlevels;
  auto nlevels2 = vars[thoID].nlevels;
  auto zaxisID = vars[thoID].zaxisID;

  if (nlevels1 != nlevels2) cdo_abort("temperature and salinity have different number of levels!");
  auto numLevels = nlevels1;

  pressure.resize(numLevels);
  cdo_zaxis_inq_levels(zaxisID, pressure.data());

  if (pin >= 0)
    for (int i = 0; i < numLevels; ++i) pressure[i] = pin;
  else
    for (int i = 0; i < numLevels; ++i) pressure[i] /= 10;

  if (Options::cdoVerbose)
  {
    cdo_print("Level Pressure");
    for (int i = 0; i < numLevels; ++i) cdo_print("%5d  %g", i + 1, pressure[i]);
  }

  int datatype = CDI_DATATYPE_FLT32;
  if (vars[thoID].dataType == CDI_DATATYPE_FLT64 && vars[saoID].dataType == CDI_DATATYPE_FLT64) datatype = CDI_DATATYPE_FLT64;

  auto vlistID2 = vlistCreate();
  vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID));

  auto tisID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);

  outputSettingFunc(vlistID2, tisID2);
  cdiDefKeyString(vlistID2, tisID2, CDI_KEY_UNITS, units.c_str());
  vlistDefVarMissval(vlistID2, tisID2, vars[thoID].missval);
  vlistDefVarDatatype(vlistID2, tisID2, datatype);

  auto saoID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);
  vlistDefVarParam(vlistID2, saoID2, cdiEncodeParam(5, 255, 255));
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_NAME, "s");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_LONGNAME, "Sea water salinity");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_STDNAME, "sea_water_salinity");
  cdiDefKeyString(vlistID2, saoID2, CDI_KEY_UNITS, "psu");
  vlistDefVarMissval(vlistID2, saoID2, vars[saoID].missval);
  vlistDefVarDatatype(vlistID2, saoID2, datatype);

  auto taxisID1 = vlistInqTaxis(vlistID);
  auto taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  auto streamID2 = cdo_open_write(1);
  cdo_def_vlist(streamID2, vlistID2);

  return IOSettings{ streamID2, vlistID2, gridsize, numLevels, taxisID1, taxisID2, tisID2, saoID2 };
}

const auto outputSetting_ADISIT = [](int p_vlistID2, int p_tisID2)
{
  vlistDefVarParam(p_vlistID2, p_tisID2, cdiEncodeParam(20, 255, 255));
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_NAME, "to");
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_LONGNAME, "Sea water temperature");
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_STDNAME, "sea_water_temperature");
};

const auto outputSetting_ADIPOT = [](int p_vlistID2, int p_tisID2)
{
  vlistDefVarParam(p_vlistID2, p_tisID2, cdiEncodeParam(2, 255, 255));
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_NAME, "tho");
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_LONGNAME, "Sea water potential temperature");
  cdiDefKeyString(p_vlistID2, p_tisID2, CDI_KEY_STDNAME, "sea_water_potential_temperature");
};

static void
check_tho_range(FieldVector &tho, int tsID, int levelID)
{
  constexpr double MIN_T = -10.0;
  constexpr double MAX_T = 40.0;
  auto mm = field_min_max(tho[levelID]);
  if (mm.min < MIN_T || mm.max > MAX_T)
    cdo_warning("Temperature in degree Celsius out of range (min=%g max=%g) [timestep:%d levelIndex:%d]!", mm.min, mm.max, tsID + 1,
                levelID + 1);
}

class Adisit : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Adisit",
    .operators = { { "adisit", 0, 0, "", AdisitHelp }, { "adipot", 0, 0, "", AdisitHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Adisit> registration = RegisterEntry<Adisit>();

private:
  int ADISIT{}, ADIPOT{};
  int thoID = -1;
  int saoID = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int operatorID{};

  size_t gridsize{};
  int numLevels{};

  int tisID2{};
  int saoID2{};

  Varray<double> pressure;

  VarList varList1;

public:
  void
  init() override
  {
    ADISIT = module.get_id("adisit");
    ADIPOT = module.get_id("adipot");

    operatorID = cdo_operator_id();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    varList1 = VarList(vlistID1);

    std::string cname_ADISIT{ "sea_water_potential_temperature" };
    std::string cname_ADIPOT{ "sea_water_temperature" };
    const std::string &cname = (operatorID == ADISIT) ? cname_ADISIT : cname_ADIPOT;

    for (auto const &var : varList1.vars)
    {
      auto code = get_code(var, cname);
      if (code == 2) { thoID = var.ID; }
      else if (code == 20 && operatorID == ADIPOT) { thoID = var.ID; }
      else if (code == 5) { saoID = var.ID; }
    }

    if (saoID == -1) cdo_abort("Sea water salinity not found!");
    if (thoID == -1) cdo_abort("%s temperature not found!", (operatorID == ADISIT) ? "Potential" : "Insitu");

    auto const &outputSettingFunc = (operatorID == ADISIT) ? outputSetting_ADISIT : outputSetting_ADIPOT;
    auto const &configResults = configureOutput(outputSettingFunc, varList1, vlistID1, thoID, saoID, pressure);

    streamID2 = configResults.streamID2;
    vlistID2 = configResults.vlistID2;
    gridsize = configResults.gridsize;
    numLevels = configResults.numLevels;
    taxisID1 = configResults.taxisID1;
    taxisID2 = configResults.taxisID2;
    tisID2 = configResults.tisID2;
    saoID2 = configResults.saoID2;
  }

  void
  run() override
  {
    FieldVector tho(numLevels), sao(numLevels), tis(numLevels);
    for (int levelID = 0; levelID < numLevels; ++levelID)
    {
      tho[levelID].resize(gridsize);
      sao[levelID].resize(gridsize);
      tis[levelID].resize(gridsize);
      tho[levelID].missval = varList1.vars[thoID].missval;
      sao[levelID].missval = varList1.vars[saoID].missval;
      tis[levelID].missval = tho[levelID].missval;
    }

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (varID == thoID) cdo_read_field(streamID1, tho[levelID]);
        if (varID == saoID) cdo_read_field(streamID1, sao[levelID]);
        if (varID == thoID) check_tho_range(tho, tsID, levelID);
      }

      auto adi_func = (operatorID == ADISIT) ? calc_adisit : calc_adipot;
      adi_func(gridsize, numLevels, pressure, tho, sao, tis);

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        field_num_mv(tis[levelID]);
        cdo_def_field(streamID2, tisID2, levelID);
        cdo_write_field(streamID2, tis[levelID]);
      }

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        field_num_mv(sao[levelID]);
        cdo_def_field(streamID2, saoID2, levelID);
        cdo_write_field(streamID2, sao[levelID]);
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
