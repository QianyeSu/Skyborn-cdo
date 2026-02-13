/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Rhopot      rhopot          potential density
*/

#include <algorithm>
#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "cdo_vlist.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "util_string.h"
#include "cdo_zaxis.h"

/*
!>
!! transformation from potential to in situ temperature
!! according to Bryden, 1973, "New polynomials for thermal expansion,
!! adiabatic temperature gradient and potential temperature of sea
!! water". Deep Sea Research and Oceanographic Abstracts. 20, 401-408
!! (GILL P.602), which gives the inverse transformation for an
!! approximate value, all terms linear in t are taken after that one
!! newton step.  for the check value 8.4678516 the accuracy is 0.2
!! mikrokelvin.
!!
*/

// compute density from insitu temperature
static double
potrho_1(const double t, const double sal, const double p)
{
  // clang-format off
  constexpr double r_a0 = 999.842594, r_a1 = 6.793952e-2, r_a2 = -9.095290e-3,
         r_a3 = 1.001685e-4, r_a4 = -1.120083e-6, r_a5 = 6.536332e-9,
         r_b0 = 8.24493e-1, r_b1 = -4.0899e-3, r_b2 = 7.6438e-5,
         r_b3 = -8.2467e-7, r_b4 = 5.3875e-9,
         r_c0 = -5.72466e-3, r_c1 = 1.0227e-4, r_c2 = -1.6546e-6,
         r_d0 = 4.8314e-4,
         r_e0 = 19652.21, r_e1 = 148.4206, r_e2 = -2.327105,
         r_e3 = 1.360477e-2, r_e4 = -5.155288e-5,
         r_f0 = 54.6746, r_f1 = -0.603459, r_f2 = 1.09987e-2,
         r_f3 = -6.1670e-5,
         r_g0 = 7.944e-2, r_g1 = 1.6483e-2, r_g2 = -5.3009e-4,
         r_h0 = 3.239908, r_h1 = 1.43713e-3, r_h2 = 1.16092e-4,
         r_h3 = -5.77905e-7,
         r_ai0 = 2.2838e-3, r_ai1 = -1.0981e-5, r_ai2 = -1.6078e-6,
         r_aj0 = 1.91075e-4,
         r_ak0 = 8.50935e-5, r_ak1 = -6.12293e-6, r_ak2 = 5.2787e-8,
         r_am0 = -9.9348e-7, r_am1 = 2.0816e-8, r_am2 = 9.1697e-10;

  double s = std::max(sal, 0.0);
  double s3h = std::sqrt(s*s*s);

  double rho = (r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))
            + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))
            + r_d0 * s*s
            + s3h * (r_c0 + t * (r_c1 + r_c2 * t)))
           / (1.
             - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))
                         + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))
                         + r_aj0 * s3h
                         + (r_ak0 + t * (r_ak1 + t * r_ak2)
                         + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)
                    + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))
                    + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))
                    + s3h * (r_g0 + t * (r_g1 + r_g2 * t))));
  // clang-format on
  return rho;
}

/*
#define N 4
int main (int argc, char *argv[])
{
  int i;
  {
    double p    = 0;
    double t[N] = {22, 25, 28, 31};
    double s[N] = {35, 35, 35, 35};
    double x[N] = {24.219, 23.343, 22.397, 21.384};
    double r[N];

    potrho_1d(t, s, p, r, N);

    for ( i = 0; i < N; ++i )
      printf("%d %5g %3g %8g %8g %8g %10.3f\n", i, p, s[i], t[i], x[i], r[i],
r[i]-x[i]);
  }

  {
    double p    = 300;
    double t[N] = {-2.140, -0.186, 1.771, 3.728};
    double s[N] = {35, 35, 35, 35};
    double x[N] = {42.191, 41.941, 41.649, 41.319};
    double r[N];

    potrho_1d(t, s, p, r, N);

    for ( i = 0; i < N; ++i )
      printf("%d %5g %3g %8g %8g %8g %10.3f\n", i, p, s[i], t[i], x[i], r[i],
r[i]-x[i]);
  }

  return 0;
}
*/

static void
calc_rhopot(size_t gridsize, size_t nlevel, Varray<double> const &pressure, FieldVector const &to, FieldVector const &sao,
            FieldVector &rho)
{
  // pressure units: hPa
  // to units:       Celsius
  // sao units:      psu

  for (size_t levelID = 0; levelID < nlevel; ++levelID)
  {
    auto const &tovec = to[levelID].vec_d;
    auto const &saovec = sao[levelID].vec_d;
    auto &rhovec = rho[levelID].vec_d;
    auto to_missval = to[levelID].missval;
    auto sao_missval = sao[levelID].missval;
    auto rho_missval = rho[levelID].missval;
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_equal(tovec[i], to_missval) || fp_is_equal(saovec[i], sao_missval))
        rhovec[i] = rho_missval;
      else
        rhovec[i] = potrho_1(tovec[i], saovec[i], pressure[levelID]);
    }
  }
}

class Rhopot : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Rhopot",
    .operators = { { "rhopot", RhopotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Rhopot> registration = RegisterEntry<Rhopot>();

  int zaxisID{};
  int toID = -1, saoID = -1, thoID = -1;
  double pin = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int numLevels{};
  int gridsize{};

  FieldVector to;
  FieldVector sao;
  FieldVector rho;

  Varray<double> pressure;

public:
  void
  init() override
  {
    if (cdo_operator_argc() == 1) pin = parameter_to_double(cdo_operator_argv(0));

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1(vlistID1);
    auto numVars = varList1.numVars();

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto code = var.code;
      if (code <= 0)
      {
        auto stdname = string_to_lower(var.stdname);

        // clang-format off
        if      (var.name == "to") code = 20;
        else if (var.name == "sao") code = 5;
        else if (var.name == "tho") code = 2;
        else if (var.name == "s") code = 5;
        else if (var.name == "t") code = 2;
        else if (stdname == "sea_water_salinity") code = 5;
        else if (stdname == "sea_water_potential_temperature") code = 2;
        // clang-format on
      }

      // clang-format off
      if      (code == 20) toID = varID;
      else if (code == 5) saoID = varID;
      else if (code == 2) thoID = varID;
      // clang-format on
    }

    if (saoID == -1) cdo_abort("Sea water salinity not found!");
    if (toID == -1 && thoID != -1)
    {
      cdo_print("Use the CDO operator 'adisit' to convert potential temperature to In-situ temperature.");
      cdo_print("Here is an example:");
      cdo_print("   cdo rhopot -adisit %s %s", cdo_get_stream_name(0), cdo_get_stream_name(1));
    }
    if (toID == -1) cdo_abort("In-situ temperature not found!");

    auto gridID = vlistGrid(vlistID1, 0);
    gridsize = vlist_check_gridsize(vlistID1);

    auto const &varSAO = varList1.vars[saoID];
    auto const &varTO = varList1.vars[toID];
    zaxisID = varSAO.zaxisID;

    if (varSAO.nlevels != varTO.nlevels) cdo_abort("temperature and salinity have different number of levels!");
    numLevels = varSAO.nlevels;

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

    to.resize(numLevels);
    sao.resize(numLevels);
    rho.resize(numLevels);

    for (int levelID = 0; levelID < numLevels; ++levelID)
    {
      to[levelID].resize(gridsize);
      sao[levelID].resize(gridsize);
      rho[levelID].resize(gridsize);
      to[levelID].missval = varTO.missval;
      sao[levelID].missval = varSAO.missval;
      rho[levelID].missval = to[levelID].missval;
    }

    int datatype = CDI_DATATYPE_FLT32;
    if (varTO.dataType == CDI_DATATYPE_FLT64 && varSAO.dataType == CDI_DATATYPE_FLT64) datatype = CDI_DATATYPE_FLT64;

    auto vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    auto varID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);
    vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(18, 255, 255));
    cdiDefKeyString(vlistID2, varID2, CDI_KEY_NAME, "rhopoto");
    cdiDefKeyString(vlistID2, varID2, CDI_KEY_LONGNAME, "Sea water potential density");
    cdiDefKeyString(vlistID2, varID2, CDI_KEY_STDNAME, "sea_water_potential_density");
    cdiDefKeyString(vlistID2, varID2, CDI_KEY_UNITS, "kg m-3");
    vlistDefVarMissval(vlistID2, varID2, rho[0].missval);
    vlistDefVarDatatype(vlistID2, varID2, datatype);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
    vlistDestroy(vlistID2);
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
        if (varID == toID) cdo_read_field(streamID1, to[levelID]);
        if (varID == saoID) cdo_read_field(streamID1, sao[levelID]);
      }

      calc_rhopot(gridsize, numLevels, pressure, to, sao, rho);

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        field_num_mv(rho[levelID]);
        cdo_def_field(streamID2, 0, levelID);
        cdo_write_field(streamID2, rho[levelID]);
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
