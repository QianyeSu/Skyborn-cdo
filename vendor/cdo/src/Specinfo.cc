/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Specinfo specinfo  Spectral information
*/

#include "process_int.h"
#include <mpim_grid.h>
#include <algorithm>

#define NTR2NSP(ntr) ((ntr + 1) * (ntr + 2))
#define NSP2NTR(nsp) ((long) ((((std::sqrt((double) (4 * nsp + 1))) - 3) / 2)))
#define NGP2NLEVEL(ngp) ((long) (log10(((double) ngp) / 80.) / std::log10(4.)))
#define NGP_ICON(nrooti, nlevel) ((long) (20 * nrooti * nrooti * ipow(4, nlevel)))
/*#define NGP_GME(ni)           ((ni+1)*(ni+1)*10)*/
#define NGP_GME(ni) (2 + ni * ni * 10)
#define NGP2NI(ngp) ((long) std::sqrt((double) ngp / 10.) - 1)

static void
fac(long nlonin, long *nlonout, int *ierr)
{
  long m = nlonin;

  while (m % 2 == 0) { m = m / 2; }
  while (m % 3 == 0) { m = m / 3; }
  while (m % 5 == 0) { m = m / 5; }

  if (m == 1)
  {
    *nlonout = nlonin;
    *ierr = 0;
  }
  else
  {
    *nlonout = nlonin + 1;
    *ierr = 1;
  }

  return;
}

static long
nlat2nlon(long nlat)
{
  if (nlat == 0) cdo_abort("nlat = 0!");

  long nlon = 2 * nlat;

  long m;
  int ierr;
  fac(nlon, &m, &ierr);
  /* adjust till fft is possible */
  while (ierr != 0)
  {
    nlon = m;
    /* correct here nlon so that nlat keeps always even */
    while (nlon % 4 != 0) nlon++;
    fac(nlon, &m, &ierr);
  }

  return nlon;
}

long
ngp2ntr(long ngp)
{
  long ntr = (long) std::lround(std::sqrt(0.25 + ngp) - 1.5);
  long nlonl = nlat_to_nlon(ntr_to_nlat_linear(ntr));
  long nlatl = nlonl / 2;

  ntr = (2 * nlatl - 1) / 2;

  return ntr;
}

static long
ipow(long i1, long i2)
{
  long i3 = 1;

  for (long i = 0; i < i2; ++i) i3 *= i1;

  return i3;
}

static constexpr long NiMax = 12;

static void
lookup_ni(long nsp, long *nroot, long *ni)
{
  long tbl2[NiMax], tbl3[NiMax], tbl5[NiMax];
  long d2 = 0, n2 = 0, d3 = 0, n3 = 0, d5 = 0, n5 = 0;

  for (long i = 0; i < NiMax; ++i)
  {
    tbl2[i] = 10 * 2 * 2 * ipow(4, (i + 1)) + 2;
    tbl3[i] = 10 * 3 * 3 * ipow(4, (i + 1)) + 2;
    tbl5[i] = 10 * 5 * 5 * ipow(4, (i + 1)) + 2;
  }

  for (long i = 0; i < NiMax; ++i)
    if (tbl2[i] >= nsp)
    {
      n2 = i;
      d2 = tbl2[n2] - nsp;
      break;
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl3[i] >= nsp)
    {
      n3 = i;
      d3 = tbl3[n3] - nsp;
      break;
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl5[i] >= nsp)
    {
      n5 = i;
      d5 = tbl5[n5] - nsp;
      break;
    }

  long d = d2;
  if (d3 < d) d = d3;
  if (d5 < d) d = d5;

  if (d == d2)
  {
    *nroot = 2;
    *ni = 2 * ipow(2, n2 + 1);
  }
  else if (d == d3)
  {
    *nroot = 3;
    *ni = 3 * ipow(2, n3 + 1);
  }
  else if (d == d5)
  {
    *nroot = 5;
    *ni = 5 * ipow(2, n5 + 1);
  }
}

static void
lookup_rl(long nsp, long *nroot, long *nlevel)
{
  long tbl2[NiMax], tbl3[NiMax], tbl5[NiMax];
  long d2 = 0, n2 = 0, d3 = 0, n3 = 0, d5 = 0, n5 = 0;

  for (long i = 0; i < NiMax; ++i)
  {
    tbl2[i] = 20 * 2 * 2 * ipow(4, (i + 1));
    tbl3[i] = 20 * 3 * 3 * ipow(4, (i + 1));
    tbl5[i] = 20 * 5 * 5 * ipow(4, (i + 1));
  }

  for (long i = 0; i < NiMax; ++i)
    if (tbl2[i] >= nsp)
    {
      n2 = i;
      d2 = tbl2[n2] - nsp;
      break;
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl3[i] >= nsp)
    {
      n3 = i;
      d3 = tbl3[n3] - nsp;
      break;
    }

  for (long i = 0; i < NiMax; ++i)
    if (tbl5[i] >= nsp)
    {
      n5 = i;
      d5 = tbl5[n5] - nsp;
      break;
    }

  long d = d2;
  if (d3 < d) d = d3;
  if (d5 < d) d = d5;

  if (d == d2)
  {
    *nroot = 2;
    *nlevel = n2 + 1;
  }
  else if (d == d3)
  {
    *nroot = 3;
    *nlevel = n3 + 1;
  }
  else if (d == d5)
  {
    *nroot = 5;
    *nlevel = n5 + 1;
  }
}

class Specinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Specinfo",
    .operators = { { "specinfo" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<Specinfo> registration = RegisterEntry<Specinfo>(module);

  char arg[128], *parg;
  std::string argument;
  struct GridSpecifications
  {
    bool nout = false;
    long ntr = 0;  // num truncations
    long nsp = 0;  // num spectral_coefficients
    long nlat = 0;
    long nlon = 0;
    long ngp = 0;  // num_grid_points
    long ni = 0;   // number intersections
    long ngp_gme = 0;

    long nlevel = 0;
    long ngp_icon = 0;
    long nrootg = 0;
    long nrooti = 0;
  };
  GridSpecifications grid_specs1;
  GridSpecifications grid_specs2;
  GridSpecifications grid_specs3;

private:
  void
  N_O(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2, GridSpecifications &p_grid_specs3, long nlon_offset = 0)
  {
    parg = &arg[1];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);

    p_grid_specs1.nlat = 2 * atoi(parg);
    p_grid_specs2.nlat = p_grid_specs1.nlat;
    p_grid_specs3.nlat = p_grid_specs1.nlat;

    p_grid_specs1.nlon = nlat2nlon(p_grid_specs1.nlat) + nlon_offset;
    p_grid_specs2.nlon = nlat2nlon(p_grid_specs2.nlat) + nlon_offset;
    p_grid_specs3.nlon = nlat2nlon(p_grid_specs3.nlat) + nlon_offset;

    p_grid_specs1.nlat = p_grid_specs1.nlon / 2;
    p_grid_specs2.nlat = p_grid_specs2.nlon / 2;
    p_grid_specs3.nlat = p_grid_specs3.nlon / 2;

    p_grid_specs1.ntr = (p_grid_specs1.nlat * 2 - 1) / 3;
    p_grid_specs2.ntr = (p_grid_specs2.nlat * 2 - 1) / 2;
    p_grid_specs3.ntr = (p_grid_specs3.nlat * 2 - 1) / 4;

    p_grid_specs1.ngp = p_grid_specs1.nlon * p_grid_specs1.nlat;
    p_grid_specs2.ngp = p_grid_specs2.nlon * p_grid_specs2.nlat;
    p_grid_specs3.ngp = p_grid_specs3.nlon * p_grid_specs3.nlat;

    p_grid_specs1.nsp = NTR2NSP(p_grid_specs1.ntr);
    p_grid_specs2.nsp = NTR2NSP(p_grid_specs2.ntr);
    p_grid_specs3.nsp = NTR2NSP(p_grid_specs3.ntr);

    lookup(p_grid_specs1);
    lookup(p_grid_specs2);
    lookup(p_grid_specs3);
  }
  void
  lookup(GridSpecifications &p_grid_specs1)
  {
    lookup_ni(p_grid_specs1.nsp, &p_grid_specs1.nrootg, &p_grid_specs1.ni);
    lookup_rl(p_grid_specs1.nsp, &p_grid_specs1.nrooti, &p_grid_specs1.nlevel);
    p_grid_specs1.nout = true;
  }

  void
  T_TL_TC(GridSpecifications &p_grid_specs, std::function<long(long)> f)
  {
    p_grid_specs.nsp = NTR2NSP(p_grid_specs.ntr);
    p_grid_specs.nlat = f(p_grid_specs.ntr);
    p_grid_specs.nlon = nlat_to_nlon(p_grid_specs.nlat);
    p_grid_specs.ngp = p_grid_specs.nlon * p_grid_specs.nlat;

    lookup(p_grid_specs);
  }
  void
  T(GridSpecifications &p_grid_specs, std::function<long(long)> f)
  {
    parg = &arg[1];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs.ntr = atoi(parg);
    T_TL_TC(p_grid_specs, f);
  }
  void
  TL_TC(GridSpecifications &p_grid_specs, std::function<long(long)> f)
  {
    parg = &arg[2];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs.ntr = atoi(parg);
    T_TL_TC(p_grid_specs, f);
  }

  void
  NI(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2)
  {
    parg = &arg[2];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs1.ni = atoi(parg);
    p_grid_specs2.ni = p_grid_specs1.ni;

    p_grid_specs1.ngp_gme = NGP_GME(p_grid_specs1.ni);
    p_grid_specs2.ngp_gme = NGP_GME(p_grid_specs2.ni);

    p_grid_specs1.ntr = ngp2ntr(p_grid_specs1.ngp_gme);
    p_grid_specs1.nsp = NTR2NSP(p_grid_specs1.ntr);

    p_grid_specs1.ntr = NSP2NTR(p_grid_specs1.nsp);
    p_grid_specs2.ntr = p_grid_specs1.ntr;

    p_grid_specs1.nlat = ntr_to_nlat(p_grid_specs1.ntr);
    p_grid_specs2.nlat = ntr_to_nlat_linear(p_grid_specs2.ntr);

    p_grid_specs1.nlon = nlat_to_nlon(p_grid_specs1.nlat);
    p_grid_specs2.nlon = nlat_to_nlon(p_grid_specs2.nlat);

    p_grid_specs1.nlat = p_grid_specs1.nlon / 2;
    p_grid_specs2.nlat = p_grid_specs2.nlon / 2;

    /* lookup_ni(p_grid_specs1.nsp, &grid_specs1.nrootg, &p_grid_specs1.ni); */
    lookup_rl(p_grid_specs1.nsp, &p_grid_specs1.nrooti, &p_grid_specs1.nlevel);

    p_grid_specs2.nrootg = p_grid_specs1.nrootg;
    p_grid_specs2.ni = p_grid_specs1.ni;
    p_grid_specs2.nrooti = p_grid_specs1.nrooti;
    p_grid_specs2.nlevel = p_grid_specs1.nlevel;

    p_grid_specs1.nout = true;
    p_grid_specs2.nout = true;
  }

  void
  NLON_NLAT(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2, GridSpecifications &p_grid_specs3)
  {
    p_grid_specs1.nlat = p_grid_specs1.nlon / 2;
    p_grid_specs2.nlat = p_grid_specs2.nlon / 2;
    p_grid_specs3.nlat = p_grid_specs3.nlon / 2;

    p_grid_specs1.ntr = (p_grid_specs1.nlat * 2 - 1) / 3;
    p_grid_specs2.ntr = (p_grid_specs2.nlat * 2 - 1) / 2;
    p_grid_specs3.ntr = (p_grid_specs3.nlat * 2 - 1) / 4;

    p_grid_specs1.ngp = p_grid_specs1.nlon * p_grid_specs1.nlat;
    p_grid_specs2.ngp = p_grid_specs2.nlon * p_grid_specs2.nlat;
    p_grid_specs3.ngp = p_grid_specs3.nlon * p_grid_specs3.nlat;

    p_grid_specs1.nsp = NTR2NSP(p_grid_specs1.ntr);
    p_grid_specs2.nsp = NTR2NSP(p_grid_specs2.ntr);
    p_grid_specs3.nsp = NTR2NSP(p_grid_specs3.ntr);

    lookup(p_grid_specs1);
    lookup(p_grid_specs2);
    lookup(p_grid_specs3);
  }

  void
  NLON(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2, GridSpecifications &p_grid_specs3)
  {
    parg = &arg[4];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs1.nlon = atoi(parg);
    p_grid_specs2.nlon = p_grid_specs1.nlon;
    p_grid_specs3.nlon = p_grid_specs1.nlon;

    p_grid_specs1.nlat = p_grid_specs1.nlon / 2;
    p_grid_specs2.nlat = p_grid_specs2.nlon / 2;
    p_grid_specs3.nlat = p_grid_specs3.nlon / 2;

    p_grid_specs1.nlon = nlat2nlon(p_grid_specs1.nlat);
    p_grid_specs2.nlon = nlat2nlon(p_grid_specs2.nlat);
    p_grid_specs3.nlon = nlat2nlon(p_grid_specs3.nlat);

    NLON_NLAT(p_grid_specs1, p_grid_specs2, p_grid_specs3);
  }

  void
  NLAT(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2, GridSpecifications &p_grid_specs3)
  {
    parg = &arg[4];
    if (*parg == '=') parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs1.nlat = atoi(parg);
    p_grid_specs2.nlat = p_grid_specs1.nlat;
    p_grid_specs3.nlat = p_grid_specs1.nlat;

    p_grid_specs1.nlon = nlat2nlon(p_grid_specs1.nlat);
    p_grid_specs2.nlon = nlat2nlon(p_grid_specs2.nlat);
    p_grid_specs3.nlon = nlat2nlon(p_grid_specs3.nlat);
    NLON_NLAT(p_grid_specs1, p_grid_specs2, p_grid_specs3);
  }

  void
  ICON(GridSpecifications &p_grid_specs1, GridSpecifications &p_grid_specs2)
  {
    parg = &arg[4];
    if (*parg != 'R') cdo_abort("Wrong parameter: %s", arg);
    parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs1.nrooti = atoi(parg);
    p_grid_specs2.nrooti = p_grid_specs1.nrooti;
    while (std::isdigit((int) *parg)) parg++;
    if (*parg != 'B') cdo_abort("Wrong parameter: %s", arg);
    parg++;
    if (!std::isdigit((int) *parg)) cdo_abort("Wrong parameter: %s", arg);
    p_grid_specs1.nlevel = atoi(parg);
    p_grid_specs2.nlevel = p_grid_specs1.nlevel;
    p_grid_specs1.ngp_icon = NGP_ICON(p_grid_specs1.nrooti, p_grid_specs1.nlevel);
    p_grid_specs2.ngp_icon = NGP_ICON(p_grid_specs1.nrooti, p_grid_specs2.nlevel);

    p_grid_specs1.ntr = ngp2ntr(p_grid_specs1.ngp_icon);
    p_grid_specs1.nsp = NTR2NSP(p_grid_specs1.ntr);
    p_grid_specs1.ntr = NSP2NTR(p_grid_specs1.nsp);
    p_grid_specs2.ntr = p_grid_specs1.ntr;

    p_grid_specs1.nlat = ntr_to_nlat(p_grid_specs1.ntr);
    p_grid_specs1.nlon = nlat_to_nlon(p_grid_specs1.nlat);
    p_grid_specs1.nlat = p_grid_specs1.nlon / 2;

    p_grid_specs2.nlat = ntr_to_nlat_linear(p_grid_specs2.ntr);
    p_grid_specs2.nlon = nlat_to_nlon(p_grid_specs2.nlat);
    p_grid_specs2.nlat = p_grid_specs2.nlon / 2;

    lookup_ni(p_grid_specs1.nsp, &p_grid_specs1.nrootg, &p_grid_specs1.ni);
    /* lookup_rl(p_grid_specs1.nsp, &p_grid_specs1.nrooti, &p_grid_specs1.nlevel);*/

    p_grid_specs2.nrootg = p_grid_specs1.nrootg;
    p_grid_specs2.ni = p_grid_specs1.ni;
    p_grid_specs2.nrooti = p_grid_specs1.nrooti;
    p_grid_specs2.nlevel = p_grid_specs1.nlevel;

    p_grid_specs1.nout = true;
    p_grid_specs2.nout = true;
  }

public:
  void
  init() override
  {
    operator_input_arg("Txx, TLxx, NLON=xx, NLAT=xx, NIxx or ICONRyyLxx");

    long len = cdo_operator_argv(0).size();

    if ((len + 1) >= 128) cdo_abort("Parameter string too large!");

    for (long i = 0; i < len; ++i) arg[i] = std::toupper(cdo_operator_argv(0)[i]);
    arg[len] = 0;

    argument = std::string(cdo_operator_argv(0));
    std::ranges::transform(argument, argument.begin(), ::toupper);
  }

  void
  run() override
  {
    if (argument.substr(0, 2) == "TL") { TL_TC(grid_specs2, ntr_to_nlat_linear); }
    else if (argument.substr(0, 2) == "TC") { TL_TC(grid_specs3, ntr_to_nlat_cubic); }
    else if (argument.substr(0, 1) == "T") { T(grid_specs1, ntr_to_nlat); }
    else if (argument.substr(0, 2) == "NI") { NI(grid_specs1, grid_specs2); }
    else if (argument.substr(0, 4) == "NLON") { NLON(grid_specs1, grid_specs2, grid_specs3); }
    else if (argument.substr(0, 4) == "NLAT") { NLAT(grid_specs1, grid_specs2, grid_specs3); }
    else if (argument.substr(0, 1) == "N") { N_O(grid_specs1, grid_specs2, grid_specs3); }
    else if (argument.substr(0, 1) == "O") { N_O(grid_specs1, grid_specs2, grid_specs3, 16); }
    else if (argument.substr(0, 4) == "ICON") { ICON(grid_specs1, grid_specs2); }
    else { cdo_abort("Unsupported parameter: %s", arg); }

    grid_specs1.nsp = NTR2NSP(grid_specs1.ntr);
    grid_specs2.nsp = NTR2NSP(grid_specs2.ntr);
    grid_specs3.nsp = NTR2NSP(grid_specs3.ntr);
    grid_specs1.ngp = grid_specs1.nlon * grid_specs1.nlat;
    grid_specs2.ngp = grid_specs2.nlon * grid_specs2.nlat;
    grid_specs3.ngp = grid_specs3.nlon * grid_specs3.nlat;
    grid_specs1.ngp_gme = NGP_GME(grid_specs1.ni);
    grid_specs2.ngp_gme = NGP_GME(grid_specs2.ni);
    grid_specs3.ngp_gme = NGP_GME(grid_specs3.ni);
    grid_specs1.ngp_icon = NGP_ICON(grid_specs1.nrooti, grid_specs1.nlevel);
    grid_specs2.ngp_icon = NGP_ICON(grid_specs2.nrooti, grid_specs2.nlevel);
    grid_specs3.ngp_icon = NGP_ICON(grid_specs3.nrooti, grid_specs3.nlevel);

    fprintf(stdout, "truncation     nsp  nlon  nlat      ngp  gme    ngp_gme  icon   ngp_icon\n");

    if (grid_specs2.nout)
      fprintf(stdout, "   TL%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", grid_specs2.ntr, grid_specs2.nsp,
              grid_specs2.nlon, grid_specs2.nlat, grid_specs2.ngp, grid_specs2.ni, grid_specs2.ngp_gme, grid_specs2.nrooti,
              grid_specs2.nlevel, grid_specs2.ngp_icon);

    if (grid_specs1.nout)
      fprintf(stdout, "   TQ%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", grid_specs1.ntr, grid_specs1.nsp,
              grid_specs1.nlon, grid_specs1.nlat, grid_specs1.ngp, grid_specs1.ni, grid_specs1.ngp_gme, grid_specs1.nrooti,
              grid_specs1.nlevel, grid_specs1.ngp_icon);

    if (grid_specs3.nout)
      fprintf(stdout, "   TC%-4ld %8ld %5ld %5ld %8ld  ni%ld %8ld  R%ldB%02ld  %8ld\n", grid_specs3.ntr, grid_specs3.nsp,
              grid_specs3.nlon, grid_specs3.nlat, grid_specs3.ngp, grid_specs3.ni, grid_specs3.ngp_gme, grid_specs3.nrooti,
              grid_specs3.nlevel, grid_specs3.ngp_icon);
  }
  void
  close() override
  {
  }
};
