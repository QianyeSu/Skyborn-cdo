/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdlib>
#include <cassert>
#include <cdi.h>

#include "cdo_options.h"
#include "field_functions.h"
#include "cdo_output.h"
#include "cdo_omp.h"
#include "interpol.h"
#include "expr.h"
#include "expr_fun.h"
#include "expr_yacc.hh"
#include "cdo_zaxis.h"

static const char *const ExIn[] = { "expr", "init" };
static const char *const tmpVarName = "_tmp_";
static int pointID = -1;
static int zonalID = -1;
static int surfaceID = -1;

enum
{
  FT_STD,
  FT_CONST,
  FT_FLD,
  FT_ZON,
  FT_VERT,
  FT_REMAP,
  FT_COORD,
  FT_1C,
  FT_2C,
  FT_0
};

// clang-format off
static double f_float(double x) { return (float) (x); }
static double f_int(double x) { return (int) (x); }
static double f_nint(double x) { return std::round(x); }
static double f_rand(double x) { (void)x; return ((double) std::rand()) / ((double) RAND_MAX); }
static double f_sqr(double x) { return x * x; }
static double f_rad(double x) { return x * M_PI / 180.0; }
static double f_deg(double x) { return x * 180.0 / M_PI; }
static double f_isMissval(double x) { (void)x; return 0.0; }

static double pt_definegrid(const ParamEntry *const p) { printf("pt_definegrid: >%s<\n", p->name.c_str()); return -1;/*return cdo_define_grid(xxx);*/ }
static double pt_grid(const ParamEntry *const p) { return p->gridID; }
static double pt_ngp(const ParamEntry *const p) { return p->ngp; }
static double pt_nlev(const ParamEntry *const p) { return p->nlev; }
static double pt_size(const ParamEntry *const p) { return p->ngp * p->nlev; }
static double pt_missval(const ParamEntry *const p) { return p->missval; }
static double ts_ctimestep(const double *const data) { return std::lround(data[CoordIndex::TIMESTEP]); }
static double ts_cdate(const double *const data) { return std::lround(data[CoordIndex::DATE]); }
static double ts_ctime(const double *const data) { return std::lround(data[CoordIndex::TIME]); }
static double ts_cdeltat(const double *const data) { return data[CoordIndex::DELTAT]; }
static double ts_cday(const double *const data) { return data[CoordIndex::DAY]; }
static double ts_cmonth(const double *const data) { return data[CoordIndex::MONTH]; }
static double ts_cyear(const double *const data) { return data[CoordIndex::YEAR]; }
static double ts_csecond(const double *const data) { return data[CoordIndex::SECOND]; }
static double ts_cminute(const double *const data) { return data[CoordIndex::MINUTE]; }
static double ts_chour(const double *const data) { return data[CoordIndex::HOUR]; }
static double ts_cdoy(const double *const data) { return data[CoordIndex::DOY]; }
static double ts_cdpy(const double *const data) { return data[CoordIndex::DPY]; }
// clang-format on

namespace
{
struct ExprFuncEntry
{
  int type{};
  int flag{};
  const std::string name{};  // function name
  void (*func)(void){};      // pointer to function
};
}  // namespace

// clang-format off
static const ExprFuncEntry funcSymbolTable[] = {
  // scalar functions   y=func(x)
  { FT_STD, 0, "abs",   (void (*)(void))((double (*)(double)) std::fabs) },  // math functions could be inlined
  { FT_STD, 0, "floor", (void (*)(void))((double (*)(double)) std::floor) },
  { FT_STD, 0, "ceil",  (void (*)(void))((double (*)(double)) std::ceil) },
  { FT_STD, 0, "sqrt",  (void (*)(void))((double (*)(double)) std::sqrt) },
  { FT_STD, 0, "exp",   (void (*)(void))((double (*)(double)) std::exp) },
  { FT_STD, 0, "erf",   (void (*)(void))((double (*)(double)) std::erf) },
  { FT_STD, 0, "log",   (void (*)(void))((double (*)(double)) std::log) },
  { FT_STD, 0, "ln",    (void (*)(void))((double (*)(double)) std::log) },
  { FT_STD, 0, "log10", (void (*)(void))((double (*)(double)) std::log10) },
  { FT_STD, 0, "sin",   (void (*)(void))((double (*)(double)) std::sin) },
  { FT_STD, 0, "cos",   (void (*)(void))((double (*)(double)) std::cos) },
  { FT_STD, 0, "tan",   (void (*)(void))((double (*)(double)) std::tan) },
  { FT_STD, 0, "sinh",  (void (*)(void))((double (*)(double)) std::sinh) },
  { FT_STD, 0, "cosh",  (void (*)(void))((double (*)(double)) std::cosh) },
  { FT_STD, 0, "tanh",  (void (*)(void))((double (*)(double)) std::tanh) },
  { FT_STD, 0, "asin",  (void (*)(void))((double (*)(double)) std::asin) },
  { FT_STD, 0, "acos",  (void (*)(void))((double (*)(double)) std::acos) },
  { FT_STD, 0, "atan",  (void (*)(void))((double (*)(double)) std::atan) },
  { FT_STD, 0, "asinh", (void (*)(void))((double (*)(double)) std::asinh) },
  { FT_STD, 0, "acosh", (void (*)(void))((double (*)(double)) std::acosh) },
  { FT_STD, 0, "atanh", (void (*)(void))((double (*)(double)) std::atanh) },
  { FT_STD, 0, "gamma", (void (*)(void))((double (*)(double)) std::tgamma) },
  // scalar functions   z=func(x,y)
  { FT_STD, 0, "mod",   (void (*)(void))((double (*)(double, double)) std::fmod) },
  { FT_STD, 0, "min",   (void (*)(void))((double (*)(double, double)) std::fmin) },
  { FT_STD, 0, "max",   (void (*)(void))((double (*)(double, double)) std::fmax) },
  { FT_STD, 0, "pow",   (void (*)(void))((double (*)(double, double)) std::pow) },
  { FT_STD, 0, "hypot", (void (*)(void))((double (*)(double, double)) std::hypot) },
  { FT_STD, 0, "atan2", (void (*)(void))((double (*)(double, double)) std::atan2) },
  // scalar functions   y=func(x)
  { FT_STD, 0, "float",     reinterpret_cast<void (*)(void)>(&f_float) },
  { FT_STD, 0, "int",       reinterpret_cast<void (*)(void)>(&f_int) },
  { FT_STD, 0, "nint",      reinterpret_cast<void (*)(void)>(&f_nint) },
  { FT_STD, 0, "rand",      reinterpret_cast<void (*)(void)>(&f_rand) },
  { FT_STD, 0, "sqr",       reinterpret_cast<void (*)(void)>(&f_sqr) },
  { FT_STD, 0, "rad",       reinterpret_cast<void (*)(void)>(&f_rad) },
  { FT_STD, 0, "deg",       reinterpret_cast<void (*)(void)>(&f_deg) },
  { FT_STD, 0, "isMissval", reinterpret_cast<void (*)(void)>(&f_isMissval) },

  // constant functions  c=func(varname)
  { FT_CONST, 0, "define_grid", reinterpret_cast<void (*)(void)>(&pt_definegrid) }, // generate CDI grid ID from name or file
  { FT_CONST, 0, "grid",        reinterpret_cast<void (*)(void)>(&pt_grid) },       // CDI grid ID
  { FT_CONST, 0, "ngp",         reinterpret_cast<void (*)(void)>(&pt_ngp) },        // number of horizontal grid points
  { FT_CONST, 0, "nlev",        reinterpret_cast<void (*)(void)>(&pt_nlev) },       // number of vertical levels
  { FT_CONST, 0, "size",        reinterpret_cast<void (*)(void)>(&pt_size) },       // ngp*nlev
  { FT_CONST, 0, "missval",     reinterpret_cast<void (*)(void)>(&pt_missval) },    // Returns the missing value of a variable

  // CDO field functions (Reduce grid to point)  varOut=fld<STAT>(varIn)
  { FT_FLD, 0, "fldmin",    reinterpret_cast<void (*)(void)>(&field_min) },
  { FT_FLD, 0, "fldmax",    reinterpret_cast<void (*)(void)>(&field_max) },
  { FT_FLD, 0, "fldrange",  reinterpret_cast<void (*)(void)>(&field_range) },
  { FT_FLD, 0, "fldsum",    reinterpret_cast<void (*)(void)>(&field_sum) },
  { FT_FLD, 1, "fldmean",   reinterpret_cast<void (*)(void)>(&field_meanw) },
  { FT_FLD, 1, "fldavg",    reinterpret_cast<void (*)(void)>(&field_avgw) },
  { FT_FLD, 1, "fldstd",    reinterpret_cast<void (*)(void)>(&field_stdw) },
  { FT_FLD, 1, "fldstd1",   reinterpret_cast<void (*)(void)>(&field_std1w) },
  { FT_FLD, 1, "fldvar",    reinterpret_cast<void (*)(void)>(&field_varw) },
  { FT_FLD, 1, "fldvar1",   reinterpret_cast<void (*)(void)>(&field_var1w) },
  { FT_FLD, 1, "fldskew",   reinterpret_cast<void (*)(void)>(&field_skew) },
  { FT_FLD, 1, "fldkurt",   reinterpret_cast<void (*)(void)>(&field_kurt) },
  { FT_FLD, 1, "fldmedian", reinterpret_cast<void (*)(void)>(&field_median) },
  { FT_FLD, 1, "fldcount",  reinterpret_cast<void (*)(void)>(&field_count) },

  // CDO zonal functions (Reduce grid to point)  varOut=zon<STAT>(varIn)
  { FT_ZON, 0, "zonmin",    reinterpret_cast<void (*)(void)>(&zonal_min) },
  { FT_ZON, 0, "zonmax",    reinterpret_cast<void (*)(void)>(&zonal_max) },
  { FT_ZON, 0, "zonrange",  reinterpret_cast<void (*)(void)>(&zonal_range) },
  { FT_ZON, 0, "zonsum",    reinterpret_cast<void (*)(void)>(&zonal_sum) },
  { FT_ZON, 0, "zonmean",   reinterpret_cast<void (*)(void)>(&zonal_mean) },
  { FT_ZON, 0, "zonavg",    reinterpret_cast<void (*)(void)>(&zonal_avg) },
  { FT_ZON, 0, "zonstd",    reinterpret_cast<void (*)(void)>(&zonal_std) },
  { FT_ZON, 0, "zonstd1",   reinterpret_cast<void (*)(void)>(&zonal_std1) },
  { FT_ZON, 0, "zonvar",    reinterpret_cast<void (*)(void)>(&zonal_var) },
  { FT_ZON, 0, "zonvar1",   reinterpret_cast<void (*)(void)>(&zonal_var1) },
  { FT_ZON, 0, "zonskew",   reinterpret_cast<void (*)(void)>(&zonal_skew) },
  { FT_ZON, 0, "zonkurt",   reinterpret_cast<void (*)(void)>(&zonal_kurt) },
  { FT_ZON, 0, "zonmedian", reinterpret_cast<void (*)(void)>(&zonal_median) },

  // CDO field functions (Reduce level to point)  varOut=vert<STAT>(varIn)
  { FT_VERT, 0, "vertmin",    reinterpret_cast<void (*)(void)>(&field_min) },
  { FT_VERT, 0, "vertmax",    reinterpret_cast<void (*)(void)>(&field_max) },
  { FT_VERT, 0, "vertrange",  reinterpret_cast<void (*)(void)>(&field_range) },
  { FT_VERT, 0, "vertsum",    reinterpret_cast<void (*)(void)>(&field_sum) },
  { FT_VERT, 1, "vertmean",   reinterpret_cast<void (*)(void)>(&field_meanw) },
  { FT_VERT, 1, "vertavg",    reinterpret_cast<void (*)(void)>(&field_avgw) },
  { FT_VERT, 1, "vertstd",    reinterpret_cast<void (*)(void)>(&field_stdw) },
  { FT_VERT, 1, "vertstd1",   reinterpret_cast<void (*)(void)>(&field_std1w) },
  { FT_VERT, 1, "vertvar",    reinterpret_cast<void (*)(void)>(&field_varw) },
  { FT_VERT, 1, "vertvar1",   reinterpret_cast<void (*)(void)>(&field_var1w) },
  { FT_VERT, 1, "vertskew",   reinterpret_cast<void (*)(void)>(&field_skew) },
  { FT_VERT, 1, "vertkurt",   reinterpret_cast<void (*)(void)>(&field_kurt) },
  { FT_VERT, 1, "vertmedian", reinterpret_cast<void (*)(void)>(&field_median) },
  // varOut=intgrid<METHOD>(varIn, gridOut)s
  { FT_REMAP, 0, "intgridbil", reinterpret_cast<void (*)(void)>(&intgrid_bil) },
  { FT_REMAP, 0, "intgridnn",  reinterpret_cast<void (*)(void)>(&intgrid_1nn) },
  // coord=func(varName)
  { FT_COORD, 0, "clon",       nullptr },
  { FT_COORD, 0, "clat",       nullptr },
  { FT_COORD, 0, "clev",       nullptr },
  { FT_COORD, 0, "clevidx",    nullptr },
  { FT_COORD, 0, "cthickness", nullptr },
  { FT_COORD, 0, "gridarea",   nullptr },
  { FT_COORD, 0, "gridweight", nullptr },
  { FT_COORD, 0, "gridindex",  nullptr },
  // c=func()
  { FT_0, 0, "ctimestep", reinterpret_cast<void (*)(void)>(&ts_ctimestep) },
  { FT_0, 0, "cdate",     reinterpret_cast<void (*)(void)>(&ts_cdate) },
  { FT_0, 0, "ctime",     reinterpret_cast<void (*)(void)>(&ts_ctime) },
  { FT_0, 0, "cdeltat",   reinterpret_cast<void (*)(void)>(&ts_cdeltat) },
  { FT_0, 0, "cday",      reinterpret_cast<void (*)(void)>(&ts_cday) },
  { FT_0, 0, "cmonth",    reinterpret_cast<void (*)(void)>(&ts_cmonth) },
  { FT_0, 0, "cyear",     reinterpret_cast<void (*)(void)>(&ts_cyear) },
  { FT_0, 0, "csecond",   reinterpret_cast<void (*)(void)>(&ts_csecond) },
  { FT_0, 0, "cminute",   reinterpret_cast<void (*)(void)>(&ts_cminute) },
  { FT_0, 0, "chour",     reinterpret_cast<void (*)(void)>(&ts_chour) },
  { FT_0, 0, "cdoy",      reinterpret_cast<void (*)(void)>(&ts_cdoy) },
  { FT_0, 0, "cdpy",      reinterpret_cast<void (*)(void)>(&ts_cdpy) },
  // varOut=sellevXXX(varIn,k)
  // varOut=sellevXXXrange(varIn,k1,k2)
  { FT_1C, 0, "sellevel",       nullptr },
  { FT_1C, 0, "sellevidx",      nullptr },
  { FT_2C, 0, "sellevelrange",  nullptr },
  { FT_2C, 0, "sellevidxrange", nullptr },
  // {FT_1C, 0, "gridindex", nullptr},
};
// clang-format on

static void
node_data_delete(nodeType *p)
{
  if (p)
  {
    if (p->param.data)
    {
      delete[] p->param.data;
      p->param.data = nullptr;
    }
  }
}

static void
node_delete(nodeType *p)
{
  if (p)
  {
    if (p->type == NodeEnum::typeVar) node_data_delete(p);
    delete p;
  }
}

static int
get_funcID(std::string const &fun)
{
  constexpr int funcTableSize = sizeof(funcSymbolTable) / sizeof(funcSymbolTable[0]);
  for (int funcID = 0; funcID < funcTableSize; funcID++)
    if (fun == funcSymbolTable[funcID].name) return funcID;

  cdo_abort("Function >%s< not available!", fun);

  return -1;
}

static constexpr bool
is_compare(int oper) noexcept
{
  return (oper == LEG || oper == GE || oper == LE || oper == EQ || oper == NE || oper == GT || oper == LT);
}

static void
param_meta_copy(ParamEntry &out, const ParamEntry &in)
{
  out.type = in.type;
  out.gridID = in.gridID;
  out.zaxisID = in.zaxisID;
  out.datatype = in.datatype;
  out.steptype = in.steptype;
  out.ngp = in.ngp;
  out.nlat = in.nlat;
  out.nlev = in.nlev;
  out.missval = in.missval;
  out.numMissVals = 0;
  out.coord = 0;
  out.hasMV = true;
  out.name = "";
  out.longname = "";
  out.units = "";
  out.data = nullptr;
}

static nodeType *
expr_con_var(int init, int oper, const nodeType *p1, const nodeType *p2)
{
  auto ngp = (p2->param.ngp > 0) ? p2->param.ngp : 1;
  auto nlev = (p2->param.nlev > 0) ? p2->param.nlev : 1;
  auto hasMV = (p2->param.numMissVals > 0);
  auto datatype = p2->param.datatype;
  auto missval2 = p2->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p2->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    auto n = ngp * nlev;
    p->param.data = new double[n];
    auto odat = p->param.data;
    const auto idat = p2->param.data;
    auto cval = p1->con().value;
    if (datatype == CDI_DATATYPE_FLT32 && is_compare(oper)) cval = (float) cval;

    oper_expr_con_var(oper, hasMV, n, missval2, odat, cval, idat);

    p->param.numMissVals = array_num_mv(n, odat, missval2);
  }

  return p;
}

static nodeType *
expr_var_con(int init, int oper, const nodeType *p1, const nodeType *p2)
{
  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto hasMV = (p1->param.numMissVals > 0);
  auto datatype = p1->param.datatype;
  auto missval1 = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    auto n = ngp * nlev;
    p->param.data = new double[n];
    auto odat = p->param.data;
    const auto idat = p1->param.data;
    auto cval = p2->con().value;
    if (datatype == CDI_DATATYPE_FLT32 && is_compare(oper)) cval = (float) cval;

    oper_expr_var_con(oper, hasMV, n, missval1, odat, idat, cval);

    p->param.numMissVals = array_num_mv(n, odat, missval1);
  }

  return p;
}

static nodeType *
expr_var_var(int init, int oper, nodeType *p1, nodeType *p2)
{
  auto px = p1;
  auto numMissVals1 = p1->param.numMissVals;
  auto numMissVals2 = p2->param.numMissVals;
  auto missval1 = p1->param.missval;
  auto missval2 = p2->param.missval;

  auto ngp1 = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto ngp2 = (p2->param.ngp > 0) ? p2->param.ngp : 1;
  auto ngp = ngp1;

  auto nlat1 = p1->param.nlat;
  auto nlat2 = p2->param.nlat;
  auto isZonal = false;

  if (ngp1 != ngp2)
  {
    if (ngp1 == 1 || ngp2 == 1)
    {
      if (ngp1 == 1)
      {
        ngp = ngp2;
        px = p2;
      }
    }
    else if (nlat1 == nlat2 && ngp1 > ngp2) { isZonal = true; }
    else
    {
      cdo_abort("%s: Number of grid points differ (%s[%zu] <-> %s[%zu])", __func__, p1->param.name, ngp1, p2->param.name, ngp2);
    }
  }

  auto nlev1 = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto nlev2 = (p2->param.nlev > 0) ? p2->param.nlev : 1;

  auto nlev = nlev1;
  if (nlev1 != nlev2)
  {
    if (nlev1 == 1 || nlev2 == 1)
    {
      if (nlev1 == 1)
      {
        nlev = nlev2;
        px = p2;
      }
    }
    else { cdo_abort("%s: Number of levels differ (%s[%zu] <-> %s[%zu])", __func__, p1->param.name, nlev1, p2->param.name, nlev2); }
  }

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, px->param);

  if (p->param.steptype == TIME_CONSTANT)
  {
    auto steptype1 = p1->param.steptype;
    auto steptype2 = p2->param.steptype;
    if (steptype1 != TIME_CONSTANT)
      p->param.steptype = steptype1;
    else if (steptype2 != TIME_CONSTANT)
      p->param.steptype = steptype2;
  }

  p->param.name = tmpVarName;
  // printf("%s %s numMissVals %zu %zu\n", p->u.var.name.c_str(), px->param.name.c_str(), numMissVals1, numMissVals2);

  if (!init)
  {
    p->param.data = new double[ngp * nlev];

    for (size_t k = 0; k < nlev; ++k)
    {
      auto loff = k * ngp;
      size_t loff1 = (nlev1 > 1) ? k * ngp1 : 0;
      size_t loff2 = (nlev2 > 1) ? k * ngp2 : 0;

      const auto idat1 = p1->param.data + loff1;
      const auto idat2 = p2->param.data + loff2;
      auto odat = p->param.data + loff;
      auto hasMV = (numMissVals1 > 0 || numMissVals2 > 0);

      if (ngp1 != ngp2)
      {
        if (isZonal)
        {
          auto nlon = ngp1 / nlat1;
          for (size_t j = 0; j < nlat1; ++j)
            oper_expr_var_con(oper, hasMV, nlon, missval1, odat + j * nlon, idat1 + j * nlon, idat2[j]);
        }
        else
        {
          (ngp2 == 1) ? oper_expr_var_con(oper, hasMV, ngp, missval1, odat, idat1, idat2[0])
                      : oper_expr_con_var(oper, hasMV, ngp, missval2, odat, idat1[0], idat2);
        }
      }
      else { oper_expr_var_var(oper, hasMV, ngp, missval1, missval2, odat, idat1, idat2); }
    }

    p->param.numMissVals = array_num_mv(ngp * nlev, p->param.data, missval1);
  }

  return p;
}

static void
ex_copy_var(int init, nodeType *p2, const nodeType *p1)
{
  auto copyConstValue1 = (p1->param.ngp == 0 && p1->param.nlev == 0);
  auto ngp1 = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev1 = (p1->param.nlev > 0) ? p1->param.nlev : 1;

  if (Options::cdoVerbose)
  {
    if (copyConstValue1)
      cdo_print("\t%s\tcopy\t%s[N%zu][L%zu] = %s", ExIn[init], p2->param.name, p2->param.ngp, p2->param.nlev, p1->param.name);
    else
      cdo_print("\t%s\tcopy\t%s[N%zu][L%zu] = %s[N%zu][L%zu]", ExIn[init], p2->param.name, p2->param.ngp, p2->param.nlev,
                p1->param.name, p1->param.ngp, p1->param.nlev);
  }

  auto ngp2 = (p2->param.ngp > 0) ? p2->param.ngp : 1;
  auto nlev2 = (p2->param.nlev > 0) ? p2->param.nlev : 1;

  if (!copyConstValue1 && ngp2 != ngp1)
    cdo_abort("%s: Number of grid points differ (%s[N%zu] = %s[N%zu])", __func__, p2->param.name, ngp2, p1->param.name, ngp1);

  if (!copyConstValue1 && nlev2 != nlev1)
    cdo_abort("%s: Number of levels differ (%s[L%zu] = %s[L%zu])", __func__, p2->param.name, nlev2, p1->param.name, nlev1);

  if (!init)
  {
    if (copyConstValue1) { std::ranges::fill_n(p2->param.data, ngp2 * nlev2, p1->param.data[0]); }
    else
    {
      array_copy(ngp2 * nlev2, p1->param.data, p2->param.data);
      p2->param.missval = p1->param.missval;
      p2->param.numMissVals = p1->param.numMissVals;
    }
  }
}

static void
ex_copy_con(int init, nodeType *p2, const nodeType *p1)
{
  auto cval = p1->con().value;

  auto ngp = p2->param.ngp;
  auto nlev = p2->param.nlev;

  if (Options::cdoVerbose)
  {
    (ngp == 0 && nlev == 0) ? cdo_print("\t%s\tcopy\t%s = %g", ExIn[init], p2->param.name, cval)
                            : cdo_print("\t%s\tcopy\t%s[N%zu][L%zu] = %g", ExIn[init], p2->param.name, ngp, nlev, cval);
  }

  if (ngp == 0 && nlev == 0)
  {
    ngp = 1;
    nlev = 1;
  }

  assert(ngp > 0);
  assert(nlev > 0);

  if (!init)
  {
    assert(p2->param.data != nullptr);
    std::ranges::fill_n(p2->param.data, ngp * nlev, cval);
  }
}

static void
ex_copy(int init, nodeType *p2, const nodeType *p1)
{
  (p1->type == NodeEnum::typeCon) ? ex_copy_con(init, p2, p1) : ex_copy_var(init, p2, p1);
}

static nodeType *
expr(int init, int oper, nodeType *p1, nodeType *p2)
{
  if (p1 == nullptr || p2 == nullptr) return nullptr;

  const char *coper = "???";

  if (Options::cdoVerbose)
  {
    switch (oper)
    {
      case LT: coper = "<"; break;
      case GT: coper = ">"; break;
      case LE: coper = "<="; break;
      case GE: coper = ">="; break;
      case NE: coper = "!="; break;
      case EQ: coper = "=="; break;
      case LEG: coper = "<=>"; break;
      case AND: coper = "&&"; break;
      case OR: coper = "||"; break;
      case '^': coper = "^"; break;
      case '+': coper = "+"; break;
      case '-': coper = "-"; break;
      case '*': coper = "*"; break;
      case '/': coper = "/"; break;
      default: cdo_abort("Internal error, expr operator %d not implemented!", oper);
    }
  }

  nodeType *p = nullptr;

  if (p1->type == NodeEnum::typeVar && p2->type == NodeEnum::typeVar)
  {
    p = expr_var_var(init, oper, p1, p2);
    if (Options::cdoVerbose)
      cdo_print("\t%s\tarith\t%s[N%zu][L%zu] = %s %s %s", ExIn[init], p->var().name, p->param.ngp, p->param.nlev, p1->var().name,
                coper, p2->var().name);
  }
  else if (p1->type == NodeEnum::typeVar && p2->type == NodeEnum::typeCon)
  {
    p = expr_var_con(init, oper, p1, p2);
    if (Options::cdoVerbose)
      cdo_print("\t%s\tarith\t%s[N%zu][L%zu] = %s %s %g", ExIn[init], p->var().name, p->param.ngp, p->param.nlev, p1->var().name,
                coper, p2->con().value);
  }
  else if (p1->type == NodeEnum::typeCon && p2->type == NodeEnum::typeVar)
  {
    p = expr_con_var(init, oper, p1, p2);
    if (Options::cdoVerbose)
      cdo_print("\t%s\tarith\t%s[N%zu][L%zu] = %g %s %s", ExIn[init], p->var().name, p->param.ngp, p->param.nlev, p1->con().value,
                coper, p2->var().name);
  }
  else if (p1->type == NodeEnum::typeCon && p2->type == NodeEnum::typeCon)
  {
    p = expr_con_con(oper, p1, p2);
    if (Options::cdoVerbose)
      cdo_print("\t%s\tarith\t%g = %g %s %g", ExIn[init], p->con().value, p1->con().value, coper, p2->con().value);
  }
  else
    cdo_abort("Internal problem!");

  if (p1->isTmpObj) node_delete(p1);
  if (p2->isTmpObj) node_delete(p2);

  return p;
}

static nodeType *
ex_fun0_con(int init, int funcID, double *data)
{
  auto &funcEntry = funcSymbolTable[funcID];
  if (funcEntry.type != FT_0) cdo_abort("Function %s not available for constant values!", funcEntry.name);

  if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s", ExIn[init], funcEntry.name);

  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  auto exprfunc = (double (*)(const double *)) funcEntry.func;
  constexpr double initDummy = -99.0;
  p->v = conNodeType(init ? initDummy : exprfunc(data));

  return p;
}

static nodeType *
ex_fun_con(int funcID, nodeType *p1)
{
  auto &funcEntry = funcSymbolTable[funcID];
  if (funcEntry.type != FT_STD) cdo_abort("Function %s not available for constant values!", funcEntry.name);

  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  auto exprfunc = (double (*)(double)) funcEntry.func;
  p->v = conNodeType(exprfunc(p1->con().value));

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static void
ex_fun_std(int funcID, size_t n, bool hasMV, double mv, const double *p1data, double *pdata)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto exprfunc = (double (*)(double)) funcEntry.func;
  if (hasMV)
  {
    if (funcEntry.name == "isMissval")
    {
      for (size_t i = 0; i < n; ++i) { pdata[i] = fp_is_equal(p1data[i], mv); }
    }
    else
    {
      for (size_t i = 0; i < n; ++i)
      {
        errno = -1;
        pdata[i] = fp_is_equal(p1data[i], mv) ? mv : exprfunc(p1data[i]);
        if (errno == EDOM || errno == ERANGE || std::isnan(pdata[i])) pdata[i] = mv;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      pdata[i] = exprfunc(p1data[i]);
      if (errno == EDOM || errno == ERANGE || std::isnan(pdata[i])) pdata[i] = mv;
    }
  }
}

void
func_expr_con_var(int funcID, bool hasMV, size_t n, double mv, double *odat, double cval, const double *idat)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto exprfunc = (double (*)(double, double)) funcEntry.func;

  if (hasMV)
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = (fp_is_equal(idat[i], mv)) ? mv : exprfunc(cval, idat[i]);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv;
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = exprfunc(cval, idat[i]);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv;
    }
  }
}

void
func_expr_var_con(int funcID, bool hasMV, size_t n, double mv, double *odat, const double *idat, double cval)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto exprfunc = (double (*)(double, double)) funcEntry.func;

  if (hasMV)
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = (fp_is_equal(idat[i], mv)) ? mv : exprfunc(idat[i], cval);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv;
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = exprfunc(idat[i], cval);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv;
    }
  }
}

void
func_expr_var_var(int funcID, bool hasMV, size_t n, double mv1, double mv2, double *odat, const double *idat1, const double *idat2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto exprfunc = (double (*)(double, double)) funcEntry.func;

  if (hasMV)
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = (fp_is_equal(idat1[i], mv1) || fp_is_equal(idat2[i], mv2)) ? mv1 : exprfunc(idat1[i], idat2[i]);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv1;
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      errno = -1;
      odat[i] = exprfunc(idat1[i], idat2[i]);
      if (errno == EDOM || errno == ERANGE || std::isnan(odat[i])) odat[i] = mv1;
    }
  }
}

static nodeType *
ex_fun_var(int init, int funcID, nodeType *p1)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto funcname = funcEntry.name;
  auto functype = funcEntry.type;
  auto funcflag = funcEntry.flag;

  auto gridID = p1->param.gridID;
  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto nlat = p1->param.nlat;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (functype == FT_CONST)
  {
    p->type = NodeEnum::typeCon;
    p->param.ngp = 0;
    p->param.nlev = 0;
    auto exprfunc = (double (*)(const ParamEntry *)) funcEntry.func;
    p->v = conNodeType(exprfunc(&p1->param));
  }
  else if (functype == FT_FLD)
  {
    p->param.gridID = pointID;
    p->param.ngp = 1;
  }
  else if (functype == FT_ZON)
  {
    if (zonalID == -1) cdo_abort("Function %s() is only available for regular 2D grids!", funcname);
    p->param.gridID = zonalID;
    p->param.ngp = nlat;
  }
  else if (functype == FT_VERT)
  {
    p->param.zaxisID = surfaceID;
    p->param.nlev = 1;
  }

  if (!init)
  {
    auto size = p->param.ngp * p->param.nlev;
    if (size == 0) size = 1;
    p->param.data = new double[size];
    double *pdata = p->param.data;
    double *p1data = p1->param.data;

    if (functype == FT_STD) { ex_fun_std(funcID, ngp * nlev, numMissVals, missval, p1data, pdata); }
    else if (functype == FT_FLD)
    {
      Field field;
      field.resize(ngp);
      if (funcflag == 1)
      {
        assert(p1->param.weight != nullptr);
        field.weightv.resize(ngp);
      }

      auto exprfunc = (double (*)(Field const &)) funcEntry.func;
      for (size_t k = 0; k < nlev; ++k)
      {
        fld_field_init(field, numMissVals, missval, ngp, p1data + k * ngp, p1->param.weight);
        pdata[k] = exprfunc(field);
      }
    }
    else if (functype == FT_ZON)
    {
      Field field1, field2;
      field1.resize(ngp);
      field2.resize(nlat);
      auto exprfunc = (void (*)(Field const &, Field &)) funcEntry.func;
      for (size_t k = 0; k < nlev; ++k)
      {
        fld_field_init(field1, numMissVals, missval, ngp, &p1data[k * ngp], nullptr);
        field1.grid = gridID;
        fld_field_init(field2, numMissVals, missval, nlat, &pdata[k * nlat], nullptr);
        exprfunc(field1, field2);
        array_copy(nlat, field2.vec_d.data(), &pdata[k * nlat]);
      }
    }
    else if (functype == FT_VERT)
    {
      Field field;
      field.resize(nlev);
      if (funcflag == 1) vert_weights(p1->param.zaxisID, nlev, field.weightv);
      auto exprfunc = (double (*)(Field const &)) funcEntry.func;
      for (size_t i = 0; i < ngp; ++i)
      {
        for (size_t k = 0; k < nlev; ++k) field.vec_d[k] = p1data[k * ngp + i];
        fld_field_init(field, numMissVals, missval, nlev, nullptr, nullptr);
        pdata[i] = exprfunc(field);
      }
    }
    else if (functype == FT_CONST) {}
    else
      cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcname);

    if (pdata) p->param.numMissVals = array_num_mv(p->param.ngp * p->param.nlev, pdata, missval);
  }

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static nodeType *
func_con_con(int funcID, nodeType *p1, nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto functype = funcEntry.type;
  if (functype != FT_STD) cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcEntry.name);

  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  auto exprfunc = (double (*)(double, double)) funcEntry.func;
  p->v = conNodeType(exprfunc(p1->con().value, p2->con().value));

  if (p1->isTmpObj) node_delete(p1);
  if (p2->isTmpObj) node_delete(p2);

  return p;
}

static nodeType *
func_con_var(int init, int funcID, const nodeType *p1, const nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto functype = funcEntry.type;
  if (functype != FT_STD) cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcEntry.name);

  auto ngp = (p2->param.ngp > 0) ? p2->param.ngp : 1;
  auto nlev = (p2->param.nlev > 0) ? p2->param.nlev : 1;
  auto hasMV = (p2->param.numMissVals > 0);
  auto missval2 = p2->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p2->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    auto n = ngp * nlev;
    p->param.data = new double[n];
    auto odat = p->param.data;
    const auto idat = p2->param.data;
    auto cval = p1->con().value;

    func_expr_con_var(funcID, hasMV, n, missval2, odat, cval, idat);

    p->param.numMissVals = array_num_mv(n, odat, missval2);
  }

  return p;
}

static nodeType *
func_var_con(int init, int funcID, const nodeType *p1, const nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto functype = funcEntry.type;
  if (functype != FT_STD) cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcEntry.name);

  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto hasMV = (p1->param.numMissVals > 0);
  auto missval1 = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    auto n = ngp * nlev;
    p->param.data = new double[n];
    auto odat = p->param.data;
    const auto idat = p1->param.data;
    auto cval = p2->con().value;

    func_expr_var_con(funcID, hasMV, n, missval1, odat, idat, cval);

    p->param.numMissVals = array_num_mv(n, odat, missval1);
  }

  return p;
}

static nodeType *
func_var_var(int init, int funcID, nodeType *p1, nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto functype = funcEntry.type;
  if (functype != FT_STD) cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcEntry.name);

  auto px = p1;
  auto numMissVals1 = p1->param.numMissVals;
  auto numMissVals2 = p2->param.numMissVals;
  auto missval1 = p1->param.missval;
  auto missval2 = p2->param.missval;

  auto ngp1 = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto ngp2 = (p2->param.ngp > 0) ? p2->param.ngp : 1;
  auto ngp = ngp1;

  auto nlat1 = p1->param.nlat;
  auto nlat2 = p2->param.nlat;
  auto isZonal = false;

  if (ngp1 != ngp2)
  {
    if (ngp1 == 1 || ngp2 == 1)
    {
      if (ngp1 == 1) { ngp = ngp2, px = p2; }
    }
    else if (nlat1 == nlat2 && ngp1 > ngp2) { isZonal = true; }
    else
    {
      cdo_abort("%s: Number of grid points differ (%s[%zu] <-> %s[%zu])", __func__, p1->param.name, ngp1, p2->param.name, ngp2);
    }
  }

  auto nlev1 = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto nlev2 = (p2->param.nlev > 0) ? p2->param.nlev : 1;

  auto nlev = nlev1;
  if (nlev1 != nlev2)
  {
    if (nlev1 == 1 || nlev2 == 1)
    {
      if (nlev1 == 1) { nlev = nlev2, px = p2; }
    }
    else { cdo_abort("%s: Number of levels differ (%s[%zu] <-> %s[%zu])", __func__, p1->param.name, nlev1, p2->param.name, nlev2); }
  }

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, px->param);

  if (p->param.steptype == TIME_CONSTANT)
  {
    auto steptype1 = p1->param.steptype;
    auto steptype2 = p2->param.steptype;
    if (steptype1 != TIME_CONSTANT) { p->param.steptype = steptype1; }
    else if (steptype2 != TIME_CONSTANT) { p->param.steptype = steptype2; }
  }

  p->param.name = tmpVarName;
  // printf("%s %s numMissVals %zu %zu\n", p->u.var.name.c_str(), px->param.name.c_str(), numMissVals1, numMissVals2);

  if (!init)
  {
    p->param.data = new double[ngp * nlev];

    for (size_t k = 0; k < nlev; ++k)
    {
      auto loff = k * ngp;
      size_t loff1 = (nlev1 > 1) ? k * ngp1 : 0;
      size_t loff2 = (nlev2 > 1) ? k * ngp2 : 0;

      const auto idat1 = p1->param.data + loff1;
      const auto idat2 = p2->param.data + loff2;
      auto odat = p->param.data + loff;
      auto hasMV = (numMissVals1 > 0 || numMissVals2 > 0);

      if (ngp1 != ngp2)
      {
        if (isZonal)
        {
          auto nlon = ngp1 / nlat1;
          for (size_t j = 0; j < nlat1; ++j)
            func_expr_var_con(funcID, hasMV, nlon, missval1, odat + j * nlon, idat1 + j * nlon, idat2[j]);
        }
        else
        {
          (ngp2 == 1) ? func_expr_var_con(funcID, hasMV, ngp, missval1, odat, idat1, idat2[0])
                      : func_expr_con_var(funcID, hasMV, ngp, missval2, odat, idat1[0], idat2);
        }
      }
      else { func_expr_var_var(funcID, hasMV, ngp, missval1, missval2, odat, idat1, idat2); }
    }

    p->param.numMissVals = array_num_mv(ngp * nlev, p->param.data, missval1);
  }

  return p;
}

static nodeType *
ex_fun(int init, int funcID, nodeType *p1)
{
  auto &funcEntry = funcSymbolTable[funcID];

  if (p1->type == NodeEnum::typeVar)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%s)", ExIn[init], funcEntry.name, p1->var().name);
    return ex_fun_var(init, funcID, p1);
  }
  else if (p1->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%g)", ExIn[init], funcEntry.name, p1->con().value);
    return ex_fun_con(funcID, p1);
  }
  else
    cdo_abort("Internal problem!");

  return nullptr;
}

static nodeType *
ex_fun2(int init, int funcID, nodeType *p1, nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];

  if (p1->type == NodeEnum::typeVar && p2->type == NodeEnum::typeVar)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%s, %s)", ExIn[init], funcEntry.name, p1->var().name, p2->var().name);
    return func_var_var(init, funcID, p1, p2);
  }
  else if (p1->type == NodeEnum::typeVar && p2->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%s, %g)", ExIn[init], funcEntry.name, p1->var().name, p2->con().value);
    return func_var_con(init, funcID, p1, p2);
  }
  else if (p1->type == NodeEnum::typeCon && p2->type == NodeEnum::typeVar)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%g, %s)", ExIn[init], funcEntry.name, p1->con().value, p2->var().name);
    return func_con_var(init, funcID, p1, p2);
  }
  else if (p1->type == NodeEnum::typeCon && p2->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tfunc\t%s (%g, %g)", ExIn[init], funcEntry.name, p1->con().value, p2->con().value);
    return func_con_con(funcID, p1, p2);
  }
  else
    cdo_abort("Internal problem!");

  return nullptr;
}

static nodeType *
ex_remap(int init, int funcID, nodeType *p1, nodeType *p2)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto funcname = funcEntry.name;
  auto functype = funcEntry.type;

  if (functype != FT_REMAP) cdo_abort("Intermal error, wrong function type (%d) for %s()!", functype, funcname);

  auto gridID = p1->param.gridID;
  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (p1->type == NodeEnum::typeVar && p2->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tremap\t%s (%s, %g)", ExIn[init], funcEntry.name, p1->var().name, p2->con().value);

    auto gridID2 = p2->con().value;
    auto ngp2 = gridInqSize(gridID2);
    p->param.gridID = gridID2;
    p->param.ngp = ngp2;
    if (!init)
    {
      auto size = ngp2 * p->param.nlev;
      p->param.data = size ? new double[size] : nullptr;
      double *pdata = p->param.data;
      double *p1data = p1->param.data;

      Field field1, field2;
      field1.resize(ngp);
      field2.resize(ngp2);
      auto exprfunc = (void (*)(Field const &, Field &)) funcEntry.func;
      for (size_t k = 0; k < nlev; ++k)
      {
        fld_field_init(field1, numMissVals, missval, ngp, &p1data[k * ngp], nullptr);
        field1.grid = gridID;
        field2.grid = gridID2;
        fld_field_init(field2, numMissVals, missval, ngp2, &pdata[k * ngp2], nullptr);
        exprfunc(field1, field2);
        array_copy(ngp2, field2.vec_d.data(), &pdata[k * ngp2]);
      }
    }
  }
  else { cdo_abort("Syntax error in call to %s(p, c), check type of parameter!", funcname); }

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static size_t
get_levidx(size_t nlev, const double *data, double value, std::string const &funcname)
{
  size_t levidx;

  for (levidx = 0; levidx < nlev; ++levidx)
    if (is_equal(data[levidx], value)) break;
  if (levidx == nlev) cdo_abort("%s(): level %g not found!", funcname, value);

  return levidx;
}

static void
get_levidxrange(size_t nlev, const double *data, double value1, double value2, std::string const &funcname, size_t &levidx1,
                size_t &levidx2)
{
  long n = nlev;
  if (data[0] <= data[nlev - 1])
  {
    long i;
    for (i = 0; i < n; ++i)
      if (data[i] >= value1) break;
    if (i == n) cdo_abort("%s(): lower level %g not found!", funcname, value1);
    levidx1 = i;

    for (i = n - 1; i >= 0; --i)
      if (data[i] <= value2) break;
    if (i < 0) cdo_abort("%s(): upper level %g not found!", funcname, value2);
    if (i < (long) levidx1) cdo_abort("%s(): level range %g to %g not found!", funcname, value1, value2);
    levidx2 = i;
  }
  else
  {
    long i;
    for (i = 0; i < n; ++i)
      if (data[i] <= value2) break;
    if (i == n) cdo_abort("%s(): upper level %g not found!", funcname, value1);
    levidx1 = i;

    for (i = n - 1; i >= 0; --i)
      if (data[i] >= value1) break;
    if (i < 0) cdo_abort("%s(): lower level %g not found!", funcname, value2);
    if (i < (long) levidx1) cdo_abort("%s(): level range %g to %g not found!", funcname, value1, value2);
    levidx2 = i;
  }
}

static nodeType *
ex_fun1c(int init, int funcID, nodeType *p1, nodeType *p2, ParseParamType &parseArg)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto funcname = funcEntry.name;
  if (p1->type != NodeEnum::typeVar) cdo_abort("1st parameter of function %s() needs to be a variable!", funcname);
  if (p1->isTmpObj) cdo_abort("Temporary objects not allowed in function %s()!", funcname);
  if (p2->type != NodeEnum::typeCon) cdo_abort("2nd parameter of function %s() needs to be a constant!", funcname);

  auto value = p2->con().value;

  if (parseArg.debug)
    cdo_print("\t%s\tfunc\t%s=%s(%s[N%zu][L%zu], %g)", ExIn[init], tmpVarName, funcname, p1->param.name, p1->param.ngp,
              p1->param.nlev, value);

  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (init)
  {
    if (p1->param.longname.size()) p->param.longname = p1->param.longname;
    if (p1->param.units.size()) p->param.units = p1->param.units;
  }

  p->param.nlev = 1;

  auto zaxisID = p1->param.zaxisID;
  auto coordID = params_get_coord_ID(parseArg, 'z', zaxisID);

  std::vector<double> data;
  double *pdata = nullptr;

  if (init)
  {
    parseArg.coords[coordID].needed = true;

    data.resize(nlev);
    pdata = data.data();
    cdo_zaxis_inq_levels(zaxisID, pdata);
  }
  else { pdata = parseArg.coords[coordID].data.data(); }

  size_t levidx = 0;
  if (funcname == "sellevidx")
  {
    auto ilevidx = std::lround(value);
    if (ilevidx < 1 || ilevidx > (long) nlev)
      cdo_abort("%s(): level index %ld out of range (range: 1-%zu)!", funcname, ilevidx, nlev);
    levidx = (size_t) ilevidx - 1;
  }
  else if (funcname == "sellevel") { levidx = get_levidx(nlev, pdata, value, funcname); }
  else
    cdo_abort("Function %s() not implemented!", funcname);

  if (init)
  {
    auto level = data[levidx];
    auto zaxisID2 = zaxisCreate(zaxisInqType(zaxisID), 1);
    zaxisDefLevels(zaxisID2, &level);
    p->param.zaxisID = zaxisID2;
  }

  if (!init)
  {
    p->param.data = new double[ngp];
    pdata = p->param.data;
    const auto p1data = p1->param.data + ngp * levidx;
    array_copy(ngp, p1data, pdata);
    if (numMissVals) numMissVals = array_num_mv(ngp, pdata, missval);
    p->param.numMissVals = numMissVals;
  }

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static nodeType *
ex_fun2c(int init, int funcID, nodeType *p1, nodeType *p2, nodeType *p3, ParseParamType &parseArg)
{
  auto &funcEntry = funcSymbolTable[funcID];
  auto funcname = funcEntry.name;
  if (p1->type != NodeEnum::typeVar) cdo_abort("Parameter of function %s() needs to be a variable!", funcname);
  if (p1->isTmpObj) cdo_abort("Temporary objects not allowed in function %s()!", funcname);
  if (p2->type != NodeEnum::typeCon) cdo_abort("2nd parameter of function %s() needs to be a constant!", funcname);
  if (p3->type != NodeEnum::typeCon) cdo_abort("3rd parameter of function %s() needs to be a constant!", funcname);
  auto value1 = p2->con().value;
  auto value2 = p3->con().value;

  if (parseArg.debug)
    cdo_print("\t%s\tfunc\t%s=%s(%s[N%zu][L%zu], %g, %g)", ExIn[init], tmpVarName, funcname, p1->param.name, p1->param.ngp,
              p1->param.nlev, value1, value2);

  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (init)
  {
    if (p1->param.longname.size()) p->param.longname = p1->param.longname;
    if (p1->param.units.size()) p->param.units = p1->param.units;
  }

  auto zaxisID = p1->param.zaxisID;
  auto coordID = params_get_coord_ID(parseArg, 'z', zaxisID);

  std::vector<double> data;
  double *pdata = nullptr;

  if (init)
  {
    parseArg.coords[coordID].needed = true;

    data.resize(nlev);
    pdata = data.data();
    cdo_zaxis_inq_levels(zaxisID, pdata);
  }
  else { pdata = parseArg.coords[coordID].data.data(); }

  size_t levidx1 = 0;
  size_t levidx2 = 0;
  if (funcname == "sellevidxrange")
  {
    if (value2 < value1) cdo_abort("%s(): first level index is greater than last level index!", funcname);
    auto ilevidx1 = std::lround(value1);
    if (ilevidx1 < 1 || ilevidx1 > (long) nlev)
      cdo_abort("%s(): level index1 %ld out of range (range: 1-%zu)!", funcname, ilevidx1, nlev);
    levidx1 = (size_t) ilevidx1 - 1;
    auto ilevidx2 = std::lround(value2);
    if (ilevidx2 < 1 || ilevidx2 > (long) nlev)
      cdo_abort("%s(): level index2 %ld out of range (range: 1-%zu)!", funcname, ilevidx2, nlev);
    levidx2 = (size_t) ilevidx2 - 1;
  }
  else if (funcname == "sellevelrange")
  {
    if (value2 < value1) cdo_abort("%s(): first level is greater than last level!", funcname);
    get_levidxrange(nlev, pdata, value1, value2, funcname, levidx1, levidx2);
  }
  else
    cdo_abort("Function %s() not implemented!", funcname);

  int nlevout = levidx2 - levidx1 + 1;
  p->param.nlev = nlevout;

  if (init)
  {
    auto zaxisID2 = zaxisCreate(zaxisInqType(zaxisID), nlevout);
    zaxisDefLevels(zaxisID2, &pdata[levidx1]);
    if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
    {
      std::vector<double> bounds(nlev);
      zaxisInqLbounds(zaxisID, bounds.data());
      zaxisDefLbounds(zaxisID2, &bounds[levidx1]);
      zaxisInqUbounds(zaxisID, bounds.data());
      zaxisDefUbounds(zaxisID2, &bounds[levidx1]);
    }
    p->param.zaxisID = zaxisID2;
  }

  if (!init)
  {
    p->param.data = new double[ngp * nlevout];
    auto paramdata = p->param.data;
    const auto p1data = p1->param.data + ngp * levidx1;
    array_copy(ngp * nlevout, p1data, paramdata);
    if (numMissVals) numMissVals = array_num_mv(ngp * nlevout, paramdata, missval);
    p->param.numMissVals = numMissVals;
  }

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static nodeType *
coord_fun(int init, int funcID, nodeType *p1, ParseParamType &parseArg)
{
  std::string funcname = funcSymbolTable[funcID].name;
  if (p1->type != NodeEnum::typeVar) cdo_abort("Parameter of function %s() needs to be a variable!", funcname);
  if (p1->isTmpObj) cdo_abort("Temporary objects not allowed in function %s()!", funcname);

  auto varName = p1->var().name;

  // clang-format off
  if      (funcname == "clon")       varName += ".x";
  else if (funcname == "clat")       varName += ".y";
  else if (funcname == "clev")       varName += ".z";
  else if (funcname == "clevidx")    varName += ".i";
  else if (funcname == "cthickness") varName += ".d";
  else if (funcname == "gridarea")   varName += ".a";
  else if (funcname == "gridweight") varName += ".w";
  else if (funcname == "gridindex")  varName += ".g";
  else cdo_abort("Implementation missing for function %s!", funcname);
  // clang-format on

  p1->var().name = varName;

  auto p = expr_run(p1, parseArg);
  p->param.hasMV = false;

  if (!init)
  {
    /*
    size_t ngp  = p1->param.ngp;
    size_t nlev = p1->param.nlev;
    p->param.data = new double[ngp * nlev];
    double *pdata  = p->param.data;
    double *p1data = p1->param.data;

    for (size_t i = 0; i < ngp*nlev; ++i) pdata[i] = p1data[i];
    */
  }

  return p;
}

static nodeType *
ex_uminus_var(int init, nodeType *p1)
{
  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    p->param.data = new double[ngp * nlev];
    double *pdata = p->param.data;
    const double *p1data = p1->param.data;

    if (numMissVals)
    {
      for (size_t i = 0; i < ngp * nlev; ++i) pdata[i] = fp_is_equal(p1data[i], missval) ? missval : -(p1data[i]);
    }
    else
    {
      for (size_t i = 0; i < ngp * nlev; ++i) pdata[i] = -(p1data[i]);
    }

    p->param.numMissVals = numMissVals;
  }

  return p;
}

static nodeType *
ex_uminus_con(nodeType *p1)
{
  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  p->v = conNodeType(-(p1->con().value));

  return p;
}

static nodeType *
ex_uminus(int init, nodeType *p1)
{
  nodeType *p = nullptr;

  if (p1->type == NodeEnum::typeVar)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tneg\t- (%s)", ExIn[init], p1->var().name);
    p = ex_uminus_var(init, p1);
  }
  else if (p1->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tneg\t- (%g)", ExIn[init], p1->con().value);
    p = ex_uminus_con(p1);
  }
  else
    cdo_abort("Internal problem!");

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static nodeType *
ex_not_var(int init, nodeType *p1)
{
  auto ngp = (p1->param.ngp > 0) ? p1->param.ngp : 1;
  auto nlev = (p1->param.nlev > 0) ? p1->param.nlev : 1;
  auto numMissVals = p1->param.numMissVals;
  auto missval = p1->param.missval;

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, p1->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    p->param.data = new double[ngp * nlev];
    double *pdata = p->param.data;
    const double *p1data = p1->param.data;

    if (numMissVals)
    {
      for (size_t i = 0; i < ngp * nlev; ++i) pdata[i] = fp_is_equal(p1data[i], missval) ? missval : unary_op_not(p1data[i]);
    }
    else
    {
      for (size_t i = 0; i < ngp * nlev; ++i) pdata[i] = unary_op_not(p1data[i]);
    }

    p->param.numMissVals = numMissVals;
  }

  return p;
}

static nodeType *
ex_not_con(const nodeType *p1)
{
  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  p->v = conNodeType(unary_op_not(p1->con().value));

  return p;
}

static nodeType *
ex_not(int init, nodeType *p1)
{
  nodeType *p = nullptr;

  if (p1->type == NodeEnum::typeVar)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tnot\t! (%s)", ExIn[init], p1->var().name);
    p = ex_not_var(init, p1);
  }
  else if (p1->type == NodeEnum::typeCon)
  {
    if (Options::cdoVerbose) cdo_print("\t%s\tnot\t! (%g)", ExIn[init], p1->con().value);
    p = ex_not_con(p1);
  }
  else
    cdo_abort("Internal problem!");

  if (p1->isTmpObj) node_delete(p1);

  return p;
}

static void
str_add_node_info(char *string, size_t stringlen, nodeType *p, const char *ext)
{
  auto len = std::strlen(string);
  if (p->type == NodeEnum::typeCon) { snprintf(string + len, stringlen - len, "%g%s", p->con().value, ext); }
  else { snprintf(string + len, stringlen - len, "%s[N%zu][L%zu]%s", p->var().name.c_str(), p->param.ngp, p->param.nlev, ext); }
}

static nodeType *
ex_ifelse(int init, nodeType *p1, nodeType *p2, nodeType *p3)
{
  if (Options::cdoVerbose)
  {
    char strbuffer[1024];
    std::snprintf(strbuffer, sizeof(strbuffer), "\t%s\tifelse\t", ExIn[init]);
    str_add_node_info(strbuffer, sizeof(strbuffer), p1, " ? ");
    str_add_node_info(strbuffer, sizeof(strbuffer), p2, " : ");
    str_add_node_info(strbuffer, sizeof(strbuffer), p3, "");
    cdo_print(strbuffer);
  }

  nodeType *p0 = nullptr;
  nodeType *pnodes[3] = { p1, p2, p3 };
  for (unsigned i = 0; i < 3; ++i)
    if (pnodes[i]->type != NodeEnum::typeCon)
    {
      p0 = pnodes[i];
      break;
    }

  if (p0 == nullptr) cdo_abort("expr?expr:expr: no data variable found!");

  auto missval = p0->param.missval;
  auto ngp = (p0->param.ngp > 0) ? p0->param.ngp : 1;
  auto nlev = (p0->param.nlev > 0) ? p0->param.nlev : 1;
  auto px = p0;

  size_t numMissVals1 = 0;
  auto missval1 = missval;
  const double *pdata1 = nullptr;
  size_t ngp1 = 1;
  size_t nlev1 = 1;

  if (p1->type == NodeEnum::typeCon) { pdata1 = &p1->con().value; }
  else
  {
    numMissVals1 = p1->param.numMissVals;
    ngp1 = (p1->param.ngp > 0) ? p1->param.ngp : 1;
    nlev1 = (p1->param.nlev > 0) ? p1->param.nlev : 1;
    missval1 = p1->param.missval;
    pdata1 = p1->param.data;

    if (ngp1 > 1 && ngp1 != ngp)
    {
      if (ngp != 1) cdo_abort("expr?expr:expr: Number of grid points differ (ngp = %zu, ngp1 = %zu)", ngp, ngp1);

      ngp = ngp1;
      px = p1;
    }

    if (nlev1 > 1 && nlev1 != nlev)
    {
      if (nlev != 1) cdo_abort("expr?expr:expr: Number of levels differ (nlev = %zu, nlev1 = %zu)", nlev, nlev1);

      nlev = nlev1;
      px = p1;
    }
  }

  auto missval2 = missval1;
  const double *pdata2 = nullptr;
  size_t ngp2 = 1;
  size_t nlev2 = 1;

  if (p2->type == NodeEnum::typeCon) { pdata2 = &p2->con().value; }
  else
  {
    ngp2 = (p2->param.ngp > 0) ? p2->param.ngp : 1;
    nlev2 = (p2->param.nlev > 0) ? p2->param.nlev : 1;
    missval2 = p2->param.missval;
    pdata2 = p2->param.data;

    if (ngp2 > 1 && ngp2 != ngp)
    {
      if (ngp != 1) cdo_abort("expr?expr:expr: Number of grid points differ (ngp = %zu, ngp2 = %zu)", ngp, ngp2);

      ngp = ngp2;
      px = p2;
    }

    if (nlev2 > 1 && nlev2 != nlev)
    {
      if (nlev != 1) cdo_abort("expr?expr:expr: Number of levels differ (nlev = %zu, nlev2 = %zu)", nlev, nlev2);

      nlev = nlev2;
      px = p2;
    }
  }

  auto missval3 = missval1;
  const double *pdata3 = nullptr;
  size_t ngp3 = 1;
  size_t nlev3 = 1;

  if (p3->type == NodeEnum::typeCon) { pdata3 = &p3->con().value; }
  else
  {
    ngp3 = (p3->param.ngp > 0) ? p3->param.ngp : 1;
    nlev3 = (p3->param.nlev > 0) ? p3->param.nlev : 1;
    missval3 = p3->param.missval;
    pdata3 = p3->param.data;

    if (ngp3 > 1 && ngp3 != ngp)
    {
      if (ngp != 1) cdo_abort("expr?expr:expr: Number of grid points differ (ngp = %zu, ngp3 = %zu)", ngp, ngp3);

      ngp = ngp3;
      px = p3;
    }

    if (nlev3 > 1 && nlev3 != nlev)
    {
      if (nlev != 1) cdo_abort("expr?expr:expr: Number of levels differ (nlev = %zu, nlev3 = %zu)", nlev, nlev3);

      nlev = nlev3;
      px = p3;
    }
  }

  auto p = new nodeType;

  p->type = NodeEnum::typeVar;
  p->isTmpObj = true;
  p->v = varNodeType(tmpVarName);
  param_meta_copy(p->param, px->param);
  p->param.name = tmpVarName;

  if (!init)
  {
    size_t numMissVals = 0;

    p->param.data = new double[ngp * nlev];

    for (size_t k = 0; k < nlev; ++k)
    {
      size_t loff1 = (nlev1 == 1) ? 0 : k * ngp1;
      size_t loff = k * ngp;
      size_t loff2 = (nlev2 == 1) ? 0 : loff;
      size_t loff3 = (nlev3 == 1) ? 0 : loff;

      const auto idat1 = pdata1 + loff1;
      const auto idat2 = pdata2 + loff2;
      const auto idat3 = pdata3 + loff3;
      auto odat = p->param.data + loff;

#ifdef _OPENMP
#pragma omp parallel for if (ngp > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < ngp; ++i)
      {
        auto ival1 = idat1[(ngp1 > 1) ? i : 0];
        auto ival2 = idat2[(ngp2 > 1) ? i : 0];
        auto ival3 = idat3[(ngp3 > 1) ? i : 0];

        if (numMissVals1 && fp_is_equal(ival1, missval1))
          odat[i] = missval1;
        else if (is_not_equal(ival1, 0.0))
          odat[i] = fp_is_equal(ival2, missval2) ? missval1 : ival2;
        else
          odat[i] = fp_is_equal(ival3, missval3) ? missval1 : ival3;
      }

      numMissVals += array_num_mv(ngp, odat, missval1);
    }

    p->param.numMissVals = numMissVals;
  }

  if (p1->isTmpObj) node_delete(p1);
  if (p2->isTmpObj) node_delete(p2);
  if (p3->isTmpObj) node_delete(p3);

  return p;
}
/*
static int
exNode(nodeType *p, ParseParamType &parseArg)
{
  if (!p) return 0;

  // node is leaf
  if (p->type == NodeEnum::typeCon || p->type == NodeEnum::typeVar || p->u.opr.nops == 0) return 0;

  // node has children
  for (int k = 0; k < p->u.opr.nops; ++k) exNode(p->u.opr.op[k], parseArg);

  return 0;
}
*/

static int
param_search_name(int nparam, std::vector<ParamEntry> const &params, std::string const &name)
{
  for (int varID = 0; varID < nparam; ++varID)
  {
    if (params[varID].isValid && params[varID].name == name) return varID;
  }

  return -1;
}

static int
param_search_name_size(int nparam, std::vector<ParamEntry> &params, std::string const &name, size_t ngp, size_t nlev)
{
  for (int varID = 0; varID < nparam; ++varID)
  {
    if (params[varID].name == name)
    {
      if ((ngp == params[varID].ngp && nlev == params[varID].nlev) || (ngp == 0 && nlev == 0)) { return varID; }

      params[varID].isValid = false;
    }
  }

  return -1;
}

static void
param_print(std::string const &vname, const ParamEntry &param, long tsID)
{
  constexpr size_t maxout = 100;
  const auto data = param.data;
  for (size_t k = 0; k < param.nlev; ++k)
    for (size_t i = 0; i < param.ngp; ++i)
    {
      if (i < maxout || i >= (param.ngp - maxout))
      {
        auto v = data[k * param.ngp + i];
        if (param.steptype == TIME_CONSTANT)
          fprintf(stdout, "   %s[lev=%zu:gp=%zu] = %g\n", vname.c_str(), k + 1, i + 1, v);
        else
          fprintf(stdout, "   %s[ts=%ld:lev=%zu:gp=%zu] = %g\n", vname.c_str(), tsID, k + 1, i + 1, v);
      }
      else if (i == maxout) { fprintf(stdout, "   .......\n"); }
    }
}

static void
add_new_constant(std::string const &varname, ParseParamType &parseArg, std::vector<ParamEntry> &params, const ParamEntry &param)
{
  auto varID = parseArg.numParams;
  if (varID >= parseArg.maxParams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxParams);

  param_meta_copy(params[varID], param);
  params[varID].type = ParamType::CONST;
  params[varID].isValid = true;
  params[varID].ngp = 0;
  params[varID].nlat = 0;
  params[varID].nlev = 0;
  params[varID].missval = -9.e33;
  params[varID].numMissVals = 0;
  params[varID].name = varname;
  parseArg.numParams++;
}

static void
add_new_param(std::string const &varname, ParseParamType &parseArg, std::vector<ParamEntry> &params, const ParamEntry &param)
{
  auto varID = parseArg.numParams;
  if (varID >= parseArg.maxParams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxParams);

  param_meta_copy(params[varID], param);
  params[varID].isValid = true;
  params[varID].hasMV = param.hasMV;
  params[varID].name = varname;
  params[varID].numMissVals = param.numMissVals;
  if (param.units.size()) params[varID].units = param.units;
  if (param.longname.size()) params[varID].longname = param.longname;
  parseArg.numParams++;
}

static nodeType *
expr_run_cmd(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;
  auto &params = parseArg.params;

  auto const &cmdName = p->cmd().cmdName;
  auto const &varName = p->cmd().varName;
  if (parseArg.debug) cdo_print("\tstatement\t\t%s(%s)", cmdName, varName);

  auto varID = param_search_name(parseArg.numParams, params, varName);
  if (varID == -1) cdo_abort("Variable %s not found, needed for statement %s(%s)!", varName, cmdName, varName);

  if (init)
  {
    if (cmdName == "remove") params[varID].remove = true;
  }
  else
  {
    if (cmdName == "print") param_print(varName, params[varID], std::lround(params[parseArg.tsID].data[0]));
  }

  return nullptr;
}

static nodeType *
expr_run_con(nodeType *p, const ParseParamType &parseArg)
{
  if (parseArg.debug) cdo_print("\tpush\tconst\t%g", p->con().value);
  return p;
}

static int
expr_run_var_grid(std::string const &vnm, int coord, ParseParamType &parseArg)
{
  auto &params = parseArg.params;

  auto varname = vnm.substr(0, vnm.size() - 2);
  auto varID = param_search_name(parseArg.numParams, params, varname);
  if (varID == -1) cdo_abort("Coordinate %c: variable >%s< not found!", coord, varname);

  auto nvarID = parseArg.numParams;
  if (nvarID >= parseArg.maxParams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxParams);

  auto coordID = params_get_coord_ID(parseArg, coord, params[varID].gridID);
  parseArg.coords[coordID].needed = true;
  auto const &units = parseArg.coords[coordID].units;
  auto const &longname = parseArg.coords[coordID].longname;

  params[nvarID].isValid = true;
  params[nvarID].coord = coord;
  params[nvarID].hasMV = false;
  params[nvarID].name = vnm;
  params[nvarID].missval = params[varID].missval;
  params[nvarID].gridID = params[varID].gridID;
  params[nvarID].zaxisID = parseArg.surfaceID;
  params[nvarID].steptype = TIME_CONSTANT;
  params[nvarID].ngp = params[varID].ngp;
  params[nvarID].nlat = params[varID].nlat;
  params[nvarID].nlev = 1;
  if (units.size()) params[nvarID].units = units;
  if (longname.size()) params[nvarID].longname = longname;
  parseArg.numParams++;

  return nvarID;
}

static int
expr_run_var_zaxis(std::string const &vnm, int coord, ParseParamType &parseArg)
{
  auto &params = parseArg.params;

  auto varname = vnm.substr(0, vnm.size() - 2);
  auto varID = param_search_name(parseArg.numParams, params, varname);
  if (varID == -1) cdo_abort("Coordinate %c: variable >%s< not found!", coord, varname);

  auto nvarID = parseArg.numParams;
  if (nvarID >= parseArg.maxParams) cdo_abort("Too many parameter (limit=%d)", parseArg.maxParams);

  auto coordID = params_get_coord_ID(parseArg, coord, params[varID].zaxisID);
  parseArg.coords[coordID].needed = true;
  auto const &units = parseArg.coords[coordID].units;
  auto const &longname = parseArg.coords[coordID].longname;

  params[nvarID].isValid = true;
  params[nvarID].coord = coord;
  params[nvarID].hasMV = false;
  params[nvarID].name = vnm;
  params[nvarID].missval = params[varID].missval;
  params[nvarID].gridID = parseArg.pointID;
  params[nvarID].zaxisID = params[varID].zaxisID;
  params[nvarID].steptype = TIME_CONSTANT;
  params[nvarID].ngp = 1;
  params[nvarID].nlev = params[varID].nlev;
  if (units.size()) params[nvarID].units = units;
  if (longname.size()) params[nvarID].longname = longname;
  parseArg.numParams++;

  return nvarID;
}

static nodeType *
expr_run_var(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;
  auto const &params = parseArg.params;

  auto const &vnm = p->var().name;
  auto varID = param_search_name(parseArg.numParams, params, vnm);
  if (varID == -1 && init)
  {
    int len = vnm.size();
    if (len > 2 && vnm[len - 2] == '.')
    {
      auto coord = vnm[len - 1];
      if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w' || coord == 'g')
      {
        varID = expr_run_var_grid(vnm, coord, parseArg);
      }
      else if (coord == 'z' || coord == 'i' || coord == 'd') { varID = expr_run_var_zaxis(vnm, coord, parseArg); }
    }
  }

  if (varID == -1) { cdo_abort("Variable >%s< not found!", vnm); }
  else if (init)
  {
    if (varID < parseArg.numVars1 && !parseArg.needed[varID]) parseArg.needed[varID] = true;
  }

  param_meta_copy(p->param, params[varID]);
  p->param.coord = params[varID].coord;
  p->param.hasMV = params[varID].hasMV;
  p->param.name = params[varID].name;
  p->param.longname = params[varID].longname;
  p->param.units = params[varID].units;
  p->isTmpObj = false;

  if (!init)
  {
    p->param.data = params[varID].data;
    p->param.numMissVals = params[varID].numMissVals;
  }

  if (parseArg.debug)
  {
    if (p->param.ngp == 0 && p->param.nlev == 0)
      cdo_print("\tpush\tvar\t%s", vnm);
    else
      cdo_print("\tpush\tvar\t%s[N%zu][L%zu]", vnm, p->param.ngp, p->param.nlev);
  }

  return p;
}

static nodeType *
expr_run_fun2(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;

  auto funcID = get_funcID(p->fun().name);
  auto const &funcEntry = funcSymbolTable[funcID];

  auto fnode1 = expr_run(p->fun().op[0], parseArg);
  auto fnode2 = expr_run(p->fun().op[1], parseArg);

  if (funcEntry.type == FT_1C) { return ex_fun1c(init, funcID, fnode1, fnode2, parseArg); }
  else if (funcEntry.type == FT_STD) { return ex_fun2(init, funcID, fnode1, fnode2); }
  else if (funcEntry.type == FT_REMAP) { return ex_remap(init, funcID, fnode1, fnode2); }
  else { cdo_abort("Syntax error in call to %s(p1, p2), check number of parameter!", p->fun().name); }

  return nullptr;
}

static nodeType *
expr_run_fun3(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;

  auto funcID = get_funcID(p->fun().name);
  auto const &funcEntry = funcSymbolTable[funcID];

  auto fnode1 = expr_run(p->fun().op[0], parseArg);
  auto fnode2 = expr_run(p->fun().op[1], parseArg);
  auto fnode3 = expr_run(p->fun().op[2], parseArg);

  if (funcEntry.type == FT_2C) { return ex_fun2c(init, funcID, fnode1, fnode2, fnode3, parseArg); }
  else { cdo_abort("Syntax error in call to %s(p1, p2, p3), check number of parameter!", p->fun().name); }

  return nullptr;
}

static nodeType *
expr_run_fun(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;
  auto const &params = parseArg.params;

  auto funcID = get_funcID(p->fun().name);
  auto &funcEntry = funcSymbolTable[funcID];

  auto numFunArgs = (p->fun().op[0]) ? p->fun().nops : 0;

  if (numFunArgs == 3) { return expr_run_fun3(p, parseArg); }
  else if (numFunArgs == 2) { return expr_run_fun2(p, parseArg); }
  else if (numFunArgs == 1)
  {
    auto fnode = expr_run(p->fun().op[0], parseArg);

    if (funcEntry.type == FT_COORD) { return coord_fun(init, funcID, fnode, parseArg); }
    else if (funcEntry.type == FT_STD || funcEntry.type == FT_FLD || funcEntry.type == FT_ZON || funcEntry.type == FT_VERT
             || funcEntry.type == FT_CONST)
    {
      if (funcEntry.flag == 1)
      {
        auto coordID = params_get_coord_ID(parseArg, 'w', fnode->param.gridID);
        if (init)
          parseArg.coords[coordID].needed = true;
        else
          fnode->param.weight = parseArg.coords[coordID].data.data();
      }

      return ex_fun(init, funcID, fnode);
    }
    else { cdo_abort("Syntax error in call to %s(p), check number of parameter!", p->fun().name); }
  }
  else if (numFunArgs == 0)
  {
    if (funcEntry.type == FT_0)
    {
      auto vartsID = parseArg.tsID;
      return ex_fun0_con(init, funcID, params[vartsID].data);
    }
    else { cdo_abort("Syntax error in call to %s(), check number of parameter!", p->fun().name); }
  }

  return nullptr;
}

static nodeType *
expr_run_opr(nodeType *p, ParseParamType &parseArg)
{
  auto init = parseArg.init;
  auto &params = parseArg.params;

  // clang-format off
  switch (p->opr().oper)
    {
    case '=':
      {
        auto rnode = expr_run(p->opr().op[1], parseArg);

        auto const &varname2 = p->opr().op[0]->var().name;

        auto isVar = (rnode && rnode->type == NodeEnum::typeVar && rnode->param.ngp > 0 && rnode->param.nlev > 0);
        if (parseArg.debug)
          {
            if (isVar)
              cdo_print("\tpop\tvar\t%s[N%zu][L%zu]", varname2, rnode->param.ngp, rnode->param.nlev);
            else
              cdo_print("\tpop\tconst\t%s", varname2);
          }

        // auto varID = param_search_name(parseArg.nparams, params, varname2);
        auto ngp = isVar ? rnode->param.ngp : 0;
        auto nlev = isVar ? rnode->param.nlev : 0;
        auto varID = param_search_name_size(parseArg.numParams, params, varname2, ngp, nlev);
        if (init)
          {
            if (varID >= 0)
              {
                // printf("  found %s\n", varname2);
                if (varID < parseArg.numVars1)
                  {
                    if (rnode->param.nlev > 0 && params[varID].nlev != rnode->param.nlev)
                      cdo_abort("The number of layers must not change (name=%s layers: in=%zu out=%zu)!", params[varID].name,
                               params[varID].nlev, rnode->param.nlev);

                    params[varID].select = true;
                    parseArg.needed[varID] = true;
                  }
                else if (params[varID].coord)
                  cdo_abort("Coordinate variable %s is read only!", varname2);
                /*
                  else
                  cdo_warning("Variable %s already defined!", varname2);
                */
              }
            else if (rnode && rnode->type == NodeEnum::typeCon)
              {
                add_new_constant(varname2, parseArg, params, rnode->param);
              }
            else if (p->opr().op[1]->type != NodeEnum::typeCon)
              {
                add_new_param(varname2, parseArg, params, rnode->param);
              }
          }
        else
          {
            if (varID < 0) cdo_abort("Variable >%s< not found!", varname2);

            if (params[varID].coord) cdo_abort("Coordinate variable %s is read only!", varname2);
            param_meta_copy(p->param, params[varID]);
            p->param.name = params[varID].name;
            p->param.data = params[varID].data;
            p->isTmpObj = false;

            ex_copy(init, p, rnode);
            params[varID].numMissVals = p->param.numMissVals;
          }

        if (rnode && rnode->isTmpObj)
          {
            node_delete(rnode);
            rnode = nullptr;
          }
        // else delete rnode;
        break;
      }
    case UMINUS: return ex_uminus(init, expr_run(p->opr().op[0], parseArg));
    case NOT:    return ex_not(init, expr_run(p->opr().op[0], parseArg));
    case '?':    return ex_ifelse(init, expr_run(p->opr().op[0], parseArg), expr_run(p->opr().op[1], parseArg),
                                  expr_run(p->opr().op[2], parseArg));
    default:     return expr(init, p->opr().oper, expr_run(p->opr().op[0], parseArg), expr_run(p->opr().op[1], parseArg));
    }
  // clang-format on

  return nullptr;
}

nodeType *
expr_run(nodeType *p, ParseParamType &parseArg)
{
  pointID = parseArg.pointID;
  zonalID = parseArg.zonalID;
  surfaceID = parseArg.surfaceID;

  if (!p) return nullptr;

  switch (p->type)
  {
    case NodeEnum::typeCmd: return expr_run_cmd(p, parseArg);
    case NodeEnum::typeCon: return expr_run_con(p, parseArg);
    case NodeEnum::typeVar: return expr_run_var(p, parseArg);
    case NodeEnum::typeFun: return expr_run_fun(p, parseArg);
    case NodeEnum::typeOpr: return expr_run_opr(p, parseArg);
    case NodeEnum::typeUndef: return nullptr;
  }

  return nullptr;
}

int
params_get_coord_ID(const ParseParamType &parseArg, int coord, int cdiID)
{
  auto ncoords = parseArg.numCoords;
  for (int coordID = 0; coordID < ncoords; ++coordID)
  {
    if (parseArg.coords[coordID].coord == coord && parseArg.coords[coordID].cdiID == cdiID) return coordID;
  }

  cdo_abort("%s: coordinate %c not found!", __func__, coord);

  return -1;
}
