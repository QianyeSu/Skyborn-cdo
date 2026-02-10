#include <algorithm>

#include "cdi.h"

#define streamDefField cdo_def_field
#define streamWriteField cdo_write_field

#include "afterburner.h"
#include "constants.h"
#include "compare.h"
#include "vertical_interp.h"

int afterDebug = 0;

static char *
FieldName(int code, const char *text)
{
  static char name[256];
  std::snprintf(name, sizeof(name), "[%3d].%s", code, text);
  return name;
}

// Free array space
template <typename T>
static void
FreeMemory(T &ptr)
{
  if (ptr)
  {
    std::free(ptr);
    ptr = nullptr;
  }
}

static void
FreeSpectral(struct Variable *vars)
{
  for (int code = MaxCodes - 1; code >= 0; --code) FreeMemory(vars[code].spectral);
}

static void
FreeFourier(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code) FreeMemory(vars[code].fourier);
}

static void
FreeHybrid(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code) FreeMemory(vars[code].hybrid);
}

static void
FreeGrid(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code) FreeMemory(vars[code].grid);
}

static void
FreeSamp(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; ++code) FreeMemory(vars[code].samp);
}

// alloc_dp - Allocate space for double array
static double *
alloc_dp(int words, const char *array_name)
{
  double *result = nullptr;

  if (words > 0)
  {
    result = (double *) std::malloc(words * sizeof(double));
    if (result == nullptr) cdo_sys_error(array_name, "No Memory!");
  }

  return result;
}

static void
IniQuaSum(double *dest, const double *restrict src, int len)
{
  for (int i = 0; i < len; ++i) dest[i] = src[i] * src[i];
}

static void
AddQuaSum(double *dest, const double *restrict src, int len)
{
  for (int i = 0; i < len; ++i) dest[i] += src[i] * src[i];
}

static void
VarQuaSum(double *Variance, const double *restrict Sum, int len, int n)
{
  double rn1 = 1.0 / (n - 1);

  for (int i = 0; i < len; ++i) Variance[i] = (Variance[i] - Sum[i] * Sum[i] * n) * rn1;

  for (int i = 0; i < len; ++i) Variance[i] = (Variance[i] > 0.0) ? std::sqrt(Variance[i]) : 0.0;
}

static void
AddVector(double *dest, const double *src, size_t len, size_t *numMissVals, double missval)
{
  if (*numMissVals)
  {
    array_add_array_mv(len, dest, src, missval);
    *numMissVals = array_num_mv(len, dest, missval);
    if (*numMissVals == 0) *numMissVals = 1;
  }
  else { array_add_array(len, dest, src); }
}

static void
Add2Vectors(double *dest, const double *restrict srcA, const double *restrict srcB, size_t len)
{
  for (size_t i = 0; i < len; ++i) dest[i] = srcA[i] + srcB[i];
}

static void
Sub2Vectors(double *dest, const double *restrict srcA, const double *restrict srcB, size_t len)
{
  for (size_t i = 0; i < len; ++i) dest[i] = srcA[i] - srcB[i];
}

static void
MultVectorScalar(double *dest, const double *restrict src, double factor, size_t len, size_t numMissVals, double missval)
{
  if (numMissVals)
  {
    for (size_t i = 0; i < len; ++i) dest[i] = (is_equal(src[i], missval)) ? missval : src[i] * factor;
  }
  else
  {
    for (size_t i = 0; i < len; ++i) dest[i] = src[i] * factor;
  }
}

static void
DivVectorIvector(double *dest, const double *restrict src, const int *samp, size_t len, size_t *numMissVals, double missval)
{
  *numMissVals = 0;

  for (size_t i = 0; i < len; ++i)
  {
    if (is_equal(src[i], missval) || samp[i] == 0)
    {
      dest[i] = missval;
      *numMissVals = *numMissVals + 1;
    }
    else
      dest[i] = src[i] / samp[i];
  }
}

void
after_gp2sp(const AfterControl &globs, struct Variable *vars, int ccode)
{
  auto *var = &vars[ccode];

  if (var->spectral == nullptr)
  {
    if (var->hybrid == nullptr)
    {
      fprintf(stderr, "%d.hybrid not found\n", ccode);
      exit(99);
    }

    if (var->fourier == nullptr)
    {
      auto fieldSize = globs.DimFC * var->hlev;
      var->fourier = alloc_dp(fieldSize, "gp2sp.fourier");
      after_GP2FC(var->hybrid, var->fourier, globs.Latitudes, globs.Longitudes, var->hlev, globs.Fouriers);
    }

    var->spectral = alloc_dp(globs.Dim3SP, "gp2sp.spectral");
    fc2sp(var->fourier, var->spectral, globs.pold.data(), var->hlev, globs.Latitudes, globs.Fouriers, globs.Truncation);
  }
}

void
after_GP2FC(double *gp, double *fc, long nlat, long nlon, long nlev, long nfc)
{
  static long ifax[10];
  static double *trig = nullptr;

  if (ifax[9] != nlon)
  {
    if (trig) std::free(trig);
    trig = (double *) std::malloc(nlon * sizeof(double));
    int status = fft_set(trig, ifax, nlon);
    if (status < 0) exit(1);
  }

  gp2fc(trig, ifax, gp, fc, nlat, nlon, nlev, nfc);
}

void
after_FC2GP(double *fc, double *gp, long nlat, long nlon, long nlev, long nfc)
{
  static long ifax[10];
  static double *trig = nullptr;

  if (ifax[9] != nlon)
  {
    if (trig) std::free(trig);
    trig = (double *) std::malloc(nlon * sizeof(double));
    int status = fft_set(trig, ifax, nlon);
    if (status < 0) exit(1);
  }

  fc2gp(trig, ifax, fc, gp, nlat, nlon, nlev, nfc);
}

// HUMTEST

static void
sh2rh(int AnalysisData, double *sphum, double *rhum, double *t, int lev, int dimgpout, const double *level, double *fullpresshybrid)
{
  /* ***************************************************** */
  /* Define constants for calculation in presence of water */
  /* ***************************************************** */
  constexpr double RGAMW = (C_RCW - C_RCPV) / C_RV;
  constexpr double RBETW = C_RLVTT / C_RV + RGAMW * C_RTT;
  double RALPW = std::log(C_RESTT) + RBETW / C_RTT + RGAMW * std::log(C_RTT);

  /* ***************************************************** */
  /* Define constants for calculation in presence of  ice  */
  /* ***************************************************** */
  /*
  double RGAMS = (C_RCS - C_RCPV) / C_RV;
  double RBETS = C_RLSTT / C_RV + RGAMS * C_RTT;
  double RALPS = std::log(C_RESTT) + RBETS / C_RTT + RGAMS * std::log(C_RTT);
  */
  const double *fullp = AnalysisData ? level : fullpresshybrid;

  /***************************************************/
  /* Diagnostics of saturation water vapour pressure */
  /* over ice makes no sense, therefore ...          */
  /* Hint of Michael Ponater                08.10.97 */
  /***************************************************/
  constexpr double RGAM = RGAMW;
  constexpr double RBET = RBETW;
  double RALP = RALPW;
  for (int lp = 0; lp < lev; lp++)
  {
    for (int i = 0; i < dimgpout; ++i)
    {
      auto lpi = lp * dimgpout + i;
      auto lfp = (1 - AnalysisData) * lpi + AnalysisData * lp;
      /*
      if (t[lpi] < C_RTT) {
        RGAM = RGAMS; RBET = RBETS; RALP = RALPS;
      } else {
        RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
      }
      */
      auto es = (std::exp(RALP - RBET / t[lpi] - RGAM * std::log(t[lpi]))) / fullp[lfp];
      // qsat = es / (1. + C_RETV * (1. - es));
      auto qsatr = (1. + C_RETV * (1. - es)) / es;
      rhum[lpi] = sphum[lpi] * 100. * qsatr;
    }
  }
}

static void
rh2sh(double *sphum, double *rhum, double *t, int lev, int dimgpout, const double *level)
{
  /* ***************************************************** */
  /* Define constants for calculation in presence of water */
  /* ***************************************************** */
  constexpr double RGAMW = (C_RCW - C_RCPV) / C_RV;
  constexpr double RBETW = C_RLVTT / C_RV + RGAMW * C_RTT;
  double RALPW = std::log(C_RESTT) + RBETW / C_RTT + RGAMW * std::log(C_RTT);

  /* ***************************************************** */
  /* Define constants for calculation in presence of  ice  */
  /* ***************************************************** */
  // double RGAMS = (C_RCS - C_RCPV) / C_RV;
  // double RBETS = C_RLSTT / C_RV + RGAMS * C_RTT;
  // double RALPS = std::log(C_RESTT) + RBETS / C_RTT + RGAMS * std::log(C_RTT);

  /***************************************************/
  /* Diagnostics of saturation water vapour pressure */
  /* over ice makes no sense, therefore ...          */
  /* Hint of Michael Ponater                08.10.97 */
  /***************************************************/

  constexpr double RGAM = RGAMW;
  constexpr double RBET = RBETW;
  double RALP = RALPW;
  for (int lp = 0; lp < lev; lp++)
  {
    for (int i = 0; i < dimgpout; ++i)
    {
      auto lpi = lp * dimgpout + i;
      /*       if (t[lpi] < C_RTT) { */
      /*          RGAM = RGAMS; RBET = RBETS; RALP = RALPS; */
      /*       }  else { */
      /*          RGAM = RGAMW; RBET = RBETW; RALP = RALPW; */
      /*       } */
      auto es = (std::exp(RALP - RBET / t[lpi] - RGAM * std::log(t[lpi]))) / level[lp];
      auto qsat = es / (1. + C_RETV * (1. - es));
      sphum[lpi] = rhum[lpi] * qsat / 100.;
    }
  }
}

void
after_FCrh2FCsh(const AfterControl &globs, struct Variable *vars)
{
  auto fieldSize = globs.DimGP * globs.NumLevelRequest;

  if (vars[RHUMIDITY].grid == nullptr) vars[RHUMIDITY].grid = alloc_dp(fieldSize, "vars[RHUMIDITY].grid");
  if (vars[TEMPERATURE].grid == nullptr) vars[TEMPERATURE].grid = alloc_dp(fieldSize, "vars[TEMPERATURE].grid");
  if (vars[HUMIDITY].grid == nullptr) vars[HUMIDITY].grid = alloc_dp(fieldSize, "vars[HUMIDITY].grid");

  after_FC2GP(vars[RHUMIDITY].fourier, vars[RHUMIDITY].grid, globs.Latitudes, globs.Longitudes, vars[RHUMIDITY].plev,
              globs.Fouriers);
  after_FC2GP(vars[TEMPERATURE].fourier, vars[TEMPERATURE].grid, globs.Latitudes, globs.Longitudes, vars[TEMPERATURE].plev,
              globs.Fouriers);

  rh2sh(vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs.NumLevelRequest, globs.DimGP, globs.LevelRequest);

  after_GP2FC(vars[HUMIDITY].grid, vars[HUMIDITY].fourier, globs.Latitudes, globs.Longitudes, vars[HUMIDITY].plev, globs.Fouriers);

  FreeMemory(vars[HUMIDITY].grid);
  FreeMemory(vars[RHUMIDITY].grid);
  FreeMemory(vars[TEMPERATURE].grid);
}

void
after_SPuv2SPdv(const AfterControl &globs, struct Variable *vars)
{
  auto numLevels = globs.NumLevelRequest;
  auto numLats = globs.Latitudes;

  auto Div = vars[DIVERGENCE].spectral;
  auto Vor = vars[VORTICITY].spectral;
  auto DivOut = Div;
  auto VorOut = Vor;
  auto fieldSize = globs.DimFC * numLevels;

  auto &u_wind = vars[U_WIND];
  auto &v_wind = vars[V_WIND];
  if (u_wind.fourier == nullptr) u_wind.fourier = alloc_dp(fieldSize, "u_wind.fourier");
  if (v_wind.fourier == nullptr) v_wind.fourier = alloc_dp(fieldSize, "v_wind.fourier");

  sp2fc(u_wind.spectral, u_wind.fourier, globs.poli.data(), numLevels, numLats, globs.Fouriers, globs.Truncation);
  sp2fc(v_wind.spectral, v_wind.fourier, globs.poli.data(), numLevels, numLats, globs.Fouriers, globs.Truncation);
  uv2dv(u_wind.fourier, v_wind.fourier, Div, Vor, globs.pol2.data(), globs.pol3.data(), numLevels, numLats, globs.Truncation);

  FreeMemory(u_wind.fourier);
  FreeMemory(v_wind.fourier);

  for (int i = 0; i < numLevels; ++i)
  {
    sp2sp(Div, globs.Truncation, DivOut, globs.Truncation);
    sp2sp(Vor, globs.Truncation, VorOut, globs.Truncation);
    Div += globs.DimSP;
    Vor += globs.DimSP;
    DivOut += globs.DimSP;
    VorOut += globs.DimSP;
  }
}

void
after_FCsh2FCrh(const AfterControl &globs, struct Variable *vars)
{
  auto fieldSize = globs.DimGP * globs.NumLevelRequest;

  if (vars[RHUMIDITY].grid == nullptr) vars[RHUMIDITY].grid = alloc_dp(fieldSize, "vars[RHUMIDITY].grid");
  if (vars[TEMPERATURE].grid == nullptr) vars[TEMPERATURE].grid = alloc_dp(fieldSize, "vars[TEMPERATURE].grid");
  if (vars[HUMIDITY].grid == nullptr) vars[HUMIDITY].grid = alloc_dp(fieldSize, "vars[HUMIDITY].grid");

  after_FC2GP(vars[HUMIDITY].fourier, vars[HUMIDITY].grid, globs.Latitudes, globs.Longitudes, vars[HUMIDITY].plev, globs.Fouriers);
  after_FC2GP(vars[TEMPERATURE].fourier, vars[TEMPERATURE].grid, globs.Latitudes, globs.Longitudes, vars[TEMPERATURE].plev,
              globs.Fouriers);

  sh2rh(globs.AnalysisData, vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs.NumLevelRequest, globs.DimGP,
        globs.LevelRequest, vars[FULL_PRESS].hybrid);

  after_GP2FC(vars[RHUMIDITY].grid, vars[RHUMIDITY].fourier, globs.Latitudes, globs.Longitudes, vars[RHUMIDITY].plev,
              globs.Fouriers);

  FreeMemory(vars[HUMIDITY].grid);
  FreeMemory(vars[RHUMIDITY].grid);
  FreeMemory(vars[TEMPERATURE].grid);
}
/* ENDE HUMTEST */

static void
CheckAnalyses(struct Variable *vars)
{
  for (int code = 0; code < 272; ++code)
    if (vars[code].needed && code != DIVERGENCE && code != VORTICITY && code != STREAM && code != U_WIND && code != HUMIDITY
        && code != VELOPOT && code != V_WIND && code != RHUMIDITY && code != GEOPOTHEIGHT && code != PS
        && vars[code].spectral == nullptr && vars[code].grid == nullptr)
    {
      afterAbort("Code  %3d not found", code);
    }
}

// Process Pressure Level data
void
after_processPL(AfterControl &globs, struct Variable *vars)
{
  globs.MeanCount++;
  globs.TermCount++;

  if (globs.MeanCount == 1)
  {
    if (globs.Debug) fprintf(stderr, "CheckAnalyses: %d %d\n", globs.TermCount, globs.MeanCount);
    CheckAnalyses(vars);
    globs.StartDate = globs.OldDate;
  }
  if (globs.TermCount > 120) globs.Debug = 0;

  /* ============================== */
  /* Computations in spectral space */
  /* ============================== */

  if (vars[TEMPERATURE].needed)
  {
    vars[TEMPERATURE].hlev = 2;
    vars[TEMPERATURE].plev = globs.NumLevelRequest;
    vars[TEMPERATURE].sfit = true;
  }

  if (vars[GEOPOTHEIGHT].comp)
  {
    vars[GEOPOTHEIGHT].hlev = 2;
    vars[GEOPOTHEIGHT].plev = globs.NumLevelRequest;
    vars[GEOPOTHEIGHT].sfit = true;
  }

  if (vars[GEOPOTHEIGHT].comp && vars[GEOPOTENTIAL].detected)
  {
    if (vars[GEOPOTHEIGHT].spectral == nullptr)
      vars[GEOPOTHEIGHT].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "GEOPOTHEIGHT.spectral");
    MultVectorScalar(vars[GEOPOTHEIGHT].spectral, vars[GEOPOTENTIAL].spectral, C_RG, globs.DimSP * globs.NumLevelRequest, 0, 0);
    vars[GEOPOTENTIAL].needed = vars[GEOPOTENTIAL].selected;
  }

  if (globs.Type == 50 && vars[HUMIDITY].needed && vars[HUMIDITY].spectral == nullptr)
  {
    vars[HUMIDITY].plev = globs.NumLevelRequest;
    vars[HUMIDITY].sfit = true;
    vars[HUMIDITY].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[HUMIDITY].spectral");
    // SPrh2SPsh();
    vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
    vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
  }

  if (vars[U_WIND].spectral && vars[V_WIND].spectral && (vars[DIVERGENCE].comp || vars[VORTICITY].comp))
  {
    vars[DIVERGENCE].hlev = vars[VORTICITY].hlev = 2;
    vars[DIVERGENCE].plev = vars[VORTICITY].plev = globs.NumLevelRequest;
    vars[DIVERGENCE].sfit = vars[VORTICITY].sfit = true;
    if (vars[DIVERGENCE].spectral == nullptr)
      vars[DIVERGENCE].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[DIVERGENCE].spectral");
    if (vars[VORTICITY].spectral == nullptr)
      vars[VORTICITY].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[VORTICITY].spectral");
    after_SPuv2SPdv(globs, vars);
  }

  if (vars[U_WIND].comp || vars[V_WIND].comp)
  {
    vars[U_WIND].hlev = vars[V_WIND].hlev = 2;
    vars[U_WIND].plev = vars[V_WIND].plev = globs.NumLevelRequest;
    vars[U_WIND].sfit = vars[V_WIND].sfit = true;
    if (vars[U_WIND].spectral == nullptr)
      vars[U_WIND].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[U_WIND].spectral");
    if (vars[V_WIND].spectral == nullptr)
      vars[V_WIND].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[V_WIND].spectral");
    dv2uv(vars[DIVERGENCE].spectral, vars[VORTICITY].spectral, vars[U_WIND].spectral, vars[V_WIND].spectral, globs.dv2uv_f1.data(),
          globs.dv2uv_f2.data(), globs.Truncation, globs.DimSP, globs.NumLevelRequest);
  }

  if (vars[VELOPOT].comp)
  {
    vars[VELOPOT].hlev = 2;
    vars[VELOPOT].plev = globs.NumLevelRequest;
    vars[VELOPOT].sfit = true;
    if (vars[VELOPOT].spectral == nullptr)
      vars[VELOPOT].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[VELOPOT].spectral");
    dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, globs.NumLevelRequest, globs.Truncation);
  }

  if (vars[STREAM].comp)
  {
    vars[STREAM].hlev = 2;
    vars[STREAM].plev = globs.NumLevelRequest;
    vars[STREAM].sfit = true;
    if (vars[STREAM].spectral == nullptr)
      vars[STREAM].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[STREAM].spectral");
    dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral, globs.NumLevelRequest, globs.Truncation);
  }

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if (globs.Type == 50)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].spectral) afterAbort("Code %d not available on spectral space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimSP;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].spectral + offset, 0);
        }
      }

    FreeSpectral(vars);
    return;
  }

  /* =============================== */
  /* Transformation to fourier space */
  /* Computations in fourier space   */
  /* =============================== */

  if (globs.Type >= 60)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].spectral)
      {
        if (vars[code].fourier == nullptr)
        {
          auto fieldSize = vars[code].plev * globs.DimFC;
          vars[code].fourier = alloc_dp(fieldSize, FieldName(code, "fourier"));
        }
        sp2fc(vars[code].spectral, vars[code].fourier, globs.poli.data(), vars[code].plev, globs.Latitudes, globs.Fouriers,
              globs.Truncation);
      }
    if (vars[U_WIND].needed && vars[U_WIND].fourier)
      scaluv(vars[U_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);
    if (vars[V_WIND].needed && vars[V_WIND].fourier)
      scaluv(vars[V_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);

    /* HUMTEST */
    if (globs.Type < 70 && vars[HUMIDITY].needed && vars[HUMIDITY].fourier == nullptr)
    {
      vars[HUMIDITY].plev = globs.NumLevelRequest;
      vars[HUMIDITY].sfit = true;
      vars[HUMIDITY].fourier = alloc_dp(globs.DimFC * globs.NumLevelRequest, "vars[HUMIDITY].fourier");

      after_FCrh2FCsh(globs, vars);

      vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    }

    if (globs.Type < 70 && vars[RHUMIDITY].needed && vars[RHUMIDITY].fourier == nullptr)
    {
      vars[RHUMIDITY].plev = globs.NumLevelRequest;
      vars[RHUMIDITY].sfit = true;
      vars[RHUMIDITY].fourier = alloc_dp(globs.DimFC * globs.NumLevelRequest, "vars[RHUMIDITY].fourier");

      after_FCsh2FCrh(globs, vars);

      vars[HUMIDITY].needed = vars[HUMIDITY].selected;
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    }
    /* ENDE HUMTEST */
  }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if (globs.Type == 60)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on fourier space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, 0);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if (globs.Type == 61)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on zonal mean!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, 0);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ============================ */
  /* Transformation to gridpoints */
  /* ============================ */

  if (vars[PS].comp && vars[LNPS].grid)
  {
    if (vars[PS].grid == nullptr) vars[PS].grid = alloc_dp(globs.DimGP, "Ps");
    for (int l = 0; l < globs.DimGP; ++l) vars[PS].grid[l] = std::exp(vars[LNPS].grid[l]);
  }

  if (globs.Type >= 70)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].fourier)
      {
        if (vars[code].grid == nullptr)
        {
          auto fieldSize = vars[code].plev * globs.DimGP;
          vars[code].grid = alloc_dp(fieldSize, FieldName(code, "grid"));
        }

        after_FC2GP(vars[code].fourier, vars[code].grid, globs.Latitudes, globs.Longitudes, vars[code].plev, globs.Fouriers);
      }
  }

  FreeFourier(vars);

  /* HUMTEST */
  /* -------------------------------- */
  /*  Computation in gridpoint space  */
  /* -------------------------------- */

  if (vars[RHUMIDITY].needed && vars[RHUMIDITY].grid == nullptr)
  {
    vars[RHUMIDITY].plev = globs.NumLevelRequest;
    vars[RHUMIDITY].sfit = true;
    vars[RHUMIDITY].grid = alloc_dp(globs.DimGP * globs.NumLevelRequest, "vars[RHUMIDITY].grid");
    sh2rh(globs.AnalysisData, vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs.NumLevelRequest, globs.DimGP,
          globs.LevelRequest, vars[FULL_PRESS].hybrid);
    vars[HUMIDITY].needed = vars[HUMIDITY].selected;
    vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
  }

  if (vars[HUMIDITY].needed && vars[HUMIDITY].grid == nullptr)
  {
    vars[HUMIDITY].plev = globs.NumLevelRequest;
    vars[HUMIDITY].sfit = true;
    vars[HUMIDITY].grid = alloc_dp(globs.DimGP * globs.NumLevelRequest, "vars[HUMIDITY].grid");
    rh2sh(vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs.NumLevelRequest, globs.DimGP,
          globs.LevelRequest);
    vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
    vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
  }
  /* HUMTEST ENDE */

  /* -------------------------- */
  /*  Computation of Means      */
  /* -------------------------- */

  if (globs.Mean)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].grid)
      {
        auto fieldSize = globs.DimGP * vars[code].plev;
        if (vars[code].mean == nullptr) vars[code].mean = alloc_dp(fieldSize, FieldName(code, "mean"));

        if (globs.MeanCount == 1)
          array_copy(fieldSize, vars[code].grid, vars[code].mean);
        else
          AddVector(vars[code].mean, vars[code].grid, fieldSize, &vars[code].numMissVals, vars[code].missval);

        if (globs.EndOfInterval)
        {
          if (vars[code].samp == nullptr)
            MultVectorScalar(vars[code].mean, vars[code].mean, 1.0 / globs.MeanCount, fieldSize, vars[code].numMissVals,
                             vars[code].missval);
          else
            DivVectorIvector(vars[code].mean, vars[code].mean, vars[code].samp, fieldSize, &vars[code].numMissVals,
                             vars[code].missval);
        }
      }

  /* -------------------------- */
  /*  Computation of Variances  */
  /* -------------------------- */

  if (globs.Mean > 1)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].mean)
      {
        auto fieldSize = globs.DimGP * vars[code].plev;
        if (vars[code].variance == nullptr) vars[code].variance = alloc_dp(fieldSize, FieldName(code, "var"));
        if (globs.MeanCount == 1)
          IniQuaSum(vars[code].variance, vars[code].grid, fieldSize);
        else
          AddQuaSum(vars[code].variance, vars[code].grid, fieldSize);

        if (globs.EndOfInterval) VarQuaSum(vars[code].variance, vars[code].mean, fieldSize, globs.MeanCount);
      }

  if (globs.Mean && !globs.EndOfInterval)
  {
    FreeGrid(vars);
    return;
  }

  /* ---------------------------------------------- */
  /*  Output of pressure level means and variances  */
  /* ---------------------------------------------- */

  if (globs.Type == 70 && globs.Mean && globs.EndOfInterval)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimGP;
          if (globs.Mean != 2)
          {
            streamDefField(globs.ostreamID, vars[code].ovarID, k);
            streamWriteField(globs.ostreamID, vars[code].mean + offset, vars[code].numMissVals);
          }
          if (globs.Mean >= 2)
          {
            streamDefField(globs.ostreamID2, vars[code].ovarID2, k);
            streamWriteField(globs.ostreamID2, vars[code].variance + offset, vars[code].numMissVals);
          }
        }
      }

    FreeSamp(vars);
    FreeGrid(vars);
    return;
  }

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if (globs.Type == 70)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimGP;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].grid + offset, vars[code].numMissVals);
        }
      }

    FreeGrid(vars);
    return;
  }
}

static void
theta(double *pthetaf, double *pthetah, const double *ph, const double *ps, const double *tf, const double *ts, int levels,
      int dimgp, int dim3gp)
{
  auto thetah = pthetah;
  auto thetaf = pthetaf;

  auto kappa = PlanetRD / C_RCPD;

  for (int h = 0; h < dimgp; h++) thetah[h] = 0.0;
  thetah += dimgp;
  for (int l = 0; l < levels - 1; ++l)
  {
    for (int h = 0; h < dimgp; h++) thetah[h] = 0.5 * (tf[h] + tf[h + dimgp]) * std::pow((ps[h] / ph[h]), kappa);

    ph += dimgp;
    tf += dimgp;
    thetah += dimgp;
  }

  array_copy(dimgp, ts, thetah);
  thetah = pthetah;
  for (int h = 0; h < dim3gp; h++) thetaf[h] = 0.5 * (thetah[h] + thetah[h + dimgp]);
}

static void
windSpeed(double *uvspeed, double *u, double *v, int dim3gp)
{
  for (int i = 0; i < dim3gp; ++i) uvspeed[i] = std::sqrt(u[i] * u[i] + v[i] * v[i]);
}

static void
Omega(double *omega_in, const double *divergence, const double *u_wind, const double *v_wind, const double *halfpress,
      const double *fullpress, double *dpsdx, const double *dpsdy, const double *vct, int dimgp, int nlev)
{
  // Compute half level part of vertical velocity

  for (int i = 0; i < dimgp; ++i) omega_in[i] = 0.0;

  for (int j = 0; j < nlev; ++j)
  {
    auto omega = omega_in + j * dimgp;
    const auto halfp = halfpress + j * dimgp;
    const auto diver = divergence + j * dimgp;
    const auto uwind = u_wind + j * dimgp;
    const auto vwind = v_wind + j * dimgp;

    auto deltaHybrid = vct[nlev + j + 2] - vct[nlev + j + 1];
    for (int i = 0; i < dimgp; ++i)
    {
      omega[i + dimgp]
          = omega[i] - diver[i] * (halfp[i + dimgp] - halfp[i]) - deltaHybrid * (uwind[i] * dpsdx[i] + vwind[i] * dpsdy[i]);
    }
  }

  // interpolate to full levels

  for (int j = 0; j < nlev; ++j)
  {
    auto omega = omega_in + j * dimgp;
    for (int i = 0; i < dimgp; ++i) omega[i] = 0.5 * (omega[i] + omega[i + dimgp]);
  }

  // compute full level part of vertical velocity

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (int j = 0; j < nlev; ++j)
  {
    auto omega = omega_in + j * dimgp;
    const auto halfp = halfpress + j * dimgp;
    const auto fullp = fullpress + j * dimgp;
    const auto uwind = u_wind + j * dimgp;
    const auto vwind = v_wind + j * dimgp;

    auto deltaHybrid = vct[nlev + j + 2] - vct[nlev + j + 1];
    if (std::fabs(deltaHybrid) > 0)
    {
      auto cterm = vct[j + 1] * vct[nlev + j + 1] - vct[j] * vct[nlev + j + 2];
      auto ctermIsValid = is_not_equal(cterm, 0.0);
      for (int i = 0; i < dimgp; ++i)
      {
        auto pterm = ctermIsValid ? cterm / (halfp[i + dimgp] - halfp[i]) * std::log(halfp[i + dimgp] / halfp[i]) : 0.0;
        omega[i] += fullp[i] * (uwind[i] * dpsdx[i] + vwind[i] * dpsdy[i]) / (halfp[i + dimgp] - halfp[i]) * (deltaHybrid + pterm);
      }
    }
  }
}

void
geopot_height_half(double *gheight, const double *ta_fl, const double *hus_fl, const double *p_hl, long ngp, long nlev)
{
  // Computes geopotential height on half level

  // gheight: geopotential height (nlev+1) (bottom level contains surface geopotential)
  // ta_fl:   air temperature on full levels
  // hus_fl:  specific humidity on full levels
  // p_hl:    pressure on half levels

  double z2log2 = 2.0 * std::log(2.0);
  double vtmp = (C_RV / PlanetRD) - 1.0;

  if (hus_fl)  // Humidity is present
  {
    for (long k = nlev; k > 1; k--)
    {
      auto offset = ngp * (k - 1);
      auto gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      const auto hus = hus_fl + offset;
      for (long i = 0; i < ngp; ++i) gh[i] = gh[i + ngp] + PlanetRD * ta[i] * (1.0 + vtmp * hus[i]) * std::log(p[i + ngp] / p[i]);
    }

    // top level
    for (long i = 0; i < ngp; ++i) gheight[i] = gheight[i + ngp] + PlanetRD * ta_fl[i] * (1.0 + vtmp * hus_fl[i]) * z2log2;
  }
  else  // No humidity
  {
    for (long k = nlev; k > 1; k--)
    {
      auto offset = ngp * (k - 1);
      auto gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      for (long i = 0; i < ngp; ++i) gh[i] = gh[i + ngp] + PlanetRD * ta[i] * std::log(p[i + ngp] / p[i]);
    }

    // top level
    for (long i = 0; i < ngp; ++i) gheight[i] = gheight[i + ngp] + PlanetRD * ta_fl[i] * z2log2;
  }

  double zrg = 1.0 / PlanetGrav;
  for (long i = 0; i < ngp * (nlev + 1); ++i) gheight[i] *= zrg;
}

void
geopot_height_full(double *gheight, const double *ta_fl, const double *hus_fl, const double *p_hl, long ngp, long nlev)
{
  // Computes geopotential height on full level

  // gheight: geopotential height (nlev+1) (bottom level contains surface geopotential)
  // ta_fl:   air temperature on full levels
  // hus_fl:  specific humidity on full levels
  // p_hl:    pressure on half levels

  double zlog2 = std::log(2.0);
  double vtmp = (C_RV / PlanetRD) - 1.0;

  if (hus_fl)  // Humidity is present
  {
    for (long k = nlev; k > 1; k--)
    {
      auto offset = ngp * (k - 1);
      auto gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      const auto hus = hus_fl + offset;
      for (long i = 0; i < ngp; ++i) gh[i] = gh[i + ngp] + PlanetRD * ta[i] * (1.0 + vtmp * hus[i]) * std::log(p[i + ngp] / p[i]);
    }

    // top level
    for (long i = 0; i < ngp; ++i) gheight[i] = gheight[i + ngp] + PlanetRD * ta_fl[i] * (1.0 + vtmp * hus_fl[i]) * zlog2;

    for (long k = 1; k < nlev; ++k)
    {
      long offset = ngp * k;
      double *gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      const auto hus = hus_fl + offset;
      for (long i = 0; i < ngp; ++i)
        gh[i] = gh[i + ngp]
                + PlanetRD * ta[i] * (1.0 + vtmp * hus[i]) * (1.0 - (p[i] / (p[i + ngp] - p[i])) * std::log(p[i + ngp] / p[i]));
    }
  }
  else  // No humidity
  {
    for (long k = nlev; k > 1; k--)
    {
      auto offset = ngp * (k - 1);
      auto gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      for (long i = 0; i < ngp; ++i) gh[i] = gh[i + ngp] + PlanetRD * ta[i] * std::log(p[i + ngp] / p[i]);
    }

    // top level
    for (long i = 0; i < ngp; ++i) gheight[i] = gheight[i + ngp] + PlanetRD * ta_fl[i] * zlog2;

    for (long k = 1; k < nlev; ++k)
    {
      auto offset = ngp * k;
      auto gh = gheight + offset;
      const auto p = p_hl + offset;
      const auto ta = ta_fl + offset;
      for (long i = 0; i < ngp; ++i)
        gh[i] = gh[i + ngp] + PlanetRD * ta[i] * (1.0 - (p[i] / (p[i + ngp] - p[i])) * std::log(p[i + ngp] / p[i]));
    }
  }

  double zrg = 1.0 / PlanetGrav;
  for (long i = 0; i < ngp * (nlev + 1); ++i) gheight[i] *= zrg;
}

constexpr double SCALESLP = 101325.0;

/* ======================================== */
/* LayerWater integral liquid water content */
/* ======================================== */

void
LayerWater(double *ww, double *ll, double pmax, double pmin, int DimGP, int HalfLevels, const double *vct)
{
  double pph[MaxLevel]{};

  int k;
  for (k = 0; k < HalfLevels; ++k) pph[k] = vct[k] + vct[k + HalfLevels] * SCALESLP;
  for (k = 0; k < HalfLevels; ++k)
    if (pph[k] > pmax) break;
  auto MaxLev = k - 1;
  for (k = HalfLevels - 1; k >= 0; k--)
    if (pph[k] < pmin) break;
  auto MinLev = k;

  std::ranges::fill_n(ll, DimGP, 0.0);

  for (k = MaxLev; k <= MinLev; ++k)
  {
    for (int i = 0; i < DimGP; ++i) ll[i] += ww[i + k * DimGP] * (pph[k + 1] - pph[k]);
  }
  for (int i = 0; i < DimGP; ++i) ll[i] /= PlanetGrav;
}

/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

void
LayerCloud(double *cc, double *ll, double pmax, double pmin, int DimGP, int HalfLevels, const double *vct)
{
  double pph[MaxLevel]{};
  constexpr double ZEPSEC = 1.0e-12;

  int k;
  for (k = 0; k < HalfLevels; ++k) pph[k] = vct[k] + vct[k + HalfLevels] * SCALESLP;
  for (k = 0; k < HalfLevels; ++k)
    if (pph[k] > pmax) break;
  auto MaxLev = k - 1;
  for (k = HalfLevels - 1; k >= 0; k--)
    if (pph[k] < pmin) break;
  auto MinLev = k;

  for (int i = 0; i < DimGP; ++i) ll[i] = 1. - cc[i + MaxLev * DimGP];

  for (k = MaxLev + 1; k <= MinLev; ++k)
  {
    for (int i = 0; i < DimGP; ++i)
      ll[i] *= (1. - std::max(cc[i + (k - 1) * DimGP], cc[i + k * DimGP])) / (1. - std::min(cc[i + (k - 1) * DimGP], 1. - ZEPSEC));
  }
  for (int i = 0; i < DimGP; ++i) ll[i] = 1. - ll[i];
}

// Grid Point Computations
void
after_EchamCompGP(const AfterControl &globs, struct Variable *vars)
{
  if (vars[GEOPOTHEIGHT].comp || vars[SLP].comp || vars[THETAF].needed || vars[HALF_PRESS].needed || vars[RHUMIDITY].comp
      || vars[OMEGA].comp || globs.Type >= 30)
  {
    if (vars[FULL_PRESS].hybrid == nullptr) vars[FULL_PRESS].hybrid = alloc_dp(globs.Dim3GP, "vars[FULL_PRESS].hybrid");

    vars[HALF_PRESS].hlev = globs.NumLevel + 1;
    vars[HALF_PRESS].plev = globs.NumLevelRequest;
    vars[HALF_PRESS].sfit = false;

    if (vars[HALF_PRESS].hybrid == nullptr)
      vars[HALF_PRESS].hybrid = alloc_dp(globs.Dim3GP + globs.DimGP, "vars[HALF_PRESS].hybrid");

    vct_to_hybrid_pressure(vars[FULL_PRESS].hybrid, vars[HALF_PRESS].hybrid, globs.vct, vars[PS_PROG].hybrid, globs.NumLevel,
                           globs.DimGP);
  }

  if (globs.unitsel > 2) FreeMemory(vars[FULL_PRESS].hybrid);

  if (vars[THETAF].needed)
  {
    auto &thetaf = vars[THETAF];
    thetaf.hlev = globs.NumLevel;
    thetaf.plev = globs.NumLevelRequest;
    thetaf.sfit = true;
    if (thetaf.hybrid == nullptr) thetaf.hybrid = alloc_dp(globs.Dim3GP, "vars[THETAF].hybrid");
    if (vars[THETAH].hybrid == nullptr) vars[THETAH].hybrid = alloc_dp(globs.Dim3GP, "vars[THETAH].hybrid");
    theta(thetaf.hybrid, vars[THETAH].hybrid, vars[HALF_PRESS].hybrid, vars[PS_PROG].hybrid, vars[TEMPERATURE].hybrid,
          vars[TS].hybrid, globs.NumLevel, globs.DimGP, globs.Dim3GP);
  }

  if (vars[GEOPOTHEIGHT].comp)
  {
    auto &geopotheight = vars[GEOPOTHEIGHT];
    geopotheight.hlev = globs.NumLevel + 1;
    geopotheight.plev = globs.NumLevelRequest;
    geopotheight.sfit = true;
    geopotheight.hybrid = alloc_dp(globs.Dim3GP + globs.DimGP, "vars[GEOPOTHEIGHT].hybrid");

    array_copy(globs.DimGP, globs.orography.data(), geopotheight.hybrid + globs.Dim3GP);
    geopot_height_half(geopotheight.hybrid, vars[TEMPERATURE].hybrid, vars[HUMIDITY].hybrid, vars[HALF_PRESS].hybrid, globs.DimGP,
                       globs.NumLevel);

    vars[HUMIDITY].needed = vars[HUMIDITY].selected;
  }
  else if (vars[GEOPOTHEIGHT].hybrid && vars[GEOPOTHEIGHT].hlev == globs.NumLevel)
  {
    auto &geopotheight = vars[GEOPOTHEIGHT];
    geopotheight.hlev = globs.NumLevel + 1;
    geopotheight.sfit = true;
    geopotheight.hybrid = (double *) std::realloc(geopotheight.hybrid, (globs.Dim3GP + globs.DimGP) * sizeof(double));
    array_copy(globs.DimGP, globs.orography.data(), geopotheight.hybrid + globs.Dim3GP);
    for (int i = 0; i < globs.DimGP; ++i) geopotheight.hybrid[globs.Dim3GP + i] /= PlanetGrav;
  }

  if (vars[DPSDX].needed || vars[DPSDY].needed)
    for (int l = 0; l < globs.DimGP; ++l)
    {
      vars[DPSDX].hybrid[l] *= vars[PS_PROG].hybrid[l];
      vars[DPSDY].hybrid[l] *= vars[PS_PROG].hybrid[l];
    }

  if (vars[OMEGA].comp)
  {
    auto &omega = vars[OMEGA];
    omega.hlev = globs.NumLevel + 1;
    omega.plev = globs.NumLevelRequest;
    omega.sfit = true;
    omega.hybrid = alloc_dp(globs.Dim3GP + globs.DimGP, "omega.hybrid");

    Omega(omega.hybrid, vars[DIVERGENCE].hybrid, vars[U_WIND].hybrid, vars[V_WIND].hybrid, vars[HALF_PRESS].hybrid,
          vars[FULL_PRESS].hybrid, vars[DPSDX].hybrid, vars[DPSDY].hybrid, globs.vct.data(), globs.DimGP, globs.NumLevel);

    vars[DPSDX].needed = vars[DPSDX].selected;
    vars[DPSDY].needed = vars[DPSDY].selected;
  }

  if (vars[WINDSPEED].comp)
  {
    auto &windspeed = vars[WINDSPEED];
    windspeed.hlev = globs.NumLevel;
    windspeed.plev = globs.NumLevelRequest;
    windspeed.sfit = true;
    windspeed.hybrid = alloc_dp(globs.Dim3GP, "windspeed.hybrid");

    windSpeed(windspeed.hybrid, vars[U_WIND].hybrid, vars[V_WIND].hybrid, globs.Dim3GP);
  }

  if (vars[RHUMIDITY].comp)
  {
    auto &rhumidity = vars[RHUMIDITY];
    rhumidity.hlev = globs.NumLevel;
    rhumidity.plev = globs.NumLevelRequest;
    rhumidity.sfit = false;
    rhumidity.hybrid = alloc_dp(globs.Dim3GP, "rhumidity.hybrid");

    sh2rh(globs.AnalysisData, vars[HUMIDITY].hybrid, rhumidity.hybrid, vars[TEMPERATURE].hybrid, globs.NumLevel, globs.DimGP,
          globs.LevelRequest, vars[FULL_PRESS].hybrid);

    vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    vars[HUMIDITY].needed = vars[HUMIDITY].selected;
  }

  if (vars[PS].comp)
  {
    auto &ps = vars[PS];
    ps.hlev = 1;
    ps.plev = 1;
    ps.sfit = true;  // ???
    ps.hybrid = alloc_dp(globs.DimGP, "ps.hybrid");
    for (int l = 0; l < globs.DimGP; ++l) ps.hybrid[l] = std::exp(vars[LNPS].hybrid[l]);
  }

  if (vars[SLP].comp)
  {
    auto &slp = vars[SLP];
    slp.hlev = 1;
    slp.plev = 1;
    slp.sfit = true;
    slp.hybrid = alloc_dp(globs.DimGP, "slp.hybrid");

    extrapolate_P(slp.hybrid, vars[HALF_PRESS].hybrid + globs.Dim3GP, vars[FULL_PRESS].hybrid + globs.Dim3GP - globs.DimGP,
                  globs.orography.data(), vars[TEMPERATURE].hybrid + globs.Dim3GP - globs.DimGP, globs.DimGP);
    vars[TEMPERATURE].needed = vars[TEMPERATURE].selected || vars[GEOPOTHEIGHT].selected;
  }

  if (vars[PRECIP].comp)
  {
    auto &precip = vars[PRECIP];
    precip.hlev = precip.plev = 1;
    precip.sfit = false;
    precip.hybrid = alloc_dp(globs.DimGP, "precip.hybrid");
    Add2Vectors(precip.hybrid, vars[142].hybrid, vars[143].hybrid, globs.DimGP);
  }

  if (vars[NET_TOP].comp)
  {
    auto &net_top = vars[NET_TOP];
    net_top.hlev = net_top.plev = 1;
    net_top.sfit = false;
    net_top.hybrid = alloc_dp(globs.DimGP, "net_top.hybrid");
    Add2Vectors(net_top.hybrid, vars[178].hybrid, vars[179].hybrid, globs.DimGP);
  }

  if (vars[NET_BOT].comp)
  {
    auto &net_bot = vars[NET_BOT];
    net_bot.hlev = net_bot.plev = 1;
    net_bot.sfit = false;
    net_bot.hybrid = alloc_dp(globs.DimGP, "net_bot.hybrid");
    Add2Vectors(net_bot.hybrid, vars[176].hybrid, vars[177].hybrid, globs.DimGP);
  }

  if (vars[NET_HEAT].comp)
  {
    auto &net_heat = vars[NET_HEAT];
    net_heat.hlev = net_heat.plev = 1;
    net_heat.sfit = false;
    net_heat.hybrid = alloc_dp(globs.DimGP, "net_heat.hybrid");
    /*
    if (Source == ECHAM5_Source)
      {
        MultVectorScalar(net_heat.hybrid, vars[218].hybrid, (-3.345e5), globs.DimGP, vars[218].numMissVals, vars[218].missval);
        Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[176].hybrid, globs.DimGP);
        Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[177].hybrid, globs.DimGP);
        Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[146].hybrid, globs.DimGP);
        Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[147].hybrid, globs.DimGP);
        Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[206].hybrid, globs.DimGP);
        Sub2Vectors(net_heat.hybrid, net_heat.hybrid, vars[208].hybrid, globs.DimGP);
        Sub2Vectors(net_heat.hybrid, net_heat.hybrid, vars[209].hybrid, globs.DimGP);
      }
    else
    */
    {
      MultVectorScalar(net_heat.hybrid, vars[218].hybrid, C_TIMES_RHOH2O, globs.DimGP, vars[218].numMissVals, vars[218].missval);
      Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[176].hybrid, globs.DimGP);
      Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[177].hybrid, globs.DimGP);
      Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[146].hybrid, globs.DimGP);
      Add2Vectors(net_heat.hybrid, net_heat.hybrid, vars[147].hybrid, globs.DimGP);
      Sub2Vectors(net_heat.hybrid, net_heat.hybrid, vars[220].hybrid, globs.DimGP);
    }
  }

  if (vars[NET_WATER].comp)
  {
    auto &net_water = vars[NET_WATER];
    net_water.hlev = net_water.plev = 1;
    net_water.sfit = false;
    net_water.hybrid = alloc_dp(globs.DimGP, "net_water.hybrid");
    Sub2Vectors(net_water.hybrid, vars[182].hybrid, vars[160].hybrid, globs.DimGP);
    Add2Vectors(net_water.hybrid, net_water.hybrid, vars[142].hybrid, globs.DimGP);
    Add2Vectors(net_water.hybrid, net_water.hybrid, vars[143].hybrid, globs.DimGP);
  }

  if (vars[LOW_WATER].comp)
  {
    vars[LOW_WATER].hlev = vars[LOW_WATER].plev = 1;
    vars[LOW_WATER].sfit = false;
    vars[LOW_WATER].hybrid = alloc_dp(globs.DimGP, "vars[LOW_WATER].hybrid");
    LayerWater(vars[222].hybrid, vars[LOW_WATER].hybrid, 75000., 101300., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[MID_WATER].comp)
  {
    vars[MID_WATER].hlev = vars[MID_WATER].plev = 1;
    vars[MID_WATER].sfit = false;
    vars[MID_WATER].hybrid = alloc_dp(globs.DimGP, "vars[MID_WATER].hybrid");
    LayerWater(vars[222].hybrid, vars[MID_WATER].hybrid, 46000., 73000., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[HIH_WATER].comp)
  {
    vars[HIH_WATER].hlev = vars[HIH_WATER].plev = 1;
    vars[HIH_WATER].sfit = false;
    vars[HIH_WATER].hybrid = alloc_dp(globs.DimGP, "vars[HIH_WATER].hybrid");
    LayerWater(vars[222].hybrid, vars[HIH_WATER].hybrid, 5000., 44000., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[ALL_WATER].comp)
  {
    vars[ALL_WATER].hlev = vars[ALL_WATER].plev = 1;
    vars[ALL_WATER].sfit = false;
    vars[ALL_WATER].hybrid = alloc_dp(globs.DimGP, "vars[ALL_WATER].hybrid");
    LayerWater(vars[222].hybrid, vars[ALL_WATER].hybrid, 5000., 101300., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[LOW_CLOUD].comp)
  {
    vars[LOW_CLOUD].hlev = vars[LOW_CLOUD].plev = 1;
    vars[LOW_CLOUD].sfit = false;
    vars[LOW_CLOUD].hybrid = alloc_dp(globs.DimGP, "vars[LOW_CLOUD].hybrid");
    LayerCloud(vars[223].hybrid, vars[LOW_CLOUD].hybrid, 75000., 101300., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[MID_CLOUD].comp)
  {
    vars[MID_CLOUD].hlev = vars[MID_CLOUD].plev = 1;
    vars[MID_CLOUD].sfit = false;
    vars[MID_CLOUD].hybrid = alloc_dp(globs.DimGP, "vars[MID_CLOUD].hybrid");
    LayerCloud(vars[223].hybrid, vars[MID_CLOUD].hybrid, 46000., 73000., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[HIH_CLOUD].comp)
  {
    vars[HIH_CLOUD].hlev = vars[HIH_CLOUD].plev = 1;
    vars[HIH_CLOUD].sfit = false;
    vars[HIH_CLOUD].hybrid = alloc_dp(globs.DimGP, "vars[HIH_CLOUD].hybrid");
    LayerCloud(vars[223].hybrid, vars[HIH_CLOUD].hybrid, 5000., 44000., globs.DimGP, globs.HalfLevels, globs.vct.data());
  }

  if (vars[SW_CLF].comp)
  {
    vars[SW_CLF].hlev = vars[SW_CLF].plev = 1;
    vars[SW_CLF].sfit = false;
    vars[SW_CLF].hybrid = alloc_dp(globs.DimGP, "SW_CLF.hybrid");
    Sub2Vectors(vars[SW_CLF].hybrid, vars[178].hybrid, vars[224].hybrid, globs.DimGP);
  }

  if (vars[SW_BOT_CLF].comp)
  {
    vars[SW_BOT_CLF].hlev = vars[SW_BOT_CLF].plev = 1;
    vars[SW_BOT_CLF].sfit = false;
    vars[SW_BOT_CLF].hybrid = alloc_dp(globs.DimGP, "vars[SW_BOT_CLF].hybrid");
    Sub2Vectors(vars[SW_BOT_CLF].hybrid, vars[176].hybrid, vars[185].hybrid, globs.DimGP);
  }

  if (vars[SW_TOP_CLF].comp)
  {
    vars[SW_TOP_CLF].hlev = vars[SW_TOP_CLF].plev = 1;
    vars[SW_TOP_CLF].sfit = false;
    vars[SW_TOP_CLF].hybrid = alloc_dp(globs.DimGP, "vars[SW_TOP_CLF].hybrid");
    Sub2Vectors(vars[SW_TOP_CLF].hybrid, vars[178].hybrid, vars[187].hybrid, globs.DimGP);
  }

  if (vars[LW_CLF].comp)
  {
    vars[LW_CLF].hlev = vars[LW_CLF].plev = 1;
    vars[LW_CLF].sfit = false;
    vars[LW_CLF].hybrid = alloc_dp(globs.DimGP, "LW_CLF.hybrid");
    Sub2Vectors(vars[LW_CLF].hybrid, vars[179].hybrid, vars[225].hybrid, globs.DimGP);
  }

  if (vars[LW_BOT_CLF].comp)
  {
    vars[LW_BOT_CLF].hlev = vars[LW_BOT_CLF].plev = 1;
    vars[LW_BOT_CLF].sfit = false;
    vars[LW_BOT_CLF].hybrid = alloc_dp(globs.DimGP, "vars[LW_BOT_CLF].hybrid");
    Sub2Vectors(vars[LW_BOT_CLF].hybrid, vars[177].hybrid, vars[186].hybrid, globs.DimGP);
  }

  if (vars[LW_TOP_CLF].comp)
  {
    vars[LW_TOP_CLF].hlev = vars[LW_TOP_CLF].plev = 1;
    vars[LW_TOP_CLF].sfit = false;
    vars[LW_TOP_CLF].hybrid = alloc_dp(globs.DimGP, "vars[LW_TOP_CLF].hybrid");
    Sub2Vectors(vars[LW_TOP_CLF].hybrid, vars[179].hybrid, vars[188].hybrid, globs.DimGP);
  }

  if (vars[NET_CLF].comp)
  {
    vars[NET_CLF].hlev = vars[NET_CLF].plev = 1;
    vars[NET_CLF].sfit = false;
    vars[NET_CLF].hybrid = alloc_dp(globs.DimGP, "NET_CLF.hybrid");
    Add2Vectors(vars[NET_CLF].hybrid, vars[178].hybrid, vars[179].hybrid, globs.DimGP);
    Sub2Vectors(vars[NET_CLF].hybrid, vars[NET_CLF].hybrid, vars[224].hybrid, globs.DimGP);
    Sub2Vectors(vars[NET_CLF].hybrid, vars[NET_CLF].hybrid, vars[225].hybrid, globs.DimGP);
  }

  if (vars[SW_ATM].comp)
  {
    vars[SW_ATM].hlev = vars[SW_ATM].plev = 1;
    vars[SW_ATM].sfit = false;
    vars[SW_ATM].hybrid = alloc_dp(globs.DimGP, "vars[SW_ATM].hybrid");
    Sub2Vectors(vars[SW_ATM].hybrid, vars[178].hybrid, vars[176].hybrid, globs.DimGP);
  }

  if (vars[LW_ATM].comp)
  {
    vars[LW_ATM].hlev = vars[LW_ATM].plev = 1;
    vars[LW_ATM].sfit = false;
    vars[LW_ATM].hybrid = alloc_dp(globs.DimGP, "vars[LW_ATM].hybrid");
    Sub2Vectors(vars[LW_ATM].hybrid, vars[179].hybrid, vars[177].hybrid, globs.DimGP);
  }

  if (vars[NET_ATM].comp)
  {
    vars[NET_ATM].hlev = vars[NET_ATM].plev = 1;
    vars[NET_ATM].sfit = false;
    vars[NET_ATM].hybrid = alloc_dp(globs.DimGP, "vars[NET_ATM].hybrid");
    Add2Vectors(vars[NET_ATM].hybrid, vars[178].hybrid, vars[179].hybrid, globs.DimGP);
    Sub2Vectors(vars[NET_ATM].hybrid, vars[NET_ATM].hybrid, vars[176].hybrid, globs.DimGP);
    Sub2Vectors(vars[NET_ATM].hybrid, vars[NET_ATM].hybrid, vars[177].hybrid, globs.DimGP);
  }

  if (vars[SURF_RUNOFF].comp)
  {
    vars[SURF_RUNOFF].hlev = vars[SURF_RUNOFF].plev = 1;
    vars[SURF_RUNOFF].sfit = false;
    vars[SURF_RUNOFF].hybrid = alloc_dp(globs.DimGP, "vars[SURF_RUNOFF].hybrid");
    Sub2Vectors(vars[SURF_RUNOFF].hybrid, vars[182].hybrid, vars[221].hybrid, globs.DimGP);
    Add2Vectors(vars[SURF_RUNOFF].hybrid, vars[SURF_RUNOFF].hybrid, vars[142].hybrid, globs.DimGP);
    Add2Vectors(vars[SURF_RUNOFF].hybrid, vars[SURF_RUNOFF].hybrid, vars[143].hybrid, globs.DimGP);
  }

  if (vars[FRESH_WATER].comp)
  {
    vars[FRESH_WATER].hlev = vars[FRESH_WATER].plev = 1;
    vars[FRESH_WATER].sfit = false;
    vars[FRESH_WATER].hybrid = alloc_dp(globs.DimGP, "vars[FRESH_WATER].hybrid");
    Add2Vectors(vars[FRESH_WATER].hybrid, vars[142].hybrid, vars[143].hybrid, globs.DimGP);
    Add2Vectors(vars[FRESH_WATER].hybrid, vars[FRESH_WATER].hybrid, vars[182].hybrid, globs.DimGP);
  }
}

static void
Derivate(double field[], double derilam[], long nlevels, long Waves, long Latitudes, Varray<double> const &derivationFactor)
{
  long i = 0;
  for (long lev = 0; lev < nlevels; lev++)
    for (long n = 0; n < Waves; ++n)
    {
      for (long l = 0; l < Latitudes; ++l)
      {
        derilam[i] = -n * field[i + Latitudes] * derivationFactor[l];
        i++;
      }
      for (long l = 0; l < Latitudes; ++l)
      {
        derilam[i] = n * field[i - Latitudes] * derivationFactor[l];
        i++;
      }
    }
}

// Process Model Level data
void
after_processML(AfterControl &globs, struct Variable *vars)
{
  int leveltype;
  size_t numMissVals;
  double *pressureLevel = nullptr;

  globs.MeanCount++;
  globs.TermCount++;

  if (globs.MeanCount == 1)
  {
    if (globs.Debug) cdo_print("TermCount = %d MeanCount = %d", globs.TermCount, globs.MeanCount);
    globs.StartDate = globs.OldDate;
  }

  if (globs.TermCount > 120) globs.Debug = 0;

  /* ============================== */
  /* Computations in spectral space */
  /* ============================== */

  if (vars[U_WIND].comp || vars[V_WIND].comp)
  {
    vars[U_WIND].hlev = vars[V_WIND].hlev = vars[DIVERGENCE].hlev;
    vars[U_WIND].plev = vars[V_WIND].plev = vars[DIVERGENCE].plev;
    vars[U_WIND].sfit = vars[V_WIND].sfit = true;
    vars[U_WIND].spectral = alloc_dp(globs.Dim3SP, "vars[U_WIND].spectral");
    vars[V_WIND].spectral = alloc_dp(globs.Dim3SP, "vars[V_WIND].spectral");

    if (vars[DIVERGENCE].spectral == nullptr) after_gp2sp(globs, vars, DIVERGENCE);
    if (vars[VORTICITY].spectral == nullptr) after_gp2sp(globs, vars, VORTICITY);

    dv2uv(vars[DIVERGENCE].spectral, vars[VORTICITY].spectral, vars[U_WIND].spectral, vars[V_WIND].spectral, globs.dv2uv_f1.data(),
          globs.dv2uv_f2.data(), globs.Truncation, globs.DimSP, vars[DIVERGENCE].hlev);
  }

  if (vars[VELOPOT].comp && globs.Type < 30)
  {
    vars[VELOPOT].hlev = vars[DIVERGENCE].hlev;
    vars[VELOPOT].plev = vars[DIVERGENCE].plev;
    vars[VELOPOT].spectral = alloc_dp(globs.Dim3SP, "vars[VELOPOT].spectral");

    if (vars[DIVERGENCE].spectral == nullptr) after_gp2sp(globs, vars, DIVERGENCE);

    dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, vars[DIVERGENCE].hlev, globs.Truncation);
  }

  if (vars[STREAM].comp && globs.Type < 30)
  {
    vars[STREAM].hlev = vars[VORTICITY].hlev;
    vars[STREAM].plev = vars[VORTICITY].plev;
    vars[STREAM].spectral = alloc_dp(globs.Dim3SP, "vars[STREAM].spectral");

    if (vars[VORTICITY].spectral == nullptr) after_gp2sp(globs, vars, VORTICITY);

    dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral, vars[VORTICITY].hlev, globs.Truncation);
  }

  if (!vars[VORTICITY].needed) FreeMemory(vars[VORTICITY].spectral);
  if (!vars[DIVERGENCE].needed) FreeMemory(vars[DIVERGENCE].spectral);

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if (globs.Type == 0)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].spectral) afterAbort("Code %d not available on spectral space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimSP;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].spectral + offset, 0);
        }
      }

    FreeSpectral(vars);
    return;
  }

  /* ------------------------------- */
  /*  Computations in fourier space  */
  /* ------------------------------- */

  if (globs.Type >= 10)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].spectral)
      {
        if (vars[code].fourier == nullptr)
        {
          auto fieldSize = vars[code].hlev * globs.DimFC;
          vars[code].fourier = alloc_dp(fieldSize, FieldName(code, "fourier"));
          sp2fc(vars[code].spectral, vars[code].fourier, globs.poli.data(), vars[code].hlev, globs.Latitudes, globs.Fouriers,
                globs.Truncation);
        }
        if (code != LNPS) FreeMemory(vars[code].spectral);
      }

    if (vars[U_WIND].needed && vars[U_WIND].fourier)
      scaluv(vars[U_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevel);
    if (vars[V_WIND].needed && vars[V_WIND].fourier)
      scaluv(vars[V_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevel);

    if (vars[DPSDX].needed)
    {
      vars[DPSDX].hlev = 1;
      vars[DPSDX].plev = 1;
      vars[DPSDX].sfit = false;
      vars[DPSDX].fourier = alloc_dp(globs.DimFC, "vars[DPSDX].fourier");
      if (vars[LNPS].fourier == nullptr) after_gp2sp(globs, vars, LNPS);
      Derivate(vars[LNPS].fourier, vars[DPSDX].fourier, 1, globs.Waves, globs.Latitudes, globs.derivationFactor);
    }
    if (vars[DPSDY].needed)
    {
      vars[DPSDY].hlev = 1;
      vars[DPSDY].plev = 1;
      vars[DPSDY].sfit = false;
      vars[DPSDY].fourier = alloc_dp(globs.DimFC, "vars[DPSDY].fourier");
      if (vars[LNPS].spectral == nullptr) after_gp2sp(globs, vars, LNPS);
      sp2fc(vars[LNPS].spectral, vars[DPSDY].fourier, globs.pdev.data(), vars[DPSDY].hlev, globs.Latitudes, globs.Fouriers,
            globs.Truncation);
    }
  }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if (globs.Type == 10)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on fourier space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, 0);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if (globs.Type == 11)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on zonal mean!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, 0);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ------------------------------ */
  /*  Transformation to gridpoints  */
  /* ------------------------------ */

  if (globs.Type >= 20)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].fourier)
      {
        if (vars[code].hybrid == nullptr)
        {
          auto fieldSize = globs.DimGP * vars[code].hlev;
          vars[code].hybrid = alloc_dp(fieldSize, FieldName(code, "hybrid"));
          after_FC2GP(vars[code].fourier, vars[code].hybrid, globs.Latitudes, globs.Longitudes, vars[code].hlev, globs.Fouriers);
        }
        FreeMemory(vars[code].fourier);
      }

    if (vars[PS_PROG].comp && vars[PS_PROG].hybrid == nullptr)
    {
      vars[PS_PROG].hybrid = alloc_dp(globs.DimGP, "PS_PROG");
      if (vars[LNPS].hybrid)
      {
        for (int l = 0; l < globs.DimGP; ++l) vars[PS_PROG].hybrid[l] = std::exp(vars[LNPS].hybrid[l]);
      }
      else if (vars[PS].hybrid)
      {
        cdo_warning("log surface pressure (code 152) not found - using surface pressure (code 134)!");
        array_copy(globs.DimGP, vars[PS].hybrid, vars[PS_PROG].hybrid);
      }
      else { afterAbort("surface pressure not found!"); }
    }
    vars[LNPS].needed = vars[LNPS].selected;

    if (globs.orography.size() == 0)
    {
      globs.orography.resize(globs.DimGP);
      if (vars[GEOPOTENTIAL].hybrid)
        array_copy(globs.DimGP, vars[GEOPOTENTIAL].hybrid, globs.orography.data());
      else
      {
        if (vars[GEOPOTENTIAL].selected || globs.Type >= 30)
        {
          cdo_warning("Orography not found - using zero orography!");
          std::ranges::fill(globs.orography, 0.0);
        }
      }
    }
    vars[GEOPOTENTIAL].needed = vars[GEOPOTENTIAL].selected;

    after_EchamCompGP(globs, vars);
  }

  FreeFourier(vars);

  if (globs.Type == 20)
  {
    /* ----------------------- */
    /*  Means on hybrid grids  */
    /* ----------------------- */

    if (globs.Mean)
    {
      for (int code = 0; code < MaxCodes; ++code)
      {
        if (vars[code].selected && vars[code].hybrid)
        {
          auto fieldSize = globs.DimGP * vars[code].hlev;

          if (vars[code].mean == nullptr) vars[code].mean = alloc_dp(fieldSize, FieldName(code, "mean"));

          if (globs.Mean > 1 && vars[code].variance == nullptr)
            vars[code].variance = alloc_dp(fieldSize, FieldName(code, "variance"));

          if (globs.MeanCount == 1)
          {
            array_copy(fieldSize, vars[code].hybrid, vars[code].mean);
            if (globs.Mean > 1) IniQuaSum(vars[code].variance, vars[code].hybrid, fieldSize);
          }
          else
          {
            AddVector(vars[code].mean, vars[code].hybrid, fieldSize, &vars[code].numMissVals, vars[code].missval);
            if (globs.Mean > 1) AddQuaSum(vars[code].variance, vars[code].hybrid, fieldSize);
          }

          if (globs.EndOfInterval)
          {
            if (vars[code].samp == nullptr)
              MultVectorScalar(vars[code].hybrid, vars[code].mean, 1.0 / globs.MeanCount, fieldSize, vars[code].numMissVals,
                               vars[code].missval);
            else
              DivVectorIvector(vars[code].hybrid, vars[code].mean, vars[code].samp, fieldSize, &vars[code].numMissVals,
                               vars[code].missval);

            if (globs.Mean > 1) VarQuaSum(vars[code].variance, vars[code].hybrid, fieldSize, globs.MeanCount);
          }
        }
      }
    }

    /* ---------------------------- */
    /* Output of hybrid level grids */
    /* ---------------------------- */

    if (globs.Mean == 0 || globs.EndOfInterval)
    {
      for (int code = 0; code < MaxCodes; ++code)
        if (vars[code].selected)
        {
          if (vars[code].hybrid == nullptr) afterAbort("Internal problem. Code %d not allocated!", code);

          int nlevels = zaxisInqSize(vars[code].ozaxisID);
          for (int k = 0; k < nlevels; ++k)
          {
            auto offset = k * globs.DimGP;
            if (globs.Mean != 2)
            {
              streamDefField(globs.ostreamID, vars[code].ovarID, k);
              streamWriteField(globs.ostreamID, vars[code].hybrid + offset, vars[code].numMissVals);
            }
            if (globs.Mean >= 2)
            {
              streamDefField(globs.ostreamID2, vars[code].ovarID2, k);
              streamWriteField(globs.ostreamID2, vars[code].variance + offset, vars[code].numMissVals);
            }
          }
        }

      FreeSamp(vars);
    }

    FreeHybrid(vars);
    return;
  }

  /* -------------------------------------- */
  /* Vertical interpolation / extrapolation */
  /* -------------------------------------- */

  if (globs.Type >= 30)
  {
    if (globs.vertIndex.size() == 0) globs.vertIndex.resize(globs.NumLevelRequest * globs.DimGP);

    if (globs.unitsel)
    {
      if (globs.p_of_height.size() == 0) globs.p_of_height.resize(globs.NumLevelRequest);
      height_to_pressure(globs.LevelRequest, globs.p_of_height.data(), globs.NumLevelRequest);
      pressureLevel = globs.p_of_height.data();
    }
    else { pressureLevel = globs.LevelRequest; }

    gen_vert_index(globs.vertIndex.data(), pressureLevel, vars[FULL_PRESS].hybrid, globs.DimGP, globs.NumLevelRequest,
                   globs.NumLevel);

    numMissVals = 0;
    if (!globs.extrapolate)
    {
      Varray<size_t> pnumMissVals(globs.NumLevelRequest);
      gen_vert_index_mv(globs.vertIndex.data(), pressureLevel, globs.DimGP, globs.NumLevelRequest, vars[PS_PROG].hybrid,
                        pnumMissVals.data());
      for (int i = 0; i < globs.NumLevelRequest; ++i) numMissVals += pnumMissVals[i];
    }

    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].hybrid)
      {
        leveltype = zaxisInqType(vars[code].izaxisID);
        if (vars[code].hlev == 1 || leveltype != ZAXIS_HYBRID || (vars[code].hlev < globs.NumLevel))
        {
          FreeMemory(vars[code].grid);
          vars[code].grid = vars[code].hybrid;
          vars[code].hybrid = nullptr;
        }
        else
        {
          if (vars[code].grid == nullptr)
          {
            auto fieldSize = globs.DimGP * globs.NumLevelRequest;
            vars[code].grid = alloc_dp(fieldSize, FieldName(code, "grid"));
          }

          if (code == TEMPERATURE)
          {
            vertical_interp_T(globs.orography.data(), vars[TEMPERATURE].hybrid, vars[TEMPERATURE].grid, vars[FULL_PRESS].hybrid,
                              vars[HALF_PRESS].hybrid, globs.vertIndex.data(), pressureLevel, globs.NumLevelRequest, globs.DimGP,
                              globs.NumLevel, vars[code].missval);
          }
          else if (code == GEOPOTHEIGHT)
          {
            if (vars[TEMPERATURE].hybrid == nullptr) afterAbort("Code  130 not found!");
            vertical_interp_Z(globs.orography.data(), vars[GEOPOTHEIGHT].hybrid, vars[GEOPOTHEIGHT].grid, vars[FULL_PRESS].hybrid,
                              vars[HALF_PRESS].hybrid, globs.vertIndex.data(), vars[TEMPERATURE].hybrid, pressureLevel,
                              globs.NumLevelRequest, globs.DimGP, globs.NumLevel, vars[code].missval);
          }
          else
          {
            int numlevel = vars[code].hlev;
            if (code == OMEGA) numlevel = globs.NumLevel;
            double *hyb_press = vars[FULL_PRESS].hybrid;
            if (numlevel == (globs.NumLevel + 1)) hyb_press = vars[HALF_PRESS].hybrid;

            vertical_interp_X(vars[code].hybrid, vars[code].grid, hyb_press, globs.vertIndex.data(), pressureLevel,
                              globs.NumLevelRequest, globs.DimGP, numlevel, vars[code].missval);
          }

          if (!globs.extrapolate) vars[code].numMissVals = numMissVals;

          if (code != TEMPERATURE) FreeMemory(vars[code].hybrid);
        }
      }
  }

  vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
  FreeHybrid(vars);
  FreeMemory(vars[HALF_PRESS].hybrid);

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if (globs.Type == 30 && globs.Mean == 0)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected && vars[code].grid)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          if (globs.Mean != 2)
          {
            auto offset = k * globs.DimGP;
            streamDefField(globs.ostreamID, vars[code].ovarID, k);
            streamWriteField(globs.ostreamID, vars[code].grid + offset, vars[code].numMissVals);
          }
        }
      }

    FreeGrid(vars);
    return;
  }

  /* ---------------------- */
  /*  Computation of Means  */
  /* ---------------------- */

  if (globs.Type >= 30 && globs.Mean)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].grid)
      {
        auto fieldSize = globs.DimGP * vars[code].plev;

        if (vars[code].mean == nullptr) vars[code].mean = alloc_dp(fieldSize, FieldName(code, "mean"));

        if (globs.MeanCount == 1)
          array_copy(fieldSize, vars[code].grid, vars[code].mean);
        else
          AddVector(vars[code].mean, vars[code].grid, fieldSize, &vars[code].numMissVals, vars[code].missval);

        if (globs.EndOfInterval)
        {
          if (vars[code].samp == nullptr)
            MultVectorScalar(vars[code].mean, vars[code].mean, 1.0 / globs.MeanCount, fieldSize, vars[code].numMissVals,
                             vars[code].missval);
          else
            DivVectorIvector(vars[code].mean, vars[code].mean, vars[code].samp, fieldSize, &vars[code].numMissVals,
                             vars[code].missval);
        }
      }

  /* -------------------------- */
  /*  Computation of Variances  */
  /* -------------------------- */

  if (globs.Type >= 30 && globs.Mean > 1)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].mean)
      {
        auto fieldSize = globs.DimGP * vars[code].plev;

        if (vars[code].variance == nullptr) vars[code].variance = alloc_dp(fieldSize, FieldName(code, "var"));

        if (globs.MeanCount == 1)
          IniQuaSum(vars[code].variance, vars[code].grid, fieldSize);
        else
          AddQuaSum(vars[code].variance, vars[code].grid, fieldSize);

        if (globs.EndOfInterval) VarQuaSum(vars[code].variance, vars[code].mean, fieldSize, globs.MeanCount);
      }

  if (globs.Mean && !globs.EndOfInterval)
  {
    FreeGrid(vars);
    return;
  }

  /* --------------------------------------------- */
  /*  Output of pressure level means and variances */
  /* --------------------------------------------- */

  if (globs.Type == 30 && globs.Mean)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected && vars[code].mean)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimGP;
          if (globs.Mean != 2)
          {
            streamDefField(globs.ostreamID, vars[code].ovarID, k);
            streamWriteField(globs.ostreamID, vars[code].mean + offset, vars[code].numMissVals);
          }
          if (globs.Mean >= 2)
          {
            streamDefField(globs.ostreamID2, vars[code].ovarID2, k);
            streamWriteField(globs.ostreamID2, vars[code].variance + offset, vars[code].numMissVals);
          }
        }
      }

    FreeSamp(vars);
    FreeGrid(vars);
    return;
  }

  /* ------------------ */
  /*  Free mean fields  */
  /* ------------------ */

  if (globs.Type >= 40 && globs.Mean)
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].mean)
      {
        FreeMemory(vars[code].variance);
        FreeMemory(vars[code].grid);
        vars[code].grid = vars[code].mean;
        vars[code].mean = nullptr;
      }

  /* --------------------------------- */
  /*  Transformation to fourier space  */
  /* --------------------------------- */

  if (globs.Type >= 40)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].grid && (vars[code].sfit || globs.Type < 70))
      {
        if (vars[code].numMissVals) afterAbort("Missing values for code %d unsupported with TYPE > 30!", code);

        if (vars[code].fourier == nullptr)
        {
          auto fieldSize = globs.DimFC * vars[code].plev;
          vars[code].fourier = alloc_dp(fieldSize, FieldName(code, "fourier"));
        }

        after_GP2FC(vars[code].grid, vars[code].fourier, globs.Latitudes, globs.Longitudes, vars[code].plev, globs.Fouriers);

        if (vars[code].sfit || globs.Type < 70) FreeMemory(vars[code].grid);
      }
  }

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].sfit || globs.Type < 70) FreeMemory(vars[code].grid);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if (globs.Type == 40)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on fourier space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, vars[code].numMissVals);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* --------------------- */
  /* Output of zonal means */
  /* --------------------- */

  if (globs.Type == 41)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on zonal mean!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, vars[code].numMissVals);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ---------------------------------- */
  /*  Transformation to spectral space  */
  /* ---------------------------------- */

  if (globs.Type >= 50)
  {
    if (vars[U_WIND].needed && vars[U_WIND].fourier)
      scaluv(vars[U_WIND].fourier, globs.coslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);
    if (vars[V_WIND].needed && vars[V_WIND].fourier)
      scaluv(vars[V_WIND].fourier, globs.coslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);

    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].fourier)
      {
        if (vars[code].spectral == nullptr)
        {
          auto fieldSize = vars[code].plev * globs.DimSP;
          vars[code].spectral = alloc_dp(fieldSize, FieldName(code, "spectral"));
        }

        fc2sp(vars[code].fourier, vars[code].spectral, globs.pold.data(), vars[code].plev, globs.Latitudes, globs.Fouriers,
              globs.Truncation);
      }

    if (vars[DIVERGENCE].needed || vars[VORTICITY].needed || vars[VELOPOT].needed || vars[STREAM].needed)
    {
      if (vars[DIVERGENCE].spectral == nullptr)
        vars[DIVERGENCE].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[DIVERGENCE].spectral");
      if (vars[VORTICITY].spectral == nullptr)
        vars[VORTICITY].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[VORTICITY].spectral");
      if ((vars[U_WIND].fourier == 0 || vars[V_WIND].fourier == 0) && globs.NumLevelRequest) afterAbort("uwind or vwind missing!");
      uv2dv(vars[U_WIND].fourier, vars[V_WIND].fourier, vars[DIVERGENCE].spectral, vars[VORTICITY].spectral, globs.pol2.data(),
            globs.pol3.data(), globs.NumLevelRequest, globs.Latitudes, globs.Truncation);
    }

    if (vars[VELOPOT].needed)
    {
      vars[VELOPOT].hlev = vars[DIVERGENCE].hlev;
      vars[VELOPOT].plev = vars[DIVERGENCE].plev;
      vars[VELOPOT].sfit = true;
      if (vars[VELOPOT].spectral == nullptr)
        vars[VELOPOT].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[VELOPOT].spectral");
      dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, globs.NumLevelRequest, globs.Truncation);
    }

    if (vars[STREAM].needed)
    {
      vars[STREAM].hlev = vars[VORTICITY].hlev;
      vars[STREAM].plev = vars[VORTICITY].plev;
      vars[STREAM].sfit = true;
      if (vars[STREAM].spectral == nullptr)
        vars[STREAM].spectral = alloc_dp(globs.DimSP * globs.NumLevelRequest, "vars[STREAM].spectral");
      dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral, globs.NumLevelRequest, globs.Truncation);
    }
  }

  for (int code = 0; code < MaxCodes; ++code)
    if (vars[code].sfit || globs.Type < 61) FreeMemory(vars[code].fourier);

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if (globs.Type == 50)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected && vars[code].spectral)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimSP;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].spectral + offset, 0);
        }
      }

    FreeSpectral(vars);
    return;
  }

  /* -------------------------------*/
  /*  Computations in fourier space */
  /* -------------------------------*/

  if (globs.Type >= 60)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].spectral)
      {
        if (vars[code].fourier == nullptr)
        {
          auto fieldSize = vars[code].plev * globs.DimFC;
          vars[code].fourier = alloc_dp(fieldSize, FieldName(code, "fourier"));
        }
        sp2fc(vars[code].spectral, vars[code].fourier, globs.poli.data(), vars[code].plev, globs.Latitudes, globs.Fouriers,
              globs.Truncation);
      }
    if (vars[U_WIND].needed && vars[U_WIND].fourier)
      scaluv(vars[U_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);
    if (vars[V_WIND].needed && vars[V_WIND].fourier)
      scaluv(vars[V_WIND].fourier, globs.rcoslat.data(), globs.Latitudes, globs.Fouriers * globs.NumLevelRequest);
  }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if (globs.Type == 60)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on fourier space!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, 0);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if (globs.Type == 61)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        if (!vars[code].fourier) afterAbort("Code %d not available on zonal mean!", code);

        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimFC;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].fourier + offset, vars[code].numMissVals);
        }
      }

    FreeFourier(vars);
    return;
  }

  /* ------------------------------ */
  /*  Transformation to gridpoints  */
  /* ------------------------------ */

  if (globs.Type >= 70)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].needed && vars[code].fourier)
      {
        auto fieldSize = vars[code].plev * globs.DimGP;
        if (vars[code].grid == nullptr) vars[code].grid = alloc_dp(fieldSize, FieldName(code, "grid"));

        after_FC2GP(vars[code].fourier, vars[code].grid, globs.Latitudes, globs.Longitudes, vars[code].plev, globs.Fouriers);
      }
  }

  FreeFourier(vars);

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if (globs.Type == 70)
  {
    for (int code = 0; code < MaxCodes; ++code)
      if (vars[code].selected)
      {
        int nlevels = zaxisInqSize(vars[code].ozaxisID);
        for (int k = 0; k < nlevels; ++k)
        {
          auto offset = k * globs.DimGP;
          streamDefField(globs.ostreamID, vars[code].ovarID, k);
          streamWriteField(globs.ostreamID, vars[code].grid + offset, vars[code].numMissVals);
        }
      }

    FreeSamp(vars);
    FreeGrid(vars);
    return;
  }
}

void
after_AnalysisAddField(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                       size_t numMissVals)
{
  int truncation;

  auto gridtype = gridInqType(gridID);
  auto leveltype = zaxisInqType(zaxisID);
  auto nlevel = zaxisInqSize(zaxisID);
  auto gridSize = gridInqSize(gridID);
  auto dataSize = gridSize * nlevel;
  auto dataOffset = gridSize * levelID;

  vars[code].numMissVals0 += numMissVals;

  if (gridtype == GRID_SPECTRAL)
  {
    vars[code].sfit = true;
    vars[code].hlev = globs->NumLevelRequest;
    vars[code].plev = globs->NumLevelRequest;
    if (nlevel > 1 && leveltype == ZAXIS_PRESSURE)
    {
      if (code != U_WIND && code != V_WIND)
      {
        auto fieldSize = globs->Dim3SP;
        if (vars[code].spectral0 == nullptr) vars[code].spectral0 = alloc_dp(fieldSize, FieldName(code, "spectral"));
        truncation = gridInqTrunc(gridID);
        sp2sp(globs->varray.data(), truncation, vars[code].spectral0 + levelID * globs->DimSP, globs->Truncation);
      }
      else
      {
        auto fieldSize = globs->Dim3SP;
        if (vars[code].spectral0 == nullptr) vars[code].spectral0 = alloc_dp(fieldSize, FieldName(code, "spectral"));
        array_copy(globs->DimSP, globs->varray.data(), vars[code].spectral0 + levelID * globs->DimSP);
      }
    }
    else
      afterAbort("Only pressure level data supported for spectral data!");
  }
  else
  {
    if (nlevel > 1 && leveltype == ZAXIS_PRESSURE)
    {
      vars[code].sfit = true;
      auto fieldSize = globs->Dim3GP;
      vars[code].hlev = globs->NumLevelRequest;
      vars[code].plev = globs->NumLevelRequest;
      if (vars[code].grid0 == nullptr) vars[code].grid0 = alloc_dp(fieldSize, FieldName(code, "grid0"));
      array_copy(globs->DimGP, globs->varray.data(), vars[code].grid0 + levelID * globs->DimGP);
    }
    else
    {
      vars[code].sfit = false;
      auto fieldSize = globs->DimGP;
      vars[code].hlev = 1;
      vars[code].plev = 1;
      if (vars[code].grid0 == nullptr) vars[code].grid0 = alloc_dp(fieldSize, FieldName(code, "grid0"));
      array_copy(globs->DimGP, globs->varray.data(), vars[code].grid0);
    }

    if (globs->Mean > 0 && (numMissVals || vars[code].samp))
    {
      if (vars[code].samp == nullptr)
      {
        vars[code].samp = (int *) std::malloc(dataSize * sizeof(int));
        for (size_t i = 0; i < dataSize; ++i) vars[code].samp[i] = globs->MeanCount0;
      }

      for (size_t i = 0; i < gridSize; ++i)
        if (is_not_equal(globs->varray[i], vars[code].missval)) vars[code].samp[i + dataOffset]++;
    }
  }
}

double *
after_get_dataptr(struct Variable *vars, int code, int gridID, int zaxisID, int levelID)
{
  auto gridtype = gridInqType(gridID);
  auto nlevel = zaxisInqSize(zaxisID);
  auto gridSize = gridInqSize(gridID);
  auto dataSize = gridSize * nlevel;
  auto dataOffset = gridSize * levelID;

  double *dataptr = nullptr;

  if (gridtype == GRID_SPECTRAL)
  {
    if (vars[code].spectral0 == nullptr) vars[code].spectral0 = alloc_dp(dataSize, FieldName(code, "spectral0"));

    dataptr = vars[code].spectral0 + dataOffset;
  }
  else
  {
    if (vars[code].hybrid0 == nullptr) vars[code].hybrid0 = alloc_dp(dataSize, FieldName(code, "hybrid0"));

    dataptr = vars[code].hybrid0 + dataOffset;
  }

  return dataptr;
}

void
after_EchamAddField(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                    size_t numMissVals)
{
  auto gridtype = gridInqType(gridID);
  auto leveltype = zaxisInqType(zaxisID);
  auto nlevel = zaxisInqSize(zaxisID);
  auto gridSize = gridInqSize(gridID);
  auto dataSize = gridSize * nlevel;
  auto dataOffset = gridSize * levelID;

  vars[code].numMissVals0 += numMissVals;

  if (gridtype == GRID_SPECTRAL)
  {
    vars[code].sfit = true;
    vars[code].hlev = nlevel;
    vars[code].plev = 1;
    if (nlevel > 1 && leveltype == ZAXIS_HYBRID) vars[code].plev = globs->NumLevelRequest;

    if (gridInqTrunc(gridID) != globs->Truncation)
      afterAbort("Resolution error. Required %d - Found %d", globs->Truncation, gridInqTrunc(gridID));
  }
  else
  {
    vars[code].hlev = nlevel;
    vars[code].plev = nlevel;
    vars[code].sfit = false;
    if (nlevel > 1 && leveltype == ZAXIS_HYBRID && (nlevel == globs->NumLevel || nlevel == globs->NumLevel + 1))
    {
      vars[code].plev = globs->NumLevelRequest;
      vars[code].sfit = true;
    }

    if (globs->Mean > 0 && (numMissVals || vars[code].samp))
    {
      if (vars[code].samp == nullptr)
      {
        vars[code].samp = (int *) std::malloc(dataSize * sizeof(int));
        for (size_t i = 0; i < dataSize; ++i) vars[code].samp[i] = globs->MeanCount0;
      }

      double *dataptr = vars[code].hybrid0 + dataOffset;
      for (size_t i = 0; i < gridSize; ++i)
        if (is_not_equal(dataptr[i], vars[code].missval)) vars[code].samp[i + dataOffset]++;
    }
  }
}

static void
MakeDependencies(struct Variable *vars, int varcode, int depcode)
{
  if (vars[varcode].needed && !vars[varcode].detected)
  {
    vars[depcode].needed = true;
    vars[varcode].comp = true;

    if (afterDebug) fprintf(stderr, "Needed code %d to compute code %d\n", depcode, varcode);

    if (vars[depcode].ivarID == -1)
    {
      if (depcode == U_WIND)
      {
        MakeDependencies(vars, U_WIND, DIVERGENCE);
        MakeDependencies(vars, U_WIND, VORTICITY);
      }
      if (depcode == V_WIND)
      {
        MakeDependencies(vars, V_WIND, DIVERGENCE);
        MakeDependencies(vars, V_WIND, VORTICITY);
      }
    }

    if (vars[varcode].ivarID == -1)
    {
      if (vars[depcode].ivarID == -1) { afterAbort("code %d undefined, needed to compute code %d", depcode, varcode); }
      else
      {
        vars[varcode].ivarID = vars[depcode].ivarID;
        vars[varcode].igridID = vars[depcode].igridID;
        vars[varcode].ogridID = vars[depcode].ogridID;
        vars[varcode].izaxisID = vars[depcode].izaxisID;
        if (varcode != GEOPOTHEIGHT) vars[varcode].ozaxisID = vars[depcode].ozaxisID;
      }
    }
  }
}

static void
CheckDependencies(struct Variable *vars, int analysisdata)
{
  MakeDependencies(vars, VELOPOT, U_WIND);
  MakeDependencies(vars, VELOPOT, V_WIND);
  MakeDependencies(vars, VELOPOT, VORTICITY);
  MakeDependencies(vars, VELOPOT, DIVERGENCE);

  MakeDependencies(vars, STREAM, U_WIND);
  MakeDependencies(vars, STREAM, V_WIND);
  MakeDependencies(vars, STREAM, VORTICITY);
  MakeDependencies(vars, STREAM, DIVERGENCE);

  MakeDependencies(vars, VORTICITY, U_WIND);
  MakeDependencies(vars, VORTICITY, V_WIND);

  MakeDependencies(vars, DIVERGENCE, U_WIND);
  MakeDependencies(vars, DIVERGENCE, V_WIND);

  MakeDependencies(vars, U_WIND, VORTICITY);
  MakeDependencies(vars, U_WIND, DIVERGENCE);
  MakeDependencies(vars, U_WIND, V_WIND);

  MakeDependencies(vars, V_WIND, VORTICITY);
  MakeDependencies(vars, V_WIND, DIVERGENCE);
  MakeDependencies(vars, V_WIND, U_WIND);

  MakeDependencies(vars, WINDSPEED, U_WIND);
  MakeDependencies(vars, WINDSPEED, V_WIND);

  if (analysisdata)
  {
    MakeDependencies(vars, RHUMIDITY, HUMIDITY);
    MakeDependencies(vars, RHUMIDITY, TEMPERATURE);
    MakeDependencies(vars, HUMIDITY, RHUMIDITY);
    MakeDependencies(vars, HUMIDITY, TEMPERATURE);
    MakeDependencies(vars, GEOPOTHEIGHT, GEOPOTENTIAL);
  }
  else
  {
    MakeDependencies(vars, THETAF, TEMPERATURE);
    MakeDependencies(vars, SLP, TEMPERATURE);
  }

  MakeDependencies(vars, SW_ATM, 176);
  MakeDependencies(vars, SW_ATM, 178);
  MakeDependencies(vars, LW_ATM, 177);
  MakeDependencies(vars, LW_ATM, 179);
  MakeDependencies(vars, NET_ATM, 176);
  MakeDependencies(vars, NET_ATM, 177);
  MakeDependencies(vars, NET_ATM, 178);
  MakeDependencies(vars, NET_ATM, 179);

  MakeDependencies(vars, SW_BOT_CLF, 176);
  MakeDependencies(vars, SW_BOT_CLF, 185);
  MakeDependencies(vars, LW_BOT_CLF, 177);
  MakeDependencies(vars, LW_BOT_CLF, 186);
  MakeDependencies(vars, SW_TOP_CLF, 178);
  MakeDependencies(vars, SW_TOP_CLF, 187);
  MakeDependencies(vars, LW_TOP_CLF, 179);
  MakeDependencies(vars, LW_TOP_CLF, 188);

  MakeDependencies(vars, ALL_WATER, 222);
  MakeDependencies(vars, LOW_WATER, 222);
  MakeDependencies(vars, MID_WATER, 222);
  MakeDependencies(vars, HIH_WATER, 222);

  MakeDependencies(vars, LOW_CLOUD, 223);
  MakeDependencies(vars, MID_CLOUD, 223);
  MakeDependencies(vars, HIH_CLOUD, 223);

  if (vars[LOW_CLOUD].comp || vars[MID_CLOUD].comp || vars[HIH_CLOUD].comp)
  {
    static int zaxisID = -999;
    if (zaxisID == -999) zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

    vars[LOW_CLOUD].izaxisID = zaxisID;
    vars[LOW_CLOUD].ozaxisID = zaxisID;
    vars[MID_CLOUD].izaxisID = zaxisID;
    vars[MID_CLOUD].ozaxisID = zaxisID;
    vars[HIH_CLOUD].izaxisID = zaxisID;
    vars[HIH_CLOUD].ozaxisID = zaxisID;
  }
}

void
after_AnalysisDependencies(struct Variable *vars, int ncodes)
{
  for (int code = 0; code < ncodes; ++code) vars[code].needed = vars[code].selected;

  MakeDependencies(vars, PS, LNPS);

  CheckDependencies(vars, 1);
}

void
after_EchamDependencies(struct Variable *vars, int ncodes, int type, int source)
{
  for (int code = 0; code < ncodes; ++code) vars[code].needed = vars[code].selected;

  for (int code = 0; code < ncodes; ++code)
    if (vars[code].detected == false) vars[code].ivarID = -1;

  if (type >= 50)
  {
    vars[U_WIND].needed |= vars[DIVERGENCE].needed;
    vars[V_WIND].needed |= vars[DIVERGENCE].needed;
    vars[U_WIND].needed |= vars[VORTICITY].needed;
    vars[V_WIND].needed |= vars[VORTICITY].needed;
    vars[U_WIND].needed |= vars[VELOPOT].needed;
    vars[V_WIND].needed |= vars[VELOPOT].needed;
    vars[U_WIND].needed |= vars[STREAM].needed;
    vars[V_WIND].needed |= vars[STREAM].needed;
  }

  if (type >= 30) vars[LNPS].needed = true;

  if (type >= 20)
  {
    MakeDependencies(vars, THETAF, LNPS);
    MakeDependencies(vars, SLP, LNPS);
    MakeDependencies(vars, SLP, GEOPOTENTIAL);
    /*
    MakeDependencies(vars, SLP, HALF_PRESS);
    MakeDependencies(vars, SLP, FULL_PRESS);
    */
    MakeDependencies(vars, RHUMIDITY, TEMPERATURE);
    MakeDependencies(vars, RHUMIDITY, HUMIDITY);
    MakeDependencies(vars, RHUMIDITY, LNPS);
    MakeDependencies(vars, GEOPOTHEIGHT, TEMPERATURE);
    MakeDependencies(vars, GEOPOTHEIGHT, HUMIDITY);
    MakeDependencies(vars, GEOPOTHEIGHT, LNPS);
    MakeDependencies(vars, OMEGA, DIVERGENCE);
    MakeDependencies(vars, OMEGA, U_WIND);
    MakeDependencies(vars, OMEGA, V_WIND);
    MakeDependencies(vars, OMEGA, LNPS);
    MakeDependencies(vars, OMEGA, DPSDX);
    MakeDependencies(vars, OMEGA, DPSDY);
  }

  MakeDependencies(vars, DPSDX, LNPS);
  MakeDependencies(vars, HALF_PRESS, LNPS);
  MakeDependencies(vars, PS, LNPS);

  MakeDependencies(vars, SURF_RUNOFF, 142);
  MakeDependencies(vars, SURF_RUNOFF, 143);
  MakeDependencies(vars, SURF_RUNOFF, 182);
  MakeDependencies(vars, SURF_RUNOFF, 221); /* snow depth change */

  MakeDependencies(vars, THETAF, TS);

  MakeDependencies(vars, FRESH_WATER, 142);
  MakeDependencies(vars, FRESH_WATER, 143);
  MakeDependencies(vars, FRESH_WATER, 182);

  MakeDependencies(vars, PRECIP, 142);
  MakeDependencies(vars, PRECIP, 143);

  if (source != ECHAM5_Source)
  {
    MakeDependencies(vars, NET_WATER, 142);
    MakeDependencies(vars, NET_WATER, 143);
    MakeDependencies(vars, NET_WATER, 160);
    MakeDependencies(vars, NET_WATER, 182);

    MakeDependencies(vars, NET_TOP, 178);
    MakeDependencies(vars, NET_TOP, 179);

    MakeDependencies(vars, NET_BOT, 176);
    MakeDependencies(vars, NET_BOT, 177);

    MakeDependencies(vars, NET_HEAT, 146);
    MakeDependencies(vars, NET_HEAT, 147);
    MakeDependencies(vars, NET_HEAT, 176);
    MakeDependencies(vars, NET_HEAT, 177);
    MakeDependencies(vars, NET_HEAT, 218);
    /*
      if ( source == ECHAM5_Source )
      {
      MakeDependencies(vars, NET_HEAT, 206);
      MakeDependencies(vars, NET_HEAT, 208);
      MakeDependencies(vars, NET_HEAT, 209);
      }
      else
    */
    {
      MakeDependencies(vars, NET_HEAT, 220);
    }
  }

  MakeDependencies(vars, SW_CLF, 178);
  MakeDependencies(vars, SW_CLF, 224);

  MakeDependencies(vars, LW_CLF, 179);
  MakeDependencies(vars, LW_CLF, 225);

  MakeDependencies(vars, NET_CLF, 178);
  MakeDependencies(vars, NET_CLF, 179);
  MakeDependencies(vars, NET_CLF, 224);
  MakeDependencies(vars, NET_CLF, 225);

  if (vars[DPSDX].needed || vars[DPSDY].needed || vars[GEOPOTHEIGHT].comp || vars[SLP].comp || vars[THETAF].needed
      || vars[HALF_PRESS].needed || vars[RHUMIDITY].comp || vars[OMEGA].comp || type >= 30)
    vars[PS_PROG].comp = true;

  CheckDependencies(vars, 0);
}

void
after_legini_setup(AfterControl &globs, struct Variable *vars)
{
  long ntr = globs.Truncation;
  long nlat = globs.Latitudes;
  long dimsp = (ntr + 1) * (ntr + 2);
  long pdim = (dimsp / 2) * nlat;

  globs.poli.resize(pdim);

  if (!globs.AnalysisData)
  {
    if (globs.Type >= 20) globs.pold.resize(pdim);
    if (vars[DPSDY].needed) globs.pdev.resize(pdim);
  }

  if ((vars[DIVERGENCE].needed || vars[VORTICITY].needed || vars[VELOPOT].needed || vars[STREAM].needed) && globs.Type > 20)
  {
    globs.pol2.resize(pdim);
    globs.pol3.resize(pdim);
  }

  if (globs.AnalysisData && (globs.Type == 70) && globs.Gaussian && !globs.Spectral)
  {
    if (globs.poli.size()) varray_free(globs.poli);
    if (globs.pol2.size()) varray_free(globs.pol2);
    if (globs.pol3.size()) varray_free(globs.pol3);
    return;
  }

  after_legini_full(ntr, nlat, globs.poli.data(), globs.pold.data(), globs.pdev.data(), globs.pol2.data(), globs.pol3.data(),
                    globs.coslat.data());

  for (long jgl = 0; jgl < nlat; ++jgl) globs.rcoslat[jgl] = 1.0 / globs.coslat[jgl];

  for (long jgl = 0; jgl < nlat; ++jgl) globs.derivationFactor[jgl] = globs.rcoslat[jgl] / PlanetRadiusDefault;
}
