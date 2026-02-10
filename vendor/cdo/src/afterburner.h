#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "process_int.h"
#include "cdo_output.h"
#include "transform.h"
#include "workerthread.h"
#include "varray.h"

constexpr int MaxLevel = 1024;
constexpr int MaxCodes = 277;

constexpr int ECHAM5_Source = 1;

struct Date
{
  int yr;
  int mo;
  int dy;
  int hr;
  int mn;
};

struct AfterControl
{
  int Verbose = 0;

  int Mean = 0;
  int MeanCount0 = 0;
  int MeanCount = 0;
  int Multi = 0;
  int Nfiles = 0;
  int TermCount = 0;

  int OutputInterval = 0;
  int EndOfInterval = 0;

  int AnalysisData = 0;  // 0 = ECHAM Data, 1 = ECMWF Spectral Analyses
  int DayIn = 0;         // day increment of infiles if Multi = true
  int Debug = false;
  bool extrapolate = true;
  bool szip = false;

  int istreamID = CDI_UNDEFID;
  CdoStreamID ostreamID = CDO_STREAM_UNDEF;
  CdoStreamID ostreamID2 = CDO_STREAM_UNDEF;
  int ivlistID = CDI_UNDEFID;
  int ovlistID = CDI_UNDEFID;
  int ovlistID2 = CDI_UNDEFID;
  int taxisID = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;

  struct Date NextDate;
  struct Date NewDate;
  struct Date OldDate;
  struct Date StartDate;

  int numHalfLevels = 0;
  int nvct = 0;
  Varray<double> vct;

  Varray<int> vertIndex;
  Varray<double> p_of_height;
  Varray<double> orography;

  int Type = 0;
  int unitsel = 0;

  int Fouriers = 0;
  int Latitudes = 0;
  int Longitudes = 0;
  int HalfLevels = 0;
  int Gaussian = 0;
  int Spectral = 0;

  int Truncation = 0;
  int Waves = 0;

  long Dim3FC = 0, Dim3SP = 0, Dim3GP = 0;
  long DimFC = 0, DimGP = 0, DimSP = 0;
  long DimSP_half = 0;

  Varray<double> poli;
  Varray<double> pold;
  Varray<double> pdev;
  Varray<double> pol2;
  Varray<double> pol3;

  Varray<double> dv2uv_f1;
  Varray<double> dv2uv_f2;

  int NumCodesRequest = 0;

  int NumLevel = 0;
  int NumLevelFound = 0;
  int NumLevelRequest = 0;
  double LevelRequest[MaxLevel];

  Varray<double> rcoslat;
  Varray<double> coslat;
  Varray<double> derivationFactor;
  Varray<double> varray;
};

struct Variable
{
  int needed0;   // var needed for process
  int needed;    // var needed for process
  int selected;  // var selected for output
  int detected;  // var detected in input
  int comp;      // compute var if selected and not detected
  int sfit;
  int hlev;
  int plev;
  int ivarID;
  int ovarID;   // 1st variable ID
  int ovarID2;  // 2nd variable ID used for variance
  int tableID;
  int igridID;
  int ogridID;
  int izaxisID;
  int ozaxisID;
  size_t numMissVals0;
  size_t numMissVals;
  double missval;
  double *spectral;
  double *spectral0;
  double *fourier;
  double *hybrid;
  double *hybrid0;
  double *height;
  double *grid;
  double *grid0;
  double *mean;
  double *variance;
  int *samp;
};

// clang-format off
#define    LOW_CLOUD   34
#define    MID_CLOUD   35
#define    HIH_CLOUD   36
#define    LOW_WATER   37  // not used ?
#define    MID_WATER   38  // not used ?
#define    HIH_WATER   39  // not used ?
#define    ALL_WATER   40  // not used ?

#define GEOPOTENTIAL  129
#define  TEMPERATURE  130
#define       U_WIND  131
#define       V_WIND  132
#define     HUMIDITY  133
#define           PS  134
#define        OMEGA  135
#define    VORTICITY  138
#define           TS  139
#define       STREAM  148
#define      VELOPOT  149
#define          SLP  151
#define         LNPS  152
#define   DIVERGENCE  155
#define GEOPOTHEIGHT  156
#define    RHUMIDITY  157

#define   SW_BOT_CLF  189  // not used ?
#define   LW_BOT_CLF  190  // not used ?
#define   SW_TOP_CLF  191  // not used ?
#define   LW_TOP_CLF  192  // not used ?

#define    WINDSPEED  259
#define       PRECIP  260
#define      NET_TOP  261
#define      NET_BOT  262
#define     NET_HEAT  263
#define    NET_WATER  264
#define       SW_CLF  265
#define       LW_CLF  266
#define      NET_CLF  267
#define       SW_ATM  268
#define       LW_ATM  269
#define      NET_ATM  270
#define  SURF_RUNOFF  271
#define        DPSDX  273
#define        DPSDY  274
#define  FRESH_WATER  275
#define      PS_PROG  276  // PS for prognostic timestep
#define   HALF_PRESS  277
#define   FULL_PRESS  278
#define       THETAH  279
#define       THETAF  280
// clang-format on

void after_gp2sp(const AfterControl &globs, struct Variable *vars, int ccode);
void after_GP2FC(double *gp, double *fc, long nlat, long nlon, long nlev, long nfc);
void after_FC2GP(double *fc, double *gp, long nlat, long nlon, long nlev, long nfc);
void after_FCrh2FCsh(const AfterControl &globs, struct Variable *vars);
void after_SPuv2SPdv(const AfterControl &globs, struct Variable *vars);
void after_FCsh2FCrh(const AfterControl &globs, struct Variable *vars);

void after_EchamCompGP(const AfterControl &globs, struct Variable *vars);
void after_processPL(AfterControl &globs, struct Variable *vars);
void after_processML(AfterControl &globs, struct Variable *vars);

void after_AnalysisAddField(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                            size_t numMissVals);
double *after_get_dataptr(struct Variable *vars, int code, int gridID, int zaxisID, int levelID);
void after_EchamAddField(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                         size_t numMissVals);

void after_AnalysisDependencies(struct Variable *vars, int ncodes);
void after_EchamDependencies(struct Variable *vars, int ncodes, int type, int source);

void after_legini_setup(AfterControl &globs, struct Variable *vars);

static WorkerThread *afterWorkerThread = nullptr;
static bool afterReadAsync = true;

template <typename... Args>
void
afterAbort(std::string const &format, Args const &...args)
{
  if (afterReadAsync && afterWorkerThread) afterWorkerThread->wait();
  cdo_abort(format, args...);
}

#endif /* AFTERBURNER_H */
