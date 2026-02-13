/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "config.h"

#include <cdi.h>
#include "cdo_options.h"
#include "cdo_timer.h"
#include "util_string.h"
#include "cdo_output.h"

#include <cstring>
#include <algorithm>

namespace cdo
{
iTimer readTimer;
iTimer writeTimer;
}  // namespace cdo

namespace cdo
{
const char *progname;
const char *Version = "Climate Data Operators version " VERSION " (https://mpimet.mpg.de/cdo)";
std::string FileSuffix;
bool stdinIsTerminal = false;
bool stdoutIsTerminal = false;
bool stderrIsTerminal = false;
}  // namespace cdo

namespace Options
{
long coresize = 0;
int numStreamWorker = 0;
int nsb = 0;  // Number of significant bits
bool benchmark = false;
bool silentMode = false;
bool test = false;
bool fast = false;
bool lazyGridLoad = false;
bool force = false;

// NetCDF4/HDF5 filter
std::string filterSpec;

int cdoShuffle = 0;
bool cdoCompress = false;
int cdoCompType = CDI_COMPRESS_NONE;
int cdoCompLevel = 0;
bool cdoInteractive = false;
bool cdoVerbose = false;
int CDO_Rusage = 0;
int cdoExitStatus = 0;
bool Timer = false;

bool CheckDatarange = false;

int CDO_flt_digits = 7;   // TODO:rename
int CDO_dbl_digits = 15;  // TODO:rename

bool Use_FFTW = true;
bool VersionInfo = true;
int CMOR_Mode = false;

bool CDO_diagnostic = false;

MemType CDO_Memtype(MemType::Native);
bool CDO_Async_Read = false;

int CDO_Reduce_Dim = false;
int CDO_Append_History = true;
bool CDO_Reset_History = false;
bool PrintFilename = false;
bool CDO_task = false;

unsigned Random_Seed = 1;

int cdoChunkType = CDI_UNDEFID;
int cdoChunkSize = CDI_UNDEFID;
int cdoChunkSizeDimT{ 0 };
int cdoChunkSizeDimZ{ 0 };
int cdoChunkSizeDimY{ 0 };
int cdoChunkSizeDimX{ 0 };
bool cdoOverwriteMode = false;
bool cdoParIO = false;
bool cdoRegulargrid = false;
std::string cdoQueryParameter;
std::vector<std::string> cdoVarnames;

size_t
cdo_num_varnames()
{
  return cdoVarnames.size();
}

bool RemapGenerateWeights{ true };

const char *cdoExpName = nullptr;
}  // namespace Options

namespace Threading
{
int ompNumMaxThreads = 1;
int ompNumUserRequestedThreads = 0;
bool cdoLockIO = false;
}  // namespace Threading

const char *
cdo_comment(void)
{
  return cdo::Version;
}

static bool
filetype_has_szip(int filetype)
{
  return (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2 || filetype == CDI_FILETYPE_NC4
          || filetype != CDI_FILETYPE_NC4C);
}

static bool
filetype_has_zip(int filetype)
{
  return (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NCZARR);
}

void
set_compression(int streamID, int filetype)
{
  if (Options::cdoCompress)
  {
    if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
    {
      Options::cdoCompType = CDI_COMPRESS_SZIP;
      Options::cdoCompLevel = 0;
    }
    else if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NCZARR)
    {
      Options::cdoCompType = CDI_COMPRESS_ZIP;
      Options::cdoCompLevel = 1;
    }
  }

  if (Options::cdoCompType != CDI_COMPRESS_NONE)
  {
    /*      streamDefShuffle(streamID, Options::cdoShuffle);*/

    streamDefCompType(streamID, Options::cdoCompType);
    streamDefCompLevel(streamID, Options::cdoCompLevel);

    if (Options::cdoCompType == CDI_COMPRESS_SZIP && !filetype_has_szip(filetype))
      cdo_warning("SZIP compression not available for non GRIB/NetCDF4 data!");

    if (Options::cdoCompType == CDI_COMPRESS_JPEG && filetype != CDI_FILETYPE_GRB2)
      cdo_warning("JPEG compression not available for non GRIB2 data!");

    if (Options::cdoCompType == CDI_COMPRESS_ZIP && !filetype_has_zip(filetype))
      cdo_warning("Deflate compression not available for non NetCDF4 data!");
  }

  if (Options::filterSpec.size() > 0) { streamDefFilter(streamID, Options::filterSpec.c_str()); }
}

static double pointSearchRadius = 180.0;  // default point search radius in degrees

// set point search radius in degrees
void
cdo_set_search_radius(double searchRadius)
{
  pointSearchRadius = searchRadius;
}

// get point search radius in degrees
double
cdo_get_search_radius(void)
{
  auto searchRadius = pointSearchRadius;
  searchRadius = std::clamp(searchRadius, 0.0, 180.0);
  return searchRadius;
}

void
cdo_print_attributes(std::FILE *fp, int cdiID, int varID, int nblanks)
{
  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int ia = 0; ia < natts; ++ia)
  {
    char attname[CDI_MAX_NAME];
    int atttype, attlen;
    cdiInqAtt(cdiID, varID, ia, attname, &atttype, &attlen);

    if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
        || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
    {
      std::vector<int> attint(attlen);
      cdiInqAttInt(cdiID, varID, attname, attlen, attint.data());
      std::fprintf(fp, "%*s", nblanks, "");
      std::fprintf(fp, "%s = ", attname);
      for (int i = 0; i < attlen; ++i)
      {
        if (i) std::fprintf(fp, ", ");
        std::fprintf(fp, "%d", attint[i]);
      }
      std::fprintf(fp, "\n");
    }
    else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
    {
      char fltstr[128];
      std::vector<double> attflt(attlen);
      cdiInqAttFlt(cdiID, varID, attname, attlen, attflt.data());
      std::fprintf(fp, "%*s", nblanks, "");
      std::fprintf(fp, "%s = ", attname);
      for (int i = 0; i < attlen; ++i)
      {
        if (i) std::fprintf(fp, ", ");
        if (atttype == CDI_DATATYPE_FLT32)
          std::fprintf(fp, "%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
        else
          std::fprintf(fp, "%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
      }
      std::fprintf(fp, "\n");
    }
    else if (atttype == CDI_DATATYPE_TXT)
    {
      std::vector<char> atttxt(attlen + 1);
      cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt.data());
      atttxt[attlen] = 0;
      std::fprintf(fp, "%*s", nblanks, "");
      std::fprintf(fp, "%s = \"%s\"\n", attname, atttxt.data());
    }
  }
}
