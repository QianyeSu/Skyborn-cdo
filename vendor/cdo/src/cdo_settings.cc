/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_FEENABLEEXCEPT
#ifndef __USE_GNU
#define __USE_GNU  // gives us feenableexcept()
#endif
#endif

#if __has_include(<execinfo.h>)
#define HAVE_EXECINFO_H 1
#include <execinfo.h>
#endif

#include <cfenv>
#include <csignal>

#include <cdi.h>

#include "cdo_settings.h"
#include "cdo_features.h"
#include "cdo_default_values.h"
#include "cdo_omp.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "mpim_grid.h"
#include "util_string.h"

void expand_filter_names(std::string &filterSpec);

namespace cdo
{

int netcdf_hdr_pad = 0;

int
evaluate_except_options(std::string const &arg)
{
  int except = -1;
  // clang-format off
  if      (arg == "DIVBYZERO")  except = FE_DIVBYZERO;
  else if (arg == "INEXACT")    except = FE_INEXACT;
  else if (arg == "INVALID")    except = FE_INVALID;
  else if (arg == "OVERFLOW")   except = FE_OVERFLOW;
  else if (arg == "UNDERFLOW")  except = FE_UNDERFLOW;
  else if (arg == "ALL_EXCEPT") except = FE_ALL_EXCEPT;
  // clang-format on
  return except;
}

#ifdef HAVE_FEENABLEEXCEPT
int
set_feenableexcept(int excepts)
{
  return feenableexcept(excepts);
}
#else
int
set_feenableexcept(int excepts)
{
  static fenv_t fenv;
  if (std::fegetenv(&fenv)) return -1;

  (void) excepts;
  int oldExcepts = -1;  // previous masks
#if defined HAVE_FENV_T___CONTROL && defined HAVE_FENV_T___MXCSR
  unsigned newExcepts = ((unsigned) excepts) & FE_ALL_EXCEPT;
  oldExcepts = (int) (fenv.__control & FE_ALL_EXCEPT);

  // unmask
  fenv.__control &= ~newExcepts;
  fenv.__mxcsr &= ~(newExcepts << 7);
#endif

  return (std::fesetenv(&fenv) ? -1 : oldExcepts);
}
#endif

void
set_cdi_options()
{
  if (cdo::dbg())
  {
    fprintf(stderr, "CMOR_Mode           = %d\n", Options::CMOR_Mode);
    fprintf(stderr, "netcdf_hdr_pad      = %d\n", netcdf_hdr_pad);
    fprintf(stderr, "\n");
  }

  if (Threading::cdoLockIO == false) cdiDefGlobal("THREADSAFE", 1);

  if (Options::CMOR_Mode) cdiDefGlobal("CMOR_MODE", Options::CMOR_Mode);  // TODO maybe reposition into effect of "cmor"
  if (Options::CDO_Reduce_Dim) cdiDefGlobal("REDUCE_DIM", Options::CDO_Reduce_Dim);
  if (netcdf_hdr_pad > 0) cdiDefGlobal("NETCDF_HDR_PAD", netcdf_hdr_pad);
}

void
set_external_proj_func(void)
{
#ifdef HAVE_CDI_PROJ_FUNCS
  proj_lonlat_to_lcc_func = proj_lonlat_to_lcc;
  proj_lcc_to_lonlat_func = proj_lcc_to_lonlat;
  proj_lonlat_to_stere_func = proj_lonlat_to_stere;
  proj_stere_to_lonlat_func = proj_stere_to_lonlat;
#endif
}

static void
stackframe()
{
#if defined HAVE_EXECINFO_H && defined HAVE_BACKTRACE
  void *callstack[32];
  auto frames = backtrace(callstack, 32);
  auto messages = backtrace_symbols(callstack, frames);

  fprintf(stderr, "[bt] Execution path:\n");
  if (messages)
  {
    for (int i = 0; i < frames; ++i) fprintf(stderr, "[bt] %s\n", messages[i]);
    std::free(messages);
  }
#endif
}

void
signal_handler(int signo)
{
  if (signo == SIGFPE)
  {
    stackframe();
    cdo_abort("floating-point exception!");
  }
}

void
set_digits(std::string const &arg)
{
  const char *carg = arg.c_str();

  char *ptr1 = 0;
  if (carg != 0 && (int) std::strlen(carg) > 0 && carg[0] != ',') Options::CDO_flt_digits = (int) std::strtol(carg, &ptr1, 10);

  if (Options::CDO_flt_digits < 1 || Options::CDO_flt_digits > 20)
    cdo_abort("Unreasonable value for float significant digits: %d", Options::CDO_flt_digits);

  if (ptr1 && *ptr1 == ',')
  {
    char *ptr2 = 0;
    Options::CDO_dbl_digits = (int) std::strtol(ptr1 + 1, &ptr2, 10);
    if (ptr2 == ptr1 + 1 || Options::CDO_dbl_digits < 1 || Options::CDO_dbl_digits > 20)
      cdo_abort("Unreasonable value for double significant digits: %d", Options::CDO_dbl_digits);
  }
}

void
set_default_filetype(std::string filetypeString)
{
  if (filetypeString.size() > 0)
  {
    std::string numbitsString;
    constexpr char delimiter = '_';

    auto pos = filetypeString.find(delimiter);
    if (pos != std::string::npos)
    {
      numbitsString = filetypeString.substr(pos + 1);
      filetypeString.resize(pos);
    }

    // clang-format off
    if      (filetypeString == "grb2")   CdoDefault::FileType = CDI_FILETYPE_GRB2;
    else if (filetypeString == "grb1")   CdoDefault::FileType = CDI_FILETYPE_GRB;
    else if (filetypeString == "grb")    CdoDefault::FileType = CDI_FILETYPE_GRB;
    else if (filetypeString == "nc2")    CdoDefault::FileType = CDI_FILETYPE_NC2;
    else if (filetypeString == "nc4c")   CdoDefault::FileType = CDI_FILETYPE_NC4C;
    else if (filetypeString == "nc4")    CdoDefault::FileType = CDI_FILETYPE_NC4;
    else if (filetypeString == "nc5")    CdoDefault::FileType = CDI_FILETYPE_NC5;
    else if (filetypeString == "nc1")    CdoDefault::FileType = CDI_FILETYPE_NC;
    else if (filetypeString == "nczarr") CdoDefault::FileType = CDI_FILETYPE_NCZARR;
    else if (filetypeString == "nc")     CdoDefault::FileType = CDI_FILETYPE_NC2;
    else if (filetypeString == "srv")    CdoDefault::FileType = CDI_FILETYPE_SRV;
    else if (filetypeString == "ext")    CdoDefault::FileType = CDI_FILETYPE_EXT;
    else if (filetypeString == "ieg")    CdoDefault::FileType = CDI_FILETYPE_IEG;
    else
    {
      cdo_warning("Unsupported filetype %s!", filetypeString);
      cdo_warning("Available filetypes: grb1/grb2/nc1/nc2/nc4/nc4c/nc5/nczarr/srv/ext/ieg");
      cdo_abort("Unsupported filetype %s!", filetypeString);
    }
    // clang-format on

    if (CdoDefault::FileType != CDI_UNDEFID && numbitsString.size() > 0) { set_default_datatype(numbitsString); }
  }
}

#undef IsBigendian
#define IsBigendian() (u_byteorder.c[sizeof(long) - 1])

void
set_default_datatype(std::string const &datatypeString)
{
  static const union
  {
    unsigned long l;
    unsigned char c[sizeof(long)];
  } u_byteorder = { 1 };
  enum
  {
    D_UINT,
    D_INT,
    D_FLT,
    D_CPX
  };
  int dtype = -1;

  auto datatypestr = datatypeString.c_str();
  auto datatype = std::tolower(*datatypestr);
  // clang-format off
  if      (datatype == 'i') { dtype = D_INT;  datatypestr++; }
  else if (datatype == 'u') { dtype = D_UINT; datatypestr++; }
  else if (datatype == 'f') { dtype = D_FLT;  datatypestr++; }
  else if (datatype == 'c') { dtype = D_CPX;  datatypestr++; }
  else if (datatype == 'p') {                 datatypestr++; }
  // clang-format on

  if (std::isdigit((int) *datatypestr))
  {
    auto nbits = atoi(datatypestr);
    datatypestr++;
    if (nbits >= 10) datatypestr++;

    if (dtype == -1)
    {
      if (nbits > 0 && nbits < 32)
        CdoDefault::DataType = nbits;
      else if (nbits == 32)
        CdoDefault::DataType = (CdoDefault::FileType == CDI_FILETYPE_GRB) ? CDI_DATATYPE_PACK32 : CDI_DATATYPE_FLT32;
      else if (nbits == 64)
        CdoDefault::DataType = CDI_DATATYPE_FLT64;
      else
      {
        cdo_warning("Unsupported number of bits %d!", nbits);
        cdo_warning("Use I8/I16/I32/F32/F64 for nc1/nc2/nc4/nc4c/nc5/nczarr; U8/U16/U32 for nc4/nc4c/nc5/nczarr; F32/F64 for "
                    "grb2/srv/ext/ieg; P1 - P24 for grb1/grb2.");
        cdo_abort("Unsupported number of bits!");
      }
    }
    else
    {
      // clang-format off
      if (dtype == D_INT)
      {
        if      (nbits ==  8) CdoDefault::DataType = CDI_DATATYPE_INT8;
        else if (nbits == 16) CdoDefault::DataType = CDI_DATATYPE_INT16;
        else if (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_INT32;
        else cdo_abort("Unsupported number of bits = %d for datatype INT!", nbits);
      }
      else if (dtype == D_UINT)
      {
        if      (nbits ==  8) CdoDefault::DataType = CDI_DATATYPE_UINT8;
        else if (nbits == 16) CdoDefault::DataType = CDI_DATATYPE_UINT16;
        else if (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_UINT32;
        else cdo_abort("Unsupported number of bits = %d for datatype UINT!", nbits);
      }
      else if (dtype == D_FLT)
      {
        if      (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_FLT32;
        else if (nbits == 64) CdoDefault::DataType = CDI_DATATYPE_FLT64;
        else cdo_abort("Unsupported number of bits = %d for datatype FLT!", nbits);
      }
      else if (dtype == D_CPX)
      {
        if      (nbits == 32) CdoDefault::DataType = CDI_DATATYPE_CPX32;
        else if (nbits == 64) CdoDefault::DataType = CDI_DATATYPE_CPX64;
        else cdo_abort("Unsupported number of bits = %d for datatype CPX!", nbits);
      }
      // clang-format on
    }
  }

  if (*datatypestr != 0)
  {
    if (*datatypestr == 'l' || *datatypestr == 'L')
    {
      if (IsBigendian()) CdoDefault::Byteorder = CDI_LITTLEENDIAN;
    }
    else if (*datatypestr == 'b' || *datatypestr == 'B')
    {
      if (!IsBigendian()) CdoDefault::Byteorder = CDI_BIGENDIAN;
    }
    else { cdo_abort("Unsupported character in number of bytes: >%s< !", datatypestr); }

    datatypestr++;
    if (*datatypestr != 0) { cdo_abort("Unsupported character in number of bytes: >%s< !", datatypestr); }
  }

  if (CdoDefault::DataType == -1) { cdo_abort("Number of bits undefined!"); }
}

void
set_filterspec(std::string const &arg)
{
  if (arg.size() > 0)
  {
    if (Options::filterSpec.size() > 0) { cdo_abort("Filter specs already defined! Only one filter specs is allowed."); }
    Options::filterSpec = string_to_lower(arg);
    expand_filter_names(Options::filterSpec);
  }
  else { cdo_abort("Filter spec missing!"); }
}

void
set_compression_type(std::string const &arg)
{
  size_t len = arg.size();

  if (arg == "szip")
  {
    Options::cdoCompType = CDI_COMPRESS_SZIP;
    Options::cdoCompLevel = 0;
  }
  else if (arg == "aec" || arg == "ccsds")
  {
    Options::cdoCompType = CDI_COMPRESS_AEC;
    Options::cdoCompLevel = 0;
  }
  else if (arg == "jpeg")
  {
    Options::cdoCompType = CDI_COMPRESS_JPEG;
    Options::cdoCompLevel = 0;
  }
  else if (arg.starts_with("zip"))
  {
    Options::cdoCompType = CDI_COMPRESS_ZIP;
    Options::cdoCompLevel = (len == 5 && arg[3] == '_' && std::isdigit(arg[4])) ? std::atoi(&arg.c_str()[4]) : 1;
  }
  else if (arg.starts_with("zstd"))
  {
    int filterIdZstd = 32015;
    int zstdLevel = (len >= 6 && len <= 7 && arg[4] == '_' && std::isdigit(arg[5])) ? std::atoi(&arg.c_str()[5]) : 1;
    if (Options::filterSpec.size() > 0) { cdo_abort("Filter specs already defined! Only one filter specs is allowed."); }
    Options::filterSpec = std::to_string(filterIdZstd) + "," + std::to_string(zstdLevel);
  }
  else { cdo_abort("Compression type '%s' unsupported!", arg); }
}

void
set_chunktype(std::string const &arg)
{
  // clang-format off
  if      ("auto"  == arg) Options::cdoChunkType = CDI_CHUNK_AUTO;
  else if ("grid"  == arg) Options::cdoChunkType = CDI_CHUNK_GRID;
  else if ("lines" == arg) Options::cdoChunkType = CDI_CHUNK_LINES;
  else cdo_abort("Chunk type '%s' unsupported!", arg);
  // clang-format on
}

void
evaluate_color_options(std::string const &arg)
{
  // clang-format off
  if      ("all"  == arg) mpmo_color_set(All);
  else if ("auto" == arg) mpmo_color_set(Auto);
  else if ("no"   == arg) mpmo_color_set(No);
  else cdo_abort("Color option <%s> unknown. Known options: auto, all, no", Yellow(arg));
  // clang-format on
}

void
setup_openMP(int numThreads)
{
#ifdef _OPENMP
  if (numThreads <= 0) numThreads = 1;
  omp_set_num_threads(numThreads);

  Threading::ompNumMaxThreads = omp_get_max_threads();
  if (omp_get_max_threads() > omp_get_num_procs())
    fprintf(stderr, "Warning: Number of OMP threads=%d is greater than number of Cores=%d!\n", omp_get_max_threads(),
            omp_get_num_procs());

  if (Threading::ompNumMaxThreads < numThreads)
    fprintf(stderr, "Warning: omp_get_max_threads() returns %d!\n", Threading::ompNumMaxThreads);

  if (cdo::dbg()) cdo::features::print_openmp_info();

  if (Options::cdoVerbose)
  {
    fprintf(stderr, " OpenMP:  num_procs=%d  max_threads=%d", omp_get_num_procs(), omp_get_max_threads());
#ifdef HAVE_OPENMP4
#ifndef __ICC
    fprintf(stderr, "  num_devices=%d", omp_get_num_devices());
#endif
#endif
    fprintf(stderr, "\n");
  }
#else
  if (numThreads > 1) fprintf(stderr, "Warning: Option -P failed, OpenMP support not compiled in!\n");
#endif
}

};  // namespace cdo
