/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>  // omp_get_num_procs()
#endif

#if __has_include(<sys/resource.h>)
#define HAVE_SYS_RESOURCE_H 1
#include <sys/resource.h>  // getrlimit
#endif

#include <cdi.h>

#include <map>

#ifdef HAVE_LIBNETCDF
#include <netcdf_meta.h>
#endif

#ifdef HAVE_HDF5_H
#include <hdf5.h>
#endif

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/xmlversion.h>
#endif

#ifdef HAVE_CURL_CURL_H
#include <curl/curl.h>
#endif

#ifdef HAVE_PROJ_H
#include <proj.h>
#endif

#include <cstdio>
#include <cstring>
#include <iostream>
#include <inttypes.h>

#include "cdo_features.h"
#include "cdo_default_values.h"
#include "cdo_rlimit.h"
#include "cdo_omp.h"
#include "cdo_options.h"
#include "mpmo_color.h"
#include "util_string.h"
#include "cpp_lib.h"
#include "lib/yac/src/yac_version.h"

#include <thread>  // std::thread::hardware_concurrency()

extern "C" size_t getMemorySize(void);

namespace cdo
{
namespace features
{

#ifdef HAVE_LIBHDF5
constexpr bool have_hdf5 = true;
#else
constexpr bool have_hdf5 = false;
#endif

#ifdef HAVE_LIBCMOR
constexpr bool have_cmor = true;
#else
constexpr bool have_cmor = false;
#endif

#ifdef HAVE_LIBMAGICS
constexpr bool have_magics = true;
#else
constexpr bool have_magics = false;
#endif

#ifdef _OPENMP
constexpr bool have_openmp = true;
#else
constexpr bool have_openmp = false;
#endif

#ifdef HAVE_LIBPROJ
constexpr bool have_proj = true;
#else
constexpr bool have_proj = false;
#endif

#ifdef HAVE_LIBPTHREAD
constexpr bool have_threads = true;
#else
constexpr bool have_threads = false;
#endif

#if __has_include(<wordexp.h>)
constexpr bool have_wordexp = true;
#else
constexpr bool have_wordexp = false;
#endif

#ifdef HIRLAM_EXTENSIONS
constexpr bool has_hirlam_extensions = true;
#else
constexpr bool has_hirlam_extensions = false;
#endif

static std::map<std::string, std::pair<std::string, bool>> configMap
    = { { "has-srv", { "SERVICE", cdiHaveFiletype(CDI_FILETYPE_SRV) } },
        { "has-ext", { "EXTRA", cdiHaveFiletype(CDI_FILETYPE_EXT) } },
        { "has-ieg", { "IEG", cdiHaveFiletype(CDI_FILETYPE_IEG) } },
        { "has-grb", { "GRIB 1", cdiHaveFiletype(CDI_FILETYPE_GRB) } },
        { "has-grb1", { "GRIB 1", cdiHaveFiletype(CDI_FILETYPE_GRB) } },
        { "has-grb2", { "GRIB 2", cdiHaveFiletype(CDI_FILETYPE_GRB2) } },
        { "has-nc", { "NetCDF", cdiHaveFiletype(CDI_FILETYPE_NC) } },
        { "has-nc2", { "NetCDF 2", cdiHaveFiletype(CDI_FILETYPE_NC2) } },
        { "has-nc4", { "NetCDF 4", cdiHaveFiletype(CDI_FILETYPE_NC4) } },
        { "has-nc4c", { "NetCDF 4 classic", cdiHaveFiletype(CDI_FILETYPE_NC4C) } },
        { "has-nc5", { "NetCDF 5", cdiHaveFiletype(CDI_FILETYPE_NC5) } },
        { "has-nczarr", { "NetCDF 4 zarr", cdiHaveFiletype(CDI_FILETYPE_NCZARR) } },
        { "has-hdf5", { "HDF5", have_hdf5 } },
        { "has-ncfilter", { "NetCDF 4 filter", cdiGetConfig(CDI_NC_HAS_FILTER) } },
        { "has-ncs3", { "NetCDF 4 s3", cdiGetConfig(CDI_NC_HAS_S3) } },
        { "has-ncdap", { "NetCDF OpenDap", cdiGetConfig(CDI_NC_HAS_DAP) } },
        { "has-cgribex", { "CGRIBEX", cdiGetConfig(CDI_HAS_CGRIBEX) } },
        { "has-cmor", { "CMOR", have_cmor } },
        { "has-magics", { "MAGICS", have_magics } },
        { "has-openmp", { "OPENMP", have_openmp } },
        { "has-proj", { "PROJ", have_proj } },
        { "has-threads", { "PTHREADS", have_threads } },
        { "has-wordexp", { "WORDEXP", have_wordexp } },
        { "has-hirlam_extensions", { "HIRLAM_EXTENSIONS", has_hirlam_extensions } } };

void
print_features()
{
  constexpr size_t gigaByte = 1024 * 1024 * 1024;
  auto fp = stdout;
  std::fprintf(fp, "Features: ");
  size_t rssCur = (size_t) cdo::get_rss_cur() / gigaByte;
  size_t memorySize = getMemorySize() / gigaByte;
  if (rssCur > 0 && rssCur < memorySize) std::fprintf(fp, "%zu/", rssCur);
  if (memorySize > 0) std::fprintf(fp, "%zuGB ", memorySize);
  auto concurrentThreads = std::thread::hardware_concurrency();
#ifdef _OPENMP
  unsigned numProcs = omp_get_num_procs();
  if (numProcs < concurrentThreads) std::fprintf(fp, "%u/", numProcs);
#endif
  std::fprintf(fp, "%uthreads", concurrentThreads);
#ifdef HAVE_OPENMP4
  auto numDevices = omp_get_num_devices();
  if (numDevices > 0) std::fprintf(fp, " %ddevices", numDevices);
#endif
  std::fprintf(fp, " c++%d", (int) ((__cplusplus - 200000) / 100));
#ifdef __FAST_MATH__
  std::fprintf(fp, "/fastmath");
#endif
#ifdef _OPENMP
  std::fprintf(fp, " OpenMP");
#if defined(HAVE_OPENMP52)
  std::fprintf(fp, "52");
#elif defined(HAVE_OPENMP51)
  std::fprintf(fp, "51");
#elif defined(HAVE_OPENMP5)
  std::fprintf(fp, "5");
#elif defined(HAVE_OPENMP45)
  std::fprintf(fp, "45");
#elif defined(HAVE_OPENMP4)
  std::fprintf(fp, "4");
#elif defined(HAVE_OPENMP3)
  std::fprintf(fp, "3");
#endif
#endif
#ifdef HAVE_CF_INTERFACE
  std::fprintf(fp, " Fortran");
#endif
#ifdef HAVE_LIBPTHREAD
  std::fprintf(fp, " pthreads");
#endif
#ifdef HAVE_LIBHDF5
  std::fprintf(fp, " HDF5");
#endif
#ifdef HAVE_LIBNETCDF
  std::fprintf(fp, " NC4");
#endif
  if (cdiGetConfig(CDI_NC_HAS_S3)) { std::fprintf(fp, "/S3"); }
  if (cdiGetConfig(CDI_NC_HAS_HDF5)) { std::fprintf(fp, "/HDF5"); }
#ifdef HAVE_NC4HDF5_THREADSAFE
  std::fprintf(fp, "/threadsafe");
#endif
  if (cdiGetConfig(CDI_NC_HAS_DAP)) { std::fprintf(fp, " dap"); }
#ifdef HAVE_LIBSZ
  std::fprintf(fp, " sz");
#endif
/*
#ifdef HAVE_LIBZ
fprintf(fp, " z");
#endif
*/
#ifdef HAVE_LIBUDUNITS2
  std::fprintf(fp, " udunits2");
#endif
#ifdef HAVE_LIBPROJ
  std::fprintf(fp, " proj");
#endif
#ifdef HAVE_LIBXML2
  std::fprintf(fp, " xml2");
#endif
#ifdef HAVE_LIBMAGICS
  std::fprintf(fp, " magics");
#endif
#ifdef HAVE_LIBDRMAA
  std::fprintf(fp, " drmaa");
#endif
#ifdef HAVE_LIBCURL
  std::fprintf(fp, " curl");
#endif
#ifdef HAVE_LIBFFTW3
  std::fprintf(fp, " fftw3");
#endif
#ifdef HAVE_LIBCMOR
  std::fprintf(fp, " cmor");
#endif
#ifdef HIRLAM_EXTENSIONS
  std::fprintf(fp, " hirlam_extensions");
#endif
#if defined(__AVX2__)
  std::fprintf(fp, " avx2");
#elif defined(__AVX__)
  std::fprintf(fp, " avx");
#elif defined(__SSE4_2__)
  std::fprintf(fp, " sse4_2");
#elif defined(__SSE4_1__)
  std::fprintf(fp, " sse4_1");
#elif defined(__SSE3__)
  std::fprintf(fp, " sse3");
#elif defined(__SSE2__)
  std::fprintf(fp, " sse2");
#endif
  std::fprintf(fp, "\n");
}

void
activate_hdf5_diag()
{
  // HDF5-DIAG messages seems to be activate by a hdf5 function call
#ifdef HAVE_LIBHDF5
#ifdef H5_VERS_MAJOR
  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
#endif
#endif
}

void
print_libraries()
{
  auto fp = stdout;
  std::fprintf(fp, "Libraries:");
  std::fprintf(fp, " yac/%s", YAC_VERSION);
#ifdef HAVE_LIBNETCDF
  std::fprintf(fp, " NetCDF");
#ifdef NC_VERSION
  std::fprintf(fp, "/%s", NC_VERSION);
#endif
#endif
#ifdef HAVE_LIBHDF5
  std::fprintf(fp, " HDF5");
#ifdef H5_VERS_MAJOR
  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  std::fprintf(fp, "/%u.%u.%u", h5l_majnum, h5l_minnum, h5l_relnum);

  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR, h5h_relnum = H5_VERS_RELEASE;
  if ((h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum) || (h5h_relnum != h5l_relnum))
    std::fprintf(fp, "(h%u.%u.%u)", h5h_majnum, h5h_minnum, h5h_relnum);
#endif
#endif
/*
#ifdef HAVE_LIBZ
{
  std::fprintf(fp, " zlib/%s", zlibVersion());
#ifdef ZLIB_VERSION
  if (std::strcmp(ZLIB_VERSION, zlibVersion()) != 0)
    std::fprintf(fp, "(h%s)", ZLIB_VERSION);
#else
  std::fprintf(fp, "(header not found)");
#endif
}
#endif
*/
#ifdef HAVE_LIBPROJ
  std::fprintf(fp, " proj");
#ifdef PROJ_VERSION_MAJOR
  std::fprintf(fp, "/%u.%u.%u", PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBCMOR
  std::fprintf(fp, " cmor");
#ifdef CMOR_VERSION_MAJOR
  std::fprintf(fp, "/%u.%u.%u", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBXML2
  std::fprintf(fp, " xml2");
#ifdef LIBXML_DOTTED_VERSION
  std::fprintf(fp, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#ifdef HAVE_LIBCURL
  {
    auto version_data = curl_version_info(CURLVERSION_NOW);
    std::fprintf(fp, " curl/%s", version_data->version);
#ifdef LIBCURL_VERSION
    if (std::strcmp(LIBCURL_VERSION, version_data->version) != 0) std::fprintf(fp, "(h%s)", LIBCURL_VERSION);
#else
    std::fprintf(fp, "(header not found)");
#endif
  }
#endif

#ifdef HAVE_LIBMAGICS
  {
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#ifdef HAVE_SYS_TYPES_H
#undef HAVE_SYS_TYPES_H
#endif
#include <magics_config.h>
#ifdef MAGICS_VERSION
    std::fprintf(fp, " magics/%s", MAGICS_VERSION);
#endif
  }
#endif

  std::fprintf(fp, "\n");
}

void
print_argument_options()
{
  std::fprintf(stdout, "Available config option:\n");
  std::fprintf(stdout, "\n");
  for (auto const &entry : configMap)
    std::fprintf(stdout, "  %-12s  whether %s is enabled\n", entry.first.c_str(), entry.second.first.c_str());
}

int
print_config(std::string const &option)
{
  int status = EXIT_SUCCESS;

  if ("all-json" == option || "all" == option)
  {
    std::cout << "{\n";
    int i = 0;
    for (auto const &entry : configMap)
    {
      if (i++) std::fprintf(stdout, ",\n");
      std::cout << "\"" << entry.first << "\":\"" << (entry.second.second ? "yes" : "no") << "\"";
    }
    std::cout << "\n}\n";
  }
  else
  {
    auto foundOption = false;
    for (auto const &entry : configMap)
    {
      if (entry.first == option)
      {
        foundOption = true;
        std::cout << (entry.second.second ? "yes" : "no") << "\n";
      }
    }

    if (!foundOption)
    {
      std::fprintf(stdout, "unknown config option: %s\n", option.c_str());
      std::fprintf(stdout, "\n");
      print_argument_options();

      status = EXIT_FAILURE;
    }
  }

  return status;
}

static auto
alignof_address(void *ptr) -> int
{
  auto n = reinterpret_cast<int64_t>(ptr);
  return (int) (n & (-n));
}

static auto
alignof_malloc_data(std::vector<int> const &tsize) -> int
{
  int align = (1 << 30);
  auto n = tsize.size();

  std::vector<double *> ptr(n);

  for (size_t i = 0; i < n; ++i)
  {
    ptr[i] = (double *) std::malloc(tsize[i] * sizeof(double));
    align = std::min(align, alignof_address(ptr[i]));
  }
  for (auto &p : ptr) std::free(p);

  return align;
}

static auto
alignof_vector_data(std::vector<int> const &tsize) -> int
{
  int align = 1 << 30;
  auto n = tsize.size();

  std::vector<std::vector<double>> ptr(n);

  for (size_t i = 0; i < n; ++i)
  {
    ptr[i].resize(tsize[i]);
    align = std::min(align, alignof_address(ptr[i].data()));
  }

  return align;
}

void
print_system_info()
{
  std::fprintf(stderr, "\n");
  std::fprintf(stderr, "CDO_Color             = %d\n", mpmo_get_color_mode());
  std::fprintf(stderr, "Options::CDO_Reset_History = %d\n", Options::CDO_Reset_History);
  std::fprintf(stderr, "CDO_FileSuffix        = %s\n", cdo::FileSuffix.c_str());
  std::fprintf(stderr, "CdoDefault::FileType  = %d\n", CdoDefault::FileType);
  std::fprintf(stderr, "CdoDefault::DataType  = %d\n", CdoDefault::DataType);
  std::fprintf(stderr, "CdoDefault::Byteorder = %d\n", CdoDefault::Byteorder);
  std::fprintf(stderr, "CdoDefault::TableID   = %d\n", CdoDefault::TableID);
  std::fprintf(stderr, "\n");

  const char *envstr;
  envstr = getenv("HOSTTYPE");
  if (envstr) std::fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if (envstr) std::fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if (envstr) std::fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if (envstr) std::fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  std::fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  std::fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  std::fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  std::fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  std::fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  std::fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  std::fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  std::fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  std::fprintf(stderr, "Predefined: __SSE2__\n");
#endif
  std::fprintf(stderr, "\n");

  std::fprintf(stderr, "sizeof(size_t)      = %zu\n", sizeof(size_t));
  {
    constexpr size_t megaByte = 1024 * 1024;
    std::vector<int> numElements = { 1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025, 1 * megaByte };
    std::fprintf(stderr, "alignof malloc data = %d\n", alignof_malloc_data(numElements));
    std::fprintf(stderr, "alignof malloc big  = %d\n", alignof_malloc_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
    std::fprintf(stderr, "alignof vector data = %d\n", alignof_vector_data(numElements));
    std::fprintf(stderr, "alignof vector big  = %d\n", alignof_vector_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
  }
  std::fprintf(stderr, "\n");

#ifdef HAVE_MMAP
  std::fprintf(stderr, "HAVE_MMAP\n");
#endif
  std::fprintf(stderr, "\n");

#ifdef _OPENACC
  std::fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
  // OPENMP3:   201107
  // OPENMP4:   201307 gcc 4.9
  // OPENMP45:  201511
#ifdef _OPENMP
  std::fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
  std::fprintf(stderr, "__cplusplus         = %ld\n", (long) __cplusplus);
#ifdef __GNUC__
  std::fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#ifdef __GNUC_MINOR__
  std::fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#ifdef __ICC
  std::fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#ifdef __STDC__
  std::fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#ifdef __STD_VERSION__
  std::fprintf(stderr, "STD VERSION         = %ld\n", (long) __STD_VERSION__);
#endif
#ifdef __STDC_VERSION__
  std::fprintf(stderr, "STDC VERSION        = %ld\n", (long) __STDC_VERSION__);
#endif
#ifdef __STD_HOSTED__
  std::fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#ifdef FLT_EVAL_METHOD
  std::fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#ifdef FP_FAST_FMA
  std::fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#ifdef __FAST_MATH__
  std::fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  std::fprintf(stderr, "\n");

#ifdef _SC_VERSION
  std::fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#ifdef _SC_ARG_MAX
  std::fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#ifdef _SC_CHILD_MAX
  std::fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#ifdef _SC_STREAM_MAX
  std::fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#ifdef _SC_OPEN_MAX
  std::fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#ifdef _SC_PAGESIZE
  std::fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  std::fprintf(stderr, "\n");

  cdo::print_rlimits();

  std::fprintf(stderr, "\n");
}

static void
print_filetypes(std::FILE *fp)
{
  // clang-format off
  const std::vector<std::pair<int, std::string>> fileTypes = {
      { CDI_FILETYPE_SRV,    "srv" },
      { CDI_FILETYPE_EXT,    "ext" },
      { CDI_FILETYPE_IEG,    "ieg" },
      { CDI_FILETYPE_GRB,    "grb1" },
      { CDI_FILETYPE_GRB2,   "grb2" },
      { CDI_FILETYPE_NC,     "nc1" },
      { CDI_FILETYPE_NC2,    "nc2" },
      { CDI_FILETYPE_NC4,    "nc4" },
      { CDI_FILETYPE_NC4C,   "nc4c" },
      { CDI_FILETYPE_NC5,    "nc5" },
      { CDI_FILETYPE_NCZARR, "nczarr" }
  };
  // clang-format on

  std::fprintf(fp, "CDI file types: ");
  set_text_color(fp, BRIGHT, GREEN);
  for (auto const &[type, name] : fileTypes)
    if (cdiHaveFiletype(type)) std::fprintf(fp, "%s ", name.c_str());
  reset_text_color(fp);
  std::fprintf(fp, "\n");
}

void
version()
{
  auto fp = stdout;
  std::fprintf(fp, "%s\n", cdo::Version);
#ifdef SYSTEM_TYPE
  std::fprintf(fp, "System: %s\n", SYSTEM_TYPE);
#endif
#ifdef CXX_COMPILER
  std::fprintf(fp, "CXX Compiler: %s\n", CXX_COMPILER);
#ifdef CXX_VERSION
  std::fprintf(fp, "CXX version : %s\n", CXX_VERSION);
#endif
  std::fprintf(fp, "CXX library :");
#ifdef HAVE_LIB_RANGES_ZIP
  std::fprintf(fp, " ranges_zip");
#endif
#ifdef HAVE_LIB_MDSPAN
  std::fprintf(fp, " mdspan");
#endif
  std::fprintf(fp, "\n");
#endif
#ifdef C_COMPILER
  std::fprintf(fp, "C Compiler: %s\n", C_COMPILER);
#ifdef C_VERSION
  std::fprintf(fp, "C version : %s\n", C_VERSION);
#endif
#endif
#ifdef F77_COMPILER
  std::fprintf(fp, "F77 Compiler: %s\n", F77_COMPILER);
#ifdef F77_VERSION
  std::fprintf(fp, "F77 version : %s\n", F77_VERSION);
#endif
#endif

  cdo::features::print_features();
  cdo::features::print_libraries();

#ifdef CDI_SIZE_TYPE
#define CDO_STRINGIFY(x) #x
#define CDO_TOSTRING(x) CDO_STRINGIFY(x)
  std::fprintf(fp, "CDI data types: SizeType=%s\n", CDO_TOSTRING(CDI_SIZE_TYPE));
#endif

  print_filetypes(fp);

  cdiPrintVersion();
  std::fprintf(fp, "\n");
}

/*
#if defined(__APPLE__) && defined(__MACH__)
#include <libproc.h>
#endif
*/
void
print_rusage()
{
#if defined HAVE_SYS_RESOURCE_H && defined RUSAGE_SELF
  struct rusage ru;
  auto status = getrusage(RUSAGE_SELF, &ru);
  if (status == 0)
  {
    double ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
    double st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

    std::fprintf(stderr, "  User time:     %.3f seconds\n", ut);
    std::fprintf(stderr, "  System time:   %.3f seconds\n", st);
    std::fprintf(stderr, "  Total time:    %.3f seconds\n", ut + st);
#if defined(__APPLE__) && defined(__MACH__)
    std::fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024.0 * 1024.0));
#else
    std::fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024.0));
#endif
    std::fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
    std::fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
    std::fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
    std::fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
    std::fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
    /*
#if defined(__APPLE__) && defined(__MACH__)
#ifdef RUSAGE_INFO_CURRENT
    pid_t pid = fork();
    rusage_info_current ruinfo;
    int rusage_ret = proc_pid_rusage(pid, RUSAGE_INFO_CURRENT, (void **)&ruinfo);
    if (rusage_ret >= 0 && ruinfo.ri_lifetime_max_phys_footprint > 0)
      {
        std::fprintf(stderr, "  Peak memory:   %.2f MBytes\n", ruinfo.ri_lifetime_max_phys_footprint / 1.0);
        std::fprintf(stderr, "  Peak memory:   %.2f MBytes\n", ruinfo.ri_lifetime_max_phys_footprint / (1024.0 * 1024.0));
      }
#endif
#endif
    */
  }
#endif
}

#ifdef _OPENMP
void
print_openmp_info()
{
  std::fprintf(stderr, "OMP num procs       = %d\n", omp_get_num_procs());
  std::fprintf(stderr, "OMP max threads     = %d\n", omp_get_max_threads());
  std::fprintf(stderr, "OMP num threads     = %d\n", omp_get_num_threads());
#ifndef HAVE_OPENMP3
  std::fprintf(stderr, "OMP thread limit    = %d\n", omp_get_thread_limit());
  omp_sched_t kind;
  int modifer;
  omp_get_schedule(&kind, &modifer);
  std::fprintf(stderr, "OMP schedule        = %d (1:static; 2:dynamic; 3:guided; 4:auto)\n", (int) kind);
#endif
#ifdef HAVE_OPENMP4
  std::fprintf(stderr, "OMP proc bind       = %d (0:false; 1:true; 2:master; 3:close; 4:spread)\n", (int) omp_get_proc_bind());
#ifndef __ICC
  std::fprintf(stderr, "OMP num devices     = %d\n", omp_get_num_devices());
#endif
#endif
}
#endif

};  // namespace features
};  // namespace cdo
