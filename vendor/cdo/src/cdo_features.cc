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

#include "cdf_config.h"

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
        { "has-ncfilter", { "NetCDF 4 filter", cdi_has_ncfilter() } },
        { "has-ncdap", { "NetCDF OpenDap", cdi_has_ncdap() } },
        { "has-cgribex", { "CGRIBEX", cdi_has_cgribex() } },
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
  fprintf(fp, "Features: ");
  size_t rssCur = (size_t) cdo::get_rss_cur() / gigaByte;
  size_t memorySize = getMemorySize() / gigaByte;
  if (rssCur > 0 && rssCur < memorySize) fprintf(fp, "%zu/", rssCur);
  if (memorySize > 0) fprintf(fp, "%zuGB ", memorySize);
  auto concurrentThreads = std::thread::hardware_concurrency();
#ifdef _OPENMP
  unsigned numProcs = omp_get_num_procs();
  if (numProcs < concurrentThreads) fprintf(fp, "%u/", numProcs);
#endif
  fprintf(fp, "%uthreads", concurrentThreads);
#ifdef HAVE_OPENMP4
  auto numDevices = omp_get_num_devices();
  if (numDevices > 0) fprintf(fp, " %ddevices", numDevices);
#endif
  fprintf(fp, " c++%d", (int) ((__cplusplus - 200000) / 100));
#ifdef __FAST_MATH__
  fprintf(fp, "/fastmath");
#endif
#ifdef _OPENMP
  fprintf(fp, " OpenMP");
#if defined(HAVE_OPENMP52)
  fprintf(fp, "52");
#elif defined(HAVE_OPENMP51)
  fprintf(fp, "51");
#elif defined(HAVE_OPENMP5)
  fprintf(fp, "5");
#elif defined(HAVE_OPENMP45)
  fprintf(fp, "45");
#elif defined(HAVE_OPENMP4)
  fprintf(fp, "4");
#elif defined(HAVE_OPENMP3)
  fprintf(fp, "3");
#endif
#endif
#ifdef HAVE_CF_INTERFACE
  fprintf(fp, " Fortran");
#endif
#ifdef HAVE_LIBPTHREAD
  fprintf(fp, " pthreads");
#endif
#ifdef HAVE_LIBHDF5
  fprintf(fp, " HDF5");
#endif
#ifdef HAVE_LIBNETCDF
  fprintf(fp, " NC4");
#endif
#ifdef HAVE_NC4S3
  fprintf(fp, "/S3");
#endif
#ifdef HAVE_NC4HDF5
  fprintf(fp, "/HDF5");
#endif
#ifdef HAVE_NC4HDF5_THREADSAFE
  fprintf(fp, "/threadsafe");
#endif
  if (cdi_has_ncdap()) fprintf(fp, " dap");
#ifdef HAVE_LIBSZ
  fprintf(fp, " sz");
#endif
/*
#ifdef HAVE_LIBZ
fprintf(fp, " z");
#endif
*/
#ifdef HAVE_LIBUDUNITS2
  fprintf(fp, " udunits2");
#endif
#ifdef HAVE_LIBPROJ
  fprintf(fp, " proj");
#endif
#ifdef HAVE_LIBXML2
  fprintf(fp, " xml2");
#endif
#ifdef HAVE_LIBMAGICS
  fprintf(fp, " magics");
#endif
#ifdef HAVE_LIBDRMAA
  fprintf(fp, " drmaa");
#endif
#ifdef HAVE_LIBCURL
  fprintf(fp, " curl");
#endif
#ifdef HAVE_LIBFFTW3
  fprintf(fp, " fftw3");
#endif
#ifdef HAVE_LIBCMOR
  fprintf(fp, " cmor");
#endif
#ifdef HIRLAM_EXTENSIONS
  fprintf(fp, " hirlam_extensions");
#endif
#if defined(__AVX2__)
  fprintf(fp, " avx2");
#elif defined(__AVX__)
  fprintf(fp, " avx");
#elif defined(__SSE4_2__)
  fprintf(fp, " sse4_2");
#elif defined(__SSE4_1__)
  fprintf(fp, " sse4_1");
#elif defined(__SSE3__)
  fprintf(fp, " sse3");
#elif defined(__SSE2__)
  fprintf(fp, " sse2");
#endif
  fprintf(fp, "\n");
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
  fprintf(fp, "Libraries:");
  fprintf(fp, " yac/%s", YAC_VERSION);
#ifdef HAVE_LIBNETCDF
  fprintf(fp, " NetCDF");
#ifdef NC_VERSION
  fprintf(fp, "/%s", NC_VERSION);
#endif
#endif
#ifdef HAVE_LIBHDF5
  fprintf(fp, " HDF5");
#ifdef H5_VERS_MAJOR
  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  fprintf(fp, "/%u.%u.%u", h5l_majnum, h5l_minnum, h5l_relnum);

  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR, h5h_relnum = H5_VERS_RELEASE;
  if ((h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum) || (h5h_relnum != h5l_relnum))
    fprintf(fp, "(h%u.%u.%u)", h5h_majnum, h5h_minnum, h5h_relnum);
#endif
#endif
/*
#ifdef HAVE_LIBZ
{
  fprintf(fp, " zlib/%s", zlibVersion());
#ifdef ZLIB_VERSION
  if (std::strcmp(ZLIB_VERSION, zlibVersion()) != 0)
    fprintf(fp, "(h%s)", ZLIB_VERSION);
#else
  fprintf(fp, "(header not found)");
#endif
}
#endif
*/
#ifdef HAVE_LIBPROJ
  fprintf(fp, " proj");
#ifdef PROJ_VERSION_MAJOR
  fprintf(fp, "/%u.%u.%u", PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBCMOR
  fprintf(fp, " cmor");
#ifdef CMOR_VERSION_MAJOR
  fprintf(fp, "/%u.%u.%u", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBXML2
  fprintf(fp, " xml2");
#ifdef LIBXML_DOTTED_VERSION
  fprintf(fp, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#ifdef HAVE_LIBCURL
  {
    auto version_data = curl_version_info(CURLVERSION_NOW);
    fprintf(fp, " curl/%s", version_data->version);
#ifdef LIBCURL_VERSION
    if (std::strcmp(LIBCURL_VERSION, version_data->version) != 0) fprintf(fp, "(h%s)", LIBCURL_VERSION);
#else
    fprintf(fp, "(header not found)");
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
    fprintf(fp, " magics/%s", MAGICS_VERSION);
#endif
  }
#endif

  fprintf(fp, "\n");
}

void
print_argument_options()
{
  fprintf(stdout, "Available config option:\n");
  fprintf(stdout, "\n");
  for (auto const &entry : configMap)
    fprintf(stdout, "  %-12s  whether %s is enabled\n", entry.first.c_str(), entry.second.first.c_str());
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
      if (i++) fprintf(stdout, ",\n");
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
      fprintf(stdout, "unknown config option: %s\n", option.c_str());
      fprintf(stdout, "\n");
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
  fprintf(stderr, "\n");
  fprintf(stderr, "CDO_Color             = %d\n", mpmo_get_color_mode());
  fprintf(stderr, "Options::CDO_Reset_History = %d\n", Options::CDO_Reset_History);
  fprintf(stderr, "CDO_FileSuffix        = %s\n", cdo::FileSuffix.c_str());
  fprintf(stderr, "CdoDefault::FileType  = %d\n", CdoDefault::FileType);
  fprintf(stderr, "CdoDefault::DataType  = %d\n", CdoDefault::DataType);
  fprintf(stderr, "CdoDefault::Byteorder = %d\n", CdoDefault::Byteorder);
  fprintf(stderr, "CdoDefault::TableID   = %d\n", CdoDefault::TableID);
  fprintf(stderr, "\n");

  const char *envstr;
  envstr = getenv("HOSTTYPE");
  if (envstr) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if (envstr) fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if (envstr) fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if (envstr) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  fprintf(stderr, "Predefined: __SSE2__\n");
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "sizeof(size_t)      = %zu\n", sizeof(size_t));
  {
    constexpr size_t megaByte = 1024 * 1024;
    std::vector<int> numElements = { 1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025, 1 * megaByte };
    fprintf(stderr, "alignof malloc data = %d\n", alignof_malloc_data(numElements));
    fprintf(stderr, "alignof malloc big  = %d\n", alignof_malloc_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
    fprintf(stderr, "alignof vector data = %d\n", alignof_vector_data(numElements));
    fprintf(stderr, "alignof vector big  = %d\n", alignof_vector_data({ 8 * megaByte, 16 * megaByte, 32 * megaByte }));
  }
  fprintf(stderr, "\n");

#ifdef HAVE_MMAP
  fprintf(stderr, "HAVE_MMAP\n");
#endif
  fprintf(stderr, "\n");

#ifdef _OPENACC
  fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
  // OPENMP3:   201107
  // OPENMP4:   201307 gcc 4.9
  // OPENMP45:  201511
#ifdef _OPENMP
  fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
  fprintf(stderr, "__cplusplus         = %ld\n", (long) __cplusplus);
#ifdef __GNUC__
  fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#ifdef __GNUC_MINOR__
  fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#ifdef __ICC
  fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#ifdef __STDC__
  fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#ifdef __STD_VERSION__
  fprintf(stderr, "STD VERSION         = %ld\n", (long) __STD_VERSION__);
#endif
#ifdef __STDC_VERSION__
  fprintf(stderr, "STDC VERSION        = %ld\n", (long) __STDC_VERSION__);
#endif
#ifdef __STD_HOSTED__
  fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#ifdef FLT_EVAL_METHOD
  fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#ifdef FP_FAST_FMA
  fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#ifdef __FAST_MATH__
  fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  fprintf(stderr, "\n");

#ifdef _SC_VERSION
  fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#ifdef _SC_ARG_MAX
  fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#ifdef _SC_CHILD_MAX
  fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#ifdef _SC_STREAM_MAX
  fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#ifdef _SC_OPEN_MAX
  fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#ifdef _SC_PAGESIZE
  fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  fprintf(stderr, "\n");

  cdo::print_rlimits();

  fprintf(stderr, "\n");
}

static void
print_filetypes(FILE *fp)
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

  fprintf(fp, "CDI file types: ");
  set_text_color(fp, BRIGHT, GREEN);
  for (auto const &[type, name] : fileTypes)
    if (cdiHaveFiletype(type)) fprintf(fp, "%s ", name.c_str());
  reset_text_color(fp);
  fprintf(fp, "\n");
}

void
version()
{
  auto fp = stdout;
  fprintf(fp, "%s\n", cdo::Version);
#ifdef SYSTEM_TYPE
  fprintf(fp, "System: %s\n", SYSTEM_TYPE);
#endif
#ifdef CXX_COMPILER
  fprintf(fp, "CXX Compiler: %s\n", CXX_COMPILER);
#ifdef CXX_VERSION
  fprintf(fp, "CXX version : %s\n", CXX_VERSION);
#endif
  fprintf(fp, "CXX library :");
#ifdef HAVE_LIB_RANGES_ZIP
  fprintf(fp, " ranges_zip");
#endif
#ifdef HAVE_LIB_MDSPAN
  fprintf(fp, " mdspan");
#endif
  fprintf(fp, "\n");
#endif
#ifdef C_COMPILER
  fprintf(fp, "C Compiler: %s\n", C_COMPILER);
#ifdef C_VERSION
  fprintf(fp, "C version : %s\n", C_VERSION);
#endif
#endif
#ifdef F77_COMPILER
  fprintf(fp, "F77 Compiler: %s\n", F77_COMPILER);
#ifdef F77_VERSION
  fprintf(fp, "F77 version : %s\n", F77_VERSION);
#endif
#endif

  cdo::features::print_features();
  cdo::features::print_libraries();

#ifdef CDI_SIZE_TYPE
#define CDO_STRINGIFY(x) #x
#define CDO_TOSTRING(x) CDO_STRINGIFY(x)
  fprintf(fp, "CDI data types: SizeType=%s\n", CDO_TOSTRING(CDI_SIZE_TYPE));
#endif

  print_filetypes(fp);

  cdiPrintVersion();
  fprintf(fp, "\n");
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

    fprintf(stderr, "  User time:     %.3f seconds\n", ut);
    fprintf(stderr, "  System time:   %.3f seconds\n", st);
    fprintf(stderr, "  Total time:    %.3f seconds\n", ut + st);
#if defined(__APPLE__) && defined(__MACH__)
    fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024.0 * 1024.0));
#else
    fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024.0));
#endif
    fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
    fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
    fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
    fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
    fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
    /*
#if defined(__APPLE__) && defined(__MACH__)
#ifdef RUSAGE_INFO_CURRENT
    pid_t pid = fork();
    rusage_info_current ruinfo;
    int rusage_ret = proc_pid_rusage(pid, RUSAGE_INFO_CURRENT, (void **)&ruinfo);
    if (rusage_ret >= 0 && ruinfo.ri_lifetime_max_phys_footprint > 0)
      {
        fprintf(stderr, "  Peak memory:   %.2f MBytes\n", ruinfo.ri_lifetime_max_phys_footprint / 1.0);
        fprintf(stderr, "  Peak memory:   %.2f MBytes\n", ruinfo.ri_lifetime_max_phys_footprint / (1024.0 * 1024.0));
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
  fprintf(stderr, "OMP num procs       = %d\n", omp_get_num_procs());
  fprintf(stderr, "OMP max threads     = %d\n", omp_get_max_threads());
  fprintf(stderr, "OMP num threads     = %d\n", omp_get_num_threads());
#ifndef HAVE_OPENMP3
  fprintf(stderr, "OMP thread limit    = %d\n", omp_get_thread_limit());
  omp_sched_t kind;
  int modifer;
  omp_get_schedule(&kind, &modifer);
  fprintf(stderr, "OMP schedule        = %d (1:static; 2:dynamic; 3:guided; 4:auto)\n", (int) kind);
#endif
#ifdef HAVE_OPENMP4
  fprintf(stderr, "OMP proc bind       = %d (0:false; 1:true; 2:master; 3:close; 4:spread)\n", (int) omp_get_proc_bind());
#ifndef __ICC
  fprintf(stderr, "OMP num devices     = %d\n", omp_get_num_devices());
#endif
#endif
}
#endif

};  // namespace features
};  // namespace cdo
