#include <csignal>
#include <cstdlib>
#include "cdo_def_options.h"
#include "cdo_getopt.h"
#include "percentiles.h"
#include "cdo_options.h"
#include "cdo_default_values.h"
#include "util_string.h"
#include "cdo_features.h"
#include "griddes.h"
#include "cdo_output.h"
#include "param_conversion.h"
#include "cdo_settings.h"
#include "cdi.h"
#include "datetime.h"
#include "table.h"
#include "mpim_grid/mpim_grid.h"
#include "grid_pointsearch.h"
#include "institution.h"
#include "cdo_zaxis.h"
#include "chunkspec.h"

#ifdef HIRLAM_EXTENSIONS
extern "C" void streamGrbDefDataScanningMode(int scanmode);
#endif

static void
set_chunkspec_parameter(std::string const &argument)
{
  auto chunkSpec = cdo::parse_chunkspec_parameter(argument);
  if (chunkSpec.t) Options::cdoChunkSizeDimT = chunkSpec.t;
  if (chunkSpec.z) Options::cdoChunkSizeDimZ = chunkSpec.z;
  if (chunkSpec.y) Options::cdoChunkSizeDimY = chunkSpec.y;
  if (chunkSpec.x) Options::cdoChunkSizeDimX = chunkSpec.x;
}

void
setup_options()
{
  CLIOptions::option("envvars")
      ->add_effect([&]() { CLIOptions::print_envvars = true; })
      ->aborts_program(true)
      ->set_category("Info")
      ->add_help("Prints the environment variables of CDO.");

  CLIOptions::option("settings")
      ->add_effect([&]() { CLIOptions::print_settings = true; })
      ->aborts_program(true)
      ->set_category("Info")
      ->add_help("Prints the settings of CDO.");

  CLIOptions::option("debug")
      ->add_effect(
          [&]()
          {
            unsigned cdoDebugLevel = 0;
            unsigned cdiDebugLevel = 0;
            cdo::parse_debug_arguments({ "1" }, cdoDebugLevel, cdiDebugLevel);

            cdiDebug(cdiDebugLevel);
            cdo::set_debug(cdoDebugLevel);
            cdo::features::version();
          })
      ->set_category("Output")
      ->add_help("Pring all available debug messages")
      ->shortform('d');

  CLIOptions::option("scoped_debug")
      ->describe_argument("comma seperated scopes")
      ->set_category("Output")
      ->on_empty_argument(
          []()
          {
            std::cerr << "No debug level given please choose: " << std::endl;
            print_debug_options();
            std::exit(EXIT_SUCCESS);
          })
      ->add_effect(
          [&](std::string const &argument)
          {
            auto tokens = split_string(argument, ",");
            if (tokens.empty())
            {
              print_debug_options();
              std::exit(EXIT_SUCCESS);
            }

            unsigned cdoDebugLevel = 0;
            unsigned cdiDebugLevel = 0;
            cdo::parse_debug_arguments(tokens, cdoDebugLevel, cdiDebugLevel);

            cdiDebug(cdiDebugLevel);
            cdo::set_debug(cdoDebugLevel);

            cdo::features::version();
          })
      ->add_help("Multiple scopes simultaneously possible. Use this option without arguments to get a list of possible scopes")
      ->shortform('D');

  CLIOptions::option("worker")
      ->describe_argument("num")
      ->add_effect([&](std::string const &argument) { Options::numStreamWorker = parameter_to_int(argument); })
      ->set_category("Multi Threading")
      ->add_help("Number of worker to decode/decompress GRIB records.");

  CLIOptions::option("precision")
      ->describe_argument("float_digits[,double_digits]")
      ->add_effect([&](std::string const &argument) { cdo::set_digits(argument); })
      ->set_category("Numeric")
      ->add_help("Precision to use in displaying floating-point data (default: 7,15).");

  CLIOptions::option("percentile")
      ->describe_argument("method")
      ->set_category("Numeric")
      ->add_effect([&](std::string const &argument) { percentile_set_method(argument); })
      ->add_help("Methods: nrank, nist, rtype8, <NumPy method (linear|lower|higher|nearest|...)>");

  CLIOptions::option("netcdf_hdr_pad")
      ->describe_argument("nbr")
      ->add_effect(
          [&](std::string const &argument)
          {
            int netcdf_hdr_pad = parameter_to_bytes(argument);
            if (netcdf_hdr_pad >= 0) cdo::netcdf_hdr_pad = netcdf_hdr_pad;
          })
      ->add_help("Pad NetCDF output header with nbr bytes.");

  CLIOptions::option("use_fftw")
      ->describe_argument("false|true")
      ->add_effect([&](std::string const &argument) { Options::Use_FFTW = (int) parameter_to_bool(argument); })
      ->add_help("Sets fftw usage.");

  CLIOptions::option("config")
      ->describe_argument("all|all-json|<specific_feature_name>")
      ->add_effect([&](std::string const &argument) { cdo::features::print_config(argument); })
      ->on_empty_argument([&]() { cdo::features::print_argument_options(); })
      ->aborts_program(true)
      ->set_category("Info")
      ->add_help("Prints all features and the enabled status.", "Use option <all> to see explicit feature names.");

  CLIOptions::option("pointsearchmethod")
      ->set_internal(true)
      ->describe_argument("<kdtree|nanoflann|spherepart|full>")
      ->set_category("Search Methods")
      ->add_effect([&](std::string const &argument) { set_pointsearch_method(argument); })
      ->add_help("Sets the point search method.");

  CLIOptions::option("gridsearchradius")
      ->describe_argument("degrees[0..180]")
      ->set_category("Search Methods")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto fval = radius_str_to_deg(argument);
            if (fval < 0 || fval > 180) cdo_abort("%s=%g out of bounds (0-180 deg)!", "gridsearchradius", fval);
            cdo_set_search_radius(fval);
          })
      ->add_help("Sets the grid search radius (0-180 deg).");

  CLIOptions::option("remap_weights")
      ->describe_argument("false|true")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto intarg = parameter_to_bool(argument);
            if (intarg != 0 && intarg != 1) cdo_abort("Unsupported value for option --remap_weights %d [false|true]", intarg);
            Options::RemapGenerateWeights = intarg;
          })
      ->add_help("Generate remap weights (default: 1).");

  CLIOptions::option("no_remap_weights")
      ->add_effect([&]() { Options::RemapGenerateWeights = 0; })
      ->add_help("Switch off generation of remap weights.");

  CLIOptions::option("enableexcept")
      ->describe_argument("except")
      ->set_category("Numeric")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto except = cdo::evaluate_except_options(argument);
            if (except < 0) cdo_abort("option --%s: unsupported argument: %s", "enableexcept", argument);
            cdo::set_feenableexcept(except);
            if (signal(SIGFPE, cdo::signal_handler) == SIG_ERR) cdo_warning("can't catch SIGFPE!");
          })
      ->add_help("Set individual floating-point traps ", "(DIVBYZERO, INEXACT, INVALID, OVERFLOW, UNDERFLOW, ALL_EXCEPT)");

  CLIOptions::option("timestat_date")
      ->describe_argument("srcdate")
      ->add_effect([&](std::string const &argument) { set_timestat_date(argument); })
      ->add_help("Target timestamp (temporal statistics): ", "first, middle, midhigh or last source timestep.");

  CLIOptions::option("ignore_time_bounds")
      ->add_effect(
          [&]()
          {
            extern bool CDO_Ignore_Time_Bounds;
            CDO_Ignore_Time_Bounds = true;
          })
      ->add_help("Ignores time bounds for time range statistics.");

  CLIOptions::option("use_time_bounds")
      ->add_effect(
          [&]()
          {
            extern bool CDO_Use_Time_Bounds;
            CDO_Use_Time_Bounds = true;
          })
      ->add_help("Enables use of timebounds.");

  CLIOptions::option("cmor")->add_effect([&]() { Options::CMOR_Mode = 1; })->add_help("CMOR conform NetCDF output.");

  CLIOptions::option("reduce_dim")->add_effect([&]() { Options::CDO_Reduce_Dim = 1; })->add_help("Reduce NetCDF dimensions.");

  CLIOptions::option("float")
      ->set_internal(true)
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Float; })
      ->set_category("Numeric")
      ->add_help("Uses single precision floats for reading data.");

  CLIOptions::option("single")
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Float; })
      ->set_category("Numeric")
      ->add_help("Uses single precision floats for reading data.");

  CLIOptions::option("double")
      ->add_effect([&]() { Options::CDO_Memtype = MemType::Double; })
      ->set_category("Numeric")
      ->add_help("Uses double precision floats for reading data.");

  CLIOptions::option("rusage")
      ->add_effect([&]() { Options::CDO_Rusage = 1; })
      ->add_help("Print information about resource utilization.")
      ->set_category("Info");

  CLIOptions::option("pedantic")->add_effect([&]() { MpMO::enable_pedantic(true); })->add_help("Warnings count as errors.");

  CLIOptions::option("eccodes")
      ->add_effect([&]() { cdiDefGlobal("ECCODES_GRIB1", true); })
      ->set_category("Format Specific")
      ->add_help("Use ecCodes to decode/encode GRIB1 messages.");

  CLIOptions::option("format")
      ->describe_argument("grb1|grb2|nc1|nc2|nc4|nc4c|nc5|nczarr|srv|ext|ieg")
      ->add_effect([&](std::string const &argument) { cdo::set_default_filetype(argument); })
      ->add_help("Format of the output file.")
      ->shortform('f');

  CLIOptions::option("history")
      ->add_effect([&]() { Options::CDO_Append_History = true; })
      ->set_category("History")
      ->add_help("Do append to NetCDF \"history\" global attribute.");

  CLIOptions::option("no_history")
      ->add_effect([&]() { Options::CDO_Append_History = false; })
      ->set_category("History")
      ->add_help("Do not append to NetCDF \"history\" global attribute.");

  CLIOptions::option("version")
      ->add_effect([&]() { cdo::features::version(); })
      ->aborts_program(true)
      ->set_category("Info")
      ->add_help("Print the version number.")
      ->shortform('V');

  CLIOptions::option("absolute_taxis")
      ->add_effect(
          [&]()
          {
            if (CdoDefault::TaxisType == TAXIS_RELATIVE)
              cdo_abort("option --%s: can't be combined with option --%s", "absolute_taxis (-a)", "relative_taxis (-r)");
            CdoDefault::TaxisType = TAXIS_ABSOLUTE;
          })
      ->add_help("Generate an absolute time axis.")
      ->shortform('a');

  CLIOptions::option("force")->add_effect([&]() { Options::force = true; })->add_help("Forcing a CDO process.");

  CLIOptions::option("fast")
      ->set_internal(true)
      ->add_effect(
          [&]()
          {
            Options::fast = true;
            Options::lazyGridLoad = true;
            cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true);
          })
      ->add_help("If available, use a faster method even if it requires more memory.");

  CLIOptions::option("lazy_grid_load")
      ->set_internal(true)
      ->describe_argument("false|true")
      ->add_effect([&](std::string const &argument) { Options::lazyGridLoad = parameter_to_bool(argument); })
      ->add_help("Enable/disable lazy grid load");

  // clang-format off
  CLIOptions::option("default_datatype")
      ->describe_argument("nbits")
      ->set_category("Numeric")
      ->add_effect([&](std::string const &argument) { cdo::set_default_datatype(argument); })
      ->add_help("Set the number of bits for the output precision",
                 "    I8|I16|I32|F32|F64     for nc1,nc2,nc4,nc4c,nc5,nczarr;",
                 "    U8|U16|U32             for nc4,nc4c,nc5;",
                 "    F32|F64                for grb2,srv,ext,ieg;",
                 "    P1 - P24               for grb1,grb2")->shortform('b');
  // clang-format on

  CLIOptions::option("check_data_range")
      ->add_effect([&]() { Options::CheckDatarange = true; })
      ->add_help("Enables checks for data overflow.")
      ->shortform('c');

  CLIOptions::option("grid")
      ->describe_argument("grid")
      ->add_effect([&](std::string const &argument) { cdo_set_grids(argument); })
      ->add_help("Set default grid name or file. Available grids: ",
                 "global_<DXY>, zonal_<DY>, r<NX>x<NY>, lon=<LON>/lat=<LAT>, F<N>, gme<NI>, hpz<ZOOM>")
      ->shortform('g');

  CLIOptions::option("institution")
      ->describe_argument("institute_name")
      ->add_effect([&](std::string const &argument) { define_institution(argument); })
      ->add_help("Sets institution name.")
      ->shortform('i');

  CLIOptions::option("chunktype")
      ->describe_argument("auto|grid|lines")
      ->set_category("Format Specific")
      ->add_effect([&](std::string const &argument) { cdo::set_chunktype(argument); })
      ->add_help("NetCDF4 chunk type (x/y dimension).")
      ->shortform('k');

  CLIOptions::option("chunksize")
      ->describe_argument("size")
      ->set_category("Format Specific")
      ->add_effect(
          [&](std::string const &argument)
          {
            int chunkSize = parameter_to_bytes(argument);
            if (chunkSize >= 0) Options::cdoChunkSize = chunkSize;
          })
      ->add_help("NetCDF4 chunk size (x/y dimension).");

  CLIOptions::option("chunkspec")
      ->describe_argument("chunkspec")
      ->set_category("Format Specific")
      ->add_effect([&](std::string const &argument) { set_chunkspec_parameter(argument); })
      ->add_help("NetCDF4 specify chunking for dimensions (x,y,z,t).");

  CLIOptions::option("copy_chunkspec")
      ->set_category("Format Specific")
      ->add_effect([&]() { cdiDefGlobal("COPY_CHUNKSPEC", true); })
      ->add_help("Copy chunk specification.");

  CLIOptions::option("remove_chunkspec")
      ->set_category("Format Specific")
      ->add_effect([&]() { cdiDefGlobal("REMOVE_CHUNKSPEC", true); })
      ->add_help("Remove chunk specification.");

  CLIOptions::option("lock_io")
      ->set_internal(true)
      ->add_effect([&]() { Threading::cdoLockIO = true; })
      ->add_help("Lock IO (sequential access).")
      ->shortform('L');

  CLIOptions::option("zaxis")
      ->describe_argument("zaxis")
      ->add_effect([&](std::string const &argument) { cdo_set_zaxes(argument); })
      ->add_help("Set default zaxis name or file.")
      ->shortform('l');

  CLIOptions::option("set_missval")
      ->describe_argument("missval")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto [success, mv] = string_to_floating<double>(argument);
            if (success)
            {
              Debug("set missval of cdi to: %f", mv);
              cdiDefMissval(mv);
            }
            else { cdo_abort("Could not convert %s to double", argument); }
          })
      ->add_help("Set the missing value of non NetCDF files (default: " + get_scientific(cdiInqMissval()) + ").")
      ->shortform('m');

  CLIOptions::option("has_missval")
      ->add_effect([&]() { cdiDefGlobal("HAVE_MISSVAL", true); })
      ->add_help("Set HAS_MISSVAL to true.")
      ->shortform('M');

  CLIOptions::option("query")
      ->set_internal(true)
      ->describe_argument("<query params>")
      ->add_effect([&](std::string const &argument) { Options::cdoQueryParameter = argument; })
      ->add_help("Set query parameter.");

  CLIOptions::option("varnames")
      ->set_internal(true)
      ->describe_argument("<varname| file>")
      ->add_effect([&](std::string const &argument) { Options::cdoVarnames = split_string(argument, ","); })
      ->add_help("Set default varnames or file.")
      ->shortform('n');

  CLIOptions::option("num_threads")
      ->describe_argument("nthreads")
      ->add_effect([&](std::string const &argument) { Threading::ompNumUserRequestedThreads = parameter_to_int(argument); })
      ->set_category("Multi Threading")
      ->add_help("Set number of OpenMP threads.")
      ->shortform('P');

  CLIOptions::option("async_read")
      ->describe_argument("true|false")
      ->add_effect(
          [&](std::string const &argument)
          {
            Options::CDO_Async_Read = parameter_to_bool(argument);
            Options::CDO_task = Options::CDO_Async_Read;
          })
      ->set_category("Multi Threading")
      ->add_help("Read input data asynchronously [default: false].",
                 "Available for the operators: diff, info, trend, detrend, Timstat");

  CLIOptions::option("p")
      ->add_effect(
          [&]()
          {
            Options::CDO_Async_Read = true;
            Options::CDO_task = true;
          })
      ->set_category("Multi Threading")
      ->add_help("Read input data asynchronously, short for '--async_read true'")
      ->shortform('p');

  CLIOptions::option("sortname")
      ->add_effect([&]() { cdiDefGlobal("SORTNAME", true); })
      ->set_category("Format Specific")
      ->add_help("Alphanumeric sorting of NetCDF parameter names.")
      ->shortform('Q');

  CLIOptions::option("seed")
      ->describe_argument("seed")
      ->set_category("Numeric")
      ->add_effect(
          [&](std::string const &argument)
          {
            int intarg = parameter_to_int(argument);
            if (intarg < 0) cdo_abort("Unsupported value for option --seed %d [>=0]", intarg);
            Options::Random_Seed = intarg;
          })
      ->add_help("Seed for a new sequence of pseudo-random numbers. <seed> must be >= 0");

  CLIOptions::option("regular")
      ->add_effect(
          [&]()
          {
            Options::cdoRegulargrid = true;
            cdiDefGlobal("REGULARGRID", true);
          })
      ->set_category("CGRIBEX")
      ->add_help("Convert GRIB1 data from global reduced to regular Gaussian grid (cgribex only).")
      ->shortform('R');

  CLIOptions::option("relative_taxis")
      ->add_effect(
          [&]()
          {
            if (CdoDefault::TaxisType == TAXIS_ABSOLUTE)
              cdo_abort("option --%s: can't be combined with option --%s", "relative_taxis (-r)", "absolute_taxis (-a)");
            CdoDefault::TaxisType = TAXIS_RELATIVE;
          })
      ->add_help("Generate a relative time axis.")
      ->shortform('r');

  CLIOptions::option("diagnostic")
      ->add_effect([&]() { Options::CDO_diagnostic = true; })
      ->add_help("Create a diagnostic output stream for the module TIMSTAT. This stream",
                 "contains the number of non missing values for each output period.")
      ->shortform('S');

  CLIOptions::option("silent")
      ->add_effect(
          [&]()
          {
            Options::silentMode = true;
            MpMO::enable_silent_mode(Options::silentMode);
          })
      ->set_category("Output")
      ->add_help("Silent mode.")
      ->shortform('s');

  CLIOptions::option("timer")->add_effect([&]() { Options::Timer = true; })->add_help("Enable timer.")->shortform('T');

  CLIOptions::option("table")
      ->describe_argument("codetab")
      ->set_category("CGRIBEX")
      ->add_effect([&](std::string const &argument) { CdoDefault::TableID = cdo::define_table(argument); })
      ->add_help("Set GRIB1 default parameter code table name or file (cgribex only).", cdo::predefined_tables(CLIOptions::padding))
      ->shortform('t');

  CLIOptions::option("sortparam")->add_effect([]() { cdiDefGlobal("SORTPARAM", true); });

  CLIOptions::option("print_filename")
      ->add_effect([]() { Options::PrintFilename = true; })
      ->add_help("Print name of all output files.");

  CLIOptions::option("verbose")
      ->add_effect(
          [&]()
          {
            Options::cdoVerbose = true;
            MpMO::enable_verbose(true);
            CLIOptions::print_envvars = true;
            gridEnableVerbose(Options::cdoVerbose);
          })
      ->add_help("Print extra details for some operators.")
      ->shortform('v');

  CLIOptions::option("disable_warnings")
      ->add_effect(
          [&]() {  // disable warning messages
            MpMO::enable_warnings(false);
            extern int _Verbose;  // CDI Warnings
            _Verbose = 0;
          })
      ->set_category("Output")
      ->add_help("Disable warning messages.")
      ->shortform('w');

  CLIOptions::option("par_io")
      ->set_internal(true)
      ->add_effect(
          [&]()
          {
            Options::cdoParIO = true;  // multi threaded I/O
          })
      ->add_help("Enables multithreaded I/O.")
      ->set_category("Multi Threading")
      ->shortform('X');

  CLIOptions::option("shuffle")
      ->add_effect([&]() { Options::cdoShuffle = true; })
      ->set_category("Compression")
      ->add_help("Specify shuffling of variable data bytes before compression (NetCDF)");

  CLIOptions::option("compress")
      ->add_effect([&]() { Options::cdoCompress = true; })
      ->set_category("Compression")
      ->add_help("Enables compression. Default = SZIP")
      ->shortform('Z');

  CLIOptions::option("filter")
      ->describe_argument("filterspec")
      ->add_effect([&](std::string const &argument) { cdo::set_filterspec(argument); })
      ->set_category("Compression")
      ->add_help("NetCDF4 filter specification")
      ->shortform('F');

  CLIOptions::option("compression_type")
      ->describe_argument("aec|jpeg|zip[_1-9]|zstd[1-19]")
      ->set_category("Compression")
      ->add_effect([&](std::string const &argument) { cdo::set_compression_type(argument); })
      ->add_help("aec         AEC compression of GRIB2 records", "jpeg        JPEG compression of GRIB2 records",
                 "zip[_1-9]   Deflate compression of NetCDF4 variables", "zstd[_1-19] Zstandard compression of NetCDF4 variables")
      ->shortform('z');

  CLIOptions::option("nsb")
      ->set_internal(true)
      ->describe_argument("1-23")
      ->add_effect([&](std::string const &argument) { Options::nsb = parameter_to_int(argument); })
      ->set_category("Numeric")
      ->add_help("Number of significant bits used for bit-rounding.");

  CLIOptions::option("show_available_options")
      ->set_internal(true)
      ->aborts_program(true)
      ->set_category("Info")
      ->add_effect([&]() { CLIOptions::print_available_options(); })
      ->add_help("Shows all available optins and prints all shortforms, only internal use for testing.");

#ifdef HIRLAM_EXTENSIONS
  CLIOptions::option("Dkext")
      ->describe_argument("debLev")
      ->set_category("Hirlam Extension")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto extDebugVal = parameter_to_int(argument);
            if (extDebugVal > 0)
            {
              extern int cdiDebugExt;
              cdoDebugExt = extDebugVal;
              cdiDebugExt = extDebugVal;
            }
          })
      ->add_help("Setting debugLevel for extensions.");

  CLIOptions::option("outputGribDataScanningMode")
      ->describe_argument("mode")
      ->set_category("Hirlam Extension")
      ->add_effect(
          [&](std::string const &argument)
          {
            auto scanningModeValue = parameter_to_int(argument);
            if (cdoDebugExt) printf("scanningModeValue=%d\n", scanningModeValue);

            if ((scanningModeValue == 0) || (scanningModeValue == 64) || (scanningModeValue == 96))
            {
              streamGrbDefDataScanningMode(scanningModeValue);  // -1: not used; allowed modes: <0,
                                                                // 64, 96>; Default is 64
            }
            else
            {
              cdo_warning("Warning: %d not in allowed modes: <0, 64, 96>; Using default: 64\n", scanningModeValue);
              streamGrbDefDataScanningMode(64);
            }
          })
      ->add_help("Setting grib scanning mode for data in output file <0, 64, 96>.", "Default is 64");
#endif  // HIRLAM_EXTENSIONS
}
