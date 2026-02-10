/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Cedrick Ansorge
          Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Filter    highpass
      Filter    lowpass
      Filter    bandpass
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#include <mutex>
static std::mutex fftwMutex;
#endif

#include "cdi.h"
#include "julian_date.h"
#include "process_int.h"
#include "param_conversion.h"
#include "cdo_fft.h"
#include "cdo_fftw3.h"
#include "cdo_options.h"
#include "datetime.h"
#include "cdo_omp.h"
#include "field_functions.h"

size_t
vars_numArrays(const FieldVector3D &vars)
{
  size_t numArrays{ 0 };
  for (size_t i = 0; i < vars[0].size(); ++i)
  {
    for (size_t j = 0; j < vars[0][i].size(); ++j) numArrays++;
  }

  return numArrays;
}

size_t
vars_allocatedMem(const FieldVector3D &vars)
{
  size_t allocatedMem{ 0 };
  auto func = [](auto &v) { return (v.size() > 0) ? v.size() * sizeof(v[0]) : 0; };
  for (size_t i = 0; i < vars[0].size(); ++i)
  {
    for (size_t j = 0; j < vars[0][i].size(); ++j) { allocatedMem += field_operation(func, vars[0][i][j]); }
  }

  return allocatedMem;
}

static void
create_fmasc(int nts, double fdata, double fmin, double fmax, std::vector<int> &fmasc)
{
  auto dimin = nts * fmin / fdata;
  auto dimax = nts * fmax / fdata;

  auto imin = (dimin < 0) ? 0 : (int) std::floor(dimin);
  auto imax = (std::ceil(dimax) > nts / 2) ? nts / 2 : (int) std::ceil(dimax);

  if (imin < 0 || imin >= nts) cdo_abort("Parameter fmin=%g: timestep %d out of bounds (1-%d)!", fmin, imin + 1, nts);
  if (imax < 0 || imax >= nts) cdo_abort("Parameter fmax=%g: timestep %d out of bounds (1-%d)!", fmax, imax + 1, nts);

  fmasc[imin] = 1;
  for (int i = imin + 1; i <= imax; ++i) fmasc[i] = fmasc[nts - i] = 1;
}

static void
filter_intrinsic(int nts, std::vector<int> const &fmasc, double *real, double *imag)
{
  auto isPowerOf2 = (nts > 0 && (nts & (nts - 1)) == 0);

  Varray<double> work_r, work_i;

  if (!isPowerOf2) work_r.resize(nts);
  if (!isPowerOf2) work_i.resize(nts);

  if (isPowerOf2)
    cdo::fft(real, imag, nts, 1);
  else
    cdo::ft_r(real, imag, nts, 1, work_r.data(), work_i.data());

  for (int i = 0; i < nts; ++i)
    if (!fmasc[i]) real[i] = imag[i] = 0;

  if (isPowerOf2)
    cdo::fft(real, imag, nts, -1);
  else
    cdo::ft_r(real, imag, nts, -1, work_r.data(), work_i.data());

  return;
}

namespace
{
struct FilterMemory
{
  Varray<double> real;
  Varray<double> imag;
#ifdef HAVE_LIBFFTW3
  fftw_complex *in_fft{};
  fftw_complex *out_fft{};
  fftw_plan p_T2S;
  fftw_plan p_S2T;
#endif
};
}  // namespace

class Filter : public Process
{
  enum
  {
    BANDPASS,
    HIGHPASS,
    LOWPASS
  };

public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Filter",
    // clang-format off
    .operators = { { "bandpass", BANDPASS, 0, FilterHelp },
                   { "highpass", HIGHPASS, 0, FilterHelp },
                   { "lowpass", LOWPASS, 0, FilterHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Filter> registration = RegisterEntry<Filter>(module);

  std::vector<std::string> tunits{};
  std::vector<int> iunits{};
  int year0{}, month0{}, day0{};
  double fdata{ 0.0 };
  TimeIncrement timeIncr0 = { 0, TimeUnits::SECONDS };
  DateTimeList dtlist{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1{};

  bool useFFTW{ false };
  int operfunc{};
  int calendar{};

public:
  void
  init() override
  {
    tunits = { "second", "minute", "hour", "day", "month", "year" };
    iunits = { 31536000, 525600, 8760, 365, 12, 1 };

    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    if (Options::Use_FFTW)
    {
#ifdef HAVE_LIBFFTW3
      if (Options::cdoVerbose) cdo_print("Using fftw3 lib");
      useFFTW = true;
#else
      if (Options::cdoVerbose) cdo_print("LIBFFTW3 support not compiled in!");
#endif
    }

    if (Options::cdoVerbose && !useFFTW) cdo_print("Using intrinsic FFT function!");

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    calendar = taxisInqCalendar(taxisID1);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector3D varsData;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

      field2D_init(varsData[tsID], varList1);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field = varsData[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
        if (field.numMissVals) cdo_abort("Missing value support for operators in module Filter not added yet!");
      }

      // get and check time increment
      if (tsID > 0)
      {
        auto vDateTime0 = dtlist.vDateTime(tsID - 1);
        auto vDateTime = dtlist.vDateTime(tsID);

        cdiDate_decode(vDateTime0.date, &year0, &month0, &day0);
        int year, month, day;
        cdiDate_decode(vDateTime.date, &year, &month, &day);

        auto julianDate0 = julianDate_encode(calendar, vDateTime0);
        auto julianDate = julianDate_encode(calendar, vDateTime);
        auto jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));

        auto timeIncr = get_time_increment(jdelta, vDateTime0.date, vDateTime.date);

        if (tsID == 1)
        {
          timeIncr0 = timeIncr;
          if (timeIncr.period == 0) cdo_abort("Time step must be different from zero!");
          if (Options::cdoVerbose) cdo_print("Time step %lld %s", timeIncr.period, tunits[(int) timeIncr.units]);
          fdata = 1.0 * iunits[(int) timeIncr.units] / timeIncr.period;
        }

        if (calendar != CALENDAR_360DAYS && calendar != CALENDAR_365DAYS && calendar != CALENDAR_366DAYS
            && timeIncr0.units < TimeUnits::MONTHS && month == 2 && day == 29 && (day0 != day || month0 != month || year0 != year))
        {
          cdo_warning("Filtering of multi-year times series doesn't works properly with a standard calendar.");
          cdo_warning("  Please delete the day %i-02-29 (cdo del29feb)", year);
        }

        if (timeIncr.period != timeIncr0.period || timeIncr.units != timeIncr0.units)
          cdo_warning("Time increment in step %d (%lld%s) differs from step 1 (%lld%s)!", tsID + 1, timeIncr.period,
                      tunits[(int) timeIncr.units], timeIncr0.period, tunits[(int) timeIncr0.units]);
      }

      tsID++;
    }

    auto nts = tsID;
    if (nts <= 1) cdo_abort("Number of time steps <= 1!");

    auto numArrays = vars_numArrays(varsData);
    auto allocatedMem = vars_allocatedMem(varsData) * nts;
    if (Options::cdoVerbose)
      cdo_print("Allocate %zu array%s over %zu steps: size=%zu Bytes", numArrays, numArrays > 1 ? "s" : "", nts, allocatedMem);

    std::vector<FilterMemory> fourierMemory(Threading::ompNumMaxThreads);

    if (useFFTW)
    {
#ifdef HAVE_LIBFFTW3
      for (auto &fm : fourierMemory)
      {
        fm.in_fft = fftw_alloc_complex(nts);
        fm.out_fft = fftw_alloc_complex(nts);
        std::scoped_lock lock(fftwMutex);
        fm.p_T2S = fftw_plan_dft_1d(nts, fm.in_fft, fm.out_fft, 1, FFTW_ESTIMATE);
        fm.p_S2T = fftw_plan_dft_1d(nts, fm.out_fft, fm.in_fft, -1, FFTW_ESTIMATE);
      }
#endif
    }
    else
    {
      for (auto &fm : fourierMemory)
      {
        fm.real.resize(nts);
        fm.imag.resize(nts);
      }
    }

    double fmin = 0.0, fmax = 0.0;
    switch (operfunc)
    {
      case BANDPASS:
      {
        operator_input_arg("lower and upper bound of frequency band");
        operator_check_argc(2);
        fmin = parameter_to_double(cdo_operator_argv(0));
        fmax = parameter_to_double(cdo_operator_argv(1));
        break;
      }
      case HIGHPASS:
      {
        operator_input_arg("lower bound of frequency pass");
        operator_check_argc(1);
        fmin = parameter_to_double(cdo_operator_argv(0));
        fmax = fdata;
        break;
      }
      case LOWPASS:
      {
        operator_input_arg("upper bound of frequency pass");
        operator_check_argc(1);
        fmin = 0.0;
        fmax = parameter_to_double(cdo_operator_argv(0));
        break;
      }
    }

    if (Options::cdoVerbose) cdo_print("fmin=%g  fmax=%g", fmin, fmax);

    std::vector<int> fmasc(nts, 0);
    create_fmasc(nts, fdata, fmin, fmax, fmasc);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        if (useFFTW)
        {
#ifdef HAVE_LIBFFTW3
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < var.gridsize; ++i)
          {
            auto ompthID = cdo_omp_get_thread_num();
            auto &fm = fourierMemory[ompthID];

            if (var.memType == MemType::Float)
              for (int t = 0; t < nts; ++t)
              {
                fm.in_fft[t][0] = varsData[t][varID][levelID].vec_f[i];
                fm.in_fft[t][1] = 0.0;
              }
            else
              for (int t = 0; t < nts; ++t)
              {
                fm.in_fft[t][0] = varsData[t][varID][levelID].vec_d[i];
                fm.in_fft[t][1] = 0.0;
              }

            filter_fftw(nts, fmasc, fm.out_fft, &fm.p_T2S, &fm.p_S2T);

            if (var.memType == MemType::Float)
              for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_f[i] = fm.in_fft[t][0] / nts;
            else
              for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_d[i] = fm.in_fft[t][0] / nts;
          }
#endif
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
          for (size_t i = 0; i < var.gridsize; ++i)
          {
            auto ompthID = cdo_omp_get_thread_num();
            auto &fm = fourierMemory[ompthID];

            if (var.memType == MemType::Float)
              for (int t = 0; t < nts; ++t) fm.real[t] = varsData[t][varID][levelID].vec_f[i];
            else
              for (int t = 0; t < nts; ++t) fm.real[t] = varsData[t][varID][levelID].vec_d[i];

            std::ranges::fill(fm.imag, 0.0);

            filter_intrinsic(nts, fmasc, fm.real.data(), fm.imag.data());

            if (var.memType == MemType::Float)
              for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_f[i] = fm.real[t];
            else
              for (int t = 0; t < nts; ++t) varsData[t][varID][levelID].vec_d[i] = fm.real[t];
          }
        }
      }
    }

#ifdef HAVE_LIBFFTW3
    if (useFFTW)
    {
      for (auto &fm : fourierMemory)
      {
        fftw_free(fm.in_fft);
        fftw_free(fm.out_fft);
        std::scoped_lock lock(fftwMutex);
        fftw_destroy_plan(fm.p_T2S);
        fftw_destroy_plan(fm.p_S2T);
      }
      fftw_cleanup();
    }
#endif

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    for (tsID = 0; tsID < nts; ++tsID)
    {
      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto &field = varsData[tsID][varID][levelID];
          if (field.hasData())
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field);
          }
        }
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
