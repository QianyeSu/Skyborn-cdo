/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Spectrum    spectrum         Spectrum
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_fft.h"
#include "field_functions.h"

static void
spectrum(int nrec, double *data, double *spectrum, double *real, double *imag, const double *window, double wssum, int detrend,
         int seg_n, int seg_l)
{
  int k;
  double sumx, sumkx;
  double a, b;
  int seg_i, offset;
  int bit;

  for (bit = seg_l; !(bit & 1); bit >>= 1)
    ;

  if (detrend == 1)
  {
    sumx = 0;
    for (k = 0; k < nrec; ++k) sumx += data[k];
    sumx /= nrec;
    for (k = 0; k < nrec; ++k) data[k] -= sumx;
  }
  else if (detrend == 2)
  {
    sumx = sumkx = 0;
    for (k = 0; k < nrec; ++k)
    {
      sumx += data[k];
      sumkx += k * data[k];
    }
    b = (sumkx - sumx * (nrec - 1) / 2.) / ((nrec + 1) * nrec * (nrec - 1) / 12.);
    a = sumx / nrec - b * (nrec - 1) / 2.;
    for (k = 0; k < nrec; ++k) data[k] -= a + b * k;
  }

  Varray<double> work_r, work_i;
  if (bit != 1)
  {
    work_r.resize(seg_l);
    work_i.resize(seg_l);
  }

  for (seg_i = 0; seg_i < seg_n; seg_i += 2)
  {
    offset = (seg_n == 1) ? 0 : (int) ((double) (nrec - seg_l) / (seg_n - 1) * seg_i);

    for (k = 0; k < seg_l; ++k) real[k] = data[offset + k];

    if (detrend == 3)
    {
      sumx = sumkx = 0;
      for (k = 0; k < seg_l; ++k)
      {
        sumx += real[k];
        sumkx += k * real[k];
      }

      b = (sumkx - sumx * (seg_l - 1) / 2.) / ((seg_l + 1) * seg_l * (seg_l - 1) / 12.);
      a = sumx / seg_l - b * (seg_l - 1) / 2.;

      for (k = 0; k < seg_l; ++k) real[k] -= a + b * k;
    }

    for (k = 0; k < seg_l; ++k) real[k] *= window[k];

    if (seg_i + 1 < seg_n)
    {
      offset = (seg_n == 1) ? 0 : (int) ((double) (nrec - seg_l) / (seg_n - 1) * (seg_i + 1));
      for (k = 0; k < seg_l; ++k) imag[k] = data[offset + k];
      if (detrend == 3)
      {
        sumx = sumkx = 0;
        for (k = 0; k < seg_l; ++k)
        {
          sumx += imag[k];
          sumkx += k * imag[k];
        }

        b = (sumkx - sumx * (seg_l - 1) / 2.) / ((seg_l + 1) * seg_l * (seg_l - 1) / 12.);
        a = sumx / seg_l - b * (seg_l - 1) / 2.;

        for (k = 0; k < seg_l; ++k) imag[k] -= a + b * k;
      }

      for (k = 0; k < seg_l; ++k) imag[k] *= window[k];
    }
    else
      for (k = 0; k < seg_l; ++k) imag[k] = 0;

    if (bit == 1)  // seg_l is a power of 2
      cdo::fft(real, imag, seg_l, 1);
    else
      cdo::ft_r(real, imag, seg_l, 1, work_r.data(), work_i.data());

    spectrum[0] += real[0] * real[0] + imag[0] * imag[0];

    for (k = 1; k < (seg_l + 1) / 2; ++k)
      spectrum[k] += real[k] * real[k] + imag[k] * imag[k] + real[seg_l - k] * real[seg_l - k] + imag[seg_l - k] * imag[seg_l - k];

    if (!(seg_l & 1)) spectrum[seg_l / 2] += real[seg_l / 2] * real[seg_l / 2] + imag[seg_l / 2] * imag[seg_l / 2];
  }

  for (k = 0; k <= seg_l / 2; ++k) spectrum[k] *= seg_l / (seg_n * wssum);

  spectrum[0] *= 2;

  if (!(seg_l & 1)) spectrum[seg_l / 2] *= 2;
}

class Spectrum : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Spectrum",
    .operators = { { "spectrum" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Spectrum> registration = RegisterEntry<Spectrum>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };

  VarList varList1;

public:
  void
  init() override
  {
    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    FieldVector3D varsData;
    std::vector<CdiDateTime> vDateTimes;

    int nalloc = 0;
    auto numVars = varList1.numVars();

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (tsID >= nalloc)
      {
        constexpr int NALLOC_INC = 1024;
        nalloc += NALLOC_INC;
        vDateTimes.resize(nalloc);
        varsData.resize(nalloc);
      }

      vDateTimes[tsID] = taxisInqVdatetime(taxisID1);

      field2D_init(varsData[tsID], varList1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];
        auto &varData = varsData[tsID][varID][levelID];
        varData.resize(var.gridsize);
        cdo_read_field(streamID1, varData);

        if (varData.numMissVals) cdo_abort("Missing values are not allowed!");
      }

      tsID++;
    }

    int nts = tsID;

    operator_input_arg("detrend type, length of segments, number of segments, window type\n\n"
                       "       detrend type: 0 - data should be used unchanged\n"
                       "                     1 - the mean of the whole time series should be subtracted\n"
                       "                     2 - the whole time series should be detrended\n"
                       "                     3 - every segment should be detrended\n\n"
                       "        window type: 0 - no data windowing\n"
                       "                     1 - Hann window\n"
                       "                     2 - Bartlett window\n"
                       "                     3 - Welch window\n");

    operator_check_argc(4);

    auto detrend = parameter_to_int(cdo_operator_argv(0));
    auto seg_l = parameter_to_int(cdo_operator_argv(1));
    auto seg_n = parameter_to_int(cdo_operator_argv(2));
    auto which_window = parameter_to_int(cdo_operator_argv(3));

    if (detrend < 0 || detrend > 3) cdo_abort("Illegal value for detrend (=%d)!", detrend);

    if (seg_l <= 2 || seg_l > nts) cdo_abort("Length must be at least 3 and at most the number of timesteps (=%d)", nts);

    if (seg_n <= 0 || seg_n > nts - seg_l + 1)
      cdo_abort("Number of segments must be positiv and not greater than %d!", nts - seg_l + 1);

    int nfreq = seg_l / 2 + 1;

    FieldVector3D varsData2(nfreq);
    for (int freq = 0; freq < nfreq; freq++) field2D_init(varsData2[freq], varList1, FIELD_VEC);

    Varray<double> array1(nts);
    Varray<double> array2(nfreq);
    Varray<double> real(seg_l);
    Varray<double> imag(seg_l);
    Varray<double> window(seg_l);

    switch (which_window)
    {
      case 0:
        for (int k = 0; k < seg_l; ++k) window[k] = 1;
        break;
      case 1:
        for (int k = 0; k < seg_l; ++k) window[k] = 1 - std::cos(2 * M_PI * (k + 1) / (seg_l + 1));
        break;
      case 2:
        for (int k = 0; k < seg_l / 2; ++k) window[k] = window[seg_l - 1 - k] = k;
        break;
      case 3:
        for (int k = 0; k < seg_l; ++k)
        {
          double temp = ((k + 1.) - 0.5 * (seg_l + 1.)) / (0.5 * (seg_l + 1.));
          window[k] = 1 - temp * temp;
        }
        break;
      default: cdo_abort("Invalid window type %d!", which_window); break;
    }

    double wssum = 0;
    for (int k = 0; k < seg_l; ++k) wssum += window[k] * window[k];

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        for (size_t i = 0; i < var.gridsize; ++i)
        {
          for (tsID = 0; tsID < nts; ++tsID) array1[tsID] = varsData[tsID][varID][levelID].vec_d[i];

          for (int freq = 0; freq < nfreq; freq++) array2[freq] = 0;

          spectrum(nts, array1.data(), array2.data(), real.data(), imag.data(), window.data(), wssum, detrend, seg_n, seg_l);

          for (int freq = 0; freq < nfreq; freq++) varsData2[freq][varID][levelID].vec_d[i] = array2[freq];
        }
      }
    }

    for (tsID = 0; tsID < nfreq; ++tsID)
    {
      taxisDefVdatetime(taxisID2, vDateTimes[0]);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto &varData2 = varsData2[tsID][varID][levelID];
          if (!varData2.empty())
          {
            varData2.numMissVals = 0;
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, varData2);
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
