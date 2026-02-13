/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Spectral   sp2gp           Spectral to gridpoint
      Spectral   sp2gpl          Spectral to gridpoint linear (sp2gp,linear)
      Spectral   gp2sp           Gridpoint to spectral
      Spectral   gp2spl          Gridpoint to spectral linear (gp2sp,linear)
      Spectral   sp2sp           Spectral to spectral
      Spectral   spcut           Cut spectral wave number
*/

#include <cdi.h>

#include "cdo_vlist.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "griddes.h"
#include "specspace.h"

static int
gp2sp_init(SP_Transformation &spTrans, int gridID1, int gridIDsp, int defaultTrunc, int (*nlat2ntr)(int))
{
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);

  long ntr = nlat2ntr(nlat);
  if (defaultTrunc > 0)
  {
    if (defaultTrunc > ntr) cdo_abort("Output trunctation=%d muss be lower than input trunctation=%d", defaultTrunc, ntr);
    ntr = defaultTrunc;
  }
  if (Options::cdoVerbose) cdo_print("trunc=%ld\n", ntr);

  if (gridIDsp != -1)
    if (ntr != gridInqTrunc(gridIDsp)) gridIDsp = -1;

  if (gridIDsp == -1)
  {
    gridIDsp = gridCreate(GRID_SPECTRAL, (ntr + 1) * (ntr + 2));
    gridDefTrunc(gridIDsp, ntr);
    gridDefComplexPacking(gridIDsp, 1);
  }

  if (gridIDsp == -1) cdo_abort("Computation of spherical harmonics failed!");

  int gridID2 = gridIDsp;

  ntr = gridInqTrunc(gridID2);
  spTrans.init(nlon, nlat, ntr, PolFlag::FC2SP);

  return gridID2;
}

static int
sp2gp_init(SP_Transformation &spTrans, int gridID1, int gridIDsp, int gridIDgp, int defaultTrunc, int (*nlat2ntr)(int),
           const char *ctype)
{
  if (defaultTrunc > 0) gridIDgp = -1;

  if (gridIDgp != -1)
  {
    long nlat = gridInqYsize(gridIDgp);
    long ntr = nlat2ntr(nlat);
    if (gridInqTrunc(gridIDsp) != ntr) gridIDgp = -1;
  }

  if (gridIDgp == -1)
  {
    int ntr = gridInqTrunc(gridIDsp);
    if (defaultTrunc > 0)
    {
      if (defaultTrunc < ntr) cdo_abort("Output trunctation=%d muss be greater than input trunctation=%d", defaultTrunc, ntr);
      ntr = defaultTrunc;
    }
    char gridname[20];
    std::snprintf(gridname, sizeof(gridname), "t%s%dgrid", ctype, ntr);
    gridIDgp = grid_from_name(gridname);
  }

  int gridID2 = gridIDgp;

  long ntr = gridInqTrunc(gridID1);
  long nlon = gridInqXsize(gridID2);
  long nlat = gridInqYsize(gridID2);
  spTrans.init(nlon, nlat, ntr, PolFlag::SP2FC);

  return gridID2;
}

class Spectral : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Spectral",
    .operators = { { "gp2sp", SpectralHelp },
                   { "gp2spl", SpectralHelp },
                   { "sp2gp", SpectralHelp },
                   { "sp2gpl", SpectralHelp },
                   { "sp2sp", SpecconvHelp },
                   { "spcut", SpecconvHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Spectral> registration = RegisterEntry<Spectral>();

private:
  int GP2SP{}, GP2SPL{}, SP2GP{}, SP2GPL{}, SP2SP{}, SPCUT{};

  int gridID1 = -1, gridID2 = -1;
  int defaultTrunc = 0;
  Varray<int> waves{};
  SP_Transformation spTrans{};
  bool dataIsUnchanged{};

  int operatorID{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  VarList varList1{};

  bool lgp2sp{};
  bool lsp2gp{};
  bool linear{};

public:
  void
  init() override
  {
    dataIsUnchanged = data_is_unchanged();

    GP2SP = module.get_id("gp2sp");
    GP2SPL = module.get_id("gp2spl");
    SP2GP = module.get_id("sp2gp");
    SP2GPL = module.get_id("sp2gpl");
    SP2SP = module.get_id("sp2sp");
    SPCUT = module.get_id("spcut");

    operatorID = cdo_operator_id();

    lgp2sp = (operatorID == GP2SP || operatorID == GP2SPL);
    lsp2gp = (operatorID == SP2GP || operatorID == SP2GPL);
    linear = (operatorID == GP2SPL || operatorID == SP2GPL);

    int (*nlat2ntr)(int) = linear ? nlat_to_ntr_linear : nlat_to_ntr;
    const char *ctype = linear ? "l" : "";

    auto paramArgc = cdo_operator_argc();
    if ((lgp2sp || lsp2gp) && paramArgc == 1)
    {
      auto const &parg = cdo_operator_argv(0);
      auto pos = parg.find('=');
      if (pos > 0 && parg.substr(0, pos) == "trunc") { defaultTrunc = parameter_to_int(parg.substr(pos + 1)); }
      else
      {
        auto type = parameter_to_word((pos > 0 && parg.substr(0, pos) == "type") ? parg.substr(pos + 1) : parg);
        if (type == "linear")
        {
          nlat2ntr = nlat_to_ntr_linear;
          ctype = "l";
        }
        else if (type == "cubic")
        {
          nlat2ntr = nlat_to_ntr_cubic;
          ctype = "c";
        }
        else if (type == "quadratic") { nlat2ntr = nlat_to_ntr; }
        else
          cdo_abort("Unsupported type: %s\n", type);
      }
    }
    else if (paramArgc > 0 && operatorID != SP2SP && operatorID != SPCUT) { cdo_abort("Too many parameters"); }

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto gridIDsp = vlist_get_first_spectral_grid(vlistID1);
    auto gridIDgp = vlist_get_first_gaussian_grid(vlistID1);

    // define output grid
    if (lgp2sp)
    {
      if (gridIDgp == -1) cdo_warning("No data on regular Gaussian grid found!");

      gridID1 = gridIDgp;
      if (gridID1 != -1) gridID2 = gp2sp_init(spTrans, gridID1, gridIDsp, defaultTrunc, nlat2ntr);
    }
    else if (lsp2gp)
    {
      if (gridIDsp == -1) cdo_warning("No spectral data found!");

      gridID1 = gridIDsp;
      if (gridID1 != -1) gridID2 = sp2gp_init(spTrans, gridID1, gridIDsp, gridIDgp, defaultTrunc, nlat2ntr, ctype);
    }
    else if (operatorID == SP2SP)
    {
      gridID1 = gridIDsp;

      operator_input_arg("truncation");
      if (gridID1 != -1)
      {
        if (!std::isdigit(cdo_operator_argv(0)[0])) cdo_abort("parameter truncation must comprise only digits [0-9]!");
        long ntr = parameter_to_int(cdo_operator_argv(0));
        long nsp = (ntr + 1) * (ntr + 2);
        gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
        gridDefTrunc(gridIDsp, ntr);
        gridDefComplexPacking(gridIDsp, 1);
      }
      else
        cdo_abort("No spectral data found!");

      gridID2 = gridIDsp;
    }
    else if (operatorID == SPCUT)
    {
      gridID1 = gridIDsp;

      operator_input_arg("wave numbers");
      if (gridID1 != -1)
      {
        long maxWaveNumbers = 1 + gridInqTrunc(gridID1);
        auto waveNumber = cdo_argv_to_intarr(cdo_get_oper_argv());
        long ncut = waveNumber.size();
        waves.resize(maxWaveNumbers);
        for (long i = 0; i < maxWaveNumbers; ++i) waves[i] = 1;
        for (long i = 0; i < ncut; ++i)
        {
          long j = waveNumber[i] - 1;
          if (j < 0 || j >= maxWaveNumbers)
            cdo_abort("wave number %ld out of range (min=1, max=%l qd)!", waveNumber[i], maxWaveNumbers);
          waves[j] = 0;
        }
      }
      else
        cdo_abort("No spectral data found!");

      gridID2 = gridIDsp;
    }

    if (gridID1 != -1) vlistChangeGrid(vlistID2, gridID1, gridID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    std::vector<bool> processVars(varList1.numVars());
    for (auto const &var : varList1.vars) { processVars[var.ID] = (gridID1 == var.gridID); }

    Varray<double> array1(varList1.gridsizeMax());
    Varray<double> array2;
    if (gridID2 != -1) array2.resize(gridInqSize(gridID2));

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        if (processVars[varID])
        {
          size_t numMissVals;
          cdo_read_field(streamID1, array1.data(), &numMissVals);
          if (numMissVals) cdo_abort("Missing values unsupported for spectral data!");

          gridID1 = varList1.vars[varID].gridID;
          // clang-format off
          if      (lgp2sp) grid2spec(spTrans, gridID1, array1, gridID2, array2);
          else if (lsp2gp) spec2grid(spTrans, gridID1, array1, gridID2, array2);
          else if (operatorID == SP2SP) spec2spec(gridID1, array1, gridID2, array2);
          else if (operatorID == SPCUT) speccut(gridID1, array1, array2, waves);
          // clang-format on

          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, array2.data(), numMissVals);
        }
        else
        {
          cdo_def_field(streamID2, varID, levelID);
          if (dataIsUnchanged) { cdo_copy_field(streamID1, streamID2); }
          else
          {
            size_t numMissVals;
            cdo_read_field(streamID1, array1.data(), &numMissVals);
            cdo_write_field(streamID2, array1.data(), numMissVals);
          }
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
