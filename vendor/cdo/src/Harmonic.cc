/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Harmonic   harmonic        Harmonic
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "varray.h"
#include "util_files.h"
#include "util_string.h"

class Harmonic : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Harmonic",
    .operators = { { "harmonic" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Harmonic> registration = RegisterEntry<Harmonic>(module);

private:
  CdiDateTime vDateTime{};

  CdoStreamID streamID1{};
  std::vector<CdoStreamID> streamIDs;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int n_out{};
  int n{};

  VarList varList1;
  VarList varList2;

public:
  void
  init() override
  {
    operator_input_arg("wave number and wave length of first harmonic in number of timesteps");

    operator_check_argc(2);

    n_out = parameter_to_int(cdo_operator_argv(0));
    n = parameter_to_int(cdo_operator_argv(1));

    if (n_out > 9) cdo_abort("Maximum number of wave numbers is 9!");

    if (n < 1 || n < 2 * n_out)
      cdo_abort("The wave length must be positive and smaller than 2 times the number of requested harmonics (=%d)!", n_out);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID2);

    auto fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    streamIDs.resize(n_out);

    for (int j = 0; j < n_out; ++j)
    {
      auto fileName = cdo_get_obase() + string_format("%1d", j + 1);
      if (fileSuffix.size() > 0) fileName += fileSuffix;

      auto streamID2 = open_write(fileName);
      cdo_def_vlist(streamID2, vlistID2);
      streamIDs[j] = streamID2;
    }
  }

  void
  run() override
  {
    Varray<double> array(varList1.gridsizeMax());

    auto numVars = varList1.numVars();

    Varray3D<double> out(n_out);
    Varray3D<double> work(2 * n_out);

    for (int j = 0; j < n_out; ++j)
    {
      out[j].resize(numVars);
      for (int varID = 0; varID < numVars; ++varID)
        out[j][varID].resize(varList1.vars[varID].nlevels * varList1.vars[varID].gridsize);
    }

    for (int j = 0; j < n_out * 2; ++j)
    {
      work[j].resize(numVars);
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        work[j][varID].resize(var.gridsize * var.nlevels, 0);
      }
    }

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      if (tsID == 0) vDateTime = taxisInqVdatetime(taxisID1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);

        if (numMissVals) cdo_abort("Missing values are not allowed!");

        auto gridsize = varList1.vars[varID].gridsize;
        auto offset = gridsize * levelID;

        for (int j = 0; j < n_out; ++j)
        {
          auto scarg = 2 * M_PI * (((j + 1) * (tsID + 1)) % n) / n;
          auto sine = std::sin(scarg);
          auto cosine = std::cos(scarg);
          for (size_t i = 0; i < gridsize; ++i)
          {
            work[j][varID][i + offset] += array[i] * sine;
            work[n_out + j][varID][i + offset] += array[i] * cosine;
          }
        }
      }

      tsID++;
    }

    auto nts = tsID;

    cdo_stream_close(streamID1);

    if (nts % n) { cdo_abort("The length of first harmonic (=%d) does not divide the number of timesteps (=%d)!", n, nts); }

    for (int j = 0; j < n_out && 2 * (j + 1) < n; ++j)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        for (int levelID = 0; levelID < var2.nlevels; ++levelID)
        {
          auto offset = var2.gridsize * levelID;
          for (size_t i = 0; i < var2.gridsize; ++i)
          {
            auto sqrwork1 = work[j][varID][i + offset] * work[j][varID][i + offset];
            auto sqrwork2 = work[n_out + j][varID][i + offset] * work[n_out + j][varID][i + offset];
            out[j][varID][i + offset] = std::sqrt(sqrwork1 + sqrwork2) * 2 / nts;
          }
        }
      }
    }

    if (2 * n_out == n)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        for (int levelID = 0; levelID < var2.nlevels; ++levelID)
        {
          auto offset = var2.gridsize * levelID;
          for (size_t i = 0; i < var2.gridsize; ++i)
            out[n_out - 1][varID][i + offset] = work[2 * n_out - 1][varID][i + offset] / nts;
        }
      }
    }

    auto nout = n_out;

    taxisDefVdatetime(taxisID2, vDateTime);
    for (int j = 0; j < nout; ++j)
    {
      auto streamID2 = streamIDs[j];
      cdo_def_timestep(streamID2, 0);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        for (int levelID = 0; levelID < var2.nlevels; ++levelID)
        {
          auto offset = var2.gridsize * levelID;
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, &out[j][varID][offset], 0);
        }
      }
    }

    for (int j = 0; j < n_out && 2 * (j + 1) < n; ++j)
    {
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        for (int levelID = 0; levelID < var2.nlevels; ++levelID)
        {
          auto offset = var2.gridsize * levelID;
          for (size_t i = 0; i < var2.gridsize; ++i)
          {
            auto work1 = work[j][varID][i + offset];
            auto work2 = work[n_out + j][varID][i + offset];
            auto tmpatan2 = std::atan2(work1, work2) * n / (j + 1) / 2 / M_PI;
            out[j][varID][i + offset] = (work1 > 0.0 || work2 > 0.0) ? tmpatan2 : var2.missval;
            if (out[j][varID][i + offset] < 0) out[j][varID][i + offset] += n / (j + 1.);
          }
        }
      }
    }

    nout = n_out;
    if (2 * n_out == n) nout -= 1;

    taxisDefVdatetime(taxisID2, vDateTime);
    for (int j = 0; j < nout; ++j)
    {
      auto streamID2 = streamIDs[j];
      cdo_def_timestep(streamID2, 1);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var2 = varList2.vars[varID];
        for (int levelID = 0; levelID < var2.nlevels; ++levelID)
        {
          auto offset = var2.gridsize * levelID;
          auto numMissVals = array_num_mv(var2.gridsize, &out[j][varID][offset], var2.missval);
          cdo_def_field(streamID2, varID, levelID);
          cdo_write_field(streamID2, &out[j][varID][offset], numMissVals);
        }
      }
    }
  }

  void
  close() override
  {
    for (auto const &streamID : streamIDs) cdo_stream_close(streamID);
  }
};
