/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

        Timstat3        varquot2test
        Timstat3        meandiff2test
*/

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "statistic.h"
#include "arithmetic.h"

constexpr int NIN = 2;
constexpr int NWORK = 6;

static void
varquot2test(double rconst, double risk, size_t gridsize, double missval, Varray2D<double> const &work, Varray<double> &out)
{
  auto missval1 = missval;
  auto missval2 = missval;

  auto varquot2test_kernel = [&](auto i, auto is_EQ)
  {
    auto temp0 = DIVMX(MULM(work[0][i], work[0][i]), work[2][i]);
    auto temp1 = DIVMX(MULM(work[3][i], work[3][i]), work[5][i]);
    auto temp2 = SUBM(work[1][i], temp0);
    auto temp3 = SUBM(work[4][i], temp1);
    auto statistic = DIVM(temp2, ADDM(temp2, MULM(rconst, temp3)));

    auto fractil1 = missval1, fractil2 = missval1;
    if (work[2][i] > 1 && work[5][i] > 1)
      cdo::beta_distr_constants((work[2][i] - 1) / 2, (work[5][i] - 1) / 2, 1 - risk, &fractil1, &fractil2);

    double result = is_EQ(statistic, missval1) ? missval1 : (statistic <= fractil1 || statistic >= fractil2);
    return result;
  };

  if (std::isnan(missval))
    for (size_t i = 0; i < gridsize; ++i) out[i] = varquot2test_kernel(i, fp_is_equal);
  else
    for (size_t i = 0; i < gridsize; ++i) out[i] = varquot2test_kernel(i, is_equal);
}

static void
meandiff2test(double rconst, double risk, size_t gridsize, double missval, Varray2D<double> const &work, Varray<double> &out)
{
  auto mul = [](auto x, auto y) { return x * y; };
  constexpr double meanFactor[] = { 1.0, -1.0 };
  constexpr double varFactor[] = { 1.0, 1.0 };
  constexpr auto factor1 = mul(mul(meanFactor[0], meanFactor[0]), varFactor[0]);
  constexpr auto factor2 = mul(mul(meanFactor[1], meanFactor[1]), varFactor[1]);

  auto missval1 = missval;
  auto missval2 = missval;

  auto meandiff2test_kernel = [&](auto i, auto is_EQ)
  {
    double temp0 = 0.0;
    double degOfFreedom = -NIN;
    auto tmp = DIVMX(MULM(work[0][i], work[0][i]), work[2][i]);
    temp0 = ADDM(temp0, DIVMX(SUBM(work[1][i], tmp), varFactor[0]));
    degOfFreedom = ADDM(degOfFreedom, work[2][i]);
    tmp = DIVMX(MULM(work[3][i], work[3][i]), work[5][i]);
    temp0 = ADDM(temp0, DIVMX(SUBM(work[4][i], tmp), varFactor[1]));
    degOfFreedom = ADDM(degOfFreedom, work[5][i]);

    if (fp_is_not_equal(temp0, missval1) && temp0 < 0) temp0 = 0;  // This is possible because of rounding errors

    auto stddevEstimator = SQRTM(DIVMX(temp0, degOfFreedom));
    auto meanEstimator = -rconst;
    meanEstimator = ADDM(meanEstimator, MULM(meanFactor[0], DIVMX(work[0][i], work[2][i])));
    meanEstimator = ADDM(meanEstimator, MULM(meanFactor[1], DIVMX(work[3][i], work[5][i])));

    double temp1 = 0.0;
    temp1 = ADDM(temp1, DIVMX(factor1, work[2][i]));
    temp1 = ADDM(temp1, DIVMX(factor2, work[5][i]));
    auto norm = SQRTM(temp1);

    auto temp2 = DIVM(DIVM(meanEstimator, norm), stddevEstimator);
    auto fractil = (degOfFreedom < 1) ? missval1 : cdo::student_t_inv(degOfFreedom, 1 - risk / 2);

    double result = (is_EQ(temp2, missval1) || is_EQ(fractil, missval1)) ? missval1 : (std::fabs(temp2) >= fractil);
    return result;
  };

  if (std::isnan(missval))
    for (size_t i = 0; i < gridsize; ++i) out[i] = meandiff2test_kernel(i, fp_is_equal);
  else
    for (size_t i = 0; i < gridsize; ++i) out[i] = meandiff2test_kernel(i, is_equal);
}

class Timstat3 : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timstat3",
    .operators = { { "meandiff2test" }, { "varquot2test" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Timstat3> registration = RegisterEntry<Timstat3>();

  int VARQUOT2TEST, MEANDIFF2TEST;
  int vlistID[NIN], vlistID2 = -1;

  CdoStreamID streamID[NIN];
  CdoStreamID streamID3;
  int taxisID1{ CDI_UNDEFID };
  int taxisID3;

  double rconst;
  double risk;

  int operatorID;

  VarList varList1;

public:
  void
  init() override
  {
    VARQUOT2TEST = module.get_id("varquot2test");
    MEANDIFF2TEST = module.get_id("meandiff2test");

    operatorID = cdo_operator_id();

    operator_input_arg("constant and risk (e.g. 0.05)");
    operator_check_argc(2);
    rconst = parameter_to_double(cdo_operator_argv(0));
    risk = parameter_to_double(cdo_operator_argv(1));

    if (rconst <= 0) cdo_abort("Constant must be positive!");
    if (risk <= 0 || risk >= 1) cdo_abort("Risk must be greater than 0 and lower than 1!");

    for (int is = 0; is < NIN; ++is) { streamID[is] = cdo_open_read(is); }
    for (int is = 0; is < NIN; ++is) { vlistID[is] = cdo_stream_inq_vlist(streamID[is]); }

    varList1 = VarList(vlistID[0]);
    for (auto &var : varList1.vars) var.memType = MemType::Double;

    for (int is = 1; is < NIN; ++is) { varList_compare(varList1, VarList(vlistID[is])); }

    auto vlistID3 = vlistDuplicate(vlistID[0]);

    taxisID1 = vlistInqTaxis(vlistID[0]);
    taxisID3 = taxisDuplicate(taxisID1);

    vlistDefTaxis(vlistID3, taxisID3);
    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    CdiDateTime vDateTime{};
    int reachedEOF[NIN]{ 0 };

    auto numVars = varList1.numVars();
    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    Field inField, outField;

    Varray4D<double> work(numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto gridsize = var.gridsize;
      auto numLevels = var.nlevels;

      work[varID].resize(numLevels);

      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        work[varID][levelID].resize(NWORK);
        for (int iw = 0; iw < NWORK; ++iw) work[varID][levelID][iw].resize(gridsize, 0);
      }
    }

    int tsID = 0;
    while (true)
    {
      int is;
      for (is = 0; is < NIN; ++is)
      {
        if (reachedEOF[is]) continue;

        auto numFields = cdo_stream_inq_timestep(streamID[is], tsID);
        if (numFields == 0)
        {
          reachedEOF[is] = 1;
          continue;
        }

        vDateTime = taxisInqVdatetime(taxisID1);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID[is]);

          if (tsID == 0 && is == 0) fieldInfoList[fieldID].set(varID, levelID);

          auto const &var = varList1.vars[varID];
          inField.init(var);
          cdo_read_field(streamID[is], inField);

          auto const &inArray = inField.vec_d;
          auto &rwork1 = work[varID][levelID][3 * is + 0];
          auto &rwork2 = work[varID][levelID][3 * is + 1];
          auto &rwork3 = work[varID][levelID][3 * is + 2];
          auto gridsize = var.gridsize;
          for (size_t i = 0; i < gridsize; ++i)
          {
            rwork1[i] += inArray[i];
            rwork2[i] += inArray[i] * inArray[i];
            rwork3[i]++;
          }
        }
      }

      for (is = 0; is < NIN; ++is)
        if (not reachedEOF[is]) break;

      if (is == NIN) break;

      tsID++;
    }

    taxisDefVdatetime(taxisID3, vDateTime);
    cdo_def_timestep(streamID3, 0);

    for (int fieldID = 0; fieldID < maxFields; ++fieldID)
    {
      auto [varID, levelID] = fieldInfoList[fieldID].get();

      auto const &var = varList1.vars[varID];
      outField.init(var);
      auto const &rwork = work[varID][levelID];

      if (operatorID == VARQUOT2TEST) { varquot2test(rconst, risk, var.gridsize, var.missval, rwork, outField.vec_d); }
      else if (operatorID == MEANDIFF2TEST) { meandiff2test(rconst, risk, var.gridsize, var.missval, rwork, outField.vec_d); }

      field_num_mv(outField);
      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, outField);
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    for (int is = 0; is < NIN; ++is) cdo_stream_close(streamID[is]);
  }
};
