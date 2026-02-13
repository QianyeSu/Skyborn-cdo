/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

        Timstat2        timcor      correlates two data files on the same grid
*/

#include <cdi.h>

#include "arithmetic.h"
#include "field.h"
#include "process_int.h"
#include "cdo_omp.h"
#include "field_functions.h"

double
calc_pvalue(double cor, size_t n)
{
  // Author: Estanislao Gavilan
  double t_stat = cor * std::sqrt((n - 2) / (1 - cor * cor));
  double pvalue = 0.5 * (1.0 + std::erf(std::fabs(t_stat / std::sqrt(2.0))));
  return pvalue;
}

// correlation in time
template <typename T1, typename T2>
void
correlation_init(Varray<T1> const &x, Varray<T2> const &y, double missval1, double missval2, bool hasMissValues, size_t gridsize,
                 Varray<size_t> &nofvals, Varray2D<double> &work)
{
  T1 xmv = missval1;
  T2 ymv = missval2;

  auto correlation_sum = [&](auto i)
  {
    double xx = x[i];
    double yy = y[i];
    work[0][i] += xx;
    work[1][i] += yy;
    work[2][i] += xx * xx;
    work[3][i] += yy * yy;
    work[4][i] += xx * yy;
    nofvals[i]++;
  };

  if (hasMissValues)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
#endif
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_not_equal(x[i], xmv) && fp_is_not_equal(y[i], ymv)) correlation_sum(i);
    }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#if _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
#endif
    for (size_t i = 0; i < gridsize; ++i) { correlation_sum(i); }
  }
}

static void
correlation_init(size_t gridsize, Field const &field1, Field const &field2, Varray<size_t> &nofvals, Varray2D<double> &work)
{
  auto hasMissValues = (field1.numMissVals > 0 || field2.numMissVals > 0);

  auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2)
  { correlation_init(v1, v2, mv1, mv2, hasMissValues, gridsize, nofvals, work); };
  field_operation2(func, field1, field2, field1.missval, field2.missval);
}

static size_t
correlation(size_t gridsize, double missval, Varray<size_t> const &nofvals, Varray2D<double> &work)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;

  size_t numMissVals = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    double cor;
    double pvalue;
    auto nvals = nofvals[i];
    if (nvals > 0)
    {
      double dnvals = nvals;
      auto temp0 = MULM(work[0][i], work[1][i]);
      auto temp1 = SUBM(work[4][i], DIVMX(temp0, dnvals));
      auto temp2 = MULM(work[0][i], work[0][i]);
      auto temp3 = MULM(work[1][i], work[1][i]);
      auto temp4 = SUBM(work[2][i], DIVMX(temp2, dnvals));
      auto temp5 = SUBM(work[3][i], DIVMX(temp3, dnvals));
      auto temp6 = MULM(temp4, temp5);

      cor = DIVM(temp1, SQRTM(temp6));
      cor = std::clamp(cor, -1.0, 1.0);

      if (fp_is_equal(cor, missval)) numMissVals++;

      pvalue = (nvals <= 2) ? missval : ((std::fabs(cor) < 1) ? calc_pvalue(cor, nvals) : 1);
    }
    else
    {
      numMissVals++;
      cor = missval;
      pvalue = missval;
    }

    work[0][i] = cor;
    work[1][i] = pvalue;
  }

  return numMissVals;
}

// covariance in time
template <typename T1, typename T2>
static void
covariance_init(Varray<T1> const &x, Varray<T2> const &y, double missval1, double missval2, bool hasMissValues, size_t gridsize,
                Varray<size_t> &nofvals, Varray2D<double> &work)
{
  T1 xmv = missval1;
  T2 ymv = missval2;

  auto covariance_sum = [&](auto i)
  {
    double xx = x[i];
    double yy = y[i];
    work[0][i] += xx;
    work[1][i] += yy;
    work[2][i] += xx * yy;
    nofvals[i]++;
  };

  if (hasMissValues)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
#endif
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_not_equal(x[i], xmv) && fp_is_not_equal(y[i], ymv)) covariance_sum(i);
    }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#if _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
#endif
    for (size_t i = 0; i < gridsize; ++i) { covariance_sum(i); }
  }
}

static void
covariance_init(size_t gridsize, Field const &field1, Field const &field2, Varray<size_t> &nofvals, Varray2D<double> &work)
{
  auto hasMissValues = (field1.numMissVals > 0 || field2.numMissVals > 0);

  auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2)
  { covariance_init(v1, v2, mv1, mv2, hasMissValues, gridsize, nofvals, work); };
  field_operation2(func, field1, field2, field1.missval, field2.missval);
}

static size_t
covariance(size_t gridsize, double missval, Varray<size_t> const &nofvals, Varray2D<double> &work)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;

  size_t numMissVals = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    double covar;
    auto nvals = nofvals[i];
    if (nvals > 0)
    {
      double dnvals = nvals;
      auto temp = DIVMX(MULM(work[0][i], work[1][i]), dnvals * dnvals);
      covar = SUBM(DIVMX(work[2][i], dnvals), temp);
      if (fp_is_equal(covar, missval)) numMissVals++;
    }
    else
    {
      numMissVals++;
      covar = missval;
    }

    work[0][i] = covar;
  }

  return numMissVals;
}

// rms in time
template <typename T1, typename T2>
static void
rmsd_init(Varray<T1> const &x, Varray<T2> const &y, double missval1, double missval2, size_t gridsize, Varray<size_t> &nofvals,
          Varray<double> &rmsd)
{
  T1 xmv = missval1;
  T2 ymv = missval2;

  for (size_t i = 0; i < gridsize; ++i)
  {
    if (fp_is_not_equal(x[i], xmv) && fp_is_not_equal(y[i], ymv))
    {
      double xx = x[i];
      double yy = y[i];
      rmsd[i] += ((xx - yy) * (xx - yy));
      nofvals[i]++;
    }
  }
}

static void
rmsd_init(size_t gridsize, Field const &field1, Field const &field2, Varray<size_t> &nofvals, Varray<double> &rmsd)
{
  auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2) { rmsd_init(v1, v2, mv1, mv2, gridsize, nofvals, rmsd); };
  field_operation2(func, field1, field2, field1.missval, field2.missval);
}

static size_t
rmsd_compute(size_t gridsize, double missval, Varray<size_t> const &nofvals, Varray<double> &rmsd)
{
  size_t numMissVals = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    if (nofvals[i] > 0) { rmsd[i] = std::sqrt(rmsd[i] / (double) nofvals[i]); }
    else
    {
      numMissVals++;
      rmsd[i] = missval;
    }
  }

  return numMissVals;
}

class Timstat2 : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Timstat2",
    // clang-format off
    .operators = { { "timcor", FieldFunc_Cor, 5, TimcorHelp },
                   { "timcovar", FieldFunc_Covar, 3, TimcovarHelp },
                   { "timrmsd", FieldFunc_Rmsd, 1, nullptr } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Timstat2> registration = RegisterEntry<Timstat2>();

  int numWork{};

  CdiDateTime vDateTime{};

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;
  int taxisID1{ CDI_UNDEFID };
  int taxisID3;

  int operfunc{};

  VarList varList1;
  VarList varList2;

  bool doWritePvalue{ false };

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);
    numWork = cdo_operator_f2(operatorID);
    auto timeIsConst = (operfunc == FieldFunc_Rmsd);

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);

    auto numVars = varList1.numVars();

    taxisID1 = vlistInqTaxis(vlistID1);
    // auto taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);

    if (timeIsConst)
      for (int varID = 0; varID < numVars; ++varID) vlistDefVarTimetype(vlistID3, varID, TIME_CONSTANT);

    auto const &var = varList1.vars[0];
    doWritePvalue = (operfunc == FieldFunc_Cor && numVars == 1 && var.nlevels == 1);

    if (doWritePvalue)
    {
      auto varID = vlistDefVar(vlistID3, var.gridID, var.zaxisID, var.timeType);
      vlistDefVarName(vlistID3, varID, "pvalue");
    }

    vlistDefNtsteps(vlistID3, 1);

    vlistDefTaxis(vlistID3, taxisID3);
    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);
  }

  void
  run() override
  {
    Field field1, field2;
    std::vector<FieldInfo> fieldInfoList(varList1.maxFields());

    auto numVars = varList1.numVars();
    Varray4D<double> work(numVars);
    Varray3D<size_t> nofvals(numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto gridsize = var.gridsize;
      auto nlevels = var.nlevels;

      work[varID].resize(nlevels);
      nofvals[varID].resize(nlevels);

      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        nofvals[varID][levelID].resize(gridsize, 0);
        work[varID][levelID].resize(numWork);
        for (int iw = 0; iw < numWork; ++iw) work[varID][levelID][iw].resize(gridsize, 0.0);
      }
    }

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      vDateTime = taxisInqVdatetime(taxisID1);

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields != numFields2) cdo_warning("Input streams have different number of fields!");

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        (void) cdo_inq_field(streamID1);
        auto [varID, levelID] = cdo_inq_field(streamID2);
        if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);

        field1.init(varList1.vars[varID]);
        field2.init(varList2.vars[varID]);

        auto gridsize = varList1.vars[varID].gridsize;

        cdo_read_field(streamID1, field1);
        cdo_read_field(streamID2, field2);

        auto &rwork = work[varID][levelID];
        auto &rnofvals = nofvals[varID][levelID];

        if (operfunc == FieldFunc_Cor) { correlation_init(gridsize, field1, field2, rnofvals, rwork); }
        else if (operfunc == FieldFunc_Covar) { covariance_init(gridsize, field1, field2, rnofvals, rwork); }
        else if (operfunc == FieldFunc_Rmsd) { rmsd_init(gridsize, field1, field2, rnofvals, rwork[0]); }
      }

      tsID++;
    }

    tsID = 0;
    taxisDefVdatetime(taxisID3, vDateTime);
    cdo_def_timestep(streamID3, tsID);

    for (int fieldID = 0; fieldID < varList1.maxFields(); ++fieldID)
    {
      auto [varID, levelID] = fieldInfoList[fieldID].get();

      auto gridsize = varList1.vars[varID].gridsize;
      auto missval = varList1.vars[varID].missval;

      auto &rwork = work[varID][levelID];
      auto const &rnofvals = nofvals[varID][levelID];

      size_t numMissVals = 0;
      if (operfunc == FieldFunc_Cor) { numMissVals = correlation(gridsize, missval, rnofvals, rwork); }
      else if (operfunc == FieldFunc_Covar) { numMissVals = covariance(gridsize, missval, rnofvals, rwork); }
      else if (operfunc == FieldFunc_Rmsd) { numMissVals = rmsd_compute(gridsize, missval, rnofvals, rwork[0]); }

      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, rwork[0].data(), numMissVals);

      if (doWritePvalue)
      {
        cdo_def_field(streamID3, 1, levelID);
        cdo_write_field(streamID3, rwork[1].data(), numMissVals);
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
