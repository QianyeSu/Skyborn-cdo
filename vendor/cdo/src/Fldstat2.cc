/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Fldstat2    fldcor         Correlation in grid space
      Fldstat2    fldcovar       Covariance in grid space
*/

#include <cdi.h>

#include "arithmetic.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "field_functions.h"

// routine corr copied from PINGO
// correclation in space
auto correlation_kernel = [](auto v1, auto mv1, auto v2, auto mv2, auto w, auto &sum0, auto &sum1, auto &sum00, auto &sum01,
                             auto &sum11, auto &wsum0, auto is_NE)
{
  if (is_NE(w, mv1) && is_NE(v1, mv1) && is_NE(v2, mv2))
  {
    sum0 += w * v1;
    sum1 += w * v2;
    sum00 += w * v1 * v1;
    sum01 += w * v1 * v2;
    sum11 += w * v2 * v2;
    wsum0 += w;
  }
};

template <typename T1, typename T2>
static double
correlation(Varray<T1> const &v1, Varray<T2> const &v2, double missval1, double missval2, size_t gridsize,
            Varray<double> const &weight)
{
  double sum0 = 0.0, sum1 = 0.0, sum00 = 0.0, sum01 = 0.0, sum11 = 0.0, wsum0 = 0.0;

  if (std::isnan(missval1) || std::isnan(missval2))
  {
    for (size_t i = 0; i < gridsize; ++i)
      correlation_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum00, sum01, sum11, wsum0, fp_is_not_equal);
  }
  else
  {
    for (size_t i = 0; i < gridsize; ++i)
      correlation_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum00, sum01, sum11, wsum0, is_not_equal);
  }

  auto is_EQ = fp_is_equal;
  auto out = is_not_equal(wsum0, 0.0)
                 ? DIVM((sum01 * wsum0 - sum0 * sum1), SQRTM((sum00 * wsum0 - sum0 * sum0) * (sum11 * wsum0 - sum1 * sum1)))
                 : missval1;

  return out;
}

static double
correlation(Field const &field1, Field const &field2, Varray<double> const &weight)
{
  auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2, size_t size)
  { return correlation(v1, v2, mv1, mv2, size, weight); };
  return field_operation2(func, field1, field2, field1.missval, field2.missval, field1.size);
}

// covariance in space
auto covariance_kernel
    = [](auto v1, auto mv1, auto v2, auto mv2, auto w, auto &sum0, auto &sum1, auto &sum01, auto &wsum0, auto is_NE)
{
  if (is_NE(w, mv1) && is_NE(v1, mv1) && is_NE(v2, mv2))
  {
    sum0 += w * v1;
    sum1 += w * v2;
    sum01 += w * v1 * v2;
    wsum0 += w;
  }
};

template <typename T1, typename T2>
static double
covariance(Varray<T1> const &v1, Varray<T2> const &v2, double missval1, double missval2, size_t gridsize,
           Varray<double> const &weight)
{
  double sum0 = 0.0, sum1 = 0.0, sum01 = 0.0, wsum0 = 0.0;

  if (std::isnan(missval1) || std::isnan(missval2))
  {
    for (size_t i = 0; i < gridsize; ++i)
      covariance_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum01, wsum0, fp_is_not_equal);
  }
  else
  {
    for (size_t i = 0; i < gridsize; ++i)
      covariance_kernel(v1[i], missval1, v2[i], missval2, weight[i], sum0, sum1, sum01, wsum0, is_not_equal);
  }

  auto out = is_not_equal(wsum0, 0.0) ? (sum01 * wsum0 - sum0 * sum1) / (wsum0 * wsum0) : missval1;

  return out;
}

static double
covariance(Field const &field1, Field const &field2, Varray<double> const &weight)
{
  auto func = [&](auto const &v1, auto const &v2, double mv1, double mv2, size_t size)
  { return covariance(v1, v2, mv1, mv2, size, weight); };
  return field_operation2(func, field1, field2, field1.missval, field2.missval, field1.size);
}

class Fldstat2 : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Fldstat2",
    .operators = { { "fldcor", FieldFunc_Cor, 0, FldcorHelp }, { "fldcovar", FieldFunc_Covar, 0, FldcovarHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Fldstat2> registration = RegisterEntry<Fldstat2>(module);
  int operfunc{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID3{};

  bool wstatus = false;
  bool needWeights = true;

  VarList varList1{};
  VarList varList2{};

  Varray<double> weight{};

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = cdo_operator_f1(operatorID);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    double slon = 0.0, slat = 0.0;
    auto gridID3 = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(gridID3, 1);
    gridDefYsize(gridID3, 1);
    gridDefXvals(gridID3, &slon);
    gridDefYvals(gridID3, &slat);

    auto numGrids = vlistNumGrids(vlistID1);

    for (int index = 0; index < numGrids; ++index) vlistChangeGridIndex(vlistID3, index, gridID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    if (needWeights) weight.resize(varList1.gridsizeMax());
  }

  void
  run() override
  {
    Field field1, field2;

    int lastgridID = -1;
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields2 == 0)
      {
        cdo_warning("Input streams have different number of time steps!");
        break;
      }

      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        auto const &var2 = varList1.vars[varID];
        field1.init(var1);
        (void) cdo_inq_field(streamID2);
        field2.init(var2);
        cdo_read_field(streamID1, field1);
        cdo_read_field(streamID2, field2);

        auto gridID = var1.gridID;
        if (needWeights && gridID != lastgridID)
        {
          lastgridID = gridID;
          wstatus = (gridcell_weights(gridID, weight) != 0);
        }
        if (wstatus && tsID == 0 && levelID == 0) cdo_warning("Using constant grid cell area weights for variable %s!", var1.name);

        double sglval = 0.0;
        if (operfunc == FieldFunc_Cor)
          sglval = correlation(field1, field2, weight);
        else if (operfunc == FieldFunc_Covar)
          sglval = covariance(field1, field2, weight);

        auto numMissVals3 = fp_is_equal(sglval, var1.missval) ? 1 : 0;

        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, &sglval, numMissVals3);
      }

      tsID++;
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
