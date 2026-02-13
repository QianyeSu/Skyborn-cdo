/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Ralf Müller

*/

/*
   This module contains the following operators:

*/

#include <atomic>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_timer.h"
#include <mpim_grid.h>
#include "grid_pointsearch.h"
#include "cdo_options.h"
#include "progress.h"
#include "cdo_omp.h"
#include "matrix_view.h"

template <typename T, typename CMP_FUNC>
T
fillmiss_kernel(int nfill, bool globgrid, long nx, long ny, long i, long j, T missval, MatrixView<T> &matrix1, CMP_FUNC is_NE)
{
  if (is_NE(matrix1[j][i], missval)) return matrix1[j][i];

  T rval = missval;
  long ir, iu, il, io;
  long k1, k2;
  double s1, s2;

  long kr = 0, ku = 0, kl = 0, ko = 0;
  double xr = 0.0, xu = 0.0, xl = 0.0, xo = 0.0;

  for (ir = i + 1; ir < nx; ir++)
    if (is_NE(matrix1[j][ir], missval))
    {
      kr = ir - i;
      xr = matrix1[j][ir];
      break;
    }

  if (globgrid && ir == nx)
  {
    for (ir = 0; ir < i; ir++)
      if (is_NE(matrix1[j][ir], missval))
      {
        kr = nx + ir - i;
        xr = matrix1[j][ir];
        break;
      }
  }

  for (il = i - 1; il >= 0; il--)
    if (is_NE(matrix1[j][il], missval))
    {
      kl = i - il;
      xl = matrix1[j][il];
      break;
    }

  if (globgrid && il == -1)
  {
    for (il = nx - 1; il > i; il--)
      if (is_NE(matrix1[j][il], missval))
      {
        kl = nx + i - il;
        xl = matrix1[j][il];
        break;
      }
  }

  for (iu = j + 1; iu < ny; iu++)
    if (is_NE(matrix1[iu][i], missval))
    {
      ku = iu - j;
      xu = matrix1[iu][i];
      break;
    }

  for (io = j - 1; io >= 0; io--)
    if (is_NE(matrix1[io][i], missval))
    {
      ko = j - io;
      xo = matrix1[io][i];
      break;
    }

  // printf("%d %d %d %d %d %d %g %g %g %g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo);

  auto kh = kl + kr;
  auto kv = ko + ku;
  // clang-format off
  if      (kh == 0) { k1 = 0; s1 = 0.0; }
  else if (kl == 0) { k1 = 1; s1 = xr; }
  else if (kr == 0) { k1 = 1; s1 = xl; }
  else              { k1 = 2; s1 = xr * kl / kh + xl * kr / kh; }

  if      (kv == 0) { k2 = 0; s2 = 0.0; }
  else if (ku == 0) { k2 = 1; s2 = xo; }
  else if (ko == 0) { k2 = 1; s2 = xu; }
  else              { k2 = 2; s2 = xu * ko / kv + xo * ku / kv; }

  auto kk = k1 + k2;
  if (kk >= nfill)
    {
      if      (kk == 0) cdo_abort("no point found!");
      else if (k1 == 0) rval = s2;
      else if (k2 == 0) rval = s1;
      else              rval = s1 * k2 / kk + s2 * k1 / kk;
    }
  else
    rval = matrix1[j][i];
  // clang-format on

  return rval;
}

template <typename T1, typename T2, typename CMP_FUNC>
void
fillmiss_x(Varray<T1> &vIn, Varray<T2> &vOut, int gridID, double mv, int nfill, CMP_FUNC is_NE)
{
  T1 missval = mv;
  long nx = gridInqXsize(gridID);
  long ny = gridInqYsize(gridID);
  auto globgrid = (bool) gridIsCircular(gridID);

  auto gridtype = gridInqType(gridID);
  if (!(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)) cdo_abort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  MatrixView<T1> matrix1(vIn.data(), ny, nx);
  MatrixView<T2> matrix2(vOut.data(), ny, nx);

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (long j = 0; j < ny; ++j)
    for (long i = 0; i < nx; ++i) { matrix2[j][i] = fillmiss_kernel(nfill, globgrid, nx, ny, i, j, missval, matrix1, is_NE); }
}

static void
fillmiss(Field &field1, Field &field2, int nfill)
{
  auto func = [&](auto &v1, auto &v2, auto grid, double mv, auto is_NE) { fillmiss_x(v1, v2, grid, mv, nfill, is_NE); };
  std::isnan(field1.missval) ? field_operation2(func, field1, field2, field1.grid, field1.missval, fp_is_not_equal)
                             : field_operation2(func, field1, field2, field1.grid, field1.missval, is_not_equal);
}

template <typename T, typename CMP_FUNC>
T
fillmiss_one_step_kernel(long nx, long ny, long i, long j, T missval, MatrixView<T> &matrix1, CMP_FUNC is_NE)
{
  if (is_NE(matrix1[j][i], missval)) return matrix1[j][i];

  T rval = missval;
  long ir, iu, il, io;
  long k1, k2;
  T s1, s2;

  long kr = 0, ku = 0, kl = 0, ko = 0;
  T xr = 0.0, xu = 0.0, xl = 0.0, xo = 0.0;

  for (ir = i + 1; ir < nx; ir++)
    if (is_NE(matrix1[j][ir], missval))
    {
      kr = ir - i;
      xr = matrix1[j][ir];
      break;
    }

  for (il = i - 1; il >= 0; il--)
    if (is_NE(matrix1[j][il], missval))
    {
      kl = i - il;
      xl = matrix1[j][il];
      break;
    }

  for (iu = j + 1; iu < ny; iu++)
    if (is_NE(matrix1[iu][i], missval))
    {
      ku = iu - j;
      xu = matrix1[iu][i];
      break;
    }

  for (io = j - 1; io >= 0; io--)
    if (is_NE(matrix1[io][i], missval))
    {
      ko = j - io;
      xo = matrix1[io][i];
      break;
    }

  auto kh = kl + kr;
  auto kv = ko + ku;
  // clang-format off
  if      (kh == 0) { s1 = 0.0; k1 = 0; }
  else if (kl == 0) { s1 = xr;  k1 = kr; }
  else if (kr == 0) { s1 = xl;  k1 = kl; }
  else              { s1 = (kl < kr) ? xl : xr;  k1 = (kl < kr) ? kl : kr; }

  if      (kv == 0) { s2 = 0.0; k2 = 0; }
  else if (ku == 0) { s2 = xo;  k2 = ko; }
  else if (ko == 0) { s2 = xu;  k2 = ku; }
  else              { s2 = (ku < ko) ? xu : xo;  k2 = (ku < ko) ? ku : ko; }

  auto kk = k1 + k2;
  if      (kk == 0) rval = matrix1[j][i];
  else if (k1 == 0) rval = s2;
  else if (k2 == 0) rval = s1;
  else              rval = (k1 <= k2) ? s1 : s2;
  // clang-format on

  return rval;
}

template <typename T1, typename T2, typename CMP_FUNC>
void
fillmiss_one_step_x(Varray<T1> &vIn, Varray<T2> &vOut, int gridID, double mv, int maxfill, CMP_FUNC is_EQ)
{
  T1 missval = mv;
  long nx = gridInqXsize(gridID);
  long ny = gridInqYsize(gridID);

  MatrixView<T1> matrix1(vIn.data(), ny, nx);
  MatrixView<T2> matrix2(vOut.data(), ny, nx);

  for (int fill_iterations = 0; fill_iterations < maxfill; fill_iterations++)
  {
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (long j = 0; j < ny; ++j)
      for (long i = 0; i < nx; ++i) matrix2[j][i] = fillmiss_one_step_kernel(nx, ny, i, j, missval, matrix1, is_EQ);

    if ((fill_iterations + 1) < maxfill)
      for (long j = 0; j < ny; ++j)
        for (long i = 0; i < nx; ++i) matrix1[j][i] = matrix2[j][i];
  }
}

static void
fillmiss_one_step(Field &field1, Field &field2, int maxfill)
{
  auto func = [&](auto &v1, auto &v2, auto grid, double mv, auto is_NE) { fillmiss_one_step_x(v1, v2, grid, mv, maxfill, is_NE); };
  std::isnan(field1.missval) ? field_operation2(func, field1, field2, field1.grid, field1.missval, fp_is_not_equal)
                             : field_operation2(func, field1, field2, field1.grid, field1.missval, is_not_equal);
}

template <typename T1, typename T2>
void
setmisstodis(Varray<T1> &vIn, Varray<T2> &vOut, int gridID, size_t numMissVals, double mv, int numNeighbors)
{
  T1 missval = mv;
  auto gridID0 = gridID;

  auto gridsize = gridInqSize(gridID);
  auto nvals = gridsize - numMissVals;
  gridID = generate_full_point_grid(gridID);

  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xvals(gridsize);
  Varray<double> yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  static bool doCheck = true;
  if (doCheck)
  {
    doCheck = false;
    check_longitude_range(xvals, "center", cdo_grid_get_units(gridID, CDI_XAXIS, "grid center lon"));
    check_latitude_range(yvals, "center", cdo_grid_get_units(gridID, CDI_YAXIS, "grid center lat"));
  }

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yvals, "grid center lat");

  std::vector<size_t> mindex(numMissVals, 1), vindex(nvals, 1);
  Varray<double> lons(nvals), lats(nvals);

  size_t nv = 0, nm = 0;
  for (size_t i = 0; i < gridsize; ++i)
  {
    vOut[i] = vIn[i];
    if (fp_is_equal(vIn[i], missval))
    {
      mindex[nm] = i;
      nm++;
    }
    else
    {
      if (nv < nvals)
      {
        lons[nv] = xvals[i];
        lats[nv] = yvals[i];
        vindex[nv] = i;
      }
      nv++;
    }
  }

  if (nv != nvals) cdo_abort("Internal problem, number of valid values differ!");

  std::vector<KnnData> knnDataList;
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(numNeighbors);

  cdo::timer timer;

  GridPointsearch gps;

  if (numMissVals)
  {
    gps.enable_extrapolation();
    grid_pointsearch_create_unstruct(gps, lons, lats, true);
  }

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds", timer.elapsed());

  cdo::Progress progress;

  timer.reset();

  std::atomic<size_t> atomicCount{ 0 };

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t i = 0; i < numMissVals; ++i)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && numMissVals > progressMinSize) progress.update((double) atomicCount / numMissVals);

    auto &knnData = knnDataList[ompthID];

    grid_search_point_unstruct(gps, PointLonLat{ xvals[mindex[i]], yvals[mindex[i]] }, knnData);

    // Compute weights if mask is false, eliminate those points
    auto numWeights = knnData.compute_weights();
    if (numWeights)
    {
      double result = 0.0;
      for (size_t n = 0; n < numWeights; ++n) result += vIn[vindex[knnData.m_indices[n]]] * knnData.m_dist[n];
      vOut[mindex[i]] = result;
    }
  }

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", timer.elapsed());

  if (gridID0 != gridID) gridDestroy(gridID);
}

static void
setmisstodis(Field &field1, Field &field2, int numNeighbors)
{
  auto func = [&](auto &v1, auto &v2, auto grid, auto numMissVals, double mv)
  { setmisstodis(v1, v2, grid, numMissVals, mv, numNeighbors); };
  field_operation2(func, field1, field2, field1.grid, field1.numMissVals, field1.missval);
}

class Fillmiss : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Fillmiss",
    .operators = { { "fillmiss", 0, 0, "nfill" },
                   { "fillmiss2", 0, 0, "nfill" },
                   { "setmisstonn", 0, 0, "", SetmissHelp },
                   { "setmisstodis", 0, 0, "numberofneighbors", SetmissHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Fillmiss> registration = RegisterEntry<Fillmiss>();

  int FILLMISS{}, FILLMISS2{}, SETMISSTONN{}, SETMISSTODIS{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int operatorID{};
  int nfill{};

  VarList varList1;

  void (*fill_method)(Field &, Field &, int) = &setmisstodis;

public:
  void
  init() override
  {
    FILLMISS = module.get_id("fillmiss");
    FILLMISS2 = module.get_id("fillmiss2");
    SETMISSTONN = module.get_id("setmisstonn");
    SETMISSTODIS = module.get_id("setmisstodis");

    operatorID = cdo_operator_id();

    // clang-format off
    if      (operatorID == FILLMISS)      fill_method = &fillmiss;
    else if (operatorID == FILLMISS2)     fill_method = &fillmiss_one_step;
    else if (operatorID == SETMISSTONN)   fill_method = &setmisstodis;
    else if (operatorID == SETMISSTODIS)  fill_method = &setmisstodis;
    // clang-format on

    nfill = (operatorID == SETMISSTODIS) ? 4 : 1;

    // Argument handling
    auto oargc = cdo_operator_argc();
    if (oargc == 1)
    {
      nfill = parameter_to_int(cdo_operator_argv(0));
      if (operatorID == FILLMISS && (nfill < 1 || nfill > 4)) cdo_abort("nfill out of range!");
    }
    else if (oargc > 1)
      cdo_abort("Too many arguments!");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Field field1, field2;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);

      cdo_def_timestep(streamID2, tsID);

      while (numFields--)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        field1.init(var1);
        cdo_read_field(streamID1, field1);

        cdo_def_field(streamID2, varID, levelID);

        if (field1.numMissVals == 0 || field1.numMissVals == var1.gridsize) { cdo_write_field(streamID2, field1); }
        else
        {
          auto gridtype = var1.gridType;
          if ((operatorID == FILLMISS || operatorID == FILLMISS2) && (gridtype == GRID_GME || gridtype == GRID_UNSTRUCTURED))
            cdo_abort("%s data unsupported!", gridNamePtr(gridtype));

          field2.init(var1);

          fill_method(field1, field2, nfill);

          field2.numMissVals = field_num_mv(field2);

          cdo_write_field(streamID2, field2);
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
