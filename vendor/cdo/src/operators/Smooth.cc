/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Smooth        smooth          Smooth grid points
      Smooth        smooth9         9 point smoothing
*/

#include <atomic>
#include <sstream>

#include <cdi.h>

#include "process_int.h"
#include "param_conversion.h"
#include "cdo_timer.h"
#include <mpim_grid.h>
#include "constants.h"  // planet radius
#include "pmlist.h"
#include "cdo_options.h"
#include "progress.h"
#include "cdo_omp.h"
#include "grid_pointsearch.h"
#include "interpol.h"

namespace
{
struct SmoothPoint
{
  double arc_radius{ 0.0 };
  double radius{ 1.0 };
  size_t maxpoints{ SIZE_MAX };
  KnnParams knnParams;

  SmoothPoint()
  {
    knnParams.weighted = WeightingMethod::linear;
    knnParams.weight0 = 0.25;
    knnParams.weightR = 0.25;
  }
};
}  // namespace

template <typename T1, typename T2>
static size_t
smooth(int gridID, double mv, Varray<T1> const &array1, Varray<T2> &array2, const SmoothPoint &spoint)
{
  T1 missval = mv;
  auto gridID0 = gridID;
  auto gridsize = gridInqSize(gridID);
  auto numNeighbors = spoint.maxpoints;
  if (numNeighbors > gridsize) numNeighbors = gridsize;

  Vmask mask(gridsize);
  for (size_t i = 0; i < gridsize; ++i) mask[i] = fp_is_not_equal(array1[i], missval);

  gridID = generate_full_point_grid(gridID);
  if (!gridHasCoordinates(gridID)) cdo_abort("Cell center coordinates missing!");

  Varray<double> xvals(gridsize), yvals(gridsize);
  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, yvals, "grid center lat");

  auto knnParams = spoint.knnParams;
  knnParams.k = numNeighbors;
  knnParams.kMin = 1;
  knnParams.searchRadius = spoint.radius;

  std::vector<KnnData> knnDataList;
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(knnParams);

  cdo::timer timer;

  GridPointsearch gps;
  gps.set_radius((spoint.arc_radius > 0.0) ? arc_to_chord_length(spoint.arc_radius) : spoint.radius);
  grid_pointsearch_create_unstruct(gps, xvals, yvals, true);

  if (Options::cdoVerbose) cdo_print("Point search created: %.2f seconds (%zu points)", timer.elapsed(), gridsize);

  cdo::Progress progress;

  timer.reset();

  size_t numWeightsMin = gridsize, numWeightsMax = 0;
  std::atomic<size_t> atomicCount{ 0 }, atomicSum{ 0 }, atomicNumMiss{ 0 };

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(shared) schedule(dynamic) reduction(min : numWeightsMin) reduction(max : numWeightsMax)
#endif
  for (size_t i = 0; i < gridsize; ++i)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && gridsize > progressMinSize) progress.update((double) atomicCount / gridsize);

    auto &knnData = knnDataList[ompthID];
    grid_search_point_smooth(gps, PointLonLat{ xvals[i], yvals[i] }, knnData);

    // Compute weights if mask is false, eliminate those points
    auto numWeights = knnData.compute_weights(mask);

    array2[i] = numWeights ? knnData.array_weights_sum(array1) : missval;
    atomicSum += numWeights;
    if (numWeights == 0) atomicNumMiss++;

    if (Options::cdoVerbose)
    {
      numWeightsMin = std::min(numWeightsMin, numWeights);
      numWeightsMax = std::max(numWeightsMax, numWeights);
    }
  }

  size_t numMissValsx = atomicNumMiss;
  size_t numPoints = atomicSum;

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds (%zu points)", timer.elapsed(), numPoints);
  if (Options::cdoVerbose) cdo_print("Min/Max points found: %zu/%zu", numWeightsMin, numWeightsMax);

  if (gridID0 != gridID) gridDestroy(gridID);

  return numMissValsx;
}

static void
smooth(Field const &field1, Field &field2, const SmoothPoint &spoint)
{
  auto func = [&](auto const &v1, auto &v2) { field2.numMissVals = smooth(field1.grid, field1.missval, v1, v2, spoint); };
  field_operation2(func, field1, field2);
}

template <typename T1, typename T2>
static size_t
smooth9(int gridID, double mv, Varray<T1> const &array1, Varray<T2> &array2)
{
  T1 missval = mv;
  auto gridsize = gridInqSize(gridID);
  auto nlon = gridInqXsize(gridID);
  auto nlat = gridInqYsize(gridID);
  auto gridIsCyclic = gridIsCircular(gridID);

  Vmask mask(gridsize);

  for (size_t i = 0; i < gridsize; ++i) mask[i] = fp_is_not_equal(missval, array1[i]);

  double avg = 0;
  double divavg = 0;
  auto smooth9_sum = [&](double sfac, size_t ij)
  {
    if (mask[ij])
    {
      avg += sfac * array1[ij];
      divavg += sfac;
    }
  };

  size_t numMissVals = 0;
  for (size_t i = 0; i < nlat; ++i)
  {
    for (size_t j = 0; j < nlon; ++j)
    {
      avg = 0;
      divavg = 0;

      if ((i == 0) || (j == 0) || (i == (nlat - 1)) || (j == (nlon - 1)))
      {
        auto ij = j + nlon * i;
        // clang-format off
        if (mask[ij])
        {
          avg += array1[ij];
          divavg += 1;
          // upper left corner
          if      (i != 0 && j != 0)                   smooth9_sum(0.3, (i - 1) * nlon + j - 1);
          else if (i != 0 && gridIsCyclic)             smooth9_sum(0.3, (i - 1) * nlon + j - 1 + nlon);
          // upper cell
          if      (i != 0)                             smooth9_sum(0.5, (i - 1) * nlon + j);
          // upper right corner
          if      (i != 0 && j != (nlon - 1))          smooth9_sum(0.3, (i - 1) * nlon + j + 1);
          else if (i != 0 && gridIsCyclic)             smooth9_sum(0.3, (i - 1) * nlon + j + 1 - nlon);
          // left cell
          if      (j != 0)                             smooth9_sum(0.5, i * nlon + j - 1);
          else if (gridIsCyclic)                       smooth9_sum(0.5, i * nlon - 1 + nlon);
          // right cell
          if      (j != (nlon - 1))                    smooth9_sum(0.5, i * nlon + j + 1);
          else if (gridIsCyclic)                       smooth9_sum(0.5, i * nlon + j + 1 - nlon);
          // lower left corner
          if      (i != (nlat - 1) && j != 0)          smooth9_sum(0.3, (i + 1) * nlon + j - 1);
          else if (i != (nlat - 1) && gridIsCyclic)    smooth9_sum(0.3, (i + 1) * nlon - 1 + nlon);
          // lower cell
          if      (i != (nlat - 1))                    smooth9_sum(0.5, (i + 1) * nlon + j);
          // lower right corner
          if      (i != (nlat - 1) && j != (nlon - 1)) smooth9_sum(0.3, (i + 1) * nlon + j + 1);
          else if (i != (nlat - 1) && gridIsCyclic)    smooth9_sum(0.3, (i + 1) * nlon + j + 1 - nlon);
        }
        // clang-format on
      }
      else if (mask[j + nlon * i])
      {
        avg += array1[j + nlon * i];
        divavg += 1;

        smooth9_sum(0.3, (i - 1) * nlon + j - 1);
        smooth9_sum(0.5, (i - 1) * nlon + j);
        smooth9_sum(0.3, (i - 1) * nlon + j + 1);
        smooth9_sum(0.5, i * nlon + j - 1);
        smooth9_sum(0.5, i * nlon + j + 1);
        smooth9_sum(0.3, (i + 1) * nlon + j - 1);
        smooth9_sum(0.5, (i + 1) * nlon + j);
        smooth9_sum(0.3, (i + 1) * nlon + j + 1);
      }

      if (std::fabs(divavg) > 0) { array2[i * nlon + j] = avg / divavg; }
      else
      {
        array2[i * nlon + j] = missval;
        numMissVals++;
      }
    }
  }

  return numMissVals;
}

static void
smooth9(Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { field2.numMissVals = smooth9(field1.grid, field1.missval, v1, v2); };
  field_operation2(func, field1, field2);
}

double
radiusDegToKm(double radiusInDeg)
{
  return radiusInDeg * (2.0 * PlanetRadiusDefault * M_PI) / (360.0 * 1000.0);
}

static void
get_parameter(int &xnsmooth, SmoothPoint &spoint)
{
  auto &knnParams = spoint.knnParams;
  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "nsmooth")     xnsmooth = parameter_to_int(value);
      else if (key == "maxpoints")   spoint.maxpoints = parameter_to_size_t(value);
      else if (key == "radius")      spoint.radius = radius_str_to_deg(value);
      else if (key == "arc_radius")  spoint.arc_radius = radius_str_to_deg(value);
      else if (key == "weighted")    knnParams.weighted = string_to_weightingMethod(parameter_to_word(value));
      else if (key == "gauss_scale") knnParams.gaussScale = parameter_to_double(value);
      else if (key == "weight0")     knnParams.weight0 = parameter_to_double(value);
      else if (key == "weightR")     knnParams.weightR = parameter_to_double(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

static void
print_parameter(SmoothPoint const &sp)
{
  auto const &kp = sp.knnParams;
  std::stringstream outbuffer;

  if (sp.arc_radius > 0.0)
    outbuffer << "arc_radius=" << sp.arc_radius << "deg(" << radiusDegToKm(sp.arc_radius) << "km)";
  else
    outbuffer << "radius=" << sp.radius << "deg(" << radiusDegToKm(sp.radius) << "km)";

  outbuffer << ", maxpoints=" << sp.maxpoints;
  outbuffer << ", weighted=" << weightingMethod_to_string(kp.weighted);
  if (kp.weighted == WeightingMethod::linear) outbuffer << ", weight0=" << kp.weight0 << ", weightR=" << kp.weightR;
  if (kp.weighted == WeightingMethod::gaussWeighted) outbuffer << ", gauss_scale=" << kp.gaussScale;

  cdo_print("%s", outbuffer.str());
}

static void
check_radius_range(double radius, const char *name)
{
  if (radius < 0.0 || radius > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", name, radius);
}

class Smooth : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Smooth",
    .operators = { { "smooth", SmoothHelp }, { "smooth9", SmoothHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Smooth> registration = RegisterEntry<Smooth>();

  int SMOOTH{}, SMOOTH9{};
  int numVars{};
  VarList varList1{};
  std::vector<bool> varIDs;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int xnsmooth = 1;
  int operatorID{};

  SmoothPoint spoint{};

public:
  void
  init() override
  {
    SMOOTH = module.get_id("smooth");
    SMOOTH9 = module.get_id("smooth9");

    operatorID = cdo_operator_id();

    if (operatorID == SMOOTH) get_parameter(xnsmooth, spoint);

    check_radius_range(spoint.radius, "radius");
    check_radius_range(spoint.arc_radius, "arc_radius");

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    varIDs.resize(numVars, false);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto gridID = var.gridID;
      auto gridtype = gridInqType(gridID);
      if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR || gridtype == GRID_PROJECTION
          || (operatorID == SMOOTH9 && gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0))
      {
        varIDs[varID] = true;
      }
      else if (operatorID == SMOOTH && gridtype == GRID_UNSTRUCTURED) { varIDs[varID] = true; }
      else { cdo_warning("Unsupported grid for variable %s", var.name); }
    }

    if (varList1.gridsizeMax() < spoint.maxpoints) spoint.maxpoints = varList1.gridsizeMax();
    if (Options::cdoVerbose && operatorID == SMOOTH) print_parameter(spoint);

    spoint.radius *= DEG2RAD;
    spoint.arc_radius *= DEG2RAD;

    streamID2 = cdo_open_write(1);
  }

  void
  run() override
  {
    Field field1, field2;

    cdo_def_vlist(streamID2, vlistID2);

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
        auto const &var = varList1.vars[varID];
        field1.init(var);
        field2.init(var);
        cdo_read_field(streamID1, field1);

        if (varIDs[varID])
        {
          for (int i = 0; i < xnsmooth; ++i)
          {
            if (operatorID == SMOOTH)
              smooth(field1, field2, spoint);
            else if (operatorID == SMOOTH9)
              smooth9(field1, field2);

            field_copy(field2, field1);
          }
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field1);
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
