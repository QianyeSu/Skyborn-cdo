/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
/*
  cdo samplegridicon,icon_grid_0007_R02B06_G.nc,icon_grid_0019_R02B05_G.nc,icon_grid_0005_R02B04_G.nc \
                     infileR02B06.nc outR02B04_mean.nc outR02B04_std.nc
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_omp.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include <mpim_grid.h>
#include "grid_pointsearch.h"
#include "verifygrid.h"
#include "field_functions.h"

constexpr int MAX_CHILDS = 9;

namespace
{
struct CellIndex
{
  long ncells;
  Varray<long> parent;  // parent cell index
  Varray<long> child;   // child cell index
  std::string filename;
};
}  // namespace

static void
copy_data_to_index(long ncells, Varray<double> const &data, long *cellindex)
{
  for (long i = 0; i < ncells; ++i) cellindex[i] = std::lround(data[i]);
}

static void
read_cellindex(std::string const &filename, CellIndex &cellindex)
{
  auto streamID = stream_open_read_locked(filename.c_str());
  auto vlistID = streamInqVlist(streamID);
  auto numGrids = vlistNumGrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < numGrids; ++index)
  {
    gridID = vlistGrid(vlistID, index);
    if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
  }

  if (gridID == -1) cdo_abort("No ICON grid found in %s!", filename);

  // int nid = CDI_UNDEFID;
  int pid = CDI_UNDEFID;
  // int cid = CDI_UNDEFID;
  auto nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
  {
    auto varname = cdo::inq_var_name(vlistID, varID);
    /*
    if (varname == "neighbor_cell_index")
      {
        nid = varID;
        break;
      }
    */
    if (varname == "parent_cell_index")
    {
      pid = varID;
      break;
    }
    /*
    if (varname == "child_cell_index")
      {
        cid = varID;
        break;
      }
    */
  }

  // if (nid == CDI_UNDEFID) cdo_abort("neighbor_cell_index not found in %s!", filename);
  // if (pid == CDI_UNDEFID) cdo_abort("parent_cell_index not found in %s!", filename);
  // if (cid == CDI_UNDEFID) cdo_abort("child_cell_index not found in %s!", filename);

  long ncells = gridInqSize(gridID);
  cellindex.ncells = ncells;

  // cellindex.neighbor.resize(3*ncells);
  cellindex.parent.resize(ncells);
  // if (cid != CDI_UNDEFID) cellindex.child.resize(MAX_CHILDS * ncells);
  Varray<double> data(ncells);

  for (long i = 0; i < ncells; ++i) cellindex.parent[i] = 0;

  auto numFields = streamInqTimestep(streamID, 0);
  for (int fieldID = 0; fieldID < numFields; ++fieldID)
  {
    int varID, levelID;
    size_t numMissVals;
    streamInqField(streamID, &varID, &levelID);
    if (varID == pid /* || varID == nid || varID == cid */)
    {
      streamReadField(streamID, data.data(), &numMissVals);
      // if (varID == pid)
      {
        if (Options::cdoVerbose) cdo_print("Read parent_cell_index");
        copy_data_to_index(ncells, data, cellindex.parent.data());
      }
      // else if ( varID == nid ) copy_data_to_index(ncells, data, cellindex.neighbor.data() + levelID * ncells);
      // else if ( varID == cid ) copy_data_to_index(ncells, data, cellindex.child.data() + levelID * ncells);
    }
  }

  // Fortran to C index
  for (long i = 0; i < ncells; ++i) cellindex.parent[i] -= 1;
  // for ( long i = 0; i < 3*ncells; ++i ) cellindex.neighbor[i] -= 1;

  streamClose(streamID);
}

static int
read_grid(std::string const &filename)
{
  auto streamID = stream_open_read_locked(filename.c_str());
  auto vlistID = streamInqVlist(streamID);
  auto numGrids = vlistNumGrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < numGrids; ++index)
  {
    gridID = vlistGrid(vlistID, index);
    if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
  }

  if (gridID == -1) cdo_abort("No ICON grid found in %s!", filename);

  auto gridID2 = gridDuplicate(gridID);

  streamClose(streamID);

  return gridID2;
}

/**
* Return the first index of element x fits.
*
* If no interval can be found return -1.

* @param *array ascending sorted list
* @param n      length of the sorted list
* @param search the element to find a position for
*/
static long
find_index(int search, long n, Varray<long> const &array)
{
  long first = 0;
  long last = n - 1;
  long middle = (first + last) / 2;

  while (first <= last)
  {
    if (array[middle] < search)
      first = middle + 1;
    else if (array[middle] == search)
    {
      for (long i = middle; i >= 0; i--)
      {
        if (array[i] == search)
          middle = i;
        else
          break;
      }
      return middle;
    }
    else
      last = middle - 1;

    middle = (first + last) / 2;
  }

  return -1;
}

namespace
{
struct SortInfo
{
  int p, i;
};
}  // namespace

static void
compute_child_from_parent(CellIndex &cellindex1, CellIndex &cellindex2)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  auto ncells1 = cellindex1.ncells;
  auto &parent1 = cellindex1.parent;

  std::vector<long> idx1(ncells1);
  for (long i = 0; i < ncells1; ++i) idx1[i] = i;
  for (long i = 1; i < ncells1; ++i)
  {
    if (parent1[i] < parent1[i - 1])
    {
      if (Options::cdoVerbose) cdo_print("Sort parent index of %s!", cellindex1.filename);
      std::vector<SortInfo> sinfo(ncells1);
      for (long j = 0; j < ncells1; ++j)
      {
        sinfo[j].p = parent1[j];
        sinfo[j].i = idx1[j];
      }
      std::ranges::sort(sinfo, {}, &SortInfo::p);
      for (long j = 0; j < ncells1; ++j)
      {
        parent1[j] = sinfo[j].p;
        idx1[j] = sinfo[j].i;
      }
      break;
    }
  }

  auto ncells2 = cellindex2.ncells;
  cellindex2.child.resize(MAX_CHILDS * ncells2);
  auto &child2 = cellindex2.child;
  for (long i = 0; i < ncells2; ++i)
  {
    for (long k = 0; k < MAX_CHILDS; ++k) child2[i * MAX_CHILDS + k] = -1;
    auto j = find_index(i, ncells1, parent1);
    if (j < 0) continue;
    for (long k = 0; k < MAX_CHILDS; ++k)
    {
      if ((size_t) (j + k) >= idx1.size() || i != parent1[j + k]) break;
      //  child2[i*MAX_CHILDS+k] = j+k;
      child2[i * MAX_CHILDS + k] = idx1[j + k];
    }
    // if ( i%10000 == 0 ) printf("%d %d %d %d %d %d\n", i, j, parent1[j], parent1[j+1], parent1[j+2], parent1[j+3]);
  }
}

static void
read_coordinates(std::string const &filename, long n, Varray<double> &lon, Varray<double> &lat)
{
  auto streamID = streamOpenRead(filename.c_str());
  auto vlistID = streamInqVlist(streamID);
  auto numGrids = vlistNumGrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < numGrids; ++index)
  {
    gridID = vlistGrid(vlistID, index);
    if (gridInqType(gridID) == GRID_UNSTRUCTURED && (long) gridInqSize(gridID) == n && gridInqNvertex(gridID) == 3) break;
  }

  if (gridID == -1) cdo_abort("No ICON grid with %ld cells found in %s!", n, filename);

  gridInqXvals(gridID, lon.data());
  gridInqYvals(gridID, lat.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, lon, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, lat, "grid center lat");

  streamClose(streamID);
}

static void
read_coordinates(std::string const &filename, long n, Varray<double> &lon, Varray<double> &lat, int nv, Varray<double> &lon_bnds,
                 Varray<double> &lat_bnds)
{
  auto streamID = streamOpenRead(filename.c_str());
  auto vlistID = streamInqVlist(streamID);
  auto numGrids = vlistNumGrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < numGrids; ++index)
  {
    gridID = vlistGrid(vlistID, index);
    if (gridInqType(gridID) == GRID_UNSTRUCTURED && (long) gridInqSize(gridID) == n && gridInqNvertex(gridID) == 3) break;
  }

  if (gridID == -1) cdo_abort("No ICON grid with %ld cells found in %s!", n, filename);

  gridInqXvals(gridID, lon.data());
  gridInqYvals(gridID, lat.data());

  // Convert lat/lon units if required
  cdo_grid_to_radian(gridID, CDI_XAXIS, lon, "grid center lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, lat, "grid center lat");

  if (nv == 3)
  {
    gridInqXbounds(gridID, lon_bnds.data());
    gridInqYbounds(gridID, lat_bnds.data());

    cdo_grid_to_radian(gridID, CDI_XAXIS, lon_bnds, "grid corner lon");
    cdo_grid_to_radian(gridID, CDI_YAXIS, lat_bnds, "grid corner lat");
  }

  streamClose(streamID);
}

static void
compute_child_from_bounds(CellIndex &cellindex2, Varray<double> &grid_center_lon2, Varray<double> &grid_center_lat2,
                          Varray<double> &grid_corner_lon2, Varray<double> &grid_corner_lat2,
                          Varray<double> const &grid_center_lon1, Varray<double> const &grid_center_lat1)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  long ncells2 = (long) grid_center_lon2.size();

  GridPointsearch gps;
  grid_pointsearch_create_unstruct(gps, grid_center_lon1, grid_center_lat1, true);

  int ncorner = 3;

  constexpr int MaxSearch = 128;
  std::vector<KnnData> knnDataList;
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(MaxSearch);

  cellindex2.child.resize(MAX_CHILDS * ncells2);
  auto &child2 = cellindex2.child;

#ifdef _OPENMP
#pragma omp parallel for if (ncells2 > 20000) default(shared) schedule(static)
#endif
  for (long cellNo2 = 0; cellNo2 < ncells2; ++cellNo2)
  {
    for (int k = 0; k < MAX_CHILDS; ++k) child2[cellNo2 * MAX_CHILDS + k] = -1;

    Varray<Point3D> cellCorners3D(4);
    set_cell_corners_3D(ncorner, &grid_corner_lon2[cellNo2 * ncorner], &grid_corner_lat2[cellNo2 * ncorner], cellCorners3D);
    cellCorners3D[ncorner] = cellCorners3D[0];

    auto coordinateToIgnore = find_coordinate_to_ignore(cellCorners3D);

    auto cval
        = (coordinateToIgnore == 1) ? cellCorners3D[0].X : ((coordinateToIgnore == 2) ? cellCorners3D[0].Y : cellCorners3D[0].Z);
    auto invertResult = (cval < 0.0);

    Varray<Point> cellCornersPlaneProjection(4);
    set_cell_corners_plane_projection(coordinateToIgnore, ncorner, cellCorners3D, cellCornersPlaneProjection);

    auto isClockwise = are_polygon_vertices_arranged_in_clockwise_order(cellCornersPlaneProjection, ncorner + 1);

    if (invertResult) isClockwise = !isClockwise;
    if (isClockwise) continue;

    auto ompthID = cdo_omp_get_thread_num();
    auto &knnData = knnDataList[ompthID];

    grid_search_point_unstruct(gps, PointLonLat{ grid_center_lon2[cellNo2], grid_center_lat2[cellNo2] }, knnData);

    int k = 0;
    double centerCoordinates[3];
    Point3D centerPoint3D;
    for (int i = 0; i < MaxSearch; ++i)
    {
      auto cellNo1 = knnData.m_indices[i];
      if (cellNo1 < SIZE_MAX)
      {
        gcLLtoXYZ(grid_center_lon1[cellNo1], grid_center_lat1[cellNo1], centerCoordinates);
        centerPoint3D.X = centerCoordinates[0];
        centerPoint3D.Y = centerCoordinates[1];
        centerPoint3D.Z = centerCoordinates[2];

        auto centerPoint2D = set_center_point_plane_projection(coordinateToIgnore, centerPoint3D);

        auto windingNumber = winding_numbers_algorithm(cellCornersPlaneProjection, ncorner + 1, centerPoint2D);
        if (windingNumber != 0)
        {
          if (k >= MAX_CHILDS) cdo_abort("Internal problem, limit of MAX_CHILDS reached (limit=9).");
          child2[cellNo2 * MAX_CHILDS + k++] = (long) cellNo1;
        }
      }
    }
  }
}

static void
compute_child_from_coordinates(const CellIndex &cellindex1, CellIndex &cellindex2)
{
  if (Options::cdoVerbose) cdo_print("%s", __func__);

  auto ncells1 = cellindex1.ncells;
  auto ncells2 = cellindex2.ncells;

  Varray<double> lon1(ncells1), lat1(ncells1), lon2(ncells2), lat2(ncells2);
  Varray<double> lon2_bnds(3 * ncells2), lat2_bnds(3 * ncells2);

  read_coordinates(cellindex1.filename, ncells1, lon1, lat1);
  read_coordinates(cellindex2.filename, ncells2, lon2, lat2, 3, lon2_bnds, lat2_bnds);

  compute_child_from_bounds(cellindex2, lon2, lat2, lon2_bnds, lat2_bnds, lon1, lat1);
}

static void
compute_child(CellIndex &cellindex1, CellIndex &cellindex2)
{
  bool lparent = true;
  auto ncells1 = cellindex1.ncells;
  auto const &parent1 = cellindex1.parent;

  long i = 0;
  for (; i < ncells1; ++i)
    if (parent1[i] >= 0) break;
  if (i == ncells1) lparent = false;

  if (lparent) { compute_child_from_parent(cellindex1, cellindex2); }
  else
  {
    compute_child_from_coordinates(cellindex1, cellindex2);
    // cdo_abort("Missing parent index of %s!", cellindex1.filename);
  }
}

static void
compute_sum(long i, long &n, double &sum, double &sumq, long kci, std::vector<CellIndex> &cellindex, Varray<double> const &array)
{
  // printf("compute: i, kci %d %d\n", i, kci);
  auto ncells2 = cellindex[kci].ncells;
  if (i < 0 || i > ncells2) cdo_abort("Child grid cell index %ld out of bounds %ld!", i, ncells2);

  for (int k = 0; k < MAX_CHILDS; ++k)
  {
    long index = cellindex[kci].child[i * MAX_CHILDS + k];
    if (index == -1) break;
    if (kci == 1)
    {
      sum += array[index];
      sumq += array[index] * array[index];
      n += 1;
    }
    else
      compute_sum(index, n, sum, sumq, kci - 1, cellindex, array);
  }
}

static void
samplegrid(double missval, long nci, std::vector<CellIndex> &cellindex, Varray<double> const &array1, Varray<double> &array2,
           Varray<double> &array3)
{
  static bool lstat = true;
  long kci = nci - 1;
  auto ncells2 = cellindex[kci].ncells;
  long nx = 0;
  double x = 0.0;
#ifdef _OPENMP
// #pragma omp parallel for default(shared)
#endif
  for (long i = 0; i < ncells2; ++i)
  {
    long n = 0;
    double sum = 0, sumq = 0;
    compute_sum(i, n, sum, sumq, kci, cellindex, array1);
    array2[i] = n ? sum / n : missval;  // mean
    double var1 = (n * n > n) ? (sumq * n - sum * sum) / (n * n - n) : missval;
    if (var1 < 0 && var1 > -1.e-5) var1 = 0;
    array3[i] = var_to_std(var1, missval);  // std1
    if (lstat && n)
    {
      nx++;
      x += n;
    }
  }
  if (Options::cdoVerbose && lstat)
  {
    lstat = false;
    cdo_print("Mean number of childs %g", nx ? x / nx : 0);
  }
}

class Samplegridicon : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Samplegridicon",
    .operators = { { "samplegridicon", 0, 0, "samplegrids" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 2, OnlyFirst },
  };
  inline static RegisterEntry<Samplegridicon> registration = RegisterEntry<Samplegridicon>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int vlistID2{ CDI_UNDEFID };

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};

  int vlistID1{ CDI_UNDEFID };
  int gridID2{};

  Varray<double> array1;
  Varray<double> array2;
  Varray<double> array3;
  std::vector<CellIndex> cellindex;

  int nsamplegrids{};
  size_t gridsize2{};

public:
  void
  init() override
  {
    nsamplegrids = cdo_operator_argc();
    if (nsamplegrids < 2) cdo_abort("Parameter missing!");

    cellindex.resize(nsamplegrids);

    for (int i = 0; i < nsamplegrids; ++i)
    {
      read_cellindex(cdo_operator_argv(i), cellindex[i]);
      cellindex[i].filename = cdo_operator_argv(i);
      if (Options::cdoVerbose) cdo_print("Found %ld grid cells in %s", cellindex[i].ncells, cellindex[i].filename);
    }

    for (int i = 0; i < nsamplegrids - 1; ++i) compute_child(cellindex[i], cellindex[i + 1]);

    gridID2 = read_grid(cdo_operator_argv(nsamplegrids - 1).c_str());

    streamID1 = cdo_open_read(0);
    vlistID1 = cdo_stream_inq_vlist(streamID1);

    VarList varList1(vlistID1);
    long gridsizeMax = varList1.gridsizeMax();
    if (Options::cdoVerbose) cdo_print("Source gridsize = %zu", gridsizeMax);
    if (gridsizeMax != cellindex[0].ncells)
      cdo_abort("Gridsize (%ld) of input stream and first grid (%ld) differ!", gridsizeMax, cellindex[0].ncells);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizeMax *= 2;
    array1 = Varray<double>(gridsizeMax);

    vlistID2 = vlistDuplicate(vlistID1);
    auto vlistID3 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);
    vlistDefTaxis(vlistID3, taxisID3);

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID = vlistGrid(vlistID1, index);
      auto gridtype = gridInqType(gridID);
      if (!(gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3))
        cdo_abort("Unsupported gridtype: %s with %d corners", gridNamePtr(gridtype), gridInqNvertex(gridID));

      vlistChangeGridIndex(vlistID2, index, gridID2);
      vlistChangeGridIndex(vlistID3, index, gridID2);
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    gridsize2 = gridInqSize(gridID2);
    if (Options::cdoVerbose) cdo_print("Target gridsize = %ld", gridsize2);
    if (vlistNumber(vlistID2) != CDI_REAL) gridsize2 *= 2;
    array2 = Varray<double>(gridsize2);
    array3 = Varray<double>(gridsize2);
  }

  void
  run() override
  {
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID2, tsID);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto missval = vlistInqVarMissval(vlistID1, varID);

        samplegrid(missval, nsamplegrids, cellindex, array1, array2, array3);

        numMissVals = varray_num_mv(gridsize2, array2, missval);
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, array2.data(), numMissVals);

        numMissVals = varray_num_mv(gridsize2, array3, missval);
        cdo_def_field(streamID3, varID, levelID);
        cdo_write_field(streamID3, array3.data(), numMissVals);
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
    vlistDestroy(vlistID2);
  }
};
