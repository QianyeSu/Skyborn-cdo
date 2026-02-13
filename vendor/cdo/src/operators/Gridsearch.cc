/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "griddes.h"

namespace
{
struct grid_type
{
  int gridID{};
  long size{};
  long numCorners{};
  Varray<double> cell_corner_lon{};
  Varray<double> cell_corner_lat{};
};

struct cellsearch_type
{
  grid_type *srcGrid{};
  grid_type *tgtGrid{};
  Varray<float> src_cell_bound_box{};
};
}  // namespace

static grid_type *
grid_new(int gridID, const char *txt)
{
  bool lgrid_destroy = false;
  auto gridtype = gridInqType(gridID);

  if (gridtype == GRID_GME)
  {
    lgrid_destroy = true;
    auto gridID_gme = gridToUnstructured(gridID, NeedCorners::Yes);
    gridCompress(gridID_gme);
    gridID = gridID_gme;
  }

  if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
  {
    lgrid_destroy = true;
    gridID = gridToCurvilinear(gridID, NeedCorners::Yes);
  }

  if (!gridHasCoordinates(gridID)) cdo_abort("%s grid corner missing!", txt);

  grid_type *grid = new grid_type;

  grid->gridID = gridID;
  grid->size = gridInqSize(grid->gridID);
  grid->numCorners = (gridInqType(grid->gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(grid->gridID) : 4;

  // printf("%s grid size %ld nv %ld\n", txt, grid->size, grid->numCorners);
  grid->cell_corner_lon.resize(grid->numCorners * grid->size);
  grid->cell_corner_lat.resize(grid->numCorners * grid->size);
  gridInqXbounds(grid->gridID, grid->cell_corner_lon.data());
  gridInqYbounds(grid->gridID, grid->cell_corner_lat.data());

  cdo_grid_to_radian(gridID, CDI_XAXIS, grid->cell_corner_lon, "grid corner lon");
  cdo_grid_to_radian(gridID, CDI_YAXIS, grid->cell_corner_lat, "grid corner lat");

  if (lgrid_destroy) gridDestroy(gridID);

  return grid;
}

static void
grid_delete(grid_type *grid)
{
  if (grid) delete grid;
}

static void
boundbox_from_corners1r(long ic, long nc, Varray<double> const &cornerLon, Varray<double> const &cornerLat, float *bound_box)
{
  auto inc = ic * nc;

  auto clat = cornerLat[inc];
  auto clon = cornerLon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for (long j = 1; j < nc; ++j)
  {
    clat = cornerLat[inc + j];
    clon = cornerLon[inc + j];

    if (clat < bound_box[0]) bound_box[0] = clat;
    if (clat > bound_box[1]) bound_box[1] = clat;
    if (clon < bound_box[2]) bound_box[2] = clon;
    if (clon > bound_box[3]) bound_box[3] = clon;
  }

  /*
  if ( std::fabs(bound_box[3] - bound_box[2]) > PI )
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }
  */
}

static void
boundbox_from_corners(long size, long nc, Varray<double> const &cornerLon, Varray<double> const &cornerLat,
                      Varray<float> &bound_box)
{
  for (long i = 0; i < size; ++i)
  {
    auto i4 = i << 2;  // *4
    auto inc = i * nc;
    auto clat = cornerLat[inc];
    auto clon = cornerLon[inc];
    bound_box[i4] = clat;
    bound_box[i4 + 1] = clat;
    bound_box[i4 + 2] = clon;
    bound_box[i4 + 3] = clon;
    for (long j = 1; j < nc; ++j)
    {
      clat = cornerLat[inc + j];
      clon = cornerLon[inc + j];
      if (clat < bound_box[i4]) bound_box[i4] = clat;
      if (clat > bound_box[i4 + 1]) bound_box[i4 + 1] = clat;
      if (clon < bound_box[i4 + 2]) bound_box[i4 + 2] = clon;
      if (clon > bound_box[i4 + 3]) bound_box[i4 + 3] = clon;
    }
  }
}

static cellsearch_type *
cellsearch_new(grid_type *srcGrid, grid_type *tgtGrid)
{
  cellsearch_type *cellsearch = new cellsearch_type;

  cellsearch->srcGrid = srcGrid;
  cellsearch->tgtGrid = tgtGrid;

  cellsearch->src_cell_bound_box.resize(4 * srcGrid->size);

  boundbox_from_corners(srcGrid->size, srcGrid->numCorners, srcGrid->cell_corner_lon, srcGrid->cell_corner_lat,
                        cellsearch->src_cell_bound_box);

  return cellsearch;
}

static void
cellsearch_delete(cellsearch_type *cellsearch)
{
  if (cellsearch) delete cellsearch;
}

static long
search_cells(cellsearch_type const *cellsearch, long tgtCellIndex, long *srchIndices)
{
  const grid_type *srcGrid = cellsearch->srcGrid;
  const grid_type *tgtGrid = cellsearch->tgtGrid;
  auto const &src_cell_bound_box = cellsearch->src_cell_bound_box;

  float tgt_cell_bound_box[4];
  boundbox_from_corners1r(tgtCellIndex, tgtGrid->numCorners, tgtGrid->cell_corner_lon, tgtGrid->cell_corner_lat,
                          tgt_cell_bound_box);

  auto bound_box_lat1 = tgt_cell_bound_box[0];
  auto bound_box_lat2 = tgt_cell_bound_box[1];
  auto bound_box_lon1 = tgt_cell_bound_box[2];
  auto bound_box_lon2 = tgt_cell_bound_box[3];

  long numSearchCells = 0;
  for (long srcCellIndex = 0; srcCellIndex < srcGrid->size; ++srcCellIndex)
  {
    auto srcCellIndexM4 = srcCellIndex << 2;
    if ((src_cell_bound_box[srcCellIndexM4 + 2] <= bound_box_lon2) && (src_cell_bound_box[srcCellIndexM4 + 3] >= bound_box_lon1))
    {
      if ((src_cell_bound_box[srcCellIndexM4] <= bound_box_lat2) && (src_cell_bound_box[srcCellIndexM4 + 1] >= bound_box_lat1))
      {
        srchIndices[numSearchCells] = srcCellIndex;
        numSearchCells++;
      }
    }
  }

  return numSearchCells;
}

static void
cell_search(int gridIDsrc, int gridIDtgt)
{
  grid_type *srcGrid = grid_new(gridIDsrc, "source");
  grid_type *tgtGrid = grid_new(gridIDtgt, "target");

  std::vector<long> srchIndices(srcGrid->size);

  cellsearch_type *cellsearch = cellsearch_new(srcGrid, tgtGrid);

  for (long tgtCellIndex = 0; tgtCellIndex < tgtGrid->size; ++tgtCellIndex)
  {
    long numSearchCells = search_cells(cellsearch, tgtCellIndex, srchIndices.data());

    if (Options::cdoVerbose && numSearchCells > 0)
    {
      printf("tgt cell %ld: found %ld src cells\n", tgtCellIndex, numSearchCells);
      for (long n = 0; n < numSearchCells; ++n) printf("   %ld: %ld\n", n + 1, srchIndices[n]);
    }
  }

  cellsearch_delete(cellsearch);
  grid_delete(srcGrid);
  grid_delete(tgtGrid);
}

class Gridsearch : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Gridsearch",
    .operators = { { "testpointsearch" }, { "testcellsearch" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 0, NoRestriction },
  };
  inline static RegisterEntry<Gridsearch> registration = RegisterEntry<Gridsearch>();
  int TESTPOINTSEARCH, TESTCELLSEARCH;
  int operatorID;

  int gridID1;
  int gridID2;

public:
  void
  init() override
  {

    TESTPOINTSEARCH = module.get_id("testpointsearch");
    TESTCELLSEARCH = module.get_id("testcellsearch");

    operatorID = cdo_operator_id();

    operator_input_arg("source and target grid description file or name");
    operator_check_argc(2);

    gridID1 = cdo_define_grid(cdo_operator_argv(0));
    gridID2 = cdo_define_grid(cdo_operator_argv(1));
  }

  void
  run() override
  {
    if (operatorID == TESTCELLSEARCH) cell_search(gridID1, gridID2);
  }

  void
  close() override
  {
  }
};
