/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <algorithm>

#include <cdi.h>

#include "cdo_omp.h"
#include "cdo_options.h"
#include "cdo_rlimit.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "util_files.h"
#include "util_string.h"

namespace
{
struct GridInfo1
{
  size_t gridsize{ 0 };
  size_t nx{ 0 }, ny{ 0 };
  int gridID{ -1 };
  int gridtype{ -1 };
  bool isReg2d{ false };
};

struct DistgridInfo
{
  std::vector<size_t> cellindex;
  double lonBounds[2] = { 999.0, -999.0 };
  double latBounds[2] = { 999.0, -999.0 };
  size_t gridsize{ 0 };
  size_t nx{ 0 }, ny{ 0 };
  size_t offset{ 0 };
  int gridID{ -1 };
};
}  // namespace

static void
calc_boundbox(size_t gridsize2, size_t nv, bool withBounds, Varray<double> const &xvals2, Varray<double> const &yvals2,
              Varray<double> const &xbounds2, Varray<double> const &ybounds2, double *lonBounds, double *latBounds)
{
  constexpr double Pi = M_PI;
  constexpr double Pi2 = 2.0 * M_PI;

  double xmin = xvals2[0];
  double xmax = xvals2[0];
  double ymin = yvals2[0];
  double ymax = yvals2[0];

  for (size_t i = 0; i < gridsize2; ++i)
  {
    auto xval = xvals2[i];
    auto yval = yvals2[i];
    if ((xval - xmax) > Pi) xval -= Pi2;
    if ((xmin - xval) > Pi) xval += Pi2;

    xmin = std::min(xmin, xval);
    xmax = std::max(xmax, xval);
    ymin = std::min(ymin, yval);
    ymax = std::max(ymax, yval);

    if (withBounds)
    {
      for (size_t k = 0; k < nv; ++k)
      {
        xval = xbounds2[i * nv + k];
        yval = ybounds2[i * nv + k];
        if ((xval - xmax) > Pi) xval -= Pi2;
        if ((xmin - xval) > Pi) xval += Pi2;
        xmin = std::min(xmin, xval);
        xmax = std::max(xmax, xval);
        ymin = std::min(ymin, yval);
        ymax = std::max(ymax, yval);
      }
    }
  }

  if (xmin < -Pi)
  {
    xmin += Pi2;
    xmax += Pi2;
  }

  lonBounds[0] = xmin;
  lonBounds[1] = xmax;
  latBounds[0] = ymin;
  latBounds[1] = ymax;
}

static void
gen_dist_grids(GridInfo1 const &gridInfo1, std::vector<DistgridInfo> &distgridInfo, size_t nxvals, size_t nyvals, size_t nxblocks,
               size_t nyblocks, size_t nsplit)
{
  auto gridID1 = gridInfo1.gridID;
  auto gridtype = gridInfo1.gridtype;
  auto isReg2d = gridInfo1.isReg2d;
  auto isUnstructured = (gridtype == GRID_UNSTRUCTURED);
  auto isCurvilinear = (gridtype == GRID_CURVILINEAR);
  auto isRegular
      = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC || gridtype == GRID_PROJECTION);
  if (!isRegular && !isCurvilinear && !isUnstructured) cdo_abort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  auto nx = gridInqXsize(gridID1);
  size_t ny = isUnstructured ? 1 : gridInqYsize(gridID1);

  auto lxcoord = (gridInqXvals(gridID1, nullptr) > 0);
  auto lycoord = (gridInqYvals(gridID1, nullptr) > 0);
  auto withBounds = (!isRegular && gridHasBounds(gridID1));

  std::vector<size_t> xlsize(nxblocks), ylsize(nyblocks);

  for (size_t ix = 0; ix < nxblocks; ++ix) xlsize[ix] = nxvals;
  if (nx % nxblocks != 0) xlsize[nxblocks - 1] = nx - (nxblocks - 1) * nxvals;
  if (Options::cdoVerbose)
    for (size_t ix = 0; ix < nxblocks; ++ix) cdo_print("xblock %zu: size=%zu", ix, xlsize[ix]);

  for (size_t iy = 0; iy < nyblocks; ++iy) ylsize[iy] = nyvals;
  if (ny % nyblocks != 0) ylsize[nyblocks - 1] = ny - (nyblocks - 1) * nyvals;
  if (Options::cdoVerbose)
    for (size_t iy = 0; iy < nyblocks; ++iy) cdo_print("yblock %zu: size=%zu", iy, ylsize[iy]);

  auto nxvmax = std::max(nxvals, xlsize[nxblocks - 1]);
  auto nyvmax = std::max(nyvals, ylsize[nyblocks - 1]);

  Varray<double> xpvals, ypvals;
  Varray<double> xvals, yvals, xbounds, ybounds;
  Varray<double> xvals2, yvals2, xbounds2, ybounds2;

  if (lxcoord)
  {
    xvals.resize(isRegular ? nx : nx * ny);
    gridInqXvals(gridID1, xvals.data());

    if (!isRegular) xvals2.resize(nxvmax * nyvmax);
  }

  if (lycoord)
  {
    yvals.resize(isRegular ? ny : nx * ny);
    gridInqYvals(gridID1, yvals.data());

    if (!isRegular) yvals2.resize(nxvmax * nyvmax);
  }

  size_t nv = 0;
  if (withBounds)
  {
    if (!isRegular)
    {
      nv = gridInqNvertex(gridID1);
      xbounds.resize(nx * ny * nv);
      ybounds.resize(nx * ny * nv);
      xbounds2.resize(nxvmax * nyvmax * nv);
      ybounds2.resize(nxvmax * nyvmax * nv);
    }

    gridInqXbounds(gridID1, xbounds.data());
    gridInqYbounds(gridID1, ybounds.data());
  }

  size_t index = 0;
  for (size_t iy = 0; iy < nyblocks; ++iy)
    for (size_t ix = 0; ix < nxblocks; ++ix)
    {
      auto offset = iy * nyvals * nx + ix * nxvals;
      auto gridsize2 = xlsize[ix] * ylsize[iy];
      if (!isReg2d) distgridInfo[index].cellindex.resize(gridsize2);
      auto &cellindex = distgridInfo[index].cellindex;

      distgridInfo[index].nx = xlsize[ix];
      distgridInfo[index].ny = ylsize[iy];
      distgridInfo[index].offset = offset;
      auto &lonBounds = distgridInfo[index].lonBounds;
      auto &latBounds = distgridInfo[index].latBounds;

      gridsize2 = 0;
      // printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
      for (size_t j = 0; j < ylsize[iy]; ++j)
        for (size_t i = 0; i < xlsize[ix]; ++i)
        {
          // printf(">> %d %d %d\n", j, i, offset + j*nx + i);
          if (!isRegular)
          {
            if (lxcoord) xvals2[gridsize2] = xvals[offset + j * nx + i];
            if (lycoord) yvals2[gridsize2] = yvals[offset + j * nx + i];
            if (lxcoord) lonBounds[0] = std::min(lonBounds[0], xvals2[gridsize2]);
            if (lxcoord) lonBounds[1] = std::max(lonBounds[1], xvals2[gridsize2]);
            if (lycoord) latBounds[0] = std::min(latBounds[0], yvals2[gridsize2]);
            if (lycoord) latBounds[1] = std::max(latBounds[1], yvals2[gridsize2]);
            if (withBounds)
            {
              for (size_t k = 0; k < nv; ++k)
              {
                xbounds2[gridsize2 * nv + k] = xbounds[(offset + j * nx + i) * nv + k];
                ybounds2[gridsize2 * nv + k] = ybounds[(offset + j * nx + i) * nv + k];
                lonBounds[0] = std::min(lonBounds[0], xbounds2[gridsize2 * nv + k]);
                lonBounds[1] = std::max(lonBounds[1], xbounds2[gridsize2 * nv + k]);
                latBounds[0] = std::min(latBounds[0], ybounds2[gridsize2 * nv + k]);
                latBounds[1] = std::max(latBounds[1], ybounds2[gridsize2 * nv + k]);
              }
            }
          }
          if (!isReg2d) cellindex[gridsize2] = offset + j * nx + i;
          gridsize2++;
        }
      // printf("gridsize2 %d\n", gridsize2);

      if (!isRegular && lxcoord && lycoord)
        calc_boundbox(gridsize2, nv, withBounds, xvals2, yvals2, xbounds2, ybounds2, lonBounds, latBounds);

      auto gridID2 = gridCreate(gridtype, gridsize2);
      if (gridtype != GRID_UNSTRUCTURED)
      {
        gridDefXsize(gridID2, xlsize[ix]);
        gridDefYsize(gridID2, ylsize[iy]);

        gridDefNP(gridID2, gridInqNP(gridID1));
      }

      if (withBounds) gridDefNvertex(gridID2, nv);

      cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);
      cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, gridID2);
      cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, gridID2);
      cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_REFERENCEURI, gridID2);
      cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_UUID, gridID2);

      grid_copy_names(gridID1, gridID2);

      if (gridtype == GRID_PROJECTION) grid_copy_mapping(gridID1, gridID2);

      if (isRegular)
      {
        if (lxcoord) gridDefXvals(gridID2, &xvals[ix * nxvals]);
        if (lycoord) gridDefYvals(gridID2, &yvals[iy * nyvals]);
      }
      else
      {
        if (lxcoord) gridDefXvals(gridID2, xvals2.data());
        if (lycoord) gridDefYvals(gridID2, yvals2.data());
        if (withBounds)
        {
          gridDefXbounds(gridID2, xbounds2.data());
          gridDefYbounds(gridID2, ybounds2.data());
        }
      }

      auto projID1 = gridInqProj(gridID1);
      if (projID1 != CDI_UNDEFID && gridInqType(projID1) == GRID_PROJECTION)
      {
        auto projID2 = gridCreate(GRID_PROJECTION, gridsize2);
        gridDefXsize(projID2, xlsize[ix]);
        gridDefYsize(projID2, ylsize[iy]);

        grid_copy_names(projID1, projID2);
        grid_copy_mapping(projID1, projID2);

        auto lxpcoord = (gridInqXvals(projID1, nullptr) > 0);
        if (lxpcoord)
        {
          if (!xpvals.size())
          {
            xpvals.resize(nx);
            gridInqXvals(projID1, xpvals.data());
          }
          gridDefXvals(projID2, &xpvals[ix * nxvals]);
        }
        auto lypcoord = (gridInqYvals(projID1, nullptr) > 0);
        if (lypcoord)
        {
          if (!ypvals.size())
          {
            ypvals.resize(ny);
            gridInqYvals(projID1, ypvals.data());
          }
          gridDefYvals(projID2, &ypvals[iy * nyvals]);
        }

        gridDefProj(gridID2, projID2);
      }

      distgridInfo[index].gridID = gridID2;
      distgridInfo[index].gridsize = gridsize2;

      index++;
      if (index > nsplit) cdo_abort("Internal problem, index exceeded bounds!");
    }

  auto numBlocks = index;

  for (size_t i = 0; i < numBlocks; ++i)
  {
    auto &lons = distgridInfo[i].lonBounds;
    auto &lats = distgridInfo[i].latBounds;
    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID1, CDI_XAXIS, 2, lons, "lon bounds");
    cdo_grid_to_degree(gridID1, CDI_YAXIS, 2, lats, "lat bounds");
  }
}

static void
dist_cells_reg2d(Field const &field1, Field &field2, DistgridInfo const &distgridInfo, size_t nlon)
{
  auto nx = distgridInfo.nx;
  auto ny = distgridInfo.ny;

  for (size_t j = 0; j < ny; ++j)
  {
    auto offset1 = distgridInfo.offset + j * nlon;
    auto offset2 = j * nx;

    auto func = [&](auto const &v1, auto &v2)
    {
      for (size_t i = 0; i < nx; ++i) { v2[offset2 + i] = v1[offset1 + i]; }
    };
    field_operation2(func, field1, field2);
  }
}

static void
dist_cells(Field const &field1, Field &field2, std::vector<size_t> const &cellIndex)
{
  auto func = [&](auto const &v1, auto &v2, auto n)
  {
    for (size_t i = 0; i < n; ++i) { v2[i] = v1[cellIndex[i]]; }
  };
  field_operation2(func, field1, field2, field2.size);
}

class Distgrid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Distgrid",
    .operators = { { "distgrid", DistgridHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, OBASE, OnlyFirst },
  };
  inline static RegisterEntry<Distgrid> registration = RegisterEntry<Distgrid>(module);
  size_t MaxBlocks = 99999;
  int gridID1;
  int gridtype = -1;
  int numGrids;

  CdoStreamID streamID1;
  std::vector<CdoStreamID> streamIDs;

  size_t nsplit;

  std::vector<int> vlistIDs;
  std::vector<GridInfo1> gridInfo1;
  VarList varList1;
  std::vector<std::vector<DistgridInfo>> distgridInfo;

public:
  void
  init() override
  {
    operator_input_arg("nxblocks, [nyblocks]");
    if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
    if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
    size_t nxblocks = parameter_to_int(cdo_operator_argv(0));
    size_t nyblocks = 1;
    if (cdo_operator_argc() == 2) nyblocks = parameter_to_int(cdo_operator_argv(1));

    if (nxblocks == 0) cdo_abort("nxblocks has to be greater than 0!");
    if (nyblocks == 0) cdo_abort("nyblocks has to be greater than 0!");

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    numGrids = vlistNumGrids(vlistID1);
    gridInfo1.resize(numGrids);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID = vlistGrid(vlistID1, index);
      gridInfo1[index].gridID = gridID;
      gridInfo1[index].gridtype = gridInqType(gridID);
      gridInfo1[index].gridsize = gridInqSize(gridID);
      gridInfo1[index].nx = gridInqXsize(gridID);
      gridInfo1[index].ny = gridInqYsize(gridID);
      gridInfo1[index].isReg2d = (gridInfo1[index].gridtype != GRID_UNSTRUCTURED);
    }

    {
      int index;
      for (index = 0; index < numGrids; ++index)
      {
        gridID1 = vlistGrid(vlistID1, index);
        gridtype = gridInqType(gridID1);
        if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED
            || gridtype == GRID_PROJECTION || (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0))
          break;
      }

      if (index == numGrids)
        cdo_abort("No Lon/Lat, Gaussian, curvilinear or generic grid found (%s data unsupported)!", gridNamePtr(gridtype));
    }

    auto isUnstructured = (gridtype == GRID_UNSTRUCTURED);

    gridID1 = vlistGrid(vlistID1, 0);
    auto gridsize = gridInqSize(gridID1);
    auto nx = gridInqXsize(gridID1);
    size_t ny = isUnstructured ? 1 : gridInqYsize(gridID1);
    for (int i = 1; i < numGrids; ++i)
    {
      gridID1 = vlistGrid(vlistID1, i);
      if (gridsize != gridInqSize(gridID1)) cdo_abort("Gridsize must not change!");
    }

    if (nxblocks > nx)
    {
      cdo_print("nxblocks (%zu) greater than nx (%zu), set to %zu!", nxblocks, nx, nx);
      nxblocks = nx;
    }
    if (nyblocks > ny)
    {
      cdo_print("nyblocks (%zu) greater than ny (%zu), set to %zu!", nyblocks, ny, ny);
      nyblocks = ny;
    }

    auto xinc = nx / nxblocks;
    auto yinc = ny / nyblocks;
    if (nx % xinc && nx % (xinc + 1) && nxblocks * (xinc + 1) <= nx) xinc++;
    if (ny % yinc && ny % (yinc + 1) && nyblocks * (yinc + 1) <= ny) yinc++;

    nsplit = nxblocks * nyblocks;
    if (nsplit > MaxBlocks) cdo_abort("Too many blocks (max = %d)!", MaxBlocks);

    cdo::set_numfiles(nsplit + 8);

    varList1 = VarList(vlistID1);

    vlistIDs = std::vector<int>(nsplit);
    streamIDs = std::vector<CdoStreamID>(nsplit);

    distgridInfo = std::vector<std::vector<DistgridInfo>>(numGrids);

    for (int i = 0; i < numGrids; ++i) distgridInfo[i].resize(nsplit);

    for (size_t index = 0; index < nsplit; ++index) vlistIDs[index] = vlistDuplicate(vlistID1);

    if (Options::cdoVerbose) cdo_print("numGrids=%d  nsplit=%zu", numGrids, nsplit);

    for (int i = 0; i < numGrids; ++i)
    {
      gen_dist_grids(gridInfo1[i], distgridInfo[i], xinc, yinc, nxblocks, nyblocks, nsplit);
      /*
      if ( Options::cdoVerbose )
        for ( size_t index = 0; index < nsplit; index++ )
          cdo_print("Block %d,  gridID %d,  gridsize %zu", index+1,
      distgridInfo[i].gridID[index], gridInqSize(grids[i].gridID[index]));
      */
      for (size_t index = 0; index < nsplit; ++index) vlistChangeGridIndex(vlistIDs[index], i, distgridInfo[i][index].gridID);
    }

    auto fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    if (Options::test && numGrids == 1 && isUnstructured)
    {
      printf("#MGF\n");
      printf("numfiles=%zu\n", nsplit);
    }

    for (size_t index = 0; index < nsplit; ++index)
    {
      auto fileName = cdo_get_obase() + string_format("%05ld", (long) index);
      if (fileSuffix.size() > 0) fileName += fileSuffix;

      streamIDs[index] = open_write(fileName);
      cdo_def_vlist(streamIDs[index], vlistIDs[index]);

      if (Options::test && numGrids == 1 && isUnstructured)
      {
        auto gridsize2 = distgridInfo[0][index].gridsize;
        auto offset = distgridInfo[0][index].offset;
        auto const &lons = distgridInfo[0][index].lonBounds;
        auto const &lats = distgridInfo[0][index].latBounds;
        printf("--- # %zu\n", index + 1);
        printf("filename=%s\n", fileName.c_str());
        printf("numcells=%zu\n", gridsize2);
        printf("offset=%zu\n", offset);
        printf("boundbox=%g/%g/%g/%g\n", lons[0], lons[1], lats[0], lats[1]);
      }
    }
  }

  void
  run() override
  {
    Field field1;
    std::vector<Field> field2vec(Threading::ompNumMaxThreads);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      for (auto const &streamID : streamIDs) cdo_def_timestep(streamID, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varIDx, levelIDx] = cdo_inq_field(streamID1);
        int varID = varIDx;  // needed for omp loop with intel icpx 2022.0.0
        int levelID = levelIDx;
        field1.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field1);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for (size_t index = 0; index < nsplit; ++index)
        {
          auto var = varList1.vars[varID];
          auto &distgrid = distgridInfo[0][index];
          var.gridID = distgrid.gridID;
          var.gridsize = distgrid.gridsize;

          auto ompthID = cdo_omp_get_thread_num();
          auto &field2 = field2vec[ompthID];
          field2.init(var);

          if (gridInfo1[0].isReg2d)
            dist_cells_reg2d(field1, field2, distgrid, gridInfo1[0].nx);
          else
            dist_cells(field1, field2, distgrid.cellindex);

          if (field1.numMissVals) field_num_mv(field2);

          cdo_def_field(streamIDs[index], varID, levelID);
          cdo_write_field(streamIDs[index], field2);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);

    for (auto const &streamID : streamIDs) cdo_stream_close(streamID);
    for (auto const &vlistID : vlistIDs) vlistDestroy(vlistID);

    for (int i = 0; i < numGrids; ++i)
      for (auto const &distgrid : distgridInfo[i]) gridDestroy(distgrid.gridID);
  }
};
