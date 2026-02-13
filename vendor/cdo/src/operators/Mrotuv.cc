/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Mrotuv      mrotuv          Forward rotation for MPIOM data
*/

#include <cdi.h>

#include <mpim_grid.h>
#include "varray.h"
#include "cdo_options.h"
#include "process_int.h"
#include "cdi_lockedIO.h"
#include "matrix_view.h"

void
rotate_uv(Varray<double> &u_i_v, Varray<double> &v_j_v, long nx, long ny, Varray<double> &lon_v, Varray<double> &lat_v,
          Varray<double> &u_lon_v, Varray<double> &v_lat_v)
{
  /*
    in      :: u_i[ny][nx], v_j[ny][nx]      vector components in i-j-direction
    in      :: lat[ny][nx], lon[ny][nx]      latitudes and longitudes
    out     :: u_lon[ny][nx], v_lat[ny][nx]  vector components in lon-lat direction
  */
  constexpr double pi = 3.14159265359;
  MatrixView<double> lon(lon_v.data(), ny, nx);
  MatrixView<double> lat(lat_v.data(), ny, nx);
  MatrixView<double> u_i(u_i_v.data(), ny, nx);
  MatrixView<double> v_j(v_j_v.data(), ny, nx);
  MatrixView<double> u_lon(u_lon_v.data(), ny, nx);
  MatrixView<double> v_lat(v_lat_v.data(), ny, nx);

  // specification whether change in sign is needed for the input arrays
  constexpr auto change_sign_u = false;
  constexpr auto change_sign_v = true;

  // initialization
  v_lat_v.assign(nx * ny, 0.0);
  u_lon_v.assign(nx * ny, 0.0);

  // rotation
  for (long j = 0; j < ny; ++j)
    for (long i = 0; i < nx; ++i)
    {
      auto ip1 = i + 1;
      auto im1 = i - 1;
      auto jp1 = j + 1;
      auto jm1 = j - 1;
      if (ip1 >= nx) ip1 = 0;  // the 0-meridian
      if (im1 < 0) im1 = nx - 1;
      if (jp1 >= ny) jp1 = j;  // treatment of the last..
      if (jm1 < 0) jm1 = j;    // .. and the fist grid-row

      // difference in latitudes
      auto dlat_i = lat[j][ip1] - lat[j][im1];
      auto dlat_j = lat[jp1][i] - lat[jm1][i];

      // difference in longitudes
      auto dlon_i = lon[j][ip1] - lon[j][im1];
      if (dlon_i > pi) dlon_i -= 2 * pi;
      if (dlon_i < (-pi)) dlon_i += 2 * pi;
      auto dlon_j = lon[jp1][i] - lon[jm1][i];
      if (dlon_j > pi) dlon_j -= 2 * pi;
      if (dlon_j < (-pi)) dlon_j += 2 * pi;

      auto lat_factor = std::cos(lat[j][i]);
      dlon_i = dlon_i * lat_factor;
      dlon_j = dlon_j * lat_factor;

      // projection by scalar product
      u_lon[j][i] = u_i[j][i] * dlon_i + v_j[j][i] * dlat_i;
      v_lat[j][i] = u_i[j][i] * dlon_j + v_j[j][i] * dlat_j;

      auto dist_i = std::sqrt(dlon_i * dlon_i + dlat_i * dlat_i);
      auto dist_j = std::sqrt(dlon_j * dlon_j + dlat_j * dlat_j);

      if (std::fabs(dist_i) > 0.0 && std::fabs(dist_j) > 0.0)
      {
        u_lon[j][i] /= dist_i;
        v_lat[j][i] /= dist_j;
      }
      else
      {
        u_lon[j][i] = 0.0;
        v_lat[j][i] = 0.0;
      }

      // velocity vector lengths
      auto absold = std::sqrt(u_i[j][i] * u_i[j][i] + v_j[j][i] * v_j[j][i]);
      auto absnew = std::sqrt(u_lon[j][i] * u_lon[j][i] + v_lat[j][i] * v_lat[j][i]);

      u_lon[j][i] *= absold;
      v_lat[j][i] *= absold;

      if (absnew > 0.0)
      {
        u_lon[j][i] /= absnew;
        v_lat[j][i] /= absnew;
      }
      else
      {
        u_lon[j][i] = 0.0;
        v_lat[j][i] = 0.0;
      }

      // change sign
      if (change_sign_u) u_lon[j][i] *= -1;
      if (change_sign_v) v_lat[j][i] *= -1;

      if (Options::cdoVerbose)
      {
        absold = std::sqrt(u_i[j][i] * u_i[j][i] + v_j[j][i] * v_j[j][i]);
        absnew = std::sqrt(u_lon[j][i] * u_lon[j][i] + v_lat[j][i] * v_lat[j][i]);

        if (i % 20 == 0 && j % 20 == 0 && absold > 0.0)
        {
          printf("(absold,absnew) %ld %ld %g %g %g %g %g %g\n", j + 1, i + 1, absold, absnew, u_i[j][i], v_j[j][i], u_lon[j][i],
                 v_lat[j][i]);

          // test orthogonality
          if ((dlon_i * dlon_j + dlat_j * dlat_i) > 0.1)
            std::fprintf(stderr, "orthogonal? %ld %ld %g\n", j + 1, i + 1, (dlon_i * dlon_j + dlat_j * dlat_i));
        }
      }
    }
}

void
p_to_uv_grid(long nlon, long nlat, Varray<double> &grid1x_v, Varray<double> &grid1y_v, Varray<double> &gridux_v,
             Varray<double> &griduy_v, Varray<double> &gridvx_v, Varray<double> &gridvy_v)
{
  MatrixView<double> grid1x(grid1x_v.data(), nlat, nlon);
  MatrixView<double> grid1y(grid1y_v.data(), nlat, nlon);
  MatrixView<double> gridux(gridux_v.data(), nlat, nlon);
  MatrixView<double> griduy(griduy_v.data(), nlat, nlon);
  MatrixView<double> gridvx(gridvx_v.data(), nlat, nlon);
  MatrixView<double> gridvy(gridvy_v.data(), nlat, nlon);

  // interpolate scalar to u points
  for (long j = 0; j < nlat; ++j)
    for (long i = 0; i < nlon; ++i)
    {
      auto ip1 = i + 1;
      if (ip1 > nlon - 1) ip1 = 0;

      gridux[j][i] = (grid1x[j][i] + grid1x[j][ip1]) * 0.5;
      if ((grid1x[j][i] > 340 && grid1x[j][ip1] < 20) || (grid1x[j][i] < 20 && grid1x[j][ip1] > 340))
      {
        gridux[j][i] += (gridux[j][i] < 180) ? 180 : -180;
      }

      griduy[j][i] = (grid1y[j][i] + grid1y[j][ip1]) * 0.5;
    }

  // interpolate scalar to v points
  for (long j = 0; j < nlat; ++j)
    for (long i = 0; i < nlon; ++i)
    {
      auto jp1 = j + 1;
      if (jp1 > nlat - 1) jp1 = nlat - 1;

      gridvx[j][i] = (grid1x[j][i] + grid1x[jp1][i]) * 0.5;
      if ((grid1x[j][i] > 340 && grid1x[jp1][i] < 20) || (grid1x[j][i] < 20 && grid1x[jp1][i] > 340))
      {
        gridvx[j][i] += (gridvx[j][i] < 180) ? 180 : -180;
      }

      gridvy[j][i] = (grid1y[j][i] + grid1y[jp1][i]) * 0.5;
    }
}

class Mrotuv : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Mrotuv",
    .operators = { { "mrotuv" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 2, NoRestriction },
  };
  inline static RegisterEntry<Mrotuv> registration = RegisterEntry<Mrotuv>();

private:
  size_t numMissVals1 = 0, numMissVals2 = 0;
  int uid = -1, vid = -1;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};

  double missval1{};
  double missval2{};

  int nlevs{};
  size_t gridsize{};

  size_t nlon{};
  size_t nlat{};

  Varray<double> grid1x{};
  Varray<double> grid1y{};
  Varray<double> gridux{};
  Varray<double> griduy{};
  Varray<double> gridvx{};
  Varray<double> gridvy{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    VarList varList1(vlistID1);

    auto numVars = varList1.numVars();
    for (int varid = 0; varid < numVars; varid++)
    {
      auto code = varList1.vars[varid].code;
      if (code == 3 || code == 131) uid = varid;
      if (code == 4 || code == 132) vid = varid;
    }

    if (uid == -1 || vid == -1)
    {
      if (numVars == 2)
      {
        uid = 0;
        vid = 1;
      }
      else
        cdo_abort("U and V not found in %s", cdo_get_stream_name(0));
    }

    auto const &varU = varList1.vars[uid];
    auto const &varV = varList1.vars[vid];
    if (varU.nlevels != varV.nlevels) cdo_abort("U and V have different number of levels!");

    auto gridID1 = varU.gridID;
    gridsize = varU.gridsize;
    if (varU.gridID != varV.gridID) cdo_abort("Input grids differ!");

    auto gridType = varU.gridType;
    if (gridType != GRID_LONLAT && gridType != GRID_GAUSSIAN && gridType != GRID_CURVILINEAR)
      cdo_abort("Grid %s unsupported!", gridNamePtr(gridType));

    if (gridType != GRID_CURVILINEAR) gridID1 = gridToCurvilinear(gridID1);

    if (gridsize != gridInqSize(gridID1)) cdo_abort("Internal problem: gridsize changed!");

    nlon = gridInqXsize(gridID1);
    nlat = gridInqYsize(gridID1);

    grid1x.resize(gridsize);
    grid1y.resize(gridsize);
    gridux.resize(gridsize);
    griduy.resize(gridsize);
    gridvx.resize(gridsize);
    gridvy.resize(gridsize);

    gridInqXvals(gridID1, grid1x.data());
    gridInqYvals(gridID1, grid1y.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID1, CDI_XAXIS, grid1x, "grid center lon");
    cdo_grid_to_degree(gridID1, CDI_YAXIS, grid1y, "grid center lat");

    p_to_uv_grid(nlon, nlat, grid1x, grid1y, gridux, griduy, gridvx, gridvy);

    auto gridIDu = gridCreate(GRID_CURVILINEAR, nlon * nlat);
    cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridIDu);
    gridDefXsize(gridIDu, nlon);
    gridDefYsize(gridIDu, nlat);
    gridDefXvals(gridIDu, gridux.data());
    gridDefYvals(gridIDu, griduy.data());

    auto gridIDv = gridCreate(GRID_CURVILINEAR, nlon * nlat);
    cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridIDv);
    gridDefXsize(gridIDv, nlon);
    gridDefYsize(gridIDv, nlat);
    gridDefXvals(gridIDv, gridvx.data());
    gridDefYvals(gridIDv, gridvy.data());

    for (size_t i = 0; i < gridsize; ++i)
    {
      grid1x[i] *= DEG2RAD;
      grid1y[i] *= DEG2RAD;
    }

    vlistClearFlag(vlistID1);
    for (int lid = 0; lid < nlevs; lid++) vlistDefFlag(vlistID1, uid, lid, true);
    auto vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistChangeVarGrid(vlistID2, 0, gridIDu);

    vlistClearFlag(vlistID1);
    for (int lid = 0; lid < nlevs; lid++) vlistDefFlag(vlistID1, vid, lid, true);
    auto vlistID3 = vlistCreate();
    cdo_vlist_copy_flag(vlistID3, vlistID1);
    vlistChangeVarGrid(vlistID3, 0, gridIDv);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID2 = cdo_open_write(1);
    streamID3 = cdo_open_write(2);

    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_vlist(streamID3, vlistID3);

    missval1 = varU.missval;
    missval2 = varV.missval;
  }

  void
  run() override
  {
    Varray<double> ufield_v(gridsize);
    Varray<double> vfield_v(gridsize);
    MatrixView<double> ufield(ufield_v.data(), nlat, nlon);
    MatrixView<double> vfield(vfield_v.data(), nlat, nlon);

    Varray2D<double> urfield(nlevs);
    for (auto &ur : urfield) { ur.resize(gridsize); }
    Varray2D<double> vrfield(nlevs);
    for (auto &vr : vrfield) { vr.resize(gridsize); }

    Varray2D<double> uhelp(nlat, Varray<double>(nlon + 2));
    Varray2D<double> vhelp(nlat, Varray<double>(nlon + 2));
    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        if (varID == uid) cdo_read_field(streamID1, urfield[levelID].data(), &numMissVals1);
        if (varID == vid) cdo_read_field(streamID1, vrfield[levelID].data(), &numMissVals2);
      }

      for (int levelID = 0; levelID < nlevs; ++levelID)
      {
        // remove missing values
        if (numMissVals1 || numMissVals2)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(urfield[levelID][i], missval1)) urfield[levelID][i] = 0.0;
            if (fp_is_equal(vrfield[levelID][i], missval2)) vrfield[levelID][i] = 0.0;
          }
        }

        // rotate
        rotate_uv(urfield[levelID], vrfield[levelID], nlon, nlat, grid1x, grid1y, ufield_v, vfield_v);

        // load to a help field
        for (size_t j = 0; j < nlat; ++j)
          for (size_t i = 0; i < nlon; ++i)
          {
            uhelp[j][i + 1] = ufield[j][i];
            vhelp[j][i + 1] = vfield[j][i];
          }

        // make help field cyclic
        for (size_t j = 0; j < nlat; ++j)
        {
          uhelp[j][0] = uhelp[j][nlon];
          uhelp[j][nlon + 1] = uhelp[j][1];
          vhelp[j][0] = vhelp[j][nlon];
          vhelp[j][nlon + 1] = vhelp[j][1];
        }

        // interpolate on u/v points
        for (size_t j = 0; j < nlat; ++j)
          for (size_t i = 0; i < nlon; ++i) { ufield[j][i] = (uhelp[j][i + 1] + uhelp[j][i + 2]) * 0.5; }

        for (size_t j = 0; j < nlat - 1; ++j)
          for (size_t i = 0; i < nlon; ++i) { vfield[j][i] = (vhelp[j][i + 1] + vhelp[j + 1][i + 1]) * 0.5; }

        for (size_t i = 0; i < nlon; ++i) { vfield[nlat - 1][i] = vhelp[nlat - 1][i + 1]; }

        cdo_def_field(streamID2, 0, levelID);
        cdo_write_field(streamID2, ufield_v.data(), numMissVals1);
        cdo_def_field(streamID3, 0, levelID);
        cdo_write_field(streamID3, vfield_v.data(), numMissVals2);
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
