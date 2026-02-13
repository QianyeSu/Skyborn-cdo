/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Mrotuvb     mrotuvb          Backward rotation for MPIOM data
*/

#include <cdi.h>

#include "varray.h"
#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "matrix_view.h"

/*
!----------------------------------------------------------------------
!
!     rotation of vectors: In ocean models with rotated grids velocity
!     vectors are given in the direction of grid lines and rows. They
!     have to be rotated in latitudinal and longitudinal direction.
!
!     Note: This routine assumes positive meridional flow for a flow
!           from grid point(i,j) to grid point(i,j+1) and positive
!           zonal flow for a flow from grid point(i,j) to point(i+1,j).
!           this is not the case for mpi-om!
!
!           If this routine is used to rotate data of mpi-om, the
!           logical change_sign_v needs to be true.
!
!  j. jungclaus: 22.01.04:
!    Note here for the coupling fields u_i,v_j are on the non-verlapping
!    (ie-2=ix) grid, furthermore, the velocity fields were previously
!    interpolated onto the scalar points!
!----------------------------------------------------------------------
*/

void
rotate_uv2(Varray<double> &u_i_v, Varray<double> &v_j_v, long nx, long ny, Varray<double> &lon_v, Varray<double> &lat_v,
           Varray<double> &u_lon_v, Varray<double> &v_lat_v)
{
  /*
    in      :: u_i[ny][nx], v_j[ny][nx]       vector components in i-j-direction
    in      :: lat[ny][nx], lon[ny][nx]       latitudes and longitudes
    out     :: u_lon[ny][nx], v_lat[ny][nx]   vector components in lon-lat direction
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

  // change sign
  if (change_sign_u)
    for (long i = 0; i < nx * ny; ++i) u_i_v[i] *= -1;

  if (change_sign_v)
    for (long i = 0; i < nx * ny; ++i) v_j_v[i] *= -1;

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

      if (std::fabs(dist_i) > 0 && std::fabs(dist_j) > 0)
      {
        u_lon[j][i] /= dist_i;
        v_lat[j][i] /= dist_j;
      }
      else
      {
        u_lon[j][i] = 0.0;
        v_lat[j][i] = 0.0;
      }

      if (Options::cdoVerbose)
      {
        // velocity vector lengths
        auto absold = std::sqrt(u_i[j][i] * u_i[j][i] + v_j[j][i] * v_j[j][i]);
        auto absnew = std::sqrt(u_lon[j][i] * u_lon[j][i] + v_lat[j][i] * v_lat[j][i]);

        if (i % 20 == 0 && j % 20 == 0 && absold > 0)
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

static void
uv_to_p_grid(size_t nlon, size_t nlat, Varray<double> &grid1x_v, Varray<double> &grid1y_v, Varray<double> &grid2x_v,
             Varray<double> &grid2y_v, Varray<double> &grid3x_v, Varray<double> &grid3y_v)
{
  MatrixView<double> grid1x(grid1x_v.data(), nlat, nlon);
  MatrixView<double> grid1y(grid1y_v.data(), nlat, nlon);
  MatrixView<double> grid2x(grid2x_v.data(), nlat, nlon);
  MatrixView<double> grid2y(grid2y_v.data(), nlat, nlon);
  MatrixView<double> grid3x(grid3x_v.data(), nlat, nlon);
  MatrixView<double> grid3y(grid3y_v.data(), nlat, nlon);

  Varray2D<double> gxhelp(nlat, Varray<double>(nlon + 2)), gyhelp(nlat, Varray<double>(nlon + 2));

  // load to a help field
  for (size_t j = 0; j < nlat; ++j)
    for (size_t i = 0; i < nlon; ++i)
    {
      gxhelp[j][i + 1] = grid1x[j][i];
      gyhelp[j][i + 1] = grid1y[j][i];
    }

  // make help field cyclic
  for (size_t j = 0; j < nlat; ++j)
  {
    gxhelp[j][0] = gxhelp[j][nlon];
    gxhelp[j][nlon + 1] = gxhelp[j][1];
    gyhelp[j][0] = gyhelp[j][nlon];
    gyhelp[j][nlon + 1] = gyhelp[j][1];
  }

  // interpolate u to scalar points
  for (size_t j = 0; j < nlat; ++j)
    for (size_t i = 0; i < nlon; ++i)
    {
      grid3x[j][i] = (gxhelp[j][i] + gxhelp[j][i + 1]) * 0.5;
      if ((gxhelp[j][i] > 340 && gxhelp[j][i + 1] < 20) || (gxhelp[j][i] < 20 && gxhelp[j][i + 1] > 340))
      {
        grid3x[j][i] += (grid3x[j][i] < 180) ? 180 : -180;
      }

      grid3y[j][i] = (gyhelp[j][i] + gyhelp[j][i + 1]) * 0.5;
    }

  // load to a help field
  for (size_t j = 0; j < nlat; ++j)
    for (size_t i = 0; i < nlon; ++i)
    {
      gxhelp[j][i + 1] = grid2x[j][i];
      gyhelp[j][i + 1] = grid2y[j][i];
    }

  // make help field cyclic
  for (size_t j = 0; j < nlat; ++j)
  {
    gxhelp[j][0] = gxhelp[j][nlon];
    gxhelp[j][nlon + 1] = gxhelp[j][1];
    gyhelp[j][0] = gyhelp[j][nlon];
    gyhelp[j][nlon + 1] = gyhelp[j][1];
  }

  // interpolate v to scalar points
  for (size_t j = 1; j < nlat - 1; ++j)
    for (size_t i = 0; i < nlon; ++i)
    {
      auto gx = (gxhelp[j][i + 1] + gxhelp[j - 1][i + 1]) * 0.5;
      if ((gxhelp[j][i + 1] > 340 && gxhelp[j - 1][i + 1] < 20) || (gxhelp[j][i + 1] < 20 && gxhelp[j - 1][i + 1] > 340))
      {
        gx += (gx < 180) ? 180 : -180;
      }

      auto gy = (gyhelp[j][i + 1] + gyhelp[j - 1][i + 1]) * 0.5;

      // printf("%d %d %g %g %g %g \n", j, i, gx, gy, grid3x[j][i], grid3y[j][i]);

      auto gx2 = (gx + grid3x[j][i]) * 0.5;
      if ((gx > 340 && grid3x[j][i] < 20) || (gx < 20 && grid3x[j][i] > 340)) { gx2 += (gx2 < 180) ? 180 : -180; }

      auto gy2 = (gy + grid3y[j][i]) * 0.5;

      grid3x[j][i] = gx2;
      grid3y[j][i] = gy2;

      // printf("%d %d %g %g %g %g \n", j, i, gx2, gy2, grid3x[j][i], grid3y[j][i]);
    }
}

class Mrotuvb : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Mrotuvb",
    .operators = { { "mrotuvb", MrotuvbHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Mrotuvb> registration = RegisterEntry<Mrotuvb>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  CdoStreamID streamID3;

  int taxisID1{ CDI_UNDEFID };
  int taxisID3;

  double missval1;
  double missval2;

  size_t nlon;
  size_t nlat;

  size_t gridsize;

  bool gpint;

  Varray<double> grid3x;
  Varray<double> grid3y;

public:
  void
  init() override
  {
    gpint = !(cdo_operator_argc() == 1 && cdo_operator_argv(0) == "noint");

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);

    VarList varList1(vlistID1);
    VarList varList2(vlistID2);

    auto const &varU = varList1.vars[0];
    auto const &varV = varList2.vars[0];

    if (varList1.numVars() > 1) cdo_abort("More than one variable found in %s", cdo_get_stream_name(0));
    if (varList2.numVars() > 1) cdo_abort("More than one variable found in %s", cdo_get_stream_name(1));

    auto gridID1 = varU.gridID;
    auto gridID2 = varV.gridID;
    gridsize = varU.gridsize;
    if (gpint && gridID1 == gridID2) cdo_abort("Input grids are the same, use parameter >noint< to disable interpolation!");
    if (!gpint && gridID1 != gridID2) cdo_abort("Input grids are not the same!");
    if (gridsize != varV.gridsize) cdo_abort("Grids have different size!");

    auto gridType = varU.gridType;
    if (gridType != GRID_LONLAT && gridType != GRID_GAUSSIAN && gridType != GRID_CURVILINEAR)
      cdo_abort("Grid %s unsupported!", gridNamePtr(gridType));

    if (gridType != GRID_CURVILINEAR) gridID1 = gridToCurvilinear(gridID1, NeedCorners::Yes);

    if (gridsize != gridInqSize(gridID1)) cdo_abort("Internal problem: gridsize changed!");

    if (varV.gridType != GRID_CURVILINEAR) gridID2 = gridToCurvilinear(gridID2, NeedCorners::Yes);

    if (gridsize != gridInqSize(gridID2)) cdo_abort("Internal problem: gridsize changed!");

    nlon = gridInqXsize(gridID1);
    nlat = gridInqYsize(gridID1);

    Varray<double> grid1x(gridsize), grid1y(gridsize);
    Varray<double> grid2x(gridsize), grid2y(gridsize);
    grid3x.resize(gridsize);
    grid3y.resize(gridsize);

    gridInqXvals(gridID1, grid1x.data());
    gridInqYvals(gridID1, grid1y.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID1, CDI_XAXIS, grid1x, "grid1 center lon");
    cdo_grid_to_degree(gridID1, CDI_YAXIS, grid1y, "grid1 center lat");

    gridInqXvals(gridID2, grid2x.data());
    gridInqYvals(gridID2, grid2y.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID2, CDI_XAXIS, grid2x, "grid2 center lon");
    cdo_grid_to_degree(gridID2, CDI_YAXIS, grid2y, "grid2 center lat");

    if (gpint) { uv_to_p_grid(nlon, nlat, grid1x, grid1y, grid2x, grid2y, grid3x, grid3y); }
    else
    {
      grid3x = grid1x;
      grid3y = grid1y;
    }

    auto gridID3 = gridCreate(GRID_CURVILINEAR, gridsize);
    cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID3);
    gridDefXsize(gridID3, nlon);
    gridDefYsize(gridID3, nlat);
    gridDefXvals(gridID3, grid3x.data());
    gridDefYvals(gridID3, grid3y.data());

    for (size_t i = 0; i < gridsize; ++i)
    {
      grid3x[i] *= DEG2RAD;
      grid3y[i] *= DEG2RAD;
    }

    auto vlistID3 = vlistDuplicate(vlistID1);
    vlistCat(vlistID3, vlistID2);

    if (varU.code == varV.code) vlistDefVarCode(vlistID3, 1, varU.code + 1);

    vlistChangeGrid(vlistID3, gridID1, gridID3);
    vlistChangeGrid(vlistID3, gridID2, gridID3);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID3, taxisID3);

    if (Options::cdoVerbose) vlistPrint(vlistID3);

    streamID3 = cdo_open_write(2);
    cdo_def_vlist(streamID3, vlistID3);

    missval1 = varU.missval;
    missval2 = varV.missval;
  }

  void
  run() override
  {
    Varray<double> urfield(gridsize);
    Varray<double> vrfield(gridsize);
    Varray<double> ufield_v(gridsize);
    Varray<double> vfield_v(gridsize);

    MatrixView<double> ufield(ufield_v.data(), nlat, nlon);
    MatrixView<double> vfield(vfield_v.data(), nlat, nlon);
    Varray<double> uhelp_v, vhelp_v;
    if (gpint)
    {
      auto gridsizex = (nlon + 2) * nlat;
      uhelp_v.resize(gridsizex);
      vhelp_v.resize(gridsizex);
    }
    MatrixView<double> uhelp(uhelp_v.data(), nlat, nlon + 2);
    MatrixView<double> vhelp(vhelp_v.data(), nlat, nlon + 2);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID3, tsID);

      auto numFields2 = cdo_stream_inq_timestep(streamID2, tsID);

      if (numFields != numFields2) cdo_abort("Input streams have different number of levels!");

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        (void) cdo_inq_field(streamID1);
        auto [varID2, levelID] = cdo_inq_field(streamID2);

        size_t numMissVals1, numMissVals2;
        cdo_read_field(streamID1, ufield_v.data(), &numMissVals1);
        cdo_read_field(streamID2, vfield_v.data(), &numMissVals2);

        // remove missing values
        if (numMissVals1 || numMissVals2)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(ufield_v[i], missval1)) ufield_v[i] = 0.0;
            if (fp_is_equal(vfield_v[i], missval2)) vfield_v[i] = 0.0;
          }
        }

        if (gpint)
        {
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

          // interpolate on pressure points
          for (size_t j = 1; j < nlat; ++j)
            for (size_t i = 0; i < nlon; ++i)
            {
              ufield[j][i] = (uhelp[j][i] + uhelp[j][i + 1]) * 0.5;
              vfield[j][i] = (vhelp[j - 1][i + 1] + vhelp[j][i + 1]) * 0.5;
            }
        }

        for (size_t i = 0; i < nlon; ++i)
        {
          ufield[0][i] = 0.0;
          vfield[0][i] = 0.0;
        }

        // rotate
        rotate_uv2(ufield_v, vfield_v, nlon, nlat, grid3x, grid3y, urfield, vrfield);

        // calc lat, lon, Auv and alpha
        /*
        {
        double lat, lon, auv, alpha;
        for ( int j = 1; j < nlat-1; j += 3 )
          for ( int i = 0; i < nlon; i += 3 )
            {
              lat = grid3y[j][i]*RAD2DEG;
              lon = grid3x[j][i]*RAD2DEG;
              auv = std::sqrt(urfield[j][i]*urfield[j][i] + vrfield[j][i]*vrfield[j][i]);
              alpha = std::atan2(vrfield[j][i], urfield[j][i]);
              alpha = 90. - alpha*RAD2DEG;

              if ( alpha <   0 ) alpha += 360.;
              if ( alpha > 360 ) alpha -= 360.;

              printf("%g %g %g %g\n", lon, lat, alpha, auv);
            }
        }
        */
        cdo_def_field(streamID3, 0, levelID);
        cdo_write_field(streamID3, urfield.data(), 0);
        cdo_def_field(streamID3, 1, levelID);
        cdo_write_field(streamID3, vrfield.data(), 0);
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
