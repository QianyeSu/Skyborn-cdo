/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Selbox     sellonlatbox    Select lon/lat box
      Selbox     selindexbox     Select index box
*/

#include <cdi.h>

#include <mpim_grid.h>
#include "grid_define.h"
#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include "selboxinfo.h"

static void
correct_xvals(long nlon, long inc, double *xvals)
{
  if (nlon > 1 && is_equal(xvals[0], xvals[(nlon - 1) * inc])) xvals[(nlon - 1) * inc] += 360;

  if (xvals[0] > xvals[(nlon - 1) * inc])
    for (long i = 0; i < nlon; ++i)
      if (xvals[i * inc] >= 180) xvals[i * inc] -= 360;

  for (long i = 0; i < nlon; ++i)
  {
    if (xvals[i * inc] < -180) xvals[i * inc] += 360;
    if (xvals[i * inc] > 360) xvals[i * inc] -= 360;
  }

  if (xvals[0] > xvals[(nlon - 1) * inc])
    for (long i = 1; i < nlon; ++i)
      if (xvals[i * inc] < xvals[(i - 1) * inc]) xvals[i * inc] += 360;
}

static void
gengridxyvals(int gridtype, int gridID1, int gridID2, long nlon, long nlat, long nlon2, long nlat2, SelboxInfo const &selboxInfo,
              bool unitsIsDegree)
{
  Varray<double> xvals1, yvals1, xvals2, yvals2;

  auto lxvals = (gridInqXvals(gridID1, nullptr) > 0);
  auto lyvals = (gridInqYvals(gridID1, nullptr) > 0);

  if (gridtype == GRID_CURVILINEAR)
  {
    if (lxvals && lyvals)
    {
      xvals1.resize(nlon * nlat);
      yvals1.resize(nlon * nlat);
      xvals2.resize(nlon2 * nlat2);
      yvals2.resize(nlon2 * nlat2);
    }
  }
  else
  {
    if (lxvals) xvals1.resize(nlon);
    if (lyvals) yvals1.resize(nlat);
    if (lxvals) xvals2.resize(nlon2);
    if (lyvals) yvals2.resize(nlat2);
  }

  auto pxvals2 = xvals2.data();
  auto pyvals2 = yvals2.data();

  if (lxvals) gridInqXvals(gridID1, xvals1.data());
  if (lyvals) gridInqYvals(gridID1, yvals1.data());

  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;

  if (gridtype == GRID_CURVILINEAR)
  {
    if (lxvals && lyvals)
      for (long ilat = lat1; ilat <= lat2; ilat++)
      {
        for (long ilon = lon21; ilon <= lon22; ilon++)
        {
          *pxvals2++ = xvals1[ilat * nlon + ilon];
          *pyvals2++ = yvals1[ilat * nlon + ilon];
        }
        for (long ilon = lon11; ilon <= lon12; ilon++)
        {
          *pxvals2++ = xvals1[ilat * nlon + ilon];
          *pyvals2++ = yvals1[ilat * nlon + ilon];
        }
      }
  }
  else
  {
    if (lxvals)
    {
      for (long i = lon21; i <= lon22; ++i) *pxvals2++ = xvals1[i];
      for (long i = lon11; i <= lon12; ++i) *pxvals2++ = xvals1[i];
      if (unitsIsDegree) correct_xvals(nlon2, 1, xvals2.data());
    }

    if (lyvals)
      for (long i = lat1; i <= lat2; ++i) *pyvals2++ = yvals1[i];
  }
  // for ( int i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
  // for ( int i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
  if (lxvals) gridDefXvals(gridID2, xvals2.data());
  if (lyvals) gridDefYvals(gridID2, yvals2.data());
}

static void
gengridXboundsCurvi(int gridID1, long nlon, long nlat, long nlon2, long nlat2, SelboxInfo const &selboxInfo,
                    Varray<double> &xbounds1, Varray<double> &xbounds2)
{
  xbounds1.resize(4 * nlon * nlat);
  xbounds2.resize(4 * nlon2 * nlat2);

  auto pxbounds2 = xbounds2.data();

  gridInqXbounds(gridID1, xbounds1.data());

  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  for (long ilat = lat1; ilat <= lat2; ilat++)
  {
    for (long ilon = 4 * lon21; ilon < 4 * (lon22 + 1); ilon++) *pxbounds2++ = xbounds1[4 * ilat * nlon + ilon];
    for (long ilon = 4 * lon11; ilon < 4 * (lon12 + 1); ilon++) *pxbounds2++ = xbounds1[4 * ilat * nlon + ilon];
  }
}

static void
gengridYboundsCurvi(int gridID1, long nlon, long nlat, long nlon2, long nlat2, SelboxInfo const &selboxInfo,
                    Varray<double> &ybounds1, Varray<double> &ybounds2)
{
  ybounds1.resize(4 * nlon * nlat);
  ybounds2.resize(4 * nlon2 * nlat2);

  auto pybounds2 = ybounds2.data();

  gridInqYbounds(gridID1, ybounds1.data());

  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  for (long ilat = lat1; ilat <= lat2; ilat++)
  {
    for (long ilon = 4 * lon21; ilon < 4 * (lon22 + 1); ilon++) *pybounds2++ = ybounds1[4 * ilat * nlon + ilon];
    for (long ilon = 4 * lon11; ilon < 4 * (lon12 + 1); ilon++) *pybounds2++ = ybounds1[4 * ilat * nlon + ilon];
  }
}

static void
gengridXboundsRect2D(int gridID1, long nlon, long nlon2, SelboxInfo const &selboxInfo, bool unitsIsDegree, Varray<double> &xbounds1,
                     Varray<double> &xbounds2)
{
  xbounds1.resize(2 * nlon);
  xbounds2.resize(2 * nlon2);

  auto pxbounds2 = xbounds2.data();

  gridInqXbounds(gridID1, xbounds1.data());

  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  for (long i = 2 * lon21; i < 2 * (lon22 + 1); ++i) *pxbounds2++ = xbounds1[i];
  for (long i = 2 * lon11; i < 2 * (lon12 + 1); ++i) *pxbounds2++ = xbounds1[i];

  if (unitsIsDegree)
  {
    correct_xvals(nlon2, 2, xbounds2.data());
    correct_xvals(nlon2, 2, xbounds2.data() + 1);
    // make sure that the first bound is less than the second
    long nx = 0;
    for (long i = 0; i < nlon2; ++i)
      if (xbounds2[2 * i] > xbounds2[2 * i + 1]) nx++;
    if (nx == nlon2 && xbounds2[0] > -180.)
      for (long i = 0; i < nlon2; ++i) xbounds2[2 * i] -= 360.;
  }
}

static void
gengridYboundsRect2D(int gridID1, long nlat, long nlat2, SelboxInfo const &selboxInfo, Varray<double> &ybounds1,
                     Varray<double> &ybounds2)
{
  ybounds1.resize(2 * nlat);
  ybounds2.resize(2 * nlat2);

  auto pybounds2 = ybounds2.data();

  gridInqYbounds(gridID1, ybounds1.data());

  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  for (long i = 2 * lat1; i < 2 * (lat2 + 1); ++i) *pybounds2++ = ybounds1[i];
}

static void
gengridXbounds(int gridID1, int gridID2, long nlon, long nlat, long nlon2, long nlat2, SelboxInfo const &selboxInfo,
               bool unitsIsDegree)
{
  Varray<double> xbounds1, xbounds2;
  auto gridtype = gridInqType(gridID1);

  if (gridtype == GRID_CURVILINEAR)
  {
    gridDefNvertex(gridID2, 4);
    gengridXboundsCurvi(gridID1, nlon, nlat, nlon2, nlat2, selboxInfo, xbounds1, xbounds2);
  }
  else
  {
    gridDefNvertex(gridID2, 2);
    gengridXboundsRect2D(gridID1, nlon, nlon2, selboxInfo, unitsIsDegree, xbounds1, xbounds2);
  }

  gridDefXbounds(gridID2, xbounds2.data());
}

static void
gengridYbounds(int gridID1, int gridID2, long nlon, long nlat, long nlon2, long nlat2, SelboxInfo const &selboxInfo)
{
  Varray<double> ybounds1, ybounds2;
  auto gridtype = gridInqType(gridID1);

  if (gridtype == GRID_CURVILINEAR)
  {
    gridDefNvertex(gridID2, 4);
    gengridYboundsCurvi(gridID1, nlon, nlat, nlon2, nlat2, selboxInfo, ybounds1, ybounds2);
  }
  else
  {
    gridDefNvertex(gridID2, 2);
    gengridYboundsRect2D(gridID1, nlat, nlat2, selboxInfo, ybounds1, ybounds2);
  }

  gridDefYbounds(gridID2, ybounds2.data());
}

static int
gengrid(int gridID1, SelboxInfo const &selboxInfo)
{
  auto const &lat1 = selboxInfo.lat1;
  auto const &lat2 = selboxInfo.lat2;
  auto const &lon11 = selboxInfo.lon11;
  auto const &lon12 = selboxInfo.lon12;
  auto const &lon21 = selboxInfo.lon21;
  auto const &lon22 = selboxInfo.lon22;
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);

  long nlon21 = lon12 - lon11 + 1;
  long nlon22 = lon22 - lon21 + 1;
  long nlon2 = nlon21 + nlon22;
  long nlat2 = lat2 - lat1 + 1;

  if (Options::cdoVerbose) cdo_print("nlon1=%ld  nlat1=%ld", nlon, nlat);
  if (Options::cdoVerbose) cdo_print("nlon2=%ld  nlat2=%ld", nlon2, nlat2);

  auto gridtype = gridInqType(gridID1);

  auto gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  if (nlon > 0) gridDefXsize(gridID2, nlon2);
  if (nlat > 0) gridDefYsize(gridID2, nlat2);
  if (nlat > 0) gridDefNP(gridID2, gridInqNP(gridID1));

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  if (gridtype == GRID_PROJECTION) grid_copy_mapping(gridID1, gridID2);

  auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);

  auto unitsIsDegree = (xunits.rfind("degree", 0) == 0);

  gengridxyvals(gridtype, gridID1, gridID2, nlon, nlat, nlon2, nlat2, selboxInfo, unitsIsDegree);

  if (gridInqXbounds(gridID1, nullptr)) gengridXbounds(gridID1, gridID2, nlon, nlat, nlon2, nlat2, selboxInfo, unitsIsDegree);
  if (gridInqYbounds(gridID1, nullptr)) gengridYbounds(gridID1, gridID2, nlon, nlat, nlon2, nlat2, selboxInfo);

  if (gridtype == GRID_CURVILINEAR && gridHasArea(gridID1))
  {
    Varray<double> areaIn(nlon * nlat), areaOut(nlon2 * nlat2);
    gridInqArea(gridID1, areaIn.data());
    auto pareaOut = areaOut.data();

    for (long ilat = lat1; ilat <= lat2; ilat++)
    {
      for (long ilon = lon21; ilon < (lon22 + 1); ilon++) *pareaOut++ = areaIn[ilat * nlon + ilon];
      for (long ilon = lon11; ilon < (lon12 + 1); ilon++) *pareaOut++ = areaIn[ilat * nlon + ilon];
    }

    gridDefArea(gridID2, areaOut.data());
  }

  auto projID1 = gridInqProj(gridID1);
  if (projID1 != CDI_UNDEFID && gridInqType(projID1) == GRID_PROJECTION)
  {
    auto projID2 = gridCreate(GRID_PROJECTION, nlon2 * nlat2);
    gridDefXsize(projID2, nlon2);
    gridDefYsize(projID2, nlat2);

    grid_copy_names(projID1, projID2);
    grid_copy_mapping(projID1, projID2);

    gengridxyvals(GRID_PROJECTION, projID1, projID2, nlon, nlat, nlon2, nlat2, selboxInfo, false);

    gridDefProj(gridID2, projID2);
  }

  return gridID2;
}

static void
copy_array_index(size_t n, const double *restrict v1, Varray<double> &v2, std::vector<int64_t> const &indices)
{
#ifdef _OPENMP
#pragma omp parallel for if (n > 99999) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) v2[i] = v1[indices[i]];
}

static void
copy_bounds_index(size_t n, size_t nv, const double *restrict v1, Varray<double> &v2, std::vector<int64_t> const &indices)
{
#ifdef _OPENMP
#pragma omp parallel for if (n > 99999) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t k = 0; k < nv; ++k) v2[i * nv + k] = v1[indices[i] * nv + k];
  }
}

int
gengridcell(int gridID1, size_t gridsize2, std::vector<int64_t> const &cellIndices)
{
  auto gridtype = gridInqType(gridID1);
  auto gridsize1 = gridInqSize(gridID1);

  auto setDimName = (gridtype == GRID_CURVILINEAR);
  if (gridtype == GRID_CURVILINEAR) gridtype = GRID_UNSTRUCTURED;

  if (gridtype != GRID_UNSTRUCTURED && !is_healpix_grid(gridID1) && gridtype != GRID_GENERIC) return -1;

  auto gridID2 = gridCreate(gridtype, gridsize2);

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  if (setDimName) cdiDefKeyString(gridID2, CDI_XAXIS, CDI_KEY_DIMNAME, "ncells");

  if (gridHasCoordinates(gridID1))
  {
    const auto *xvals1 = gridInqXvalsPtr(gridID1);
    const auto *yvals1 = gridInqYvalsPtr(gridID1);

    Varray<double> xyvals2(gridsize2);

    copy_array_index(gridsize2, xvals1, xyvals2, cellIndices);
    gridDefXvals(gridID2, xyvals2.data());

    copy_array_index(gridsize2, yvals1, xyvals2, cellIndices);
    gridDefYvals(gridID2, xyvals2.data());
  }

  if (gridHasBounds(gridID1))
  {
    const auto *xbounds1 = gridInqXboundsPtr(gridID1);
    const auto *ybounds1 = gridInqYboundsPtr(gridID1);

    auto nv = gridInqNvertex(gridID1);
    gridDefNvertex(gridID2, nv);

    Varray<double> xybounds2(nv * gridsize2);

    copy_bounds_index(gridsize2, nv, xbounds1, xybounds2, cellIndices);
    gridDefXbounds(gridID2, xybounds2.data());

    copy_bounds_index(gridsize2, nv, ybounds1, xybounds2, cellIndices);
    gridDefYbounds(gridID2, xybounds2.data());
  }

  if (gridHasArea(gridID1))
  {
    Varray<double> areaIn(gridsize1), areaOut(gridsize2);
    gridInqArea(gridID1, areaIn.data());

    for (size_t i = 0; i < gridsize2; ++i) areaOut[i] = areaIn[cellIndices[i]];

    gridDefArea(gridID2, areaOut.data());
  }

  return gridID2;
}

SelboxInfo
gen_lonlat_selbox_reg2d(int gridID, double xlon1, double xlon2, double xlat1, double xlat2)
{
  SelboxInfo selboxInfo;
  auto &lat1 = selboxInfo.lat1;
  auto &lat2 = selboxInfo.lat2;
  auto &lon11 = selboxInfo.lon11;
  auto &lon12 = selboxInfo.lon12;
  auto &lon21 = selboxInfo.lon21;
  auto &lon22 = selboxInfo.lon22;
  lat1 = 0;
  lat2 = 0;
  lon11 = 1;
  lon12 = 0;
  lon21 = 0;
  lon22 = 0;

  if (is_not_equal(xlon1, xlon2))
  {
    xlon2 -= 360 * std::floor((xlon2 - xlon1) / 360);
    if (is_equal(xlon1, xlon2)) xlon2 += 360;
  }
  else { xlon2 += 0.00001; }

  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  Varray<double> xvals(nlon), yvals(nlat);

  if (nlon > 0)
  {
    gridInqXvals(gridID, xvals.data());
    cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");

    xlon2 -= 360 * std::floor((xlon1 - xvals[0]) / 360);
    xlon1 -= 360 * std::floor((xlon1 - xvals[0]) / 360);
  }

  if (nlat > 0)
  {
    gridInqYvals(gridID, yvals.data());
    cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");
  }

  if (nlon > 0)
  {
    // while ( nlon == 1 || (xvals[nlon-1] - xvals[0]) >= 360 ) nlon--;

    for (lon21 = 0; lon21 < nlon && xvals[lon21] < xlon1; lon21++);
    for (lon22 = lon21; lon22 < nlon && xvals[lon22] < xlon2; lon22++);

    if (lon22 >= nlon || xvals[lon22] > xlon2) lon22--;

    xlon1 -= 360;
    xlon2 -= 360;

    for (lon11 = 0; xvals[lon11] < xlon1; lon11++);
    for (lon12 = lon11; lon12 < nlon && xvals[lon12] < xlon2; lon12++);

    // lon12--;
    if (lon12 >= nlon || xvals[lon12] > xlon2) lon12--;
    if (lon21 < nlon && lon12 >= 0 && is_equal(xvals[lon12], xvals[lon21])) lon12--;

    if (lon12 - lon11 + 1 + lon22 - lon21 + 1 <= 0) cdo_abort("Longitudinal dimension is too small!");
  }

  if (nlat > 0)
  {
    if (yvals[0] > yvals[nlat - 1])
    {
      if (xlat1 > xlat2)
      {
        for (lat1 = 0; lat1 < nlat && yvals[lat1] > xlat1; lat1++);
        for (lat2 = nlat - 1; lat2 && yvals[lat2] < xlat2; lat2--);
      }
      else
      {
        for (lat1 = 0; lat1 < nlat && yvals[lat1] > xlat2; lat1++);
        for (lat2 = nlat - 1; lat2 && yvals[lat2] < xlat1; lat2--);
      }
    }
    else
    {
      if (xlat1 < xlat2)
      {
        for (lat1 = 0; lat1 < nlat && yvals[lat1] < xlat1; lat1++);
        for (lat2 = nlat - 1; lat2 && yvals[lat2] > xlat2; lat2--);
      }
      else
      {
        for (lat1 = 0; lat1 < nlat && yvals[lat1] < xlat2; lat1++);
        for (lat2 = nlat - 1; lat2 && yvals[lat2] > xlat1; lat2--);
      }
    }

    if (lat2 - lat1 + 1 <= 0) cdo_abort("Latitudinal dimension is too small!");
  }

  return selboxInfo;
}

static SelboxInfo
gen_lonlat_selbox_curv(int gridID, double xlon1, double xlon2, double xlat1, double xlat2)
{
  SelboxInfo selboxInfo;
  auto &lat1 = selboxInfo.lat1;
  auto &lat2 = selboxInfo.lat2;
  auto &lon11 = selboxInfo.lon11;
  auto &lon12 = selboxInfo.lon12;
  auto &lon21 = selboxInfo.lon21;
  auto &lon22 = selboxInfo.lon22;
  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);
  size_t gridsize = nlon * nlat;

  auto grid_is_circular = gridIsCircular(gridID);

  Varray<double> xvals(gridsize);
  Varray<double> yvals(gridsize);

  gridInqXvals(gridID, xvals.data());
  gridInqYvals(gridID, yvals.data());

  // Convert lat/lon units if required
  cdo_grid_to_degree(gridID, CDI_XAXIS, xvals, "grid center lon");
  cdo_grid_to_degree(gridID, CDI_YAXIS, yvals, "grid center lat");

  if (xlon1 > xlon2) cdo_abort("The second longitude have to be greater than the first one!");
  /*
  if ( xlon2 < 180. )
    {
      xlon1 += 360;
      xlon2 += 360;
    }
  */
  if (xlat1 > xlat2) std::swap(xlat1, xlat2);

  lat1 = nlat - 1;
  lat2 = 0;
  lon11 = 0;
  lon12 = -1;
  // lon11 = nlon-1;
  // lon12 = 0;
  lon21 = nlon - 1;
  lon22 = 0;

  bool lp2 = false;
  double xfirst, xlast, ylast;

  if (grid_is_circular)
  {
    for (long ilat = 0; ilat < nlat; ilat++)
    {
      xlast = xvals[ilat * nlon + nlon - 1];
      ylast = yvals[ilat * nlon + nlon - 1];
      if (ylast >= xlat1 && ylast <= xlat2)
        if (xlon1 <= xlast && xlon2 > xlast && (xlon2 - xlon1) < 360)
        {
          lon11 = nlon - 1;
          lon12 = 0;
          lp2 = true;
          break;
        }
    }
  }

  for (long ilat = 0; ilat < nlat; ilat++)
    for (long ilon = 0; ilon < nlon; ilon++)
    {
      auto xval = xvals[ilat * nlon + ilon];
      auto yval = yvals[ilat * nlon + ilon];
      if (yval >= xlat1 && yval <= xlat2)
      {
        if (lp2)
        {
          xfirst = xvals[ilat * nlon];
          if (xfirst < xlon1) xfirst = xlon1;

          xlast = xvals[ilat * nlon + nlon - 1];
          if (xlast > xlon2) xlast = xlon2;

          // printf("%g %g %g %g %g %g\n", yval, xval, xlon1, xlast,  xfirst, xlon2);
          if (xval >= xlon1 && xval <= xlast)
          {
            lon21 = std::min(lon21, ilon);
            lon22 = std::max(lon22, ilon);
            lat1 = std::min(lat1, ilat);
            lat2 = std::max(lat2, ilat);
          }
          else if (xval >= xfirst && xval <= xlon2)
          {
            lon11 = std::min(lon11, ilon);
            lon12 = std::max(lon12, ilon);
            lat1 = std::min(lat1, ilat);
            lat2 = std::max(lat2, ilat);
          }
        }
        else
        {
          if (((xval >= xlon1 && xval <= xlon2) || (xval - 360 >= xlon1 && xval - 360 <= xlon2)
               || (xval + 360 >= xlon1 && xval + 360 <= xlon2)))
          {
            lon21 = std::min(lon21, ilon);
            lon22 = std::max(lon22, ilon);
            lat1 = std::min(lat1, ilat);
            lat2 = std::max(lat2, ilat);
          }
        }
        /*
        if ( xval >= xlon1 && xval <= xlon2 )
          {
            lon21 = std::min(lon21, ilon);
            lon22 = std::max(lon22, ilon);
            lat1 = std::min(lat1, ilat);
            lat2 = std::max(lat2, ilat);
          }
        else if ( xval >= xlon1-360 && xval <= xlon2-360 )
          {
            lon11 = std::min(lon11, ilon);
            lon12 = std::max(lon12, ilon);
            lat1 = std::min(lat1, ilat);
            lat2 = std::max(lat2, ilat);
          }
        */
      }
    }

  // while ( lon12 >= lon21 ) lon12--;
  // if ( lon12 <= lon11 ) { lon11 = nlon-1; lon12 = 0; }
  if (lon12 == 0 && lon11 > 0)
  {
    lon11 = 0;
    lon12 = -1;
  }
  /*
  printf("lon21, lon22, lon11, lon12  idx: %ld %ld %ld %ld  lon: %g %g %g %g\n",
         lon21, lon22, lon11, lon12, xvals[lon21], xvals[lon22], xvals[lon11], xvals[lon12]);
  */
  if (lat2 - lat1 + 1 <= 0) cdo_abort("Latitudinal dimension is too small!");

  return selboxInfo;
}

void
getlonlatparams(int argcOffset, double &xlon1, double &xlon2, double &xlat1, double &xlat2)
{
  bool lset = false;

  auto nargc = cdo_operator_argc() - argcOffset;
  if (nargc == 1)
  {
    auto gridname = cdo_operator_argv(argcOffset + 0);
    if (gridname == "europe")
    {
      xlon1 = -30;
      xlon2 = 60;
      xlat1 = 30;
      xlat2 = 80;
      lset = true;
    }
  }

  if (!lset)
  {
    operator_check_argc(argcOffset + 4);

    xlon1 = parameter_to_double(cdo_operator_argv(argcOffset + 0));
    xlon2 = parameter_to_double(cdo_operator_argv(argcOffset + 1));
    xlat1 = parameter_to_double(cdo_operator_argv(argcOffset + 2));
    xlat2 = parameter_to_double(cdo_operator_argv(argcOffset + 3));
  }
}

SelboxInfo
gen_lonlat_selbox(int argcOffset, int gridID)
{
  double xlon1 = 0, xlon2 = 0, xlat1 = 0, xlat2 = 0;
  getlonlatparams(argcOffset, xlon1, xlon2, xlat1, xlat2);

  auto gridtype = gridInqType(gridID);

  if (gridtype == GRID_CURVILINEAR)
    return gen_lonlat_selbox_curv(gridID, xlon1, xlon2, xlat1, xlat2);
  else
    return gen_lonlat_selbox_reg2d(gridID, xlon1, xlon2, xlat1, xlat2);
}

static SelboxInfo
selbox_lonlat_grid(int gridID1)
{
  auto selboxInfo = gen_lonlat_selbox(0, gridID1);
  selboxInfo.gridID1 = gridID1;
  selboxInfo.gridtype = gridInqType(gridID1);

  selboxInfo.gridID2 = gengrid(gridID1, selboxInfo);
  return selboxInfo;
}

static SelboxInfo
selbox_cell_grid(int gridID1)
{
  SelboxInfo selboxInfo;
  selboxInfo.gridID1 = gridID1;
  selboxInfo.gridtype = gridInqType(gridID1);

  auto &gridsize2 = selboxInfo.nvals;
  auto &cellIndices = selboxInfo.cellIndices;

  int argcOffset = 0;
  operator_check_argc(argcOffset + 4);

  auto xlon1 = parameter_to_double(cdo_operator_argv(argcOffset + 0));
  auto xlon2 = parameter_to_double(cdo_operator_argv(argcOffset + 1));
  auto xlat1 = parameter_to_double(cdo_operator_argv(argcOffset + 2));
  auto xlat2 = parameter_to_double(cdo_operator_argv(argcOffset + 3));

  if (xlon1 >= xlon2) std::swap(xlon1, xlon2);
  if (xlat1 >= xlat2) std::swap(xlat1, xlat2);

  auto gridtype = gridInqType(gridID1);
  auto gridsize1 = gridInqSize(gridID1);

  if (gridtype != GRID_UNSTRUCTURED && !is_healpix_grid(gridID1)) cdo_abort("Internal problem, wrong grid type!");

  auto gridID0 = gridID1;
  gridID1 = generate_full_grid(gridID1);
  if (!gridHasCoordinates(gridID1)) cdo_abort("Cell center coordinates missing!");

  {
    Varray<double> xvals(gridsize1), yvals(gridsize1);
    gridInqXvals(gridID1, xvals.data());
    gridInqYvals(gridID1, yvals.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID1, CDI_XAXIS, xvals, "grid center lon");
    cdo_grid_to_degree(gridID1, CDI_YAXIS, yvals, "grid center lat");

    // find gridsize2
    long maxcell = 0;
    long nvals = 0;
    for (size_t i = 0; i < gridsize1; ++i)
    {
      auto xval = xvals[i];
      auto yval = yvals[i];
      if (yval >= xlat1 && yval <= xlat2)
        if ((xval >= xlon1 && xval <= xlon2) || (xval + 360 >= xlon1 && xval + 360 <= xlon2)
            || (xval - 360 >= xlon1 && xval - 360 <= xlon2))
        {
          nvals++;
          if (nvals > maxcell)
          {
            constexpr long cellinc = 4096;
            maxcell += cellinc;
            cellIndices.resize(maxcell);
          }
          cellIndices[nvals - 1] = i;
        }
    }

    if (nvals == 0) cdo_abort("No grid points found!");

    gridsize2 = nvals;
  }

  selboxInfo.gridID2 = gengridcell(gridID1, gridsize2, cellIndices);
  if (gridID0 != gridID1) gridDestroy(gridID1);

  return selboxInfo;
}

static void
checkBounds(const char *txt, long minVal, long maxVal, long &lval)
{
  if (lval < minVal)
  {
    cdo_warning("%s index out of range, set to %ld!", txt, minVal);
    lval = minVal;
  }

  if (lval > maxVal)
  {
    cdo_warning("%s index out of range, set to %ld!", txt, maxVal);
    lval = maxVal;
  }
}

SelboxInfo
gen_index_selbox(int argcOffset, int gridID)
{
  SelboxInfo selboxInfo;
  auto &lat1 = selboxInfo.lat1;
  auto &lat2 = selboxInfo.lat2;
  auto &lon11 = selboxInfo.lon11;
  auto &lon12 = selboxInfo.lon12;
  auto &lon21 = selboxInfo.lon21;
  auto &lon22 = selboxInfo.lon22;

  operator_check_argc(argcOffset + 4);

  lon11 = parameter_to_int(cdo_operator_argv(argcOffset + 0));
  lon12 = parameter_to_int(cdo_operator_argv(argcOffset + 1));
  lat1 = parameter_to_int(cdo_operator_argv(argcOffset + 2));
  lat2 = parameter_to_int(cdo_operator_argv(argcOffset + 3));

  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  if (lon11 < 0) lon11 = nlon + lon11 + 1;
  if (lon12 < 0) lon12 = nlon + lon12 + 1;
  if (lat1 < 0) lat1 = nlat + lat1 + 1;
  if (lat2 < 0) lat2 = nlat + lat2 + 1;

  if (lat1 > lat2) std::swap(lat1, lat2);

  checkBounds("First latitude", 1, nlat, lat1);
  checkBounds("Last latitude", 1, nlat, lat2);
  checkBounds("First longitude", 1, nlon, lon11);
  checkBounds("Last longitude", 1, nlon, lon12);

  lon11--;
  lon12--;
  lat1--;
  lat2--;

  if (lon11 > lon12)
  {
    lon21 = lon11;
    lon22 = nlon - 1;
    lon11 = 0;
  }
  else
  {
    if (lon12 > nlon - 1)
    {
      lon21 = lon11;
      lon22 = nlon - 1;
      lon11 = 0;
      lon12 = 0;
    }
    else
    {
      lon21 = 0;
      lon22 = -1;
    }
  }

  return selboxInfo;
}

static SelboxInfo
selbox_index_grid(int gridID)
{
  auto selboxInfo = gen_index_selbox(0, gridID);
  selboxInfo.gridID1 = gridID;
  selboxInfo.gridtype = gridInqType(gridID);

  if (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_LCC)
    selboxInfo.gridID2 = cdo_define_subgrid_grid(gridID, selboxInfo.lon11, selboxInfo.lon12, selboxInfo.lat1, selboxInfo.lat2);
  else
    selboxInfo.gridID2 = gengrid(gridID, selboxInfo);

  return selboxInfo;
}

template <typename T1, typename T2>
static void
window_box(int nwpv, const T1 *array1, int gridID, T2 *array2, long lat1, long lat2, long lon11, long lon12, long lon21, long lon22)
{
  long nlon = gridInqXsize(gridID);

  if (nwpv == 2)
  {
    for (long ilat = lat1; ilat <= lat2; ilat++)
    {
      for (long ilon = lon21; ilon <= lon22; ilon++)
      {
        *array2++ = array1[ilat * nlon * 2 + ilon * 2];
        *array2++ = array1[ilat * nlon * 2 + ilon * 2 + 1];
      }
      for (long ilon = lon11; ilon <= lon12; ilon++)
      {
        *array2++ = array1[ilat * nlon * 2 + ilon * 2];
        *array2++ = array1[ilat * nlon * 2 + ilon * 2 + 1];
      }
    }
  }
  else
  {
    for (long ilat = lat1; ilat <= lat2; ilat++)
    {
      for (long ilon = lon21; ilon <= lon22; ilon++) *array2++ = array1[ilat * nlon + ilon];
      for (long ilon = lon11; ilon <= lon12; ilon++) *array2++ = array1[ilat * nlon + ilon];
    }
  }
}

static void
window_box(Field const &field1, Field &field2, long lat1, long lat2, long lon11, long lon12, long lon21, long lon22)
{
  auto func = [&](auto const &v1, auto &v2)
  { window_box(field1.nwpv, v1.data(), field1.grid, v2.data(), lat1, lat2, lon11, lon12, lon21, lon22); };
  field_operation2(func, field1, field2);
}

template <typename T1, typename T2>
static void
window_cell(int nwpv, Varray<T1> const &v1, Varray<T2> &v2, long n, const std::vector<int64_t> &indices)
{
  if (nwpv == 2)
  {
    for (long i = 0; i < n; ++i)
    {
      v2[i * 2] = v1[indices[i] * 2];
      v2[i * 2 + 1] = v1[indices[i] * 2 + 1];
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > 99999) default(shared) schedule(static)
#endif
    for (long i = 0; i < n; ++i) v2[i] = v1[indices[i]];
  }
}

void
window_cell(Field const &field1, Field &field2, std::vector<int64_t> const &cellIndices)
{
  auto func = [&](auto &v1, auto &v2) { window_cell(field1.nwpv, v1, v2, field2.gridsize, cellIndices); };
  field_operation2(func, field1, field2);
}

static std::vector<SelboxInfo>
get_selboxInfoList(int vlistID1, int vlistID2, bool operIndexBox)
{
  std::vector<SelboxInfo> selboxInfoList;

  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    if (gridInqSize(gridID1) == 1) continue;

    auto gridtype = gridInqType(gridID1);
    auto projtype = gridInqProjType(gridID1);

    auto isReg2dGeoGrid = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR);
    auto projHasGeoCoords = (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_RLL);
    //                     || (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_LCC)
    //                     || (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_STERE);
    if (isReg2dGeoGrid || projHasGeoCoords || (operIndexBox && (gridtype == GRID_GENERIC || gridtype == GRID_PROJECTION))
        || (!operIndexBox && (gridtype == GRID_UNSTRUCTURED || is_healpix_grid(gridID1))))
    {
      if (operIndexBox) { selboxInfoList.push_back(selbox_index_grid(gridID1)); }
      else
      {
        selboxInfoList.push_back((gridtype == GRID_UNSTRUCTURED || is_healpix_grid(gridID1)) ? selbox_cell_grid(gridID1)
                                                                                             : selbox_lonlat_grid(gridID1));
      }

      vlistChangeGridIndex(vlistID2, index, selboxInfoList.back().gridID2);
    }
    else
    {
      if (gridInqSize(gridID1) > 2) cdo_warning("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
  }

  if (Options::cdoVerbose)
  {
    for (auto const &sb : selboxInfoList)
      if (sb.gridtype != GRID_UNSTRUCTURED && !is_healpix_grid(sb.gridID1))
      {
        cdo_print("box1 - idx1,idx2,idy1,idy2: %ld,%ld,%ld,%ld", sb.lon21 + 1, sb.lon22 + 1, sb.lat1 + 1, sb.lat2 + 1);
        cdo_print("box2 - idx1,idx2,idy1,idy2: %ld,%ld,%ld,%ld", sb.lon11 + 1, sb.lon12 + 1, sb.lat1 + 1, sb.lat2 + 1);
      }
  }

  return selboxInfoList;
}

static std::vector<bool>
get_processVars(VarList const &varList, std::vector<SelboxInfo> const &selboxInfoList)
{
  auto numVars = varList.numVars();

  std::vector<bool> processVars(numVars, false);

  int varID;
  for (auto const &sb : selboxInfoList)
  {
    for (varID = 0; varID < numVars; ++varID)
      if (sb.gridID1 == varList.vars[varID].gridID) processVars[varID] = true;
  }

  for (varID = 0; varID < numVars; ++varID)
    if (processVars[varID]) break;

  if (varID >= numVars) cdo_abort("No processable variable found!");

  return processVars;
}

static const SelboxInfo &
select_selboxInfo(int gridID, std::vector<SelboxInfo> const &selboxInfoList)
{
  for (auto const &selboxInfo : selboxInfoList)
    if (gridID == selboxInfo.gridID1) return selboxInfo;

  cdo_abort("Internal problem, grid not found!");

  return selboxInfoList[0];
}

class Selbox : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selbox",
    .operators = { { "sellonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude", SelboxHelp },
                   { "selindexbox", 0, 0, "index of first and last longitude and index of first and last latitude", SelboxHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Selbox> registration = RegisterEntry<Selbox>(module);

  int SELLONLATBOX, SELINDEXBOX;
  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;

  int vlistID2{ CDI_UNDEFID };
  int operatorID;

  std::vector<SelboxInfo> selboxInfoList;
  std::vector<bool> processVars;

public:
  void
  init() override
  {
    SELLONLATBOX = module.get_id("sellonlatbox");
    SELINDEXBOX = module.get_id("selindexbox");

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    selboxInfoList = get_selboxInfoList(vlistID1, vlistID2, (operatorID == SELINDEXBOX));

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);

    processVars = get_processVars(varList1, selboxInfoList);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
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

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var1 = varList1.vars[varID];
        field1.init(var1);
        cdo_read_field(streamID1, field1);

        cdo_def_field(streamID2, varID, levelID);

        if (processVars[varID])
        {
          field2.init(varList2.vars[varID]);

          auto const &sb = select_selboxInfo(var1.gridID, selboxInfoList);

          if (operatorID == SELLONLATBOX && (sb.gridtype == GRID_UNSTRUCTURED || is_healpix_grid(sb.gridID1)))
            window_cell(field1, field2, sb.cellIndices);
          else
            window_box(field1, field2, sb.lat1, sb.lat2, sb.lon11, sb.lon12, sb.lon21, sb.lon22);

          if (field1.numMissVals) field_num_mv(field2);

          cdo_write_field(streamID2, field2);
        }
        else { cdo_write_field(streamID2, field1); }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
