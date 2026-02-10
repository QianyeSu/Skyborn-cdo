/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include <cdi.h>

#include <utility>

#include "process_int.h"
#include "param_conversion.h"
#include "util_string.h"
#include "pmlist.h"
#include <mpim_grid.h>

constexpr double UndefValue = DBL_MAX;

namespace
{
struct Parameter
{
  long east = 0;
  long west = 0;
  long south = 0;
  long north = 0;
  double value = UndefValue;
};

struct HaloInfo
{
  Parameter param;
  int gridID1 = -1;
  int gridID2 = -1;
};
}  // namespace

static HaloInfo
gen_tripolar_grid(int gridID1)
{
  auto nlon1 = gridInqXsize(gridID1);
  auto nlat1 = gridInqYsize(gridID1);

  auto nlon2 = nlon1;
  auto nlat2 = nlat1 + 2;

  auto gridtype = gridInqType(gridID1);

  auto gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  // auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  // auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);

  if (gridHasCoordinates(gridID1))
  {
    if (gridtype == GRID_CURVILINEAR)
    {
      Varray<double> xvals1(nlon1 * nlat1), yvals1(nlon1 * nlat1);
      Varray<double> xvals2(nlon2 * nlat2), yvals2(nlon2 * nlat2);

      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      for (size_t ilat = 0; ilat < nlat1; ilat++)
      {
        for (size_t ilon = 0; ilon < nlon1; ilon++)
        {
          xvals2[(ilat + 2) * nlon1 + ilon] = xvals1[ilat * nlon1 + ilon];
          yvals2[(ilat + 2) * nlon1 + ilon] = yvals1[ilat * nlon1 + ilon];
        }
      }

      for (size_t ilon = 0; ilon < nlon1; ilon++)
      {
        auto ilonr = nlon1 - ilon - 1;
        xvals2[1 * nlon1 + ilon] = xvals2[2 * nlon1 + ilonr];  // syncronise line 2 with line 3
        xvals2[0 * nlon1 + ilon] = xvals2[3 * nlon1 + ilonr];  // syncronise line 1 with line 4
        yvals2[1 * nlon1 + ilon] = yvals2[2 * nlon1 + ilonr];  // syncronise line 2 with line 3
        yvals2[0 * nlon1 + ilon] = yvals2[3 * nlon1 + ilonr];  // syncronise line 1 with line 4
      }

      gridDefXvals(gridID2, xvals2.data());
      gridDefYvals(gridID2, yvals2.data());
    }
  }

  if (gridHasBounds(gridID1))
  {
    if (gridtype == GRID_CURVILINEAR)
    {
      Varray<double> xbounds1(4 * nlon1 * nlat1), ybounds1(4 * nlon1 * nlat1);
      Varray<double> xbounds2(4 * nlon2 * nlat2), ybounds2(4 * nlon2 * nlat2);

      gridInqXbounds(gridID1, xbounds1.data());
      gridInqYbounds(gridID1, ybounds1.data());

      if (gridtype == GRID_CURVILINEAR)
      {
        gridDefNvertex(gridID2, 4);

        auto nlon4 = 4 * nlon1;
        for (size_t ilat = 0; ilat < nlat1; ilat++)
        {
          for (size_t ilon = 0; ilon < nlon4; ilon++)
          {
            xbounds2[(ilat + 2) * nlon4 + ilon] = xbounds1[ilat * nlon4 + ilon];
            ybounds2[(ilat + 2) * nlon4 + ilon] = ybounds1[ilat * nlon4 + ilon];
          }
        }

        for (size_t ilon = 0; ilon < nlon1; ilon++)
        {
          auto ilonr = nlon1 - ilon - 1;
          for (size_t k = 0; k < 4; ++k)
          {
            auto kr = 3 - k;
            xbounds2[1 * nlon4 + 4 * ilon + k] = xbounds2[2 * nlon4 + 4 * ilonr + kr];
            xbounds2[0 * nlon4 + 4 * ilon + k] = xbounds2[3 * nlon4 + 4 * ilonr + kr];
            ybounds2[1 * nlon4 + 4 * ilon + k] = ybounds2[2 * nlon4 + 4 * ilonr + kr];
            ybounds2[0 * nlon4 + 4 * ilon + k] = ybounds2[3 * nlon4 + 4 * ilonr + kr];
          }
        }
        /*
        for (size_t ilon = 0; ilon < 4*nlon1; ilon++)
          {
            auto ilonr = nlon4 - ilon - 1;
            xbounds2[1*nlon4 + ilon] = xbounds2[2*nlon4 + ilonr]; xbounds2[0*nlon4 + ilon] = xbounds2[3*nlon4 + ilonr];
            ybounds2[1*nlon4 + ilon] = ybounds2[2*nlon4 + ilonr]; ybounds2[0*nlon4 + ilon] = ybounds2[3*nlon4 + ilonr];
          }
        */
      }

      gridDefXbounds(gridID2, xbounds2.data());
      gridDefYbounds(gridID2, ybounds2.data());
    }
  }

  HaloInfo haloInfo;
  haloInfo.gridID1 = gridID1;
  haloInfo.gridID2 = gridID2;

  return haloInfo;
}

static int
gen_regular_grid(int gridID1, Parameter &haloParam)
{
  double cpi2 = M_PI * 2;

  auto isCircularGrid = gridIsCircular(gridID1);

  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);

  long nlon2 = nlon1 + haloParam.east + haloParam.west;
  long nlat2 = nlat1 + haloParam.south + haloParam.north;

  long lonMinIdx = 0;
  long lonMaxIdx = nlon1;
  if (haloParam.east < 0) lonMinIdx = -haloParam.east;
  if (haloParam.west < 0) lonMaxIdx += haloParam.west;
  long latMinIdx = 0;
  long latMaxIdx = nlat1;
  if (haloParam.south < 0) latMinIdx = -haloParam.south;
  if (haloParam.north < 0) latMaxIdx += haloParam.north;

  if (Options::cdoVerbose)
    printf("nlon1=%ld/nlat1=%ld nlon2=%ld/nlat2=%ld east=%ld/west=%ld/south=%ld/north=%ld\n", nlon1, nlat1, nlon2, nlat2,
           haloParam.east, haloParam.west, haloParam.south, haloParam.north);

  auto gridtype = gridInqType(gridID1);
  auto isCurvilinearGrid = (gridtype == GRID_CURVILINEAR);

  auto xinc = isCurvilinearGrid ? 0.0 : gridInqXinc(gridID1);
  auto yinc = isCurvilinearGrid ? 0.0 : gridInqYinc(gridID1);

  auto gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  cdiCopyKey(gridID1, CDI_GLOBAL, CDI_KEY_DATATYPE, gridID2);

  grid_copy_names(gridID1, gridID2);

  auto xunits = cdo::inq_key_string(gridID1, CDI_XAXIS, CDI_KEY_UNITS);
  // auto yunits = cdo::inq_key_string(gridID1, CDI_YAXIS, CDI_KEY_UNITS);

  if (xunits.rfind("degree", 0) == 0) cpi2 *= RAD2DEG;

  if (gridHasCoordinates(gridID1))
  {
    Varray<double> xvals1(isCurvilinearGrid ? nlon1 * nlat1 : nlon1);
    Varray<double> yvals1(isCurvilinearGrid ? nlon1 * nlat1 : nlat1);
    Varray<double> xvals2(isCurvilinearGrid ? nlon2 * nlat2 : nlon2);
    Varray<double> yvals2(isCurvilinearGrid ? nlon2 * nlat2 : nlat2);

    auto pxvals2 = xvals2.data();
    auto pyvals2 = yvals2.data();

    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    if (isCurvilinearGrid)
    {
      if (yvals1[0] > yvals1[nlon1 * nlat1 - 1]) std::swap(haloParam.south, haloParam.north);

      if (haloParam.south > 0)
      {
        pxvals2 += nlon2 * haloParam.south;
        pyvals2 += nlon2 * haloParam.south;
      }

      for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++)
      {
        const auto pxvals1 = &xvals1[ilat * nlon1];
        const auto pyvals1 = &yvals1[ilat * nlon1];
        if (isCircularGrid)
        {
          for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *pxvals2++ = pxvals1[ilon];
          for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *pyvals2++ = pyvals1[ilon];
        }
        else
        {
          for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *pxvals2++ = pxvals1[0];
          for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *pyvals2++ = pyvals1[0];
        }

        for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *pxvals2++ = pxvals1[ilon];
        for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *pyvals2++ = pyvals1[ilon];

        if (isCircularGrid)
        {
          for (long ilon = 0; ilon < haloParam.west; ilon++) *pxvals2++ = pxvals1[ilon];
          for (long ilon = 0; ilon < haloParam.west; ilon++) *pyvals2++ = pyvals1[ilon];
        }
        else
        {
          for (long ilon = 0; ilon < haloParam.west; ilon++) *pxvals2++ = pxvals1[nlon1 - 1];
          for (long ilon = 0; ilon < haloParam.west; ilon++) *pyvals2++ = pyvals1[nlon1 - 1];
        }
      }
      for (long ilat = 0; ilat < haloParam.north; ilat++)
      {
        long offset = nlon2 * (nlat1 + haloParam.south - 1);
        for (long ilon = 0; ilon < nlon2; ilon++) *pxvals2++ = xvals2[offset + ilon];
        for (long ilon = 0; ilon < nlon2; ilon++) *pyvals2++ = yvals2[offset + ilon];
      }
      for (long ilat = 0; ilat < haloParam.south; ilat++)
      {
        long offset = nlon2 * haloParam.south;
        for (long ilon = 0; ilon < nlon2; ilon++) xvals2[nlon2 * ilat + ilon] = xvals2[offset + ilon];
        for (long ilon = 0; ilon < nlon2; ilon++) yvals2[nlon2 * ilat + ilon] = yvals2[offset + ilon];
      }
    }
    else
    {
      if (yvals1[0] > yvals1[nlat1 - 1]) std::swap(haloParam.south, haloParam.north);

      if (isCircularGrid)
      {
        // clang-format off
        for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *pxvals2++ = xvals1[ilon] - cpi2;
        for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)         *pxvals2++ = xvals1[ilon];
        for (long ilon = 0; ilon < haloParam.west; ilon++)             *pxvals2++ = xvals1[ilon] + cpi2;
        // clang-format on
      }
      else
      {
        // clang-format off
        for (long ilon = haloParam.east; ilon > 0; ilon--)     *pxvals2++ = xvals1[0] - ilon * xinc;
        for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *pxvals2++ = xvals1[ilon];
        for (long ilon = 1; ilon <= haloParam.west; ilon++)    *pxvals2++ = xvals1[nlon1 - 1] + ilon * xinc;
        // clang-format on
      }

      // clang-format off
      for (long ilat = haloParam.south; ilat > 0; ilat--)    *pyvals2++ = yvals1[0] - ilat * yinc;
      for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++) *pyvals2++ = yvals1[ilat];
      for (long ilat = 1; ilat <= haloParam.north; ilat++)   *pyvals2++ = yvals1[nlat1 - 1] + ilat * yinc;
      // clang-format on
    }

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }

  if (isCircularGrid && gridHasBounds(gridID1))
  {
    auto nv = isCurvilinearGrid ? 4 : 2;
    gridDefNvertex(gridID2, nv);

    Varray<double> xbounds1(isCurvilinearGrid ? nv * nlon1 * nlat1 : nv * nlon1);
    Varray<double> ybounds1(isCurvilinearGrid ? nv * nlon1 * nlat1 : nv * nlat1);
    Varray<double> xbounds2(isCurvilinearGrid ? nv * nlon2 * nlat2 : nv * nlon2);
    Varray<double> ybounds2(isCurvilinearGrid ? nv * nlon2 * nlat2 : nv * nlat2);

    double *pxbounds2 = xbounds2.data();
    double *pybounds2 = ybounds2.data();

    gridInqXbounds(gridID1, xbounds1.data());
    gridInqYbounds(gridID1, ybounds1.data());

    if (isCurvilinearGrid)
    {
      if (isCircularGrid)
      {
        for (long ilat = 0; ilat < nlat1; ilat++)
        {
          const auto pxbounds1 = &xbounds1[nv * ilat * nlon1];
          const auto pybounds1 = &ybounds1[nv * ilat * nlon1];
          for (long ilon = nv * (nlon1 - haloParam.east); ilon < nv * nlon1; ilon++) *pxbounds2++ = pxbounds1[ilon];
          for (long ilon = nv * (nlon1 - haloParam.east); ilon < nv * nlon1; ilon++) *pybounds2++ = pybounds1[ilon];

          for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++) *pxbounds2++ = pxbounds1[ilon];
          for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++) *pybounds2++ = pybounds1[ilon];

          for (long ilon = 0; ilon < nv * haloParam.west; ilon++) *pxbounds2++ = pxbounds1[ilon];
          for (long ilon = 0; ilon < nv * haloParam.west; ilon++) *pybounds2++ = pybounds1[ilon];
        }
      }
    }
    else
    {
      if (isCircularGrid)
      {
        // clang-format off
        for (long ilon = nv * (nlon1 - haloParam.east); ilon < nv * nlon1; ilon++) *pxbounds2++ = xbounds1[ilon] - cpi2;
        for (long ilon = nv * lonMinIdx; ilon < nv * lonMaxIdx; ilon++)           *pxbounds2++ = xbounds1[ilon];
        for (long ilon = 0; ilon < nv * haloParam.west; ilon++)                    *pxbounds2++ = xbounds1[ilon] + cpi2;
        // clang-format on
      }

      for (long ilat = 0; ilat < nv * nlat2; ++ilat) ybounds2[ilat] = ybounds1[ilat];
    }

    gridDefXbounds(gridID2, xbounds2.data());
    gridDefYbounds(gridID2, ybounds2.data());
  }

  return gridID2;
}

static Parameter
get_parameter(void)
{
  Parameter haloParam;

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
      if      (key == "east")   haloParam.east = parameter_to_long(value);
      else if (key == "west")   haloParam.west = parameter_to_long(value);
      else if (key == "south")  haloParam.south = parameter_to_long(value);
      else if (key == "north")  haloParam.north = parameter_to_long(value);
      else if (key == "value")  haloParam.value = parameter_to_double(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return haloParam;
}

static Parameter
get_params()
{
  Parameter haloParam;

  if (cdo_operator_argc() == 2 && string_is_int(cdo_operator_argv(0)) && string_is_int(cdo_operator_argv(1)))
  {
    haloParam.east = parameter_to_long(cdo_operator_argv(0));
    haloParam.west = parameter_to_long(cdo_operator_argv(1));
  }
  else { haloParam = get_parameter(); }

  return haloParam;
}

static HaloInfo
gen_index_grid(int gridID1, Parameter &haloParam)
{
  auto isCircularGrid = gridIsCircular(gridID1);
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);

  if (isCircularGrid && haloParam.east > nlon) cdo_abort("east halo out of range (max=%ld).", nlon);
  if (isCircularGrid && haloParam.west > nlon) cdo_abort("west halo out of range (max=%ld).", nlon);
  if (haloParam.east < 0 && -haloParam.east > nlon - 1) cdo_abort("negative east halo out of range (max=%ld).", -nlon + 1);
  if (haloParam.west < 0 && -haloParam.west > nlon - 1) cdo_abort("negative west halo out of range (max=%ld).", -nlon + 1);
  if (haloParam.south < 0 && -haloParam.south > nlat - 1) cdo_abort("negative south halo out of range (max=%ld).", -nlat + 1);
  if (haloParam.north < 0 && -haloParam.north > nlat - 1) cdo_abort("negative north halo out of range (max=%ld).", -nlat + 1);
  if (haloParam.east + haloParam.west < -nlon + 1)
    cdo_abort("sum of negative east and west halo out of range (max=%ld).", -nlon + 1);
  if (haloParam.south + haloParam.north < -nlat + 1)
    cdo_abort("sum of negative south and north halo out of range (max=%ld).", -nlat + 1);

  auto gridID2 = gen_regular_grid(gridID1, haloParam);

  HaloInfo haloInfo;
  haloInfo.param = haloParam;
  haloInfo.gridID1 = gridID1;
  haloInfo.gridID2 = gridID2;

  return haloInfo;
}

static bool
regular_halo(Varray<double> const &array1, int gridID1, Varray<double> &array2, const Parameter &haloParam, double missval)
{
  auto recalcNumMiss = false;

  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);
  long nlon2 = nlon1 + haloParam.east + haloParam.west;

  auto isCircularGrid = gridIsCircular(gridID1);
  auto fillValueDefined = is_not_equal(haloParam.value, UndefValue);
  auto fillValue = fillValueDefined ? haloParam.value : missval;
  auto useFillValue = (!isCircularGrid || fillValueDefined);
  if (!fillValueDefined && (haloParam.east > 0 || haloParam.west > 0 || haloParam.south > 0 || haloParam.north > 0))
    recalcNumMiss = true;

  long lonMinIdx = 0;
  long lonMaxIdx = nlon1;
  if (haloParam.east < 0) lonMinIdx = -haloParam.east;
  if (haloParam.west < 0) lonMaxIdx += haloParam.west;
  long latMinIdx = 0;
  long latMaxIdx = nlat1;
  if (haloParam.south < 0) latMinIdx = -haloParam.south;
  if (haloParam.north < 0) latMaxIdx += haloParam.north;

  auto parray2 = array2.data();

  for (long i = 0; i < nlon2 * haloParam.south; i++) *parray2++ = fillValue;

  for (long ilat = latMinIdx; ilat < latMaxIdx; ilat++)
  {
    const auto parray1 = &array1[ilat * nlon1];
    if (useFillValue)
    {
      // clang-format off
      for (long ilon = 0; ilon < haloParam.east; ilon++)    *parray2++ = fillValue;
      for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++) *parray2++ = parray1[ilon];
      for (long ilon = 0; ilon < haloParam.west; ilon++)    *parray2++ = fillValue;
      // clang-format on
    }
    else
    {
      // clang-format off
      for (long ilon = nlon1 - haloParam.east; ilon < nlon1; ilon++) *parray2++ = parray1[ilon];
      for (long ilon = lonMinIdx; ilon < lonMaxIdx; ilon++)          *parray2++ = parray1[ilon];
      for (long ilon = 0; ilon < haloParam.west; ilon++)             *parray2++ = parray1[ilon];
      // clang-format on
    }
  }

  for (long i = 0; i < nlon2 * haloParam.north; i++) *parray2++ = fillValue;

  return recalcNumMiss;
}

static void
tripolar_halo(Varray<double> const &array1, int gridID1, Varray<double> &array2)
{
  auto nlon = gridInqXsize(gridID1);
  auto nlat = gridInqYsize(gridID1);

  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++) array2[(ilat + 2) * nlon + ilon] = array1[ilat * nlon + ilon];

  for (size_t ilon = 0; ilon < nlon; ilon++)
  {
    size_t ilonr = nlon - ilon - 1;
    array2[1 * nlon + ilon] = array2[2 * nlon + ilonr];  // syncronise line 2 with line 3
    array2[0 * nlon + ilon] = array2[3 * nlon + ilonr];  // syncronise line 1 with line 4
  }
}

static std::vector<HaloInfo>
get_haloInfoList(int vlistID1, int vlistID2, Parameter &haloParam, bool operSetHalo)
{
  std::vector<HaloInfo> haloInfoList;

  auto numGrids = vlistNumGrids(vlistID1);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID1 = vlistGrid(vlistID1, index);
    if (gridInqSize(gridID1) == 1) continue;

    auto gridtype = gridInqType(gridID1);

    auto isReg2dGeoGrid = (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR);
    auto isGenericGrid = (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0);
    if (isReg2dGeoGrid || isGenericGrid)
    {
      haloInfoList.push_back(operSetHalo ? gen_index_grid(gridID1, haloParam) : gen_tripolar_grid(gridID1));
      vlistChangeGridIndex(vlistID2, index, haloInfoList.back().gridID2);
    }
    else
    {
      if (gridInqSize(gridID1) > 2) cdo_warning("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
  }

  if (Options::cdoVerbose)
  {
    /*  for (auto const &sb : haloInfoList)
    if (sb.gridtype != GRID_UNSTRUCTURED && !is_healpix_grid(sb.gridID1))
      {
        cdo_print("box1 - idx1,idx2,idy1,idy2: %ld,%ld,%ld,%ld", sb.lon21 + 1, sb.lon22 + 1, sb.lat1 + 1, sb.lat2 + 1);
        cdo_print("box2 - idx1,idx2,idy1,idy2: %ld,%ld,%ld,%ld", sb.lon11 + 1, sb.lon12 + 1, sb.lat1 + 1, sb.lat2 + 1);
      }
      */
  }

  return haloInfoList;
}

static std::vector<bool>
get_processVars(VarList const &varList, std::vector<HaloInfo> const &haloInfoList)
{
  auto numVars = varList.numVars();

  std::vector<bool> processVars(numVars, false);

  int varID;
  for (auto const &haloInfo : haloInfoList)
  {
    for (varID = 0; varID < numVars; ++varID)
      if (haloInfo.gridID1 == varList.vars[varID].gridID) processVars[varID] = true;
  }

  for (varID = 0; varID < numVars; ++varID)
    if (processVars[varID]) break;

  if (varID >= numVars) cdo_abort("No processable variable found!");

  return processVars;
}

static const HaloInfo &
select_haloInfo(int gridID, std::vector<HaloInfo> const &haloInfoList)
{
  for (auto const &haloInfo : haloInfoList)
    if (gridID == haloInfo.gridID1) return haloInfo;

  cdo_abort("Internal problem, grid not found!");

  return haloInfoList[0];
}

class Sethalo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Sethalo",
    .operators = { { "sethalo", SethaloHelp }, { "tpnhalo", SethaloHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Sethalo> registration = RegisterEntry<Sethalo>(module);

private:
  int SETHALO{};
  int operatorID{};

  Parameter haloParam;
  std::vector<HaloInfo> haloInfoList;

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1;
  VarList varList2;

  std::vector<bool> processVars;

public:
  void
  init() override
  {
    SETHALO = module.get_id("sethalo");

    operatorID = cdo_operator_id();

    haloParam = get_params();

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    haloInfoList = get_haloInfoList(vlistID1, vlistID2, haloParam, (operatorID == SETHALO));

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    for (auto &var : varList1.vars) var.memType = MemType::Double;
    for (auto &var : varList2.vars) var.memType = MemType::Double;

    processVars = get_processVars(varList1, haloInfoList);

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
          auto recalcNumMiss = false;
          auto missval = var1.missval;

          field2.init(varList2.vars[varID]);

          auto const &haloInfo = select_haloInfo(var1.gridID, haloInfoList);

          if (operatorID == SETHALO)
            recalcNumMiss = regular_halo(field1.vec_d, haloInfo.gridID1, field2.vec_d, haloParam, missval);
          else
            tripolar_halo(field1.vec_d, haloInfo.gridID1, field2.vec_d);

          if (field1.numMissVals || recalcNumMiss) field_num_mv(field2);

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
  }
};
