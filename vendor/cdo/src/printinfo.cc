#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "cdo_options.h"
#include "util_string.h"
#include "cdo_cdi_wrapper.h"
#include "cdi_uuid.h"
#include "mpmo_color.h"
#include "datetime.h"
#include "compare.h"
#include "printinfo.h"

#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

std::string
date_to_string(CdiDate date)
{
  char cstr[32];
  std::snprintf(cstr, sizeof(cstr), DATE_FORMAT, date.year, date.month, date.day);

  return std::string(cstr);
}

std::string
time_to_string(CdiTime time)
{
  static bool readEnv = true;
  static int msDigitsNum = 0;
  if (readEnv)
  {
    readEnv = false;
    auto envString = getenv_string("CDO_MS_DIGITS");
    if (envString.size())
    {
      int ival = std::stoi(envString);
      if (ival > 0) msDigitsNum = ival;
      if (ival > 3) msDigitsNum = 3;
    }
  }

  char cstr[32];
  if (msDigitsNum)
    std::snprintf(cstr, sizeof(cstr), "%2.2d:%2.2d:%0*.*f", time.hour, time.minute, msDigitsNum + 3, msDigitsNum,
                  time.second + time.ms / 1000.0);
  else
    std::snprintf(cstr, sizeof(cstr), TIME_FORMAT, time.hour, time.minute, time.second);

  return std::string(cstr);
}

std::string
datetime_to_string(CdiDateTime dt)
{
  return date_to_string(dt.date) + "T" + time_to_string(dt.time);
}

const char *
comptype_to_name(int compType)
{
  switch (compType)
  {
    case CDI_COMPRESS_SZIP: return "szip";
    case CDI_COMPRESS_AEC: return "aec";
    case CDI_COMPRESS_ZIP: return "zip";
    case CDI_COMPRESS_JPEG: return "jpeg";
    case CDI_COMPRESS_FILTER: return "filter";
  }
  return " ";
}

void
print_filetype(CdoStreamID streamID, int vlistID)
{
  auto filetype = cdo_inq_filetype(streamID);

  auto filetypestr = cdo::filetype_to_cstr(filetype);
  if (filetypestr == nullptr || *filetypestr == 0)
    printf("  unsupported filetype %d", filetype);
  else
    printf("%s", filetypestr);

  if (filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG)
  {
    switch (cdo_inq_byteorder(streamID))
    {
      case CDI_BIGENDIAN: printf("  BIGENDIAN"); break;
      case CDI_LITTLEENDIAN: printf("  LITTLEENDIAN"); break;
      default: printf("  byteorder: %d undefined", cdo_inq_byteorder(streamID)); break;
    }
  }

  int nvars = vlistNvars(vlistID);
  constexpr int comps[] = { CDI_COMPRESS_SZIP, CDI_COMPRESS_AEC, CDI_COMPRESS_ZIP, CDI_COMPRESS_JPEG, CDI_COMPRESS_FILTER };
  int kk = 0;
  for (size_t k = 0; k < sizeof(comps) / sizeof(int); ++k)
    for (int varID = 0; varID < nvars; ++varID)
    {
      auto comptype = vlistInqVarCompType(vlistID, varID);
      if (comptype == comps[k])
      {
        printf("%c%s", (kk++ == 0) ? ' ' : '/', comptype_to_name(comptype));
        break;
      }
    }

  printf("\n");
}

static void
print_xvals(int gridID, int dig)
{
  auto xsize = gridInqXsize(gridID);
  if (xsize && gridInqXvals(gridID, NULL))
  {
    auto xname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
    auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
    auto xfirst = gridInqXval(gridID, 0);
    auto xlast = gridInqXval(gridID, xsize - 1);
    auto xinc = gridInqXinc(gridID);
    auto gridType = gridInqType(gridID);

    if (gridType == GRID_GAUSSIAN_REDUCED)
    {
      if (xsize == 2) std::fprintf(stdout, "%33s : %.*g to %.*g [%s]\n", xname.c_str(), dig, xfirst, dig, xlast, xunits.c_str());
    }
    else
    {
      std::fprintf(stdout, "%33s : %.*g", xname.c_str(), dig, xfirst);
      if (xsize > 1)
      {
        std::fprintf(stdout, " to %.*g", dig, xlast);
        if (is_not_equal(xinc, 0.0)) std::fprintf(stdout, " by %.*g", dig, xinc);
      }
      if (xunits.size()) std::fprintf(stdout, " [%s]", xunits.c_str());
      if (gridIsCircular(gridID)) std::fprintf(stdout, "  circular");
      std::fprintf(stdout, "\n");
    }
  }
}

static void
print_yvals(int gridID, int dig)
{
  auto ysize = gridInqYsize(gridID);
  if (ysize && gridInqYvals(gridID, NULL))
  {
    auto yname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
    auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
    auto yfirst = gridInqYval(gridID, 0);
    auto ylast = gridInqYval(gridID, ysize - 1);
    auto yinc = gridInqYinc(gridID);
    std::fprintf(stdout, "%33s : %.*g", yname.c_str(), dig, yfirst);
    if (ysize > 1)
    {
      auto gridtype = gridInqType(gridID);
      std::fprintf(stdout, " to %.*g", dig, ylast);
      if (is_not_equal(yinc, 0.0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED)
        std::fprintf(stdout, " by %.*g", dig, yinc);
    }
    if (yunits.size()) std::fprintf(stdout, " [%s]", yunits.c_str());
    std::fprintf(stdout, "\n");
  }
}

static double
calc_curvi_xinc(int gridID)
{
  auto xinc = 0.0;
  if (gridInqType(gridID) != GRID_CURVILINEAR) return xinc;
  const auto *xvals2D = gridInqXvalsPtr(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);
  if (xsize > 1)
  {
    Varray<double> xvals(xsize);
    for (size_t i = 0; i < xsize; ++i) xvals[i] = xvals2D[i];
    xinc = std::fabs(xvals[xsize - 1] - xvals[0]) / (xsize - 1);
    for (size_t i = 1; i < xsize; ++i)
      if (std::fabs(std::fabs(xvals[i - 1] - xvals[i]) - xinc) > 0.005 * xinc)
      {
        xinc = 0.0;
        break;
      }
    if (is_not_equal(xinc, 0.0))
    {
      for (size_t i = 1; i < ysize; ++i)
        if (is_not_equal(xvals2D[i * xsize], xvals2D[0]) || is_not_equal(xvals2D[(i + 1) * xsize - 1], xvals2D[xsize - 1]))
        {
          xinc = 0.0;
          break;
        }
    }
  }

  return xinc;
}

static double
calc_curvi_yinc(int gridID)
{
  auto yinc = 0.0;
  if (gridInqType(gridID) != GRID_CURVILINEAR) return yinc;
  const auto *yvals2D = gridInqYvalsPtr(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);
  if (ysize > 1)
  {
    Varray<double> yvals(ysize);
    for (size_t i = 0; i < ysize; ++i) yvals[i] = yvals2D[i * xsize];
    yinc = std::fabs(yvals[ysize - 1] - yvals[0]) / (ysize - 1);
    for (size_t i = 1; i < ysize; ++i)
      if (std::fabs(std::fabs(yvals[i - 1] - yvals[i]) - yinc) > 0.005 * yinc)
      {
        yinc = 0.0;
        break;
      }
    if (is_not_equal(yinc, 0.0))
    {
      for (size_t i = 1; i < xsize; ++i)
        if (is_not_equal(yvals2D[i], yvals2D[0]) || is_not_equal(yvals2D[(ysize - 1) * xsize + i], yvals2D[(ysize - 1) * xsize]))
        {
          yinc = 0.0;
          break;
        }
    }
  }

  return yinc;
}

static void
print_xyvals2D(int gridID, int dig)
{
  if (gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL))
  {
    auto xname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
    auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
    auto yname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
    auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
    auto gridsize = gridInqSize(gridID);
    const auto *xvals2D = gridInqXvalsPtr(gridID);
    const auto *yvals2D = gridInqYvalsPtr(gridID);
    auto xmm = varray_min_max(gridsize, xvals2D);
    auto ymm = varray_min_max(gridsize, yvals2D);
    auto gridtype = gridInqType(gridID);
    auto xinc = (gridtype == GRID_CURVILINEAR) ? calc_curvi_xinc(gridID) : 0.0;
    auto yinc = (gridtype == GRID_CURVILINEAR) ? calc_curvi_yinc(gridID) : 0.0;

    std::fprintf(stdout, "%33s : %.*g", xname.c_str(), dig, xmm.min);
    if (gridsize > 1) std::fprintf(stdout, " to %.*g", dig, xmm.max);
    if (is_not_equal(xinc, 0.0)) std::fprintf(stdout, " by %.*g", dig, xinc);
    if (xunits.size()) std::fprintf(stdout, " [%s]", xunits.c_str());
    if (gridIsCircular(gridID)) std::fprintf(stdout, "  circular");
    std::fprintf(stdout, "\n");
    std::fprintf(stdout, "%33s : %.*g", yname.c_str(), dig, ymm.min);
    if (gridsize > 1) std::fprintf(stdout, " to %.*g", dig, ymm.max);
    if (is_not_equal(yinc, 0.0)) std::fprintf(stdout, " by %.*g", dig, yinc);
    if (yunits.size()) std::fprintf(stdout, " [%s]", yunits.c_str());
    std::fprintf(stdout, "\n");
  }
}

static void
printGridNumPoints(int gridtype, int gridID, size_t gridsize, size_t xsize, size_t ysize)
{
  std::fprintf(stdout, "points=%zu", gridsize);
  if (gridtype == GRID_GAUSSIAN_REDUCED)
    std::fprintf(stdout, "  nlat=%zu", ysize);
  else if (xsize && ysize)
    std::fprintf(stdout, " (%zux%zu)", xsize, ysize);

  auto numLPE = gridInqNP(gridID);
  if (numLPE > 0)
  {
    if (gridtype == GRID_GAUSSIAN) std::fprintf(stdout, "  F%d", numLPE);
    if (gridtype == GRID_GAUSSIAN_REDUCED) std::fprintf(stdout, "  N%d", numLPE);
  }
  reset_text_color(stdout);

  std::fprintf(stdout, "\n");
}

static void
print_projection_parameter(int gridID, std::string gmapname)
{
  if (gmapname.empty()) gmapname = "undefined";
  set_text_color(stdout, BLUE);
  std::fprintf(stdout, "         %24s", "mapping");
  reset_text_color(stdout);
  std::fprintf(stdout, " : ");
  set_text_color(stdout, GREEN);
  std::fprintf(stdout, "%s", gmapname.c_str());

  if (gmapname == "healpix")
  {
    auto nside = cdo::inq_att_int(gridID, CDI_GLOBAL, "healpix_nside");
    auto order = cdo::inq_att_string(gridID, CDI_GLOBAL, "healpix_order");
    std::fprintf(stdout, "  nside=%d", nside);
    std::fprintf(stdout, "  order=%s", order.size() ? order.c_str() : "undefined");
  }

  std::fprintf(stdout, "\n");
  reset_text_color(stdout);

  if (grid_has_proj_params(gridID))
  {
    auto projParams = grid_get_proj_params(gridID);
    if (projParams.size()) std::fprintf(stdout, "         %24s : %s\n", "proj_params", projParams.c_str());
  }
}

static void
print_healpix_parameter(int gridID, size_t gridSize)
{
  set_text_color(stdout, BLUE);
  std::fprintf(stdout, "         %24s", "parameter");
  reset_text_color(stdout);
  std::fprintf(stdout, " : ");
  set_text_color(stdout, GREEN);

  auto refinementLevel = cdo::inq_att_int(gridID, CDI_GLOBAL, "refinement_level");
  auto order = cdo::inq_att_string(gridID, CDI_GLOBAL, "indexing_scheme");
  std::fprintf(stdout, "refinement_level=%d", refinementLevel);
  std::fprintf(stdout, "  indexing_scheme=%s", order.size() ? order.c_str() : "undefined");

  if (gridSize == gridInqIndices(gridID, nullptr)) { std::fprintf(stdout, "  cell_index=true"); }

  std::fprintf(stdout, "\n");
  reset_text_color(stdout);
}

static void
printGridInfoKernel(int gridID, int index, int lproj)
{
  auto dig = Options::CDO_flt_digits;
  auto gridType = gridInqType(gridID);

  if (lproj && gridType != GRID_PROJECTION)
    std::fprintf(stderr, "Internal problem (%s): sub grid not equal GRID_PROJECTION!\n", __func__);

  auto trunc = gridInqTrunc(gridID);
  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);

  if (!lproj)
  {
    std::fprintf(stdout, "  %4d : ", index + 1);
    set_text_color(stdout, BLUE);
    std::fprintf(stdout, "%-24s", gridNamePtr(gridType));
    reset_text_color(stdout);
    std::fprintf(stdout, " : ");
  }

  if (gridType == GRID_LONLAT || gridType == GRID_PROJECTION || gridType == GRID_GENERIC || gridType == GRID_CHARXY
      || gridType == GRID_GAUSSIAN || gridType == GRID_GAUSSIAN_REDUCED || gridType == GRID_HEALPIX)
  {
    if (!lproj)
    {
      set_text_color(stdout, GREEN);
      printGridNumPoints(gridType, gridID, gridsize, xsize, ysize);
    }

    if (gridType == GRID_HEALPIX) { print_healpix_parameter(gridID, gridsize); }
    else
    {
      auto gmapname = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME);
      if (gridType == GRID_PROJECTION || gmapname.size()) { print_projection_parameter(gridID, gmapname); }
    }

    print_xvals(gridID, dig);
    print_yvals(gridID, dig);

    if (gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL))
    {
      std::fprintf(stdout, "%33s :", "available");
      if (gridType == GRID_GAUSSIAN_REDUCED && gridInqXvals(gridID, NULL)) std::fprintf(stdout, " xvals");
      // clang-format off
      if      (gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL)) std::fprintf(stdout, " cellbounds");
      else if (gridInqXbounds(gridID, NULL)) std::fprintf(stdout, " xbounds");
      else if (gridInqYbounds(gridID, NULL)) std::fprintf(stdout, " ybounds");
      // clang-format on
      if (gridHasArea(gridID)) std::fprintf(stdout, " area");
      if (gridInqMask(gridID, NULL)) std::fprintf(stdout, " mask");
      std::fprintf(stdout, "\n");
    }
  }
  else if (gridType == GRID_SPECTRAL)
  {
    set_text_color(stdout, GREEN);
    std::fprintf(stdout, "points=%zu  nsp=%zu  T%d", gridsize, gridsize / 2, trunc);
    if (gridInqComplexPacking(gridID)) std::fprintf(stdout, "  complexPacking");
    reset_text_color(stdout);
    std::fprintf(stdout, "\n");
  }
  else if (gridType == GRID_FOURIER)
  {
    set_text_color(stdout, GREEN);
    std::fprintf(stdout, "points=%zu  nfc=%zu  T%d\n", gridsize, gridsize / 2, trunc);
    reset_text_color(stdout);
  }
  else if (gridType == GRID_GME)
  {
    int nd, ni, ni2, ni3;
    gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
    set_text_color(stdout, GREEN);
    std::fprintf(stdout, "points=%zu  nd=%d  ni=%d\n", gridsize, nd, ni);
    reset_text_color(stdout);
  }
  else if (gridType == GRID_CURVILINEAR || gridType == GRID_UNSTRUCTURED)
  {
    set_text_color(stdout, GREEN);
    if (gridType == GRID_CURVILINEAR) { std::fprintf(stdout, "points=%zu (%zux%zu)", gridsize, xsize, ysize); }
    else { std::fprintf(stdout, "points=%zu", gridsize); }

    if (gridType == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0) std::fprintf(stdout, "  nvertex=%d", gridInqNvertex(gridID));
    reset_text_color(stdout);

    std::fprintf(stdout, "\n");

    if (gridType == GRID_UNSTRUCTURED)
    {
      int number = 0;
      cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
      int position = 0;
      cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
      if (number > 0) std::fprintf(stdout, "%33s : number=%d  position=%d\n", "grid", number, position);

      int length = 0;
      if (CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
      {
        char referenceLink[8192];
        cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, referenceLink, &length);
        std::fprintf(stdout, "%33s : %s\n", "uri", referenceLink);
      }
    }

    print_xyvals2D(gridID, dig);
  }
  else  // if ( gridtype == GRID_GENERIC )
  {
    set_text_color(stdout, GREEN);
    if (ysize == 0)
      std::fprintf(stdout, "points=%zu\n", gridsize);
    else
      std::fprintf(stdout, "points=%zu (%zux%zu)\n", gridsize, xsize, ysize);
    reset_text_color(stdout);
  }

  if (gridType == GRID_CURVILINEAR || gridType == GRID_UNSTRUCTURED)
  {
    if (gridHasArea(gridID) || gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL))
    {
      std::fprintf(stdout, "%33s :", "available");
      if (gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL)) std::fprintf(stdout, " cellbounds");
      if (gridHasArea(gridID)) std::fprintf(stdout, " area");
      if (gridInqMask(gridID, NULL)) std::fprintf(stdout, " mask");
      std::fprintf(stdout, "\n");
    }
  }

  unsigned char uuid[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  auto status = cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
  if (status == CDI_NOERR && !cdiUUIDIsNull(uuid))
  {
    char uuidStr[uuidNumHexChars + 1] = { 0 };
    if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) std::fprintf(stdout, "%33s : %s\n", "uuid", uuidStr);
  }

  if (Options::cdoVerbose)
  {
    int datatype;
    cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
    std::fprintf(stdout, "%33s : %s\n", "datatype", cdo::datatype_to_cstr(datatype));
    std::fprintf(stdout, "%33s : %d\n", "gridID", gridID);
  }
}

void
print_grid_info(int vlistID)
{
  auto numGrids = vlistNumGrids(vlistID);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(vlistID, index);
    printGridInfoKernel(gridID, index, false);
    auto projID = gridInqProj(gridID);
    if (projID != CDI_UNDEFID) printGridInfoKernel(projID, index, true);
  }
}

static void
printZaxisBoundsInfo(int zaxisID, int dig, int levelsize, double zinc, std::string const &zunits)
{
  auto level1 = zaxisInqLbound(zaxisID, 0);
  auto level2 = zaxisInqUbound(zaxisID, 0);
  if (!(levelsize == 1 && is_equal(level1, level2) && std::fabs(level1) <= 0))
  {
    std::fprintf(stdout, "%33s : ", "bounds");
    std::fprintf(stdout, "%.*g-%.*g", dig, level1, dig, level2);
    if (levelsize > 1)
    {
      level1 = zaxisInqLbound(zaxisID, levelsize - 1);
      level2 = zaxisInqUbound(zaxisID, levelsize - 1);
      std::fprintf(stdout, " to %.*g-%.*g", dig, level1, dig, level2);
      if (is_not_equal(zinc, 0)) std::fprintf(stdout, " by %.*g", dig, zinc);
    }
    if (zunits.size()) std::fprintf(stdout, " [%s]", zunits.c_str());
    std::fprintf(stdout, "\n");
  }
}

static bool
zaxisTypeIsSingleLayer(int zaxistype)
{
  switch (zaxistype)
  {
    case ZAXIS_MEANSEA:
    case ZAXIS_TROPOPAUSE:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_ATMOSPHERE:
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_LAKE_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM_TA:
    case ZAXIS_SEDIMENT_BOTTOM_TW:
    case ZAXIS_SURFACE: return true;
  }

  return false;
}

static void
printZaxisLevelInfo(int levelsize, int zaxisID, int zaxistype, double &zinc, int dig, std::string const &zname,
                    std::string const &zunits)
{
  Varray<double> levels(levelsize);
  zaxisInqLevels(zaxisID, levels.data());

  if (!(zaxisTypeIsSingleLayer(zaxistype) && levelsize == 1 && std::fabs(levels[0]) <= 0))
  {
    auto zfirst = levels[0];
    auto zlast = levels[levelsize - 1];
    if (levelsize > 2)
    {
      zinc = (levels[levelsize - 1] - levels[0]) / (levelsize - 1);
      for (int levelID = 2; levelID < levelsize; ++levelID)
        if (std::fabs(std::fabs(levels[levelID] - levels[levelID - 1]) - zinc) > 0.001 * zinc)
        {
          zinc = 0;
          break;
        }
    }

    std::fprintf(stdout, "%33s : %.*g", zname.c_str(), dig, zfirst);
    if (levelsize > 1)
    {
      std::fprintf(stdout, " to %.*g", dig, zlast);
      if (is_not_equal(zinc, 0.0)) std::fprintf(stdout, " by %.*g", dig, zinc);
    }
    if (zunits.size()) std::fprintf(stdout, " [%s]", zunits.c_str());
    std::fprintf(stdout, "\n");
  }
}

static void
printZaxisHybridInfo(int zaxisID)
{
  auto psname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_PSNAME);
  int vctsize = zaxisInqVctSize(zaxisID);
  if (vctsize || psname.size())
  {
    std::fprintf(stdout, "%33s :", "available");
    if (vctsize) std::fprintf(stdout, " vct");
    if (psname.size()) std::fprintf(stdout, "  ps: %s", psname.c_str());
    std::fprintf(stdout, "\n");
  }
}

static void
printZaxisGenericInfo(int ltype, int zaxistype, const char *zaxisname)
{
  if (zaxistype == ZAXIS_GENERIC && ltype != 0) { std::fprintf(stdout, "%-12s (ltype=%3d)", zaxisname, ltype); }
  else { std::fprintf(stdout, "%-24s", zaxisname); }
}

static void
printZaxisReferenceInfo(int zaxisID)
{
  int number = 0;
  // cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, &number)
  // if (number > 0)
  {
    std::fprintf(stdout, "%33s : ", "zaxis");
    std::fprintf(stdout, "number=%d\n", number);
  }

  unsigned char uuid[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  auto status = cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
  if (status == CDI_NOERR && !cdiUUIDIsNull(uuid))
  {
    char uuidStr[uuidNumHexChars + 1] = { 0 };
    if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) std::fprintf(stdout, "%33s : %s\n", "uuid", uuidStr);
  }
}

void
print_zaxis_info(int vlistID)
{
  auto dig = Options::CDO_flt_digits;
  char zaxisname[CDI_MAX_NAME];

  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    auto zaxistype = zaxisInqType(zaxisID);
    auto levelsize = zaxisInqSize(zaxisID);
    int ltype = 0, ltype2 = -1;
    cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
    cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, &ltype2);

    zaxisName(zaxistype, zaxisname);
    auto zname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
    auto zunits = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
    if (zunits.size() > 12) zunits.resize(12);

    std::fprintf(stdout, "  %4d : ", vlistZaxisIndex(vlistID, zaxisID) + 1);
    set_text_color(stdout, BLUE);
    printZaxisGenericInfo(ltype, zaxistype, zaxisname);

    reset_text_color(stdout);

    std::fprintf(stdout, " :");

    set_text_color(stdout, GREEN);
    std::fprintf(stdout, " levels=%d", levelsize);
    int zscalar = (levelsize == 1) ? zaxisInqScalar(zaxisID) : false;
    if (zscalar) std::fprintf(stdout, "  scalar");
    reset_text_color(stdout);
    std::fprintf(stdout, "\n");

    double zinc = 0.0;
    if (zaxisInqLevels(zaxisID, NULL)) printZaxisLevelInfo(levelsize, zaxisID, zaxistype, zinc, dig, zname, zunits);

    if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
      printZaxisBoundsInfo(zaxisID, dig, levelsize, zinc, zunits);

    if (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) printZaxisHybridInfo(zaxisID);

    if (zaxistype == ZAXIS_REFERENCE) printZaxisReferenceInfo(zaxisID);

    if (ltype != ltype2 && ltype2 != -1) std::fprintf(stdout, "%33s : %d\n", "typeOfSecondFixedSurface", ltype2);

    if (Options::cdoVerbose)
    {
      int datatype;
      cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
      std::fprintf(stdout, "%33s : %s\n", "datatype", cdo::datatype_to_cstr(datatype));
      std::fprintf(stdout, "%33s : %d\n", "zaxisID", zaxisID);
    }
  }
}

void
print_subtype_info(int vlistID)
{
  auto nsubtypes = vlistNsubtypes(vlistID);
  for (int index = 0; index < nsubtypes; ++index)
  {
    auto subtypeID = vlistSubtype(vlistID, index);
    auto subtypesize = subtypeInqSize(subtypeID);
    // subtypePrint(subtypeID);
    std::fprintf(stdout, "  %4d : %-24s :", vlistSubtypeIndex(vlistID, subtypeID) + 1, "tiles");
    std::fprintf(stdout, " ntiles=%d", subtypesize);
    std::fprintf(stdout, "\n");
  }
}

static int
printDateTime(int ntimeout, CdiDateTime vDateTime)
{
  if (ntimeout == 4)
  {
    ntimeout = 0;
    std::fprintf(stdout, "\n");
  }

  std::fprintf(stdout, " %s %s", date_to_string(vDateTime.date).c_str(), time_to_string(vDateTime.time).c_str());

  return ++ntimeout;
}

constexpr int NumTimestep = 60;

static int
printDot(int ndotout, int *nfact, int *ncout)
{
  constexpr int MaxDots = 80;

  // printf("ncout %d %d %d\n",*ncout, (*ncout)%(*nfact), *nfact);
  if ((*ncout) % (*nfact) == 0)
  {
    if (ndotout == MaxDots)
    {
      *ncout = 0;
      ndotout = 0;
      std::fprintf(stdout, "\n   ");
      (*nfact) *= 10;
    }

    std::fprintf(stdout, ".");
    fflush(stdout);
    ndotout++;
  }

  (*ncout)++;

  return ndotout;
}

void
print_timesteps(CdoStreamID streamID, int taxisID, int verbose)
{
  struct datetime
  {
    CdiDateTime vDateTime{};
    struct datetime *next = nullptr;
  };
  struct datetime vdatetime[NumTimestep];
  struct datetime *next_vdatetime = vdatetime;

  for (int i = 0; i < NumTimestep - 1; ++i) vdatetime[i].next = &vdatetime[i + 1];
  vdatetime[NumTimestep - 1].next = &vdatetime[0];

  int ntimeout = 0;
  int ndotout = 0;
  int nvdatetime = 0;
  int ncout = 0;
  int nfact = 1;
  int tsID = 0;

  DateTimeList dtlist;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(streamID, tsID);
    if (numFields == 0) break;

    dtlist.taxis_inq_timestep(taxisID, 0);
    auto vDateTime = dtlist.vDateTime(0);

    if (verbose || tsID < NumTimestep) { ntimeout = printDateTime(ntimeout, vDateTime); }
    else
    {
      if (tsID == 2 * NumTimestep) std::fprintf(stdout, "\n   ");
      if (tsID >= 2 * NumTimestep) ndotout = printDot(ndotout, &nfact, &ncout);

      if (nvdatetime < NumTimestep)
      {
        vdatetime[nvdatetime].vDateTime = vDateTime;
        nvdatetime++;
      }
      else
      {
        next_vdatetime->vDateTime = vDateTime;
        next_vdatetime = next_vdatetime->next;
      }
    }

    tsID++;
  }

  if (nvdatetime)
  {
    std::fprintf(stdout, "\n");

    ntimeout = 0;
    int toff = 0;
    if (tsID > 2 * NumTimestep)
    {
      toff = tsID % 4;
      if (toff > 0) toff = 4 - toff;
      for (int i = 0; i < toff; ++i) next_vdatetime = next_vdatetime->next;
    }
    for (int i = toff; i < nvdatetime; ++i)
    {
      auto vDateTime = next_vdatetime->vDateTime;
      ntimeout = printDateTime(ntimeout, vDateTime);
      next_vdatetime = next_vdatetime->next;
    }
  }
}
