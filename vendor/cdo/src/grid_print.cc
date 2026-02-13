/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include "cdi_uuid.h"

#include <cstring>
#include <string>

#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"
#include "compare.h"
#include "util_string.h"

constexpr int MaxLen = 120;

static void
printDblsPrefixAutoBrk(std::FILE *fp, int dig, std::string const &prefix, size_t n, std::vector<double> const &values)
{
  int nbyte0 = prefix.size();
  fputs(prefix.c_str(), fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
  {
    if (nbyte > MaxLen)
    {
      std::fprintf(fp, "\n%*s", nbyte0, "");
      nbyte = nbyte0;
    }
    nbyte += std::fprintf(fp, "%.*g ", dig, values[i]);
  }
  fputs("\n", fp);
}

template <typename T>
static void
printIntsPrefixAutoBrk(std::FILE *fp, std::string const &prefix, size_t n, std::vector<T> const &values)
{
  int nbyte0 = prefix.size();
  fputs(prefix.c_str(), fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
  {
    if (nbyte > MaxLen)
    {
      std::fprintf(fp, "\n%*s", nbyte0, "");
      nbyte = nbyte0;
    }
    nbyte += std::fprintf(fp, "%lld ", static_cast<long long>(values[i]));
  }
  fputs("\n", fp);
}

static void
print_bounds(std::FILE *fp, int dig, std::string const &prefix, size_t n, size_t nvertex, std::vector<double> const &bounds)
{
  fputs(prefix.c_str(), fp);
  for (size_t i = 0; i < n; ++i)
  {
    if (i > 0) std::fprintf(fp, "\n%*s", (int) prefix.size(), "");
    for (size_t iv = 0; iv < nvertex; iv++) std::fprintf(fp, "%.*g ", dig, bounds[i * nvertex + iv]);
  }
  fputs("\n", fp);
}

static void
print_mask(std::FILE *fp, std::string const &prefix, size_t n, std::vector<int> const &mask)
{
  int nbyte0 = (int) prefix.size();
  fputs(prefix.c_str(), fp);
  auto nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
  {
    if (nbyte > MaxLen)
    {
      std::fprintf(fp, "\n%*s", nbyte0, "");
      nbyte = nbyte0;
    }
    nbyte += (size_t) std::fprintf(fp, "%d ", mask[i]);
  }
  fputs("\n", fp);
}

static void
print_attribute_txt(std::FILE *fp, int cdiID, int varID, int attlen, char const (&attname)[CDI_MAX_NAME + 1])
{
  std::vector<char> atttxt(attlen + 1);
  cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt.data());
  atttxt[attlen] = 0;
  if (std::strchr(atttxt.data(), '"'))
    std::fprintf(fp, "%s = '%s'\n", attname, atttxt.data());
  else
    std::fprintf(fp, "%s = \"%s\"\n", attname, atttxt.data());
}

static void
print_attribute_int(std::FILE *fp, int cdiID, int varID, int attlen, char const (&attname)[CDI_MAX_NAME + 1])
{
  std::vector<int> attint(attlen);
  cdiInqAttInt(cdiID, varID, attname, attlen, attint.data());
  std::fprintf(fp, "%s =", attname);
  for (int i = 0; i < attlen; ++i) std::fprintf(fp, " %d", attint[i]);
  std::fprintf(fp, "\n");
}

static void
print_attribute_flt(std::FILE *fp, int cdiID, int varID, int atttype, int attlen, char const (&attname)[CDI_MAX_NAME + 1])
{
  char fltstr[128];
  std::vector<double> attflt(attlen);
  cdiInqAttFlt(cdiID, varID, attname, attlen, attflt.data());
  std::fprintf(fp, "%s =", attname);
  if (atttype == CDI_DATATYPE_FLT32)
    for (int i = 0; i < attlen; ++i)
      std::fprintf(fp, " %sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
  else
    for (int i = 0; i < attlen; ++i)
      std::fprintf(fp, " %s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
  std::fprintf(fp, "\n");
}

static void
grid_print_attributes(std::FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  char attname[CDI_MAX_NAME + 1];

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int iatt = 0; iatt < natts; ++iatt)
  {
    int atttype, attlen;
    cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);
    if (attlen == 0) continue;

    if (cdo_cmpstr(attname, "grid_mapping_name")) continue;

    if (atttype == CDI_DATATYPE_TXT) { print_attribute_txt(fp, cdiID, varID, attlen, attname); }
    else if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
             || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
    {
      print_attribute_int(fp, cdiID, varID, attlen, attname);
    }
    else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
    {
      print_attribute_flt(fp, cdiID, varID, atttype, attlen, attname);
    }
  }
}

static void
print_xaxis(int gridID, std::FILE *fp, std::string const &xdimname)
{
  auto xname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_NAME);
  auto xlongname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_LONGNAME);
  auto xunits = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_UNITS);
  if (xname.size()) std::fprintf(fp, "xname     = %s\n", xname.c_str());
  if (xdimname.size() && xdimname != xname) std::fprintf(fp, "xdimname  = %s\n", xdimname.c_str());
  if (xlongname.size()) std::fprintf(fp, "xlongname = \"%s\"\n", xlongname.c_str());
  if (xunits.size()) std::fprintf(fp, "xunits    = \"%s\"\n", xunits.c_str());
}

static void
print_yaxis(int gridID, std::FILE *fp, std::string const &ydimname)
{
  auto yname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_NAME);
  auto ylongname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_LONGNAME);
  auto yunits = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);
  if (yname.size()) std::fprintf(fp, "yname     = %s\n", yname.c_str());
  if (ydimname.size() && ydimname != yname) std::fprintf(fp, "ydimname  = %s\n", ydimname.c_str());
  if (ylongname.size()) std::fprintf(fp, "ylongname = \"%s\"\n", ylongname.c_str());
  if (yunits.size()) std::fprintf(fp, "yunits    = \"%s\"\n", yunits.c_str());
}

static void
printf_unstructured_keys(int gridID, std::FILE *&fp)
{
  int number = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
  if (number > 0)
  {
    std::fprintf(fp, "number    = %d\n", number);
    int position = 0;
    cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
    if (position >= 0) std::fprintf(fp, "position  = %d\n", position);
  }

  int length = 0;
  if (CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
  {
    char referenceLink[8192];
    cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, referenceLink, &length);
    std::fprintf(fp, "uri       = %s\n", referenceLink);
  }
}

static void
print_xvals(int gridID, int opt, std::FILE *fp, SizeType nxvals, int type, int dig)
{
  double xfirst = 0.0, xinc = 0.0;

  if (type == GRID_LONLAT || type == GRID_GAUSSIAN || type == GRID_PROJECTION || type == GRID_GENERIC)
  {
    xfirst = gridInqXval(gridID, 0);
    xinc = gridInqXinc(gridID);
  }

  if (is_not_equal(xinc, 0.0) && opt)
  {
    std::fprintf(fp, "xfirst    = %.*g\n", dig, xfirst);
    std::fprintf(fp, "xinc      = %.*g\n", dig, xinc);
  }
  else
  {
    std::vector<double> xvals(nxvals);
    gridInqXvals(gridID, xvals.data());
    printDblsPrefixAutoBrk(fp, dig, "xvals     = ", nxvals, xvals);
  }
}

static void
print_yvals(int gridID, int opt, std::FILE *fp, SizeType nyvals, int type, int dig)
{
  double yfirst = 0.0, yinc = 0.0;

  if (type == GRID_LONLAT || type == GRID_GENERIC || type == GRID_PROJECTION || type == GRID_GENERIC)
  {
    yfirst = gridInqYval(gridID, 0);
    yinc = gridInqYinc(gridID);
  }

  if (is_not_equal(yinc, 0.0) && opt)
  {
    std::fprintf(fp, "yfirst    = %.*g\n", dig, yfirst);
    std::fprintf(fp, "yinc      = %.*g\n", dig, yinc);
  }
  else
  {
    std::vector<double> yvals(nyvals);
    gridInqYvals(gridID, yvals.data());
    printDblsPrefixAutoBrk(fp, dig, "yvals     = ", nyvals, yvals);
  }
}

static void
print_xcvals(int gridID, std::FILE *fp, SizeType xsize, int xstrlen)
{
  char **xcvals = new char *[xsize];
  for (size_t i = 0; i < xsize; ++i) xcvals[i] = new char[xstrlen + 1];
  gridInqXCvals(gridID, xcvals);
  for (size_t i = 0; i < xsize; ++i) xcvals[i][xstrlen] = 0;
  for (size_t i = 0; i < xsize; ++i)
    for (int k = xstrlen - 1; k; k--)
    {
      if (xcvals[i][k] == ' ')
        xcvals[i][k] = 0;
      else
        break;
    }

  std::fprintf(fp, "xcvals    = \"%.*s\"", xstrlen, xcvals[0]);
  for (size_t i = 1; i < xsize; ++i) std::fprintf(fp, ", \"%.*s\"", xstrlen, xcvals[i]);
  std::fprintf(fp, "\n");

  for (size_t i = 0; i < xsize; ++i) delete[] xcvals[i];
  delete[] xcvals;
}

static void
print_ycvals(int gridID, std::FILE *fp, SizeType ysize, int ystrlen)
{
  char **ycvals = new char *[ysize];
  for (size_t i = 0; i < ysize; ++i) ycvals[i] = new char[ystrlen + 1];
  gridInqYCvals(gridID, ycvals);
  for (size_t i = 0; i < ysize; ++i) ycvals[i][ystrlen] = 0;
  for (size_t i = 0; i < ysize; ++i)
    for (int k = ystrlen - 1; k; k--)
    {
      if (ycvals[i][k] == ' ')
        ycvals[i][k] = 0;
      else
        break;
    }

  std::fprintf(fp, "ycvals    = \"%.*s\"", ystrlen, ycvals[0]);
  for (size_t i = 1; i < ysize; ++i) std::fprintf(fp, ", \"%.*s\"", ystrlen, ycvals[i]);
  std::fprintf(fp, "\n");

  for (size_t i = 0; i < ysize; ++i) delete[] ycvals[i];
  delete[] ycvals;
}

static void
print_projection(int gridID, std::FILE *fp)
{
  auto gridMapping = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME);
  auto gridMappingName = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME);
  if (gridMapping.size()) std::fprintf(fp, "grid_mapping = %s\n", gridMapping.c_str());
  if (gridMappingName.size()) std::fprintf(fp, "grid_mapping_name = %s\n", gridMappingName.c_str());
  grid_print_attributes(fp, gridID);
}

static void
print_healpix(int gridID, std::FILE *fp, SizeType gridsize)
{
  grid_print_attributes(fp, gridID);
  if (gridsize == gridInqIndices(gridID, nullptr))
  {
    std::vector<int64_t> indices(gridsize);
    gridInqIndices(gridID, indices.data());
    printIntsPrefixAutoBrk(fp, "cell_index = ", gridsize, indices);
  }
}

static void
grid_print_kernel(int gridID, int opt, std::FILE *fp)
{
  auto nxvals = gridInqXvals(gridID, nullptr);
  auto nyvals = gridInqYvals(gridID, nullptr);
  auto nxbounds = gridInqXbounds(gridID, nullptr);
  auto nybounds = gridInqYbounds(gridID, nullptr);

  auto gridsize = gridInqSize(gridID);
  auto xsize = gridInqXsize(gridID);
  auto ysize = gridInqYsize(gridID);
  auto type = gridInqType(gridID);
  auto nvertex = gridInqNvertex(gridID);
  auto xstrlen = gridInqXIsc(gridID);
  auto ystrlen = gridInqYIsc(gridID);
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);

  int dig = (datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

  std::fprintf(fp, "gridtype  = %s\n", gridNamePtr(type));
  std::fprintf(fp, "gridsize  = %zu\n", gridsize);
  if (datatype == CDI_DATATYPE_FLT32) std::fprintf(fp, "datatype  = float\n");

  if (type != GRID_GME)
  {
    if (type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER)
    {
      if (xsize > 0) std::fprintf(fp, "xsize     = %zu\n", xsize);
      if (ysize > 0) std::fprintf(fp, "ysize     = %zu\n", ysize);
    }

    auto xdimname = cdo::inq_key_string(gridID, CDI_XAXIS, CDI_KEY_DIMNAME);
    if (nxvals > 0 || xstrlen) { print_xaxis(gridID, fp, xdimname); }
    else if (xsize > 0 && xdimname.size()) { std::fprintf(fp, "xdimname  = %s\n", xdimname.c_str()); }

    auto ydimname = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_DIMNAME);
    if (nyvals > 0 || ystrlen) { print_yaxis(gridID, fp, ydimname); }
    else if (ysize > 0 && ydimname.size()) { std::fprintf(fp, "ydimname  = %s\n", ydimname.c_str()); }

    if (type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR)
    {
      auto vdimName = cdo::inq_key_string(gridID, CDI_GLOBAL, CDI_KEY_VDIMNAME);
      if (vdimName.size()) std::fprintf(fp, "vdimname  = %s\n", vdimName.c_str());
    }
    if (type == GRID_UNSTRUCTURED && nvertex > 0) std::fprintf(fp, "nvertex   = %d\n", nvertex);
  }

  switch (type)
  {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_TRAJECTORY:
    case GRID_CHARXY:
    {
      if (type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED) std::fprintf(fp, "numLPE    = %d\n", gridInqNP(gridID));

      size_t xdim, ydim;
      if (type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED)
      {
        xdim = gridsize;
        ydim = gridsize;
      }
      else if (type == GRID_GAUSSIAN_REDUCED)
      {
        xdim = 2;
        ydim = ysize;
      }
      else
      {
        xdim = xsize;
        ydim = ysize;
      }

      if (type == GRID_UNSTRUCTURED) printf_unstructured_keys(gridID, fp);

      if (nxvals > 0) print_xvals(gridID, opt, fp, nxvals, type, dig);

      if (nxbounds)
      {
        std::vector<double> xbounds(nxbounds);
        gridInqXbounds(gridID, xbounds.data());
        print_bounds(fp, dig, "xbounds   = ", xdim, nvertex, xbounds);
      }

      if (xstrlen) print_xcvals(gridID, fp, xsize, xstrlen);

      if (nyvals > 0) print_yvals(gridID, opt, fp, nyvals, type, dig);

      if (nybounds)
      {
        std::vector<double> ybounds(nybounds);
        gridInqYbounds(gridID, ybounds.data());
        print_bounds(fp, dig, "ybounds   = ", ydim, nvertex, ybounds);
      }

      if (ystrlen) print_ycvals(gridID, fp, ysize, ystrlen);

      if (gridHasArea(gridID))
      {
        std::vector<double> area(gridsize);
        gridInqArea(gridID, area.data());
        printDblsPrefixAutoBrk(fp, dig, "area      = ", gridsize, area);
      }

      if (type == GRID_GAUSSIAN_REDUCED)
      {
        std::vector<int> reducedPoints(ysize);
        gridInqReducedPoints(gridID, reducedPoints.data());
        printIntsPrefixAutoBrk(fp, "reducedPoints = ", (size_t) (ysize > 0 ? ysize : 0), reducedPoints);
      }

#ifdef HIRLAM_EXTENSIONS
      {
        int scanningMode = 0;
        cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_SCANNINGMODE, &scanningMode);
        std::fprintf(fp, "scanningMode = %d\n", scanningMode);
      }
#endif

      if (type == GRID_PROJECTION) print_projection(gridID, fp);

      break;
    }
    case GRID_HEALPIX:
    {
      print_healpix(gridID, fp, gridsize);
      break;
    }
    case GRID_SPECTRAL:
    {
      std::fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
      std::fprintf(fp, "complexpacking = %d\n", gridInqComplexPacking(gridID));
      break;
    }
    case GRID_FOURIER:
    {
      std::fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
      break;
    }
    case GRID_GME:
    {
      int nd, ni, ni2, ni3;
      gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
      std::fprintf(fp, "ni        = %d\n", ni);
      break;
    }
    default:
    {
      std::fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
      break;
    }
  }

  unsigned char uuid[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
  if (!cdiUUIDIsNull(uuid))
  {
    char uuidStr[uuidNumHexChars + 1] = { 0 };
    if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) std::fprintf(fp, "uuid      = %s\n", uuidStr);
  }

  if (gridInqMask(gridID, nullptr))
  {
    std::vector<int> mask(gridsize);
    gridInqMask(gridID, mask.data());
    print_mask(fp, "mask      = ", (size_t) (gridsize > 0 ? gridsize : 0), mask);
  }

  auto projID = gridInqProj(gridID);
  if (projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION) grid_print_kernel(projID, opt, fp);
}

void
cdo_print_griddes(int gridID, int opt)
{
  grid_print_kernel(gridID, opt, stdout);
}
