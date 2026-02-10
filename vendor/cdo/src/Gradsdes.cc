/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Gradsdes   gradsdes       GrADS data descriptor file
*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* _FILE_OFFSET_BITS influence off_t */
#endif

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"
#include <mpim_grid.h>
#include "util_fileextensions.h"
#include "util_files.h"
#include "cdo_zaxis.h"

/*
  Values output into the grib map file:

  Header:

  hipnt info:  0 - version number (1)
               1 - number of times in file
               2 - number of records per time
               3 - Grid type
                 255 - user defined grid.  descriptor
                       describes grid exactly; one record
                       per grid.
                  29 - Predefined grid set 29 and 30.
                       Two records per grid.

  hfpnt info:  None

  Info:

  intpnt info (for each mapped grib record) :
                 0 - position of start of data in file
                 1 - position of start of bit map in file
                 2 - number of bits per data element

  fltpnt info :
                 0 - decimal scale factor for this record
                 1 - binary scale factor
                 2 - reference value

*/

namespace
{
struct gaindx
{
  int type;      /* Indexing file type             */
  int hinum;     /* Number of ints in header       */
  int hfnum;     /* Number of floats in header     */
  int intnum;    /* Number of index ints (long)    */
  int fltnum;    /* Number of index floats         */
  int *hipnt;    /* Pointer to header int values   */
  float *hfpnt;  /* Pointer to header float values */
  int *intpnt;   /* Pointer to int index values    */
  float *fltpnt; /* Pointer to float index values  */
};

struct gaindxb
{
  int bignum;    /* Number of off_t values */
  off_t *bigpnt; /* Pointer to off_t values */
};
}  // namespace

/* Byte swap requested number of 4 byte elements */

static void
gabswp(void *r, int cnt)
{
  char *ch1 = (char *) r;
  char *ch2 = ch1 + 1;
  char *ch3 = ch2 + 1;
  char *ch4 = ch3 + 1;
  for (int i = 0; i < cnt; ++i)
  {
    char cc1 = *ch1;
    char cc2 = *ch2;
    *ch1 = *ch4;
    *ch2 = *ch3;
    *ch3 = cc2;
    *ch4 = cc1;
    ch1 += 4;
    ch2 += 4;
    ch3 += 4;
    ch4 += 4;
  }
}

/*
 * convert an IBM float to single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 */
static float
ibm2flt(unsigned char *ibm)
{
  int positive = ((ibm[0] & 0x80) == 0);
  long mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
  int power = (int) (ibm[0] & 0x7f) - 64;
  unsigned int abspower = (power > 0) ? power : -power;

  /* calc exp */
  double exp = 16.0;
  double value = 1.0;
  while (abspower)
  {
    if (abspower & 1) value *= exp;
    exp = exp * exp;
    abspower >>= 1;
  }

  if (power < 0) value = 1.0 / value;
  value = value * mant / 16777216.0;
  if (positive == 0) value = -value;

  return (float) value;
}

/*
 * convert a float to an IBM single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 *
 * doesn't handle subnormal numbers
 */
static int
flt2ibm(float x, unsigned char *ibm)
{
  int sign, exp, i;

  if (std::fabs((double) x) <= 0)
  {
    ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
    return 0;
  }

  /* sign bit */
  if (x < 0.0)
  {
    sign = 128;
    x = -x;
  }
  else
    sign = 0;

  double mant = frexp((double) x, &exp);

  /* round up by adding 2**-24 */
  /* mant = mant + 1.0/16777216.0; */

  if (mant >= 1.0)
  {
    mant = 0.5;
    exp++;
  }
  while (exp & 3)
  {
    mant *= 0.5;
    exp++;
  }

  exp = exp / 4 + 64;

  if (exp < 0)
  {
    fprintf(stderr, "underflow in flt2ibm\n");
    ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
    return 0;
  }
  if (exp > 127)
  {
    fprintf(stderr, "overflow in flt2ibm\n");
    ibm[0] = sign | 127;
    ibm[1] = ibm[2] = ibm[3] = 255;
    return -1;
  }

  /* normal number */

  ibm[0] = sign | exp;

  mant = mant * 256.0;
  i = (int) std::floor(mant);
  mant = mant - i;
  ibm[1] = i;

  mant = mant * 256.0;
  i = (int) std::floor(mant);
  mant = mant - i;
  ibm[2] = i;

  ibm[3] = (int) std::floor(mant * 256.0);

  return 0;
}

#define GET_UINT4(a, b, c, d) ((int) ((a << 24) + (b << 16) + (c << 8) + (d)))
#define Put1Byte(buf, cnt, ival) (buf[cnt++] = (ival))
#define Put2Byte(buf, cnt, ival) ((buf[cnt++] = (ival) >> 8), (buf[cnt++] = (ival)))
#define Put4Byte(buf, cnt, ival) \
  ((buf[cnt++] = (ival) >> 24), (buf[cnt++] = (ival) >> 16), (buf[cnt++] = (ival) >> 8), (buf[cnt++] = (ival)))

#define PutInt(buf, cnt, ival) ((ival < 0) ? Put4Byte(buf, cnt, 0x7fffffffU - ival + 1) : Put4Byte(buf, cnt, ival))

static void
dumpmap()
{
  unsigned char urec[4];
  unsigned char vermap;
  unsigned char mrec[512];
  int swpflg = 0;
  int nrecords = 0;
  struct gaindx indx;
  struct gaindxb indxb;
  [[maybe_unused]] size_t nbytes;

  indx.hipnt = nullptr;
  indx.hfpnt = nullptr;
  indx.intpnt = nullptr;
  indx.fltpnt = nullptr;
  indxb.bigpnt = nullptr;
  indxb.bignum = 0;

  auto mapfp = std::fopen(cdo_get_stream_name(0), "r");
  if (mapfp == nullptr) cdo_abort("Open failed on %s", cdo_get_stream_name(0));

  /* check the version number */

  fseek(mapfp, 1, 0);
  nbytes = fread(&vermap, sizeof(unsigned char), 1, mapfp);

  if (vermap == 0) vermap = 1;

  printf("gribmap version = %d\n", vermap);

  if (vermap == 2)
  {
    fseek(mapfp, 2, 0);

    nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
    indx.hinum = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);

    nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
    indx.hfnum = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);

    nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
    indx.intnum = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);

    nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
    indx.fltnum = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);

    nbytes = fread(mrec, sizeof(unsigned char), 7, mapfp);

    if (indx.hinum > 0)
    {
      indx.hipnt = new int[indx.hinum];
      for (int i = 0; i < indx.hinum; ++i)
      {
        nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
        indx.hipnt[i] = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);
      }
    }
    if (indx.hfnum > 0)
    {
      indx.hfpnt = new float[indx.hfnum];
      nbytes = fread(indx.hfpnt, sizeof(float), indx.hfnum, mapfp);
    }
    if (indx.intnum > 0)
    {
      indx.intpnt = new int[indx.intnum];
      for (int i = 0; i < indx.intnum; ++i)
      {
        nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
        indx.intpnt[i] = GET_UINT4(mrec[0], mrec[1], mrec[2], mrec[3]);
        if (indx.intpnt[i] < 0) indx.intpnt[i] = 0x7fffffffU - indx.intpnt[i] + 1;
      }
    }
    if (indx.fltnum > 0)
    {
      indx.fltpnt = new float[indx.fltnum];
      for (int i = 0; i < indx.fltnum; ++i)
      {
        nbytes = fread(urec, sizeof(unsigned char), 4, mapfp);
        indx.fltpnt[i] = ibm2flt(urec);
      }
    }
  }
  else
  {
    fseek(mapfp, 0, 0);
    nbytes = fread(&indx, sizeof(struct gaindx), 1, mapfp);

    if (indx.type >> 24 > 0) swpflg = 1;
    if (swpflg) printf("swap endian!\n");
    if (swpflg) gabswp((float *) &indx.type, 5);

    if (indx.hinum > 0)
    {
      indx.hipnt = new int[indx.hinum];
      nbytes = fread(indx.hipnt, sizeof(int), indx.hinum, mapfp);
      if (swpflg) gabswp((float *) (indx.hipnt), indx.hinum);
    }
    if (indx.hfnum > 0)
    {
      indx.hfpnt = new float[indx.hfnum];
      nbytes = fread(indx.hfpnt, sizeof(float), indx.hfnum, mapfp);
      if (swpflg) gabswp(indx.hfpnt, indx.hfnum);
    }

    if (indx.intnum > 0)
    {
      indx.intpnt = new int[indx.intnum];
      nbytes = fread(indx.intpnt, sizeof(int), indx.intnum, mapfp);
      if (swpflg) gabswp((float *) (indx.intpnt), indx.intnum);
    }
    if (indx.fltnum > 0)
    {
      indx.fltpnt = new float[indx.fltnum];
      nbytes = fread(indx.fltpnt, sizeof(float), indx.fltnum, mapfp);
      if (swpflg) gabswp(indx.fltpnt, indx.fltnum);
    }

    if (indx.hipnt[0] == 4)
    {
      indxb.bignum = indx.hipnt[4];
      if (indxb.bignum > 0)
      {
        indxb.bigpnt = new off_t[indxb.bignum];
        nbytes = fread(indxb.bigpnt, sizeof(off_t), indxb.bignum, mapfp);
        if (swpflg) gabswp(indxb.bigpnt, indxb.bignum);
      }
    }
  }

  std::fclose(mapfp);

  printf("hinum: %d\n", indx.hinum);
  for (int i = 0; i < indx.hinum; ++i) printf("%3d %5d\n", i + 1, indx.hipnt[i]);

  printf("\n");
  printf("hfnum: %d\n", indx.hfnum);
  for (int i = 0; i < indx.hfnum; ++i) printf("%3d %g\n", i + 1, indx.hfpnt[i]);

  printf("\n");

  nrecords = indx.hipnt[1] * indx.hipnt[2];

  if (indx.intnum == indx.fltnum)
  {
    printf("num: %d\n", indx.intnum);
    for (int i = 0; i < indx.intnum / 3; ++i)
      printf("%3d %8d %6d %4d %8g %10g %8g\n", i + 1, indx.intpnt[i * 3], indx.intpnt[i * 3 + 1], indx.intpnt[i * 3 + 2],
             indx.fltpnt[i * 3], indx.fltpnt[i * 3 + 1], indx.fltpnt[i * 3 + 2]);
  }
  else if (indx.intnum == nrecords && indx.fltnum == nrecords * 3 && indxb.bignum == nrecords * 2)
  {
    printf("nrecords: %d\n", nrecords);
    for (int i = 0; i < nrecords; ++i)
      printf("%3d %8zd %6zd %4d %8g %10g %8g\n", i + 1, (size_t) indxb.bigpnt[i * 2], (size_t) indxb.bigpnt[i * 2 + 1],
             indx.intpnt[i], indx.fltpnt[i * 3], indx.fltpnt[i * 3 + 1], indx.fltpnt[i * 3 + 2]);
  }
  else
  {
    printf("intnum: %d\n", indx.intnum);
    for (int i = 0; i < indx.intnum; ++i) printf("%3d %d\n", i + 1, indx.intpnt[i]);

    printf("\n");
    printf("fltnum: %d\n", indx.fltnum);
    for (int i = 0; i < indx.fltnum; ++i) printf("%3d %g\n", i + 1, indx.fltpnt[i]);

    printf("\n");
    printf("bignum: %d\n", indxb.bignum);
    for (int i = 0; i < indxb.bignum; ++i) printf("%3d %zd\n", i + 1, (size_t) indxb.bigpnt[i]);
  }
}

static void
ctl_xydef(FILE *gdp, int gridID, bool *yrev)
{
  int j;

  *yrev = false;

  int xsize = gridInqXsize(gridID);
  int ysize = gridInqYsize(gridID);

  int gridtype = gridInqType(gridID);
  int projtype = gridInqProjType(gridID);

  /* XDEF */

  if (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_LCC)
  {
    double xmin = 1.e10, xmax = -1.e10, ymin = 1.e10, ymax = -1.e10;
    double xrange, yrange;
    int nx = 0, ny = 0, ni;
    double inc[] = { 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001 };

    auto xinc = gridInqXinc(gridID);
    auto yinc = gridInqYinc(gridID);

    struct CDI_GridProjParams gpp;
    gridInqParamsLCC(gridID, &gpp);

    fprintf(gdp, "PDEF %d %d LCCR %g %g 1 1 %g %g %g %g %g\n", xsize, ysize, gpp.xval_0, gpp.yval_0, gpp.lat_1, gpp.lat_2,
            gpp.lon_0, xinc, yinc);

    gridID = gridToCurvilinear(gridID);
    Varray<double> xvals(xsize * ysize), yvals(xsize * ysize);
    gridInqXvals(gridID, xvals.data());
    gridInqYvals(gridID, yvals.data());
    for (int i = 0; i < xsize * ysize; ++i)
    {
      if (xvals[i] > 180) xvals[i] -= 360;
      if (xvals[i] < xmin) xmin = xvals[i];
      if (xvals[i] > xmax) xmax = xvals[i];
      if (yvals[i] < ymin) ymin = yvals[i];
      if (yvals[i] > ymax) ymax = yvals[i];
    }

    double xfirst = ((int) (xmin - 0.0));
    double yfirst = ((int) (ymin - 0.0));
    xrange = ((int) (xmax + 1.5)) - xfirst;
    yrange = ((int) (ymax + 1.5)) - yfirst;

    ni = sizeof(inc) / sizeof(inc[0]);
    for (int i = 0; i < ni; ++i)
    {
      xinc = yinc = inc[i];
      nx = 1 + (int) (xrange / xinc);
      ny = 1 + (int) (yrange / yinc);

      if (nx > 1.5 * xsize && ny > 1.5 * ysize) break;
    }

    fprintf(gdp, "XDEF %d LINEAR %f %f\n", nx, xfirst, xinc);
    fprintf(gdp, "YDEF %d LINEAR %f %f\n", ny, yfirst, yinc);

    fprintf(gdp, "* XDEF 3600 LINEAR -179.95 0.1\n");
    fprintf(gdp, "* YDEF 1800 LINEAR  -89.95 0.1\n");
  }
  else
  {
    auto xfirst = gridInqXval(gridID, 0);
    auto xinc = gridInqXinc(gridID);
    if (is_equal(xinc, 0) && gridInqXvals(gridID, nullptr))
    {
      Varray<double> xvals(xsize);
      gridInqXvals(gridID, xvals.data());
      fprintf(gdp, "XDEF %d LEVELS ", xsize);
      j = 0;
      for (int i = 0; i < xsize; ++i)
      {
        fprintf(gdp, "%7.3f ", xvals[i]);
        j++;
        if (j == 6)
        {
          fprintf(gdp, "\n");
          j = 0;
          if (i != xsize - 1) fprintf(gdp, "               ");
        }
      }
      if (j) fprintf(gdp, "\n");
    }
    else
    {
      if (is_equal(xinc, 0)) xinc = 360.0 / xsize;
      fprintf(gdp, "XDEF %d LINEAR %f %f\n", xsize, xfirst, xinc);
    }
  }

  /* YDEF */

  if (!(gridtype == GRID_PROJECTION && projtype == CDI_PROJ_LCC))
  {
    auto yfirst = gridInqYval(gridID, 0);
    auto yinc = gridInqYinc(gridID);
    if (gridtype == GRID_GAUSSIAN) yinc = 0;

    if (is_equal(yinc, 0) && gridInqYvals(gridID, nullptr))
    {
      Varray<double> yvals(ysize);
      gridInqYvals(gridID, yvals.data());
      fprintf(gdp, "YDEF %d LEVELS ", ysize);
      j = 0;
      if (yvals[0] > yvals[ysize - 1])
      {
        *yrev = true;
        for (int i = ysize - 1; i >= 0; i--)
        {
          fprintf(gdp, "%7.3f ", yvals[i]);
          j++;
          if (j == 6)
          {
            fprintf(gdp, "\n");
            j = 0;
            if (i != 0) fprintf(gdp, "               ");
          }
        }
      }
      else
      {
        for (int i = 0; i < ysize; ++i)
        {
          fprintf(gdp, "%7.3f ", yvals[i]);
          j++;
          if (j == 6)
          {
            fprintf(gdp, "\n");
            j = 0;
            if (i != ysize - 1) fprintf(gdp, "               ");
          }
        }
      }

      if (j) fprintf(gdp, "\n");
    }
    else
    {
      if (is_equal(yinc, 0)) yinc = 180.0 / ysize;
      if (yinc < 0)
      {
        *yrev = true;
        fprintf(gdp, "YDEF %d LINEAR %f %f\n", ysize, yfirst + yinc * (ysize - 1), -yinc);
      }
      else
        fprintf(gdp, "YDEF %d LINEAR %f %f\n", ysize, yfirst, yinc);
    }
  }
}

static void
ctl_zdef(FILE *gdp, int vlistID, bool *zrev)
{
  int j;
  int zaxisIDmax = -1;
  int zaxisID, nlev;
  double levinc = 0;

  *zrev = false;
  auto numZaxes = vlistNumZaxis(vlistID);

  int nlevmax = 0;
  for (int index = 0; index < numZaxes; ++index)
  {
    zaxisID = vlistZaxis(vlistID, index);
    nlev = zaxisInqSize(zaxisID);
    if (nlev > nlevmax)
    {
      nlevmax = nlev;
      zaxisIDmax = zaxisID;
    }
  }

  Varray<double> levels(nlevmax);
  cdo_zaxis_inq_levels(zaxisIDmax, levels.data());
  bool lplev = (zaxisInqType(zaxisIDmax) == ZAXIS_PRESSURE);
  double level0 = levels[0];
  if (nlevmax > 1)
  {
    if (levels[0] < levels[1] && zaxisInqType(zaxisIDmax) != ZAXIS_HYBRID) *zrev = true;

    levinc = levels[1] - levels[0];

    if (is_equal(levinc, 1)) *zrev = false;

    for (int i = 1; i < nlevmax; ++i)
    {
      if (is_not_equal(levinc, (levels[i] - levels[i - 1])))
      {
        levinc = 0;
        break;
      }
    }
  }

  if (is_not_equal(levinc, 0))
    fprintf(gdp, "ZDEF %d LINEAR %g %g\n", nlevmax, level0, levinc);
  else
  {
    fprintf(gdp, "ZDEF %d LEVELS ", nlevmax);
    j = 0;
    /* zrev not needed !!!
    if ( *zrev )
      {
        for (i = nlevmax-1; i >= 0; i--)
          {
            if ( lplev ) fprintf(gdp, "%g ", levels[i]/100);
            else         fprintf(gdp, "%d ", (int) levels[i]);
            j++;
            if ( j == 10 )
              {
                fprintf(gdp, "\n");
                j = 0;
                if ( i != 0 ) fprintf(gdp, "               ");
              }
          }
      }
    else
    */
    {
      for (int i = 0; i < nlevmax; ++i)
      {
        if (lplev)
          fprintf(gdp, "%g ", levels[i] / 100);
        else
          fprintf(gdp, "%g ", levels[i]);
        j++;
        if (j == 10)
        {
          fprintf(gdp, "\n");
          j = 0;
          if (i != (nlevmax - 1)) fprintf(gdp, "               ");
        }
      }
    }
    if (j) fprintf(gdp, "\n");
  }
}

static void
ctl_options(FILE *gdp, bool yrev, bool zrev, bool sequential, bool bigendian, bool littleendian, bool flt64, bool cal365day)
{
  /* if ( filetype == CDI_FILETYPE_GRB ) zrev = false; */

  if (yrev || zrev || sequential || bigendian || littleendian || flt64)
  {
    fprintf(gdp, "OPTIONS");
    if (yrev) fprintf(gdp, " yrev");
    if (zrev) fprintf(gdp, " zrev");
    if (sequential) fprintf(gdp, " sequential");
    if (bigendian) fprintf(gdp, " big_endian");
    if (littleendian) fprintf(gdp, " little_endian");
    if (flt64) fprintf(gdp, " flt64");
    if (cal365day) fprintf(gdp, " 365_day_calendar");
    fprintf(gdp, "\n");
  }
}

static void
ctl_undef(FILE *gdp, double missval)
{
  fprintf(gdp, "UNDEF  %g\n", missval);
}

static void
ctl_vars(FILE *gdp, int filetype, int vlistID, VarList const &varList, int nvarsout, int *vars)
{
  char varname[CDI_MAX_NAME];

  fprintf(gdp, "VARS  %d\n", nvarsout);

  int numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    if (vars[varID] == true)
    {
      auto const &var = varList.vars[varID];
      int zaxisID = var.zaxisID;
      int ltype = 0;
      cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
      int nlev = zaxisInqSize(zaxisID);

      std::strcpy(varname, var.name.c_str());
      int len = (int) std::strlen(varname);
      int i;
      for (i = 0; i < len; ++i)
        if (varname[i] == '-') break;

      if (i < len)
        for (int j = i; j < len; ++j) varname[j] = varname[j + 1];

      fprintf(gdp, "%-15s", varname);

      if (nlev == 1) nlev = 0;

      fprintf(gdp, "  %3d", nlev);

      if (filetype == CDI_FILETYPE_GRB)
      {
        /*
        if      ( ltype == ZAXIS_SURFACE )  ltype = 1;
        else if ( ltype == ZAXIS_PRESSURE ) ltype = 99;
        else if ( nlev == 1 )  ltype = 1;
        else ltype = 99;
        */
        fprintf(gdp, "  %d,%d", var.code, ltype);
      }
      else if (filetype == CDI_FILETYPE_NC)
      {
        int xyz = vlistInqVarXYZ(vlistID, varID);

        fprintf(gdp, "  ");
        if (!var.isConstant) fprintf(gdp, "t,");
        if (xyz == 321)
        {
          if (nlev > 0) fprintf(gdp, "z,");
          fprintf(gdp, "y,x");
        }
        else if (xyz == 312)
        {
          if (nlev > 0) fprintf(gdp, "z,");
          fprintf(gdp, "x,y");
        }
        else if (xyz == 231)
        {
          fprintf(gdp, "y,");
          if (nlev > 0) fprintf(gdp, "z,");
          fprintf(gdp, "x");
        }
        else if (xyz == 132)
        {
          fprintf(gdp, "x,");
          if (nlev > 0) fprintf(gdp, "z,");
          fprintf(gdp, "y");
        }
        else
        {
          if (nlev > 0) fprintf(gdp, "z,");
          fprintf(gdp, "y,x");
        }
      }
      else
        fprintf(gdp, "  99");

      if (var.longname.size())
        fprintf(gdp, "  %s", var.longname.c_str());
      else
        fprintf(gdp, "  %s", varname);

      if (var.units.size()) fprintf(gdp, "  [%s]", var.units.c_str());

      fprintf(gdp, "\n");
    }
  }

  fprintf(gdp, "ENDVARS\n");
}

static void
write_map_grib1(const char *ctlfile, int map_version, int nrecords, int *intnum, float *fltnum, off_t *bignum)
{
  struct gaindx indx;
  struct gaindxb indxb;

  memset(&indx, 0, sizeof(struct gaindx));

  auto mapfp = std::fopen(ctlfile, "w");
  if (mapfp == nullptr) cdo_abort("Open failed on %s", ctlfile);

  indx.type = map_version;
  indx.hfnum = 0;
  if (map_version == 4)
  {
    indx.hinum = 5;
    indx.intnum = nrecords;
    indxb.bignum = 2 * nrecords;
  }
  else
  {
    indx.hinum = 4;
    indx.intnum = 3 * nrecords;
    indxb.bignum = 0;
  }
  indx.fltnum = 3 * nrecords;

  indx.hipnt = nullptr;
  indx.hfpnt = nullptr;
  indx.intpnt = nullptr;
  indx.fltpnt = nullptr;
  indxb.bigpnt = nullptr;

  int hinum[5];
  hinum[0] = map_version;
  hinum[1] = 1;
  hinum[2] = nrecords;
  hinum[3] = 255;
  hinum[4] = indxb.bignum;

  if (map_version == 2)
  {
    int nb, bcnt, rc;
    float fdum;
    unsigned char ibmfloat[4];

    /* calculate the size of the ver==1 index file */

    nb = 2 + (indx.hinum * 4) + /* version in byte 2, then 4 ints with number of each data type */
         indx.hinum * sizeof(int) + indx.hfnum * sizeof(int) + indx.intnum * sizeof(int) + indx.fltnum * sizeof(float);

    /* add additional info */

    nb += 7;     /* base time (+ sec)  for compatibility with earlier version 2 maps */
    nb += 8 * 4; /* grvals for time <-> grid conversion */

    unsigned char *map = new unsigned char[nb];

    bcnt = 0;
    Put1Byte(map, bcnt, 0);
    Put1Byte(map, bcnt, map_version);

    Put4Byte(map, bcnt, indx.hinum);
    Put4Byte(map, bcnt, indx.hfnum);
    Put4Byte(map, bcnt, indx.intnum);
    Put4Byte(map, bcnt, indx.fltnum);

    Put2Byte(map, bcnt, 0); /* initial year   */
    Put1Byte(map, bcnt, 0); /* initial month  */
    Put1Byte(map, bcnt, 0); /* initial day    */
    Put1Byte(map, bcnt, 0); /* initial hour   */
    Put1Byte(map, bcnt, 0); /* initial minute */
    Put1Byte(map, bcnt, 0); /* initial second */

    if (indx.hinum)
      for (int i = 0; i < indx.hinum; ++i) Put4Byte(map, bcnt, hinum[i]);

    if (indx.hfnum) { /* blank for now */ }

    for (int i = 0; i < nrecords; ++i)
    {
      PutInt(map, bcnt, (int) bignum[i * 2]);
      PutInt(map, bcnt, (int) bignum[i * 2 + 1]);
      PutInt(map, bcnt, intnum[i]);
    }

    for (int i = 0; i < indx.fltnum; ++i)
    {
      fdum = fltnum[i];
      rc = flt2ibm(fdum, ibmfloat);
      if (rc < 0) cdo_abort("overflow in IBM float conversion");
      for (int j = 0; j < 4; ++j) map[bcnt++] = ibmfloat[j];
    }

    /* write out the factors for converting from grid to absolute time */

    for (int i = 0; i < 8; ++i)
    {
      fdum = 0;
      rc = flt2ibm(fdum, ibmfloat);
      if (rc < 0) cdo_abort("overflow in IBM float conversion");
      for (int j = 0; j < 4; ++j) map[bcnt++] = ibmfloat[j];
    }

    fwrite(map, 1, bcnt, mapfp);

    delete[] map;
  }
  else
  {
    fwrite(&indx, sizeof(struct gaindx), 1, mapfp);
    if (indx.hinum > 0) fwrite(hinum, sizeof(int), indx.hinum, mapfp);
    if (map_version == 1)
    {
      std::vector<int> intnumbuf(indx.intnum);
      for (int i = 0; i < nrecords; ++i)
      {
        intnumbuf[i * 3 + 0] = (int) bignum[i * 2];
        intnumbuf[i * 3 + 1] = (int) bignum[i * 2 + 1];
        intnumbuf[i * 3 + 2] = intnum[i];
      }
      if (indx.intnum > 0) fwrite(intnumbuf.data(), sizeof(int), indx.intnum, mapfp);
      if (indx.fltnum > 0) fwrite(fltnum, sizeof(float), indx.fltnum, mapfp);
    }
    else
    {
      if (indx.intnum > 0) fwrite(intnum, sizeof(int), indx.intnum, mapfp);
      if (indx.fltnum > 0) fwrite(fltnum, sizeof(float), indx.fltnum, mapfp);
      if (indxb.bignum > 0) fwrite(bignum, sizeof(off_t), indxb.bignum, mapfp);
    }
  }

  std::fclose(mapfp);
}

/*
static
void write_map_grib2(const char *ctlfile, int map_version, int nrecords, int
*intnum, float *fltnum, off_t *bignum)
{
}
*/

class Gradsdes : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Gradsdes",
    .operators = { { "gradsdes", GradsdesHelp }, { "dumpmap", GradsdesHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, FilesOnly },
  };
  inline static RegisterEntry<Gradsdes> registration = RegisterEntry<Gradsdes>(module);

  int GRADSDES, DUMPMAP;
  int gridID = -1;
  int gridtype = -1;
  int index;
  char *idxfile = nullptr;
  bool yrev = false;
  bool zrev = false;
  int xyheader = 0;
  int nrecords = 0;
  bool bigendian = false, littleendian = false;
  bool sequential = false;
  bool flt64 = false;
  char Time[30], Incr[12] = { "1mn" };
  std::vector<std::string> IncrKey = { "mn", "hr", "dy", "mo", "yr" };
  int isd, imn, ihh, iyy, imm, idd;
  int isds = 0, imns = 0, ihhs = 0, iyys = 0, imms = 0, idds = 0;
  int imn0 = 0, ihh0 = 0, iyy0 = 0, imm0 = 0, idd0 = 0;
  int idmn, idhh, idmm, idyy, iddd;
  int dt = 1, iik = 0, mdt = 0;
  int prec;
  int map_version = 2;
  int maxrecs = 0;
  int monavg = -1;
  std::vector<int> intnum;
  std::vector<float> fltnum;
  std::vector<off_t> bignum;
  Varray<double> array;
  std::vector<std::string> cmons = { "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec" };

  CdoStreamID streamID;

  int operatorID;
  int taxisID;
  int vlistID;

  int filetype;
  int nvarsout = 0, numFieldsout = 0;

  VarList varList;

  std::vector<int> vars;
  std::vector<int> recoffset;

  bool cal365day;

  FILE *gdp;
  char *ctlfile;
  const char *datfile;

  void
  _run(int &tsID)
  {
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);

      int ms;
      if (tsID == 0)
      {
        cdiDate_decode(vDateTime.date, &iyys, &imms, &idds);
        cdiTime_decode(vDateTime.time, &ihhs, &imns, &isds, &ms);

        if (imms < 1 || imms > 12) imms = 1;

        ihh0 = ihhs;
        imn0 = imns;
        iyy0 = iyys;
        imm0 = imms;
        idd0 = idds;
      }

      if (tsID == 1)
      {
        cdiDate_decode(vDateTime.date, &iyy, &imm, &idd);
        cdiTime_decode(vDateTime.time, &ihh, &imn, &isd, &ms);

        idmn = imn - imns;
        idhh = ihh - ihhs;
        iddd = idd - idds;
        idmm = imm - imms;
        idyy = iyy - iyys;

        if (idmn != 0) { dt = idmn + (idhh + (iddd + (idmm * 30 + idyy * 12) * 30) * 24) * 60; }
        else if (idhh != 0)
        {
          dt = idhh + (iddd + (idmm + idyy * 12) * 30) * 24;
          iik = 1;
        }
        else if (iddd != 0)
        {
          dt = iddd + (idmm + idyy * 12) * 30;
          iik = 2;
        }
        else if (idmm != 0)
        {
          dt = idmm + idyy * 12;
          iik = 3;
        }
        else if (idyy != 0)
        {
          dt = idyy;
          iik = 4;
        }

        if (dt <= 0) dt = 1;
      }

      if (tsID > 0 && tsID < 6 && iik != 3 && (monavg == true || monavg == -1))
      {
        cdiDate_decode(vDateTime.date, &iyy, &imm, &idd);
        cdiTime_decode(vDateTime.time, &ihh, &imn, &isd, &ms);

        idmn = imn - imns;
        idhh = ihh - ihhs;
        iddd = idd - idds;
        idmm = imm - imms;
        idyy = iyy - iyys;

        if (iddd < 0) iddd *= -1;
        if (idyy > 0) idmm += idyy * 12;

        if (/*idmn == 0 && idhh == 0 &&*/ (iddd == 0 || iddd == 1 || idd > 27) && idmm > 0 && (mdt == 0 || idmm == mdt))
        {
          mdt = idmm;
          monavg = true;
        }
        else { monavg = false; }
        /*
        printf("monavg %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n",
               tsID, monavg, mdt, imm , imms, idmm, iyy, iyys, idyy, idd, idds, iddd);
        */
        imns = imn;
        ihhs = ihh;
        idds = idd;
        imms = imm;
        iyys = iyy;
      }

      if (filetype == CDI_FILETYPE_GRB)
      {
        nrecords += numFieldsout;
        if (nrecords >= maxrecs)
        {
          maxrecs = nrecords;
          intnum.resize(1 * maxrecs);
          fltnum.resize(3 * maxrecs);
          bignum.resize(2 * maxrecs);
        }

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID);
          if (vars[varID] == true)
          {
            size_t numMissVals;
            cdo_read_field(streamID, array.data(), &numMissVals);

            index = (tsID * numFieldsout + recoffset[varID] + levelID);

            cdo_inq_grib_info(streamID, &intnum[index], &fltnum[index * 3], &bignum[index * 2]);

            if (map_version != 4)
            {
              long checksize = (long) bignum[index * 2] + (long) varList.gridsizeMax() * intnum[index] / 8;
              if (checksize < 0L || checksize > 2147483647L)
              {
                nrecords -= numFieldsout;
                cdo_warning("File size limit reached for GrADS GRIB map_version=%d! Only the first %d time "
                            "steps (2GB) are processed.",
                            map_version, tsID);
                return;
              }
            }
          }
        }
      }

      tsID++;
    }
  }

public:
  void
  init() override
  {
    GRADSDES = module.get_id("gradsdes");
    DUMPMAP = module.get_id("dumpmap");

    (void) (GRADSDES);  // unused

    operatorID = cdo_operator_id();

    datfile = cdo_get_stream_name(0);
    auto len = std::strlen(datfile);
    ctlfile = new char[len + 10];
    std::strcpy(ctlfile, datfile);

    if (operatorID == DUMPMAP)
    {
      dumpmap();
      return;
    }

    if (cdo_operator_argc() > 1) cdo_abort("Too many arguments!");

    if (cdo_operator_argc() == 1)
    {
      map_version = parameter_to_int(cdo_operator_argv(0));
      if (map_version != 1 && map_version != 2 && map_version != 4) cdo_abort("map_version=%d unsupported!", map_version);
    }
    else
    {
      if (FileUtils::size(cdo_get_stream_name(0)) > 2147483647L) map_version = 4;
    }

    if (Options::cdoVerbose) cdo_print("GrADS GRIB map version: %d", map_version);

    if (map_version == 4 && sizeof(off_t) != 8)
      cdo_abort("GrADS GRIB map version %d requires size of off_t to be 8! The size of off_t is %ld.", map_version, sizeof(off_t));

    streamID = cdo_open_read(0);

    vlistID = cdo_stream_inq_vlist(streamID);
    varList = VarList(vlistID);

    auto numVars = varList.numVars();
    auto numSteps = varList.numSteps();
    auto numGrids = varList.numGrids();

    filetype = cdo_inq_filetype(streamID);
    auto byteorder = cdo_inq_byteorder(streamID);

    if (filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C
        || filetype == CDI_FILETYPE_NC5 || filetype == CDI_FILETYPE_NCZARR)
      filetype = CDI_FILETYPE_NC;

    if (filetype != CDI_FILETYPE_SRV && filetype != CDI_FILETYPE_EXT && filetype != CDI_FILETYPE_IEG
        && filetype != CDI_FILETYPE_GRB)
    {
      if (filetype == CDI_FILETYPE_NC)
        // cdo_abort("Unsupported file format: NetCDF");
        ;
      else if (filetype == CDI_FILETYPE_GRB2)
        // cdo_abort("Unsupported file format: GRIB2");
        ;
      else
        cdo_abort("Unsupported file format!");
    }

    // find the first lonlat or Gaussian grid
    for (index = 0; index < numGrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      gridtype = gridInqType(gridID);
      auto projtype = gridInqProjType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_LCC)
          || (gridtype == GRID_PROJECTION && projtype == CDI_PROJ_RLL))
        break;
    }

    if (index == numGrids) cdo_abort("No Lon/Lat, Gaussian or Lambert grid found (%s data unsupported)!", gridNamePtr(gridtype));

    // select all variables with used gridID
    vars = std::vector<int>(numVars);
    recoffset = std::vector<int>(numVars);
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList.vars[varID];
      if (var.gridID == gridID)
      {
        if (filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG)
        {
          prec = var.dataType;
          if (prec == CDI_DATATYPE_FLT64) flt64 = true;
        }
        vars[varID] = true;
        recoffset[varID] = numFieldsout;
        nvarsout++;
        numFieldsout += zaxisInqSize(var.zaxisID);
        if (numSteps != 1 && numSteps != 0 && var.isConstant)
          cdo_abort("Unsupported GrADS record structure! Variable %d has only 1 time step.", var.code);
      }
      else
      {
        cdo_print("Unsupported grid type >%s<, skipped variable %s!", gridNamePtr(gridInqType(var.gridID)), var.name);
        vars[varID] = false;
      }
    }

    if (filetype != CDI_FILETYPE_GRB && numVars != nvarsout) cdo_abort("Too many different grids!");

    if (filetype == CDI_FILETYPE_SRV)
    {
      xyheader = 40;
      if (flt64) xyheader = 72;
      sequential = true;
      if (byteorder == CDI_BIGENDIAN) bigendian = true;
      if (byteorder == CDI_LITTLEENDIAN) littleendian = true;
    }

    if (filetype == CDI_FILETYPE_EXT)
    {
      xyheader = 24;
      if (flt64) xyheader = 40;
      sequential = true;
      if (byteorder == CDI_BIGENDIAN) bigendian = true;
      if (byteorder == CDI_LITTLEENDIAN) littleendian = true;
    }

    if (filetype == CDI_FILETYPE_IEG)
    {
      xyheader = 644;
      if (flt64) xyheader = 1048;
      sequential = true;
      if (byteorder == CDI_BIGENDIAN) bigendian = true;
      if (byteorder == CDI_LITTLEENDIAN) littleendian = true;
    }

    // ctl file name
    repl_filetypeext(ctlfile, filetypeext(filetype), ".ctl");

    // open ctl file
    gdp = std::fopen(ctlfile, "w");
    if (gdp == nullptr) cdo_abort("Open failed on %s", ctlfile);

    // VERSION
    fprintf(gdp, "* Generated by CDO operator gradsdes\n");
    fprintf(gdp, "*\n");

    // DSET
    if (datfile[0] == '/')
      fprintf(gdp, "DSET  %s\n", datfile);
    else
    {
      datfile = strrchr(datfile, '/');
      if (datfile == 0)
        datfile = cdo_get_stream_name(0);
      else
        datfile++;
      fprintf(gdp, "DSET  ^%s\n", datfile);
    }

    /*
     * DTYPE Print file type
     * INDEX Print filename of the control/index file .ctl/.idx
     */
    if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
    {
      idxfile = strdup(ctlfile);
      char *pidxfile = idxfile;

      // print GRIB[12] file type
      // generate the index file
      if (filetype == CDI_FILETYPE_GRB)
      {
        fprintf(gdp, "DTYPE  GRIB\n");
        repl_filetypeext(idxfile, ".ctl", ".gmp");
      }
      else if (filetype == CDI_FILETYPE_GRB2)
      {
        fprintf(gdp, "DTYPE  GRIB2\n");
        repl_filetypeext(pidxfile, ".ctl", ".idx");
      }

      // print file name of index file
      if (datfile[0] == '/')
        fprintf(gdp, "INDEX  %s\n", pidxfile);
      else
      {
        pidxfile = strrchr(pidxfile, '/');
        if (pidxfile == 0)
          pidxfile = idxfile;
        else
          pidxfile++;
        fprintf(gdp, "INDEX  ^%s\n", pidxfile);
      }

      array.resize(varList.gridsizeMax());
    }
    else if (filetype == CDI_FILETYPE_NC) { fprintf(gdp, "DTYPE  NetCDF\n"); }

    // XYHEADER
    if (xyheader) fprintf(gdp, "XYHEADER  %d\n", xyheader);

    // TIME

    taxisID = vlistInqTaxis(vlistID);

    cal365day = (taxisInqCalendar(taxisID) == CALENDAR_365DAYS);
  }

  void
  run() override
  {
    if (operatorID == DUMPMAP) return;

    int tsID = 0;

    _run(tsID);

    // XYDEF
    ctl_xydef(gdp, gridID, &yrev);

    // ZDEF
    ctl_zdef(gdp, vlistID, &zrev);

    // TDEF

    if (monavg == true)
    {
      dt = mdt;
      iik = 3;
      if (idd0 > 28)
      {
        /* int iddx = idd0; */
        idd0 = 1;
        cdo_print("Reset start date to %02d:%02dZ%02d%s%04d", ihh0, imn0, idd0, cmons[imm0 - 1], iyy0);
      }
    }

    std::snprintf(Time, sizeof(Time), "%02d:%02dZ%02d%s%04d", ihh0, imn0, idd0, cmons[imm0 - 1].c_str(), iyy0);
    std::snprintf(Incr, sizeof(Incr), "%d%s", dt, IncrKey[iik].c_str());

    fprintf(gdp, "TDEF %d LINEAR %s %s\n", tsID, Time, Incr);

    // TITLE

    int xsize = gridInqXsize(gridID);
    int ysize = gridInqYsize(gridID);

    int res = (gridtype == GRID_GAUSSIAN) ? nlat_to_ntr(ysize) : 0;

    if (res)
      fprintf(gdp, "TITLE  %s  T%d grid\n", datfile, res);
    else
      fprintf(gdp, "TITLE  %s  %dx%d grid\n", datfile, xsize, ysize);

    // OPTIONS
    ctl_options(gdp, yrev, zrev, sequential, bigendian, littleendian, flt64, cal365day);

    // UNDEF
    ctl_undef(gdp, varList.vars[0].missval);

    // VARS
    ctl_vars(gdp, filetype, vlistID, varList, nvarsout, vars.data());

    // INDEX file
    if (filetype == CDI_FILETYPE_GRB)
    {
      write_map_grib1(idxfile, map_version, nrecords, intnum.data(), fltnum.data(), bignum.data());
    }
    if (filetype == CDI_FILETYPE_GRB2)
    {
      cdo_abort("The fileformat GRIB2 is not fully supported yet for the gradsdes operator.\n"
                "The .ctl file %s was generated. You can add the necessary .idx file by running\n\tgribmap -i %s",
                ctlfile, ctlfile);
      // write_map_grib2(idxfile, map_version, nrecords, intnum, fltnum, bignum);
    }
  }

  void
  close() override
  {
    if (operatorID == DUMPMAP) return;

    cdo_stream_close(streamID);

    if (ctlfile) delete[] ctlfile;
    if (idxfile) delete[] idxfile;
  }
};
