/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#include <cdi.h>

#include <cstring>
#include <string>

#include "cdi_uuid.h"
#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"

constexpr int MaxLen = 120;

static void
printDblsPrefixAutoBrk(std::FILE *fp, int dig, std::string const &prefix, size_t n, std::vector<double> const &vals,
                       size_t extbreak)
{
  int nbyte0 = (int) prefix.size();
  fputs(prefix.c_str(), fp);
  int nbyte = nbyte0;
  for (size_t i = 0; i < n; ++i)
  {
    if (nbyte > MaxLen || (i && i == extbreak))
    {
      std::fprintf(fp, "\n%*s", nbyte0, "");
      nbyte = nbyte0;
    }
    nbyte += std::fprintf(fp, "%.*g ", dig, vals[i]);
  }
  fputs("\n", fp);
}

static void
zaxis_print_kernel(int zaxisID, std::FILE *fp)
{
  auto type = zaxisInqType(zaxisID);
  auto nlevels = zaxisInqSize(zaxisID);
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
  auto nvals = (size_t) zaxisInqLevels(zaxisID, nullptr);

  auto dig = (datatype == CDI_DATATYPE_FLT64) ? Options::CDO_dbl_digits : Options::CDO_flt_digits;

  std::fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  std::fprintf(fp, "size      = %d\n", nlevels);
  // clang-format off
  if      (datatype == CDI_DATATYPE_FLT32) std::fprintf(fp, "datatype  = float\n");
  else if (datatype == CDI_DATATYPE_INT32) std::fprintf(fp, "datatype  = int\n");
  else if (datatype == CDI_DATATYPE_INT16) std::fprintf(fp, "datatype  = short\n");
  // clang-format on

  if (nlevels == 1 && zaxisInqScalar(zaxisID)) std::fprintf(fp, "scalar    = true\n");

  auto zname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_NAME);
  auto zlongname = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME);
  auto zunits = cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS);
  if (zname.size()) std::fprintf(fp, "name      = %s\n", zname.c_str());
  if (zlongname.size()) std::fprintf(fp, "longname  = \"%s\"\n", zlongname.c_str());
  if (zunits.size()) std::fprintf(fp, "units     = \"%s\"\n", zunits.c_str());

  std::vector<double> vals;
  if (nvals) vals.resize(nvals);

  if (nvals)
  {
    zaxisInqLevels(zaxisID, vals.data());
    printDblsPrefixAutoBrk(fp, dig, "levels    = ", nvals, vals, 0);
  }
  else if (type == ZAXIS_CHAR)
  {
    auto clen = zaxisInqCLen(zaxisID);
    char **cvals = nullptr;
    zaxisInqCVals(zaxisID, &cvals);
    std::fprintf(fp, "levels    = \n");
    for (int i = 0; i < nlevels; ++i)
    {
      std::fprintf(fp, "     [%2d] = %.*s\n", i, clen, cvals[i]);
      std::free(cvals[i]);
    }
    if (cvals) std::free(cvals);
  }

  if (zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
  {
    zaxisInqLbounds(zaxisID, vals.data());
    printDblsPrefixAutoBrk(fp, dig, "lbounds   = ", nvals, vals, 0);

    zaxisInqUbounds(zaxisID, vals.data());
    printDblsPrefixAutoBrk(fp, dig, "ubounds   = ", nvals, vals, 0);
  }

  if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
  {
    auto vctsize = zaxisInqVctSize(zaxisID);
    if (vctsize)
    {
      std::fprintf(fp, "vctsize   = %d\n", vctsize);
      std::vector<double> vct(vctsize);
      zaxisInqVct(zaxisID, vct.data());
      printDblsPrefixAutoBrk(fp, dig, "vct       = ", vctsize, vct, vctsize / 2);
    }
  }

  if (type == ZAXIS_REFERENCE)
  {
    unsigned char uuid[CDI_UUID_SIZE] = { 0 };
    int length = CDI_UUID_SIZE;
    cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
    if (!cdiUUIDIsNull(uuid))
    {
      char uuidStr[uuidNumHexChars + 1] = { 0 };
      if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) std::fprintf(fp, "uuid      = %s\n", uuidStr);
    }
  }

  cdo_print_attributes(fp, zaxisID, CDI_GLOBAL, 0);
}

void
cdoPrintZaxis(int zaxisID)
{
  zaxis_print_kernel(zaxisID, stdout);
}
