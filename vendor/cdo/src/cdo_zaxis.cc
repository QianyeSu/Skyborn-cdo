/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstring>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

#include <cdi.h>

#include "c_wrapper.h"
#include "param_conversion.h"
#include "cdo_options.h"
#include "cdo_zaxis.h"
#include "util_string.h"
#include "parse_literals.h"
#include "pmlist.h"
#include "cdo_output.h"
#include "varray.h"
#include "cdi_lockedIO.h"

namespace
{
struct ZaxisDesciption
{
  Varray<double> vals;
  Varray<double> lbounds;
  Varray<double> ubounds;
  Varray<double> vct;
  size_t vctsize;
  int type = CDI_UNDEFID;
  int datatype = CDI_UNDEFID;
  size_t size = 0;
  bool scalar = false;
  std::string name;
  std::string longname;
  std::string units;
};
}  // namespace

int
zaxisDefine(ZaxisDesciption zaxis)
{
  if (zaxis.type == CDI_UNDEFID) cdo_abort("zaxistype undefined!");
  if (zaxis.size == 0) cdo_abort("zaxis size undefined!");

  auto zaxisID = zaxisCreate(zaxis.type, (int) zaxis.size);

  if (zaxis.size == 1 && zaxis.scalar) zaxisDefScalar(zaxisID);

  if (zaxis.datatype != CDI_UNDEFID) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, zaxis.datatype);

  if (zaxis.vals.size()) zaxisDefLevels(zaxisID, zaxis.vals.data());
  if (zaxis.lbounds.size()) zaxisDefLbounds(zaxisID, zaxis.lbounds.data());
  if (zaxis.ubounds.size()) zaxisDefUbounds(zaxisID, zaxis.ubounds.data());

  if (zaxis.name.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zaxis.name.c_str());
  if (zaxis.longname.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, zaxis.longname.c_str());
  if (zaxis.units.size()) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zaxis.units.c_str());

  if (zaxis.type == ZAXIS_HYBRID || zaxis.type == ZAXIS_HYBRID_HALF)
  {
    if (zaxis.vctsize && zaxis.vct.size())
      zaxisDefVct(zaxisID, (int) zaxis.vctsize, zaxis.vct.data());
    else
      cdo_warning("vct undefined!");
  }

  return zaxisID;
}

namespace
{
struct KVMap
{
  KeyValues *kv = nullptr;
  bool isValid{ false };
};
}  // namespace

static void
zaxis_read_data(std::vector<KVMap> const &kvmap, ZaxisDesciption &zaxis, size_t &natts, std::string const &dname)
{
  for (size_t ik = 0; ik < kvmap.size(); ++ik)
  {
    if (!kvmap[ik].isValid) continue;

    auto kv = kvmap[ik].kv;
    auto const &key = kv->key;
    size_t nvalues = kv->nvalues;
    if (nvalues == 0) continue;
    auto const &values = kv->values;
    auto const &value = kv->values[0];

    // clang-format off
    if (key == "zaxistype")
    {
      auto zaxistype = parameter_to_word(value);
      if      (zaxistype == "pressure") zaxis.type = ZAXIS_PRESSURE;
      else if (zaxistype == "hybrid_half") zaxis.type = ZAXIS_HYBRID_HALF;
      else if (zaxistype == "hybrid") zaxis.type = ZAXIS_HYBRID;
      else if (zaxistype == "height") zaxis.type = ZAXIS_HEIGHT;
      else if (zaxistype == "depth_below_sea") zaxis.type = ZAXIS_DEPTH_BELOW_SEA;
      else if (zaxistype == "depth_below_land") zaxis.type = ZAXIS_DEPTH_BELOW_LAND;
      else if (zaxistype == "isentropic") zaxis.type = ZAXIS_ISENTROPIC;
      else if (zaxistype == "surface") zaxis.type = ZAXIS_SURFACE;
      else if (zaxistype == "generic") zaxis.type = ZAXIS_GENERIC;
      else if (zaxistype == "generalized_height") zaxis.type = ZAXIS_REFERENCE;
      else cdo_abort("Invalid zaxis type: %s (zaxis description file: %s)", zaxistype, dname);
    }
    else if (key == "datatype")
    {
      auto datatype = parameter_to_word(value);
      if      (datatype == "double") zaxis.datatype = CDI_DATATYPE_FLT64;
      else if (datatype == "float")  zaxis.datatype = CDI_DATATYPE_FLT32;
      else if (datatype == "int")    zaxis.datatype = CDI_DATATYPE_INT32;
      else if (datatype == "short")  zaxis.datatype = CDI_DATATYPE_INT16;
      else cdo_abort("Invalid datatype: %s (zaxis description file: %s)", datatype, dname);
    }
    else if (key == "size")     zaxis.size = parameter_to_int(value);
    else if (key == "scalar")   zaxis.scalar = parameter_to_bool(value);
    else if (key == "vctsize")  zaxis.vctsize = parameter_to_int(value);
    else if (key == "name")     zaxis.name = parameter_to_word(value);
    else if (key == "units")    zaxis.units = parameter_to_word(value);
    else if (key == "longname") zaxis.longname = value;
    else if (key == "levels")
    {
      if (zaxis.size == 0) cdo_abort("size undefined (zaxis description file: %s)!", dname);
      if (zaxis.size != nvalues) cdo_abort("size=%zu and number of levels=%zu differ!", zaxis.size, nvalues);
      zaxis.vals.resize(zaxis.size);
      for (size_t i = 0; i < zaxis.size; ++i) zaxis.vals[i] = parameter_to_double(values[i]);
    }
    else if (key == "lbounds")
    {
      if (zaxis.size == 0) cdo_abort("size undefined (zaxis description file: %s)!", dname);
      if (zaxis.size != nvalues) cdo_abort("size=%zu and number of lbounds=%zu differ!", zaxis.size, nvalues);
      zaxis.lbounds.resize(zaxis.size);
      for (size_t i = 0; i < zaxis.size; ++i) zaxis.lbounds[i] = parameter_to_double(values[i]);
    }
    else if (key == "ubounds")
    {
      if (zaxis.size == 0) cdo_abort("size undefined (zaxis description file: %s)!", dname);
      if (zaxis.size != nvalues) cdo_abort("size=%zu and number of ubounds=%zu differ!", zaxis.size, nvalues);
      zaxis.ubounds.resize(zaxis.size);
      for (size_t i = 0; i < zaxis.size; ++i) zaxis.ubounds[i] = parameter_to_double(values[i]);
    }
    else if (key == "vct")
    {
      if (zaxis.vctsize == 0) cdo_abort("vctsize undefined (zaxis description file: %s)!", dname);
      if (zaxis.vctsize != nvalues) cdo_abort("vctsize=%zu and size of vct=%zu differ!", zaxis.vctsize, nvalues);
      zaxis.vct.resize(zaxis.vctsize);
      for (size_t i = 0; i < zaxis.vctsize; ++i) zaxis.vct[i] = parameter_to_double(values[i]);
    }
    else
    {
      natts = ik;
      break;
    }
    // clang-format on
  }
}

static void
zaxis_read_attributes(size_t natts, std::vector<KVMap> const &kvmap, int zaxisID)
{
  const std::vector<std::string> reservedKeys
      = { "zaxistype", "size", "scalar", "vctsize", "name", "units", "longname", "levels", "lbounds", "ubounds", "vct" };
  const char *attkey0 = nullptr;

  for (size_t ik = natts; ik < kvmap.size(); ++ik)
  {
    if (!kvmap[ik].isValid) continue;

    const auto kv = kvmap[ik].kv;
    auto const &key = kv->key;
    size_t nvalues = kv->nvalues;
    if (nvalues == 0) continue;
    auto const &values = kv->values;
    auto const &value = kv->values[0];

    if (ik == natts)
      attkey0 = key.c_str();
    else
    {
      for (auto const &reservedKey : reservedKeys)
        if (key == reservedKey)
          cdo_abort("Found reserved keyword >%s< in attribute names! Check name or position of >%s<.", key, attkey0);
    }

    auto dtype = literals_find_datatype(nvalues, values);

    if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
    {
      std::vector<int> ivals(nvalues);
      for (size_t i = 0; i < nvalues; ++i) ivals[i] = literal_to_int(values[i]);
      cdiDefAttInt(zaxisID, CDI_GLOBAL, key.c_str(), dtype, nvalues, ivals.data());
    }
    else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)
    {
      Varray<double> dvals(nvalues);
      for (size_t i = 0; i < nvalues; ++i) dvals[i] = literal_to_double(values[i]);
      cdiDefAttFlt(zaxisID, CDI_GLOBAL, key.c_str(), dtype, nvalues, dvals.data());
    }
    else
    {
      auto len = (int) value.size();
      cdiDefAttTxt(zaxisID, CDI_GLOBAL, key.c_str(), len, value.c_str());
    }
  }
}

int
zaxis_from_file(std::FILE *zfp, std::string const &filename)
{
  PMList pmlist;
  pmlist.read_namelist(zfp, filename);
  if (pmlist.size() == 0) return -1;
  KVList &kvlist = pmlist.front();

  auto nkv = kvlist.size();
  if (nkv == 0) return -1;

  std::vector<KVMap> kvmap(nkv);
  for (size_t i = 0; i < nkv; ++i) kvmap[i].isValid = false;

  size_t ik = 0;
  const std::string firstKey = "zaxistype";
  for (auto &kv : kvlist)
  {
    if (ik == 0 && kv.key != firstKey) cdo_abort("First zaxis description keyword must be >%s< (found: %s)!", firstKey, kv.key);

    if (kv.nvalues == 0) { cdo_warning("Z-axis description keyword %s has no values, skipped!", kv.key); }
    else
    {
      kvmap[ik].isValid = true;
      kvmap[ik].kv = &kv;
    }
    ik++;
  }

  ZaxisDesciption zaxis;

  size_t natts = 0;
  zaxis_read_data(kvmap, zaxis, natts, filename);

  int zaxisID = (zaxis.type == CDI_UNDEFID) ? CDI_UNDEFID : zaxisDefine(zaxis);
  if (zaxisID != CDI_UNDEFID && natts > 0) zaxis_read_attributes(natts, kvmap, zaxisID);

  return zaxisID;
}

static void
gen_zaxis_height(ZaxisDesciption &zaxis, std::string const &zaxisname)
{
  int zaxistype = ZAXIS_HEIGHT;

  if (Options::CMOR_Mode) zaxis.scalar = true;

  double value = 0.0;
  std::vector<char> unitsString(zaxisname.size(), 0);
  auto numVals = std::sscanf(zaxisname.c_str(), "height_%lf%s", &value, unitsString.data());
  // printf("%s: %d %g%s\n", zaxisname.c_str(), numVals, value, unitsString.data());

  if (numVals != 2) return;

  std::string units = unitsString.data();

  zaxis.type = zaxistype;
  zaxis.size = 1;
  // zaxis.scalar = true;
  zaxis.vals.resize(1);
  zaxis.vals[0] = value;
  zaxis.units = units;

  auto len = units.size();
  if (len > 2 && units[len - 2] == '_' && units[len - 1] == 's')
  {
    zaxis.units.resize(len - 2);
    zaxis.scalar = true;
  }
}

int
zaxis_from_name(std::string const &zaxisname)
{
  int zaxisID = CDI_UNDEFID;

  auto zaxisNameLC = string_to_lower(zaxisname);

  ZaxisDesciption zaxis;

  if (zaxisNameLC == "surface")  // surface
  {
    zaxis.type = ZAXIS_SURFACE;
    zaxis.size = 1;
    zaxis.vals.resize(zaxis.size);
    zaxis.vals[0] = 0;
  }
  else if (zaxisNameLC.starts_with("height")) { gen_zaxis_height(zaxis, zaxisNameLC); }

  if (zaxis.type != CDI_UNDEFID) zaxisID = zaxisDefine(zaxis);

  return zaxisID;
}

int
cdo_define_zaxis(std::string zaxisfile)
{
  int zaxisID = CDI_UNDEFID;
  bool isreg = false;

  auto len = zaxisfile.size();
  int zaxisNumber = 1;
  if (len > 2 && zaxisfile[len - 2] == ':' && std::isdigit(zaxisfile[len - 1]))
  {
    zaxisNumber = zaxisfile[len - 1] - '0';
    zaxisfile.pop_back();
    zaxisfile.pop_back();
  }

  auto fileno = open(zaxisfile.c_str(), O_RDONLY);
  if (fileno >= 0)
  {
    struct stat filestat;
    if (fstat(fileno, &filestat) == 0) isreg = S_ISREG(filestat.st_mode);
  }

  if (fileno == -1 || !isreg)
  {
    if (isreg) close(fileno);
    zaxisID = zaxis_from_name(zaxisfile);
    if (zaxisID == CDI_UNDEFID) cdo_abort("Open failed on %s!", zaxisfile);
  }
  else
  {
    char buffer[4];
    if (read(fileno, buffer, 4) != 4) cdo_sys_error("Read zaxis from %s failed!", zaxisfile);

    close(fileno);

    if (zaxisID == CDI_UNDEFID)
    {
      Debug(cdoDebug, "Zaxis from CDI file");
      open_lock();
      auto streamID = streamOpenRead(zaxisfile.c_str());
      open_unlock();
      if (streamID >= 0)
      {
        auto vlistID = streamInqVlist(streamID);
        auto numZaxes = vlistNumZaxis(vlistID);
        if (zaxisNumber < 1 || zaxisNumber > numZaxes) cdo_abort("Zaxis number %d not available in %s!", zaxisNumber, zaxisfile);
        zaxisID = vlistZaxis(vlistID, zaxisNumber - 1);
        streamClose(streamID);
      }
    }

    if (zaxisID == CDI_UNDEFID)
    {
      Debug(cdoDebug, "zaxis from ASCII file");
      auto fobj = c_fopen(zaxisfile, "r");
      if (fobj.get() != nullptr) { zaxisID = zaxis_from_file(fobj.get(), zaxisfile); }
    }

    if (zaxisID == CDI_UNDEFID) cdo_abort("Invalid zaxis description file %s!", zaxisfile);
  }

  return zaxisID;
}

void
cdo_set_zaxes(std::string const &zaxisarg)
{
  auto zaxisNames = split_string(zaxisarg, ",");
  for (auto const &name : zaxisNames) cdo_define_zaxis(name);
}

static int
ztype2ltype(int zaxistype)
{
  int ltype = CDI_UNDEFID;

  // clang-format off
  if      (zaxistype == ZAXIS_SURFACE         )  ltype =   1;
  else if (zaxistype == ZAXIS_PRESSURE        )  ltype = 100;
  else if (zaxistype == ZAXIS_ALTITUDE        )  ltype = 103;
  else if (zaxistype == ZAXIS_HEIGHT          )  ltype = 105;
  else if (zaxistype == ZAXIS_SIGMA           )  ltype = 107;
  else if (zaxistype == ZAXIS_HYBRID          )  ltype = 109;
  else if (zaxistype == ZAXIS_HYBRID_HALF     )  ltype = 109;
  else if (zaxistype == ZAXIS_DEPTH_BELOW_LAND)  ltype = 111;
  else if (zaxistype == ZAXIS_ISENTROPIC      )  ltype = 113;
  else if (zaxistype == ZAXIS_DEPTH_BELOW_SEA )  ltype = 160;
  // clang-format on

  return ltype;
}

int
zaxis_to_ltype(int zaxisID)
{
  int ltype = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
  if (ltype <= 0) ltype = ztype2ltype(zaxisInqType(zaxisID));

  return ltype;
}

void
gen_layer_bounds(int nlev, Varray<double> const &levels, Varray<double> &lbounds, Varray<double> &ubounds)
{
  if (nlev == 1)
  {
    lbounds[0] = 0.0;
    ubounds[0] = 1.0;
  }
  else
  {
    lbounds[0] = levels[0];
    ubounds[nlev - 1] = levels[nlev - 1];
    for (int i = 0; i < nlev - 1; ++i)
    {
      auto bound = 0.5 * (levels[i] + levels[i + 1]);
      lbounds[i + 1] = bound;
      ubounds[i] = bound;
    }
  }
}

int
get_layer_thickness(bool useWeights, bool genBounds, int index, int zaxisID, int nlev, Varray<double> &thickness,
                    Varray<double> &weights)
{
  int status = 0;
  Varray<double> levels(nlev), lbounds(nlev, 0.0), ubounds(nlev, 1.0);

  cdo_zaxis_inq_levels(zaxisID, levels.data());
  if (genBounds)
  {
    status = 2;
    gen_layer_bounds(nlev, levels, lbounds, ubounds);
  }
  else if (useWeights && zaxisInqLbounds(zaxisID, nullptr) && zaxisInqUbounds(zaxisID, nullptr))
  {
    status = 1;
    zaxisInqLbounds(zaxisID, lbounds.data());
    zaxisInqUbounds(zaxisID, ubounds.data());
  }

  for (int i = 0; i < nlev; ++i) thickness[i] = std::fabs(ubounds[i] - lbounds[i]);

  auto layerSum = varray_sum(nlev, thickness);
  varray_copy(nlev, thickness, weights);
  varray_divc(nlev, weights, (layerSum / nlev));
  auto weightSum = varray_sum(nlev, weights);

  if (Options::cdoVerbose)
  {
    cdo_print("zaxisID=%d  nlev=%d  layersum=%g  weightsum=%g", index, nlev, layerSum, weightSum);
    printf("         level     bounds   thickness  weight\n");
    for (int i = 0; i < nlev; ++i)
      printf("   %3d  %6g  %6g/%-6g  %6g  %6g\n", i + 1, levels[i], lbounds[i], ubounds[i], thickness[i], weights[i]);
  }

  return status;
}
