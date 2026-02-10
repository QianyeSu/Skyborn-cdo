#include <cdi.h>

#include <cstdlib>
#include <cmath>
#include <vector>

// #include "cdi_int.h"
extern "C"
{
  void cdiDefTableID(int tableID);
  void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
  void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);
}

#include "cdo_default_values.h"
#include "cdo_cdi_wrapper.h"
#include "cdo_output.h"

namespace cdo
{

const char *
filetype_to_cstr(int filetype)
{
  switch (filetype)
  {
      // clang-format off
    case CDI_FILETYPE_GRB:    return "GRIB";
    case CDI_FILETYPE_GRB2:   return "GRIB2";
    case CDI_FILETYPE_NC:     return "NetCDF";
    case CDI_FILETYPE_NC2:    return "NetCDF2";
    case CDI_FILETYPE_NC4:    return "NetCDF4";
    case CDI_FILETYPE_NC4C:   return "NetCDF4 classic";
    case CDI_FILETYPE_NC5:    return "NetCDF5";
    case CDI_FILETYPE_NCZARR: return "NCZarr";
    case CDI_FILETYPE_SRV:    return "SERVICE";
    case CDI_FILETYPE_EXT:    return "EXTRA";
    case CDI_FILETYPE_IEG:    return "IEG";
    default:                  return "";
      // clang-format on
  }
}

const char *
datatype_to_cstr(int datatype)
{
  static char cstr[20]{};
  if (datatype > 0 && datatype <= 32) std::snprintf(cstr, sizeof(cstr), "P%d", datatype);

  // clang-format off
  if      (datatype == CDI_DATATYPE_PACK  ) return "P0";
  else if (datatype > 0 && datatype <= 32 ) return cstr;
  else if (datatype == CDI_DATATYPE_CPX32 ) return "C32";
  else if (datatype == CDI_DATATYPE_CPX64 ) return "C64";
  else if (datatype == CDI_DATATYPE_FLT32 ) return "F32";
  else if (datatype == CDI_DATATYPE_FLT64 ) return "F64";
  else if (datatype == CDI_DATATYPE_INT8  ) return "I8";
  else if (datatype == CDI_DATATYPE_INT16 ) return "I16";
  else if (datatype == CDI_DATATYPE_INT32 ) return "I32";
  else if (datatype == CDI_DATATYPE_UINT8 ) return "U8";
  else if (datatype == CDI_DATATYPE_UINT16) return "U16";
  else if (datatype == CDI_DATATYPE_UINT32) return "U32";
  else                                      return "";
  // clang-format on
}

int
str_to_datatype(std::string const &datatypeStr)
{
  if (datatypeStr.size() > 1)
  {
    auto ilen = atoi(datatypeStr.c_str() + 1);
    // clang-format off
    if      (datatypeStr == "P0")     return CDI_DATATYPE_PACK;
    else if (datatypeStr[0] == 'P' && ilen > 0 && ilen <= 32) return ilen;
    else if (datatypeStr == "C32")    return CDI_DATATYPE_CPX32;
    else if (datatypeStr == "C64")    return CDI_DATATYPE_CPX64;
    else if (datatypeStr == "F32")    return CDI_DATATYPE_FLT32;
    else if (datatypeStr == "F64")    return CDI_DATATYPE_FLT64;
    else if (datatypeStr == "I8")     return CDI_DATATYPE_INT8;
    else if (datatypeStr == "I16")    return CDI_DATATYPE_INT16;
    else if (datatypeStr == "I32")    return CDI_DATATYPE_INT32;
    else if (datatypeStr == "U8")     return CDI_DATATYPE_UINT8;
    else if (datatypeStr == "U16")    return CDI_DATATYPE_UINT16;
    else if (datatypeStr == "U32")    return CDI_DATATYPE_UINT32;
    else if (datatypeStr == "real")   return CDI_DATATYPE_FLT32;
    else if (datatypeStr == "double") return CDI_DATATYPE_FLT64;
    // clang-format on
  }

  return -1;
}

const char *
get_steptype_name(int tsteptype)
{
  // clang-format off
  if      (tsteptype == TSTEP_INSTANT)  return "instant";
  else if (tsteptype == TSTEP_INSTANT2) return "instant";
  else if (tsteptype == TSTEP_INSTANT3) return "instant";
  else if (tsteptype == TSTEP_MIN)      return "min";
  else if (tsteptype == TSTEP_MAX)      return "max";
  else if (tsteptype == TSTEP_AVG)      return "avg";
  else if (tsteptype == TSTEP_ACCUM)    return "accum";
  else if (tsteptype == TSTEP_RANGE)    return "range";
  else if (tsteptype == TSTEP_DIFF)     return "diff";
  else if (tsteptype == TSTEP_SUM)      return "sum";
  // clang-format on
  return "unknown";
}

}  // namespace cdo

int
cdo_taxis_create(int taxisType)
{
  if (CdoDefault::TaxisType != CDI_UNDEFID) taxisType = CdoDefault::TaxisType;
  return taxisCreate(taxisType);
}

void
cdo_taxis_copy_timestep(int taxisIDdes, int taxisIDsrc)
{
  taxisCopyTimestep(taxisIDdes, taxisIDsrc);
}

void
cdo_def_table_id(int tableID)
{
  cdiDefTableID(tableID);
}

void
grid_gen_xvals(int xsize, double xfirst, double xlast, double xinc, double *xvals)
{
  gridGenXvals(xsize, xfirst, xlast, xinc, xvals);
}

void
grid_gen_yvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals)
{
  gridGenYvals(gridtype, ysize, yfirst, ylast, yinc, yvals);
}

size_t
cdo_vlist_gridsizemax(int vlistID)
{
  size_t gridsizemax = 0;

  auto numGrids = vlistNumGrids(vlistID);
  for (int index = 0; index < numGrids; index++)
  {
    auto gridID = vlistGrid(vlistID, index);
    auto gridsize = gridInqSize(gridID);
    if (gridsize > gridsizemax) gridsizemax = gridsize;
  }

  return gridsizemax;
}

namespace cdo
{

int
inq_att_int(int cdiID, int varID, std::string const &attname)
{
  int attint = -1;
  cdiInqAttInt(cdiID, varID, attname.c_str(), 1, &attint);
  return attint;
}

std::string
inq_att_string(int cdiID, int varID, std::string const &attname)
{
  int attlen = cdiInqAttLen(cdiID, varID, attname.c_str());
  std::vector<char> atttxt(1, 0);
  if (attlen > 0)
  {
    atttxt.resize(attlen + 1);
    cdiInqAttTxt(cdiID, varID, attname.c_str(), attlen, atttxt.data());
    atttxt[attlen] = 0;
  }

  return std::string(atttxt.data());
}

int
inq_key_int(int cdiID, int varID, int key)
{
  int intValue{ 0 };
  cdiInqKeyInt(cdiID, varID, key, &intValue);
  return intValue;
}

std::string
inq_key_string(int cdiID, int varID, int key)
{
  char cstr[CDI_MAX_NAME]{};
  int length = CDI_MAX_NAME;
  cdiInqKeyString(cdiID, varID, key, cstr, &length);

  return std::string(cstr);
}

std::string
inq_var_name(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME]{};
  vlistInqVarName(vlistID, varID, cstr);
  return std::string(cstr);
}

std::string
inq_var_longname(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME]{};
  vlistInqVarLongname(vlistID, varID, cstr);
  return std::string(cstr);
}

std::string
inq_var_units(int vlistID, int varID)
{
  char cstr[CDI_MAX_NAME]{};
  vlistInqVarUnits(vlistID, varID, cstr);
  return std::string(cstr);
}

HpParams
get_healpix_params(int gridID)
{
  auto gridType = gridInqType(gridID);
  auto projection = "healpix";
  auto refinementLevel = cdo::inq_att_int(gridID, CDI_GLOBAL, "refinement_level");
  auto nside = cdo::inq_att_int(gridID, CDI_GLOBAL, "healpix_nside");
  auto paramName = (gridType == GRID_HEALPIX) ? "refinement_level" : "healpix_nside";
  if (refinementLevel == -1 && nside == -1) { cdo_abort("%s mapping parameter %s missing!", projection, paramName); }

  if (refinementLevel == -1) refinementLevel = std::log2(nside);
  if (nside == -1) nside = std::lround(std::pow(2.0, refinementLevel));

  auto indexingSchemeName = (gridType == GRID_HEALPIX) ? "indexing_scheme" : "healpix_order";
  auto order = cdo::inq_att_string(gridID, CDI_GLOBAL, indexingSchemeName);
  if (order.empty()) { cdo_abort("%s mapping parameter %s missing!", projection, indexingSchemeName); }

  if (nside < 1) cdo_abort("%s mapping parameter %s < 1!", projection, paramName);

  auto hpOrder = hp_get_order(order);
  if (hpOrder == HpOrder::Undef) cdo_abort("%s mapping parameter %s=%s unsupported!", projection, indexingSchemeName, order);

  return HpParams{ nside, hpOrder };
}

std::string
get_chunkspec_string(int vlistID, int varID)
{
  std::vector<std::tuple<int, char>> chunkList = { { CDI_KEY_CHUNKSIZE_DIMX, 'x' },
                                                   { CDI_KEY_CHUNKSIZE_DIMY, 'y' },
                                                   { CDI_KEY_CHUNKSIZE_DIMZ, 'z' },
                                                   { CDI_KEY_CHUNKSIZE_DIMT, 't' } };
  std::string chunkSpecString;
  int chunkSize;
  int numDims = 0;
  for (auto [key, dim] : chunkList)
  {
    if (cdiInqKeyInt(vlistID, varID, key, &chunkSize) == 0)
    {
      if (numDims++) chunkSpecString += ",";
      chunkSpecString += dim;
      chunkSpecString += "=" + std::to_string(chunkSize);
    }
  }

  return chunkSpecString;
}
/*
void
delete_chunks_dimT(int vlistID)
{
  constexpr int dimT{ 0 };
  constexpr int chunkKeys[] = { CDI_KEY_CHUNKSIZE_DIMT, CDI_KEY_CHUNKSIZE_DIMZ, CDI_KEY_CHUNKSIZE_DIMY, CDI_KEY_CHUNKSIZE_DIMX };
  auto numVars = vlistNvars(vlistID);
  for (int varID = 0; varID < numVars; ++varID)
  {
    if (cdo::inq_key_int(vlistID, varID, chunkKeys[dimT]) != 0)
    {
      cdiDeleteKey(vlistID, varID, chunkKeys[dimT]);
      for (int i = 1; i < 4; ++i)
      {
        if (cdo::inq_key_int(vlistID, varID, chunkKeys[i]) != 0) cdiDeleteKey(vlistID, varID, chunkKeys[i]);
      }
    }
  }
}
*/
}  // namespace cdo
