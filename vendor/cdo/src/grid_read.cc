/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>
#include "cdi_uuid.h"

#include "cdo_output.h"
#include "param_conversion.h"
#include "griddes.h"
#include "parse_literals.h"
#include "pmlist.h"

namespace
{
struct KVMap
{
  KeyValues *kv = nullptr;
  bool isValid{ false };
};
}  // namespace

static void
grid_set_gridtype(GridDesciption &grid, std::string const &value, const char *dname, int ik, int &iproj)
{
  auto gridtype = parameter_to_word(value);

  if (grid.type != CDI_UNDEFID)
  {
    if (gridtype == "projection") { iproj = ik; }
    return;
  }

  if (gridtype == "lonlat") { grid.type = GRID_LONLAT; }
  else if (gridtype == "latlon") { grid.type = GRID_LONLAT; }
  else if (gridtype == "gaussian") { grid.type = GRID_GAUSSIAN; }
  else if (gridtype == "gaussian_reduced") { grid.type = GRID_GAUSSIAN_REDUCED; }
  else if (gridtype == "curvilinear") { grid.type = GRID_CURVILINEAR; }
  else if (gridtype == "unstructured") { grid.type = GRID_UNSTRUCTURED; }
  else if (gridtype == "cell") { grid.type = GRID_UNSTRUCTURED; }
  else if (gridtype == "spectral") { grid.type = GRID_SPECTRAL; }
  else if (gridtype == "gme") { grid.type = GRID_GME; }
  else if (gridtype == "projection") { grid.type = GRID_PROJECTION; }
  else if (gridtype == "generic") { grid.type = GRID_GENERIC; }
  else if (gridtype == "characterXY") { grid.type = GRID_CHARXY; }
  else if (gridtype == "healpix") { grid.type = GRID_HEALPIX; }
  else { cdo_abort("Invalid gridtype: %s (grid description file: %s)", gridtype, dname); }

  if (grid.type == GRID_LONLAT || grid.type == GRID_GAUSSIAN || grid.type == GRID_GAUSSIAN_REDUCED || grid.type == GRID_PROJECTION)
  {
    grid.nvertex = 2;
  }
  else if (grid.type == GRID_CURVILINEAR) { grid.nvertex = 4; }
}

static void
grid_set_scanningMode(GridDesciption &grid, std::string const &value)
{
  grid.scanningMode = 64;
  auto scmode = parameter_to_int(value);
  // -1: not used; allowed modes: <0, 64, 96>; Default is 64
  if (scmode == 0 || scmode == 64 || scmode == 96) { grid.scanningMode = scmode; }
  else { cdo_warning("Warning: %d not in allowed modes: <0, 64, 96>; Using default: 64\n", scmode); }
}

static void
grid_set_datatype(GridDesciption &grid, std::string const &value, const char *dname)
{
  auto datatype = parameter_to_word(value);
  if (datatype == "double") { grid.datatype = CDI_DATATYPE_FLT64; }
  else if (datatype == "float") { grid.datatype = CDI_DATATYPE_FLT32; }
  else { cdo_abort("Invalid datatype: %s (grid description file: %s)", datatype, dname); }
}

static void
grid_set_ni(GridDesciption &grid, std::string const &value)
{
  grid.ni = parameter_to_int(value);
  grid.nd = 10;
}

static void
grid_set_xcvals(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = grid.xsize;
  if (size == 0) cdo_abort("xsize or gridsize undefined (grid description file: %s)!", dname);
  if (size != nvalues)
    cdo_abort("Number of xcvals=%zu and size of xcvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

  grid.xcvals = new char *[size];
  size_t xstrlen = 64;
  for (size_t i = 0; i < size; ++i) xstrlen = std::max(xstrlen, values[i].size());
  for (size_t i = 0; i < size; ++i) grid.xcvals[i] = new char[xstrlen + 1];
  for (size_t i = 0; i < size; ++i) std::strcpy(grid.xcvals[i], parameter_to_word(values[i].c_str()));
}

static void
grid_set_ycvals(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = grid.ysize;
  if (size == 0) cdo_abort("ysize or gridsize undefined (grid description file: %s)!", dname);
  if (size != nvalues)
    cdo_abort("Number of ycvals=%zu and size of ycvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

  grid.ycvals = new char *[size];
  size_t ystrlen = 64;
  for (size_t i = 0; i < size; ++i) ystrlen = std::max(ystrlen, values[i].size());
  for (size_t i = 0; i < size; ++i) grid.ycvals[i] = new char[ystrlen + 1];
  for (size_t i = 0; i < size; ++i) std::strcpy(grid.ycvals[i], parameter_to_word(values[i].c_str()));
}

static void
grid_set_xvals(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.xsize;
  if (size == 0) cdo_abort("xsize or gridsize undefined (grid description file: %s)!", dname);
  if (size != nvalues)
    cdo_abort("Number of xvals=%zu and size of xvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

  grid.xvals.resize(size);
  for (size_t i = 0; i < size; ++i) grid.xvals[i] = parameter_to_double(values[i]);
}

static void
grid_set_yvals(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.ysize;
  if (size == 0) cdo_abort("ysize or gridsize undefined (grid description file: %s)!", dname);
  if (size != nvalues)
    cdo_abort("Number of yvals=%zu and size of yvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

  grid.yvals.resize(size);
  for (size_t i = 0; i < size; ++i) grid.yvals[i] = parameter_to_double(values[i]);
}

static void
grid_set_xbounds(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.xsize;
  if (size == 0) cdo_abort("xsize or gridsize undefined (grid description file: %s)!", dname);
  if (grid.nvertex == 0) cdo_abort("nvertex undefined (grid description file: %s)!", dname);
  if (grid.nvertex * size != nvalues)
    cdo_abort("Number of xbounds=%zu and size of xbounds=%zu differ (grid description file: %s)!", nvalues, grid.nvertex * size,
              dname);

  grid.xbounds.resize(grid.nvertex * size);
  for (size_t i = 0; i < grid.nvertex * size; ++i) grid.xbounds[i] = parameter_to_double(values[i]);
}

static void
grid_set_ybounds(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.ysize;
  if (size == 0) cdo_abort("ysize or gridsize undefined (grid description file: %s)!", dname);
  if (grid.nvertex == 0) cdo_abort("nvertex undefined (grid description file: %s)!", dname);
  if (grid.nvertex * size != nvalues)
    cdo_abort("Number of ybounds=%zu and size of ybounds=%zu differ (grid description file: %s)!", nvalues, grid.nvertex * size,
              dname);

  grid.ybounds.resize(grid.nvertex * size);
  for (size_t i = 0; i < grid.nvertex * size; ++i) grid.ybounds[i] = parameter_to_double(values[i]);
}

static void
grid_set_gridlatlon(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  if (grid.size == 0) grid.size = grid.xsize * grid.ysize;
  if (grid.size == 0) cdo_abort("gridsize undefined (grid description file: %s)!", dname);
  if (grid.size * 2 != nvalues)
    cdo_abort("Number of gridlonlat values=%zu and size of grid=%zu differ (grid description file: %s)!", nvalues, grid.size * 2,
              dname);
  grid.xvals.resize(grid.size);
  grid.yvals.resize(grid.size);
  for (size_t i = 0; i < grid.size; ++i)
  {
    grid.yvals[i] = parameter_to_double(values[2 * i]);
    grid.xvals[i] = parameter_to_double(values[2 * i + 1]);
  }
}

static void
grid_set_mask(GridDesciption &grid, std::vector<std::string> const &values, const char *dname, size_t nvalues)
{
  auto size = grid.size;
  if (grid.size == 0) cdo_abort("gridsize undefined (grid description file: %s)!", dname);
  if (size != nvalues)
    cdo_abort("Number of mask values=%zu and size of grid=%zu differ (grid description file: %s)!", nvalues, size, dname);
  grid.mask.resize(size);
  size_t count = 0;
  for (size_t i = 0; i < size; ++i)
  {
    grid.mask[i] = parameter_to_int(values[i]);
    if (grid.mask[i] == 1) count++;
  }
  if (count == size)
  {
    grid.mask.clear();
    grid.mask.shrink_to_fit();
  }
}

static void
grid_set_reducedPoints(GridDesciption &grid, std::vector<std::string> const &values, const char *dname)
{
  auto size = grid.ysize;
  if (size == 0) cdo_abort("ysize undefined (grid description file: %s)!", dname);
  grid.reducedPoints.resize(size);
  for (size_t i = 0; i < size; ++i) grid.reducedPoints[i] = parameter_to_int(values[i]);
}

static void
grid_set_indices(GridDesciption &grid, std::vector<std::string> const &values, const char *dname)
{
  auto size = grid.size;
  if (size == 0) cdo_abort("grid size undefined (grid description file: %s)!", dname);
  grid.indices.resize(size);
  for (size_t i = 0; i < size; ++i) grid.indices[i] = parameter_to_long(values[i]);
}

static void
grid_read_data(int ikv, int nkv, std::vector<KVMap> const &kvmap, GridDesciption &grid, int &iproj, int &igmap, const char *dname)
{
  for (int ik = ikv; ik < nkv; ++ik)
  {
    if (!kvmap[ik].isValid) continue;

    const auto kv = kvmap[ik].kv;
    auto const &key = kv->key;
    size_t nvalues = kv->nvalues;
    if (nvalues == 0) continue;
    auto const &values = kv->values;
    auto const &value = kv->values[0];

    // clang-format off
    if      (key == "gridtype")   { grid_set_gridtype(grid, value, dname, ik, iproj); }
    else if (key == "datatype")   { grid_set_datatype(grid, value, dname); }
    else if (key == "gridsize")   { grid.size = parameter_to_size_t(value); }
    else if (key == "xsize")      { grid.xsize = parameter_to_size_t(value); }
    else if (key == "nlon")       { grid.xsize = parameter_to_size_t(value); }
    else if (key == "ysize")      { grid.ysize = parameter_to_size_t(value); }
    else if (key == "nlat")       { grid.ysize = parameter_to_size_t(value); }
    else if (key == "truncation") { grid.ntr = parameter_to_int(value); }
    else if (key == "numLPE")     { grid.numLPE = parameter_to_int(value); }
    else if (key == "np")         { grid.numLPE = parameter_to_int(value); } // np: obsolete
    else if (key == "nvertex")    { grid.nvertex = parameter_to_int(value); }
    else if (key == "complexpacking") { grid.lcomplex = parameter_to_int(value); }
    else if (key == "ni")         { grid_set_ni(grid, value); }
    else if (key == "position")   { grid.position = parameter_to_int(value); }
    else if (key == "number")     { grid.number = parameter_to_int(value); }
    else if (key == "scanningMode") { grid_set_scanningMode(grid, value); }
    else if (key == "xname")      { grid.xname = parameter_to_word(value); }
    else if (key == "yname")      { grid.yname = parameter_to_word(value); }
    else if (key == "xdimname")   { grid.xdimname = parameter_to_word(value); }
    else if (key == "ydimname")   { grid.ydimname = parameter_to_word(value); }
    else if (key == "vdimname")   { grid.vdimname = parameter_to_word(value); }
    else if (key == "xlongname")  { grid.xlongname = value; }
    else if (key == "ylongname")  { grid.ylongname = value; }
    else if (key == "xunits")     { grid.xunits = value; }
    else if (key == "yunits")     { grid.yunits = value; }
    else if (key == "path")       { grid.path = value; }
    else if (key == "uuid")       { cdiStr2UUID(value.c_str(), grid.uuid); }
    else if (key == "xfirst")     { grid.xfirst = parameter_to_double(value); }
    else if (key == "yfirst")     { grid.yfirst = parameter_to_double(value); }
    else if (key == "xlast")      { grid.xlast = parameter_to_double(value); }
    else if (key == "ylast")      { grid.ylast = parameter_to_double(value); }
    else if (key == "xinc")       { grid.xinc = parameter_to_double(value); }
    else if (key == "yinc")       { grid.yinc = parameter_to_double(value); }
    else if (key == "a")          { grid.a = parameter_to_double(value); }
    else if (key == "refinement_level")  { grid.refinementLevel = parameter_to_int(value); }
    else if (key == "indexing_scheme")   { grid.healpixOrder = value; }
    else if (key == "xcvals")     { grid_set_xcvals(grid, values,  dname, nvalues); }
    else if (key == "ycvals")     { grid_set_ycvals(grid, values, dname, nvalues); }
    else if (key == "xvals")      { grid_set_xvals(grid, values, dname, nvalues); }
    else if (key == "yvals")      { grid_set_yvals(grid, values, dname, nvalues); }
    else if (key == "xbounds")    { grid_set_xbounds(grid, values, dname, nvalues); }
    else if (key == "ybounds")    { grid_set_ybounds(grid, values, dname, nvalues); }
    else if (key == "gridlatlon") { grid_set_gridlatlon(grid, values, dname, nvalues); }
    else if (key == "mask")       { grid_set_mask(grid, values, dname, nvalues); }
    else if (key == "reducedPoints") { grid_set_reducedPoints(grid, values, dname); }
    else if (key == "cell_index")    { grid_set_indices(grid, values, dname); }
    else if (key == "grid_mapping_name") { igmap = ik;  break; }
    else if (key == "grid_mapping")      { igmap = ik;  break; }
    else { cdo_abort("Invalid keyword >%s< (grid description file: %s)", key, dname); }
    // clang-format on
  }
}

static void
grid_read_mapping(int igmap, int nkv, std::vector<KVMap> const &kvmap, int gridID)
{
  auto hasGridmapVarname = false;
  for (int ik = igmap; ik < nkv; ++ik)
  {
    if (!kvmap[ik].isValid) continue;

    const auto kv = kvmap[ik].kv;
    auto const &key = kv->key;
    size_t nvalues = kv->nvalues;
    if (nvalues == 0) continue;
    auto const &values = kv->values;
    auto const &value = kv->values[0];

    if (key == "grid_mapping")
    {
      hasGridmapVarname = true;
      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, value.c_str());
      continue;
    }

    if (key == "grid_mapping_name")
    {
      if (!hasGridmapVarname)
      {
        cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME,
                        (value == "rotated_latitude_longitude") ? "rotated_pole" : "crs");
      }

      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, value.c_str());
    }

    auto dtype = literals_find_datatype(nvalues, values);
    if (dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32)
    {
      std::vector<int> ivals(nvalues);
      for (size_t i = 0; i < nvalues; ++i) ivals[i] = literal_to_int(values[i]);
      cdiDefAttInt(gridID, CDI_GLOBAL, key.c_str(), dtype, nvalues, ivals.data());
    }
    else if (dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64)
    {
      std::vector<double> dvals(nvalues);
      for (size_t i = 0; i < nvalues; ++i) dvals[i] = literal_to_double(values[i]);
      cdiDefAttFlt(gridID, CDI_GLOBAL, key.c_str(), dtype, nvalues, dvals.data());
    }
    else
    {
      auto len = (int) value.size();
      cdiDefAttTxt(gridID, CDI_GLOBAL, key.c_str(), len, value.c_str());
    }
  }
}

int
grid_read(std::FILE *gfp, const char *dname)
{
  PMList pmlist;
  pmlist.read_namelist(gfp, dname);
  if (pmlist.size() == 0) return -1;
  KVList &kvlist = pmlist.front();

  auto nkv = static_cast<int>(kvlist.size());
  if (nkv == 0) return -1;

  std::vector<KVMap> kvmap(nkv);
  for (int i = 0; i < nkv; ++i) kvmap[i].isValid = false;

  int ik = 0;
  const std::string firstKey = "gridtype";
  for (auto &kv : kvlist)
  {
    if (ik == 0 && kv.key != firstKey) cdo_abort("First grid description keyword must be >%s< (found: %s)!", firstKey, kv.key);

    if (kv.nvalues == 0) { cdo_warning("Grid description keyword %s has no values, skipped!", kv.key); }
    else
    {
      kvmap[ik].isValid = true;
      kvmap[ik].kv = &kv;
    }
    ik++;
  }

  int iproj = 0;
  int igmap = 0;
  GridDesciption grid;
  grid_read_data(0, nkv, kvmap, grid, iproj, igmap, dname);

  auto gridID = (grid.type == CDI_UNDEFID) ? CDI_UNDEFID : grid_define(grid);

  if (gridID != CDI_UNDEFID)
  {
    auto gridprojID = gridID;

    if (iproj > 0)
    {
      GridDesciption proj;
      grid_read_data(iproj, nkv, kvmap, proj, iproj, igmap, dname);

      auto projID = (proj.type == CDI_UNDEFID) ? CDI_UNDEFID : grid_define(proj);
      if (projID != CDI_UNDEFID)
      {
        gridDefProj(gridID, projID);
        gridprojID = projID;
      }
    }

    if (igmap > 0) grid_read_mapping(igmap, nkv, kvmap, gridprojID);
  }

  return gridID;
}
