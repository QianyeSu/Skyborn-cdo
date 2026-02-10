/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <time.h>

#include <cdi.h>

#include "process_int.h"
#include "griddes.h"
#include "commandline.h"
#include "cdo_default_values.h"

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

static constexpr int nvars_ml = 4;
#ifdef HAVE_LIBNETCDF
static const char strfiletype_ml[] = "Initial file spectral";
#endif

struct VAR
{
  int gridtype;
  int zaxistype;
  int code;
  std::string name;
  std::string longname;
  std::string units;
  int gridID;
  int zaxisID;
  size_t gridsize;
  int nlev;
  double *ptr;
};

struct ATTS
{
  int naint;
  int naflt;
  int natxt;
  char *atxtname[1024];
  char *atxtentry[1024];
};

static void
iniatts(ATTS *atts)
{
  atts->naint = 0;
  atts->naflt = 0;
  atts->natxt = 0;
}

static void
inivar(VAR &var, int gridtype, int zaxistype, int code, std::string const &name, std::string const &longname,
       std::string const &units)
{
  var.gridtype = gridtype;
  var.zaxistype = zaxistype;
  var.code = code;
  var.name = name;
  var.longname = longname;
  var.units = units;
}

#ifdef HAVE_LIBNETCDF
static void
inivars_ml(std::vector<VAR> &vars)
{
  vars.resize(nvars_ml + 1);

  inivar(vars[0], GRID_GAUSSIAN, ZAXIS_HYBRID, 133, "Q", "specific humidity", "kg/kg");
  inivar(vars[1], GRID_SPECTRAL, ZAXIS_HYBRID, 138, "SVO", "vorticity", "1/s");
  inivar(vars[2], GRID_SPECTRAL, ZAXIS_HYBRID, 155, "SD", "divergence", "1/s");
  inivar(vars[3], GRID_SPECTRAL, ZAXIS_HYBRID, 130, "STP", "temperature", "K");
  // Don't change the order (lsp must be the last one)!
  inivar(vars[4], GRID_SPECTRAL, ZAXIS_SURFACE, 152, "LSP", "log surface pressure", "");
}

static void
nce(int istat)
{
  // This routine provides a simple interface to NetCDF error message routine.

  if (istat != NC_NOERR) cdo_abort(nc_strerror(istat));
}
#endif

static int
import_e5ml(const char *filename, std::vector<VAR> &vars)
{
  int nvars = 0;
#ifdef HAVE_LIBNETCDF
  // open file and check file type
  auto nc_file_id = cdo_cdf_openread(filename);

  char filetype[256];
  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  size_t attlen;
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if (std::strcmp(filetype, strfiletype_ml) != 0) return 0;

  inivars_ml(vars);

  // read dimensions

  int nc_dim_id;
  nce(nc_inq_dimid(nc_file_id, "lon", &nc_dim_id));
  size_t dimlen;
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  auto nlon = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "lat", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  auto nlat = (int) dimlen;

  auto gridIDgp = gridCreate(GRID_GAUSSIAN, nlon * nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  nce(nc_inq_dimid(nc_file_id, "nsp", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  auto nsp = (int) dimlen;

  auto gridIDsp = gridCreate(GRID_SPECTRAL, nsp * 2);
  gridDefComplexPacking(gridIDsp, 1);

  nce(nc_inq_dimid(nc_file_id, "nlev", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  auto nlev = (int) dimlen;
  auto nlevp1 = nlev + 1;
  auto nvct = nlevp1 * 2;

  auto zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  auto zaxisIDml = zaxisCreate(ZAXIS_HYBRID, nlev);

  {
    Varray<double> levs(nlev);
    for (int i = 0; i < nlev; ++i) levs[i] = i + 1;
    zaxisDefLevels(zaxisIDml, levs.data());
  }

  // read variables

  int nc_var_id;
  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));

  {
    auto xvals = new double[nlon];
    auto yvals = new double[nlat];

    nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

    nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
    nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

    gridDefXvals(gridIDgp, xvals);
    gridDefYvals(gridIDgp, yvals);

    delete[] xvals;
    delete[] yvals;
  }

  auto vct = new double[nvct];

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct + nlevp1));

  zaxisDefVct(zaxisIDml, 2 * nlevp1, vct);
  delete[] vct;

  size_t start[3], count[3];
  for (int iv = 0; iv < nvars_ml; iv++)
  {
    size_t nvals = 0;

    auto gridtype = vars[iv].gridtype;

    if (gridtype == GRID_GAUSSIAN)
    {
      vars[iv].gridID = gridIDgp;
      nvals += nlon * nlat;
    }
    else
    {
      vars[iv].gridID = gridIDsp;
      nvals += nsp * 2;
    }

    vars[iv].zaxisID = zaxisIDml;
    vars[iv].gridsize = nvals;
    vars[iv].nlev = nlev;

    vars[iv].ptr = new double[nlev * nvals];

    for (int i = 0; i < nlev; ++i)
    {
      if (gridtype == GRID_GAUSSIAN)
      {
        start[0] = 0;
        start[1] = i;
        start[2] = 0;
        count[0] = nlat;
        count[1] = 1;
        count[2] = nlon;
      }
      else
      {
        start[0] = 0;
        start[1] = 0;
        start[2] = i;
        count[0] = nsp;
        count[1] = 2;
        count[2] = 1;
      }

      nce(nc_inq_varid(nc_file_id, vars[iv].name.c_str(), &nc_var_id));
      nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, vars[iv].ptr + i * nvals));
    }
  }

  // read lsp

  vars[nvars_ml].gridID = gridIDsp;
  vars[nvars_ml].zaxisID = zaxisIDsfc;
  vars[nvars_ml].gridsize = nsp * 2;
  vars[nvars_ml].nlev = 1;

  start[0] = 0;
  start[1] = 0;
  start[2] = nlev;
  count[0] = nsp;
  count[1] = 2;
  count[2] = 1;

  vars[nvars_ml].ptr = new double[nsp * 2];

  nce(nc_inq_varid(nc_file_id, "STP", &nc_var_id));
  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, vars[nvars_ml].ptr));

  // close input file
  cdo_cdf_close(nc_file_id);

  nvars = nvars_ml + 1;

#else
  cdo_abort("NetCDF support not compiled in!");
#endif

  return nvars;
}

static void
export_e5ml(const char *filename, const std::vector<VAR> &vars, int nvars, int vdate, int vtime, int ntr)
{
#ifdef HAVE_LIBNETCDF

  auto date_and_time_in_sec = time(nullptr);

  char timestr[30];
  timestr[0] = 0;

  if (date_and_time_in_sec != -1)
  {
    auto date_and_time = localtime(&date_and_time_in_sec);
    (void) strftime(timestr, sizeof(timestr), "%d/%m/%Y %H:%M", date_and_time);
  }

  const char *username = getenv("LOGNAME");
  if (username == nullptr)
  {
    username = getenv("USER");
    if (username == nullptr) username = "unknown";
  }

  int n2 = 2;
  int lon = 0;
  int lat = 0;
  int nsp = 0;
  int nlev = 0;
  int nlevp1 = 0;
  int nvclev = 0;
  int gridIDgp = -1, zaxisIDml = -1;
  for (int varid = 0; varid < nvars; ++varid)
  {
    auto gridtype = vars[varid].gridtype;
    auto zaxistype = vars[varid].zaxistype;

    if (gridtype == GRID_GAUSSIAN && lat == 0)
    {
      gridIDgp = vars[varid].gridID;
      lon = gridInqXsize(vars[varid].gridID);
      lat = gridInqYsize(vars[varid].gridID);
    }
    else if (gridtype == GRID_SPECTRAL && nsp == 0)
    {
      nsp = gridInqSize(vars[varid].gridID);
      nsp = nsp / 2;
    }

    if (zaxistype == ZAXIS_HYBRID && nlev == 0)
    {
      zaxisIDml = vars[varid].zaxisID;
      nlev = zaxisInqSize(vars[varid].zaxisID);
      nlevp1 = nlev + 1;
      nvclev = nlev + 1;
    }
  }

  if (lat == 0) cdo_abort("Gaussian grid not found!");
  if (nsp == 0) cdo_abort("Spectral data not found!");
  if (nlev == 0) cdo_abort("Hybrid level not found!");

  size_t nlon = lon;
  size_t nlat = lat;

  size_t data_size = nlon + nlat + 2 * nvclev + 2 * nsp * 2 * nlev + nsp * 2 * nlevp1 + nlon * nlat * nlev;

  int writemode = NC_CLOBBER;
  if (data_size * 8 > 2147000000)
  {
#if defined(NC_64BIT_OFFSET)
    writemode = NC_CLOBBER | NC_64BIT_OFFSET;
#else
    cdoWarning("Datasize > 2GB and NC_64BIT_OFFSET not available!");
#endif
  }

  // create file
  int nc_file_id;
  nce(nc_create(filename, writemode, &nc_file_id));

  char atttext[1024];
  std::strcpy(atttext, "IEEE");
  size_t attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "source_type", attlen, atttext));

  std::strcpy(atttext, cdo::command_line());
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "history", attlen, atttext));

  std::strcpy(atttext, username);
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "user", attlen, atttext));

  std::strcpy(atttext, timestr);
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "created", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_1", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_2", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_3", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_4", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_5", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_6", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_7", attlen, atttext));

  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_8", attlen, atttext));

  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "fdate", NC_INT, 1, &vdate));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "ftime", NC_INT, 1, &vtime));

  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "vdate", NC_INT, 1, &vdate));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "vtime", NC_INT, 1, &vtime));

  // attint = 31;
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_n", NC_INT, 1, &ntr));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_m", NC_INT, 1, &ntr));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_k", NC_INT, 1, &ntr));

  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "file_type", std::strlen(strfiletype_ml), strfiletype_ml));

  int lat_dimid;
  nce(nc_def_dim(nc_file_id, "lat", lat, &lat_dimid));

  int lon_dimid;
  nce(nc_def_dim(nc_file_id, "lon", lon, &lon_dimid));

  int nlev_dimid, nlevp1_dimid;
  nce(nc_def_dim(nc_file_id, "nlev", nlev, &nlev_dimid));
  nce(nc_def_dim(nc_file_id, "nlevp1", nlevp1, &nlevp1_dimid));

  int nsp_dimid;
  nce(nc_def_dim(nc_file_id, "nsp", nsp, &nsp_dimid));

  int nvclev_dimid;
  nce(nc_def_dim(nc_file_id, "nvclev", nvclev, &nvclev_dimid));

  int n2_dimid;
  nce(nc_def_dim(nc_file_id, "n2", n2, &n2_dimid));

  nce(nc_enddef(nc_file_id));

  // define gaussian grid

  auto xvals = new double[nlon];
  auto yvals = new double[nlat];

  gridInqXvals(gridIDgp, xvals);
  gridInqYvals(gridIDgp, yvals);

  nce(nc_redef(nc_file_id));
  int nc_var_id;
  nce(nc_def_var(nc_file_id, "lat", NC_DOUBLE, 1, &lat_dimid, &nc_var_id));
  std::strcpy(atttext, "Gaussian latitude");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  std::strcpy(atttext, "degrees_N");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, yvals));

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "lon", NC_DOUBLE, 1, &lon_dimid, &nc_var_id));
  std::strcpy(atttext, "longitude");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  std::strcpy(atttext, "degrees_E");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, xvals));

  delete[] xvals;
  delete[] yvals;

  // define model level

  // nvct = nvclev*2;

  auto vct = zaxisInqVctPtr(zaxisIDml);

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_a", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  std::strcpy(atttext, "vertical-coordinate parameter set A");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_b", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  std::strcpy(atttext, "vertical-coordinate parameter set B");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  std::strcpy(atttext, "");
  attlen = std::strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct + nlevp1));

  delete[] vct;

  int lspid = -1;
  int nc_stpid = -1;

  size_t start[3], count[3];
  for (int varid = 0; varid < nvars; varid++)
  {
    size_t nvals = 0;

    auto code = vars[varid].code;
    auto gridtype = vars[varid].gridtype;

    auto ilev = zaxisInqSize(vars[varid].zaxisID);

    if (ilev == 1)
    {
      if (code == 152)
      {
        lspid = varid;
        if (gridtype != GRID_SPECTRAL) cdo_abort("%s has wrong gridtype!", vars[varid].name);
      }
      continue;
    }

    if (nlev != ilev) cdo_abort("Unexpected number of level %d!", ilev);

    int dimidsp[9];
    if (gridtype == GRID_GAUSSIAN)
    {
      nvals = nlon * nlat;

      dimidsp[0] = lat_dimid;
      dimidsp[1] = nlev_dimid;
      dimidsp[2] = lon_dimid;
    }
    else if (gridtype == GRID_SPECTRAL)
    {
      nvals = nsp * 2;

      dimidsp[0] = nsp_dimid;
      dimidsp[1] = n2_dimid;

      if (vars[varid].name == "STP" || vars[varid].name == "T")
        dimidsp[2] = nlevp1_dimid;
      else
        dimidsp[2] = nlev_dimid;
    }
    else
      cdo_abort("Unsupported grid!");

    nce(nc_redef(nc_file_id));
    nce(nc_def_var(nc_file_id, vars[varid].name.c_str(), NC_DOUBLE, 3, dimidsp, &nc_var_id));
    if (vars[varid].longname.size())
      nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", vars[varid].longname.size(), vars[varid].longname.c_str()));
    if (vars[varid].units.size())
      nce(nc_put_att_text(nc_file_id, nc_var_id, "units", vars[varid].units.size(), vars[varid].units.c_str()));
    nce(nc_enddef(nc_file_id));

    if (dimidsp[2] == nlevp1_dimid) nc_stpid = nc_var_id;

    for (int i = 0; i < nlev; ++i)
    {
      if (gridtype == GRID_GAUSSIAN)
      {
        start[0] = 0;
        start[1] = i;
        start[2] = 0;
        count[0] = nlat;
        count[1] = 1;
        count[2] = nlon;
      }
      else
      {
        start[0] = 0;
        start[1] = 0;
        start[2] = i;
        count[0] = nsp;
        count[1] = 2;
        count[2] = 1;
      }

      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr + i * nvals));
    }
  }

  if (lspid == -1) cdo_abort("LSP not found!");
  if (nc_stpid == -1) cdo_abort("STP not found!");

  // write lsp
  start[0] = 0;
  start[1] = 0;
  start[2] = nlev;
  count[0] = nsp;
  count[1] = 2;
  count[2] = 1;

  nce(nc_put_vara_double(nc_file_id, nc_stpid, start, count, vars[lspid].ptr));

  // close input file
  nce(nc_close(nc_file_id));

#else
  cdo_abort("NetCDF support not compiled in!");
#endif
}

class Echam5ini : public Process
{
  using Process::Process;

public:
  void
  ex_e5ml()
  {
    auto streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto taxisID = vlistInqTaxis(vlistID1);

    VarList varList1(vlistID1);
    auto numVars = varList1.numVars();

    std::vector<VAR> vars(numVars);

    int ntr = 0;
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      auto code = var1.code;
      auto name = var1.name;
      auto longname = var1.longname;
      auto units = var1.units;

      if (code < 0) code = 0;
      if (name.substr(0, 3) == "var")
      {
        if (code > 0)
        {
          if (code == 133)
          {
            name = "Q";
            longname = "specific humidity";
            units = "kg/kg";
          }
          if (code == 138)
          {
            name = "SVO";
            longname = "vorticity";
            units = "1/s";
          }
          if (code == 155)
          {
            name = "SD";
            longname = "divergence";
            units = "1/s";
          }
          if (code == 130)
          {
            name = "STP";
            longname = "temperature";
            units = "K";
          }
          if (code == 152)
          {
            name = "LSP";
            longname = "log surface pressure";
          }
        }
      }
      else if (name.substr(0, 3) == "LSP")
        code = 152;

      if (var1.zaxisType == GRID_SPECTRAL && ntr == 0) ntr = gridInqTrunc(var1.gridID);

      auto zaxisType = var1.zaxisType;
      if (zaxisType == ZAXIS_HYBRID && var1.nlevels == 1) zaxisType = ZAXIS_SURFACE;

      inivar(vars[varID], var1.gridType, zaxisType, code, name, longname, units);

      vars[varID].gridID = var1.gridID;
      vars[varID].zaxisID = var1.zaxisID;
      vars[varID].gridsize = var1.gridsize;
      vars[varID].nlev = var1.nlevels;
      vars[varID].ptr = new double[var1.nlevels * var1.gridsize];
    }

    auto numFields = cdo_stream_inq_timestep(streamID1, 0);
    auto vDateTime = taxisInqVdatetime(taxisID);

    auto vdate = cdiDate_get(vDateTime.date);
    auto vtime = cdiTime_get(vDateTime.time);
    if (vdate == 0)
    {
      vdate = 19890101;
      vtime = 120000;
    }

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);

      auto gridsize = vars[varID].gridsize;
      size_t numMissVals;
      cdo_read_field(streamID1, vars[varID].ptr + levelID * gridsize, &numMissVals);
    }

    cdo_stream_close(streamID1);

    export_e5ml(cdo_get_stream_name(1), vars, numVars, vdate, vtime, ntr);
  }

  void
  im_e5ml()
  {
    ATTS atts;
    iniatts(&atts);

    std::vector<VAR> vars;
    auto numVars = import_e5ml(cdo_get_stream_name(0), vars);
    if (numVars == 0) cdo_abort("Unsupported file type!");

    auto vlistID2 = vlistCreate();
    vlistDefNtsteps(vlistID2, 0);

    for (int iv = 0; iv < numVars; iv++)
    {
      int varID = vlistDefVar(vlistID2, vars[iv].gridID, vars[iv].zaxisID, TIME_CONSTANT);
      if (vars[iv].code > 0) vlistDefVarCode(vlistID2, varID, vars[iv].code);
      if (vars[iv].name.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, vars[iv].name.c_str());
      if (vars[iv].longname.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, vars[iv].longname.c_str());
      if (vars[iv].units.size()) cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, vars[iv].units.c_str());
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }

    for (int iatt = 0; iatt < atts.natxt; ++iatt)
    {
      // printf("%s: %s\n", atts.atxtname[iatt], atts.atxtentry[iatt]);
      cdiDefAttTxt(vlistID2, CDI_GLOBAL, atts.atxtname[iatt], (int) std::strlen(atts.atxtentry[iatt]) + 1, atts.atxtentry[iatt]);
    }

    auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID);

    if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = CDI_FILETYPE_NC;

    auto streamID2 = cdo_open_write(1);

    cdo_def_vlist(streamID2, vlistID2);

    int tsID = 0;
    cdo_def_timestep(streamID2, tsID);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = vars[varID].gridsize;
      auto nlev = vars[varID].nlev;

      for (int levelID = 0; levelID < nlev; ++levelID)
      {
        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, vars[varID].ptr + levelID * gridsize, 0);
      }
    }

    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }

  void
  init() override
  {
  }
  void
  close() override
  {
    // vlistDestroy(vlistID2);
  }
};

class Echam5ini_import : public Echam5ini
{
public:
  using Echam5ini::Echam5ini;
  inline static CdoModule module = {
    .name = "Echam5ini",
    .operators = { { "import_e5ml" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Echam5ini_import> registration = RegisterEntry<Echam5ini_import>(module);

public:
  void
  run() override
  {
    im_e5ml();
  }
};

class Echam5ini_export : public Echam5ini
{
public:
  using Echam5ini::Echam5ini;
  inline static CdoModule module = {
    .name = "Echam5ini",
    .operators = { { "import_e5ml" }, { "export_e5ml" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, OnlyFirst },
  };
  inline static RegisterEntry<Echam5ini_export> registration = RegisterEntry<Echam5ini_export>(module);

public:
  void
  run() override
  {
    ex_e5ml();
  }
};
