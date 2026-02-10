/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cassert>

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "printinfo.h"
#include "gaussian_latitudes.h"

extern "C"
{
#include "lib/gradsdes/gradsdes.h"
}

static void
get_dim_vals(dsets_t *pfi, double *vals, int dimlen, int dim)
{
  gadouble (*conv)(gadouble *, gadouble);

  assert(dimlen == pfi->dnum[dim]);

  if (pfi->linear[dim] == 0)
  {
    for (int i = 0; i < dimlen; ++i) vals[i] = pfi->grvals[dim][i + 1];
  }
  else if (pfi->linear[dim] == 1)
  {
    conv = pfi->gr2ab[dim];
    gadouble *cvals = pfi->grvals[dim];
    for (int i = 0; i < dimlen; ++i) vals[i] = conv(cvals, i + 1);
  }
}

static void
rev_vals(double *vals, int n)
{
  for (int i = 0; i < n / 2; ++i)
  {
    double dum = vals[i];
    vals[i] = vals[n - 1 - i];
    vals[n - 1 - i] = dum;
  }
}

static int
define_grid(dsets_t *pfi)
{
  int nx = pfi->dnum[0];
  int ny = pfi->dnum[1];

  Varray<double> xvals(nx);
  Varray<double> yvals(ny);

  get_dim_vals(pfi, xvals.data(), nx, 0);
  get_dim_vals(pfi, yvals.data(), ny, 1);

  if (pfi->yrflg) rev_vals(yvals.data(), ny);

  bool isGaussLat = (pfi->linear[1] == 0) ? is_gaussian_latitudes((size_t) ny, yvals.data()) : false;

  int gridtype = isGaussLat ? GRID_GAUSSIAN : GRID_LONLAT;

  int gridID = gridCreate(gridtype, nx * ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);

  gridDefXvals(gridID, xvals.data());
  gridDefYvals(gridID, yvals.data());

  return gridID;
}

static int
define_level(dsets_t *pfi, int nlev)
{
  int zaxisID = -1;
  int nz = pfi->dnum[2];

  if (nz)
  {
    Varray<double> zvals(nz);

    get_dim_vals(pfi, zvals.data(), nz, 2);

    if (nz == 1 && is_equal(zvals[0], 0))
      zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
    else
    {
      if (nlev > 0 && nlev < nz) nz = nlev;
      if (pfi->zrflg) rev_vals(zvals.data(), nz);
      zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
    }
    zaxisDefLevels(zaxisID, zvals.data());
  }
  else
  {
    double level = 0;
    nz = 1;
    zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
    zaxisDefLevels(zaxisID, &level);
  }

  return zaxisID;
}

class Importbinary : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Importbinary",
    .operators = { { "import_binary", ImportbinaryHelp } },
    .aliases = { { "import_grads", "import_binary" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Importbinary> registration = RegisterEntry<Importbinary>(module);

  size_t numMissVals = 0;
  int told{}, fnum{};
  int tmin = 0, tmax = 0;
  char *ch = nullptr;
  int flag{};
  struct dt dtim{}, dtimi{};
  double sfclevel = 0;

  int taxisID{};
  CdoStreamID streamID{};

  int vlistID{};
  int gridID{};
  int zaxisID{};

  size_t gridsize{};
  dsets_t pfi{};

  int numFields{};
  size_t recsize{};
  Varray<char> rec{};
  Varray<double> array{};

  int recoffset{};
  Varray<int> var_dfrm{};

  std::vector<FieldInfo> fieldInfoList{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    dsets_init(&pfi);

    int status = read_gradsdes((char *) cdo_get_stream_name(0), &pfi);
    if (Options::cdoVerbose) fprintf(stderr, "status %d\n", status);
    // if (status) cdo_abort("Open failed on %s!", pfi.name);
    if (status) cdo_abort("Open failed!");

    numFields = pfi.trecs;
    int nvars = pfi.vnum;
    struct gavar *pvar = pfi.pvar1;

    if (nvars == 0) cdo_abort("No variables found!");

    gridID = define_grid(&pfi);
    zaxisID = define_level(&pfi, 0);

    auto zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
    zaxisDefLevels(zaxisIDsfc, &sfclevel);

    vlistID = vlistCreate();

    Varray<int> var_zaxisID(nvars);
    var_dfrm = Varray<int>(numFields);
    fieldInfoList = std::vector<FieldInfo>(numFields);

    int fieldID = 0;
    for (int ivar = 0; ivar < nvars; ++ivar)
    {
      int varID = -1;
      /*
      if ( Options::cdoVerbose )
        fprintf(stderr, "1:%s 2:%s %d %d %d %d 3:%s %d \n",
                pvar->abbrv, pvar->longnm, pvar->offset, pvar->recoff,
      pvar->levels, pvar->nvardims, pvar->varnm, pvar->var_t);
      */
      int nlevels = pvar->levels;

      if (nlevels == 0)
      {
        nlevels = 1;
        varID = vlistDefVar(vlistID, gridID, zaxisIDsfc, TIME_VARYING);
      }
      else
      {
        if (nlevels > zaxisInqSize(zaxisID))
          cdo_abort("Variable %s has too many number of levels!", pvar->abbrv);
        else if (nlevels < zaxisInqSize(zaxisID))
        {
          int vid, zid = -1, nlev;
          for (vid = 0; vid < ivar; ++vid)
          {
            zid = var_zaxisID[vid];
            nlev = zaxisInqSize(zid);
            if (nlev == nlevels) break;
          }

          if (vid == ivar) zid = define_level(&pfi, nlevels);
          varID = vlistDefVar(vlistID, gridID, zid, TIME_VARYING);
        }
        else
          varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      }

      var_zaxisID[varID] = vlistInqVarZaxis(vlistID, varID);

      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, pvar->abbrv);
      {
        char *longname = pvar->varnm;
        int len = (int) std::strlen(longname);
        if (longname[0] == '\'' && longname[len - 1] == '\'')
        {
          longname[len - 1] = 0;
          longname++;
        }
        if (longname[0] == '\t') longname++;
        cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname);
      }

      double missval = pfi.undef;
      int datatype = CDI_DATATYPE_FLT32;

      if (pvar->dfrm == 1)
      {
        datatype = CDI_DATATYPE_UINT8;
        if (missval < 0 || missval > 255) missval = 255;
      }
      else if (pvar->dfrm == 2)
      {
        datatype = CDI_DATATYPE_UINT16;
        if (missval < 0 || missval > 65535) missval = 65535;
      }
      else if (pvar->dfrm == -2)
      {
        datatype = CDI_DATATYPE_INT16;
        if (missval < -32768 || missval > 32767) missval = -32768;
      }
      else if (pvar->dfrm == 4)
      {
        datatype = CDI_DATATYPE_INT32;
        if (missval < -2147483648 || missval > 2147483647) missval = -2147483648;
      }
      else if (pfi.flt64)
        datatype = CDI_DATATYPE_FLT64;

      vlistDefVarDatatype(vlistID, varID, datatype);
      vlistDefVarMissval(vlistID, varID, missval);

      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        if (fieldID >= numFields) cdo_abort("Internal problem with number of records!");
        fieldInfoList[fieldID].set(varID, levelID);
        var_dfrm[fieldID] = pvar->dfrm;
        fieldID++;
      }

      pvar++;
    }

    gr2t(pfi.grvals[3], (gadouble) 1, &dtim);
    CdiDateTime rDateTime{};
    rDateTime.date = cdiDate_encode(dtim.yr, dtim.mo, dtim.dy);
    rDateTime.time = cdiTime_encode(dtim.hr, dtim.mn, 0, 0);

    auto calendar = CALENDAR_STANDARD;
    taxisID = cdo_taxis_create(TAXIS_RELATIVE);
    taxisDefCalendar(taxisID, calendar);
    taxisDefRdatetime(taxisID, rDateTime);
    vlistDefTaxis(vlistID, taxisID);

    streamID = cdo_open_write(1);
    cdo_def_vlist(streamID, vlistID);

    gridsize = pfi.dnum[0] * pfi.dnum[1];
    recoffset = pfi.xyhdr * (pfi.flt64 ? 8 : 4);
    if (pfi.seqflg) recoffset += 4;

    // recsize = pfi.gsiz*4;
    recsize = pfi.gsiz * 8;

    rec = Varray<char>(recsize);
    array = Varray<double>(gridsize);
  }
  void
  run() override
  {
    pfi.infile = nullptr;
    int tcur = 0;
    int e = 1;
    while (1)
    {  // loop over all times for this ensemble
      if (pfi.tmplat)
      {
        // make sure no file is open
        if (pfi.infile != nullptr)
        {
          std::fclose(pfi.infile);
          pfi.infile = nullptr;
        }
        // advance to first valid time step for this ensemble
        if (tcur == 0)
        {
          told = 0;
          tcur = 1;
          while (pfi.fnums[tcur - 1] == -1) tcur++;
        }
        else
        {  // tcur!=0
          told = pfi.fnums[tcur - 1];
          // increment time step until fnums changes
          while (told == pfi.fnums[tcur - 1] && tcur <= pfi.dnum[3])
          {
            tcur++;
            if (tcur > pfi.dnum[3]) break;
          }
        }

        // make sure we haven't advanced past end of time axis
        if (tcur > pfi.dnum[3]) break;

        // check if we're past all valid time steps for this ensemble
        if ((told != -1) && (pfi.fnums[tcur - 1] == -1)) break;

        /* Find the range of t indexes that have the same fnums value.
           These are the times that are contained in this particular file */
        tmin = tcur;
        tmax = tcur - 1;
        fnum = pfi.fnums[tcur - 1];
        if (fnum != -1)
        {
          while (fnum == pfi.fnums[tmax])
          {
            tmax++;
            if (tmax == pfi.dnum[3]) break;
          }
          gr2t(pfi.grvals[3], (gadouble) tcur, &dtim);
          gr2t(pfi.grvals[3], (gadouble) 1, &dtimi);
          ch = gafndt(pfi.name, &dtim, &dtimi, pfi.abvals[3], pfi.pchsub1, nullptr, tcur, e, &flag);
          if (ch == nullptr) cdo_abort("Couldn't determine data file name for e=%d t=%d!", e, tcur);
        }
      }
      else
      {
        // Data set is not templated
        ch = pfi.name;
        tmin = 1;
        tmax = pfi.dnum[3];
      }

      // Open this file and position to start of first record
      if (Options::cdoVerbose) cdo_print("Opening file: %s", ch);
      pfi.infile = std::fopen(ch, "rb");
      if (pfi.infile == nullptr)
      {
        if (pfi.tmplat)
        {
          cdo_warning("Could not open file: %s", ch);
          break;
        }
        else { cdo_abort("Could not open file: %s", ch); }
      }

      // file header
      if (pfi.fhdr > 0) fseeko(pfi.infile, pfi.fhdr, SEEK_SET);

      // Get file size
      /*
      fseeko(pfi.infile,0L,2);
      flen = ftello(pfi.infile);

      printf("flen %d tsiz %d\n", flen, pfi.tsiz);

      fseeko (pfi.infile,0,0);
      */
      for (int tsID = tmin - 1; tsID < tmax; ++tsID)
      {
        gr2t(pfi.grvals[3], (gadouble) (tsID + 1), &dtim);
        CdiDateTime vDateTime{};
        vDateTime.date = cdiDate_encode(dtim.yr, dtim.mo, dtim.dy);
        vDateTime.time = cdiTime_encode(dtim.hr, dtim.mn, 0, 0);
        if (Options::cdoVerbose) cdo_print(" Reading timestep: %3d %s", tsID + 1, datetime_to_string(vDateTime));

        taxisDefVdatetime(taxisID, vDateTime);
        cdo_def_timestep(streamID, tsID);

        for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          // record size depends on data type
          if (var_dfrm[fieldID] == 1) { recsize = pfi.gsiz; }
          else if ((var_dfrm[fieldID] == 2) || (var_dfrm[fieldID] == -2)) { recsize = pfi.gsiz * 2; }
          else { recsize = pfi.flt64 ? pfi.gsiz * 8 : pfi.gsiz * 4; }

          size_t rc = std::fread(rec.data(), 1, recsize, pfi.infile);
          if (rc < recsize) cdo_abort("I/O error reading record=%d of timestep=%d!", fieldID + 1, tsID + 1);

          char *cdata = &rec[recoffset];

          // convert
          if (var_dfrm[fieldID] == 1)
          {
            const unsigned char *carray = (const unsigned char *) cdata;
            for (size_t i = 0; i < gridsize; ++i) array[i] = (double) carray[i];
          }
          else if (var_dfrm[fieldID] == 2)
          {
            if (pfi.bswap) gabswp2(cdata, gridsize);
            const unsigned short *sarray = (const unsigned short *) cdata;
            for (size_t i = 0; i < gridsize; ++i) array[i] = (double) sarray[i];
          }
          else if (var_dfrm[fieldID] == -2)
          {
            if (pfi.bswap) gabswp2(cdata, gridsize);
            const short *sarray = (const short *) cdata;
            for (size_t i = 0; i < gridsize; ++i) array[i] = (double) sarray[i];
          }
          else if (var_dfrm[fieldID] == 4)
          {
            if (pfi.bswap) gabswp(cdata, gridsize);
            const int *iarray = (const int *) cdata;
            for (size_t i = 0; i < gridsize; ++i) array[i] = (double) iarray[i];
          }
          else
          {
            if (pfi.flt64)
            {
              if (pfi.bswap) cdo_abort("Byte swap not implemented for 64-bit floats!");
              const double *darray = (const double *) cdata;
              for (size_t i = 0; i < gridsize; ++i) array[i] = darray[i];
            }
            else
            {
              if (pfi.bswap) gabswp(cdata, gridsize);
              const float *farray = (const float *) cdata;
              for (size_t i = 0; i < gridsize; ++i) array[i] = (double) farray[i];
            }
          }

          double fmin = 1.e99;
          double fmax = -1.e99;
          numMissVals = 0;
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (array[i] > pfi.ulow && array[i] < pfi.uhi)
            {
              array[i] = pfi.undef;
              numMissVals++;
            }
            else if (std::isnan(array[i]))
            {
              array[i] = pfi.undef;
              numMissVals++;
            }
            else
            {
              fmin = std::min(fmin, array[i]);
              fmax = std::max(fmax, array[i]);
            }
          }

          auto [varID, levelID] = fieldInfoList[fieldID].get();
          cdo_def_field(streamID, varID, levelID);
          cdo_write_field(streamID, array.data(), numMissVals);
        }
      }

      // break out if not templating
      if (!pfi.tmplat) break;

    }  // end of while (1) loop
  }

  void
  close() override
  {
    process_def_var_num(vlistNvars(vlistID));

    cdo_stream_close(streamID);
    vlistDestroy(vlistID);

    if (pfi.infile) std::fclose(pfi.infile);
  }
};
