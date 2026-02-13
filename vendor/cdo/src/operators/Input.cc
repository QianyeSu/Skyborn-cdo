/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Input     input          ASCII input
      Input     inputsrv       SERVICE input
      Input     inputext       EXTRA input
*/

#include <cdi.h>

#include <mpim_grid.h>
#include "process_int.h"
#include "griddes.h"
#include "cdo_zaxis.h"
#include "grid_read_pingo.h"
#include "cdo_default_values.h"

static size_t
input_iarray(size_t numVals, int *array)
{
  size_t readVals = 0;

  for (size_t i = 0; i < numVals; ++i)
  {
    if (scanf("%d", &array[i]) != 1) break;
    readVals++;
  }

  return readVals;
}

static int
read_record(int numFields, int numLevels, size_t &gridsize0, int &code, int &date, int &time, Varray<double> &array)
{
  code = -1;
  date = 10101;
  time = 0;

  if (numFields == 0) array.resize(gridsize0 * numLevels);

  cdo_print("Enter all %zu elements of timestep %d!", gridsize0 * numLevels, numFields + 1);

  auto rval = input_darray(stdin, gridsize0 * numLevels, array);

  if (numFields > 0 && rval == 0) return -1;

  if (rval != gridsize0 * numLevels) cdo_abort("Too few input elements (%zu of %zu)!", rval, gridsize0 * numLevels);

  if (feof(stdin)) return -1;

  return 0;
}

static int
read_record_ext(int numFields, size_t &gridsize0, int &code, int &date, int &time, Varray<double> &array, int &gridID,
                double &dlevel)
{
  cdo_print("Enter header (date,code,level,gridsize) of record %d (or EOF(=^D))!", numFields + 1);

  int ihead[4];
  auto rval = input_iarray(4, ihead);
  if (feof(stdin) && numFields == 0) cdo_abort("Too few header elements (%d of %d)!", rval, 4);
  if (feof(stdin)) return -1;
  if (rval != 4) cdo_abort("Invalid header input!");

  date = ihead[0];
  code = ihead[1];
  int level = ihead[2];
  size_t gridsize = ihead[3];

  time = 0;

  if (numFields == 0)
  {
    dlevel = level;
    gridsize0 = gridsize;
    array.resize(gridsize);
    gridID = gridCreate(GRID_GENERIC, gridsize);
  }
  else
  {
    if (gridsize != gridsize0) cdo_abort("Gridsize must not change!");
  }

  cdo_print("Enter all %zu elements of record %d!", gridsize, numFields + 1);

  rval = input_darray(stdin, gridsize, array);
  if (rval != gridsize) cdo_abort("Invalid data input!");

  return 0;
}

static int
read_record_srv(int numFields, size_t &gridsize0, int &code, int &date, int &time, Varray<double> &array, int &gridID,
                double &dlevel)
{
  cdo_print("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2) of record %d (or EOF(=^D))!", numFields + 1);

  int ihead[8];
  auto rval = input_iarray(8, ihead);
  if (feof(stdin) && numFields == 0) cdo_abort("Too few header elements (%d of %d)!", rval, 8);
  if (feof(stdin)) return -1;
  if (rval != 8) cdo_abort("Invalid header input!");

  code = ihead[0];
  int level = ihead[1];
  date = ihead[2];
  time = ihead[3];
  int nlon = ihead[4];
  int nlat = ihead[5];

  size_t gridsize = nlon * nlat;

  if (numFields == 0)
  {
    dlevel = level;
    gridsize0 = gridsize;
    array.resize(gridsize);
    gridID = gridCreate(GRID_GENERIC, gridsize);
    gridDefXsize(gridID, nlon);
    gridDefYsize(gridID, nlat);
  }
  else
  {
    if (gridsize != gridsize0) cdo_abort("Gridsize must not change!");
  }

  cdo_print("Enter all %zu elements of record %d!", gridsize, numFields + 1);

  rval = input_darray(stdin, gridsize, array);
  if (rval != gridsize) cdo_abort("Invalid data input!");

  return 0;
}

class Input : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Input",
    // clang-format off
    .operators = { { "input", InputHelp },
                   { "inputsrv", InputHelp },
                   { "inputext", InputHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 0, 1, NoRestriction },
  };
  inline static RegisterEntry<Input> registration = RegisterEntry<Input>();

  int INPUT{}, INPUTSRV{}, INPUTEXT{};

  int operatorID{};

  int outputFiletype = CDI_FILETYPE_GRB;
  CdoStreamID streamID = CDO_STREAM_UNDEF;
  int vlistID = -1;

  int gridID = -1;
  int zaxisID = -1;

public:
  void
  init() override
  {
    INPUT = module.get_id("input");
    INPUTSRV = module.get_id("inputsrv");
    INPUTEXT = module.get_id("inputext");

    operatorID = cdo_operator_id();

    if (operatorID == INPUT)
    {
      operator_input_arg("grid description file or name");
      operator_check_argc((cdo_operator_argc() == 1) ? 1 : 2);

      gridID = cdo_define_grid(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) zaxisID = cdo_define_zaxis(cdo_operator_argv(1));
    }
    else { operator_check_argc(0); }

    if (operatorID == INPUT) { outputFiletype = cdo_filetype(); }
    else if (operatorID == INPUTEXT)
    {
      outputFiletype = (CdoDefault::FileType == CDI_UNDEFID) ? CDI_FILETYPE_EXT : CdoDefault::FileType;
    }
    else if (operatorID == INPUTSRV)
    {
      outputFiletype = (CdoDefault::FileType == CDI_UNDEFID) ? CDI_FILETYPE_SRV : CdoDefault::FileType;
    }
  }

  void
  run() override
  {
    size_t gridsize0 = (gridID == -1) ? 0 : gridInqSize(gridID);
    int numLevels = (zaxisID == -1) ? 1 : zaxisInqSize(zaxisID);

    int taxisID = 0;
    int code = 0, date = 0, time = 0;
    double missval = 0;
    double dlevel = 0;
    Varray<double> array;

    int numFields = 0;
    int tsID = 0;
    while (true)
    {
      if (operatorID == INPUT)
      {
        if (-1 == read_record(numFields, numLevels, gridsize0, code, date, time, array)) break;
      }
      else if (operatorID == INPUTEXT)
      {
        if (-1 == read_record_ext(numFields, gridsize0, code, date, time, array, gridID, dlevel)) break;
      }
      else if (operatorID == INPUTSRV)
      {
        if (-1 == read_record_srv(numFields, gridsize0, code, date, time, array, gridID, dlevel)) break;
      }

      if (numFields == 0)
      {
        if (zaxisID == -1)
        {
          zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
          zaxisDefLevels(zaxisID, &dlevel);
        }

        vlistID = vlistCreate();
        int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
        vlistDefVarParam(vlistID, varID, cdiEncodeParam(code, 255, 255));

        missval = vlistInqVarMissval(vlistID, varID);

        taxisID = cdo_taxis_create(TAXIS_RELATIVE);
        vlistDefTaxis(vlistID, taxisID);

        streamID = cdo_open_write(0, outputFiletype);
        cdo_def_vlist(streamID, vlistID);
      }

      CdiDateTime vDateTime{};
      vDateTime.date = cdiDate_set(date);
      vDateTime.time = cdiTime_set(time);
      taxisDefVdatetime(taxisID, vDateTime);
      cdo_def_timestep(streamID, tsID);

      constexpr int varID = 0;
      for (int levelID = 0; levelID < numLevels; ++levelID)
      {
        auto offset = gridsize0 * levelID;
        auto numMissVals = array_num_mv(gridsize0, &array[offset], missval);
        cdo_def_field(streamID, varID, levelID);
        cdo_write_field(streamID, &array[offset], numMissVals);
      }

      numFields++;
      tsID++;
    }
  }

  void
  close() override
  {
    if (streamID)
    {
      cdo_stream_close(streamID);
      vlistDestroy(vlistID);
    }
  }
};
