/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

*/

#include "cdi.h"
#include "julian_date.h"

#include "process_int.h"
#include "interpol.h"
#include "varray.h"
#include "printinfo.h"

static int
read_nextpos(std::FILE *fp, int calendar, JulianDate &julianDate, double &xpos, double &ypos)
{
  xpos = 0.0;
  ypos = 0.0;

  int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0, ms = 0;
  int stat = std::fscanf(fp, "%d-%d-%d %d:%d:%d %lf %lf", &year, &month, &day, &hour, &minute, &second, &xpos, &ypos);

  CdiDateTime dateTime{};
  dateTime.date = cdiDate_set(10101);
  if (stat != EOF)
  {
    dateTime.date = cdiDate_encode(year, month, day);
    dateTime.time = cdiTime_encode(hour, minute, second, ms);
  }

  julianDate = julianDate_encode(calendar, dateTime);

  return stat;
}

class Intgridtraj : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intgridtraj",
    .operators = { { "intgridtraj" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Intgridtraj> registration = RegisterEntry<Intgridtraj>();

  size_t numMissVals = 0;
  double xpos{}, ypos{};
  int calendar = CALENDAR_STANDARD;

  CdoStreamID streamID1 = CDO_STREAM_UNDEF;
  CdoStreamID streamID2 = CDO_STREAM_UNDEF;

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };
  int gridID2{};
  int numVars{};

  VarList varList1{};

  Varray2D<double> vardata1;
  Varray2D<double> vardata2;

  JulianDate julianDate{};

  std::FILE *fp{};

public:
  void
  init() override
  {
    operator_input_arg("filename with grid trajectories");
    operator_check_argc(1);

    auto posfile = cdo_operator_argv(0).c_str();
    fp = std::fopen(posfile, "r");
    if (fp == nullptr) cdo_abort("Open failed on %s!", posfile);

    read_nextpos(fp, calendar, julianDate, xpos, ypos);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);
    numVars = varList1.numVars();

    vardata1.resize(numVars);
    vardata2.resize(numVars);
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto gridsize = varList1.vars[varID].gridsize;
      auto nlevels = varList1.vars[varID].nlevels;
      vardata1[varID].resize(gridsize * nlevels);
      vardata2[varID].resize(gridsize * nlevels);
    }

    gridID2 = gridCreate(GRID_TRAJECTORY, 1);
    gridDefXsize(gridID2, 1);
    gridDefYsize(gridID2, 1);
    gridDefXvals(gridID2, &xpos);
    gridDefYvals(gridID2, &ypos);

    vlistID2 = vlistDuplicate(vlistID1);

    auto numGrids = vlistNumGrids(vlistID1);
    for (int index = 0; index < numGrids; ++index)
    {
      auto gridID1 = vlistGrid(vlistID1, index);

      if (gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
        cdo_abort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID1)));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);
  }

  void
  run() override
  {
    Field field1, field2;
    field1.resize(varList1.gridsizeMax());
    field2.resize(1);

    auto maxFields = varList1.maxFields();
    std::vector<FieldInfo> fieldInfoList(maxFields);

    int tsID = 0;
    auto numFields = cdo_stream_inq_timestep(streamID1, tsID++);
    auto julianDate1 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));
    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(streamID1);
      auto gridsize = varList1.vars[varID].gridsize;
      auto offset = gridsize * levelID;
      auto single1 = &vardata1[varID][offset];
      cdo_read_field(streamID1, single1, &numMissVals);
      if (numMissVals) cdo_abort("Missing values unsupported for this operator!");
    }

    int tsIDo = 0;
    while (julianDate_to_seconds(julianDate1) <= julianDate_to_seconds(julianDate))
    {
      numFields = cdo_stream_inq_timestep(streamID1, tsID++);
      if (numFields == 0) break;
      auto julianDate2 = julianDate_encode(calendar, taxisInqVdatetime(taxisID1));

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        fieldInfoList[fieldID].set(varID, levelID);

        auto gridsize = varList1.vars[varID].gridsize;
        auto offset = gridsize * levelID;
        auto single2 = &vardata2[varID][offset];
        cdo_read_field(streamID1, single2, &numMissVals);
        if (numMissVals) cdo_abort("Missing values unsupported for this operator!");
      }

      while (julianDate_to_seconds(julianDate) < julianDate_to_seconds(julianDate2))
      {
        if (julianDate_to_seconds(julianDate) >= julianDate_to_seconds(julianDate1)
            && julianDate_to_seconds(julianDate) < julianDate_to_seconds(julianDate2))
        {
          if (streamID2 == CDO_STREAM_UNDEF)
          {
            streamID2 = cdo_open_write(1);
            cdo_def_vlist(streamID2, vlistID2);
          }

          taxisDefVdatetime(taxisID2, julianDate_decode(calendar, julianDate));
          cdo_def_timestep(streamID2, tsIDo++);

          auto deltat = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1));
          auto fac1 = julianDate_to_seconds(julianDate_sub(julianDate2, julianDate)) / deltat;
          auto fac2 = julianDate_to_seconds(julianDate_sub(julianDate, julianDate1)) / deltat;
          // printf("      %f %f %f %f %f\n", julianDate_to_seconds(julianDate), julianDate_to_seconds(julianDate1),
          // julianDate_to_seconds(julianDate2), fac1, fac2);
          for (int fieldID = 0; fieldID < numFields; ++fieldID)
          {
            auto [varID, levelID] = fieldInfoList[fieldID].get();

            auto &var1 = varList1.vars[varID];
            auto gridsize = var1.gridsize;
            auto offset = gridsize * levelID;
            auto single1 = &vardata1[varID][offset];
            auto single2 = &vardata2[varID][offset];

            for (size_t i = 0; i < gridsize; ++i) field1.vec_d[i] = single1[i] * fac1 + single2[i] * fac2;

            field1.grid = var1.gridID;
            field1.missval = var1.missval;
            field1.numMissVals = numMissVals;
            field2.grid = gridID2;
            field2.numMissVals = 0;

            intgrid_bil(field1, field2);

            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field2);
          }
        }

        if (read_nextpos(fp, calendar, julianDate, xpos, ypos) == EOF) break;
        gridDefXvals(gridID2, &xpos);
        gridDefYvals(gridID2, &ypos);
      }

      julianDate1 = julianDate2;
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto vardatap = vardata1[varID];
        vardata1[varID] = vardata2[varID];
        vardata2[varID] = vardatap;
      }
    }

    if (tsIDo == 0)
    {
      auto dt = julianDate_decode(calendar, julianDate);
      cdo_warning("Date/time %s %s not found!", date_to_string(dt.date), time_to_string(dt.time));
    }
  }

  void
  close() override
  {
    std::fclose(fp);
    if (streamID2 != CDO_STREAM_UNDEF) cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
