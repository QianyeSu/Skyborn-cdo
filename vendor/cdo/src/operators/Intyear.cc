/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Intyear    intyear         Year interpolation
*/

#include <cdi.h>

#include "cdo_rlimit.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"
#include "util_string.h"

static size_t
intlin_year(double fac1, double fac2, size_t gridsize, Varray<double> const &array1, Varray<double> const &array2,
            Varray<double> &array3, bool withMissval, double missval1, double missval2)
{
  size_t numMissVals3 = 0;

  if (withMissval)
  {
    for (size_t i = 0; i < gridsize; ++i)
    {
      if (fp_is_not_equal(array1[i], missval1) && fp_is_not_equal(array2[i], missval2))
      {
        array3[i] = array1[i] * fac1 + array2[i] * fac2;
      }
      else
      {
        array3[i] = missval1;
        numMissVals3++;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < gridsize; ++i) array3[i] = array1[i] * fac1 + array2[i] * fac2;
  }

  return numMissVals3;
}

class Intyear : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Intyear",
    .operators = { { "intyear", IntyearHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, OBASE, NoRestriction },
  };
  inline static RegisterEntry<Intyear> registration = RegisterEntry<Intyear>();
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3;

  std::vector<int> iyears;
  int numYears{};

  VarList varList1{};
  VarList varList2{};

  std::vector<CdoStreamID> streamIDs;

public:
  void
  init() override
  {
    operator_input_arg("years");

    iyears = cdo_argv_to_intarr(cdo_get_oper_argv());
    numYears = iyears.size();

    cdo::set_numfiles(numYears + 8);

    streamIDs.resize(numYears);

    streamID1 = cdo_open_read(0);
    streamID2 = cdo_open_read(1);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto vlistID2 = cdo_stream_inq_vlist(streamID2);
    auto vlistID3 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);
    varList2 = VarList(vlistID2);
    varList_compare(varList1, varList2);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = vlistInqTaxis(vlistID2);
    taxisID3 = taxisDuplicate(taxisID1);
    if (taxisHasBounds(taxisID3)) taxisDeleteBounds(taxisID3);
    vlistDefTaxis(vlistID3, taxisID3);

    auto fileSuffix = FileUtils::gen_suffix(cdo_inq_filetype(streamID1), vlistID1, cdo_get_stream_name(0));

    for (int iy = 0; iy < numYears; iy++)
    {
      auto fileName = cdo_get_obase() + string_format("%04d", iyears[iy]);
      if (fileSuffix.size() > 0) fileName += fileSuffix;

      streamIDs[iy] = open_write(fileName);
      cdo_def_vlist(streamIDs[iy], vlistID3);
    }
  }

  void
  run() override
  {
    Varray<double> array1(varList1.gridsizeMax());
    Varray<double> array2(varList1.gridsizeMax());
    Varray<double> array3(varList1.gridsizeMax());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;
      numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) cdo_abort("Too few timesteps in second inputfile!");

      auto vDateTime1 = taxisInqVdatetime(taxisID1);
      auto vDateTime2 = taxisInqVdatetime(taxisID2);
      auto year1 = vDateTime1.date.year;
      auto year2 = vDateTime2.date.year;

      for (int iy = 0; iy < numYears; iy++)
      {
        if (iyears[iy] < year1 || iyears[iy] > year2)
          cdo_abort("Year %d out of bounds (first year %d; last year %d)!", iyears[iy], year1, year2);

        auto vDateTime3 = vDateTime1;
        vDateTime3.date.year = iyears[iy];
        taxisDefVdatetime(taxisID3, vDateTime3);
        cdo_def_timestep(streamIDs[iy], tsID);
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        (void) cdo_inq_field(streamID1);
        auto [varID, levelID] = cdo_inq_field(streamID2);

        size_t numMissVals1, numMissVals2;
        cdo_read_field(streamID1, array1.data(), &numMissVals1);
        cdo_read_field(streamID2, array2.data(), &numMissVals2);

        auto gridsize = varList1.vars[varID].gridsize;
        auto missval1 = varList1.vars[varID].missval;
        auto missval2 = varList2.vars[varID].missval;
        auto withMissval = (numMissVals1 || numMissVals2);

        for (int iy = 0; iy < numYears; iy++)
        {
          auto fac1 = ((double) year2 - iyears[iy]) / (year2 - year1);
          auto fac2 = ((double) iyears[iy] - year1) / (year2 - year1);

          auto numMissVals3 = intlin_year(fac1, fac2, gridsize, array1, array2, array3, withMissval, missval1, missval2);

          cdo_def_field(streamIDs[iy], varID, levelID);
          cdo_write_field(streamIDs[iy], array3.data(), numMissVals3);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    for (auto const &streamID : streamIDs) cdo_stream_close(streamID);

    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
