/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Mastrfu    mastrfu         Mass stream function
*/

#include <cdi.h>

#include "varray.h"
#include "process_int.h"
#include <mpim_grid.h>
#include "cdo_zaxis.h"

static void
mastrfu(int gridID, int zaxisID, Varray2D<double> const &field1, Varray2D<double> &field2, size_t numMissVals, double missval)
{
  auto fact = 4.0 * std::atan(1.0) * 6371000.0 / 9.81;

  auto nlat = gridInqSize(gridID);
  auto nlev = zaxisInqSize(zaxisID);
  Varray<double> phi(nlat), cosphi(nlat), plevel(nlev);

  cdo_zaxis_inq_levels(zaxisID, plevel.data());

  auto plevel1 = plevel[0];
  auto plevelN = plevel[nlev - 1];
  if (plevel1 < plevelN) cdo_abort("The 3d pressure level data is upside down! Use the operator invertlev to invert the levels.");

  gridInqYvals(gridID, phi.data());

  auto units = cdo::inq_key_string(gridID, CDI_YAXIS, CDI_KEY_UNITS);

  if (units.rfind("degree", 0) == 0)
    for (size_t ilat = 0; ilat < nlat; ilat++) phi[ilat] *= DEG2RAD;

  for (size_t ilat = 0; ilat < nlat; ilat++) phi[ilat] = std::sin(phi[ilat]);

  for (size_t ilat = 0; ilat < nlat; ilat++) cosphi[ilat] = std::sqrt(1.0 - phi[ilat] * phi[ilat]);

  for (int ilev = 0; ilev < nlev; ilev++)
    for (size_t ilat = 0; ilat < nlat; ilat++) field2[ilev][ilat] = 0.0;

  if (numMissVals == 0)
  {
    for (int ilev = nlev - 1; ilev >= 0; ilev--)
      for (int n = ilev; n < nlev - 1; ++n)
        for (size_t ilat = 0; ilat < nlat; ilat++)
        {
          field2[ilev][ilat] += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
        }
  }
  else
  {
    for (size_t ilat = 0; ilat < nlat; ilat++)
      for (int ilev = nlev - 1; ilev >= 0; ilev--)
        for (int n = ilev; n < nlev - 1; ++n)
        {
          if (fp_is_equal(field1[n][ilat], missval) || fp_is_equal(field1[n + 1][ilat], missval))
          {
            field2[ilev][ilat] = missval;
            break;
          }
          else
            field2[ilev][ilat] += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
        }
  }
}

class Mastrfu : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Mastrfu",
    .operators = { { "mastrfu", MastrfuHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Mastrfu> registration = RegisterEntry<Mastrfu>();

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int zaxisID{};

  int gridID{};
  size_t nlat{};

  double missval{};

public:
  void
  init() override
  {
    streamID1 = cdo_open_read(0);

    operator_check_argc(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    VarList varList1(vlistID1);

    auto numVars = varList1.numVars();
    if (numVars != 1) cdo_abort("This operator works only with one variable!");

    auto const &var0 = varList1.vars[0];
    if (var0.code > 0 && var0.code != 132) cdo_warning("Unexpected code %d!", var0.code);

    missval = var0.missval;
    gridID = var0.gridID;
    zaxisID = var0.zaxisID;

    if (var0.zaxisType != ZAXIS_PRESSURE && var0.zaxisType != ZAXIS_GENERIC)
    {
      cdo_warning("Unexpected vertical grid %s!", cdo::inq_key_string(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME));
    }

    if (gridInqXsize(var0.gridID) > 1) cdo_abort("Grid must be a zonal mean!");

    nlat = gridInqSize(var0.gridID);

    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    vlistDefVarCode(vlistID2, 0, 272);
    cdiDefKeyString(vlistID2, 0, CDI_KEY_NAME, "mastrfu");
    cdiDefKeyString(vlistID2, 0, CDI_KEY_LONGNAME, "mass stream function");
    cdiDefKeyString(vlistID2, 0, CDI_KEY_UNITS, "kg/s");
    vlistDefVarDatatype(vlistID2, 0, CDI_DATATYPE_FLT32);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto numLevels = zaxisInqSize(zaxisID);
    Varray2D<double> array1(numLevels, Varray<double>(nlat));
    Varray2D<double> array2(numLevels, Varray<double>(nlat));

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      size_t numMissVals = 0;
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        size_t numMissVals1;
        cdo_read_field(streamID1, array1[levelID].data(), &numMissVals1);
        numMissVals += numMissVals1;
      }

      mastrfu(gridID, zaxisID, array1, array2, numMissVals, missval);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        int varID = 0;
        int levelID = fieldID;
        cdo_def_field(streamID2, varID, levelID);
        numMissVals = array_num_mv(nlat, array2[levelID].data(), missval);
        cdo_write_field(streamID2, array2[levelID].data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
  }
};
