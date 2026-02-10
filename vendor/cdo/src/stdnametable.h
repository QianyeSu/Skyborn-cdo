/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef STDNAMETABLE_H
#define STDNAMETABLE_H

#include <string>

enum stdnameid
{
  air_pressure,
  pressure_thickness,
  surface_geopotential,
  geopotential,
  air_temperature,
  specific_humidity,
  surface_air_pressure,
  air_density,
  air_pressure_at_sea_level,
  geopotential_height,
  geometric_height_at_full_level_center,
  geometric_height_at_half_level_center
};

int var_echamcode(int varid);
const char *var_name(int varid);
const char *var_stdname(int varid);
const char *var_units(int varid);

int stdname_to_echamcode(std::string const &stdname);

struct GribCodes
{
  int geopot = 0;
  int ta = 0;
  int hus = 0;
  int ps = 0;
  int lsp = 0;
  int gheight = 0;
  int wind = 0;
  int uwind = 0;
  int vwind = 0;
};

GribCodes echam_gribcodes();
GribCodes wmo_gribcodes();
GribCodes hirlam_harmonie_gribcodes();

#endif
