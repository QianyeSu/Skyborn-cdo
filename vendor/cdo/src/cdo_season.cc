#include "cdo_season.h"
#include "cdo_output.h"
#include "cdo_options.h"
#include "util_string.h"

#include <stdlib.h>

const char *seasonNamesDec[4] = { "DJF", "MAM", "JJA", "SON" };
const char *seasonNamesJan[4] = { "JFM", "AMJ", "JAS", "OND" };

SeasonStart
get_season_start(void)
{
  static auto seasonStart = SeasonStart::DEC;
  static auto doEnvRead = true;

  if (doEnvRead)
    {
      doEnvRead = false;

      auto envString = getenv_string("CDO_SEASON_START");
      if (envString.size())
        {
          // clang-format off
          if      (envString == "DEC") seasonStart = SeasonStart::DEC;
          else if (envString == "JAN") seasonStart = SeasonStart::JAN;
          // clang-format on

          if (Options::cdoVerbose)
            {
              // clang-format off
              if      (seasonStart == SeasonStart::DEC) cdo_print("Set SEASON_START to December");
              else if (seasonStart == SeasonStart::JAN) cdo_print("Set SEASON_START to January");
              // clang-format on
            }
        }
    }

  return seasonStart;
}

const char **
get_season_name(void)
{
  return (get_season_start() == SeasonStart::DEC) ? seasonNamesDec : seasonNamesJan;
}

int
month_to_season(int month)
{
  if (month < 0 || month > 16) cdo_abort("Month %d out of range!", month);

  auto zmonth = (get_season_start() == SeasonStart::DEC) ? (month % 12) : (month - 1);
  auto season = (month <= 12) ? zmonth / 3 : month - 13;

  if (season < 0 || season > 3) cdo_abort("Season %d out of range!", season + 1);

  return season;
}
