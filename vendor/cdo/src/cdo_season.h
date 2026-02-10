#ifndef CDO_SEASON_H
#define CDO_SEASON_H

enum struct SeasonStart
{
  DEC,
  JAN
};

SeasonStart get_season_start(void);
const char **get_season_name(void);
int month_to_season(int month);

#endif
