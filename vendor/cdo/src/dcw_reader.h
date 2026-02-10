#ifndef DCW_READER_H
#define DCW_READER_H

#include <string>
#include <vector>

struct Region
{
  double west = 180.0, east = -180.0, south = 90.0, north = -90.0;
};

struct DCW_Country  // Information per country
{
  char continent[4] = { 0 };  // 2-char continent code (EU, NA, SA, AF, AU, AN)
  char code[4] = { 0 };       // 2-char country code ISO 3166-1 (e.g., NO, US)
  char name[80] = { 0 };      // Full name of the country
};

struct DCW_State  // Information per state
{
  char country[4] = { 0 };  // 2-char country code ISO 3166-1 (e.g., BR, US)
  char code[4] = { 0 };     // 2/3-char state codes for US, Canada, China, Argentina, Australia, Brazil, Russia (e.g., TX)
  char name[80] = { 0 };    // Full name of the state
};

struct DCW_Lists
{
  std::vector<DCW_Country> countries;
  std::vector<DCW_State> states;
};

int dcw_load_lists(DCW_Lists &dcw_lists);
void dcw_print_path();
void dcw_print_countries(const DCW_Lists &dcw_lists);
void dcw_print_states(const DCW_Lists &dcw_lists);
void dcw_print_polygons(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList);
void dcw_sort_countries(DCW_Lists &dcw_lists);
int dcw_get_region(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList, Region &region);
std::vector<std::string> dcw_expand_code_list(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList);
int dcw_get_lonlat(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList, std::vector<double> &lon,
                   std::vector<double> &lat);

#endif
