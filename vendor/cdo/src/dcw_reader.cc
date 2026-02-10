// g++ -Wall -I/opt/local/include read_dcw.cc -L/opt/local/lib -lnetcdf

// modified code from gmt_dcw.c

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#include "dcw_reader.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <algorithm>  // sort

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

#define DCW_SITE "ftp://ftp.soest.hawaii.edu/gmt"

struct DCW_Country_State  // Information per country with state
{
  char country[4];  // 2/3-char country code ISO 3166-1 (e.g. BR, US) for countries with states
};

static bool
dcw_get_path(const char *name, const char *suffix, char *path)
{
  static bool readEnv = true;
  static char *DCW_Dir = nullptr;
  bool found = false;

  if (readEnv)
    {
      DCW_Dir = getenv("DIR_DCW");
      readEnv = false;
    }

  if (DCW_Dir)
    {
      std::snprintf(path, PATH_MAX, "%s/%s%s", DCW_Dir, name, suffix);
      if (access(path, R_OK) == 0) found = true;
    }
  else
    {
      fprintf(stderr, "Environment variable DIR_DCW not set!\n");
      return found;
    }

  if (!found)
    {
      printf("Unable to find or open the Digital Chart of the World\n");
      printf("Perhaps you did not install this file in DIR_DCW?\n");
      printf("Use your package manager to install package dcw-gmt.\n");
      printf("Alternatively, get the latest dcw-gmt-<version>.tar.gz or dcw-gmt-<version>.zip from the %s.\n", DCW_SITE);
    }

  return found;
}

int
dcw_load_lists(DCW_Lists &dcw_lists)
{
  unsigned int dim[3];
  // Open and read list of countries and states and return via two struct and one char arrays plus dimensions in dim
  size_t n_alloc = 300;
  unsigned int k, n;
  char line[BUFSIZ] = { 0 };
  // DCW_Country_State *Country_State = NULL;

  char path[PATH_MAX] = { 0 };
  if (!dcw_get_path("dcw-countries", ".txt", path)) return -1;

  // Get countries first
  auto fp = std::fopen(path, "r");
  if (fp == nullptr)
    {
      fprintf(stderr, "Unable to open file %s [permission trouble?]\n", path);
      return -1;
    }

  auto &countries = dcw_lists.countries;
  countries.resize(n_alloc);
  k = 0;
  while (fgets(line, BUFSIZ, fp))
    {
      if (line[0] == '#') continue;  // Skip comments
      std::sscanf(line, "%s %s %[^\n]", countries[k].continent, countries[k].code, countries[k].name);
      k++;
      if (k == n_alloc)
        {
          n_alloc += 100;
          countries.resize(n_alloc);
        }
    }

  std::fclose(fp);

  dim[0] = k;  // Number of countries

  countries.resize(k);

  // Get states
  if (!dcw_get_path("dcw-states", ".txt", path)) { return -1; }

  fp = std::fopen(path, "r");
  if (fp == nullptr)
    {
      fprintf(stderr, "Unable to open file %s [permission trouble?]\n", path);
      return -1;
    }

  auto &states = dcw_lists.states;
  states.resize(n_alloc);
  k = 0;
  n = 1;
  while (fgets(line, BUFSIZ, fp))
    {
      if (line[0] == '#') continue;  // Skip comments
      std::sscanf(line, "%s %s %[^\n]", states[k].country, states[k].code, states[k].name);
      if (k && std::strcmp(states[k].country, states[k - 1].country)) n++;  // New country with states
      k++;
      if (k == n_alloc)
        {
          n_alloc += 100;
          states.resize(n_alloc);
        }
    }

  std::fclose(fp);

  dim[1] = k;  // Number of states
  states.resize(k);

  // Get list of countries with states

  dim[2] = n;  // Number of countries with states
  /*
  if (0)  // Wants list returned
  {
    Country_State = (DCW_Country_State *) malloc(n * sizeof(struct DCW_Country_State));
    std::memcpy(Country_State[0].country, states[0].country, 4);
    for (k = n = 1; k < dim[1]; ++k)
      {
        if (std::strcmp(states[k].country, states[k - 1].country)) std::memcpy(Country_State[n++].country, states[k].country, 4);
      }
    // *CS = Country_State;
  }
  */

  // fprintf(stderr, "# DCW: Found %u countries, %u countries with states, and %u states\n", dim[0], dim[2], dim[1]);

  return 0;
}

static int
dcw_find_country(std::string const &code, const std::vector<DCW_Country> &list)
{
  int low = 0, high = (int) list.size() - 1;

  while (low <= high)
    {
      auto midpoint = low + (high - low) / 2;

      // check to see if value is equal to item in array
      auto way = std::strcmp(code.c_str(), list[midpoint].code);
      if (way == 0) return midpoint;

      if (way < 0)
        high = midpoint - 1;
      else
        low = midpoint + 1;
    }

  // item was not found
  return -1;
}

static int
dcw_find_state(std::string const &code, const std::vector<DCW_State> &list)
{
  auto country = code.substr(0, 2);
  auto state = code.substr(2);
  for (size_t i = 0, n = list.size(); i < n; ++i)
    {
      if (country == list[i].country && state == list[i].code) return i;
    }

  // item was not found
  return -1;
}

#ifdef HAVE_LIBNETCDF
static int
nc_get_minmax(int ncid, std::string const &name, double &minval, double &maxval)
{
  int varid;
  if (nc_inq_varid(ncid, name.c_str(), &varid)) return 1;
  if (nc_get_att_double(ncid, varid, "min", &minval)) return 1;
  if (nc_get_att_double(ncid, varid, "max", &maxval)) return 1;

  return 0;
}
#endif

static int
dcw_open_nc()
{
  char path[PATH_MAX] = { 0 };
  if (!dcw_get_path("dcw-gmt", ".nc", path)) return 1;

  int ncid = -1;
#ifdef HAVE_LIBNETCDF
  auto status = nc_open(path, NC_NOWRITE, &ncid);
  if (status)
    {
      fprintf(stderr, "Cannot open file %s!\n", path);
      return -1;
    }
#else
  fprintf(stderr, "dcw_open_nc failed: NetCDF support not compiled in!\n");
#endif

  return ncid;
}

int
dcw_get_region(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList, Region &region)
{
  auto ncid = dcw_open_nc();
  if (ncid < 0) return 1;

  int numEntries = 0;
  for (auto &code : codeList)
    {
      auto ks = (code.size() > 2) ? dcw_find_state(code, dcw_lists.states) : dcw_find_country(code, dcw_lists.countries);
      if (ks == -1)
        {
          fprintf(stderr, "No country code matching <%s>!\n", code.c_str());
          return 3;
        }

      // auto &country = dcw_lists.countries[ks];
      // auto is_Antarctica = (!std::strncmp(country.code, "AQ", 2U));

      auto xname = code + "_lon";
      auto yname = code + "_lat";

      double west = 0.0, east = 0.0, south = 0.0, north = 0.0;
#ifdef HAVE_LIBNETCDF
      if (nc_get_minmax(ncid, xname, west, east)) continue;
      if (nc_get_minmax(ncid, yname, south, north)) continue;
#endif
      // if (west <= 0 && east >=0) fprintf(stderr, "%s: %g %g\n", code.c_str(), west, east);
      // if (west >= 180 && east >=180) fprintf(stderr, "%s: %g %g\n", code.c_str(), west, east);
      if (west >= 180.0 && east > 180.0)
        {
          west -= 360.0;
          east -= 360.0;
        }
      if (numEntries > 0)
        {
          if (east <= 0.0 && region.west >= 0.0)
            {
              west += 360.0;
              east += 360.0;
            }
          else if (east >= 180.0 && region.east <= 0.0)
            {
              region.west += 360.0;
              region.east += 360.0;
            }
        }
      region.west = std::min(region.west, west);
      region.south = std::min(region.south, south);
      region.east = std::max(region.east, east);
      region.north = std::max(region.north, north);

      numEntries++;
    }

#ifdef HAVE_LIBNETCDF
  nc_close(ncid);
#endif

  if (numEntries == 0)
    {
      printf("Empty code list!\n");
      return 2;
    }

  return 0;
}
/*
#ifdef HAVE_LIBNETCDF
static int
nc_get_data(int ncid, const char *name, std::vector<double> &data)
{
  int varid;
  double scale, minval;
  if (nc_inq_varid(ncid, name, &varid)) return 1;
  if (nc_get_att_double(ncid, varid, "min", &minval)) return 1;
  if (nc_get_att_double(ncid, varid, "scale", &scale)) return 1;

  nc_type xtype;
  int nvdims, nvatts, dimids[9];
  if (nc_inq_var(ncid, varid, NULL, &xtype, &nvdims, dimids, &nvatts)) return 1;

  if (nvdims != 1) return 1;

  size_t np = 0;
  if (nc_inq_dimlen(ncid, dimids[0], &np)) return 1;

  std::vector<unsigned short> shortData(np);
  if (nc_get_var_ushort(ncid, varid, shortData.data())) return 1;

  data.resize(np);
  scale = 1.0 / scale;
  for (size_t i = 0; i < np; ++i)  // Unpack
    data[i] = (shortData[i] == 65535U) ? 0.0 : shortData[i] * scale + minval;

  return 0;
}
#endif
*/
#ifdef HAVE_LIBNETCDF
static int
nc_get_lonlat(int ncid, std::string const &xname, std::string const &yname, std::vector<double> &x, std::vector<double> &y)
{
  int xvarid;
  double xscale, xmin;
  if (nc_inq_varid(ncid, xname.c_str(), &xvarid)) return 1;
  if (nc_get_att_double(ncid, xvarid, "min", &xmin)) return 1;
  if (nc_get_att_double(ncid, xvarid, "scale", &xscale)) return 1;

  int yvarid;
  double yscale, ymin;
  if (nc_inq_varid(ncid, yname.c_str(), &yvarid)) return 1;
  if (nc_get_att_double(ncid, yvarid, "min", &ymin)) return 1;
  if (nc_get_att_double(ncid, yvarid, "scale", &yscale)) return 1;

  nc_type xtype;
  int nvdims, nvatts, dimids[9];
  if (nc_inq_var(ncid, xvarid, NULL, &xtype, &nvdims, dimids, &nvatts)) return 1;

  if (nvdims != 1) return 1;

  size_t np = 0;
  if (nc_inq_dimlen(ncid, dimids[0], &np)) return 1;

  std::vector<unsigned short> xshort(np), yshort(np);
  if (nc_get_var_ushort(ncid, xvarid, xshort.data())) return 1;
  if (nc_get_var_ushort(ncid, yvarid, yshort.data())) return 1;

  x.resize(np);
  y.resize(np);
  xscale = 1.0 / xscale;
  yscale = 1.0 / yscale;
  for (size_t i = 0; i < np; ++i)  // Unpack
    {
      x[i] = (xshort[i] == 65535U) ? 0.0 : xshort[i] * xscale + xmin;
      y[i] = (xshort[i] == 65535U) ? 0.0 : yshort[i] * yscale + ymin;
      //  use ^ xshort to check for undefined values !!!
    }

  return 0;
}
#endif

int
dcw_get_lonlat(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList, std::vector<double> &lon,
               std::vector<double> &lat)
{
  auto ncid = dcw_open_nc();
  if (ncid < 0) return 1;

  size_t numVals = 0;
  int numEntries = 0;
  for (auto &code : codeList)
    {
      auto ks = (code.size() > 2) ? dcw_find_state(code, dcw_lists.states) : dcw_find_country(code, dcw_lists.countries);
      if (ks == -1)
        {
          fprintf(stderr, "No country code matching <%s>!\n", code.c_str());
          return 3;
        }

      // auto &country = dcw_lists.countries[ks];
      // auto is_Antarctica = (!std::strncmp(country.code, "AQ", 2U));

      auto xname = code + "_lon";
      auto yname = code + "_lat";

      std::vector<double> x, y;
#ifdef HAVE_LIBNETCDF
      if (nc_get_lonlat(ncid, xname, yname, x, y)) continue;
#endif
      auto offset = lon.size();
      numVals += x.size();
      lon.resize(numVals);
      lat.resize(numVals);
      for (size_t i = 0, n = x.size(); i < n; ++i) lon[offset + i] = x[i];
      for (size_t i = 0, n = y.size(); i < n; ++i) lat[offset + i] = y[i];

      numEntries++;
    }

#ifdef HAVE_LIBNETCDF
  nc_close(ncid);
#endif

  if (numEntries == 0)
    {
      printf("Empty code list!\n");
      return 2;
    }

  return 0;
}

static int
dcw_print_lonlat(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList)
{
  auto ncid = dcw_open_nc();
  if (ncid < 0) return 1;

  int numEntries = 0;
  for (auto &code : codeList)
    {
      auto ks = (code.size() > 2) ? dcw_find_state(code, dcw_lists.states) : dcw_find_country(code, dcw_lists.countries);
      if (ks == -1)
        {
          printf("No country code matching <%s>!\n", code.c_str());
          return 3;
        }

      auto &country = dcw_lists.countries[ks];
      // auto is_Antarctica = (!std::strncmp(country.code, "AQ", 2U));

      auto xname = code + "_lon";
      auto yname = code + "_lat";

      std::vector<double> x, y;
#ifdef HAVE_LIBNETCDF
      if (nc_get_lonlat(ncid, xname, yname, x, y)) continue;
#endif
      size_t nseg = 0;
      auto n = x.size();
      for (size_t i = 0; i < n; ++i)
        {
          if (x[i] == 0.0 && y[i] == 0.0)
            printf("> %s  %s  Segment %zu\n", country.code, country.name, nseg++);
          else
            printf("%.12g  %.12g\n", x[i], y[i]);
        }

      numEntries++;
    }

#ifdef HAVE_LIBNETCDF
  nc_close(ncid);
#endif

  if (numEntries == 0)
    {
      printf("Empty code list!\n");
      return 2;
    }

  return 0;
}

static bool
compare_code(const DCW_Country &a, const DCW_Country &b)
{
  return (std::strcmp(a.code, b.code) < 0);
}

std::vector<std::string>
dcw_expand_code_list(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList)
{
  std::vector<std::string> codeListExpand;
  codeListExpand.reserve(codeList.size());

  for (auto &code : codeList)
    {
      if (code.size() == 2 || code.size() == 4 || code.size() == 5) { codeListExpand.push_back(code); }
      else if (code.size() == 3 && code[0] == '=')
        {
          int n = 0;
          for (auto &country : dcw_lists.countries)
            {
              if (std::strncmp(country.continent, &code[1], 2)) continue;  // Not this one
              codeListExpand.push_back(country.code);
              n++;
            }

          if (n == 0) fprintf(stderr, "Wrong DCW country or continent code <%s>\n", code.c_str() + 1);
        }
      else { fprintf(stderr, "Wrong DCW country or continent code <%s>\n", code.c_str()); }
    }

  return codeListExpand;
}

void
dcw_print_path()
{
  char path[PATH_MAX] = { 0 };
  dcw_get_path("", "", path);
  printf("%s\n", path);
}

void
dcw_print_countries(const DCW_Lists &dcw_lists)
{
  printf("# DCW countries\n");
  printf("# List of continent-code country-code country-name\n");

  for (auto &country : dcw_lists.countries) printf("%s %s %s\n", country.continent, country.code, country.name);
}

void
dcw_print_states(const DCW_Lists &dcw_lists)
{
  printf("# DCW states\n");
  printf("# List of country-code state-code state-name\n");

  for (auto &state : dcw_lists.states) printf("%s %s %s\n", state.country, state.code, state.name);
}

void
dcw_print_polygons(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList)
{
  if (dcw_print_lonlat(dcw_lists, codeList)) fprintf(stderr, "dcw_print_lonlat failed!\n");
}

void
dcw_sort_countries(DCW_Lists &dcw_lists)
{
  std::ranges::sort(dcw_lists.countries, compare_code);
}

/*
static void
analyze_segments(std::vector<double> const &lon, std::vector<double> const &lat)
{
  size_t nseg = 0;
  auto n = lon.size();
  size_t sbegin = n;
  size_t send = 0;
  for (size_t i = 0; i < n; ++i)
    {
      if (lon[i] == 0.0 && lat[i] == 0.0)
        {
          sbegin = i + 1;
        }
      else if ((i == (n - 1)) || (i < (n - 1) && lon[i + 1] == 0.0 && lat[i + 1] == 0.0))
        {
          double xmin = lon[sbegin], xmax = lon[sbegin], ymin = lat[sbegin], ymax = lat[sbegin];
          send = i;
          for (size_t k = sbegin + 1; k <= send; ++k)
            {
              xmin = std::min(xmin, lon[k]);
              xmax = std::max(xmax, lon[k]);
              ymin = std::min(ymin, lat[k]);
              ymax = std::max(ymax, lat[k]);
            }

          if (xmin < 180.0 && xmax > 180.0)
            {
              fprintf(stderr, "> Segment %zu\n", nseg++);
              fprintf(stderr, "> xmin, xmax, ymin, ymax %g %g %g %g\n", xmin, xmax, ymin, ymax);
            }
        }
    }
}

static void
dump_lonlat(std::vector<double> const &lon, std::vector<double> const &lat)
{
  size_t nseg = 0;
  auto n = lon.size();
  for (size_t i = 0; i < n; ++i)
    {
      if (lon[i] == 0.0 && lat[i] == 0.0)
        printf("> Segment %zu\n", nseg++);
      else
        printf("%.12g  %.12g\n", lon[i], lat[i]);
    }
}

static void
print_lonlat(const DCW_Lists &dcw_lists, std::vector<std::string> const &codeList)
{
  std::vector<double> lon, lat;
  if (dcw_get_lonlat(dcw_lists, codeList, lon, lat)) fprintf(stderr, "dcw_get_lonlat failed!\n");

  analyze_segments(lon, lat);
  dump_lonlat(lon, lat);
}

int main(void)
{
  DCW_Lists dcw_lists;
  if (dcw_load_lists(dcw_lists))
    {
      fprintf(stderr, "dcw_load_lists failed!\n");
      return -1;
    }

  // dcw_print_path();
  // dcw_print_countries(dcw_lists.countries);
  // dcw_print_states(dcw_lists.states);

  dcw_sort_countries(dcw_lists);

  std::vector<std::string> codeList;
  // codeList.push_back("DE");
  // codeList.push_back("FR");
  // codeList.push_back("PL");
  // codeList.push_back("PT");
  codeList.push_back("=EU");
  // codeList.push_back("=AF");
  // codeList.push_back("=AS");
  // codeList.push_back("=NA");
  // codeList.push_back("=OC");
  // codeList.push_back("=SA");

  printf("# Digital Chart of the World\n");
  printf("# Region for country:");
  for (auto &code : codeList) printf(" %s", code.c_str());
  printf("\n");

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  Region region;
  if (dcw_get_region(dcw_lists, codeList, region)) fprintf(stderr, "dcw_get_region() failed!\n");

  printf("#   West=%g  East=%g  South=%g  North=%g\n", region.west, region.east, region.south, region.north);
  printf("#\n");

  // print_lonlat(dcw_lists, codeList);
  dcw_print_polygons(dcw_lists, codeList);

  return 0;
}
*/
