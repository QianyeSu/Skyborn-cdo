/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Modali Kameswarrao

*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>

#include <cstdlib>
#include "c_wrapper.h"
#include "process_int.h"
#include "cdo_options.h"
#include <mpim_grid.h>

#ifdef HAVE_LIBMAGICS

#include <magics_api.h>
#include "magics_template_parser.h"
#include "results_template_parser.h"
#include "util_string.h"
#include "printinfo.h"

#define DBG 0

/***** ADDED for handling plots with  defined lat lon min max *****/
/*** LAT_MIN,LAT_MAX, LON_MIN,LON_MAX ****/
/**** lat_min,lat_max,lon_min,lon_max ****/

/****
subpage_lower_left_latitude
subpage_lower_left_longitude
subpage_upper_right_latitude
subpage_upper_right_longitude
****/

static const char *contour_params[]
    = { "min",    "max",       "count",      "interval", "list",    "colour",  "thickness", "style",     "RGB",
        "device", "step_freq", "file_split", "lat_min",  "lat_max", "lon_min", "lon_max",   "projection" };
static int contour_param_count = sizeof(contour_params) / sizeof(char *);

static const char *shaded_params[]
    = { "min",          "max",    "count",     "interval",   "list",    "colour_min", "colour_max", "colour_table", "RGB",
        "colour_triad", "device", "step_freq", "file_split", "lat_min", "lat_max",    "lon_min",    "lon_max",      "projection" };
static int shaded_param_count = sizeof(shaded_params) / sizeof(char *);

static const char *grfill_params[]
    = { "min",          "max",        "count",   "interval",     "list",      "colour_min", "colour_max",
        "colour_table", "resolution", "RGB",     "colour_triad", "device",    "step_freq",  "file_split",
        "lat_min",      "lat_max",    "lon_min", "lon_max",      "projection" };
static int grfill_param_count = sizeof(grfill_params) / sizeof(char *);

static const char *STD_COLOUR_TABLE[] = { "red",
                                          "green",
                                          "blue",
                                          "yellow",
                                          "cyan",
                                          "magenta",
                                          "black",
                                          "avocado",
                                          "beige",
                                          "brick",
                                          "brown",
                                          "burgundy",
                                          "charcoal",
                                          "chestnut",
                                          "coral",
                                          "cream",
                                          "evergreen",
                                          "gold",
                                          "grey",
                                          "khaki",
                                          "kellygreen",
                                          "lavender",
                                          "mustard",
                                          "navy",
                                          "ochre",
                                          "olive",
                                          "peach",
                                          "pink",
                                          "rose",
                                          "rust",
                                          "sky",
                                          "tan",
                                          "tangerine",
                                          "turquoise",
                                          "violet",
                                          "reddishpurple",
                                          "purplered",
                                          "purplishred",
                                          "orangishred",
                                          "redorange",
                                          "reddishorange",
                                          "orange",
                                          "yellowishorange",
                                          "orangeyellow",
                                          "orangishyellow",
                                          "greenishyellow",
                                          "yellowgreen",
                                          "yellowishgreen",
                                          "bluishgreen",
                                          "bluegreen",
                                          "greenishblue",
                                          "purplishblue",
                                          "bluepurple",
                                          "bluishpurple",
                                          "purple",
                                          "white" };

static char **USR_COLOUR_TABLE = nullptr;

static int STD_COLOUR_COUNT = sizeof(STD_COLOUR_TABLE) / sizeof(char *);
static int USR_COLOUR_COUNT = 0;

static const char *STYLE_TABLE[] = { "SOLID", "DASH", "DOT", "CHAIN_DASH", "CHAIN_DOT" };
static int STYLE_COUNT = sizeof(STYLE_TABLE) / sizeof(char *);

static const char *DEVICE_TABLE[] = { "PS", "EPS", "PDF", "PNG", "GIF", "GIF_ANIMATION", "JPEG", "SVG", "KML" };
static int DEVICE_COUNT = sizeof(DEVICE_TABLE) / sizeof(char *);

/*char *PROJECTION_TABLE[] = { "cylindrical", "polar_stereographic",
 * "polar_north", "geos", "meteosat", "meteosat_57E", "goes_east", "lambert",
 * "EPSG3857", "goode", "collignon", "mollweide", "robinson", "bonne", "google",
 * "efas", "EPSG4326", "lambert_north_atlantic", "mercator", "cartesian",
 * "taylor", "tephigram" };
 */
/** The following projections are having some issues to be clarified with
 * Magics++ **/

static const char *PROJECTION_TABLE[] = { "cylindrical",
                                          "polar_stereographic",
                                          "polar_north",
                                          "geos",
                                          "meteosat",
                                          "meteosat_57E",
                                          "lambert",
                                          "EPSG3857",
                                          "goode",
                                          "collignon",
                                          "mollweide",
                                          "robinson",
                                          "bonne",
                                          "google",
                                          "efas",
                                          "EPSG4326",
                                          "lambert_north_atlantic",
                                          "mercator" };

static int PROJECTION_COUNT = sizeof(PROJECTION_TABLE) / sizeof(char *);

int ANIM_FLAG = 0, STEP_FREQ = 0;  // '0' for static images like jpeg,ps, etc.. , '1' for animation formats

int checkcolour(char *colour_in);
int ReadColourTable(char *filepath);
int checkstyle(char *style_in);
int checkdevice(char *device_in);
int checkprojection(char *projection_in);

// Magics default values
static int COUNT = 10, isRGB = false, THICKNESS = 1, NUM_LEVELS = 0, FILE_SPLIT = false;
static double YMIN = 1.0e+200, YMAX = -1.0e+200, INTERVAL = 8.0, RESOLUTION = 10.0;
static double *LEV_LIST = nullptr;
static double LAT_MIN = 1.0e+200, LAT_MAX = -1.e+200;
static double LON_MIN = 1.0e+200, LON_MAX = -1.e+200;
const char *COLOUR = nullptr, *COLOUR_MIN = nullptr, *COLOUR_MAX = nullptr, *STYLE = nullptr;
const char *DEVICE = nullptr, *COLOUR_TRIAD = nullptr, *PROJECTION = nullptr;

static void
init_MAGICS()
{
  setenv("MAGPLUS_QUIET", "1", 1);  // To suppress magics messages

  mag_open();
  // Some standard parameters affectng the magics environment, moved from the xml file  ** begin **
  mag_setc("page_id_line", "off");
  mag_setc("output_name_first_page_number", "off");
  if (FILE_SPLIT == true) mag_setc("output_ps_split", "on");
}

static void
quit_MAGICS()
{
  mag_close();
  if (DBG) std::fprintf(stderr, "Exiting From MAGICS\n");
}

int
checkcolour(char *colour_in)
{
  float rgb_values[3];
  char temp[256];

  auto ref = colour_in;

  if (isRGB)
  {
    if (std::strchr(colour_in, ';') == nullptr || std::strstr(colour_in, "RGB(") == nullptr)
    {
      cdo_warning("Found 'RGB=true',Specify Colour in 'RGB(r;g;b)' ( where r,g,b in [0.0,1.0] ) format!");
      return 1;
    }

    int n = std::strlen(colour_in);

    if (DBG) std::fprintf(stdout, "  count %d  original colour %s RGB %d\n", n, colour_in, isRGB);

    int i;
    for (i = 0; i < n - 1; ++i)
    {
      if (i > 3) temp[i - 4] = *colour_in;
      colour_in++;
    }
    temp[i - 4] = '\0';

    if (DBG) std::fprintf(stdout, "  count %d  modified color %s \n", (int) std::strlen(temp), temp);

    const auto splitStrings = split_string(temp, ";");

    if (splitStrings.size() != 3)
    {
      cdo_warning(" Colour specified in Improper format!");
      return 1;
    }

    for (int k = 0; k < 3; ++k) rgb_values[k] = std::stof(splitStrings[k]);

    if (rgb_values[0] + rgb_values[1] + rgb_values[2] > 3.0f || rgb_values[0] + rgb_values[1] + rgb_values[2] < 0.0f)
    {
      cdo_warning(" RGB Colour specified with Improper values!");
      return 1;
    }
  }
  else
  {
    if (std::strchr(colour_in, ';') != nullptr || std::strstr(colour_in, "RGB(") != nullptr)
    {
      cdo_warning("Found Colour with 'RGB(r;g;b)' format, set parameter RGB='true' !");
      return 1;
    }

    cstr_to_lower(colour_in);
    for (int i = 0; i < STD_COLOUR_COUNT; ++i)
    {
      if (cdo_cmpstr(STD_COLOUR_TABLE[i], colour_in)) return 0;
    }
    cdo_warning("Specified Colour not in Standard colour list, resetting to blue(default colour)!");
    return 1;
  }

  if (DBG) cdo_warning("Colour %s verified!", ref);

  return 0;
}

int
ReadColourTable(char *filepath)
{
  auto fobj = c_fopen(filepath, "r");
  if (!fobj.get())
  {
    std::fprintf(stdout, "File Not available!");
    return 1;
  }

  int num_colors = 0;
  std::fscanf(fobj.get(), "%d", &num_colors);

  if (DBG) std::fprintf(stderr, "Num Colours %d\n", num_colors);

  if (!num_colors)
  {
    cdo_warning("No colours found in File, proceeding with Standard Colour table!");
    return 1;
  }

  USR_COLOUR_COUNT = 0;
  USR_COLOUR_TABLE = (char **) std::malloc(num_colors * sizeof(char *));
  char **temp_table = (char **) std::malloc(num_colors * sizeof(char *));

  for (int i = 0; i < num_colors; ++i)
  {
    temp_table[i] = (char *) std::malloc(256 * sizeof(char));
    std::fscanf(fobj.get(), "%s", temp_table[i]);
    if (DBG) std::fprintf(stdout, "%s\n", temp_table[i]);
  }

  char orig_char = ';', rep_char = ',';
  for (int i = 0; i < num_colors; ++i)
  {
    if (DBG) std::fprintf(stdout, "%s \n", temp_table[i]);

    if (!checkcolour(temp_table[i]))
    {
      if (isRGB) cstr_replace_char(temp_table[i], orig_char, rep_char);  // replace ';' in RGB format to ','

      if (DBG) std::fprintf(stdout, "Before appending %s\n", temp_table[i]);

      USR_COLOUR_TABLE[USR_COLOUR_COUNT] = strdup(temp_table[i]);

      // std::strcpy( USR_COLOUR_TABLE[ USR_COLOUR_COUNT ], temp_table[i] );
      USR_COLOUR_COUNT++;

      if (DBG) std::fprintf(stdout, "After appending %s\n", temp_table[i]);
    }
  }

  if (USR_COLOUR_COUNT < num_colors) cdo_warning(" Discarding improper format colours and continuing!");

  for (int i = 0; i < num_colors; ++i) std::free(temp_table[i]);
  std::free(temp_table);

  return 0;
}

int
checkstyle(char *style_in)
{
  auto found = false;
  cstr_to_upper(style_in);
  for (int i = 0; i < STYLE_COUNT; ++i)
  {
    if (DBG) std::fprintf(stderr, "Input %s ref %s\n", style_in, STYLE_TABLE[i]);

    if (cdo_cmpstr(STYLE_TABLE[i], style_in))
    {
      found = true;
      STYLE = style_in;
      return 0;
    }
  }

  if (!found) cdo_warning(" Style specified with Improper value!");

  return 1;
}

int
checkdevice(char *device_in)
{
  auto found = false;
  cstr_to_upper(device_in);
  for (int i = 0; i < DEVICE_COUNT; ++i)
  {
    if (DBG) std::fprintf(stderr, "Input %s ref %s\n", device_in, DEVICE_TABLE[i]);

    if (cdo_cmpstr(DEVICE_TABLE[i], device_in))
    {
      found = true;

      DEVICE = device_in;
      if (DBG) std::fprintf(stderr, "DEVICE %s\n", DEVICE);

      if (cdo_cmpstr("GIF_ANIMATION", device_in) || cdo_cmpstr("KML", device_in))
      {
        ANIM_FLAG = 1;
        STEP_FREQ = 1;
      }
      return 0;
    }
  }

  if (!found) cdo_warning(" Device specified with Improper value!");

  return 1;
}

int
checkprojection(char *projection_in)
{
  auto found = false;

  // cstr_to_upper( projection_in );

  for (int i = 0; i < PROJECTION_COUNT; ++i)
  {
    if (DBG) std::fprintf(stderr, "Input %s ref %s\n", projection_in, PROJECTION_TABLE[i]);

    if (cdo_cmpstr(PROJECTION_TABLE[i], projection_in))
    {
      found = true;
      PROJECTION = projection_in;
      return 0;
    }
  }

  if (!found)
  {
    cdo_warning(" Projection specified with Improper value!");
    cdo_warning(" Specify one of the following:");
    cdo_warning(" cylindrical polar_stereographic polar_north geos meteosat "
                "meteosat_57E geos_east lambert EPSG3857 goode collignon "
                "mollweide robinson bonne google efas EPSG4326 "
                "lambert_north_atlantic mercator cartesian taylor tephigram");
  }

  return 1;
}
#endif

class Magplot : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Magplot",
    .operators = { { "contour", MagplotHelp }, { "shaded", MagplotHelp }, { "grfill", MagplotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Magplot> registration = RegisterEntry<Magplot>();
  int CONTOUR, SHADED, GRFILL;
  CdoStreamID streamID;

  int operatorID;

  int taxisID;
  int vlistID;

  bool lregular = false;

  size_t gridsize;

  size_t nlon;
  size_t nlat;

  int nparam;
  std::vector<std::string> pnames;

  Varray<double> grid_center_lat;
  Varray<double> grid_center_lon;

#ifdef HAVE_LIBMAGICS
  void
  verifyPlotParameters(int num_param, std::vector<std::string> const &param_names, int opID)
  {
    bool halt_flag = false;
    int param_count = 0;
    const char **params = nullptr;
    char *temp_str;
    const char orig_char = ';', rep_char = ',';

    /*
      char *contour_params[] = {"ymin","ymax","count","interval","list","colour","thickness","style"};
      char *shaded_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colour_table","step_freq"};
      char *grfill_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colour_table","resolution"};
    */

    for (int i = 0; i < num_param; ++i)
    {
      auto found = false;
      auto syntax = true;
      const auto splitStrings = split_string(param_names[i], "=");

      if (DBG) std::fprintf(stderr, "Verifying params!\n");

      if (splitStrings.size() > 1)
      {
        auto const &key = splitStrings[0];
        auto value = strdup(splitStrings[1].c_str());

        if (opID == CONTOUR)
        {
          param_count = contour_param_count;
          params = contour_params;
        }
        else if (opID == SHADED)
        {
          param_count = shaded_param_count;
          params = shaded_params;
        }
        else if (opID == GRFILL)
        {
          param_count = grfill_param_count;
          params = grfill_params;
        }

        for (int j = 0; j < param_count; ++j)
        {
          if (key == params[j])
          {
            found = true;
            if (key == "colour" || key == "style" || key == "colour_min" || key == "colour_max" || key == "RGB"
                || key == "colour_triad" || key == "device" || key == "file_split" || key == "projection")
            {
              if (string_is_float(value)) { syntax = false; }
              else
              {
                if (key == "RGB" || key == "file_split")
                {
                  temp_str = strdup(value);
                  cstr_to_lower(temp_str);
                  if (std::strcmp(temp_str, "true") && std::strcmp(temp_str, "false")) { syntax = false; }
                  else
                  {
                    if (key == "RGB")
                      isRGB = cdo_cmpstr(temp_str, "true");
                    else if (key == "file_split")
                      FILE_SPLIT = cdo_cmpstr(temp_str, "true");
                  }
                }
                else if (key == "style")
                {
                  if (checkstyle(value)) syntax = false;
                }
                else if (key == "colour" || key == "colour_min" || key == "colour_max")
                {

                  if (checkcolour(value)) { syntax = false; }
                  else
                  {
                    if (key == "colour")
                    {
                      temp_str = strdup(value);
                      if (!isRGB)
                        cstr_to_lower(temp_str);
                      else
                      {
                        cstr_to_upper(temp_str);
                        cstr_replace_char(temp_str, orig_char, rep_char);  // replace ';' in RGB format to ','
                      }
                      COLOUR = temp_str;
                      if (DBG) std::fprintf(stderr, "COLOUR %s\n", COLOUR);
                    }
                    if (key == "colour_min")
                    {
                      temp_str = strdup(value);
                      if (!isRGB)
                        cstr_to_lower(temp_str);
                      else
                      {
                        cstr_to_upper(temp_str);
                        cstr_replace_char(temp_str, orig_char, rep_char);  // replace ';' in RGB format to ','
                      }
                      COLOUR_MIN = temp_str;
                      if (DBG) std::fprintf(stderr, "COLOUR %s\n", COLOUR_MIN);
                    }
                    if (key == "colour_max")
                    {
                      temp_str = strdup(value);
                      if (!isRGB)
                        cstr_to_lower(temp_str);
                      else
                      {
                        cstr_to_upper(temp_str);
                        cstr_replace_char(temp_str, orig_char, rep_char);  // replace ';' in RGB format to ','
                      }
                      COLOUR_MAX = temp_str;
                      if (DBG) std::fprintf(stderr, "COLOUR %s\n", COLOUR_MAX);
                    }
                  }
                }
                else if (key == "device")
                {
                  if (checkdevice(value)) syntax = false;
                }
                else if (key == "colour_triad")
                {
                  temp_str = strdup(value);
                  cstr_to_upper(temp_str);
                  if (std::strcmp(temp_str, "CW") && std::strcmp(temp_str, "ACW"))
                    syntax = false;
                  else
                  {
                    if (DBG) std::fprintf(stderr, "TRIAD check  %s!\n", temp_str);
                    COLOUR_TRIAD = cdo_cmpstr(temp_str, "CW") ? "clockwise" : "anti_clockwise";
                  }
                }
                else if (key == "projection")
                {
                  if (checkprojection(value)) syntax = false;
                }
              }
            }

            if (key == "min" || key == "max" || key == "lat_min" || key == "lat_max" || key == "lon_min" || key == "lon_max"
                || key == "count" || key == "interval" || key == "thickness" || key == "resolution" || key == "step_freq")
            {
              if (!string_is_float(value))
                syntax = false;
              else
              {
                if (key == "min") YMIN = atof(value);
                if (key == "max") YMAX = atof(value);
                if (key == "count") COUNT = atoi(value);
                if (key == "interval") INTERVAL = atof(value);
                if (key == "thickness") THICKNESS = atoi(value);
                if (key == "resolution") RESOLUTION = atoi(value);
                if (key == "step_freq") STEP_FREQ = atoi(value);
                if (key == "lat_min") LAT_MIN = atof(value);
                if (key == "lat_max") LAT_MAX = atof(value);
                if (key == "lon_min") LON_MIN = atof(value);
                if (key == "lon_max") LON_MAX = atof(value);
              }
            }

            if (key == "colour_table")
            {
              auto fobj = c_fopen(value, "r");
              if (fobj.get() == nullptr)
              {
                std::fprintf(stderr, "Input Color Table File not found in specified path '%s'\n", value);
                halt_flag = true;
              }
              else { ReadColourTable(value); }
            }

            if (key == "list")
            {
              const auto splitStrings2 = split_string(value, ";");
              if (!splitStrings2.size()) { syntax = false; }
              else
              {
                for (size_t k = 0, n = splitStrings2.size(); k < n; ++k)
                {
                  if (!string_is_float(splitStrings2[k])) syntax = false;
                }
                if (syntax == true)
                {
                  NUM_LEVELS = (int) splitStrings2.size();
                  LEV_LIST = (double *) std::malloc(NUM_LEVELS * sizeof(double));
                  for (int k = 0; k < NUM_LEVELS; ++k) LEV_LIST[k] = std::stof(splitStrings2[k]);
                }
              }
            }
          } /*** if(key == params[j])  ***/
        }   /*** Loop over param count ***/

        // if (value) std::free(value); // value is use e.g. for DEVICE
      } /*** (splitStrings.size() > 1) ***/
      else { syntax = false; }

      if (!found)
      {
        halt_flag = true;
        std::fprintf(stderr, "Invalid parameter  '%s'\n", param_names[i].c_str());
      }
      if (found && syntax == false)
      {
        halt_flag = true;
        std::fprintf(stderr, "Invalid parameter specification  '%s'\n", param_names[i].c_str());
      }
    } /*** Loop over params ****/

    if (halt_flag) std::exit(0);
  }
#endif

#ifdef HAVE_LIBMAGICS
  void
  magplot(const char *plotfile, int p_operatorID, std::string const &varname, std::string const &units, long p_nlon, long p_nlat,
          Varray<double> &p_grid_center_lon, Varray<double> &p_grid_center_lat, Varray<double> &p_array, size_t numMissVals,
          double missval, int p_nparam, std::vector<std::string> const &params, std::string const &datetimeStr, bool p_lregular)

  {
    double dlon = 0, dlat = 0;
    char plotfilename[4096];
    char *titlename;
    char tempname[256];

    if (DBG)
    {
      std::fprintf(stderr, "Num params %d\n", p_nparam);

      for (int i = 0; i < p_nparam; ++i) std::fprintf(stderr, "Param %s\n", params[i].c_str());
      fflush(stderr);

      for (int i = 0; i < p_nparam; ++i)
      {
        const auto splitStrings = split_string(params[i], "=");
        auto const &key = splitStrings[0];

        if (key == "min") std::fprintf(stderr, "Min Val %g\n", YMIN);
        if (key == "max") std::fprintf(stderr, "Max Val %g\n", YMAX);
        // if (key == "resolution") std::fprintf(stderr,"RESOLUTION %g\n",RESOLUTION );
        if (key == "colour") std::fprintf(stderr, "COLOUR %s\n", COLOUR);
        if (key == "colour_min") std::fprintf(stderr, "COLOUR %s\n", COLOUR_MIN);
        if (key == "colour_max") std::fprintf(stderr, "COLOUR %s\n", COLOUR_MAX);
        if (key == "interval") std::fprintf(stderr, "INTERVAL %f\n", INTERVAL);
        if (key == "count") std::fprintf(stderr, "COUNT %d\n", COUNT);

        if (key == "list")
          for (int j = 0; j < NUM_LEVELS; ++j) std::fprintf(stderr, "LIST %f\n", LEV_LIST[j]);

        if (key == "thickness") std::fprintf(stderr, "THICKNESS %d\n", THICKNESS);
        if (key == "style") std::fprintf(stderr, "STYLE %s\n", STYLE);
        if (key == "device") std::fprintf(stderr, "DEVICE %s\n", DEVICE);
        if (key == "step_freq") std::fprintf(stderr, "STEP_FREQ %d\n", STEP_FREQ);
        if (key == "lat_min") std::fprintf(stderr, "Lat Min Val %g\n", LAT_MIN);
        if (key == "lat_max") std::fprintf(stderr, "Lat Max Val %g\n", LAT_MAX);
        if (key == "lon_min") std::fprintf(stderr, "Lon Min Val %g\n", LON_MIN);
        if (key == "lon_max") std::fprintf(stderr, "Lon Max Val %g\n", LON_MAX);
        if (key == "projection") std::fprintf(stderr, "PROJECTION %s\n", PROJECTION);
      }
    }

    if (p_lregular)
    {
      if (p_nlon > 1)
      {
        for (int i = 1; i < p_nlon; ++i) dlon += (p_grid_center_lon[i] - p_grid_center_lon[i - 1]);
        dlon /= (p_nlon - 1);
      }
      if (p_nlat > 1)
      {
        for (int i = 1; i < p_nlat; ++i) dlat += (p_grid_center_lat[p_nlon * i] - p_grid_center_lat[p_nlon * (i - 1)]);
        dlat /= (p_nlat - 1);
      }
    }

    std::snprintf(plotfilename, sizeof(plotfilename), "%s [%s] %s", varname.c_str(), units.c_str(), datetimeStr.c_str());
    titlename = strdup(plotfilename);
    std::snprintf(plotfilename, sizeof(plotfilename), "%s_%s", plotfile, varname.c_str());

    if (Options::cdoVerbose) cdo_print("p_nlon: %zu  p_nlat: %zu", p_nlon, p_nlat);
    if (Options::cdoVerbose) cdo_print("dlon: %g  dlat: %g", dlon, dlat);

    auto mm = numMissVals ? varray_min_max_mv(p_nlon * p_nlat, p_array, missval) : varray_min_max(p_nlon * p_nlat, p_array);

    if (Options::cdoVerbose) cdo_print("min: %g  max: %g", mm.min, mm.max);
    if (Options::cdoVerbose) cdo_print("input_field_organization: %s", p_lregular ? "REGULAR" : "NONREGULAR");

    mag_setc("output_name", plotfilename);
    mag_new("page");

    // Set the input data arrays to magics++

    mag_set2r("input_field", p_array.data(), p_nlon, p_nlat);

    mag_setr("input_field_suppress_below", mm.min);
    mag_setr("input_field_suppress_above", mm.max);

    if (lregular)
    {
      mag_setc("input_field_organization", "REGULAR");
      // mag_setc("input_field_organization", "GAUSSIAN");

      mag_setr("input_field_initial_latitude", p_grid_center_lat[0]);
      mag_setr("input_field_latitude_step", dlat);

      mag_setr("input_field_initial_longitude", p_grid_center_lon[0]);
      mag_setr("input_field_longitude_step", dlon);
    }
    else
    {
      mag_setc("input_field_organization", "NONREGULAR");

      mag_set2r("input_field_latitudes", p_grid_center_lat.data(), p_nlon, p_nlat);
      mag_set2r("input_field_longitudes", p_grid_center_lon.data(), p_nlon, p_nlat);
    }

    /* magics_template_parser( magics_node ); */
    /* results_template_parser(results_node, varname.c_str() ); */

    /* set up the coastline attributes */
    /* mag_setc ("map_coastline_colour", "khaki"); */
    /* mag_setc ("map_grid_colour",      "grey");  */

    // Parameters common to all operators
    if (DEVICE) mag_setc("output_format", DEVICE);

    if (PROJECTION) mag_setc("subpage_map_projection", PROJECTION);

    mag_seti("map_label_latitude_frequency", 2);
    mag_seti("map_label_longitude_frequency", 2);
    /*mag_setr ("map_label_height",0.5);*/
    mag_setr("map_label_height", 0.4);

    /* define the contouring parameters */
    if (p_operatorID == SHADED)
    {
      mag_setc("contour", "off");
      mag_setc("contour_shade", "on");
      mag_setc("contour_shade_method", "area_fill");
      mag_setc("contour_label", "off");

      if (LAT_MIN < 1.0e+200) mag_setr("subpage_lower_left_latitude", LAT_MIN);
      if (LON_MIN < 1.0e+200) mag_setr("subpage_lower_left_longitude", LON_MIN);
      if (LAT_MAX > -1.0e+200) mag_setr("subpage_upper_right_latitude", LAT_MAX);
      if (LON_MAX > -1.0e+200) mag_setr("subpage_upper_right_longitude", LON_MAX);

      if (YMIN < 1.0e+200)
      {
        mag_setr("contour_shade_min_level", YMIN);
        mag_setr("contour_min_level", YMIN);
      }

      if (YMAX > -1.0e+200)
      {
        mag_setr("contour_shade_max_level", YMAX);
        mag_setr("contour_max_level", YMAX);
      }

      if (COLOUR_MAX) mag_setc("contour_shade_max_level_colour", COLOUR_MAX);
      if (COLOUR_MIN) mag_setc("contour_shade_min_level_colour", COLOUR_MIN);

      if (is_not_equal(INTERVAL, 8.0f))
      {
        mag_setc("contour_level_selection_type", "INTERVAL");
        mag_setr("contour_interval", INTERVAL);
      }

      if (COUNT != 10)
      {
        mag_setc("contour_level_selection_type", "COUNT");
        mag_seti("contour_level_count", COUNT);
      }

      if (NUM_LEVELS)
      {
        mag_setc("contour_level_selection_type", "LEVEL_LIST");
        mag_set1r("contour_level_list", LEV_LIST, NUM_LEVELS);
      }

      if (USR_COLOUR_COUNT)
      {
        mag_setc("contour_shade_colour_method", "LIST");
        mag_set1c("contour_shade_colour_list", (const char **) USR_COLOUR_TABLE, USR_COLOUR_COUNT);
      }

      if (COLOUR_TRIAD) { mag_setc("contour_shade_colour_direction", COLOUR_TRIAD); }

      // Adjust Set The page slightly to fit the legend
      mag_setr("subpage_x_length", 24.);
      mag_setr("subpage_y_length", 30.);

      // Legend Settings
      mag_setc("legend", "on");
      mag_setc("legend_display_type", "continuous");
      mag_setc("legend_entry_plot_direction", "column");
      mag_setc("legend_box_mode", "positional");
      mag_setr("legend_box_x_position", 26.5);
      mag_setr("legend_box_y_position", 0.39);
      mag_setr("legend_box_x_length", 2.0);
      mag_setr("legend_box_y_length", 12.69);

      if (DBG)
      {
        mag_enqc("output_name", (char *) &tempname);
        std::fprintf(stderr, " SHADED Done %s!\n", tempname);
        std::fprintf(stderr, " SHADED Done!\n");
      }
    }
    else if (p_operatorID == CONTOUR)
    {
      mag_setc("contour", "on");
      mag_setc("contour_shade", "off");
      mag_setc("contour_label", "on");
      mag_setc("contour_highlight", "off");

      if (LAT_MIN < 1.0e+200) mag_setr("subpage_lower_left_latitude", LAT_MIN);
      if (LON_MIN < 1.0e+200) mag_setr("subpage_lower_left_longitude", LON_MIN);
      if (LAT_MAX > -1.0e+200) mag_setr("subpage_upper_right_latitude", LAT_MAX);
      if (LON_MAX > -1.0e+200) mag_setr("subpage_upper_right_longitude", LON_MAX);

      if (YMIN < 1.0e+200) mag_setr("contour_min_level", YMIN);
      if (YMAX > -1.0e+200) mag_setr("contour_max_level", YMAX);

      if (COLOUR) mag_setc("contour_line_colour", COLOUR);

      if (is_not_equal(INTERVAL, 8.0f))
      {
        mag_setc("contour_level_selection_type", "INTERVAL");
        mag_setr("contour_interval", INTERVAL);
      }

      if (COUNT != 10)
      {
        mag_setc("contour_level_selection_type", "COUNT");
        mag_seti("contour_level_count", COUNT);
      }

      if (NUM_LEVELS)
      {
        mag_setc("contour_level_selection_type", "LEVEL_LIST");
        mag_set1r("contour_level_list", LEV_LIST, NUM_LEVELS);
      }

      if (THICKNESS != 1) mag_seti("contour_line_thickness", THICKNESS);

      if (STYLE) mag_setc("contour_line_style", STYLE);

      // Adjust Set The page slightly to fit the legend
      mag_setr("subpage_x_length", 24.);
      mag_setr("subpage_y_length", 30.);

      if (DBG) std::fprintf(stderr, " CONTOUR Done!\n");
    }
    else if (p_operatorID == GRFILL)
    {
      mag_setc("contour", "off");
      mag_setc("contour_shade", "on");

      // mag_setc ( "contour_shade_technique", "cell_shading" );
      mag_setc("contour_shade_technique", "grid_shading");

      mag_setc("contour_shade_method", "area_fill");
      mag_setc("contour_label", "off");

      if (LAT_MIN < 1.0e+200) mag_setr("subpage_lower_left_latitude", LAT_MIN);
      if (LON_MIN < 1.0e+200) mag_setr("subpage_lower_left_longitude", LON_MIN);
      if (LAT_MAX > -1.0e+200) mag_setr("subpage_upper_right_latitude", LAT_MAX);
      if (LON_MAX > -1.0e+200) mag_setr("subpage_upper_right_longitude", LON_MAX);

      if (YMIN < 1.0e+200)
      {
        mag_setr("contour_shade_min_level", YMIN);
        mag_setr("contour_min_level", YMIN);
      }

      if (YMAX > -1.0e+200)
      {
        mag_setr("contour_shade_max_level", YMAX);
        mag_setr("contour_max_level", YMAX);
      }

      // if( YMIN < 1.0e+200  ) mag_setr( "contour_shade_min_level", YMIN );
      // if( YMAX > -1.0e+200 ) mag_setr( "contour_shade_max_level", YMAX );

      if (COLOUR_MIN) mag_setc("contour_shade_min_level_colour", COLOUR_MIN);
      if (COLOUR_MAX) mag_setc("contour_shade_max_level_colour", COLOUR_MAX);

      if (is_not_equal(INTERVAL, 8.0f))
      {
        mag_setc("contour_level_selection_type", "INTERVAL");
        mag_setr("contour_interval", INTERVAL);
      }

      if (COUNT != 10)
      {
        mag_setc("contour_level_selection_type", "COUNT");
        mag_seti("contour_level_count", COUNT);
      }

      if (NUM_LEVELS)
      {
        mag_setc("contour_level_selection_type", "LEVEL_LIST");
        mag_set1r("contour_level_list", LEV_LIST, NUM_LEVELS);
      }

      if (USR_COLOUR_COUNT)
      {
        mag_setc("contour_shade_colour_method", "LIST");
        mag_set1c("contour_shade_colour_list", (const char **) USR_COLOUR_TABLE, USR_COLOUR_COUNT);
      }
      /*
      if( is_not_equal(RESOLUTION, 10.0f) )
        mag_setr( "contour_shade_cell_resolution", RESOLUTION );
      */
      if (COLOUR_TRIAD) mag_setc("contour_shade_colour_direction", COLOUR_TRIAD);

      // Adjust Set The page slightly to fit the legend
      mag_setr("subpage_x_length", 24.);
      mag_setr("subpage_y_length", 30.);

      // Legend Settings
      mag_setc("legend", "on");
      mag_setc("legend_display_type", "continuous");
      mag_setc("legend_entry_plot_direction", "column");
      mag_setc("legend_box_mode", "positional");
      mag_setr("legend_box_x_position", 26.5);
      mag_setr("legend_box_y_position", 0.39);
      mag_setr("legend_box_x_length", 2.0);
      mag_setr("legend_box_y_length", 12.69);

      if (DBG) std::fprintf(stderr, " GrFILL Done!\n");
    }

    // plot the title text and the coastlines
    mag_cont();
    mag_coast();

    mag_set1c("text_lines", (const char **) &titlename, 1);
    mag_setc("text_colour", "black");

    /*
      mag_setr("text_font_size", 0.6);
      mag_setc("text_mode", "positional");
      mag_setr("text_box_x_position", 1.5);
      mag_setr("text_box_y_position", 16.5);
      mag_setr("text_box_x_length", 20.);
      mag_setr("text_box_y_length", 2.5);
      mag_setc("text_border", "off");
    */

    mag_setc("text_justification", "left");
    mag_text();

    if (LEV_LIST) std::free(LEV_LIST);
  }

#endif

public:
  void
  init() override
  {

#ifdef HAVE_LIBMAGICS
    nparam = cdo_operator_argc();
    pnames = cdo_get_oper_argv();

    CONTOUR = module.get_id("contour");
    SHADED = module.get_id("shaded");
    GRFILL = module.get_id("grfill");

    operatorID = cdo_operator_id();

    if (nparam)
    {
      if (DBG)
        for (int i = 0; i < nparam; ++i) std::fprintf(stderr, "Param %d is %s\n", i + 1, pnames[i].c_str());

      verifyPlotParameters(nparam, pnames, operatorID);
    }

    streamID = cdo_open_read(0);

    vlistID = cdo_stream_inq_vlist(streamID);
    taxisID = vlistInqTaxis(vlistID);

    auto gridID = vlistInqVarGrid(vlistID, 0);
    // int zaxisID = vlistInqVarZaxis(vlistID, 0);

    const auto gridtype = gridInqType(gridID);
    if (gridtype == GRID_GME) cdo_abort("GME grid unsupported!");
    if (gridtype == GRID_UNSTRUCTURED) cdo_abort("Unstructured grids unsupported!");
    if (gridtype == GRID_GENERIC) cdo_abort("Generic grids unsupported!");

    if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) lregular = true;

    if (!lregular) cdo_abort("Curvilinear grids unsupported!");

    if (gridtype != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, NeedCorners::Yes);

    gridsize = gridInqSize(gridID);
    nlon = gridInqXsize(gridID);
    nlat = gridInqYsize(gridID);

    grid_center_lat.resize(gridsize);
    grid_center_lon.resize(gridsize);
    gridInqXvals(gridID, grid_center_lon.data());
    gridInqYvals(gridID, grid_center_lat.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID, CDI_XAXIS, grid_center_lon, "grid center lon");
    cdo_grid_to_degree(gridID, CDI_YAXIS, grid_center_lat, "grid center lat");

    // HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR
    /*
       init_XML_template_parser( Filename );
       updatemagics_and_results_nodes( );
    */

    init_MAGICS();

#endif
  }
  void
  run() override
  {
#ifdef HAVE_LIBMAGICS
    Varray<double> array(gridsize);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;

      if (ANIM_FLAG)
      {
        if (numFields > 1)
        {
          cdo_warning("File has more than one variable! Animation creation not possible!");
          break;
        }
        if (tsID % STEP_FREQ)
        {
          tsID++;
          continue;
        }
      }
      else
      {
        if (STEP_FREQ)
        {
          if (tsID % STEP_FREQ)
          {
            tsID++;
            cdo_warning("NOT PLOTTING STEP %d!", tsID);
            continue;
          }
        }
        else
        {
          if (tsID)
          {
            cdo_warning("File variables have values at more than one time step! Images created for first time step!");
            cdo_warning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!");
            cdo_warning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection "
                        "operators!");
            break;
          }
        }
      }

      const auto datetimeStr = datetime_to_string(taxisInqVdatetime(taxisID));
      if (DBG) std::fprintf(stderr, "Date/Time %s\n", datetimeStr.c_str());

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID);
        size_t numMissVals;
        cdo_read_field(streamID, array.data(), &numMissVals);

        const auto missval = vlistInqVarMissval(vlistID, varID);
        if (numMissVals) cdo_set_nan(missval, gridsize, array.data());

        auto varname = cdo::inq_var_name(vlistID, varID);
        auto units = cdo::inq_var_units(vlistID, varID);

        if (operatorID == SHADED || operatorID == CONTOUR || operatorID == GRFILL)
        {
          if (DBG)
          {
            // clang-format off
            if      (operatorID == SHADED)  std::fprintf(stderr, " Creating SHADED PLOT for %s\n", varname.c_str());
            else if (operatorID == CONTOUR) std::fprintf(stderr, " Creating CONTOUR PLOT for %s\n", varname.c_str());
            else if (operatorID == GRFILL)  std::fprintf(stderr, " Creating GRFILL PLOT for %s\n", varname.c_str());
            // clang-format on
          }

          if (DBG) std::fprintf(stderr, "Plot %d\n", varID);
          magplot(cdo_get_stream_name(1), operatorID, varname, units, nlon, nlat, grid_center_lon, grid_center_lat, array,
                  numMissVals, missval, nparam, pnames, datetimeStr, lregular);
        }
        else
          std::fprintf(stderr, "operator not implemented\n");
      }

      if (DBG) std::fprintf(stderr, "TimeStep %d\n", tsID);

      tsID++;
      /*
      if ( !STEP_FREQ  && tsID )
        {
          cdo_warning("File variables have values at more than one time step! Images created for first time step!");
          cdo_warning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!");
          cdo_warning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection
      operators!"); break;
        }
      else
        {
           tsID++;
           if (DBG) std::fprintf(stderr, "TimeStep %d\n", tsID);
        }
      */
    }

    if (ANIM_FLAG && FILE_SPLIT) cdo_warning("File split parameter ignored!");

#else

    cdo_abort("MAGICS support not compiled in!");

#endif
  }
  void
  close() override
  {
#ifdef HAVE_LIBMAGICS

    quit_MAGICS();

    cdo_stream_close(streamID);

    // quit_XML_template_parser( );
#endif
  }
};
