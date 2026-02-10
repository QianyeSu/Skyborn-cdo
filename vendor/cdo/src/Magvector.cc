/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Modali Kameswarrao

*/

#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>

#include "process_int.h"
#include <mpim_grid.h>
#include "printinfo.h"

#ifdef HAVE_LIBMAGICS

#include <magics_api.h>

#include "magics_template_parser.h"
#include "results_template_parser.h"
#include "util_string.h"

#define DBG 0

const char *vector_params[] = { "thin_fac", "unit_vec", "device", "step_freq" };
int vector_param_count = sizeof(vector_params) / sizeof(char *);

// Default Magics Values
double THIN_FAC = 2.0, UNIT_VEC = 25.0;
extern int ANIM_FLAG, STEP_FREQ;

int checkdevice(char *device_in);

extern const char *DEVICE;

static void
init_MAGICS()

{

  setenv("MAGPLUS_QUIET", "1", 1); /* To suppress magics messages */
  mag_open();

  /* Some standard parameters affectng the magics environment, moved from the xml file  ** begin ** */
  mag_setc("page_id_line", "off");
}

static void
quit_MAGICS()
{
  mag_close();
  if (DBG) fprintf(stdout, "Exiting From MAGICS\n");
}

#endif

#ifdef HAVE_LIBMAGICS
void
magvector(const char *plotfile, long nlon, long nlat, Varray<double> &grid_center_lon, Varray<double> &grid_center_lat,
          Varray<double> &uarray, Varray<double> &varray, int nparam, std::vector<std::string> &params,
          std::string const &datetimeStr)

{
  long i;
  double dlon = 0, dlat = 0;
  char plotfilename[4096];
  char *temp_str = nullptr;

  if (uarray.empty() && varray.empty())
  {
    fprintf(stderr, " No Velocity Components in input file, cannot creaate Vector PLOT!\n");
    return;
  }

  if (uarray.empty() || varray.empty())
  {
    fprintf(stderr, " Found only one Velocity Component in input file, cannot create Vector PLOT!\n");
    return;
  }

  if (DBG)
  {
    fprintf(stderr, "Num params %zu\n", params.size());

    for (i = 0; i < nparam; ++i) fprintf(stderr, "Param %s\n", params[i].c_str());
    fflush(stderr);
  }

  for (i = 0; i < nparam; ++i)
  {
    auto splitStrings = split_string(params[i], "=");
    auto const &key = splitStrings[0];
    auto const &value = splitStrings[1];

    if (key == "thin_fac")
    {
      THIN_FAC = std::stof(value);
      if (DBG) fprintf(stderr, "THIN FACTOR %g\n", THIN_FAC);
    }

    if (key == "unit_vec")
    {
      UNIT_VEC = std::stof(value);
      if (DBG) fprintf(stderr, "UNIT VECTOR %g\n", UNIT_VEC);
    }

    if (key == "device")
    {
      temp_str = strdup(value.c_str());
      cstr_to_upper(temp_str);
      DEVICE = temp_str;
      if (DBG) fprintf(stderr, "DEVICE %s\n", DEVICE);

      mag_setc("output_format", DEVICE);
    }

    if (key == "step_freq")
    {
      STEP_FREQ = std::stoi(value);
      if (DBG) fprintf(stderr, "STEP FREQ %d\n", STEP_FREQ);
    }
  }

  if (nlon > 1)
  {
    for (i = 1; i < nlon; ++i) dlon += (grid_center_lon[i] - grid_center_lon[i - 1]);
    dlon /= (nlon - 1);
  }

  if (nlat > 1)
  {
    for (i = 1; i < nlat; ++i) dlat += (grid_center_lat[nlon * i] - grid_center_lat[nlon * (i - 1)]);
    dlat /= (nlat - 1);
  }

  /* magics_template_parser( magics_node ); */

  /* results_template_parser(results_node, varname ); */

  std::snprintf(plotfilename, sizeof(plotfilename), "Velocity Vectors %s", datetimeStr.c_str());
  char *titlename = strdup(plotfilename);
  std::snprintf(plotfilename, sizeof(plotfilename), "%s", plotfile);

  mag_setc("output_name", plotfilename);
  mag_new("page");

  /* Set the input data */
  mag_setr("input_field_initial_latitude", grid_center_lat[0]);
  mag_setr("input_field_latitude_step", dlat);

  mag_setr("input_field_initial_longitude", grid_center_lon[0]);
  mag_setr("input_field_longitude_step", dlon);

  mag_set2r("input_wind_u_component", uarray.data(), nlon, nlat);
  mag_set2r("input_wind_v_component", varray.data(), nlon, nlat);

  mag_seti("map_label_latitude_frequency", 2);
  mag_seti("map_label_longitude_frequency", 2);
  /*mag_setr ("map_label_height",0.5);*/
  mag_setr("map_label_height", 0.4);

  // if (operatorID == VECTOR)
  {
    /* Magics functions for performing vector operation */
    /*
      mag_setc("wind_legend_only", "on" );
      mag_setc("wind_legend_text", "on" );
    */

    mag_setc("legend", "on");
    mag_setc("wind_flag_cross_boundary", "on");
    mag_seti("wind_arrow_thickness", 1);
    mag_coast();

    if (is_not_equal(THIN_FAC, 2.0f)) mag_setr("wind_thinning_factor", THIN_FAC);

    /*wind_arrow_unit_velocity */
    if (is_not_equal(UNIT_VEC, 25.0f)) mag_setr("wind_arrow_unit_velocity", UNIT_VEC);

    mag_wind();

    mag_set1c("text_lines", (const char **) &titlename, 1);
    mag_setc("text_colour", "black");
    mag_setc("text_justification", "centre");
    mag_text();
  }

  std::free(titlename);
}

void
verify_vector_parameters(int num_param, std::vector<std::string> &param_names, bool vectorMode)
{

  int i, j;
  auto halt_flag = false;
  int param_count = 0;
  const char **params = nullptr;

  // char  *vector_params[] = {"min","max","count","interval","list","colour","thickness","style","RGB"};

  for (i = 0; i < num_param; ++i)
  {
    auto found = false;
    auto syntax = true;
    auto splitStrings = split_string(param_names[i], "=");

    if (DBG) fprintf(stderr, "Verifying params!\n");

    if (splitStrings.size() > 1)
    {
      auto const &key = splitStrings[0];
      auto const &value = splitStrings[1];
      if (vectorMode)
      {
        param_count = vector_param_count;
        params = vector_params;
      }

      for (j = 0; j < param_count; ++j)
      {
        if (key == params[j])
        {
          found = true;

          if (key == "thin_fac" || key == "unit_vec" || key == "step_freq")
          {
            if (!string_is_float(value)) syntax = false;
          }

          if (key == "device")
          {
            if (string_is_float(value))
              syntax = false;
            else
            {
              if (DBG) fprintf(stderr, "Parameter value '%s'\n", value.c_str());
              char *deviceCstr = strdup(value.c_str());
              if (checkdevice(deviceCstr)) syntax = false;

              // Vector not supported in google earth format
              if (value == "KML" || value == "kml")
              {
                syntax = false;
                if (DBG) fprintf(stderr, "Parameter value '%s'\n", value.c_str());
              }
            }
          }
        }
      }
    }
    else { syntax = false; }

    if (!found)
    {
      halt_flag = true;
      fprintf(stderr, "Invalid parameter  '%s'\n", param_names[i].c_str());
    }
    if (found && !syntax)
    {
      halt_flag = true;
      fprintf(stderr, "Invalid parameter specification  '%s'\n", param_names[i].c_str());
    }
  }

  if (halt_flag) exit(0);
}
#endif

class Magvector : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Magvector",
    .operators = { { "vector", MagvectorHelp }, { "stream", MagvectorHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Magvector> registration = RegisterEntry<Magvector>(module);

  int VECTOR, STREAM;
  CdoStreamID streamID;
  int taxisID;
  int vlistID;
  int operatorID;

  size_t gridsize;

  int nlon;
  int nlat;

  int nparam;

  std::vector<std::string> pnames;

  Varray<double> uarray;
  Varray<double> varray;
  Varray<double> grid_center_lat;
  Varray<double> grid_center_lon;

  int found = 0;

public:
  void
  init() override
  {
#ifdef HAVE_LIBMAGICS
    nparam = cdo_operator_argc();
    pnames = cdo_get_oper_argv();

    VECTOR = module.get_id("vector");
    STREAM = module.get_id("stream");

    operatorID = cdo_operator_id();

    if (nparam)
    {
      if (DBG)
        for (int i = 0; i < nparam; ++i) fprintf(stderr, "Param %d is %s!\n", i + 1, pnames[i].c_str());

      verify_vector_parameters(nparam, pnames, (operatorID == VECTOR));
    }

    streamID = cdo_open_read(0);

    vlistID = cdo_stream_inq_vlist(streamID);
    taxisID = vlistInqTaxis(vlistID);

    auto gridID = vlistInqVarGrid(vlistID, 0);
    // int zaxisID = vlistInqVarZaxis(vlistID, 0);

    if (gridInqType(gridID) == GRID_GME) cdo_abort("GME grid unsupported!");
    if (gridInqType(gridID) == GRID_UNSTRUCTURED) cdo_abort("Unstructured grid unsupported!");

    if (gridInqType(gridID) != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, NeedCorners::Yes);

    gridsize = gridInqSize(gridID);

    nlon = gridInqXsize(gridID);
    nlat = gridInqYsize(gridID);
    // int nlev     = zaxisInqSize(zaxisID);

    uarray = Varray<double>(gridsize);
    varray = Varray<double>(gridsize);
    grid_center_lat = Varray<double>(gridsize);
    grid_center_lon = Varray<double>(gridsize);

    gridInqYvals(gridID, grid_center_lat.data());
    gridInqXvals(gridID, grid_center_lon.data());

    // Convert lat/lon units if required
    cdo_grid_to_degree(gridID, CDI_XAXIS, grid_center_lon, "grid center lon");
    cdo_grid_to_degree(gridID, CDI_YAXIS, grid_center_lat, "grid center lat");

    /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
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
    int tsID = 0;

    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;

      if (ANIM_FLAG)
      {
        if (tsID % STEP_FREQ)
        {
          tsID++;
          continue;
        }
      }
      else
      {
        if (!STEP_FREQ && tsID)
        {
          cdo_warning("File has values at more than one time step! Image created for first time step!!!");
          break;
        }
      }

      auto datetimeStr = datetime_to_string(taxisInqVdatetime(taxisID));

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID);

        auto varname = cdo::inq_var_name(vlistID, varID);

        if (operatorID == VECTOR)
        {
          if (varname == "var131" || varname == "u")  // U Velocity as per GRIB is var131, as per NC 'u'
          {
            if (DBG) fprintf(stderr, "Found U VEL in Varname %s\n", varname.c_str());
            size_t numMissVals;
            cdo_read_field(streamID, uarray.data(), &numMissVals);
            if (numMissVals) cdo_set_nan(vlistInqVarMissval(vlistID, varID), gridsize, uarray.data());
            found++;
          }
          if (varname == "var132" || varname == "v")  // V Velocity as per GRIB  is var132, as per NC 'v'
          {
            if (DBG) fprintf(stderr, "Found V VEL in Varname %s\n", varname.c_str());
            size_t numMissVals;
            cdo_read_field(streamID, varray.data(), &numMissVals);
            if (numMissVals) cdo_set_nan(vlistInqVarMissval(vlistID, varID), gridsize, varray.data());
            found++;
          }
          if (found == 2) break;
        }
        else if (operatorID == STREAM)
          fprintf(stderr, " Stream Operator Un-Supported!\n");
        else
          fprintf(stderr, " Operator Un-Supported!\n");
      }

      if (operatorID == VECTOR)
      {
        if (found == 2)
        {
          if (DBG) fprintf(stderr, "Found Both U & V VEL, Creating vector fields! \n");
          magvector(cdo_get_stream_name(1), nlon, nlat, grid_center_lon, grid_center_lat, uarray, varray, nparam, pnames,
                    datetimeStr);
        }
        else if (found == 1)
        {
          fprintf(stderr, "Found only one Velocity Component in input file, cannot create Vector PLOT!\n");
          break;
        }
        else if (found == 0)
        {
          fprintf(stderr, "No Velocity Components in input file, cannot create Vector PLOT!\n");
          break;
        }
      }

      tsID++;

      /*
      if( ANIM_FLAG )
        tsID++;
      else
        {
           cdo_warning("File has values at more than one time step! Image created for first time step!!!");
           if( STEP_FREQ > 1 ) cdo_warning("Step frequency parameter ignored!!!"); break;
        }
      */
    }
#else
    cdo_abort("MAGICS support not compiled in!");
#endif
  }

  void
  close() override
  {
#ifdef HAVE_LIBMAGICS
    cdo_stream_close(streamID);
    // quit_XML_template_parser();
    quit_MAGICS();
#endif
  }
};
