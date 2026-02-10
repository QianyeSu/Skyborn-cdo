#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "cdo_magics_mapper.h"
#include "magics_template_parser.h"

#define PARAM_COUNT sizeof(mapper) / sizeof(CdoMagicsMapper)

// extern int _set_magics_parameter_value(const char *param_name, const char *param_type, const char *param_value)

int Set_magics_param_CCOLS(const char *user_name, const char *param_value);
int Reset_magics_param_CCOLS(const char *user_name);

int Set_magics_param_CLEVS(const char *user_name, const char *param_value);
int Reset_magics_param_CLEVS(const char *user_name);

int Set_magics_param_CTABLE(const char *user_name, const char *param_value);
int Reset_magics_param_CTABLE(const char *user_name);

/* Define an array of Mapper structures to sort. */

struct CdoMagicsMapper
{
  const char *cdo_name = nullptr;
  const char *magics_name = nullptr;
  const char *magics_type = nullptr;
  // Function to Update the Corresponding Magics parameters
  int (*Set_magics_param)(const char *user_name, const char *param_value){};
  // Function to Reset the Corresponding Magics parameters
  int (*Reset_magics_param)(const char *user_name){};
};

CdoMagicsMapper mapper[] =

    { { "clevs", "contour_level_list", "floatarray", &Set_magics_param_CLEVS, &Reset_magics_param_CLEVS },
      { "ccols", "contour_param_2", "intarray", &Set_magics_param_CCOLS, &Reset_magics_param_CCOLS },
      { "color_table", "contour_param_3", "intarray", &Set_magics_param_CTABLE, &Reset_magics_param_CTABLE } };

int
Set_magics_param_CCOLS(const char *user_name, const char *param_value)

{
  if (user_name == nullptr) return 1;
  printf("Setting the CCOLS magics params \n");

  _set_magics_parameter_value("contour_shade_colour_method", "string", "list");
  _set_magics_parameter_value("contour_shade_colour_list", "stringarray", param_value);
#if 0
#endif
  return 0;
}

int
Reset_magics_param_CCOLS(const char *user_name)
{
  (void) user_name;
  printf("Re-Setting the CCOLS magics params \n");
  return 0;
}

int
Set_magics_param_CLEVS(const char *user_name, const char *param_value)
{
  if (user_name == nullptr) return 1;

  _set_magics_parameter_value("contour_level_selection_type", "string", "level_list");
  _set_magics_parameter_value("contour_level_list", "floatarray", param_value);

  return 0;
}

int
Reset_magics_param_CLEVS(const char *user_name)
{
  (void) user_name;
  _set_magics_parameter_value("contour_level_selection_type", "string", "count");
  printf("Re-Setting the CLEVS magics params \n");
  return 0;
}

int
Set_magics_param_CTABLE(const char *user_name, const char *param_value)
{
  (void) param_value;

  if (user_name == nullptr) return 1;
  printf("Setting the CTABLE magics params \n");
#if 0
  _set_magics_parameter_value( "contour_level_list", "floatarray", param_value );
#endif
  return 0;
}

int
Reset_magics_param_CTABLE(const char *user_name)
{
  (void) user_name;
  printf("Re-Setting the CTABLE magics params \n");
  return 0;
}

// This is the comparison function used for sorting and searching.

int
Compare(const void *p1, const void *p2)
{
  return std::strcmp(((CdoMagicsMapper *) p1)->cdo_name, ((CdoMagicsMapper *) p2)->cdo_name);
}

// Do the lookup into the sorted array.

// int get_magics_parameter_info( const char *user_name, char **magics_name, char **magics_type )

int
get_magics_parameter_info(const char *user_name, char *param_value)
{
  static int once = 1;
  int ret_flag = 0;

  if (once)
  {
    std::qsort(mapper, PARAM_COUNT, sizeof(CdoMagicsMapper), Compare);
    once = 0;
  }

  CdoMagicsMapper target;
  target.cdo_name = (char *) user_name;
  CdoMagicsMapper *result = (CdoMagicsMapper *) bsearch(&target, mapper, PARAM_COUNT, sizeof(CdoMagicsMapper), Compare);
  if (result)
  {
    result->Set_magics_param(result->cdo_name, param_value);

    // magics_name = result->magics_name;
    // magics_type = result->magics_type;
  }
  else
  {
    // Call the Reset functions of all the features to Reset the magics params to default in the calling function
    ret_flag = 1;
  }

  return ret_flag;
}
