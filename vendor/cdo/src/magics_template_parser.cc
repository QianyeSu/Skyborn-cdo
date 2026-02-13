#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include "cdo_output.h"
#include "compare.h"
#include "magics_template_parser.h"
#include "util_string.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#ifdef HAVE_LIBMAGICS
#include "magics_api.h"
#endif

#define DBG 0

#ifdef HAVE_LIBXML2
extern void *magics_node;
#endif

/* Recursive function that sets the Magics parameters from the XML structure */

#ifdef HAVE_LIBXML2
int
magics_template_parser(void *node)
{
  xmlNode *a_node = (xmlNode *) node;
  int param_set_flag;
  xmlNode *cur_node = nullptr;
  const char *param_name, *param_type, *param_value;

  if (a_node == nullptr) return 0;

#if 0
    std::fprintf( stdout,"Parsing the magics Node \n");
#endif

  if (cdo_cmpstr((const char *) a_node->name, "magics"))
  {
    const char *value = (const char *) xmlGetProp(a_node, (const xmlChar *) "version");

    if (value)
    {
      if (DBG) printf("Version %s \n", value);

      if (atof(value) > 3.0f) { return 1; }
    }
  }

  for (cur_node = a_node->children; cur_node; cur_node = cur_node->next)
  {
    param_name = nullptr;
    param_type = nullptr;
    param_value = nullptr;

    if (cur_node->type == XML_ELEMENT_NODE)
    {

      if (DBG) printf("Node Name: %s \n", cur_node->name);

#if 0
            std::fprintf( stdout,"Node Name: %s \n", cur_node->name );
#endif

      if (cur_node->properties == nullptr)
      {
        if (cur_node->children == nullptr) { printf("NO ATTRIBUTES!!!\n"); }
      }
      else
      {

        param_name = (const char *) xmlGetProp(cur_node, (const xmlChar *) "parameter");
        param_type = (const char *) xmlGetProp(cur_node, (const xmlChar *) "type");
        param_value = (const char *) xmlGetProp(cur_node, (const xmlChar *) "value");
#if 0
    		printf( "\t\tAttr name: %s Type: %s Value: %s \n", param_name,param_type,param_value);
#endif

        param_set_flag = _set_magics_parameter_value(param_name, param_type, param_value);

        if (param_set_flag) printf(" Error in Setting the Parameter %s\n", param_name);
      }
    }
  }

  return 0;
}
#else

int
magics_template_parser(void *node)
{
  (void) node;
  cdo_abort("XML2 support not compiled in!");
  return 0;
}

#endif

#ifdef HAVE_LIBMAGICS
int
_set_magics_parameter_value(const char *param_name, const char *param_type, const char *param_value)
#else
int
_set_magics_parameter_value(const char *, const char *, const char *)
#endif
{
  int ret_flag = 0;
#ifdef HAVE_LIBMAGICS
  std::string sep_char = ",";
  const char search_char = ';';

  if (param_name == nullptr)
  {
    ret_flag = 1;
    return ret_flag;
  }

  if (param_value == nullptr) ret_flag = 2;

  // MAGICS++ ENV RELATED PARAMETERS
  if (cdo_cmpstr(param_type, "environvar"))
  {
    if (cdo_cmpstr(param_name, "quiet_option"))
    {
      if (cdo_cmpstr(param_value, "off") || cdo_cmpstr(param_value, "OFF"))
      {
#if 0
              printf( "Quiet Option %s \n", param_value );
#endif
        if (!unsetenv("MAGPLUS_QUIET"))
        {
          if (DBG) std::fprintf(stderr, "Quiet Option %s is un-set successfully!!! \n", param_value);
        }
        else
          std::fprintf(stderr, "Quiet Option %s COULDN'T be UNSET!!!\n", param_value);
      }

      if (cdo_cmpstr(param_value, "on") || cdo_cmpstr(param_value, "ON"))
      {
#if 0
              printf( "Quiet Option %s \n", param_value );
#endif
        if (!setenv("MAGPLUS_QUIET", "1", 1))
        {
          if (DBG) std::fprintf(stderr, "Quiet Option %s is set successfully!!! \n", param_value);
        }
        else
          std::fprintf(stderr, "Quiet Option %s COULDN'T be SET!!!\n", param_value);
      }
    }
  }

  // MAGICS++ FLOAT TYPE PARAMETERS
  else if (cdo_cmpstr(param_type, "float")) { mag_setr(param_name, atof(param_value)); }

  // MAGICS++ FLOAT ARRAY  TYPE    PARAMETERS
  else if (cdo_cmpstr(param_type, "floatarray"))
  {

#if 0
      std::fprintf(stderr, "param_name : %s\tparam_value: %s\n", param_name, param_value);
#endif
    if (std::strchr(param_value, ';')) sep_char = ';';
    const auto splitStrings = split_string(param_value, sep_char);
    if (splitStrings.size())
    {
      std::vector<double> float_param_list(splitStrings.size());
      for (int i = 0, n = splitStrings.size(); i < n; ++i)
      {
#if 0
	      std::fprintf(stderr, "%d %d %s\n", i, (int)splitStrings.size(), splitStrings[i].c_str());
#endif
        float_param_list[i] = std::stod(splitStrings[i]);
      }
      mag_set1r(param_name, float_param_list.data(), (int) splitStrings.size());
    }
  }

  // MAGICS++ INT TYPE    PARAMETERS
  else if (cdo_cmpstr(param_type, "int")) { mag_seti(param_name, atoi(param_value)); }

  // MAGICS++ INT ARRAY  TYPE    PARAMETERS
  else if (cdo_cmpstr(param_type, "intarray"))
  {
    if (std::strchr(param_value, ';')) sep_char = ';';
    const auto splitStrings = split_string(param_value, sep_char);
    if (splitStrings.size())
    {
      std::vector<int> int_param_list(splitStrings.size());
      for (int i = 0, n = splitStrings.size(); i < n; ++i) { int_param_list[i] = std::stoi(splitStrings[i]); }
      mag_set1i(param_name, int_param_list.data(), (int) splitStrings.size());
    }
  }

  // MAGICS++ STRING TYPE    PARAMETERS
  else if (cdo_cmpstr(param_type, "string")) { mag_setc(param_name, param_value); }

  // MAGICS++ STRINGARRAY  TYPE    PARAMETERS
  else if (cdo_cmpstr(param_type, "stringarray"))
  {
    if (DBG) std::fprintf(stderr, "Input strarr is %s  Sep char is %s Search char is %c\n", param_value, sep_char.c_str(), search_char);
    if (std::strstr(param_value, ";")) sep_char = ';';

    if (DBG) std::fprintf(stderr, "Input strarr is %s  Sep char is %s\n", param_value, sep_char.c_str());
    const auto splitStrings = split_string(param_value, sep_char);

    if (DBG)
      std::fprintf(stderr, "Input strarr is %s split str count is %d Sep char is %s\n", param_value, (int) splitStrings.size(),
              sep_char.c_str());

    const char **split_str = (const char **) std::malloc(splitStrings.size() * sizeof(char *));
    for (size_t k = 0, n = splitStrings.size(); k < n; ++k) split_str[k] = splitStrings[k].c_str();
    mag_set1c(param_name, split_str, (int) splitStrings.size());
    std::free(split_str);
  }
  else
  {
    ret_flag = 3;
    std::fprintf(stderr, "Unknown Parameter Type\n");
  }
#endif

  return ret_flag;
}
