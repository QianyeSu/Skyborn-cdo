#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_output.h"
#include "compare.h"
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#define DBG_MSG 0

// extern int get_magics_parameter_info( const char *user_name, char **magics_name, char **magics_type );

#ifdef HAVE_LIBXML2
extern int get_magics_parameter_info(const char *user_name, char *param_value);
#endif

// Recursive function that sets the results parameters from the XML structure
#ifdef HAVE_LIBXML2
int
results_template_parser(void *node, const char *varname)
{
  xmlNode *a_node = (xmlNode *) node;
  xmlNode *cur_node = nullptr;
  xmlAttrPtr attr = nullptr;
  xmlChar *param_value;

  if (a_node == nullptr) return 1;

  if (cdo_cmpstr((const char *) a_node->name, "results"))
  {
    const char *value = (const char *) xmlGetProp(a_node, (const xmlChar *) "version");

    if (value)
    {
      if (DBG_MSG) printf("Version %s \n", value);

      if (atof(value) > 3.0f) return 1;
    }
  }

  for (cur_node = a_node->children; cur_node; cur_node = cur_node->next)
  {
#if 0
      xmlChar *param_name = nullptr;
      char *param_type = nullptr;
#endif
    param_value = nullptr;

    if (cur_node->type == XML_ELEMENT_NODE)
    {

      if (DBG_MSG) printf("Node Name: %s \n", cur_node->name);

      if (cur_node->properties == nullptr)
      {
        if (cur_node->children == nullptr) printf("NO ATTRIBUTES!!!\n");
      }
      else
      {
        // Loop Over the attributes and get the corresponding  Magics Parameter name and type, set the value
#if 0
	      printf( "Finding varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif

        if (std::strcmp(varname, (const char *) xmlGetProp(cur_node, (xmlChar *) "name")) == 0)
        {
#if 0
	          printf( "Found varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif
          for (attr = cur_node->properties; attr; attr = attr->next)
          {
            if (attr != nullptr)
            {
              param_value = xmlNodeGetContent(attr->children);

              // if( !get_magics_parameter_info( attr->name, &magics_param_name, &param_type ) )
              if (!get_magics_parameter_info((const char *) attr->name, (char *) param_value))
              {
#if 0
			      printf("Done corresponding Magics Parameter found!\n");
			      printf("Setting corresponding Magics Parameter %s and type %s!\n",magics_param_name, param_type );
	  		      fprintf(stderr, "param_value: %s\n", param_value);
			      _set_magics_parameter_value( magics_param_name, param_type, param_value );
#endif
              }
              else
              {
#if 0
			      printf("No corresponding Magics Parameter found!\n");
#endif
              }
            }
          }

          break;
        }
        else
        {
          fprintf(stderr, "Var Name not matching resetting Magics Params!\n");
          // Call the Reset functions of all the features to Reset the magics params to default
        }
      }
    }
  }

  return 0;
}
#else

int
results_template_parser(void *node, const char *varname)
{
  (void) node;
  (void) varname;
  cdo_abort("XML2 support not compiled in!");
  return 0;
}
#endif  // HAVE_LIBXML2
