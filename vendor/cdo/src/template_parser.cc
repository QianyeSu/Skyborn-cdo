#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_output.h"
#include "compare.h"
#include "template_parser.h"
#include "magics_template_parser.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
static xmlNode *root_node;
static xmlDoc *param_doc;
#endif

#define DBG_MSG 0

void *magics_node = nullptr;
void *results_node = nullptr;

#ifdef HAVE_LIBXML2
int
init_XML_template_parser(char *Filename)
{
  param_doc = xmlReadFile(Filename, nullptr, 0);
  if (param_doc == nullptr)
  {
    printf("Error: Could not parse the file \"%s\"\n", Filename);
    return (1);
  }
  else
  {
    fprintf(stderr, "XML file %s being parsed \n", Filename);
    root_node = xmlDocGetRootElement(param_doc);
  }

  return 0;
}
#else
int
init_XML_template_parser(char *Filename)
{
  (void) Filename;
  cdo_abort("XML2 support not compiled in!");
  return -1;
}
#endif

int
updatemagics_and_results_nodes(void)
{
#ifdef HAVE_LIBXML2
  xmlNode *cur_node = nullptr;

  if (root_node == nullptr)
  {
    printf("Invalid Root Node\n");
    return 0;
  }

  for (cur_node = root_node->children; cur_node; cur_node = cur_node->next)
  {
    if (cur_node->type == XML_ELEMENT_NODE)
    {
#if DBG_MSG
      fprintf(stdout, "Node Name: %s \n", cur_node->name);
#endif
      if (cdo_cmpstr((const char *) cur_node->name, "magics"))
      {
        magics_node = (void *) cur_node;
#if DBG_MSG
        fprintf(stdout, "Node Name: %s \n", cur_node->name);
#endif
      }

      if (cdo_cmpstr((const char *) cur_node->name, "results"))
      {
        results_node = (void *) cur_node;
#if DBG_MSG
        fprintf(stdout, "Node Name: %s \n", cur_node->name);
#endif
      }
    }
  }
#else

  cdo_abort("XML2 support not compiled in!");

#endif

  return 0;
}

int
quit_XML_template_parser(void)
{
#ifdef HAVE_LIBXML2
  xmlFreeDoc(param_doc);
  xmlCleanupParser();
  if (param_doc == nullptr) printf("Cleaned XML parser\n");
#if DBG_MSG
  fprintf(stdout, "Cleaned XML parser\n");
#endif
#else

  cdo_abort("XML2 support not compiled in!");

#endif

  return 0;
}
