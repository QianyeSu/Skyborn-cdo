#ifndef CDO_NODE_ATTACH_EXCEPTION_H
#define CDO_NODE_ATTACH_EXCEPTION_H

#include <vector>
#include <string>
#include <memory>
#include "node.h"
#include "cdo_exception.h"

struct NodeAttachException : public CdoException
{
  std::vector<std::string>::const_iterator iter;
  NodeAttachException(std::shared_ptr<Node> p_node, std::string const &p_msg, std::string const &p_file = "?",
                      std::string const &p_line = "?")
      : CdoException(p_msg, p_file, p_line), iter(p_node->iter){};
  NodeAttachException(std::vector<std::string>::const_iterator p_iter, std::string const &p_msg, std::string const &p_file = "?",
                      std::string const &p_line = "?")
      : CdoException(p_msg, p_file, p_line), iter(p_iter){};
};

#endif
