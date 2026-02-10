#ifndef CDO_NODE_HPP
#define CDO_NODE_HPP

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include "cdo_module.h"

static std::string errmsg_node_to_many_inputs = "To many inputs";
static std::string errmsg_node_no_output = "Operator has no output, cannot be used with pipes unless used first";
static std::string errmsg_node_unassigned = "Could not be assigned, leftover input";
static std::string errmsg_node_file_to_file = "Attempted to attach file to file";
static std::string errmsg_node_only_accepts_files = "Operator cannot be piped into an operator that takes only files";
static std::string errmsg_node_not_in_first_position = "This operator can't be combined with other operators!";

class Node
{
public:
  enum NodeType
  {
    OPERATOR,
    INFILE,
    OUTFILE,
    IN_MEM_BUFFER,
    OUT_MEM_BUFFER
  };

  std::vector<std::string>::const_iterator iter;
  const std::string oper;
  const std::string arguments;
  module_constraints constraints;
  std::vector<std::shared_ptr<Node>> children = {};

  Node(std::vector<std::string>::const_iterator p_iter, std::string const &p_operName, std::string const &p_args,
       module_constraints p_constraints);
  Node(std::string const &p_operName, std::string const &p_args, module_constraints p_constraints);
  Node(std::vector<std::string>::const_iterator p_iter, NodeType nodeType);
  Node(std::string p_filename, NodeType nodeType);
  Node(int ncid, NodeType type);

  explicit Node(Node *p_nodePtr);
  std::shared_ptr<Node> copy();

  // Ready to be returned and process
  bool has_missing_input();
  // Done in terms of beein on the stack
  bool is_done();
  bool is_temporary_leaf();
  bool is_leaf();

  void add_leaf(std::shared_ptr<Node> &p_newNode);
  void append(std::vector<std::shared_ptr<Node>> &p_node);
  void append(std::shared_ptr<Node> &p_node);

  bool
  isInFile()
  {
    return type == INFILE;
  }
  bool
  isOperator()
  {
    return type == OPERATOR;
  }
  int
  numMaxChildren()
  {
    return constraints.streamInCnt;
  }
  int
  numOut()
  {
    return constraints.streamOutCnt;
  }
  PositionRestrictions
  get_restriction()
  {
    return constraints.pos_restriction;
  }

  const NodeType type;
  int ncid = -1;
  std::string to_string();

  bool
  has_required_inputs()
  {
    Debug(CDO_NODE, "Checking required inputs");
    bool variable_done = (constraints.streamInCnt == -1 && children.size() >= 1);
    bool node_done = (is_done() || variable_done);

    bool has_required_inputs = (!isInFile() && !node_done);
    return has_required_inputs;
  }
};

#endif
