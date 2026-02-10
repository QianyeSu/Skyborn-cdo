#include "node.h"
#include "cdo_output.h"
#include "cdo_node_attach_exception.h"

Node::Node(std::vector<std::string>::const_iterator p_iter, std::string const &p_oper, std::string const &args,
           module_constraints p_constraints)
    : iter(p_iter), oper(p_oper), arguments(args), constraints(p_constraints), type(OPERATOR)
{
}

Node::Node(std::vector<std::string>::const_iterator p_iter, NodeType nodeType = INFILE)
    : iter(p_iter), oper(*p_iter), arguments(""),
      constraints({ (nodeType != INFILE), (nodeType == INFILE), PositionRestrictions::NoRestriction }), type(nodeType)
{
}

Node::Node(Node *node_ptr)
    : iter(node_ptr->iter), oper(node_ptr->oper), arguments(node_ptr->arguments), constraints(node_ptr->constraints),
      type(node_ptr->type)
{
}

Node::Node(int p_ncid, NodeType p_type)
    : oper("file"), arguments(""),
      constraints({ p_type != IN_MEM_BUFFER, p_type == IN_MEM_BUFFER, PositionRestrictions::NoRestriction }), type(p_type),
      ncid(p_ncid)
{
}

Node::Node(std::string const &p_operName, std::string const &p_args, module_constraints p_constraints)
    : oper(p_operName), arguments(p_args), constraints(p_constraints), type(OPERATOR)
{
}

Node::Node(std::string p_filename, NodeType nodeType)
    : oper(p_filename), arguments(""),
      constraints({ (nodeType != INFILE), (nodeType == INFILE), PositionRestrictions::NoRestriction }), type(nodeType)
{
}

bool
Node::has_missing_input()
{
  if (children.size() == 0 && constraints.streamInCnt != 0) return true;
  if (INFILE || children.size() == (size_t) constraints.streamInCnt || constraints.streamInCnt == -1) return false;
  return true;
}

// Returns True for files and operators that have the required number of inputs
// Variable input operators always return false
// 0 Input operators always return true becaus children.size() == 0
bool
Node::is_done()
{
  bool done = false;

  if (constraints.streamInCnt == -1) { done = false; }  // varibale inputs always false
  else if (type == INFILE) { done = true; }             // files always true
  else if ((int) children.size() == constraints.streamInCnt) { done = true; }
  Debug(CDO_NODE, "%s is done: %s", oper, done ? "true" : "false");

  return done;
}

void
Node::append(std::shared_ptr<Node> &node)
{
  if (type == OPERATOR and node->constraints.pos_restriction == PositionRestrictions::OnlyFirst)
  {
    throw NodeAttachException(node, errmsg_node_not_in_first_position);
  }
  Debug(CDO_NODE, "appending  %s to %s", node->oper, oper);
  if (type == INFILE && (node->type == INFILE || node->type == OUTFILE))
  {
    throw NodeAttachException(node, errmsg_node_file_to_file);
  }
  if (type == OUTFILE && is_done()) { throw NodeAttachException(node, errmsg_node_unassigned); }
  if (constraints.streamInCnt >= 0 && (int) children.size() == constraints.streamInCnt)
  {
    throw NodeAttachException(iter, errmsg_node_to_many_inputs);
  }

  if (node->numOut() == 0) { throw NodeAttachException(node, errmsg_node_no_output); }

  if (node->type == OPERATOR && get_restriction() == PositionRestrictions::FilesOnly)
  {
    throw NodeAttachException(node, errmsg_node_only_accepts_files);
  }
  children.push_back(node);
}

void
Node::append(std::vector<std::shared_ptr<Node>> &n)
{
  for (auto &x : n) append(x);
}

bool
Node::is_temporary_leaf()
{
  return (children.size() != (size_t) constraints.streamInCnt) && type == OPERATOR && (constraints.streamInCnt != 0);
}

void
Node::add_leaf(std::shared_ptr<Node> &new_node)
{
  Debug(CDO_NODE, "add_leaf of node: %s", new_node->oper);
  if (is_temporary_leaf())
  {
    Debug(CDO_NODE, "adding leaf to %s", oper);
    append(new_node);
  }
  else
  {
    for (auto &c : children) { c->add_leaf(new_node); }
  }
}

std::shared_ptr<Node>
Node::copy()
{
  auto copiedNode = std::make_shared<Node>(this);
  copiedNode->children.clear();
  for (auto &child : children) { copiedNode->children.push_back(child->copy()); }
  return copiedNode;
}

std::string
Node::to_string()
{
  std::string r = oper;
  if (!arguments.empty()) { r += "," + arguments; }
  if (children.size() > 0) r += " [";
  for (auto &c : children) { r += " " + c->to_string(); }
  if (children.size() > 0) r += " ]";
  return r;
}
