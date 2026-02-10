#include <functional>
#include <vector>
#include <stack>

#include "node.h"
#include "parser.h"
#include "cdo_syntax_error.h"
#include "cdo_node_attach_exception.h"
#include "factory.h"
#include "cdo_output.h"

namespace Parser
{
#define debug_parser(...) Debug(PARSER, Util::tab() + " " + __VA_ARGS__)
#define THROW(type, iter, msg) throw type(iter, msg, __FILE__, std::to_string(__LINE__))

class Parser;
using ARGV_ITERATOR = std::vector<std::string>::const_iterator;
using PARSER_STACK = std::stack<Parser>;

static int staticParserCounter = 0;

namespace Util
{

void
extract_name_and_argument(std::string const &command, std::string &operatorName, std::string &operatorArgument)
{
  constexpr char delimiter = ',';

  const size_t start = (command[0] == '-') ? 1 : 0;

  size_t len = command.find(delimiter);
  if (len == std::string::npos)
  {
    len = command.size();
    operatorArgument = "";
  }
  else { operatorArgument = command.substr(len + 1, std::string::npos); }

  operatorName = command.substr(start, len - start);
}

std::string
tab()
{
  return "|" + std::string(staticParserCounter, '\t');
}

std::vector<std::string>
generate_tokens(std::string const &p_oper)
{
  std::vector<std::string> result = {};

  auto end = p_oper.find(' ');
  auto start = 0;
  while (end != std::string::npos)
  {
    auto oper = p_oper.substr(start, end - start);
    result.push_back(oper);
    start = end + 1;
    end = p_oper.find(' ', start);
  }
  auto oper = p_oper.substr(start, end - start);
  result.push_back(oper);

  return result;
}

std::string
result_to_string(std::vector<std::shared_ptr<Node>> p_roots, std::string p_text)
{
  for (auto &x : p_roots) { p_text += x->to_string() + " | "; }
  return Green(p_text);
}

template <typename T>
void
append(std::shared_ptr<T> &parent, std::shared_ptr<T> &child)
{
  try
  {
    parent->append(child);
  }
  catch (NodeAttachException &e)
  {
    THROW(InternalCdoSyntaxError, e.iter, e.what());
  }
}

}  // namespace Util

// Factory Funcs!

// Factory Funcs!
static std::shared_ptr<Node>
create_operator_node(ARGV_ITERATOR &p_curentArgument, bool is_first = false)
{
  debug_parser("Creating new operator node: %s", *p_curentArgument);

  std::string operatorName = "";
  std::string operatorArguments = "";

  Util::extract_name_and_argument(*p_curentArgument, operatorName, operatorArguments);

  auto moduleIterator = Factory::find(
      operatorName, [&]() { THROW(InternalCdoSyntaxError, p_curentArgument, Factory::err_msg_oper_not_found(operatorName)); });

  auto mod = Factory::get_module(moduleIterator);
  auto newNode = std::make_shared<Node>(p_curentArgument, operatorName, operatorArguments, mod.constraints);
  return newNode;
}
static std::shared_ptr<Node>
create_first_operator_node(ARGV_ITERATOR &p_curentArgument)
{
  return create_operator_node(p_curentArgument, true);
}

static std::shared_ptr<Node>
create_file_node(ARGV_ITERATOR &p_curentArgument)
{
  debug_parser("Creating new file node, %s", *p_curentArgument);
  auto newNode = std::make_shared<Node>(p_curentArgument, Node::NodeType::INFILE);
  return newNode;
}

static std::shared_ptr<Node>
create_node(ARGV_ITERATOR &p_curentArgument)
{
  return ((*p_curentArgument)[0] == '-') ? create_operator_node(p_curentArgument) : create_file_node(p_curentArgument);
}

class Parser
{
private:
  // TODO comment why vector
  std::vector<std::shared_ptr<Node>> stack = {};
  std::vector<std::shared_ptr<Node>> roots = {};

public:
  ARGV_ITERATOR start;
  int cntVariableInputs = 0;
  bool handlingApply = false;

  // Regualr Parser
  Parser(ARGV_ITERATOR &p_start) : start(p_start) {}
  // Apply Constructor
  Parser(ARGV_ITERATOR &p_start, bool p_handlingApply) : start(p_start), handlingApply(p_handlingApply) {}
  // TODO: comment
  Parser(std::shared_ptr<Node> first_operator, ARGV_ITERATOR &p_start) : stack({ first_operator }), start(p_start)
  {
    if (first_operator->constraints.streamInCnt == -1) { cntVariableInputs++; }
  }

  std::string
  to_string()
  {
    return "    roots: " + Yellow(Util::result_to_string(roots, "")) + ", stack: " + Yellow(Util::result_to_string(stack, ""));
  }
  void
  add_root(std::shared_ptr<Node> p_root)
  {
    roots.push_back(p_root);
  }
  std::vector<std::shared_ptr<Node>> &
  get_roots()
  {
    return roots;
  }
  void
  push(std::shared_ptr<Node> &node)
  {
    bool undecidable_from_here_on
        = !stack.empty() && (stack.back()->constraints.streamInCnt == -1) && (stack.back()->children.size() > 1);
    if (cntVariableInputs > 1 && undecidable_from_here_on) { THROW(InternalCdoSyntaxError, node->iter, errmsg_multiple_variable); }
    if (node->constraints.streamInCnt == -1) { cntVariableInputs++; }
    debug_parser("pushing new node: %s", node->oper);
    stack.push_back(node);
  }
  void
  pop()
  {
    auto const &node = stack.back();
    debug_parser("poping node: %s", node->oper);
    if (node->constraints.streamInCnt == -1) { cntVariableInputs--; }
    debug_parser("poping node: %s", node->oper);
    stack.pop_back();
  }

  bool
  finished()
  {
    for (auto &n : stack)
    {
      if (n->constraints.streamInCnt == -1)
      {
        if (n->children.size() == 0) { return false; }
      }
      else
      {
        if ((size_t) n->constraints.streamInCnt != n->children.size()) { return false; }
      }
    }
    return true;
  }

  bool
  empty()
  {
    return stack.empty();
  }
  size_t
  size()
  {
    return stack.size();
  }

  std::shared_ptr<Node> &
  top()
  {
    return stack.back();
  }
};

static void
pop_parser(PARSER_STACK &stack)
{
  debug_parser("%s", Yellow("Poping Parser: \n" + stack.top().to_string()));
  stack.pop();
}

static Parser
Apply(ARGV_ITERATOR &cur_arg)
{
  debug_parser("creating apply subgroup");
  return Parser(cur_arg, true);
}

static Parser
Subgroup(ARGV_ITERATOR &cur_arg)
{
  debug_parser("creating normal subgroup");
  return Parser(cur_arg, false);
}
// Functions for each case TODO: better comment
static void
handle_node(Parser &parser, ARGV_ITERATOR &cur_arg)
{
  debug_parser("handling Node %s", *(cur_arg));
  auto node = create_node(cur_arg);
  if (!parser.empty())
  {
    std::shared_ptr<Node> &parent = parser.top();
    debug_parser("adding %s as leaf to %s", node->oper, parent->oper);
    Util::append(parent, node);
    if (parent->is_done()) parser.pop();
  }
  else
  {
    debug_parser("stack empty: adding to root: %s", Yellow(node->oper));
    parser.add_root(node);
  }
  parser.push(node);

  debug_parser("removing finished from stack");
  while (!parser.empty() && parser.top()->is_done()) { parser.pop(); }
}
/*  triggered on ':' */
static void
handle_apply(ARGV_ITERATOR &p_cur_arg, PARSER_STACK &stack)
{
  debug_parser("handling apply");
  if (stack.top().get_roots().empty()) { THROW(InternalCdoSyntaxError, --(p_cur_arg), errmsg_apply_no_inputs); }
  stack.push(Apply(p_cur_arg));
  staticParserCounter = stack.size();
}

static void
handle_apply_end(PARSER_STACK &p_parser_stack, ARGV_ITERATOR &p_cur_arg)
{
  /* Developers NOTE:
   * When we are here the two groups that apply needs to work should be already done.
   * That means that the position n and n -1 on the stach should be the two groups.
   */
  debug_parser(Red("handling apply end"));
  auto roots = p_parser_stack.top().get_roots();
  if (roots.empty()) { THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_apply_missing_argument); }

  pop_parser(p_parser_stack);
  if (p_parser_stack.top().get_roots().empty()) { THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_apply_missing_argument); }
  auto &to_be_applied = p_parser_stack.top().get_roots()[0];
  debug_parser("to_be_applied %s:", Util::result_to_string({ to_be_applied }, "to_be_applied: "));

  // creating new subgroup that later is used to return results from apply
  auto result = Subgroup(p_cur_arg);

  if (to_be_applied->constraints.streamInCnt != 1)
  {
    THROW(InternalCdoSyntaxError, to_be_applied->iter, errmsg_only_1_to_1_operators);
  }

  for (auto &r : roots)
  {
    debug_parser("copy: %s", to_be_applied->oper);
    auto new_root = to_be_applied->copy();
    debug_parser("add: %s to %s", r->oper, to_be_applied->oper);
    new_root->add_leaf(r);
    debug_parser("adding new root %s to root", new_root->oper);
    result.add_root(new_root);
  }

  debug_parser("%s", p_parser_stack.top().to_string());
  debug_parser("%s", Util::result_to_string(result.get_roots(), "result: "));
  pop_parser(p_parser_stack);
  p_parser_stack.push(result);
  staticParserCounter = p_parser_stack.size();
}

/*  triggered on '-apply,' */
/* This function works like this because the two apply versions only really differ in their syntax.
 * Both need two groups where the first group gets copied as many times as the second has roots.
 * For the older syntax we need to get the first group from the argument which we do with this function.
 * We then just push that group on the stack as if it was part of the rest of argv.
 */
static void
handle_old_apply(ARGV_ITERATOR &p_cur_arg, ARGV_ITERATOR end, PARSER_STACK &parser_stack)
{
  if (p_cur_arg + 1 == end || (*(p_cur_arg + 1))[0] != '[')
  {
    THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_apply_requires_bracket);
  }
  debug_parser("%s", Red("Handling old Apply"));
  std::string currentArgv = (*p_cur_arg);
  const auto pos = currentArgv.find(',');
  if (pos == std::string::npos) { THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_apply_missing_argument); }

  auto parameter = currentArgv.substr(pos + 1);
  auto tokens = Util::generate_tokens(parameter);
  ARGV_ITERATOR iterCur = tokens.begin();
  ARGV_ITERATOR iterEnd = tokens.end();

  parser_stack.push(Subgroup(p_cur_arg));
  staticParserCounter = parser_stack.size();
  auto &parser = parser_stack.top();
  while (iterCur != iterEnd)
  {
    try
    {
      auto node = create_operator_node(iterCur);
      parser.add_root(node);  // TODO: why root???
    }
    catch (InternalCdoSyntaxError &e)
    {
      THROW(InternalCdoSyntaxError, p_cur_arg, e.what());
    }
    iterCur++;
  }
  debug_parser("Result tokenizer: %s", Yellow(parser.to_string()));
  handle_apply(++p_cur_arg, parser_stack);
}

/* triggered on '[' */
static void
handle_sub_group(ARGV_ITERATOR &p_cur_arg, PARSER_STACK &stack)
{
  debug_parser("handling sub group start");
  stack.push(Subgroup(p_cur_arg));
  staticParserCounter = stack.size();
}

static void
handle_sub_group_end(ARGV_ITERATOR &p_cur_arg, PARSER_STACK &p_parser_stack)
{
  if (p_parser_stack.size() == 1) { THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_missing_sub_group); }

  if (p_parser_stack.top().handlingApply)
  {
    staticParserCounter++;
    handle_apply_end(p_parser_stack, p_cur_arg);
    staticParserCounter--;
    /* Since handle appy end takes the two needed groups of the stack and replaces them with its results we can
     * simply carry on after executing the required steps for the apply handling
     * */
  }

  debug_parser(Red("handling sub group end"));
  auto finished_group = p_parser_stack.top();
  pop_parser(p_parser_stack);

  auto &cur_parser = p_parser_stack.top();
  if (finished_group.get_roots().empty()) { THROW(InternalCdoSyntaxError, p_cur_arg, errmsg_empty_subgroup); }
  if (cur_parser.empty())
  {
    // passing on all root to roots of next parser in case of variable input operators, apply and too many brackets
    // This also allows to ignore double brackets e.g -merge [ [ -topo -topo ] ]
    for (auto const &node : finished_group.get_roots()) { cur_parser.add_root(node); }
  }
  else
  {
    debug_parser("adding to %s", cur_parser.top()->oper);

    if (cur_parser.top()->children.size() != 0) { THROW(InternalCdoSyntaxError, cur_parser.top()->iter++, errmsg_mixed_input); }
    for (auto &node : finished_group.get_roots())
    {
      if (node->has_required_inputs()) { THROW(InternalCdoSyntaxError, node->iter, errmsg_missing_inputs); }
      else { Util::append(cur_parser.top(), node); }
    }

    cur_parser.cntVariableInputs--;
    cur_parser.pop();
  }
  staticParserCounter = p_parser_stack.size();
}

static void
iterate(PARSER_STACK &parser_stack, ARGV_ITERATOR &cur_arg, const ARGV_ITERATOR &end)
{
  for (auto &arg = cur_arg; arg != end; arg++)
  {
    debug_parser("current arg: %s", Green(*arg));
    if (*arg == "[") { handle_sub_group(arg, parser_stack); }
    else if (*arg == "]") { handle_sub_group_end(arg, parser_stack); }
    else if (*arg == ":") { handle_apply(arg, parser_stack); }
    else if ((*arg).find("-apply,") == 0) { handle_old_apply(arg, end, parser_stack); }
    else if ((*arg).find("-apply") == 0) { THROW(InternalCdoSyntaxError, arg, errmsg_apply_missing_argument); }
    else { handle_node(parser_stack.top(), cur_arg); }
  }
}

std::vector<std::shared_ptr<Node>>
run(std::vector<std::string> &p_argv)
{
  std::shared_ptr<Node> first_operator;
  std::stack<Parser> parser_stack;

  ARGV_ITERATOR cur_arg = p_argv.begin();
  ARGV_ITERATOR end = p_argv.end();

  std::string first_argv = *cur_arg;
  if (first_argv.starts_with("-apply,")) { THROW(InternalCdoSyntaxError, cur_arg, errmsg_apply_in_first_pos); }

  first_operator = create_first_operator_node(cur_arg);

  // TODO remove abs replace with numOut < 0 : 1 else numOut
  int numOut = first_operator->numOut();
  if (numOut < 0) numOut = 1;
  ARGV_ITERATOR without_out_files = p_argv.end() - numOut;

  cur_arg++;

  if (std::distance(cur_arg, without_out_files) < 0) { THROW(InternalCdoSyntaxError, cur_arg, errmsg_missing_outputs); }

  parser_stack.push({ first_operator, cur_arg });

  // MAIN LOOP
  iterate(parser_stack, cur_arg, without_out_files);

  if (parser_stack.size() > 1) { THROW(InternalCdoSyntaxError, parser_stack.top().start, errmsg_bracket_not_closed); }

  if (parser_stack.top().get_roots().size() > 0)
  {
    const ARGV_ITERATOR &target_iterator = parser_stack.top().get_roots()[0]->iter;
    THROW(InternalCdoSyntaxError, target_iterator, errmsg_unprocessed_inputs);
  }
  if (parser_stack.top().finished() == false)
  {
    THROW(InternalCdoSyntaxError, parser_stack.top().top()->iter, errmsg_missing_inputs);
  }

  // all thats left should be files or obase for the output files
  std::vector<std::shared_ptr<Node>> node_structure = {};
  if (first_operator->constraints.streamOutCnt != 0)
  {
    while (cur_arg != end)
    {
      if (*cur_arg == "[" || *cur_arg == "]" || *cur_arg == ":" || (*cur_arg).find("-apply") == 0)
      {
        THROW(InternalCdoSyntaxError, cur_arg, errmsg_keyword_output);
      }
      else
      {
        bool isOutFile = true;
        debug_parser("creating new out file: %s", *cur_arg);
        auto outFileNode = std::make_shared<Node>(cur_arg, Node::NodeType::OUTFILE);
        node_structure.push_back(outFileNode);
        Util::append(outFileNode, first_operator);
        cur_arg++;
      }
    }
  }
  else { node_structure = { first_operator }; }

  return node_structure;
}
std::vector<std::shared_ptr<Node>>
_parse(std::vector<std::string> p_argv)
{
  std::vector<std::shared_ptr<Node>> res = {};
  try
  {
    res = run(p_argv);
    Debug("%s", Blue(Util::result_to_string(res)));
    return res;
  }
  catch (InternalCdoSyntaxError &e)
  {
    Debug("Error");
    throw CdoSyntaxError(e, p_argv);
  }
  Debug("we got here");
  return res;
}

std::vector<std::shared_ptr<Node>>
parse(std::vector<std::string> p_argv, const char *(*context)(void) )
{
  std::vector<std::shared_ptr<Node>> res = {};
  try
  {
    return _parse(p_argv);
  }
  catch (const CdoSyntaxError &e)
  {
    cdo_abort("\n%s", e.what());
  }

  catch (const MissingOutFileException &e)
  {
    std::string errLine = Red("\033[4m" + p_argv[0] + "\033[0m") + " ";
    auto prompt = std::string(context());
    for (auto it = p_argv.begin() + 1; it < p_argv.end(); it++) { errLine += *it + " "; }
    cdo_abort("%s %s", e.what());
  }
  catch (...)
  {
    cdo_abort("Unhandled exception");
  }
  return res;
}

}  // namespace Parser
