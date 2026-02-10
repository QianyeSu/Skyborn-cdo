#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <stack>
#include "node.h"
#include "cdo_syntax_error.h"

namespace Parser
{

static std::string apply_help
    = "        This feature allows to prepend simple cdo chains to other chains\n"
      "        Apply syntax:\n"
      "            (1) [ chain : file1 file2 file ]             (Recommended Syntax)\n"
      "            (2) -apply,<chain> [ file1 file2 file_n ]    (Old Syntax) \n"
      "\n"
      "        For example the call:\n"
      "            -merge [ -select,name=topo : *.grb ]\n"
      "        would merge all grib files in the folder after selecting the variable topo from them\n"
      "\n"
      "        The example:\n"
      "            \"-merge [ -addc,1 -mulc,2 : -add file1 file2 -subc,1 file3 file4 ] out\""
      "        would result in:\n"
      "            -merge -addc,1 -mulc,2 -add file1 file2 -addc,1 -mulc,2 -subc,3 file3 -addc,2 -mulc,23 file4 out\n"
      "\n"
      "        In combination with the subgroup (see --argument_groups) feature this allows rather complex calls\n"
      "            -merge [ [ -addc,1 : *1991.grb ] -merge [ -mulc,23 : *1990.grb ] -add file3 file4 ] outfile \n";

static std::string subgroup_help
    = "        This feature allow to use multiple operators with variable number of inputs\n"
      "        Notes:\n"
      "            When a bracket is closed it is no longer possible to add aditional inputs.\n"
      "            When a bracket is closed another variable input operator can be used without brackets\n"
      "\n"
      "        Where it is normally not possible to chain multiple operators of that kind, with subgroups a arbitrary number can "
      "be chained\n"
      "            -merge -merge file1 operator file1 > error:\n"
      "                it cannot be decided which inputs belong to which operator\n"
      "            -merge [ -merge file1 operator ] file1 > success:\n"
      "\n"
      "        With the brackets it is possible to have multiple variable inputs as inputs for another variable input\n"
      "            -merge [ -merge [ *.grb ] -merge [ *.nc ] ] out\n"
      "        In combination with the apply (see --apply) feature this allows rather complex calls\n"
      "            -merge [ [ -addc,1 : *1991.grb ] -merge [ -mulc,23 : *1990.grb ] -add file3 file4 ] outfile \n";

//'Regular' parser messages
static std::string errmsg_multiple_variable
    = "Operator cannot be assigned.\n       Reason:\n         Multiple variable input operators used.\n         Use subgroups via "
      "[ ] to clarify relations (help: --argument_groups).\n";
static std::string errmsg_missing_outputs = "Missing outputs";
static std::string errmsg_missing_inputs = "Missing inputs";
static std::string errmsg_unprocessed_inputs
    = "Operator cannot be assigned.\n       Reason:\n         No Operators with missing input left.\n";
static std::string errmsg_keyword_output = "Keywords cannot be used as file names";

// Subgroup errors
static std::string errmsg_mixed_input = "Mixing of normal inputs and subgroups is not allowed";
static std::string errmsg_missing_sub_group = "Closing bracket without open subgroup";
static std::string errmsg_empty_subgroup = "Empty Subgroup";
static std::string errmsg_bracket_not_closed = "Bracket not closed";
static std::string errmsg_malformed_subgroup = "Malformed Subgroup";

// Apply error messages
static std::string errmsg_only_1_to_1_operators = "Only operators with a single in and output allowed";
static std::string errmsg_apply_missing_argument = "Missing arguments";
static std::string errmsg_apply_multiple_roots = "Apply can only process chains with a single in and out put";
static std::string errmsg_apply_requires_bracket = "Apply requires brackets";
static std::string errmsg_apply_no_inputs = "Apply content has no available free inputs";
static std::string errmsg_apply_in_first_pos = "Apply can not be in first position";

std::vector<std::shared_ptr<Node>> run(std::vector<std::string> &p_argv);
std::vector<std::shared_ptr<Node>> parse(std::vector<std::string> p_argv, const char *(*context)(void) );
std::vector<std::shared_ptr<Node>> _parse(std::vector<std::string> p_argv);

namespace Util
{
void extract_name_and_argument(std::string const &command, std::string &operatorName, std::string &operatorArgument);

std::string result_to_string(std::vector<std::shared_ptr<Node>> p_roots, std::string p_text = "returning: ");

std::string build_err_msg(std::vector<std::string> &p_argv, const std::vector<std::string>::const_iterator &iter,
                          std::string const &prompt, int cdo_abort_prompt_spacing = 10);
}  // namespace Util

struct MissingOutFileException : public std::invalid_argument
{
  explicit MissingOutFileException(std::string const &p_msg) : std::invalid_argument(p_msg) {}
};

}  // namespace Parser

#endif
