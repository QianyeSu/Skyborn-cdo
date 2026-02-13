#ifndef MODULE_INFO_H
#define MODULE_INFO_H

#include <string>
#include <map>
#include <vector>

#include "factory.h"
#include "operator_help.h"  // for CdoHelp

const std::string s_obase = "obase";
const std::string s_arbIn = "arbitrary";
const std::string s_filesOnly = "filesOnly";
const std::string s_onlyFirst = "onlyFirst";
const std::string s_noOutput = "noOutput";

struct ModListOptions
{
  bool printAll = false;
  bool operInfoRequested = false;
  std::map<const std::string, int> opt
      = { { s_obase, false }, { s_arbIn, false }, { s_filesOnly, false }, { s_onlyFirst, false }, { s_noOutput, false } };

  bool requested(std::string const &name);
  bool mod_info_requested();
  bool parse_request(std::string const &requestString);
};

std::string get_operator_description(std::string const &p_current_op_name, std::vector<std::string> const &p_help);
void operator_print_list(ModListOptions &p_modListOpt);
std::pair<std::map<std::string, std::vector<std::string>>, std::map<std::string, std::pair<std::string, std::vector<std::string>>>>
create_help_sections(const CdoHelp &p_help);

namespace Modules
{
void print_help(std::string const &p_operatorName);
void print_help(Factory::OperatorMap::iterator &it);
}  // namespace Modules

#endif
