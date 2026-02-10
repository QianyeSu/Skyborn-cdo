#include "cdo_module.h"
#include "cdo_output.h"

#include <algorithm>

oper_t::oper_t() : help(default_help) {}

oper_t::oper_t(const char *_name, int _f1, int _f2, const char *_enter)
    : name(_name), f1(_f1), f2(_f2), enter(_enter), help(default_help)
{
}

oper_t::oper_t(const char *_name, int _f1, int _f2) : name(_name), f1(_f1), f2(_f2), enter(nullptr), help(default_help) {}

oper_t::oper_t(const char *_name, int _f1, int _f2, const char *_enter, const CdoHelp &p_help)
    : name(_name), f1(_f1), f2(_f2), enter(_enter), help(p_help)
{
}

oper_t::oper_t(const char *_name, int _f1, int _f2, const CdoHelp &p_help)
    : name(_name), f1(_f1), f2(_f2), enter(nullptr), help(p_help)
{
}

oper_t::oper_t(const char *_name) : name(_name), help(default_help) {}

oper_t::oper_t(const char *_name, const CdoHelp &p_help) : name(_name), enter(nullptr), help(p_help) {}

int
CdoModule::get_id(std::string const &oper) const
{
  for (size_t i = 0; i < operators.size(); i++)
    {
      if (oper == operators[i].name) { return i; }
    }
  return -1;
}

std::string
CdoModule::toString() const
{
  std::string inp = (constraints.streamInCnt >= 0) ? std::to_string(get_stream_in_cnt()) : "Arbitrary";
  std::string out = (constraints.streamOutCnt >= 0) ? std::to_string(get_stream_out_cnt()) : "Output base";
  std::string restriction = "none";

  if (get_pos_restriction() == OnlyFirst) restriction = "Can only be the first operator";
  if (get_pos_restriction() == FilesOnly)
    (restriction != Green("none")) ? restriction = Yellow("Can only use files as input.")
                                   : restriction += Yellow(", Can only use files as input");

  std::string desc = "Input: " + inp + ", Ouput: " + out + ", Restricton: " + restriction;

  return desc;
}

int
CdoModule::get_stream_in_cnt() const
{
  return constraints.streamInCnt;
}

int
CdoModule::get_stream_out_cnt() const
{
  return constraints.streamOutCnt;
}

int
CdoModule::get_number() const
{
  return number;
}

int
CdoModule::get_mode() const
{
  return mode;
}

int
CdoModule::get_pos_restriction() const
{
  return constraints.pos_restriction;
}

int
CdoModule::is_alias(std::string const &subject) const
{
  for (size_t i = 0; i < aliases.size(); i++)
    {
      if (aliases[i].alias == subject) { return i; }
    }
  return -1;
}

std::vector<oper_t>::const_iterator
CdoModule::find_operator(std::string const &p_operatorName) const
{
  auto res
      = std::find_if(begin(operators), end(operators), [&p_operatorName](const oper_t &o) { return o.name == p_operatorName; });
  return res;
}

const CdoHelp &
CdoModule::get_help(std::string const &p_operatorName) const
{
  auto oper_iter = find_operator(p_operatorName);
  if (oper_iter == operators.end()) cdo_abort("Help for %s in module %s not found", p_operatorName, name);

  return oper_iter->help;
}
