#include "module_info.h"
#include "mpmo_color.h"
#include "util_string.h"
#include "factory.h"
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

typedef std::function<bool(const CdoModule &mod)> ModuleQuery;

bool
ModListOptions::requested(std::string const &name)
{
  return opt[name];
}

bool
ModListOptions::mod_info_requested()
{
  return (operInfoRequested || printAll || requested(s_obase) || requested(s_arbIn) || requested(s_filesOnly)
          || requested(s_onlyFirst) || requested(s_arbIn) || requested(s_noOutput));
}

bool
ModListOptions::parse_request(std::string const &requestString)
{
  auto all = true;
  const auto splitString = split_string(requestString, ",");

  if (requestString.size() > 0)
  {
    all = false;
    for (const auto &s : splitString)
    {
      auto it = Factory::get().find(s);
      if (it != Factory::get().end())
      {
        auto &module = Factory::get_module(it);
        operInfoRequested = true;
        std::cerr << s << ": " << module.toString() << "\n";
      }
      else
      {
        if (opt.find(s) != opt.end()) { opt[s] = 1; }
        else
        {
          std::cerr << "option " << s << " not found" << "\n";
          return false;
        }
      }
    }
  }
  printAll = all;

  return true;
}

std::pair<std::map<std::string, std::vector<std::string>>, std::map<std::string, std::pair<std::string, std::vector<std::string>>>>
create_help_sections(const CdoHelp &p_help)
{
  std::map<std::string, std::vector<std::string>> sections = {};
  std::map<std::string, std::vector<std::string>> operators = {};
  std::map<std::string, std::string> oper_synopsis = {};

  std::string key;
  std::string operator_name;
  bool operator_section_active = true;
  bool synopsis_active = true;

  for (auto &line : p_help)
  {
    bool is_all_caps = std::all_of(begin(line), end(line), [](char c) { return isupper(c); });
    if (line.compare("OPERATORS") == 0)
    {
      operator_section_active = true;
      synopsis_active = false;
    }
    else if (line.compare("SYNOPSIS") == 0)
    {
      key = line;
      operator_section_active = false;
      synopsis_active = true;
    }

    else if (is_all_caps == true)
    {
      key = line;
      sections[key] = std::vector<std::string>();
      operator_section_active = false;
      synopsis_active = false;
    }
    else if (operator_section_active == true)
    {
      constexpr size_t op_name_padding = 4;
      if (line.find_first_not_of(' ') == op_name_padding)
      {
        int first_space_after_name = line.find_first_of(' ', op_name_padding);
        int name_length = first_space_after_name - op_name_padding;
        operator_name = line.substr(op_name_padding, name_length);

        operators[operator_name] = std::vector<std::string>{ std::string(line) };
      }
      else { operators[operator_name].push_back(std::string(line)); }
    }
    else if (synopsis_active == true)
    {
      constexpr size_t op_name_padding = 4;
      std::size_t general_op_desc_pos = line.find("<operator>");
      if (general_op_desc_pos != std::string::npos)
      {
        std::size_t len_generic = std::string("<operator>").size();
        auto cleaned_line = std::string(line).erase(0, general_op_desc_pos + len_generic);
        sections["SYNOPSIS"] = { cleaned_line };
      }
      else if (line.find_first_not_of(' ') == op_name_padding)
      {
        int first_space_after_name = line.find_first_of("[, ", op_name_padding);
        int name_length = first_space_after_name - op_name_padding;
        operator_name = line.substr(op_name_padding, name_length);
        oper_synopsis[operator_name] = std::string(line);
      }
    }

    else { sections[key].push_back(std::string(line)); }
  }

  std::map<std::string, std::pair<std::string, std::vector<std::string>>> oper_syn_map;
  for (auto const &op : operators)
  {
    std::string syn = "";
    auto it = oper_synopsis.find(op.first);
    if (it != oper_synopsis.end())
    {
      syn = oper_synopsis[op.first];
      oper_synopsis.erase(it);
    }
    else { syn = "    " + op.first + " " + sections["SYNOPSIS"][0]; }
    oper_syn_map[op.first] = std::make_pair(syn, op.second);
  }
  for (auto const &syn : oper_synopsis) { oper_syn_map[syn.first].first = syn.second; }

  return std::make_pair(sections, oper_syn_map);
}

std::string
get_operator_description(std::string const &p_current_op_name, const CdoHelp &p_help)
{
  std::string description = "";
  if (p_help.empty()) return description;

  // search for operator section

  auto it = std::find_if(begin(p_help), end(p_help), [&](auto const &l) { return l.find("OPERATORS") != std::string::npos; });
  // if no operator section is found
  if (it == end(p_help))
  {
    std::string name_section = std::string(p_help[0]);
    it = std::find_if(begin(p_help), end(p_help),
                      [&](auto const &l) { return l.find("    " + p_current_op_name) != std::string::npos; });

    if (it != end(p_help))
    {
      name_section += *it;
      description = name_section.substr(name_section.find_first_of('-') + 2, name_section.size());
    }
  }
  else
  {
    /* don't delete, needed for new operator_help.cc!
    auto it2 = std::find_if(it + 1, end(p_help), [&](auto const &l) { return l.compare("    " + p_current_op_name) == 0; });
    if (it2 != p_help.end())
    {
      for (int i = 0; i < (int) p_help.size() - 1; ++i)
      {
        if (p_help[i].compare("    " + p_current_op_name) == 0)
        {
          std::string line = std::string(p_help[i + 1]);
          description = line.substr(line.find_first_not_of(" \t"));
          break;
        }
      }
    }
    else
    */
    {
      it = std::find_if(++it, end(p_help),
                        [&](auto const &l) { return l.find("    " + p_current_op_name + " ") != std::string::npos; });
      if (it != p_help.end())
      {
        std::string line = std::string(*it);
        auto pos = line.find("    " + p_current_op_name + " ");
        if (pos != std::string::npos)
        {
          auto op_name_start = line.find_first_not_of(" \t");

          description = line.substr(line.find_first_not_of(" \t", op_name_start + p_current_op_name.size()), line.size());
        }
      }
    }
  }

  return description;
}

// helper function for setting the spacing in operator_print_list
static std::string
get_spacing_for(int p_space, std::string const &str)
{
  std::string spacing = "";
  for (int i = str.size(); i <= p_space; ++i) spacing += " ";
  return spacing;
}

static std::string
operatorGetShortInfoString(std::string &current_op_name, const CdoModule &p_module)
{
  std::string shortInfo = current_op_name;
  int alias_index = p_module.is_alias(current_op_name);
  if (-1 != alias_index)
  {
    shortInfo += std::string(get_spacing_for(16, current_op_name) + "--> " + p_module.aliases[alias_index].original);
  }
  else if (!p_module.get_help(current_op_name).empty())
  {
    // add spaceing and saving output line to the output list
    const auto description = get_operator_description(current_op_name, p_module.get_help(current_op_name));
    shortInfo += get_spacing_for(16, current_op_name) + description;
  }
  std::string in_out_info
      = "(" + std::to_string(p_module.get_stream_in_cnt()) + "|" + std::to_string(p_module.get_stream_out_cnt()) + ")";
  shortInfo += get_spacing_for(90, shortInfo) + in_out_info;
  return shortInfo;
}

void
operator_print_list(std::function<bool(const CdoModule &)> selectionCriteria)
{
  std::vector<std::string> output_list;

  for (auto &current_op_name : Factory::get_sorted_operator_name_list())
  {
    const CdoModule &current_module = Factory::get_module(current_op_name);
    if (selectionCriteria(current_module)) { output_list.push_back(operatorGetShortInfoString(current_op_name, current_module)); }
  }
  // print generated output list
  for (std::string const &str : output_list) { std::cout << str << std::endl; }
}

void
operator_print_list(ModListOptions &p_opt)
{
  set_text_color(stderr, GREEN);

  if (p_opt.printAll == true)
  {
    operator_print_list([](const CdoModule &) { return true; });
  }
  else
  {

    ModuleQuery defaultModuleQuery = [](const CdoModule &) -> bool { return false; };
    ModuleQuery runquestDefaultModuleQuery = [](const CdoModule &) -> bool { return true; };

    // clang-format off
    ModuleQuery hasObase  = p_opt.requested(s_obase)     ? [](const CdoModule &mod) -> bool { return mod.get_stream_out_cnt() == -1;        } : defaultModuleQuery;
    ModuleQuery hasNoOut  = p_opt.requested(s_noOutput)  ? [](const CdoModule &mod) -> bool { return mod.get_stream_out_cnt() ==  0;        } : defaultModuleQuery;
    ModuleQuery hasArb    = p_opt.requested(s_arbIn)     ? [](const CdoModule &mod) -> bool { return mod.get_stream_in_cnt()  == -1;        } : defaultModuleQuery;
    ModuleQuery filesOnly = p_opt.requested(s_filesOnly) ? [](const CdoModule &mod) -> bool { return mod.get_pos_restriction() == FilesOnly; } : defaultModuleQuery;
    ModuleQuery onlyFirst = p_opt.requested(s_onlyFirst) ? [](const CdoModule &mod) -> bool { return mod.get_pos_restriction() == OnlyFirst; } : defaultModuleQuery;
    // clang-format on

    operator_print_list([&](const CdoModule &mod)
                        { return (hasObase(mod) || hasArb(mod) || hasNoOut(mod) || filesOnly(mod) || onlyFirst(mod)); });
  }

  reset_text_color(stderr);

  return;
}

std::vector<std::string>
get_no_output_operator_list()
{
  std::vector<std::string> names;
  auto &factory = Factory::get();
  for (auto &factory_entry : factory)
  {
    auto const &module = Factory::get_module(factory_entry.first);
    if (module.mode == 1 && module.constraints.streamOutCnt == 0) { names.push_back(factory_entry.first); }
  }
  std::sort(names.begin(), names.end());

  return names;
}

void
operatorPrintAll(void)
{
  int number_of_chars = 0;
  std::string tab = "   ";
  int tab_width = tab.size();
  // using a set because it sorts the operators alphabetically on its own
  std::vector<std::string> sorted_operator_names = Factory::get_sorted_operator_name_list();

  std::cout << tab;
  for (auto const &operatorName : sorted_operator_names)
  {
    if (number_of_chars > 85)
    {
      number_of_chars = tab_width;
      std::cerr << std::endl << tab;
    }

    std::cerr << " " << operatorName;
    number_of_chars += 1 + operatorName.size();
  }

  std::cerr << std::endl;
}

void
operator_print_list(bool print_no_output)
{
  std::vector<std::string> output_list = print_no_output ? get_no_output_operator_list() : Factory::get_sorted_operator_name_list();

  auto list_length = output_list.size();

  // help variables

  for (size_t out_list_idx = 0; out_list_idx < list_length; out_list_idx++)
  {
    const std::string current_op_name = output_list[out_list_idx];
    auto &current_module = Factory::get_module(current_op_name);
    if (current_module.is_alias(current_op_name) != -1)
    {
      output_list[out_list_idx] += get_spacing_for(16, current_op_name) + "--> " + Factory::get_original(current_op_name);
    }
    else if (current_module.get_help(current_op_name).empty())
    {
      // add spaceing and saving output line to the output list
      auto description = get_operator_description(current_op_name, current_module.get_help(current_op_name));
      output_list[out_list_idx] += get_spacing_for(16, current_op_name) + description;
    }
    std::string in_out_info = " (" + std::to_string(current_module.constraints.streamOutCnt) + "|"
                              + std::to_string(current_module.constraints.streamOutCnt) + ")";
    output_list[out_list_idx] += get_spacing_for(90, output_list[out_list_idx]) + in_out_info;
  }
  // print generated output list
  for (std::string const &str : output_list) { std::cout << str << std::endl; }
}
namespace Modules
{
void
print_help(Factory::OperatorMap::iterator &it)
{
  const CdoHelp &help = Factory::get_help(it);
  if (help.empty())
    std::fprintf(stderr, "No help available for this operator!\n");
  else
  {
    for (auto &line : help)
    {
      constexpr std::array<std::string_view, 9> headers
          = { "NAME", "SYNOPSIS", "DESCRIPTION", "OPERATORS", "NAMELIST", "PARAMETER", "ENVIRONMENT", "NOTE", "EXAMPLES" };

      auto useBold = (color_enabled() && std::ranges::find(headers, line) != headers.end());
      if (useBold) set_text_color(stdout, BRIGHT);
      std::cout << line << "\n";
      if (useBold) reset_text_color(stdout);
    }
  }
}

void
print_help(std::string const &p_operator_name)
{
  auto it
      = Factory::find(p_operator_name, [&p_operator_name]() { cdo_abort("%s", Factory::err_msg_oper_not_found(p_operator_name)); });
  print_help(it);
}
}  // namespace Modules
