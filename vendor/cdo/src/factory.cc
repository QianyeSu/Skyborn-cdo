/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#include <map>
#include <string>
#include <algorithm>
#include <functional>

#include "factory.h"
#include "util_string.h"

namespace Factory
{
std::string
err_msg_oper_not_found(std::string const &operatorName)
{
  // Checking if the operatorname is an existing file name
  auto fp = std::fopen(operatorName.c_str(), "r");

// Sicnce std::format is incomplete, we sadly cannot use it yet
// When the time comes: please remove the ifdef and keep the std::format variant
#ifdef FINALLY_STD_FORMAT
  std::string err_msg;

  if (fp)
  {
    std::fclose(fp);
    err_msg = std::format("Operator missing, {} is a file on disk!", operatorName);
  }
  else
  {
    // Operator is no filename
    // Checking for similar operators
    err_msg = std::format("Operator >{}< not found!\n"
                          "Similar operators are:\n{}",
                          operatorName.c_str(), Factory::find_similar_operators(operatorName));
  }
  return err_msg;
#else
  char err_msg[1024]{ 0 };
  if (fp)
  {
    std::fclose(fp);
    std::snprintf(err_msg, sizeof(err_msg), "Operator missing, %s is a file on disk!", operatorName.c_str());
  }
  else
  {
    // Operator is no filename
    // Checking for similar operators
    auto similar = Factory::find_similar_operators(operatorName);
    std::snprintf(err_msg, sizeof(err_msg), "Operator >%s< not found!\nSimilar operators are:\n%s", operatorName.c_str(),
                  similar.c_str());
  }
  return std::string(err_msg);
#endif
}

OperatorMap &
get()
{
  static OperatorMap factory;
  return factory;
}

OperatorMap::iterator
find(std::string const &p_operName)
{
  return find(p_operName, [&p_operName]() { cdo_abort("Operator >%s< not found!", p_operName); });
}

OperatorMap::iterator
find(std::string const &p_operName, std::function<void()> p_onError)
{
  auto &operator_map = get();
  auto it = operator_map.find(p_operName);
  if (it == operator_map.end()) { p_onError(); }
  return it;
}

bool
exists(const std::string arg)
{
  std::string name = arg;
  if (arg[0] == '-') { name = arg.substr(1); }
  return get().find(name) != get().end();
}

const CdoModule &
get_module(std::string const &p_operName)
{
  auto it = find(p_operName);
  return it->second.module;
}
const CdoModule &
get_module(const OperatorMap::iterator &it)
{
  return it->second.module;
}

ModuleConstructor
get_constructor(std::string const &p_operName)
{
  auto it = find(p_operName);
  return it->second.constructor;
}

ModuleConstructor
get_constructor(const OperatorMap::iterator it)
{
  return it->second.constructor;
}

const CdoHelp &
get_help(std::string const &p_operName)
{
  auto it = find(p_operName, [&p_operName]() { cdo_abort("%s", err_msg_oper_not_found(p_operName)); });
  return get_help(it);
};

const CdoHelp &
get_help(OperatorMap::iterator p_it)
{
  auto &mod = get_module(p_it);
  auto &help_iter = mod.get_help(p_it->first);
  return help_iter;
}

/**
 * @param a pointer to a string/substring
 * @param b pointer to a string/substring
 * @param alen length of string a
 * @param blen length of string b
 * @retval true if a is similar to b
 * @retval false if a is not similar to b
 *
 * Recursive function for finding substrings of a operator name that match other operators.
 */
static bool
similar(const char *a, const char *b, unsigned long alen, unsigned long blen)
{
  if (alen > 2 && blen > 2 && std::strstr(b, a)) return true;

  while (*a && *b && *a == *b)
  {
    a++;
    b++;
  }
  if (!*a && !*b) return true;

  //  printf("%d %d %s %s\n", alen, blen, a, b);

  if (alen >= 2 && blen >= 1 && *a && similar(a + 1, b, alen - 2, blen - 1)) return true;

  if (alen >= 1 && blen >= 2 && *b && similar(a, b + 1, alen - 1, blen - 2)) return true;

  return false;
}

/**
 * @param original string tested for similarity to \p other
 * @param other string that \p original will be compared to
 * @retval true if original and other are similar
 * @retval false if not
 *
 * Wrapper function for #similar() to parse c++ strings to c strings
 */
static bool
similar(std::string const &original, std::string const &other)
{
  return (similar(original.c_str(), other.c_str(), original.size(), other.size()));
}

/***
 * function for finding similar operator names for the given string
 * @param operatorName operator name to find similar operators for
 * @returns A string with all found names. The string is seqmented into lines
 * with a max length of 75 characters
 */
std::string
find_similar_operators(std::string const &operatorName)
{
  std::string found_similar_operators = "";
  size_t lines = 1;
  constexpr size_t line_length = 105;

  if (operatorName != "")
  {
    // Searching for similar operator names in operator to module map
    for (auto const &str : Factory::get())
    {
      if (similar(string_to_lower(operatorName), str.first))
      {
        if (found_similar_operators.size() + str.first.size() > lines * line_length)
        {
          found_similar_operators += "\n";
          lines++;
        }
        found_similar_operators += str.first;
        found_similar_operators += " ";
      }
    }
  }

  if (found_similar_operators.size() == 0) { found_similar_operators = "(not found)"; }
  return found_similar_operators;
}

/***
 * Creates a sorted vector with all operator names and alisases excluding all modules that are marked as internal
 * @return a sorted std::vector containing all operator names and aliases
 * excluding all operators which modules are marked as internal
 */
std::vector<std::string>
get_sorted_operator_name_list()
{
  std::vector<std::string> names;

  auto &factory = Factory::get();

  for (auto const &factory_entry : factory)
  {
    auto &module = factory_entry.second.module;
    if (module.mode == 1) { names.push_back(factory_entry.first); }
  }

  std::ranges::sort(names);

  return names;
}

std::string
get_original(std::string const &operator_name)
{
  auto module = Factory::get_module(operator_name);
  auto index = module.is_alias(operator_name);
  if (index != -1) return module.aliases[index].original;
  return operator_name;
}

/***
 * Prints all operator names and their short descriptions
 * Aliases are listed and point to their original operator name.
 * If the module is not documented the description is empty
 * If a module has only one operator the short module description is listed
 * If the operator is not documented the description is empty
 */

OperatorMap::iterator
find_module(std::string const &operator_name)
{
  return Factory::find(operator_name);
}

const CdoHelp
get_module_help(std::string const &module_name)
{
  CdoHelp operator_names = {};
  auto &modules = Factory::get();

  std::string lower_name = "";
  for (auto c : module_name) { lower_name += tolower(c); }

  for (auto &registered_mod : modules)
  {
    const CdoModule &mod = registered_mod.second.module;
    std::string lower_registered_name = "";
    for (auto c : mod.name) { lower_registered_name += tolower(c); }
    if (lower_registered_name == lower_name) { return mod.operators.begin()->help; }
  }
  return operator_names;
}
};  // namespace Factory
