#include <string>
#include <any>
#include <functional>

#include "oper_args.h"
#include "util_string.h"

bool
ArgumentHandler::check(std::string key)
{
  if (arguments.handlers.find(key) == arguments.handlers.end())
  {
    cdo_abort("Unkown Option %s", key);
    return false;
  }
  if (keyValuePairs.find(key) == keyValuePairs.end())
  {
    if (arguments.handlers[key].required) { cdo_abort("Argument >%s< is required!", key); }
    return false;
  };
  return true;
}

int
ArgumentHandler::parse(std::vector<std::string> const &argv)
{
  /*  this function assumes input in the form of
   *  value_name1=10,2,3,4,value_name2=201,23
   **/

  size_t equalPos = argv[0].find('=');
  if (equalPos == std::string::npos)
  {
    std::fprintf(stderr, "missing '=' in key/value string: >%s<\n", argv[0].c_str());
    return -1;
  }

  for (std::string const &arg : argv)
  {
    auto current = keyValuePairs.end();
    equalPos = arg.find('=');
    if (equalPos != std::string::npos)
    {
      auto key = arg.substr(0, equalPos);
      auto success = found_keys.insert(key);
      auto current_arg = arguments.handlers[key];

      for (auto &ew : current_arg.exclusive_with)
      {
        std::cout << key << " is exclsive with: " << ew << std::endl;
        if (found_keys.find(ew) != found_keys.end())
        {
          cdo_abort("%s can not be combined with any of %s", key, cdo_argv_to_string(current_arg.exclusive_with));
        }
      }

      if (success.second == false) { cdo_abort("Error while creating argument parser: duplicate key <%s>", key); }
      keyValuePairs[key] = {};
      current = keyValuePairs.find(key);

      std::string value = Util::String::trim(arg.substr(equalPos + 1));
      if (value.empty()) { cdo_abort("%s has no value", arg); }
      {  // value vector gets moved at the end of scope

        if (!value.empty()) { current->second.push_back(value); }
      }
    }
    else { current->second.push_back(arg); }
  }

  return 0;
}

OperArg
optional(std::string const &key, std::function<std::any(const std::string)> p_func, std::string mut_exclusive)
{
  return OperArg(key, p_func, { mut_exclusive });
}

OperArg
optional(std::string const &key, std::function<std::any(const std::string)> p_func, std::vector<std::string> const &mut_exclusive)
{
  return OperArg(key, p_func, mut_exclusive);
}

OperArg
required(std::string const &key, std::function<std::any(const std::string)> p_func, std::string mut_exclusive)
{
  auto operarg = OperArg(key, p_func, { mut_exclusive });
  operarg.required = true;
  return operarg;
}

OperArg
required(std::string const &key, std::function<std::any(const std::string)> p_func, std::vector<std::string> const &mut_exclusive)
{
  auto operarg = OperArg(key, p_func, mut_exclusive);
  operarg.required = true;
  return operarg;
}
