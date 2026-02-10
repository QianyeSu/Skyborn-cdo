#ifndef OPER_ARGS_H
#define OPER_ARGS_H

#include <string>
#include <map>
#include <any>
#include <set>
#include <functional>

#include "cdo_output.h"

class OperArg
{
public:
  OperArg(std::string const &p_key, std::function<std::any(const std::string)> p_func,
          std::vector<std::string> const &p_exclusive_with)
      : key(p_key), func([p_func](const std::string l_key) { return p_func(l_key); }), exclusive_with(p_exclusive_with)
  {
  }

  const std::string key;
  bool required = false;
  std::function<std::any(const std::string)> func;
  std::vector<std::string> exclusive_with;
  OperArg() {};
};

struct Arguments
{
  std::map<std::string, OperArg> handlers;
  template <typename... OperArgTypes>
  Arguments(std::initializer_list<OperArg> &&a_args)
  {
    for (auto &a : a_args) { handlers.emplace(a.key, std::move(a)); }
  }
};

struct ArgumentHandler
{
private:
  Arguments arguments;

  std::map<const std::string, std::vector<std::string>> keyValuePairs;
  std::set<std::string> found_keys;

public:
  ArgumentHandler(Arguments &args) : arguments(args) {}
  ArgumentHandler(const Arguments &args) : arguments(args) {}
  ArgumentHandler() : arguments({}) {}

  bool check(std::string key);
  int parse(std::vector<std::string> const &argv);

  // Templates
  template <typename T>
  bool
  get(std::string const &key, T &value)
  {

    Debug(false, "getting arg %s", key);
    if (not check(key)) return false;
    std::any return_val = arguments.handlers[key].func(keyValuePairs[key][0]);
    try
      {
        value = std::any_cast<T>(return_val);
      }
    catch (std::bad_any_cast &e)
      {
        cdo_abort("Mismatch while getting argument for %s: requested type was %s, actual type was %s", key, typeid(T).name(),
                  return_val.type().name());
      }
    return true;
  }

  template <class T, template <class, class Allocator = std::allocator<T>> class V>
  bool
  get(std::string const &key, V<T> &values)
  {

    Debug(false, "getting arg %s", key);
    if (not check(key)) return false;
    for (auto &s : keyValuePairs[key])
      {
        T val;
        std::any return_val = arguments.handlers[key].func(s);
        try
          {
            val = std::any_cast<T>(return_val);
          }
        catch (std::bad_any_cast &e)
          {
            cdo_abort("Mismatch while getting argument for %s: requested type was %s, actual type was %s", key, typeid(T).name(),
                      return_val.type().name());
          }
        values.push_back(val);
      }
    return true;
  }
};

OperArg optional(std::string key, std::function<std::any(const std::string)> p_func, std::string mut_exclusive);
OperArg optional(std::string key, std::function<std::any(const std::string)> p_func, std::vector<std::string> mut_exclusive = {});
OperArg required(std::string key, std::function<std::any(const std::string)> p_func, std::string mut_exclusive);
OperArg required(std::string key, std::function<std::any(const std::string)> p_func, std::vector<std::string> mut_exclusive = {});

#endif
