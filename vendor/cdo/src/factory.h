/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef FACTORY_H
#define FACTORY_H

#include <vector>
#include <map>
#include <memory>
#include <functional>

#include "cdo_module.h"
#include "cdo_output.h"
#include "process.h"

namespace Factory
{
// string -> pair(CdoModule,std::function)
typedef std::function<std::shared_ptr<Process>(int, std::string const &, std::vector<std::string> const &)> ModuleConstructor;

struct FactoryEntry
{
  const CdoModule &module;
  Factory::ModuleConstructor constructor;
  ArgumentHandler argHandlers;

  FactoryEntry(const CdoModule &mod, Factory::ModuleConstructor con, ArgumentHandler &p_argHandlers)
      : module(mod), constructor(con), argHandlers(p_argHandlers)
  {
  }
};

typedef std::map<std::string, FactoryEntry> OperatorMap;

std::string err_msg_oper_not_found(std::string const &operatorname);

std::string find_similar_operators(std::string const &operatorName);

const CdoHelp get_module_help(std::string const &module_name);

std::string get_original(std::string const &operatorName);

std::vector<std::string> get_sorted_operator_name_list();

bool exists(const std::string arg);
OperatorMap &get();  // Factory::get()

OperatorMap::iterator find_module(std::string const &operatorName);
OperatorMap::iterator find(std::string const &p_opername);
OperatorMap::iterator find(std::string const &p_opername, std::function<void()> p_onError);

const CdoModule &get_module(std::string const &p_operName);
const CdoModule &get_module(const OperatorMap::iterator &it);

ModuleConstructor get_constructor(std::string const &p_operName);
ModuleConstructor get_constructor(const OperatorMap::iterator it);

const CdoHelp &get_help(std::string const &p_operName);
const CdoHelp &get_help(OperatorMap::iterator p_it);
};  // namespace Factory

template <typename T>
struct RegisterEntry
{
  Factory::ModuleConstructor
  create_constructor(const CdoModule &mod)
  {
    return
        [&mod](int p_ID, std::string const &p_operName, std::vector<std::string> const &p_operatorArguments) -> std::shared_ptr<T> {
          Debug(FACTORY, "Creating process via factory function, %d = ID, %s = name, %s = mod_name", p_ID, p_operName, mod.name);
          auto new_process = std::make_shared<T>(p_ID, p_operName, p_operatorArguments, mod);
          return new_process;
        };
  }
  void
  register_operator(const CdoModule &mod, std::string const &p_oper_name, ArgumentHandler &arghandler)
  {
    Factory::get().insert(std::make_pair(p_oper_name, Factory::FactoryEntry(mod, create_constructor(mod), arghandler)));
  }

public:
  explicit RegisterEntry(ArgumentHandler arghandler = ArgumentHandler())
  {
    for (auto &oper : T::module.operators) { register_operator(T::module, oper.name,arghandler); }
    for (auto &alias : T::module.aliases) { register_operator(T::module, alias.alias,arghandler); }
  };
};
#endif
