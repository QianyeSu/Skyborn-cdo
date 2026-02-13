#ifndef CDO_GETOPT_H
#define CDO_GETOPT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <map>
#include <vector>
#include <string>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <algorithm>

#include "cdo_output.h"

static std::string DEFAULT_CATEGORY = "Options";

struct cdo_option_argument
{
  std::string description = std::string();
  bool is_optional = false;
  std::string errmsg = std::string();
  bool allow_keyword = false;
  std::string default_value = "";
  std::string value = std::string();
  std::function<void(std::string argument)> on_empty_argument = [](std::string s = "") { cdo_abort("missing argument for %s", s); };
};

struct CliOption
{
  std::function<void(std::string argument)> effect = [](std::string)
  {
    printf("ERROR effect has not been set. Called empty effect");
    std::exit(0);
  };
  bool hasArgument = false;
  bool abortOnExecute = false;
  cdo_option_argument argument = {};
  std::vector<std::string> description = {};
  std::string name = {};
  std::string shortName = {};
  bool isInternal = false;
  std::string linkedEnvironmentVariable = {};
  std::string category = DEFAULT_CATEGORY;

  CliOption *
  set_category(std::string const &p_category)
  {
    category = p_category;
    return this;
  }

  CliOption *
  on_empty_argument(const std::function<void(std::string p_argument)> p_on_empty_argument)
  {
    argument.on_empty_argument = p_on_empty_argument;
    return this;
  }

  CliOption *
  on_empty_argument(const std::function<void()> p_on_empty_argument)
  {
    argument.on_empty_argument = [p_on_empty_argument](std::string) { p_on_empty_argument(); };
    return this;
  }

  CliOption *
  add_effect(const std::function<void(std::string argument)> p_effect)
  {
    hasArgument = true;
    effect = p_effect;
    return this;
  }

  CliOption *
  add_effect(const std::function<void(void)> p_effect)
  {
    effect = [p_effect](std::string) { p_effect(); };
    return this;
  }

  CliOption *
  set_internal(bool p_is_internal)
  {
    isInternal = p_is_internal;
    return this;
  }

  CliOption *
  add_default(std::string const &p_default_value)
  {
    argument.default_value = p_default_value;
    return this;
  }

  CliOption *
  accepts_keyword(bool enable)
  {
    argument.allow_keyword = enable;
    return this;
  }

  CliOption *
  describe_argument(std::string const &desc)
  {
    hasArgument = true;
    argument.description = desc;
    return this;
  }
  CliOption *
  set_argument_optional(bool p_isOptional)
  {
    argument.is_optional = p_isOptional;
    return this;
  }

  CliOption *
  add_help(std::string desc...)
  {
    description.push_back(desc);
    return this;
  }

  template <typename... T>
  CliOption *
  add_help(std::string const &arg, T... args)
  {
    add_help(arg);
    add_help(args...);
    return this;
  }

  std::string
  get_arg()
  {
    return argument.value;
  }

  void
  execute()
  {
    if (effect == nullptr) { cdo_abort("effect not set"); }
    effect(get_arg());
  }

  CliOption *
  aborts_program(bool aborts)
  {
    abortOnExecute = aborts;
    return this;
  }
  const CliOption *shortform(char p_shortform);

  CliOption *
  add_env_var(std::string const &p_envVarName)
  {
    linkedEnvironmentVariable = p_envVarName;
    return this;
  }

  void
  evaluate_env_var()
  {
    const char *env_var = getenv(linkedEnvironmentVariable.c_str());
    if (env_var != nullptr)
    {
      if (hasArgument)
      {
        if (*env_var == '\0') { cdo_abort("Error: %s defined but has no value", linkedEnvironmentVariable); }
        effect(std::string(env_var));
      }
      else { effect(std::string()); }
    }
  }
};

class CLIOptions
{
private:
  static std::map<std::string, std::unique_ptr<CliOption>> optionMap;
  static std::map<std::string, std::unique_ptr<CliOption>> envvarMap;
  static std::function<bool(std::string)> is_keyword;

  static int handle_argument(cdo_option_argument &argument, const std::vector<std::string> &p_argv, int i);

public:
  static void set_keyword_detection(std::function<bool(const std::string &)> f);
  static const std::map<std::string, std::unique_ptr<CliOption>> &
  get_options()
  {
    return optionMap;
  }
  static const std::map<std::string, std::unique_ptr<CliOption>> &
  get_envvars()
  {
    return envvarMap;
  }
  static std::string pad_size_terminal(char padWith, std::string const &sourround = std::string());
  const static int EXIT_REQUESTED;
  const static int ABORT_REQUESTED;
  const static int padding;
  static bool print_settings;
  static bool print_envvars;

  static int parse(std::vector<std::string> p_argv);
  static void get_env_vars();
  static std::unique_ptr<CliOption> &option(std::string const &p_name);
  static std::unique_ptr<CliOption> &envvar(std::string const &p_name);
  static CliOption *option_from_envvar(std::string const &p_envvarName);
  static void
  shortform(std::unique_ptr<CliOption> &&long_form, std::string shortform)
  {
    if (optionMap.find(shortform) != optionMap.end()) { cdo_abort("option shortform already exists: %s", shortform); }
    optionMap[shortform] = std::move(long_form);
  }

  static cdo_option_argument &
  get_argument(const std::map<std::string, std::unique_ptr<CliOption>>::iterator &it)
  {
    return it->second->argument;
  }
  static void print_available_options();
  static std::string print_options_help(std::string const &p_category = DEFAULT_CATEGORY);
  static std::string print_envvar_help();
  static void print_registry(const std::map<std::string, std::unique_ptr<CliOption>> &p_registry);
  static std::string print_option(const std::string &p_optionName);
  static std::string print_option(const std::unique_ptr<CliOption> &iter);
};
#endif /* _CDO_GETOPT_H */
