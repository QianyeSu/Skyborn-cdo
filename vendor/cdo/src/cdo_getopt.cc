/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_getopt.h"

#include <cstdio>
#include <sstream>
#include <cstdlib>

#include <vector>
#include <map>

#define CLIOP_DBG false

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include <windows.h>
#else
#include <sys/ioctl.h>
#endif
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>

int get_width()
{
  int terminal_width = 120;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  if (GetConsoleScreenBufferInfo(GetStdHandle(STD_ERROR_HANDLE), &csbi))
    terminal_width = csbi.srWindow.Right - csbi.srWindow.Left + 1;
#else
  struct stat statbuf;
  fstat(2, &statbuf);
  if (S_ISCHR(statbuf.st_mode))
  {
#if defined(__APPLE__)
    // macOS: use default terminal width
#elif defined(__linux__) || defined(__unix__) || defined(_POSIX_VERSION)
    struct winsize w;
    ioctl(STDERR_FILENO, TIOCGWINSZ, &w);
    terminal_width = w.ws_col;
#endif
  }
#endif
  return terminal_width;
}

std::string
CLIOptions::pad_size_terminal(char padWith, std::string const &sourround)
{
  int terminal_width = get_width();

  std::string line;
  if (sourround.empty())
  {
    line = std::string(terminal_width, padWith);
  }
  else
  {
    int front_pad = 3;
    int num_spaces = 2;
    int space_left = terminal_width - front_pad - num_spaces - sourround.size();
    if (space_left > terminal_width)
    {
      space_left = 1;
    }
    if (space_left < 0)
    {
      space_left = 0;
    }

    line = std::string(front_pad, padWith);
    line += " " + sourround + " ";
    line += std::string(space_left, padWith);
  }
  return line + "\n";
}

std::map<std::string, std::unique_ptr<cdo_option_2>> CLIOptions::optionMap;

std::map<std::string, std::unique_ptr<cdo_option_2>> CLIOptions::envvarMap;

const int CLIOptions::EXIT_REQUESTED = -2;
const int CLIOptions::ABORT_REQUESTED = -1;
bool CLIOptions::print_settings = false;
bool CLIOptions::print_envvars = false;
const int CLIOptions::padding = 46;

int CLIOptions::parse(std::vector<std::string> p_argv)
{
  int retval = p_argv.size();
  for (size_t i = 1, n = p_argv.size(); i < n; ++i)
  {
    Debug(CLIOP_DBG, "Checking: %s", p_argv[i]);
    const std::map<std::string, std::unique_ptr<cdo_option_2>>::iterator &it = optionMap.find(p_argv[i]);

    if (it == optionMap.end())
    {
      std::string arg = p_argv[i];
      bool isLongFrom = arg.size() >= 2 && arg[0] == '-' && arg[1] == '-';
      bool isShortForm = arg.size() == 2 && arg[0] == '-';
      if (isLongFrom || isShortForm)
      {
        cdo_abort("Option %s not found", p_argv[i]);
      }
      retval = i;
      break;
    }
    if (it->second->hasArgument)
    {
      auto &argument = it->second->argument;
      i++;

      bool out_of_bounds = i >= p_argv.size();
      if (out_of_bounds || optionMap.find(p_argv[i]) != optionMap.end())
      {
        if (argument.default_value.size() > 0)
        {
          argument.value = argument.default_value;
        }
        else if (!argument.is_optional)
        {
          argument.on_empty_argument(p_argv[i - 1]);
          retval = EXIT_REQUESTED;
          break;
        }
      }
      else
      {
        argument.value = p_argv[i];
      }
    }

    Debug(CLIOP_DBG, "executing option %s", (*it).first);
    it->second->execute();

    if (it->second->abortOnExecute)
    {
      retval = EXIT_REQUESTED;
      break;
    }
    if (i >= p_argv.size() - 1)
    {
      // TODO: missing err msg
      retval = ABORT_REQUESTED;
      break;
    }
  }
  if (print_settings)
    print_registry(optionMap);
  if (print_envvars)
    print_registry(envvarMap);

  return retval;
}

const cdo_option_2 *
cdo_option_2::shortform(char p_shortform)
{
  std::string shortform_key = "-";
  shortform_key += p_shortform;
  CLIOptions::shortform(std::make_unique<cdo_option_2>(*this), shortform_key);
  shortName = shortform_key;
  return this;
}

void CLIOptions::get_env_vars()
{
  for (auto &map_entry : envvarMap)
  {
    auto &setting_ptr = map_entry.second;
    const char *envVarValue = getenv(setting_ptr->name.c_str());
    if (envVarValue)
    {
      if (!setting_ptr->hasArgument || *envVarValue)
      {
        Debug(CLIOP_DBG, "Executing envvar %s", setting_ptr->name);
        setting_ptr->argument.value = envVarValue ? std::string(envVarValue) : std::string();
        setting_ptr->execute();
      }
    }
  }
}

void CLIOptions::print_registry(const std::map<std::string, std::unique_ptr<cdo_option_2>> &p_registry)
{
  for (auto const &it : p_registry)
  {
    auto &arg = it.second->argument;
    if (arg.value.size() > 0)
      fprintf(stderr, "%s = %s\n", it.first.c_str(), arg.value.c_str());
  }
}

std::unique_ptr<cdo_option_2> &
CLIOptions::envvar(std::string const &p_name)
{
  if (envvarMap.find(p_name) == envvarMap.end())
  {
    envvarMap[p_name] = std::make_unique<cdo_option_2>();
    envvarMap[p_name]->name = p_name;
  }
  else
  {
    cdo_abort("Environment Variable already registered!");
  }
  return envvarMap[p_name];
}

std::unique_ptr<cdo_option_2> &
CLIOptions::option(std::string const &p_name)
{
  std::string name = "--" + p_name;
  Debug(CLIOP_DBG, "registering key: %s", name);
  if (optionMap.find(name) != optionMap.end())
  {
    cdo_abort("option name already exists: %s", name);
  }

  optionMap[name] = std::make_unique<cdo_option_2>();
  optionMap[name]->name = name;

  return optionMap[name];
}

std::string
CLIOptions::print_envvar_help()
{
  std::stringstream helpString;
  for (auto const &it : envvarMap)
  {
    if (!it.second->isInternal)
    {
      auto len0 = helpString.str().size();
      helpString << std::string(4, ' ');
      helpString << it.second->name;
      if (it.second->hasArgument)
        helpString << " <" + it.second->argument.description + "> ";
      int spaceLeft = padding - helpString.str().size() + len0;
      if (spaceLeft < 0)
        spaceLeft = 0;
      for (auto const &line : it.second->description)
      {
        helpString << std::string(spaceLeft, ' ') + line + "\n";
      }
    }
  }
  return helpString.str();
}

std::string
CLIOptions::print_options_help(std::string const &p_category)
{
  std::stringstream help;
  for (auto &iter : optionMap)
  {
    if (iter.first.size() != 2 && !iter.second->isInternal && iter.second->category == p_category)
    {
      auto &option = iter.second;
      std::string helpString = "    ";
      if (!option->shortName.empty())
      {
        helpString += option->shortName + ", ";
      }
      else
      {
        helpString += "    ";
      }
      helpString += option->name + " ";
      if (option->hasArgument)
        helpString += " <" + option->argument.description + "> ";
      if (option->hasArgument && option->argument.description.empty())
      {
        std::cerr << "error: help argument of " << option->name << " has no description!" << std::endl;
        exit(0);
      }
      int spaceLeft = padding - helpString.size();
      if (spaceLeft <= 0)
      {
        helpString += "\n";
        spaceLeft = padding;
      }
      for (auto const &line : option->description)
      {
        helpString += std::string(spaceLeft, ' ') + line + "\n";
        spaceLeft = padding;
      }
      if (option->description.empty())
        helpString += "\n";
      help << helpString;
    }
  }
  return help.str();
}

#include <algorithm>

cdo_option_2 *
CLIOptions::option_from_envvar(std::string const &p_envvarName)
{
  if (envvarMap.find(p_envvarName) == envvarMap.end())
  {
    cdo_abort("Error envvar %s does not exist", p_envvarName);
  }
  std::string ENVVAR_SUFFIX = "CDO_";
  std::string optionName = p_envvarName.substr(ENVVAR_SUFFIX.size(), p_envvarName.size());

  std::ranges::transform(optionName, optionName.begin(), ::tolower);

  optionName = "--" + optionName;

  if (optionMap.find(optionName) != optionMap.end())
  {
    cdo_abort("Error autogenerated name %s for envvar option %s does already exist!", optionName, p_envvarName);
  }
  auto newOption = std::make_unique<cdo_option_2>();
  auto &envVar = envvarMap[p_envvarName];
  newOption->name = optionName;
  newOption->hasArgument = envVar->hasArgument;
  newOption->description = envVar->description;
  newOption->argument.description = envVar->argument.description;
  newOption->description = {"This option is generated from " + p_envvarName + " and will overwrite it.",
                            "See help of corresponding environment variable."};
  newOption->effect = envVar->effect;
  newOption->isInternal = envVar->isInternal;
  newOption->argument.default_value = envVar->argument.default_value;
  optionMap[optionName] = std::move(newOption);
  return optionMap[optionName].get();
}

void CLIOptions::print_available_options()
{
  for (auto const &iter : optionMap)
  {
    std::cerr << iter.first << std::endl;
  }
  std::cerr << "_---------------------------------_" << std::endl;
  for (auto &iter : optionMap)
  {
    if (iter.first.size() == 2)
      std::cerr << iter.second->shortName[1] << (iter.second->hasArgument ? ":" : "");
  }
}
