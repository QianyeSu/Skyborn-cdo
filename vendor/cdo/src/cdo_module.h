/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef CDO_MODULE_H
#define CDO_MODULE_H

#include <string>
#include <vector>

#include "operator_help.h"

// Obase uses input name as base name for files e.g 'test' gets used as test_001 test_002 which are created inside the operator
#define OBASE -1
#define INTERNAL 0
#define EXPOSED 1

struct Alias
{
  Alias(std::string const &_alias, std::string const &_original) : alias(_alias), original(_original) {}
  std::string alias;
  std::string original;
};

enum PositionRestrictions
{
  NoRestriction = 0,
  FilesOnly = 1,
  OnlyFirst = 2
};

struct module_constraints
{
  short streamInCnt;   // Number of input streams
  short streamOutCnt;  // Number of output streams
  PositionRestrictions pos_restriction = NoRestriction;
};

class oper_t
{
  inline static const CdoHelp default_help = {};

public:
  std::string name;
  int f1 = 0;
  int f2 = 0;
  const char *enter = nullptr;
  const CdoHelp &help = default_help;
  oper_t();
  oper_t(const char *_name);
  oper_t(const char *_name, int _f1, int _f2, const char *_enter);
  oper_t(const char *_name, int _f1, int _f2);
  oper_t(const char *_name, int _f1, int _f2, const CdoHelp &p_help);
  oper_t(const char *_name, int _f1, int _f2, const char *_enter, const CdoHelp &p_help);
  oper_t(const char *_name, const CdoHelp &p_help);
  oper_t(const char *_name, int _f1, int _f2, const CdoHelp &&p_help) = delete;
  oper_t(const char *_name, int _f1, int _f2, const char *_enter, const CdoHelp &&p_help) = delete;
  oper_t(const char *_name, const CdoHelp &&p_help) = delete;

  oper_t(int _f1, int _f2, const char *_name, const char *_enter);
};

#include "oper_args.h"

struct CdoModule
{
public:
  std::string name;
  std::vector<oper_t> operators;  // Operator names
  std::vector<Alias> aliases;
  short mode;    // Module mode: 0:intern 1:extern
  short number;  // Allowed number type
  module_constraints constraints;
  Arguments arguments = {};

  std::string toString() const;

  int get_id(std::string const &oper) const;
  int get_stream_in_cnt() const;
  int get_stream_out_cnt() const;
  int get_number() const;
  int get_mode() const;
  int get_pos_restriction() const;

  // returns the alias id or if not found1
  int is_alias(std::string const &subject) const;

  std::vector<oper_t>::const_iterator find_operator(std::string const &p_operatorName) const;
  const CdoHelp &get_help(std::string const &p_operatorName) const;
};

#endif
