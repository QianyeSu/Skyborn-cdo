/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef CDO_EXPR_H
#define CDO_EXPR_H

#include <cstdio>
#include <cstdarg>
#include <cassert>
#include <string>
#include <vector>
#include <variant>

#include "varray.h"

extern int CDO_parser_errorno;

enum CoordIndex
{
  TIMESTEP = 0,
  DATE,
  TIME,
  DELTAT,
  DAY,
  MONTH,
  YEAR,
  SECOND,
  MINUTE,
  HOUR,
  CALENDAR,
  DOY,  // day of year
  DPY,  // days per year
  LEN
};

enum struct NodeEnum
{
  typeUndef,
  typeCon,
  typeVar,
  typeFun,
  typeOpr,
  typeCmd
};

// constants
struct conNodeType
{
  double value;  // value of constant

  conNodeType() : value(0.0) {}  // Default ctor is needed to initialize variant object
  explicit conNodeType(double _value) : value(_value) {}
};

// variables
struct varNodeType
{
  std::string name;  // variable name

  explicit varNodeType(const char *_name) : name(_name) {}
};

// commands
struct cmdNodeType
{
  std::string cmdName;  // command name
  std::string varName;  // variable name

  cmdNodeType(const char *cname, const char *vname) : cmdName(cname), varName(vname) {}
};

// functions
struct funNodeType
{
  std::string name;                   // function name
  int nops;                           // number of operands
  std::vector<struct nodeType *> op;  // operands

  funNodeType(std::string const &fname, int _nops, va_list args) : name(fname), nops(_nops)
  {
    op.resize(nops);
    for (int i = 0; i < nops; i++) op[i] = va_arg(args, nodeType *);
  }
};

// operators
struct oprNodeType
{
  int oper;                           // operator
  int nops;                           // number of operands
  std::vector<struct nodeType *> op;  // operands

  oprNodeType(int _oper, int _nops, va_list args) : oper(_oper), nops(_nops)
  {
    op.resize(nops);
    for (int i = 0; i < nops; i++) op[i] = va_arg(args, nodeType *);
  }
};

enum struct ParamType
{
  UNDEFINED,
  VAR,
  CONST
};

// parameter
struct ParamEntry
{
  ParamType type = ParamType::UNDEFINED;
  bool isValid = false;
  bool select = false;
  bool remove = false;
  bool hasMV = false;
  int coord = 0;
  int gridID = -1;
  int zaxisID = -1;
  int datatype = -1;
  int steptype = -1;
  size_t ngp = 0;
  size_t nlat = 0;
  size_t nlev = 0;
  size_t numMissVals = 0;
  std::string name;
  std::string longname;
  std::string stdname;
  std::string units;
  double *data = nullptr;
  double *weight = nullptr;
  double missval = 0.0;
};

// clang-format off
struct nodeType
{
  ParamEntry param;
  NodeEnum type{ NodeEnum::typeUndef };  // type of node
  bool isTmpObj = false;
  std::variant<conNodeType, varNodeType, cmdNodeType, funNodeType, oprNodeType> v;

  auto &con() const { assert(std::holds_alternative<conNodeType>(v)); return std::get<conNodeType>(v); }
  auto &var()       { assert(std::holds_alternative<varNodeType>(v)); return std::get<varNodeType>(v); }
  auto &var() const { assert(std::holds_alternative<varNodeType>(v)); return std::get<varNodeType>(v); }
  auto &cmd() const { assert(std::holds_alternative<cmdNodeType>(v)); return std::get<cmdNodeType>(v); }
  auto &fun() const { assert(std::holds_alternative<funNodeType>(v)); return std::get<funNodeType>(v); }
  auto &opr() const { assert(std::holds_alternative<oprNodeType>(v)); return std::get<oprNodeType>(v); }
};
// clang-format on

struct CoordType
{
  Varray<double> data;
  std::string units;
  std::string longname;
  size_t size;
  int coord;
  int cdiID;
  bool needed;
};

struct ParseParamType
{
  std::vector<ParamEntry> params;
  std::vector<CoordType> coords;
  std::vector<bool> needed;
  int maxParams;
  int numParams;
  int numVars1;
  int numCoords;
  int maxCoords;
  int tsID;
  int pointID;
  int zonalID;
  int surfaceID;
  bool init;
  bool debug;
};

typedef union
{
  double cvalue;   // constant value
  char *varnm;     // variable name
  char *fname;     // function name
  nodeType *nPtr;  // node pointer
} yysType;

#define YYSTYPE yysType
#define YY_EXTRA_TYPE ParseParamType *

#define YY_DECL int yylex(YYSTYPE *yylval_param, void *yyscanner)
YY_DECL;

int yyparse(ParseParamType &parseArg, void *);
void yyerror(const ParseParamType &parseArg, void *scanner, const char *errstr);

int yylex_init(void **);
int yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);

nodeType *expr_run(nodeType *p, ParseParamType &parseArg);
int params_get_coord_ID(const ParseParamType &parseArg, int coord, int cdiID);

#endif
