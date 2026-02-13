/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Math       abs             Absolute value
      Math       sqr             Square
      Math       sqrt            Square root
      Math       exp             Exponential
      Math       ln              Natural logarithm
      Math       log10           Base 10 logarithm
      Math       sin             Sine
      Math       cos             Cosine
      Math       tan             Tangent
      Math       asin            Arc sine
      Math       acos            Arc cosine
      Math       atan            Arc tangent
      Math       pow             Power
      Math       reci            Reciprocal
*/

#include <cstdlib>
#include <cdi.h>

#include "arithmetic.h"
#include "cdo_options.h"
#include "process_int.h"
#include "param_conversion.h"

template <typename T>
static void
check_out_of_range(Varray<T> &v, size_t &numMissVals, double missval_, double rmin_, double rmax_)
{
  T missval = static_cast<T>(missval_);
  T rmin = static_cast<T>(rmin_);
  T rmax = static_cast<T>(rmax_);

  if (numMissVals)
  {
    for (size_t i = 0, n = v.size(); i < n; ++i)
      if (fp_is_not_equal(v[i], missval) && (v[i] < rmin || v[i] > rmax))
      {
        v[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0, n = v.size(); i < n; ++i)
      if (v[i] < rmin || v[i] > rmax)
      {
        v[i] = missval;
        numMissVals++;
      }
  }
}

static void
check_out_of_range(Field &field, double rmin, double rmax)
{
  auto func = [&](auto &v) { check_out_of_range(v, field.numMissVals, field.missval, rmin, rmax); };
  field_operation(func, field);
  field_num_mv(field);
}

template <typename T>
static void
check_lower_range(Varray<T> &v, size_t &numMissVals, double missval_, double rmin_)
{
  T missval = static_cast<T>(missval_);
  T rmin = static_cast<T>(rmin_);

  if (numMissVals)
  {
    for (size_t i = 0, n = v.size(); i < n; ++i)
      if (fp_is_not_equal(v[i], missval) && v[i] < rmin)
      {
        v[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0, n = v.size(); i < n; ++i)
      if (v[i] < rmin)
      {
        v[i] = missval;
        numMissVals++;
      }
  }
}

static void
check_lower_range(Field &field, double rmin)
{
  auto func = [&](auto &v) { check_lower_range(v, field.numMissVals, field.missval, rmin); };
  field_operation(func, field);
  field_num_mv(field);
}

template <typename T, class UnaryOperation>
void
math_field_transform(Varray<T> &v, size_t numMissVals, double missval_, UnaryOperation unary_op)
{
  T missval = static_cast<T>(missval_);

  if (numMissVals)
    for (size_t i = 0, n = v.size(); i < n; ++i) { v[i] = fp_is_equal(v[i], missval) ? missval : unary_op(v[i]); }
  else
    for (size_t i = 0, n = v.size(); i < n; ++i) { v[i] = unary_op(v[i]); }
}

template <class UnaryOperation>
static void
math_field_transform(Field &field, UnaryOperation unary_op)
{
  auto func = [&](auto &v) { math_field_transform(v, field.numMissVals, field.missval, unary_op); };
  field_operation(func, field);
}

template <typename T>
void
math_field_sqrt(Varray<T> &v, double missval)
{
  T missval1 = static_cast<T>(missval);
  auto is_EQ = fp_is_equal;
  for (size_t i = 0, n = v.size(); i < n; ++i) { v[i] = SQRTM(v[i]); }
}

static void
math_field_sqrt(Field &field)
{
  auto func = [&](auto &v) { math_field_sqrt(v, field.missval); };
  field_operation(func, field);
}

template <typename T>
void
math_field_pow(Varray<T> &v, double missval_, double rc)
{
  T missval = static_cast<T>(missval_);
  for (size_t i = 0, n = v.size(); i < n; ++i) { v[i] = fp_is_equal(v[i], missval) ? missval : std::pow(v[i], rc); }
}

static void
math_field_pow(Field &field, double rc)
{
  auto func = [&](auto &v) { math_field_pow(v, field.missval, rc); };
  field_operation(func, field);
}

template <typename T>
void
math_field_rand(Varray<T> &v, double missval_)
{
  T missval = static_cast<T>(missval_);
  for (size_t i = 0, n = v.size(); i < n; ++i)
  {
    v[i] = fp_is_equal(v[i], missval) ? missval : ((double) std::rand()) / ((double) RAND_MAX);
  }
}

static void
math_field_rand(Field &field)
{
  auto func = [&](auto &v) { math_field_rand(v, field.missval); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_sqr_cplx(Varray<T> &vc, size_t n)
{
  for (size_t i = 0; i < n; ++i)
  {
    vc[i * 2] = vc[i * 2] * vc[i * 2] + vc[i * 2 + 1] * vc[i * 2 + 1];
    vc[i * 2 + 1] = 0.0;
  }
}

static void
math_field_sqr_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_sqr_cplx(v, field.gridsize); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_sqrt_cplx(Varray<T> &vc, size_t n, double missval)
{
  auto is_EQ = fp_is_equal;
  T missval1 = static_cast<T>(missval);
  T missval2 = static_cast<T>(missval);
  auto rsqrt2 = 1.0 / std::sqrt(2.0);
  for (size_t i = 0; i < n; ++i)
  {
    double abs = SQRTM(ADDM(MULM(vc[2 * i], vc[2 * i]), MULM(vc[2 * i + 1], vc[2 * i + 1])));
    vc[i * 2] = MULM(rsqrt2, SQRTM(ADDM(vc[i * 2], abs)));
    vc[i * 2 + 1] = MULM(rsqrt2, DIVM(vc[2 * i + 1], SQRTM(ADDM(vc[2 * i], abs))));
  }
}

static void
math_field_sqrt_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_sqrt_cplx(v, field.gridsize, field.missval); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_conj_cplx(Varray<T> &vc, size_t n)
{
  for (size_t i = 0; i < n; ++i)
  {
    // vc[i * 2] = vc[i * 2];
    vc[i * 2 + 1] = -vc[i * 2 + 1];
  }
}

static void
math_field_conj_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_conj_cplx(v, field.gridsize); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_abs_cplx(Varray<T> &vc, size_t n, double missval)
{
  auto is_EQ = fp_is_equal;
  T missval1 = static_cast<T>(missval);
  T missval2 = static_cast<T>(missval);
  for (size_t i = 0; i < n; ++i) { vc[i] = SQRTM(ADDM(MULM(vc[2 * i], vc[2 * i]), MULM(vc[2 * i + 1], vc[2 * i + 1]))); }
}

static void
math_field_abs_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_abs_cplx(v, field.gridsize, field.missval); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_arg_cplx(Varray<T> &vc, size_t n, double missval_)
{
  T missval = static_cast<T>(missval_);
  for (size_t i = 0; i < n; ++i)
  {
    vc[i]
        = (fp_is_equal(vc[2 * i], missval) || fp_is_equal(vc[2 * i + 1], missval)) ? missval : std::atan2(vc[2 * i + 1], vc[2 * i]);
  }
}

static void
math_field_arg_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_arg_cplx(v, field.gridsize, field.missval); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_re_cplx(Varray<T> &vc, size_t n)
{
  for (size_t i = 0; i < n; ++i) { vc[i] = vc[i * 2]; }
}

static void
math_field_re_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_re_cplx(v, field.gridsize); };
  field_operation(func, field);
}

template <typename T>
static void
math_field_im_cplx(Varray<T> &vc, size_t n)
{
  for (size_t i = 0; i < n; ++i) { vc[i] = vc[i * 2 + 1]; }
}

static void
math_field_im_cplx(Field &field)
{
  auto func = [&](auto &v) { math_field_im_cplx(v, field.gridsize); };
  field_operation(func, field);
}

enum struct Oper
{
  Abs,
  Int,
  Nint,
  Sqr,
  Sqrt,
  Exp,
  Ln,
  Log10,
  Sin,
  Cos,
  Tan,
  Asin,
  Acos,
  Atan,
  Pow,
  Rand,
  Reci,
  Not,
  Conj,
  Re,
  Im,
  Arg
};

class Math : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Math",
    .operators = { { "abs", (int) Oper::Abs, 0, MathHelp },   { "int", (int) Oper::Int, 0, MathHelp },
                   { "nint", (int) Oper::Nint, 0, MathHelp }, { "sqr", (int) Oper::Sqr, 0, MathHelp },
                   { "sqrt", (int) Oper::Sqrt, 0, MathHelp }, { "exp", (int) Oper::Exp, 0, MathHelp },
                   { "ln", (int) Oper::Ln, 0, MathHelp },     { "log10", (int) Oper::Log10, 0, MathHelp },
                   { "sin", (int) Oper::Sin, 0, MathHelp },   { "cos", (int) Oper::Cos, 0, MathHelp },
                   { "tan", (int) Oper::Tan, 0, MathHelp },   { "asin", (int) Oper::Asin, 0, MathHelp },
                   { "acos", (int) Oper::Acos, 0, MathHelp }, { "atan", (int) Oper::Atan, 0, MathHelp },
                   { "pow", (int) Oper::Pow, 0, MathHelp },   { "rand", (int) Oper::Rand, 0, MathHelp },
                   { "reci", (int) Oper::Reci, 0, MathHelp }, { "not", (int) Oper::Not, 0, MathHelp },
                   { "conj", (int) Oper::Conj, 0, MathHelp }, { "re", (int) Oper::Re, 0, MathHelp },
                   { "im", (int) Oper::Im, 0, MathHelp },     { "arg", (int) Oper::Arg, 0, MathHelp } },
    .aliases = { { "log", "ln" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Math> registration = RegisterEntry<Math>();

private:
  Oper operfunc{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  VarList varList1{};

  double rc = 0.0;

public:
  void
  init() override
  {
    auto operatorID = cdo_operator_id();
    operfunc = (Oper) cdo_operator_f1(operatorID);

    if (operfunc == Oper::Pow)
    {
      operator_input_arg("value");
      rc = parameter_to_double(cdo_operator_argv(0));
    }
    else { operator_check_argc(0); }

    if (operfunc == Oper::Rand) std::srand(Options::Random_Seed);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    varList1 = VarList(vlistID1);

    if (operfunc == Oper::Re || operfunc == Oper::Im || operfunc == Oper::Abs || operfunc == Oper::Arg)
    {
      auto numVars = varList1.numVars();
      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var1 = varList1.vars[varID];
        if (var1.dataType == CDI_DATATYPE_CPX32) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
        if (var1.dataType == CDI_DATATYPE_CPX64) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
      }
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        auto number = var.nwpv;
        if (number == CDI_REAL)
        {
          // clang-format off
          switch (operfunc)
          {
            case Oper::Abs:   math_field_transform(field, unary_op_abs); break;
            case Oper::Int:   math_field_transform(field, unary_op_int); break;
            case Oper::Nint:  math_field_transform(field, unary_op_nint); break;
            case Oper::Sqr:   math_field_transform(field, unary_op_sqr); break;
            case Oper::Sqrt:  math_field_sqrt(field); break;
            case Oper::Exp:   math_field_transform(field, unary_op_exp); break;
            case Oper::Ln:    check_lower_range(field, -1);
                              math_field_transform(field, unary_op_log); break;
            case Oper::Log10: check_lower_range(field, -1);
                              math_field_transform(field, unary_op_log10); break;
            case Oper::Sin:   math_field_transform(field, unary_op_sin); break;
            case Oper::Cos:   math_field_transform(field, unary_op_cos); break;
            case Oper::Tan:   math_field_transform(field, unary_op_tan); break;
            case Oper::Asin:  check_out_of_range(field, -1, 1);
                              math_field_transform(field, unary_op_asin); break;
            case Oper::Acos:  check_out_of_range(field, -1, 1);
                              math_field_transform(field, unary_op_acos); break;
            case Oper::Atan:  math_field_transform(field, unary_op_atan); break;
            case Oper::Pow:   math_field_pow(field, rc); break;
            case Oper::Rand:  math_field_rand(field); break;
            case Oper::Reci:  math_field_transform(field, unary_op_reci); break;
            case Oper::Not:   math_field_transform(field, unary_op_not); break;
            case Oper::Re:
            case Oper::Arg:   math_field_transform(field, unary_op_nop); break;
            default: cdo_abort("Operator not implemented for real data!"); break;
          }
          // clang-format on

          field_num_mv(field);
        }
        else
        {
          // clang-format off
          switch (operfunc)
          {
            case Oper::Sqr:   math_field_sqr_cplx(field); break;
            case Oper::Sqrt:  math_field_sqrt_cplx(field); break;
            case Oper::Conj:  math_field_conj_cplx(field); break;
            case Oper::Re:    math_field_re_cplx(field); break;
            case Oper::Im:    math_field_im_cplx(field); break;
            case Oper::Abs:   math_field_abs_cplx(field); break;
            case Oper::Arg:   math_field_arg_cplx(field); break;
            default: cdo_abort("Fields with complex numbers are not supported by this operator!"); break;
          }
          // clang-format on

          field.numMissVals = 0;
        }

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
