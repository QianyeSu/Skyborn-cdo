/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "arithmetic.h"
#include "cdo_output.h"
#include "field_functions.h"

template <typename T>
static void
fieldc_mul_complex(size_t len, size_t numMissVals, Varray<T> &v, double missval, const double rconst[2])
{
  (void) numMissVals;
  // z1 x z2 = (x1x2 - y1y2) + i(x1y2 + x2y1)
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
  {
    double v1r = v[2 * i];
    double v1i = v[2 * i + 1];
    v[2 * i] = SUBM(MULM(v1r, rconst[0]), MULM(v1i, rconst[1]));
    v[2 * i + 1] = ADDM(MULM(v1r, rconst[1]), MULM(v1i, rconst[0]));
  }
}

static void
fieldc_mul_complex(Field &field, const double rconst[2])
{
  auto func = [&](auto &v, auto size, auto numMissVals, double mv) { fieldc_mul_complex(size, numMissVals, v, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_div_complex(size_t len, size_t numMissVals, Varray<T> &v, double missval, const double rconst[2])
{
  (void) numMissVals;
  // z1 / z2 = (x1x2 + y1y2) / (x2x2 + y2y2) + i (y1x2 - x1y2) / (x2x2 + y2y2)
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
  {
    double v1r = v[2 * i];
    double v1i = v[2 * i + 1];
    auto denominator = ADDM(MULM(rconst[0], rconst[0]), MULM(rconst[1], rconst[1]));
    v[2 * i] = DIVM(ADDM(MULM(v1r, rconst[0]), MULM(v1i, rconst[1])), denominator);
    v[2 * i + 1] = DIVM(SUBM(MULM(v1i, rconst[0]), MULM(v1r, rconst[1])), denominator);
  }
}

static void
fieldc_div_complex(Field &field, const double rconst[2])
{
  auto func = [&](auto &v, auto size, auto numMissVals, double mv) { fieldc_div_complex(size, numMissVals, v, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_add_complex(size_t len, size_t numMissVals, Varray<T> &v, double missval, const double rconst[2])
{
  (void) numMissVals;
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
  {
    v[2 * i] = ADDM(v[2 * i], rconst[0]);
    v[2 * i + 1] = ADDM(v[2 * i + 1], rconst[1]);
  }
}

static void
fieldc_add_complex(Field &field, const double rconst[2])
{
  auto func = [&](auto &v, auto size, auto numMissVals, double mv) { fieldc_add_complex(size, numMissVals, v, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_sub_complex(size_t len, size_t numMissVals, Varray<T> &v, double missval, const double rconst[2])
{
  (void) numMissVals;
  auto is_EQ = fp_is_equal;
  auto missval1 = missval;
  auto missval2 = missval;
  for (size_t i = 0; i < len; ++i)
  {
    v[2 * i] = SUBM(v[2 * i], rconst[0]);
    v[2 * i + 1] = SUBM(v[2 * i + 1], rconst[1]);
  }
}

static void
fieldc_sub_complex(Field &field, const double rconst[2])
{
  auto func = [&](auto &v, auto size, auto numMissVals, double mv) { fieldc_sub_complex(size, numMissVals, v, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

void
fieldc_function_complex(Field &field, const double rconst[2], int function)
{
  switch (function)
  {
    case FieldFunc_Add: fieldc_add_complex(field, rconst); break;
    case FieldFunc_Sub: fieldc_sub_complex(field, rconst); break;
    case FieldFunc_Mul: fieldc_mul_complex(field, rconst); break;
    case FieldFunc_Div: fieldc_div_complex(field, rconst); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
  }
}
