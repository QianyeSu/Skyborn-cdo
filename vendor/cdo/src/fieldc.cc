/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "arithmetic.h"
#include "cdo_output.h"
#include "field_functions.h"

template <typename T>
static void
fieldc_mul(Varray<T> &v, size_t len, size_t numMissVals, double missval, double rconst)
{
  auto is_EQ = fp_is_equal;
  T missval1 = missval;
  T missval2 = missval;

  if (numMissVals)
  {
    for (size_t i = 0; i < len; ++i) v[i] = MULM(v[i], rconst);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) v[i] *= rconst;
  }
}

void
fieldc_mul(Field &field, double rconst)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { fieldc_mul(v, n, numMissVals, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static size_t
fieldc_div(Varray<T> &v, size_t len, size_t numMissVals, double missval, double rconst)
{
  auto is_EQ = fp_is_equal;
  T missval1 = missval;
  // T missval2 = missval;

  if (numMissVals || is_equal(rconst, 0))
  {
    for (size_t i = 0; i < len; ++i) v[i] = DIVMX(v[i], rconst);

    if (is_equal(rconst, 0)) numMissVals = len;
  }
  else
  {
    for (size_t i = 0; i < len; ++i) v[i] /= rconst;
  }

  return numMissVals;
}

void
fieldc_div(Field &field, double rconst)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { return fieldc_div(v, n, numMissVals, mv, rconst); };
  field.numMissVals = field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_add(Varray<T> &v, size_t len, size_t numMissVals, double missval, double rconst)
{
  auto is_EQ = fp_is_equal;
  T missval1 = missval;
  T missval2 = missval;

  if (numMissVals)
  {
    for (size_t i = 0; i < len; ++i) v[i] = ADDM(v[i], rconst);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) v[i] += rconst;
  }
}

void
fieldc_add(Field &field, double rconst)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { fieldc_add(v, n, numMissVals, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

void
fieldc_sub(Field &field, double rconst)
{
  fieldc_add(field, -rconst);
}

template <typename T>
static void
fieldc_min(Varray<T> &v, size_t len, size_t numMissVals, double missval, double constVal)
{
  auto is_NE = fp_is_not_equal;
  T missval1 = missval;
  T rconst = constVal;

  if (numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (is_NE(v[i], missval1) && v[i] > rconst) v[i] = rconst;
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
      if (v[i] > rconst) v[i] = rconst;
  }
}

void
fieldc_min(Field &field, double rconst)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { fieldc_min(v, n, numMissVals, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_max(Varray<T> &v, size_t len, size_t numMissVals, double missval, double constVal)
{
  auto is_NE = fp_is_not_equal;
  T missval1 = missval;
  T rconst = constVal;

  if (numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (is_NE(v[i], missval1) && v[i] < rconst) v[i] = rconst;
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
      if (v[i] < rconst) v[i] = rconst;
  }
}

void
fieldc_max(Field &field, double rconst)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { fieldc_max(v, n, numMissVals, mv, rconst); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

template <typename T>
static void
fieldc_mod(Varray<T> &v, size_t len, size_t numMissVals, double missval, double divisor)
{
  (void) numMissVals;
  auto is_EQ = fp_is_equal;
  T missval1 = missval;

  for (size_t i = 0; i < len; ++i) { v[i] = is_EQ(v[i], missval1) ? missval1 : std::fmod((double) v[i], divisor); }
}

void
fieldc_mod(Field &field, double divisor)
{
  auto func = [&](auto &v, auto n, auto numMissVals, double mv) { fieldc_mod(v, n, numMissVals, mv, divisor); };
  field_operation(func, field, field.size, field.numMissVals, field.missval);
}

void
fieldc_function(Field &field, double rconst, int function)
{
  switch (function)
  {
    case FieldFunc_Add: fieldc_add(field, rconst); break;
    case FieldFunc_Sub: fieldc_sub(field, rconst); break;
    case FieldFunc_Mul: fieldc_mul(field, rconst); break;
    case FieldFunc_Div: fieldc_div(field, rconst); break;
    case FieldFunc_Min: fieldc_min(field, rconst); break;
    case FieldFunc_Max: fieldc_max(field, rconst); break;
    case FieldFunc_Mod: fieldc_mod(field, rconst); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
  }
}
