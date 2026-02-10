/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "arithmetic.h"
#include "cdo_output.h"
#include "cdo_omp.h"
#include "field_functions.h"

// clang-format off

// arithmetic
const auto arith_add   = [](double a, double b) { return a + b; };
const auto arith_sub   = [](double a, double b) { return a - b; };
const auto arith_mul   = [](double a, double b) { return a * b; };
const auto arith_div   = [](double a, double b) { return a / b; };
const auto arith_min   = [](auto a, auto b) { return (b > a) ? a : b; };
const auto arith_max   = [](auto a, auto b) { return (b < a) ? a : b; };
const auto arith_sumq  = [](double a, double b) { return a + b * b; };
const auto arith_atan2 = [](double a, double b) { return std::atan2(a, b); };

// arithmetic with missing values
const auto arith_add_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = (is_EQ(a, mv_a) || is_EQ(b, mv_b)) ? mv_a : arith_add(a, b); };
const auto arith_sub_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = (is_EQ(a, mv_a) || is_EQ(b, mv_b)) ? mv_a : arith_sub(a, b); };
const auto arith_mul_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = (is_EQ(a, 0) || is_EQ(b, 0)) ? 0 : ((is_EQ(a, mv_a) || is_EQ(b, mv_b)) ? mv_a : arith_mul(a, b)); };
const auto arith_div_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = (is_EQ(a, mv_a) || is_EQ(b, mv_b) || is_EQ(b, 0)) ? mv_a : arith_div(a, b); };
const auto arith_min_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = is_EQ(b, mv_b) ? a : (is_EQ(a, mv_a) ? b : arith_min(a, b)); };
const auto arith_max_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { a = is_EQ(b, mv_b) ? a : (is_EQ(a, mv_a) ? b : arith_max(a, b)); };
const auto arith_sum_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                          { if (!is_EQ(b, mv_b)) a = (is_EQ(a, mv_a) ? b : arith_add(a, b)); };
const auto arith_sumq_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                           { if (!is_EQ(b, mv_b)) a = (is_EQ(a, mv_a) ? arith_mul(b, b) : arith_sumq(a, b)); };
const auto arith_atan2_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                            { a = (is_EQ(a, mv_a) || is_EQ(b, mv_b)) ? mv_a : arith_atan2(a, b); };
const auto arith_vinit_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                            { a = !is_EQ(b, mv_b); (void)mv_a; };
const auto arith_vincr_mv = [](auto &a, auto mv_a, auto b, auto mv_b, auto is_EQ)
                            { if (!is_EQ(b, mv_b)) a = a + 1; (void)mv_a; };
// clang-format on

template <typename T1, typename T2>
static void
varray2_div(Varray<T1> &v1, Varray<T2> const &v2, size_t n, double missval1)
{
  assert(n > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(n <= v1.size());
  assert(n <= v2.size());

#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) { v1[i] = is_equal(v2[i], 0.0) ? missval1 : arith_div(v1[i], v2[i]); }
}

template <typename T1, typename T2, typename FUNC>
static void
varray2_arith(Varray<T1> &v1, Varray<T2> const &v2, size_t n, FUNC func)
{
  assert(n > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(n <= v1.size());
  assert(n <= v2.size());

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) { v1[i] = func(v1[i], v2[i]); }
}

auto varray2_func = [](auto &v1, auto const &v2, size_t n, auto arith_func) { varray2_arith(v1, v2, n, arith_func); };

template <typename T1, typename T2, typename FUNC>
static void
varray2_arith_mv(Varray<T1> &v1, Varray<T2> const &v2, size_t n, double missval1, double missval2, FUNC func)
{
  assert(n > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(n <= v1.size());
  assert(n <= v2.size());

  T1 mv1 = missval1;
  T2 mv2 = missval2;

  if (std::isnan(missval1) || std::isnan(missval2))
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n; ++i) { func(v1[i], mv1, v2[i], mv2, fp_is_equal); }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < n; ++i) { func(v1[i], mv1, v2[i], mv2, is_equal); }
  }
}

auto varray2_func_mv = [](auto &v1, auto const &v2, size_t n, double mv1, double mv2, auto arith_mv_func)
{ varray2_arith_mv(v1, v2, n, mv1, mv2, arith_mv_func); };

// init valid values
void
field2_vinit(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_vinit_mv);
  field1.numMissVals = field2.numMissVals;
}

// increment valid values
void
field2_vincr(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_vincr_mv);
  field1.numMissVals = field2.numMissVals;
}

// init valid values
template <typename T1, typename T2>
static void
field2_vinit(Varray<T1> &v1, Varray<T2> const &v2, size_t len, double mv, int vinit)
{
  auto missval = static_cast<T2>(mv);
  for (size_t i = 0; i < len; ++i) v1[i] = (fp_is_equal(v2[i], missval)) ? 0 : vinit;
}

void
field2_vinit(Field &field1, Field const &field2, int vinit)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);
  /*
  auto func = [](auto &v1, auto &v2, auto n, double mv, int vinit) {
    auto missval = static_cast<T2>(mv);
    for (size_t i = 0; i < n; ++i) v1[i] = (fp_is_equal(v2[i], missval)) ? 0 : vinit;
  };
  */
  auto func = [&](auto &v1, auto const &v2, size_t n, double mv) { field2_vinit(v1, v2, n, mv, vinit); };
  field_operation2(func, field1, field2, field2.size, field2.missval);

  field1.numMissVals = field2.numMissVals;
}

// increment valid values
template <typename T1, typename T2>
static void
field2_vincr(Varray<T1> &v1, Varray<T2> const &v2, size_t len, double mv, int vincr)
{
  auto missval = static_cast<T2>(mv);
  for (size_t i = 0; i < len; ++i)
    if (fp_is_not_equal(v2[i], missval)) v1[i] += vincr;
}

void
field2_vincr(Field &field1, Field const &field2, int vincr)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  auto func = [&](auto &v1, auto const &v2, size_t n, double mv) { field2_vincr(v1, v2, n, mv, vincr); };
  field_operation2(func, field1, field2, field2.size, field2.missval);

  field1.numMissVals = field2.numMissVals;
}

void
field2_add(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_add_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_add); }
}

void
field2_sum(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_sum_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_add); }
}

void
field2_sumw(Field &field1, Field const &field2, double w)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i)
      if (fp_is_not_equal(array2[i], missval2))
      {
        if (fp_is_not_equal(array1[i], missval1))
          array1[i] += w * array2[i];
        else
          array1[i] = w * array2[i];
      }

    field_num_mv(field1);
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i) array1[i] += w * array2[i];
  }
}

/*
 * Compute the occurrence of values in field, if they do not equal refval.
 * This can be used to compute the lengths of multiple periods in a timeseries.
 * Missing field values are handled like refval, i.e. they stop a running period.
 * If there is missing data in the occurence field, missing fields values do not
 * change anything (they do not start a non-period by setting occurrence to zero).
 */
void
field2_sumtr(Field &occur, Field const &field, double refval)
{
  auto omissval = occur.missval;
  auto fmissval = field.missval;
  auto &oarray = occur.vec_d;
  auto const &farray = field.vec_d;

  auto len = occur.size;
  if (len != field.size) cdo_abort("Fields have different size (%s)", __func__);

  if (occur.numMissVals || field.numMissVals)
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i)
      if (fp_is_not_equal(farray[i], fmissval))
      {
        if (fp_is_not_equal(oarray[i], omissval))
          oarray[i] = (fp_is_equal(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
        else
          oarray[i] = (fp_is_equal(farray[i], refval)) ? 0.0 : 1.0;
      }
      else
      {
        if (fp_is_not_equal(oarray[i], omissval)) oarray[i] = 0.0;
      }

    field_num_mv(occur);
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i) oarray[i] = (fp_is_equal(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
  }
}

void
field2_sumq(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_sumq_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_sumq); }
}

void
field2_sumsumq(Field &field1, Field &field2, Field const &field3)
{
  field2_sumq(field2, field3);
  field2_sum(field1, field3);
}

void
field2_sumqw(Field &field1, Field const &field2, double w)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i)
      if (fp_is_not_equal(array2[i], missval2))
      {
        if (fp_is_not_equal(array1[i], missval1))
          array1[i] += w * array2[i] * array2[i];
        else
          array1[i] = w * array2[i] * array2[i];
      }

    field_num_mv(field1);
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i) array1[i] += w * array2[i] * array2[i];
  }
}

void
field2_sub(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_sub_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_sub); }
}

void
field2_mul(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_mul_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_mul); }
}

void
field2_div(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_div_mv);
    field_num_mv(field1);
  }
  else
  {
    auto func = [](auto &v1, auto const &v2, size_t n, double mv) { varray2_div(v1, v2, n, mv); };
    field_operation2(func, field1, field2, field1.size, field1.missval);
    field_num_mv(field1);
  }
}

void
field2_atan2(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_atan2_mv);
  field_num_mv(field1);
}

void
field2_set_miss(Field &field1, Field const &field2)
{
  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals)
  {
    auto missval1 = field1.missval;
    auto &array1 = field1.vec_d;
    auto const &array2 = field2.vec_d;

    for (size_t i = 0; i < len; ++i) array1[i] = fp_is_equal(array1[i], missval1) ? array2[i] : array1[i];

    field_num_mv(field1);
  }
}

void
field2_min(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_min_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_min); }
}

void
field2_max(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    field_operation2(varray2_func_mv, field1, field2, field2.size, field1.missval, field2.missval, arith_max_mv);
    field_num_mv(field1);
  }
  else { field_operation2(varray2_func, field1, field2, field2.size, arith_max); }
}

void
field2_maxmin(Field &field1, Field &field2, Field const &field3)
{
  field2_min(field2, field3);
  field2_max(field1, field3);
}

auto field2_set_index = [](auto &v1, auto &v2, auto v3, auto idx)
{
  v2 = v3;
  v1 = idx;
};

template <typename T>
void
field2_minidx(size_t numMissVals3, size_t len, double mv3, Field &field1, Field &field2, Varray<T> const &v3, int idx)
{
  T missval3 = mv3;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &v1 = field1.vec_d;
  auto &v2 = field2.vec_d;

  if (field2.numMissVals || numMissVals3)
  {
    for (size_t i = 0; i < len; ++i)
    {
      if (fp_is_equal(v3[i], missval3))
      {
        if (fp_is_equal(v2[i], missval2)) v1[i] = missval1;
      }
      else if (v3[i] < v2[i] || fp_is_equal(v2[i], missval2)) { field2_set_index(v1[i], v2[i], v3[i], idx); }
    }

    field_num_mv(field1);
    field_num_mv(field2);
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
    {
      if (v3[i] < v2[i]) field2_set_index(v1[i], v2[i], v3[i], idx);
    }
  }
}

void
field2_minidx(Field &field1, Field &field2, Field const &field3, int idx)
{
  if (field1.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  if (field2.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  auto func = [&](auto const &v3) { field2_minidx(field3.numMissVals, field3.size, field3.missval, field1, field2, v3, idx); };
  field_operation(func, field3);
}

template <typename T>
void
field2_maxidx(size_t numMissVals3, size_t len, double mv3, Field &field1, Field &field2, Varray<T> const &v3, int idx)
{
  T missval3 = mv3;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &v1 = field1.vec_d;
  auto &v2 = field2.vec_d;

  if (field2.numMissVals || numMissVals3)
  {
    for (size_t i = 0; i < len; ++i)
    {
      if (fp_is_equal(v3[i], missval3))
      {
        if (fp_is_equal(v2[i], missval2)) v1[i] = missval1;
      }
      else if (v3[i] > v2[i] || fp_is_equal(v2[i], missval2)) { field2_set_index(v1[i], v2[i], v3[i], idx); }
    }

    field_num_mv(field1);
    field_num_mv(field2);
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
    {
      if (v3[i] > v2[i]) field2_set_index(v1[i], v2[i], v3[i], idx);
    }
  }
}

void
field2_maxidx(Field &field1, Field &field2, Field const &field3, int idx)
{
  if (field1.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  if (field2.size != field3.size) cdo_abort("Fields have different size (%s)", __func__);
  auto func = [&](auto const &v3) { field2_maxidx(field3.numMissVals, field3.size, field3.missval, field1, field2, v3, idx); };
  field_operation(func, field3);
}

static size_t
field_set_numMissVals(size_t len, Varray<double> &v, double missval)
{
  size_t numMissVals = 0;

  if (std::isnan(missval))
  {
    for (size_t i = 0; i < len; ++i)
      if (fp_is_equal(v[i], missval) || v[i] < 0)
      {
        v[i] = missval;
        numMissVals++;
      }
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
      if (is_equal(v[i], missval) || v[i] < 0)
      {
        v[i] = missval;
        numMissVals++;
      }
  }

  return numMissVals;
}

void
field2_var(Field &field1, Field const &field2, Field const &field3, int divisor)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;
  auto const &array3 = field3.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    auto is_EQ = fp_is_equal;
    for (size_t i = 0; i < len; ++i)
    {
      auto temp = DIVM(MULM(array1[i], array1[i]), array3[i]);
      array1[i] = DIVM(SUBM(array2[i], temp), array3[i] - divisor);
      if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
    }
  }
  else
  {
    auto is_EQ = is_equal;
    for (size_t i = 0; i < len; ++i)
    {
      auto temp = DIVX(MUL(array1[i], array1[i]), array3[i]);
      array1[i] = DIVX(SUB(array2[i], temp), array3[i] - divisor);
      if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
    }
  }

  field1.numMissVals = field_set_numMissVals(len, array1, missval1);
}

void
field2_std(Field &field1, Field const &field2, Field const &field3, int divisor)
{
  auto missval1 = field1.missval;
  auto &array1 = field1.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  field2_var(field1, field2, field3, divisor);

  size_t numMissVals = 0;
  for (size_t i = 0; i < len; ++i)
    if (fp_is_equal(array1[i], missval1) || array1[i] < 0)
    {
      array1[i] = missval1;
      numMissVals++;
    }
    else { array1[i] = is_not_equal(array1[i], 0) ? std::sqrt(array1[i]) : 0; }
  field1.numMissVals = numMissVals;
}

void
fieldc_var(Field &field1, Field const &field2, int numSets, int divisor)
{
  auto nsetx = numSets - divisor;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (nsetx == 0)
  {
    for (size_t i = 0; i < len; ++i) array1[i] = missval1;
  }
  else if (field1.numMissVals || field2.numMissVals)
  {
    auto is_EQ = fp_is_equal;
    for (size_t i = 0; i < len; ++i)
    {
      auto temp = MULM(array1[i], array1[i]) / numSets;
      array1[i] = SUBM(array2[i], temp) / nsetx;
      if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
    }
  }
  else
  {
    for (size_t i = 0; i < len; ++i)
    {
      auto temp = MUL(array1[i], array1[i]) / numSets;
      array1[i] = SUB(array2[i], temp) / nsetx;
      if (array1[i] < 0 && array1[i] > -1.e-5) array1[i] = 0;
    }
  }

  field1.numMissVals = field_set_numMissVals(len, array1, missval1);
}

void
fieldc_std(Field &field1, Field const &field2, int numSets, int divisor)
{
  auto missval1 = field1.missval;
  auto &array1 = field1.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  fieldc_var(field1, field2, numSets, divisor);

  size_t numMissVals = 0;
  for (size_t i = 0; i < len; ++i)
    if (fp_is_equal(array1[i], missval1) || array1[i] < 0)
    {
      array1[i] = missval1;
      numMissVals++;
    }
    else { array1[i] = is_not_equal(array1[i], 0) ? std::sqrt(array1[i]) : 0; }

  field1.numMissVals = numMissVals;
}

void
field2_moq(Field &field1, Field const &field2)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field2.numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (fp_is_not_equal(array2[i], missval2))
        array1[i] = array2[i] * array2[i];
      else
        array1[i] = missval1;

    field_num_mv(field1);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) array1[i] = array2[i] * array2[i];
  }
}

void
field2_moqw(Field &field1, Field const &field2, double w)
{
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  auto len = field1.size;
  if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field2.numMissVals)
  {
    for (size_t i = 0; i < len; ++i)
      if (fp_is_not_equal(array2[i], missval2))
        array1[i] = w * array2[i] * array2[i];
      else
        array1[i] = missval1;

    field_num_mv(field1);
  }
  else
  {
    for (size_t i = 0; i < len; ++i) array1[i] = w * array2[i] * array2[i];
  }
}

/**
 * Counts the number of nonumMissValsing values. The result of the operation is computed according to the following rules:
 *
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    miss
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */
template <typename T1, typename T2>
static void
varray2_count_mv(Varray<T1> &v1, Varray<T2> const &v2, size_t len, double mv1, double mv2)
{
  assert(len > 0);
  assert(v1.size() > 0);
  assert(v2.size() > 0);
  assert(len <= v1.size());
  assert(len <= v2.size());

  for (size_t i = 0; i < len; ++i)
    if (fp_is_not_equal(v2[i], mv2))
    {
      if (fp_is_not_equal(v1[i], mv1))
        v1[i] += 1.0;
      else
        v1[i] = 1.0;
    }
}

template <typename T>
static void
varray_count(Varray<T> &v, size_t len)
{
  assert(len > 0);

  for (size_t i = 0; i < len; ++i) v[i] += 1.0;
}

void
field2_count(Field &field1, Field const &field2)
{
  if (field1.size != field2.size) cdo_abort("Fields have different size (%s)", __func__);

  if (field1.numMissVals || field2.numMissVals)
  {
    auto func = [](auto &v1, auto const &v2, size_t n, double mv1, double mv2) { varray2_count_mv(v1, v2, n, mv1, mv2); };
    field_operation2(func, field1, field2, field1.size, field1.missval, field2.missval);
    field_num_mv(field1);
  }
  else
  {
    auto func = [](auto &v, size_t n) { varray_count(v, n); };
    field_operation(func, field1, field1.size);
  }
}

void
field2_function(Field &field1, Field const &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Add:     field2_add(field1, field2);   break;
    case FieldFunc_Min:     field2_min(field1, field2);   break;
    case FieldFunc_Max:     field2_max(field1, field2);   break;
    case FieldFunc_Sum:     field2_sum(field1, field2);   break;
    case FieldFunc_Mean:    field2_sum(field1, field2);   break;
    case FieldFunc_Avg:     field2_add(field1, field2);   break;
    case FieldFunc_Sub:     field2_sub(field1, field2);   break;
    case FieldFunc_Mul:     field2_mul(field1, field2);   break;
    case FieldFunc_Div:     field2_div(field1, field2);   break;
    case FieldFunc_Atan2:   field2_atan2(field1, field2); break;
    case FieldFunc_Setmiss: field2_set_miss(field1, field2); break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}
