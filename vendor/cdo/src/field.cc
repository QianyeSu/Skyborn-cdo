/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include <cassert>
#include <algorithm>

#include "arithmetic.h"
#include "percentiles.h"
#include "varray.h"
#include "field_functions.h"
#include "cdo_output.h"
#include "cdo_omp.h"

void
Field::init(const CdoVar &var)
{
  fpeRaised = 0;
  nwpv = var.nwpv;
  grid = var.gridID;
  gridsize = var.gridsize;
  numMissVals = 0;
  missval = var.missval;
  memType = var.memType;
  size = var.gridsize * var.nwpv;
  m_count = size;
  if (memType == MemType::Float)
    varrayResize(vec_f, size);
  else
    varrayResize(vec_d, size);
}

void
Field::resize(size_t count)
{
  memType = MemType::Double;
  m_count = count;
  varrayResize(vec_d, m_count);
  if (!size) size = m_count;
}

void
Field::resize(size_t count, double value)
{
  memType = MemType::Double;
  m_count = count;
  varrayResizeInit(vec_d, m_count, value);
  if (!size) size = m_count;
}

void
Field::resizef(size_t count)
{
  memType = MemType::Float;
  m_count = count;
  varrayResize(vec_f, m_count);
  if (!size) size = m_count;
}

void
Field::resizef(size_t count, const float value)
{
  memType = MemType::Float;
  m_count = count;
  varrayResizeInit(vec_f, m_count, value);
  if (!size) size = m_count;
}

bool
Field::empty() const
{
  return m_count == 0;
}

void
Field::check_gridsize() const
{
  if (size == 0) fprintf(stderr, "Internal problem, size of field not set!\n");
  if (size > m_count) fprintf(stderr, "Internal problem, size of field is greater than allocated size of field!\n");
}

void
Field3D::init(const CdoVar &var)
{
  nlevels = var.nlevels;
  nwpv = var.nwpv;
  grid = var.gridID;
  gridsize = var.gridsize;
  missval = var.missval;
  memType = var.memType;
  size = var.nlevels * var.gridsize * var.nwpv;
  if (memType == MemType::Float)
    varrayResize(vec_f, size);
  else
    varrayResize(vec_d, size);
}

void
field_fill(Field &field, double value)
{
  field.check_gridsize();

  auto func = [&](auto &v, auto n) { std::fill(v.begin(), v.begin() + n, value); };
  field_operation(func, field, field.size);
}

void
field_ncopy(size_t n, Field const &fieldIn, Field &fieldOut)
{
  if (n > fieldIn.size) cdo_abort("Source field to small (%s)", __func__);
  if (n > fieldOut.size) cdo_abort("Target field to small (%s)", __func__);

  fieldOut.numMissVals = fieldIn.numMissVals;

  auto func = [&](auto const &v1, auto &v2) { varray_copy(n, v1, v2); };
  field_operation2(func, fieldIn, fieldOut);
}

void
field_copy(Field const &fieldIn, Field &fieldOut)
{
  field_ncopy(fieldIn.size, fieldIn, fieldOut);
}

void
field_copy(const Field3D &fieldIn, Field3D &fieldOut)
{
  if (fieldIn.size > fieldOut.size) cdo_abort("Target field to small (%s)", __func__);

  fieldOut.numMissVals = fieldIn.numMissVals;

  auto func = [&](auto &v1, auto &v2) { std::copy(v1.begin(), v1.end(), v2.begin()); };
  field_operation2(func, fieldIn, fieldOut);
}

void
field_copy(const Field3D &fieldIn, int levelID, Field &fieldOut)
{
  auto size = fieldIn.gridsize * fieldIn.nwpv;
  auto offset = levelID * size;
  auto func = [&](auto &v1, auto &v2) { std::copy(v1.begin() + offset, v1.begin() + offset + size, v2.begin()); };
  field_operation2(func, fieldIn, fieldOut);
}

void
field_add(Field &field1, const Field3D &field2, int levelID)
{
  auto size = field1.gridsize * field1.nwpv;
  auto offset = levelID * size;
  auto func = [&](auto &v1, auto const &v2)
  {
    for (size_t i = 0; i < size; ++i) v1[i] += v2[offset + i];
  };
  field_operation2(func, field1, field2);
}

// functor that returns true if value is equal to the value of the constructor parameter provided
class valueDblIsEqual
{
  double _missval;

public:
  explicit valueDblIsEqual(double missval) : _missval(missval) {}
  bool
  operator()(double value) const
  {
    return fp_is_equal(value, _missval);
  }
};

// functor that returns true if value is equal to the value of the constructor parameter provided
class valueIsEqual
{
  double _missval;

public:
  explicit valueIsEqual(double missval) : _missval(missval) {}
  bool
  operator()(double value) const
  {
    return is_equal(value, _missval);
  }
};

size_t
field_num_NANs(Field const &field)
{
  auto func = [&](auto const &v, auto n)
  {
    size_t numNANs = 0;
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : numNANs)
#endif
    for (size_t i = 0; i < n; ++i) { numNANs += std::isnan(v[i]); }
    return numNANs;
  };
  return field_operation(func, field, field.size);
}

size_t
field_num_mv(Field &field)
{
  auto func = [](auto const &v, auto n, auto mv) { return varray_num_mv(n, v, mv); };
  field.numMissVals = field_operation(func, field, field.size, field.missval);
  return field.numMissVals;
}

MinMax
field_min_max(Field const &field)
{
  auto func = [](auto const &v) { return varray_min_max(v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_min_max_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field);
}

double
field_min(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_min(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_min_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_max(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_max(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_max_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_range(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_range(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_range_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_sum(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_sum(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_sum_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_mean(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_mean(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_mean_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_meanw(Field const &field)
{
  auto func = [](auto const &v, auto n, auto mv, auto const &w) { return varray_weighted_mean(n, v, w, mv); };
  auto func_mv = [](auto const &v, auto n, auto mv, auto const &w) { return varray_weighted_mean_mv(n, v, w, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval, field.weightv)
                           : field_operation(func, field, field.size, field.missval, field.weightv);
}

double
field_avg(Field const &field)
{
  auto func = [](auto const &v, auto n) { return varray_mean(n, v); };
  auto func_mv = [](auto const &v, auto n, auto mv) { return varray_avg_mv(n, v, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval) : field_operation(func, field, field.size);
}

double
field_avgw(Field const &field)
{
  auto func = [](auto const &v, auto n, auto mv, auto const &w) { return varray_weighted_mean(n, v, w, mv); };
  auto func_mv = [](auto const &v, auto n, auto mv, auto const &w) { return varray_weighted_avg_mv(n, v, w, mv); };
  return field.numMissVals ? field_operation(func_mv, field, field.size, field.missval, field.weightv)
                           : field_operation(func, field, field.size, field.missval, field.weightv);
}

double
field_var(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_var(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_var1(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_var_1(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_skew(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_skew(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_kurt(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_kurt(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_median(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_median(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_count(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv) { return varray_count(n, v, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

double
field_varw(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv, auto const &w)
  { return varray_weighted_var(n, v, w, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval, field.weightv);
}

double
field_var1w(Field const &field)
{
  auto func = [](auto const &v, auto n, auto numMissVals, auto mv, auto const &w)
  { return varray_weighted_var_1(n, v, w, numMissVals, mv); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval, field.weightv);
}

double
var_to_std(double rvar, double missval)
{
  if (fp_is_equal(rvar, missval) || rvar < 0) return missval;

  return is_not_equal(rvar, 0) ? std::sqrt(rvar) : 0;
}

double
field_std(Field const &field)
{
  return var_to_std(field_var(field), field.missval);
}

double
field_std1(Field const &field)
{
  return var_to_std(field_var1(field), field.missval);
}

double
field_stdw(Field const &field)
{
  return var_to_std(field_varw(field), field.missval);
}

double
field_std1w(Field const &field)
{
  return var_to_std(field_var1w(field), field.missval);
}

void
field_rms(Field const &field, Field const &field2, Field &field3)
{
  size_t rnumMissVals = 0;
  auto grid1 = field.grid;
  //  size_t numMissVals1   = field.numMissVals;
  const auto array1 = field.vec_d.data();
  auto grid2 = field2.grid;
  //  size_t numMissVals2   = field2.numMissVals;
  const auto array2 = field2.vec_d.data();
  auto missval1 = field.missval;
  auto missval2 = field2.missval;
  auto const &w = field.weightv;
  auto rsum = 0.0, rsumw = 0.0;
  auto is_NE = fp_is_not_equal;
  auto is_EQ = fp_is_equal;

  auto len = gridInqSize(grid1);
  if (len != gridInqSize(grid2)) cdo_abort("fields have different size!");

  // if ( numMissVals1 )
  {
    for (size_t i = 0; i < len; ++i)
      if (is_NE(w[i], missval1))
      {
        rsum = ADDM(rsum, MULM(w[i], MULM(SUBM(array2[i], array1[i]), SUBM(array2[i], array1[i]))));
        rsumw = ADDM(rsumw, w[i]);
      }
  }
  /*
else
  {
    for ( i = 0; i < len; i++ )
      {
        rsum  += w[i] * array1[i];
        rsumw += w[i];
      }
  }
  */

  auto ravg = SQRTM(DIVM(rsum, rsumw));

  if (is_EQ(ravg, missval1)) rnumMissVals++;

  field3.vec_d[0] = ravg;
  field3.numMissVals = rnumMissVals;
}

template <typename T>
double
array_pctl(size_t len, Varray<T> &v, size_t numMissVals, double mv, double pn)
{
  T missval = mv;
  double pctl = missval;

  if (len != numMissVals)
  {
    if (numMissVals)
    {
      Varray<T> v2(len);

      size_t j = 0;
      for (size_t i = 0; i < len; ++i)
        if (fp_is_not_equal(v[i], missval)) v2[j++] = v[i];

      if (numMissVals != len - j)
        cdo_warning("Internal problem, inconsistent number of missing values (numMissVals: exprected=%zu found=%zu!)", numMissVals,
                    len - j);

      pctl = percentile(v2.data(), j, pn);
    }
    else { pctl = percentile(v.data(), len, pn); }
  }

  return pctl;
}

double
field_pctl(Field &field, double pn)
{
  auto func = [&](auto &v, auto n, auto numMissVals, auto mv) { return array_pctl(n, v, numMissVals, mv, pn); };
  return field_operation(func, field, field.size, field.numMissVals, field.missval);
}

static int
compare_double(const void *const a, const void *const b)
{
  const auto *const x = static_cast<const double *>(a);
  const auto *const y = static_cast<const double *>(b);
  return (*x < *y) ? -1 : (*x > *y);
}

double
field_rank(Field &field)
{
  auto res = 0.0;
  // Using first value as reference (observation)
  auto val = field.vec_d[0];
  const auto array = &field.vec_d[1];
  auto len = field.size - 1;

  if (field.numMissVals) return field.missval;

  std::qsort(array, len, sizeof(double), compare_double);

  if (val > array[len - 1])
    res = (double) len;
  else
    for (size_t j = 0; j < len; ++j)
      if (array[j] >= val)
      {
        res = (double) j;
        break;
      }

  return res;
}

double
field_function(Field const &field, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Min:    return field_min(field);
    case FieldFunc_Max:    return field_max(field);
    case FieldFunc_Range:  return field_range(field);
    case FieldFunc_Sum:    return field_sum(field);
    case FieldFunc_Mean:   return field_mean(field);
    case FieldFunc_Avg:    return field_avg(field);
    case FieldFunc_Std:    return field_std(field);
    case FieldFunc_Std1:   return field_std1(field);
    case FieldFunc_Var:    return field_var(field);
    case FieldFunc_Var1:   return field_var1(field);
    case FieldFunc_Meanw:  return field_meanw(field);
    case FieldFunc_Avgw:   return field_avgw(field);
    case FieldFunc_Stdw:   return field_stdw(field);
    case FieldFunc_Std1w:  return field_std1w(field);
    case FieldFunc_Varw:   return field_varw(field);
    case FieldFunc_Var1w:  return field_var1w(field);
    case FieldFunc_Skew:   return field_skew(field);
    case FieldFunc_Kurt:   return field_kurt(field);
    case FieldFunc_Median: return field_median(field);
    case FieldFunc_Count:  return field_count(field);
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
  return 0.0;
}
