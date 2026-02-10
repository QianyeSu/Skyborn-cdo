/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "percentiles.h"
#include "field_functions.h"
#include "cdo_output.h"

using funcType1 = double(size_t, Varray<double> const &);
using funcTypeMV1 = double(size_t, Varray<double> const &, double);
using funcType2 = double(size_t, Varray<double> const &, Varray<double> const &, double);
using funcTypeMV2 = double(size_t, Varray<double> const &, Varray<double> const &, double);
using funcType3 = double(size_t, Varray<double> const &, Varray<double> const &, size_t, double);
using funcType4 = double(size_t, Varray<double> const &, size_t, double);

template <typename T>
static void
varray_copy_meridional(size_t i, size_t nx, size_t ny, Varray<T> const &v1, Varray<double> &v2)
{
  for (size_t j = 0; j < ny; ++j) v2[j] = v1[j * nx + i];
}

static void
meridional_kernel_1(Field const &field1, Field &field2, funcType1 func, funcTypeMV1 funcMV)
{
  size_t rnumMissVals = 0;
  auto numMissVals = field1.numMissVals;
  auto missval = field1.missval;
  auto nx = gridInqXsize(field1.grid);
  auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  for (size_t i = 0; i < nx; ++i)
  {
    auto func_copy_meridional = [&](auto const &v1) { varray_copy_meridional(i, nx, ny, v1, v); };
    field_operation(func_copy_meridional, field1);

    auto result = numMissVals ? funcMV(ny, v, missval) : func(ny, v);
    if (fp_is_equal(result, missval)) rnumMissVals++;
    field2.vec_d[i] = result;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_kernel_2(Field const &field1, Field &field2, funcType2 func, funcTypeMV2 funcMV)
{
  size_t rnumMissVals = 0;
  auto numMissVals = field1.numMissVals;
  auto missval = field1.missval;
  auto nx = gridInqXsize(field1.grid);
  auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny), w(ny);

  for (size_t i = 0; i < nx; ++i)
  {
    varray_copy_meridional(i, nx, ny, field1.weightv, w);
    auto func_copy_meridional = [&](auto const &v1) { varray_copy_meridional(i, nx, ny, v1, v); };
    field_operation(func_copy_meridional, field1);

    auto result = numMissVals ? funcMV(ny, v, w, missval) : func(ny, v, w, missval);
    if (fp_is_equal(result, missval)) rnumMissVals++;
    field2.vec_d[i] = result;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_kernel_3(Field const &field1, Field &field2, funcType3 func)
{
  size_t rnumMissVals = 0;
  auto numMissVals = field1.numMissVals;
  auto missval = field1.missval;
  auto nx = gridInqXsize(field1.grid);
  auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny), w(ny);

  for (size_t i = 0; i < nx; ++i)
  {
    varray_copy_meridional(i, nx, ny, field1.weightv, w);
    auto func_copy_meridional = [&](auto const &v1) { varray_copy_meridional(i, nx, ny, v1, v); };
    field_operation(func_copy_meridional, field1);

    auto result = func(ny, v, w, numMissVals, missval);
    if (fp_is_equal(result, missval)) rnumMissVals++;
    field2.vec_d[i] = result;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_kernel_4(Field const &field1, Field &field2, funcType4 func)
{
  size_t rnumMissVals = 0;
  auto numMissVals = field1.numMissVals;
  auto missval = field1.missval;
  auto nx = gridInqXsize(field1.grid);
  auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  for (size_t i = 0; i < nx; ++i)
  {
    auto func_copy_meridional = [&](auto const &v1) { varray_copy_meridional(i, nx, ny, v1, v); };
    field_operation(func_copy_meridional, field1);

    auto numMissval = ny - varray_count(ny, v, numMissVals, missval);
    auto result = func(ny, v, numMissval, missval);
    if (fp_is_equal(result, missval)) rnumMissVals++;
    field2.vec_d[i] = result;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_min(Field const &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_min, varray_min_mv);
}

static void
meridional_max(Field const &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_max, varray_max_mv);
}

static void
meridional_range(Field const &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_range, varray_range_mv);
}

static void
meridional_sum(Field const &field1, Field &field2)
{
  meridional_kernel_1(field1, field2, varray_sum, varray_sum_mv);
}

static void
meridional_meanw(Field const &field1, Field &field2)
{
  meridional_kernel_2(field1, field2, varray_weighted_mean, varray_weighted_mean_mv);
}

static void
meridional_avgw(Field const &field1, Field &field2)
{
  meridional_kernel_2(field1, field2, varray_weighted_mean, varray_weighted_avg_mv);
}

static void
meridional_varw(Field const &field1, Field &field2)
{
  meridional_kernel_3(field1, field2, varray_weighted_var);
}

static void
meridional_var1w(Field const &field1, Field &field2)
{
  meridional_kernel_3(field1, field2, varray_weighted_var_1);
}

static void
meridional_stdw(Field const &field1, Field &field2)
{
  size_t rnumMissVals = 0;
  auto missval = field1.missval;

  auto nx = gridInqXsize(field1.grid);

  meridional_varw(field1, field2);

  for (size_t i = 0; i < nx; ++i)
  {
    auto rstd = var_to_std(field2.vec_d[i], missval);
    if (fp_is_equal(rstd, missval)) rnumMissVals++;
    field2.vec_d[i] = rstd;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_std1w(Field const &field1, Field &field2)
{
  size_t rnumMissVals = 0;
  auto missval = field1.missval;

  auto nx = gridInqXsize(field1.grid);

  meridional_var1w(field1, field2);

  for (size_t i = 0; i < nx; ++i)
  {
    auto rstd = var_to_std(field2.vec_d[i], missval);
    if (fp_is_equal(rstd, missval)) rnumMissVals++;
    field2.vec_d[i] = rstd;
  }

  field2.numMissVals = rnumMissVals;
}

static void
meridional_skew(Field const &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_skew);
}

static void
meridional_kurt(Field const &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_kurt);
}

static void
meridional_median(Field const &field1, Field &field2)
{
  meridional_kernel_4(field1, field2, varray_median);
}

void
meridional_pctl(Field const &field1, Field &field2, double pn)
{
  size_t rnumMissVals = 0;
  auto missval = field1.missval;

  auto nx = gridInqXsize(field1.grid);
  auto ny = gridInqYsize(field1.grid);

  Varray<double> v(ny);

  if (field1.numMissVals)
  {
    for (size_t i = 0; i < nx; ++i)
    {
      size_t k = 0;
      auto func = [&](auto &v1)
      {
        for (size_t j = 0; j < ny; ++j)
          if (fp_is_not_equal(v1[j * nx + i], missval)) v[k++] = v1[j * nx + i];
      };
      field_operation(func, field1);

      if (k > 0) { field2.vec_d[i] = percentile(v.data(), k, pn); }
      else
      {
        field2.vec_d[i] = missval;
        rnumMissVals++;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < nx; ++i)
    {
      if (ny > 0)
      {
        auto func_copy_meridional = [&](auto const &v1) { varray_copy_meridional(i, nx, ny, v1, v); };
        field_operation(func_copy_meridional, field1);

        field2.vec_d[i] = percentile(v.data(), ny, pn);
      }
      else
      {
        field2.vec_d[i] = missval;
        rnumMissVals++;
      }
    }
  }

  field2.numMissVals = rnumMissVals;
}

void
meridional_function(Field const &field1, Field &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Min:    return meridional_min(field1, field2);
    case FieldFunc_Max:    return meridional_max(field1, field2);
    case FieldFunc_Range:  return meridional_range(field1, field2);
    case FieldFunc_Sum:    return meridional_sum(field1, field2);
    case FieldFunc_Meanw:  return meridional_meanw(field1, field2);
    case FieldFunc_Avgw:   return meridional_avgw(field1, field2);
    case FieldFunc_Stdw:   return meridional_stdw(field1, field2);
    case FieldFunc_Std1w:  return meridional_std1w(field1, field2);
    case FieldFunc_Varw:   return meridional_varw(field1, field2);
    case FieldFunc_Var1w:  return meridional_var1w(field1, field2);
    case FieldFunc_Skew:   return meridional_skew(field1, field2);
    case FieldFunc_Kurt:   return meridional_kurt(field1, field2);
    case FieldFunc_Median: return meridional_median(field1, field2);
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}
