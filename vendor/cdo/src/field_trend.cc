#include "field_trend.h"
#include "compare.h"
#include "cdo_omp.h"
#include "arithmetic.h"

template <typename T>
static void
calc_trend_sum(FieldVector3D &work, bool hasMissvals, size_t len, Varray<T> const &varray, double mv, double zj, int varID,
               int levelID)
{
  T missval = mv;
  auto &sumj = work[0][varID][levelID].vec_d;
  auto &sumjj = work[1][varID][levelID].vec_d;
  auto &sumjx = work[2][varID][levelID].vec_d;
  auto &sumx = work[3][varID][levelID].vec_d;
  auto &zn = work[4][varID][levelID].vec_d;

  auto trend_sum = [&](auto i, double value)
  {
    sumj[i] += zj;
    sumjj[i] += zj * zj;
    sumjx[i] += zj * value;
    sumx[i] += value;
    zn[i]++;
  };

  auto trend_sum_mv = [&](auto i, T value, auto is_NE)
  {
    if (is_NE(value, missval)) trend_sum(i, value);
  };

  if (hasMissvals)
  {
    if (std::isnan(missval))
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) { trend_sum_mv(i, varray[i], fp_is_not_equal); }
    else
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < len; ++i) { trend_sum_mv(i, varray[i], is_not_equal); }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i) { trend_sum(i, varray[i]); }
  }
}

void
calc_trend_sum(FieldVector3D &work, Field const &field, double zj, int varID, int levelID)
{
  auto hasMissvals = (field.numMissVals > 0);
  auto func = [&](auto const &v) { calc_trend_sum(work, hasMissvals, field.size, v, field.missval, zj, varID, levelID); };
  field_operation(func, field);
}

template <typename T>
static void
sub_trend(double zj, Varray<T> &v1, Varray<double> const &v2, Varray<double> const &v3, size_t len, double mv)
{
  auto missval1 = mv;
  auto missval2 = mv;

  auto sub_kernel = [&](auto is_EQ)
  {
#ifdef _OPENMP
#pragma omp parallel for if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
    for (size_t i = 0; i < len; ++i)
    {
      auto tmp = (is_EQ(v2[i], mv) || is_EQ(v3[i], mv)) ? mv : (v2[i] + v3[i] * zj);
      v1[i] = SUBM(v1[i], tmp);
    }
  };

  std::isnan(mv) ? sub_kernel(fp_is_equal) : sub_kernel(is_equal);
}

void
sub_trend(double zj, Field &field1, Field const &field2, Field const &field3)
{
  auto func = [&](auto &v1) { sub_trend(zj, v1, field2.vec_d, field3.vec_d, field1.size, field1.missval); };
  field_operation(func, field1);
}

void
calc_trend_param(const FieldVector3D &work, Field &paramA, Field &paramB, int varID, int levelID)
{
  auto gridsize = paramA.size;
  auto missval1 = paramA.missval;
  auto missval2 = paramA.missval;

  auto const &sumj = work[0][varID][levelID].vec_d;
  auto const &sumjj = work[1][varID][levelID].vec_d;
  auto const &sumjx = work[2][varID][levelID].vec_d;
  auto const &sumx = work[3][varID][levelID].vec_d;
  auto const &zn = work[4][varID][levelID].vec_d;

  auto trend_kernel = [&](auto i, auto is_EQ)
  {
    auto temp1 = SUBM(sumjx[i], DIVMX(MULM(sumj[i], sumx[i]), zn[i]));
    auto temp2 = SUBM(sumjj[i], DIVMX(MULM(sumj[i], sumj[i]), zn[i]));
    auto temp3 = DIVM(temp1, temp2);

    paramA.vec_d[i] = SUBM(DIVMX(sumx[i], zn[i]), MULM(DIVMX(sumj[i], zn[i]), temp3));
    paramB.vec_d[i] = temp3;
  };

  if (std::isnan(missval1))
    for (size_t i = 0; i < gridsize; ++i) trend_kernel(i, fp_is_equal);
  else
    for (size_t i = 0; i < gridsize; ++i) trend_kernel(i, is_equal);
}
