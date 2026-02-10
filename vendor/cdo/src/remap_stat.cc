/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_output.h"
#include "remap.h"

static void
check_links(size_t numLinks, size_t gridSize, Varray<size_t> const &cellIndices)
{
  cdo_print("  Number of sparse matrix entries: %zu", numLinks);

  std::vector<size_t> count(gridSize, 0);

#ifdef HAVE_OPENMP4
// #pragma omp simd -> wrong result with clang
#endif
  for (size_t i = 0; i < numLinks; ++i) count[cellIndices[i]]++;

  size_t imin = SIZE_MAX;
  size_t imax = 0;
  for (size_t i = 0; i < gridSize; ++i)
  {
    if (count[i] > 0)
    {
      imin = std::min(imin, count[i]);
      imax = std::max(imax, count[i]);
    }
  }

  size_t idiff = (imax - imin) / 10 + 1;
  size_t icount = 0;
  for (size_t i = 0; i < gridSize; ++i)
    if (count[i] > 0) icount++;

  cdo_print("  Number of target cells participating in remap: %zu", icount);

  if (icount)
  {
    cdo_print("  Min no of entries/cell: %zu", imin);
    cdo_print("  Max no of entries/cell: %zu", imax);

    imax = imin + idiff;
    for (size_t n = 0; n < 10; ++n)
    {
      icount = 0;
      for (size_t i = 0; i < gridSize; ++i)
        if (count[i] >= imin && count[i] < imax) icount++;

      if (icount) cdo_print("  Num of cells with entries between %zu - %zu:  %zu", imin, imax - 1, icount);

      imin = imin + idiff;
      imax = imax + idiff;
    }
  }
}

template <typename T1, typename T2>
static void
remap_stat(Varray<T1> const &array1, Varray<T2> const &array2, double missval, int remapOrder, RemapGrid &srcGrid,
           RemapGrid &tgtGrid, RemapVars &rv)
{
  T1 mv1 = missval;
  T2 mv2 = missval;

  cdo_print("%s order mapping from grid1 (%zu) to grid2 (%zu):", (remapOrder == 2) ? "Second" : "First", srcGrid.size,
            tgtGrid.size);
  cdo_print("----------------------------------------------");

  auto mmm = varray_min_max_mean_mv(array1, srcGrid.size, mv1);
  cdo_print("  Grid1 min,mean,max: %g %g %g", mmm.min, mmm.mean, mmm.max);

  mmm = varray_min_max_mean_mv(array2, tgtGrid.size, mv2);
  cdo_print("  Grid2 min,mean,max: %g %g %g", mmm.min, mmm.mean, mmm.max);

  // Conservation Test

  if (srcGrid.cellArea.size())
  {
    double sum = 0.0;
    for (size_t i = 0; i < srcGrid.size; ++i)
      if (fp_is_not_equal(array1[i], mv1)) sum += array1[i] * srcGrid.cellArea[i] * srcGrid.cellFrac[i];
    cdo_print("  Grid1 Integral: %g", sum);

    sum = 0;
    for (size_t i = 0; i < tgtGrid.size; ++i)
      if (fp_is_not_equal(array2[i], mv2)) sum += array2[i] * tgtGrid.cellArea[i] * tgtGrid.cellFrac[i];
    cdo_print("  Grid2 Integral: %g", sum);
  }

  cdo_print("  Number of weights: %zu", rv.numWeights);

  if (rv.numLinks > 0) check_links(rv.numLinks, tgtGrid.size, rv.tgtCellIndices);
}

namespace remap
{

void
stat(int remapOrder, RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv, Field const &field1, Field const &field2)
{
  auto func = [&](auto const &v1, auto const &v2, double mv1) { remap_stat(v1, v2, mv1, remapOrder, srcGrid, tgtGrid, rv); };
  field_operation2(func, field1, field2, field1.missval);
}

};  // namespace remap
