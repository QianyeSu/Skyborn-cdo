/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include "cpp_lib.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_timer.h"
#include "cdo_omp.h"
#include <mpim_grid.h>
#include "remap.h"
#include "remap_store_link.h"
#include "progress.h"

// bicubic interpolation

static void
bicubic_set_weights(double xfrac, double yfrac, double (&weights)[4][4])
{
  auto xfrac1 = xfrac * xfrac * (xfrac - 1.0);
  auto xfrac2 = xfrac * (xfrac - 1.0) * (xfrac - 1.0);
  auto xfrac3 = xfrac * xfrac * (3.0 - 2.0 * xfrac);
  auto yfrac1 = yfrac * yfrac * (yfrac - 1.0);
  auto yfrac2 = yfrac * (yfrac - 1.0) * (yfrac - 1.0);
  auto yfrac3 = yfrac * yfrac * (3.0 - 2.0 * yfrac);
  // clang-format off
  weights[0][0] = (1.0 - yfrac3) * (1.0 - xfrac3);
  weights[1][0] = (1.0 - yfrac3) *        xfrac3;
  weights[2][0] =        yfrac3  *        xfrac3;
  weights[3][0] =        yfrac3  * (1.0 - xfrac3);
  weights[0][1] = (1.0 - yfrac3) *        xfrac2;
  weights[1][1] = (1.0 - yfrac3) *        xfrac1;
  weights[2][1] =        yfrac3  *        xfrac1;
  weights[3][1] =        yfrac3  *        xfrac2;
  weights[0][2] =        yfrac2  * (1.0 - xfrac3);
  weights[1][2] =        yfrac2  *        xfrac3;
  weights[2][2] =        yfrac1  *        xfrac3;
  weights[3][2] =        yfrac1  * (1.0 - xfrac3);
  weights[0][3] =        yfrac2  *        xfrac2;
  weights[1][3] =        yfrac2  *        xfrac1;
  weights[2][3] =        yfrac1  *        xfrac1;
  weights[3][3] =        yfrac1  *        xfrac2;
  // clang-format on
}

int num_src_points(Vmask const &mask, const size_t (&indices)[4], double (&lats)[4]);

static void
renormalize_weights(const double (&lats)[4], double (&weights)[4][4])
{
  // sum of weights for normalization
  auto sumWeights = std::fabs(lats[0]) + std::fabs(lats[1]) + std::fabs(lats[2]) + std::fabs(lats[3]);
  for (int i = 0; i < 4; ++i) weights[i][0] = std::fabs(lats[i]) / sumWeights;
  for (int i = 0; i < 4; ++i) weights[i][1] = 0.0;
  for (int i = 0; i < 4; ++i) weights[i][2] = 0.0;
  for (int i = 0; i < 4; ++i) weights[i][3] = 0.0;
}
/*
#ifdef HAVE_LIB_RANGES_ZIP
#include <ranges>
static void
bicubic_sort_weights_by_index_zip(size_t (&indices)[4], double (&weights)[4][4])
{
  auto r = std::views::zip(indices, weights);
  std::sort(r.begin(), r.end, [](auto a, auto b) {
    auto [ai, aw] = a;
    auto [bi, bw] = b;
    return (ai < bi);
  });
}
#endif
*/
static void
bicubic_sort_weights_by_index(size_t (&indices)[4], double (&weights)[4][4])
{
  constexpr size_t numWeights = 4;

  struct IndexWeightX
  {
    size_t index;
    double weight[4];
  };

  std::array<IndexWeightX, numWeights> indexWeights;

  for (size_t i = 0; i < numWeights; ++i)
  {
    indexWeights[i].index = indices[i];
    for (int k = 0; k < 4; ++k) indexWeights[i].weight[k] = weights[i][k];
  }

  std::ranges::sort(indexWeights, {}, &IndexWeightX::index);

  for (size_t i = 0; i < numWeights; ++i)
  {
    indices[i] = indexWeights[i].index;
    for (int k = 0; k < 4; ++k) weights[i][k] = indexWeights[i].weight[k];
  }
}

static void
bicubic_sort_weights(size_t (&indices)[4], double (&weights)[4][4])
{
  constexpr size_t numWeights = 4;
  if (is_sorted_list(numWeights, indices)) return;

  // #ifdef HAVE_LIB_RANGES_ZIP
  //   bicubic_sort_weights_by_index_zip(indices, weights);
  // #else
  bicubic_sort_weights_by_index(indices, weights);
  // #endif
}

static void
bicubic_warning()
{
  static auto printWarning = true;
  if (Options::cdoVerbose || printWarning)
  {
    printWarning = false;
    cdo_warning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
  }
}

// -----------------------------------------------------------------------
// This routine computes the weights for a bicubic interpolation.
// -----------------------------------------------------------------------
void
remap_bicubic_weights(RemapSearch &rsearch, RemapVars &rv)
{
  auto const *srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (srcGrid->rank != 2) cdo_abort("Can't do bicubic interpolation if the source grid is not a regular 2D grid!");

  cdo::timer timer;
  cdo::Progress progress;

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks4> weightLinks(tgtGridSize);
  weight_links_4_alloc(tgtGridSize, weightLinks);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    weightLinks[tgtCellIndex].nlinks = 0;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    auto pointLL = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

    SquareCorners squareCorners;
    double weights[4][4];  //  bicubic weights for four corners

    // Find nearest square of grid points on source grid
    auto searchResult = remap_search_square(rsearch, pointLL, squareCorners);

    // Check to see if points are mask points
    if (searchResult > 0) searchResult = remap_check_mask_indices(squareCorners.indices, srcGrid->mask);

    // If point found, find local xfrac/yfrac coordinates for weights
    if (searchResult > 0)
    {
      tgtGrid->cellFrac[tgtCellIndex] = 1.0;

      auto [xfrac, yfrac] = remap_find_weights(pointLL, squareCorners);
      if (xfrac >= 0.0 && yfrac >= 0.0)
      {
        // Successfully found xfrac, yfrac - compute weights
        bicubic_set_weights(xfrac, yfrac, weights);
        store_weightlinks_bicubic(squareCorners.indices, weights, tgtCellIndex, weightLinks);
      }
      else
      {
        bicubic_warning();
        searchResult = -1;
      }
    }

    // Search for bicubic failed - use a distance-weighted average instead
    // (this is typically near the pole) Distance was stored in srcLats!
    if (searchResult < 0)
    {
      if (num_src_points(srcGrid->mask, squareCorners.indices, squareCorners.lats) > 0)
      {
        tgtGrid->cellFrac[tgtCellIndex] = 1.0;
        renormalize_weights(squareCorners.lats, weights);
        store_weightlinks_bicubic(squareCorners.indices, weights, tgtCellIndex, weightLinks);
      }
    }
  }

  weight_links_4_to_remap_links(tgtGridSize, weightLinks, rv);

  rv.numLinksPerValue = 4;

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, timer.elapsed());
}  // remap_bicubic_weights

// -----------------------------------------------------------------------
// This routine computes and apply the weights for a bicubic interpolation.
// -----------------------------------------------------------------------
template <typename T>
static T
bicubic_remap(Varray<T> const &srcArray, const double (&wgt)[4][4], const size_t (&ind)[4], const RemapGradients &gradients)
{
  auto const &glat = gradients.lat;
  auto const &glon = gradients.lon;
  auto const &glatlon = gradients.latLon;

  double tgtPoint = 0.0;
  for (int i = 0; i < 4; ++i)
    tgtPoint += srcArray[ind[i]] * wgt[i][0] + glat[ind[i]] * wgt[i][1] + glon[ind[i]] * wgt[i][2] + glatlon[ind[i]] * wgt[i][3];

  return tgtPoint;
}

template <typename T1, typename T2>
static void
remap_bicubic(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double srcMissval, size_t numMissVals, RemapSearch &rsearch)
{
  T1 missval = srcMissval;
  auto srcGrid = rsearch.srcGrid;
  const auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (srcGrid->rank != 2) cdo_abort("Can't do bicubic interpolation if the source grid is not a regular 2D grid!");

  cdo::timer timer;
  cdo::Progress progress;

  auto tgtGridSize = tgtGrid->size;
  auto srcGridSize = srcGrid->size;

  Vmask srcGridMask(srcGridSize, 1);
  if (numMissVals) remap_set_mask(srcArray, srcGridSize, numMissVals, srcMissval, srcGridMask);

  // Compute mappings from source to target grid

  RemapGradients gradients(srcGrid->size);
  remap::gradients(srcArray, *srcGrid, srcGridMask, gradients);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    auto &tgtValue = tgtArray[tgtCellIndex];
    tgtValue = missval;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    auto pointLL = remapgrid_get_lonlat(tgtGrid, tgtCellIndex);

    SquareCorners squareCorners;
    double weights[4][4];  //  bicubic weights for four corners

    // Find nearest square of grid points on source grid
    auto searchResult = remap_search_square(rsearch, pointLL, squareCorners);

    // Check to see if points are mask points
    if (searchResult > 0) searchResult = remap_check_mask_indices(squareCorners.indices, srcGridMask);

    // If point found, find local xfrac/yfrac coordinates for weights
    if (searchResult > 0)
    {
      auto [xfrac, yfrac] = remap_find_weights(pointLL, squareCorners);
      if (xfrac >= 0.0 && yfrac >= 0.0)
      {
        // Successfully found xfrac, yfrac - compute weights
        bicubic_set_weights(xfrac, yfrac, weights);
        bicubic_sort_weights(squareCorners.indices, weights);
        tgtValue = bicubic_remap(srcArray, weights, squareCorners.indices, gradients);
      }
      else
      {
        bicubic_warning();
        searchResult = -1;
      }
    }

    // Search for bicubic failed - use a distance-weighted average instead
    // (this is typically near the pole) Distance was stored in srcLats!
    if (searchResult < 0)
    {
      if (srcGridMask.size() == 0 || num_src_points(srcGridMask, squareCorners.indices, squareCorners.lats) > 0)
      {
        renormalize_weights(squareCorners.lats, weights);
        bicubic_sort_weights(squareCorners.indices, weights);
        tgtValue = bicubic_remap(srcArray, weights, squareCorners.indices, gradients);
      }
    }
  }

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, timer.elapsed());
}  // remap_bicubic

void
remap_bicubic(RemapSearch &remapSearch, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { remap_bicubic(v1, v2, field1.missval, field1.numMissVals, remapSearch); };
  field_operation2(func, field1, field2);
}
