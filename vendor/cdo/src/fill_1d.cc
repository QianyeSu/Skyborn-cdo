/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "fill_1d.h"
#include "interpol.h"

FillMethod
string_to_fillmethod(std::string const &methodStr)
{
  // clang-format off
  if      (methodStr == "nearest")  return FillMethod::Nearest;
  else if (methodStr == "linear")   return FillMethod::Linear;
  else if (methodStr == "forward")  return FillMethod::Forward;
  else if (methodStr == "backward") return FillMethod::Backward;
  // clang-format on

  return FillMethod::Undefined;
}

static void
nearest_neighbour_with_lowest_delta(Varray<double> const &xValues, Varray<double> &yValues, int firstIndex, int lastIndex,
                                    int limit)
{
  auto startIndex = firstIndex + 1;
  auto endIndex = (limit > 0) ? std::min(lastIndex, startIndex + limit) : lastIndex;
  for (int k = startIndex; k < endIndex; ++k)
  {
    // nearest_neighbour
    auto delta1 = xValues[k] - xValues[firstIndex];
    auto delta2 = xValues[lastIndex] - xValues[k];
    yValues[k] = yValues[(delta1 <= delta2) ? firstIndex : lastIndex];
  }
}

constexpr int UndefIndex = -1;

void
fill_1d_nearest(int numValues, Varray<double> const &timeValues, Varray<double> &dataValues, double missval, int limit, int maxGaps)
{
  int numGaps = 0;
  int firstIndex = UndefIndex;
  int lastIndex = UndefIndex;
  for (int i = 0; i < numValues; ++i)
  {
    // if is_missing_val
    if (fp_is_equal(dataValues[i], missval))
    {
      // search for next non missing value
      for (int k = i + 1; k < numValues; ++k)
      {
        if (fp_is_not_equal(dataValues[k], missval))
        {
          lastIndex = k;
          break;
        }
      }

      if (firstIndex == UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = i;
        if (limit > 0) startIndex = std::max(startIndex, lastIndex - limit);
        for (int k = startIndex; k < lastIndex; ++k) dataValues[k] = dataValues[lastIndex];
        i = lastIndex;           // advance iterator
        firstIndex = lastIndex;  // firstIndex is now the first ever found
        lastIndex = UndefIndex;  // reset current is the first ever found
      }
      else if (firstIndex != UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        nearest_neighbour_with_lowest_delta(timeValues, dataValues, firstIndex, lastIndex, limit);
        i = lastIndex;
        firstIndex = lastIndex;
        lastIndex = UndefIndex;
      }
      else if (firstIndex != UndefIndex && lastIndex == UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = i;
        auto endIndex = (limit > 0) ? std::min(numValues, startIndex + limit) : numValues;
        for (int k = startIndex; k < endIndex; ++k) dataValues[k] = dataValues[firstIndex];
        break;
      }
      else if (firstIndex == UndefIndex && lastIndex == UndefIndex) { break; }
    }
    else { firstIndex = i; }
  }
}

static void
fill_linear(Varray<double> const &xValues, Varray<double> &yValues, int firstIndex, int lastIndex, int limit)
{
  auto startIndex = firstIndex + 1;
  auto endIndex = (limit > 0) ? std::min(lastIndex, startIndex + limit) : lastIndex;
  for (int k = startIndex; k < endIndex; ++k)
  {
    yValues[k] = intlin(xValues[k], yValues[firstIndex], xValues[firstIndex], yValues[lastIndex], xValues[lastIndex]);
  }
}

void
fill_1d_linear(int numValues, Varray<double> const &timeValues, Varray<double> &dataValues, double missval, int limit, int maxGaps)
{
  int numGaps = 0;
  int firstIndex = UndefIndex;
  int lastIndex = UndefIndex;
  for (int i = 0; i < numValues; ++i)
  {
    // if is_missing_val
    if (fp_is_equal(dataValues[i], missval))
    {
      // search for next non missing value
      for (int k = i + 1; k < numValues; ++k)
      {
        if (fp_is_not_equal(dataValues[k], missval))
        {
          lastIndex = k;
          break;
        }
      }

      if (firstIndex == UndefIndex && lastIndex != UndefIndex)
      {
        i = lastIndex;           // advance iterator
        firstIndex = lastIndex;  // firstIndex is now the first ever found
        lastIndex = UndefIndex;  // reset current is the first ever found
      }
      else if (firstIndex != UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        fill_linear(timeValues, dataValues, firstIndex, lastIndex, limit);
        i = lastIndex;
        firstIndex = lastIndex;
        lastIndex = UndefIndex;
      }
      else if (lastIndex == UndefIndex) { break; }
    }
    else { firstIndex = i; }
  }
}

void
fill_1d_forward(int numValues, Varray<double> &dataValues, double missval, int limit, int maxGaps)
{
  int numGaps = 0;
  int firstIndex = UndefIndex;
  int lastIndex = UndefIndex;
  for (int i = 0; i < numValues; ++i)
  {
    // if is_missing_val
    if (fp_is_equal(dataValues[i], missval))
    {
      // search for next non missing value
      for (int k = i + 1; k < numValues; ++k)
      {
        if (fp_is_not_equal(dataValues[k], missval))
        {
          lastIndex = k;
          break;
        }
      }

      if (firstIndex == UndefIndex && lastIndex != UndefIndex)
      {
        i = lastIndex;           // advance iterator
        firstIndex = lastIndex;  // firstIndex is now the first ever found
        lastIndex = UndefIndex;  // reset current is the first ever found
        continue;
      }
      else if (firstIndex != UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = firstIndex + 1;
        auto endIndex = (limit > 0) ? std::min(lastIndex, startIndex + limit) : lastIndex;
        for (int k = startIndex; k < endIndex; ++k) dataValues[k] = dataValues[firstIndex];
        i = lastIndex;
        firstIndex = lastIndex;
        lastIndex = UndefIndex;
      }
      else if (firstIndex != UndefIndex && lastIndex == UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = i;
        auto endIndex = (limit > 0) ? std::min(numValues, startIndex + limit) : numValues;
        for (int k = startIndex; k < endIndex; ++k) dataValues[k] = dataValues[firstIndex];
        break;
      }
      else if (firstIndex == UndefIndex && lastIndex == UndefIndex) { break; }
    }
    else { firstIndex = i; }
  }
}

void
fill_1d_backward(int numValues, Varray<double> &dataValues, double missval, int limit, int maxGaps)
{
  int numGaps = 0;
  int firstIndex = UndefIndex;
  int lastIndex = UndefIndex;
  for (int i = 0; i < numValues; ++i)
  {
    // if is_missing_val
    if (fp_is_equal(dataValues[i], missval))
    {
      // search for next non missing value
      for (int k = i + 1; k < numValues; ++k)
      {
        if (fp_is_not_equal(dataValues[k], missval))
        {
          lastIndex = k;
          break;
        }
      }

      if (firstIndex == UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = i;
        if (limit > 0) startIndex = std::max(startIndex, lastIndex - limit);
        for (int k = startIndex; k < lastIndex; ++k) dataValues[k] = dataValues[lastIndex];
        i = lastIndex;           // advance iterator
        firstIndex = lastIndex;  // firstIndex is now the first ever found
        lastIndex = UndefIndex;  // reset current is the first ever found
      }
      else if (firstIndex != UndefIndex && lastIndex != UndefIndex)
      {
        if (maxGaps > 0 && numGaps >= maxGaps) break;
        numGaps++;
        auto startIndex = firstIndex + 1;
        if (limit > 0) startIndex = std::max(startIndex, lastIndex - limit);
        for (int k = startIndex; k < lastIndex; ++k) dataValues[k] = dataValues[lastIndex];
        i = lastIndex;
        firstIndex = lastIndex;
        lastIndex = UndefIndex;
      }
      else if (lastIndex == UndefIndex) { break; }
    }
    else { firstIndex = i; }
  }
}
