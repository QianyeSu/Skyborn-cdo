/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_output.h"
#include "cdo_omp.h"
#include "remap.h"

namespace remap
{

template <typename T>
void
gradients(Varray<T> const &array, RemapGrid &grid, Vmask const &mask, RemapGradients &gradients)
{
  if (grid.rank != 2) cdo_abort("Internal problem (remap_gradients), grid rank = %d!", grid.rank);

  auto gridSize = grid.size;
  auto nx = grid.dims[0];
  auto ny = grid.dims[1];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t n = 0; n < gridSize; ++n)
  {
    if (mask[n])
    {
      // clang-format off
      auto delew = 0.5;
      auto delns = 0.5;

      auto j = n / nx + 1;
      auto i = n - (j - 1) * nx + 1;

      auto ip1 = i + 1;
      auto im1 = i - 1;
      auto jp1 = j + 1;
      auto jm1 = j - 1;

      if (ip1 > nx) ip1 = ip1 - nx;
      if (im1 < 1)  im1 = nx;
      if (jp1 > ny) { jp1 = j; delns = 1.0; }
      if (jm1 < 1)  { jm1 = j; delns = 1.0; }

      auto in = (jp1 - 1) * nx + i - 1;
      auto is = (jm1 - 1) * nx + i - 1;
      auto ie = (j - 1) * nx + ip1 - 1;
      auto iw = (j - 1) * nx + im1 - 1;

      auto ine = (jp1 - 1) * nx + ip1 - 1;
      auto inw = (jp1 - 1) * nx + im1 - 1;
      auto ise = (jm1 - 1) * nx + ip1 - 1;
      auto isw = (jm1 - 1) * nx + im1 - 1;

      // Compute i-gradient
      if (!mask[ie]) { ie = n; delew = 1.0; }
      if (!mask[iw]) { iw = n; delew = 1.0; }

      gradients.lat[n] = delew * (array[ie] - array[iw]);

      // Compute j-gradient
      if (!mask[in]) { in = n; delns = 1.0; }
      if (!mask[is]) { is = n; delns = 1.0; }

      gradients.lon[n] = delns * (array[in] - array[is]);
      // clang-format on

      // Compute ij-gradient
      delew = 0.5;
      delns = (jp1 == j || jm1 == j) ? 1.0 : 0.5;

      if (!mask[ine])
      {
        if (in != n)
        {
          ine = in;
          delew = 1.0;
        }
        else if (ie != n)
        {
          ine = ie;
          inw = iw;
          if (inw == n) delew = 1.0;
          delns = 1.0;
        }
        else
        {
          ine = n;
          inw = iw;
          delew = 1.0;
          delns = 1.0;
        }
      }

      if (!mask[inw])
      {
        if (in != n)
        {
          inw = in;
          delew = 1.0;
        }
        else if (iw != n)
        {
          inw = iw;
          ine = ie;
          if (ie == n) delew = 1.0;
          delns = 1.0;
        }
        else
        {
          inw = n;
          ine = ie;
          delew = 1.0;
          delns = 1.0;
        }
      }

      auto gradLatZero = delew * (array[ine] - array[inw]);

      if (!mask[ise])
      {
        if (is != n)
        {
          ise = is;
          delew = 1.0;
        }
        else if (ie != n)
        {
          ise = ie;
          isw = iw;
          if (isw == n) delew = 1.0;
          delns = 1.0;
        }
        else
        {
          ise = n;
          isw = iw;
          delew = 1.0;
          delns = 1.0;
        }
      }

      if (!mask[isw])
      {
        if (is != n)
        {
          isw = is;
          delew = 1.0;
        }
        else if (iw != n)
        {
          isw = iw;
          ise = ie;
          if (ie == n) delew = 1.0;
          delns = 1.0;
        }
        else
        {
          isw = n;
          ise = ie;
          delew = 1.0;
          delns = 1.0;
        }
      }

      auto gradLonZero = delew * (array[ise] - array[isw]);
      gradients.latLon[n] = delns * (gradLatZero - gradLonZero);
    }
    else
    {
      gradients.lat[n] = 0.0;
      gradients.lon[n] = 0.0;
      gradients.latLon[n] = 0.0;
    }
  }
}  // remap_gradients

// Explicit instantiation
template void gradients(Varray<float> const &array, RemapGrid &grid, Vmask const &mask, RemapGradients &gradients);
template void gradients(Varray<double> const &array, RemapGrid &grid, Vmask const &mask, RemapGradients &gradients);

void
gradients(Field const &field, RemapGrid &grid, RemapGradients &gradients)
{
  auto len = grid.size;
  Vmask mask(len);
#ifdef _OPENMP
#pragma omp parallel for simd if (len > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < len; ++i) mask[i] = (grid.mask[i] > 0);

  auto func = [&](auto const &v) { remap::gradients(v, grid, mask, gradients); };
  field_operation(func, field);
}

};  // namespace remap
