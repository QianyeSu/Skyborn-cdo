/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef DATARANGELIST_H
#define DATARANGELIST_H

#include <cdi.h>
#include <cstddef>
#include "cdo_output.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "varray.h"

struct Datarange
{
  size_t gridsize = 0;
  double missval = 0.0;
  double addoffset = 0.0;
  double scalefactor = 1.0;
  int datatype = 0;
  bool checkDatarange = false;

  template <typename T>
  void
  check_datarange(T *array, size_t numMissVals)
  {
    const auto mm = numMissVals ? varray_min_max_mv(gridsize, array, (T) missval) : varray_min_max(gridsize, array);
    const auto numberOfValues = numMissVals ? mm.n : gridsize;

    if (numberOfValues > 0)
      {
        auto smin = (mm.min - addoffset) / scalefactor;
        auto smax = (mm.max - addoffset) / scalefactor;

        if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
            || datatype == CDI_DATATYPE_UINT16)
          {
            smin = (int) std::lround(smin);
            smax = (int) std::lround(smax);
          }

        double vmin = 0.0, vmax = 0.0;
        // clang-format off
        if      (datatype == CDI_DATATYPE_INT8  ) { vmin =        -128.0; vmax =        127.0; }
        else if (datatype == CDI_DATATYPE_UINT8 ) { vmin =           0.0; vmax =        255.0; }
        else if (datatype == CDI_DATATYPE_INT16 ) { vmin =      -32768.0; vmax =      32767.0; }
        else if (datatype == CDI_DATATYPE_UINT16) { vmin =           0.0; vmax =      65535.0; }
        else if (datatype == CDI_DATATYPE_INT32 ) { vmin = -2147483648.0; vmax = 2147483647.0; }
        else if (datatype == CDI_DATATYPE_UINT32) { vmin =           0.0; vmax = 4294967295.0; }
        else if (datatype == CDI_DATATYPE_FLT32 ) { vmin =  -3.40282e+38; vmax =  3.40282e+38; }
        else                                      { vmin =      -1.e+300; vmax =      1.e+300; }
        // clang-format on

        if (smin < vmin || smax > vmax)
          cdo_warning("Some data values (min=%g max=%g) are outside the\n"
                      "    valid range (%g - %g) of the used output precision!\n"
                      "    Use the CDO option%s -b F64 to increase the output precision.",
                      smin, smax, vmin, vmax, (datatype == CDI_DATATYPE_FLT32) ? "" : " -b F32 or");
      }
  }
};

#endif
