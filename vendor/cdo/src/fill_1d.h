/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FILL_1D_H
#define FILL_1D_H

#include "varray.h"

enum struct FillMethod
{
  Undefined,
  Nearest,
  Linear,
  Forward,
  Backward
};

FillMethod string_to_fillmethod(std::string const &methodStr);

void fill_1d_nearest(int numValues, Varray<double> const &timeValues, Varray<double> &dataValues, double missval, int limit,
                     int maxGaps);
void fill_1d_linear(int numValues, Varray<double> const &timeValues, Varray<double> &dataValues, double missval, int limit,
                    int maxGaps);
void fill_1d_forward(int numValues, Varray<double> &dataValues, double missval, int limit, int maxGaps);
void fill_1d_backward(int numValues, Varray<double> &dataValues, double missval, int limit, int maxGaps);

#endif
