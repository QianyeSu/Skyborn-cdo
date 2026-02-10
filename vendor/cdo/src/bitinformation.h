#ifndef BITINFORMATION_H
#define BITINFORMATION_H

#include <cstddef>

constexpr int NBITS = 32;  // Number of bits in type `float`

struct MutualInformation
{
  double M[NBITS] = {};
};

void signed_exponent(float *A, size_t n);
MutualInformation bitinformation(float *A, size_t n);

#endif
