#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <cstdint>
#include <cstdbool>
#include <cmath>
#include <vector>

#ifdef CDO
// BitInformation v0.5.1 (https://github.com/milankl/BitInformation.jl)
// Converted from Julia to C++ by Uwe Schulzweida

#include "cdo_options.h"
#include "cdo_omp.h"
#endif

#include "statistic.h"  // cdo::normal_inv()
#include "bitinformation.h"

struct BitpairCounters
{
  int C[NBITS][2][2] = {};
};

/*
  p₁ = binom_confidence(n::Int,c::Real)

  Returns the probability `p₁` of successes in the binomial distribution (p=1/2) of
  `n` trials with confidence `c`.

  # Example
  At c=0.95, i.e. 95% confidence, n=1000 tosses of a coin will yield not more than
  ```julia
  julia> p₁ = BitInformation.binom_confidence(1000,0.95)
  0.5309897516152281
  ```
  about 53.1% heads (or tails).
*/
static double
binom_confidence(size_t n, double c)
{
  double v = 1.0 - (1.0 - c) * 0.5;
  double p = 0.5 + cdo::normal_inv(v) / (2.0 * std::sqrt(n));
  // printf("p = %g %g %g \n", p, cdo::normal_inv(v), v);
  return std::min(1.0, p);  // cap probability at 1 (only important for n small)
}

/*
  entropy2(p1, p2)

  Compute the entropy of a collection of probabilities `p1` and `p2`,
  The entropy is scaled by `1/log(2.0)`.
  Elements with probability 0 or 1 add 0 to the entropy.
*/

static inline double
xlogx(double x)
{
  return x * std::log(x);
}

static inline double
entropy2_kernel(double p1, double p2)
{
  return -(xlogx(p1) + xlogx(p2));
}

static double
entropy2(double p1, double p2)
{
  constexpr double base = 2.0;
  return entropy2_kernel(p1, p2) / std::log(base);
}

/*
  Hf = binom_free_entropy(n::Int,c::Real)

  Returns the free entropy `Hf` associated with `binom_confidence`.
*/
static double
binom_free_entropy(size_t n, double c)
{
  double p = binom_confidence(n, c);
  return 1.0 - entropy2(p, 1.0 - p);
}

/*
  set_zero_insignificant!(H::Vector,nelements::Int,confidence::Real)

  Remove binary information in the vector `H` of entropies that is insignificantly
  different from a random 50/50 by setting it to zero.
*/
static void
set_zero_insignificant(double *H, size_t nelements, double confidence)
{
  double Hfree = binom_free_entropy(nelements, confidence);             // free entropy of random 50/50 at trial size
                                                                        // get chance p for 1 (or 0) from binom distr
  for (int i = 0; i < NBITS; ++i) H[i] = (H[i] <= Hfree) ? 0.0 : H[i];  // set zero what's insignificant
}

/*
  bitpair_count(a::T,b::T,C::Array{Int,3}) where {T<:Integer}

  Update counter array C of size nbits x 2 x 2 for every 00|01|10|11-bitpairing in a,b.
  `nbits` is the number of bits in type `T`.
*/
static void
bitpair_count_kernel(uint32_t a, uint32_t b, BitpairCounters &BC)
{
  uint32_t mask = 1;                    // start with least significant bit
  for (uint32_t i = 0; i < NBITS; ++i)  // loop from least to most significant bit
  {
    uint32_t j = ((a & mask) >> i);  // isolate that bit in a,b
    uint32_t k = ((b & mask) >> i);  // and move to 0x0 or 0x1s
    BC.C[NBITS - i - 1][j][k] += 1;  // to be used as index j,k to increase counter C
    mask <<= 1;                      // shift mask to get the next significant bit
  }
}

/*
  C = bitpair_count(A::AbstractArray{T},B::AbstractArray{T}) where {T<:Union{Integer,AbstractFloat}}

  Returns counter array C of size NBITS x 2 x 2 for every 00|01|10|11-bitpairing in elements of A,B.
*/
static BitpairCounters
bitpair_count(float *A, float *B, size_t n)
{
  // @assert size(A) == size(B) "Size of A=$(size(A)) does not match size of B=$(size(B))"

  // reinterpret arrays A,B as UInt (no mem allocation required)
  uint32_t *Auint = (uint32_t *) A;
  uint32_t *Buint = (uint32_t *) B;

#ifdef CDO
  std::vector<BitpairCounters> BC(Threading::ompNumMaxThreads);
#else
  std::vector<BitpairCounters> BC(1);
#endif

  // loop over all elements in A,B pairwise, inner loop (within bitpair_count): bit positions
  // note this is faster than swapping inner & outer loop
#ifdef CDO
#ifdef _OPENMP
#pragma omp parallel for if (n > 999999) default(shared) schedule(static)
#endif
#endif
  for (size_t i = 0; i < n; ++i)
  {
#ifdef CDO
    auto ompthID = cdo_omp_get_thread_num();
#else
    auto ompthID = 0;
#endif
    bitpair_count_kernel(Auint[i], Buint[i], BC[ompthID]);  // count the bits and update counter array C
  }

#ifdef CDO
#ifdef _OPENMP
  for (int m = 1; m < Threading::ompNumMaxThreads; ++m)
  {
    for (int k = 0; k < NBITS; ++k)
      for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) BC[0].C[k][j][i] += BC[m].C[k][j][i];
  }
#endif
#endif

  return BC[0];
}

/*
  Mutual information from the joint probability mass function p
  of two variables X,Y. p is an nx x ny array which elements represent
  the probabilities of observing x~X and y~Y.
*/
static double
mutual_information_kernel(double p[2][2])
{
  //@assert sum(p) ≈ one(T)     "Entries in p have to sum to 1"
  //@assert all(p .>= zero(T))  "Entries in p have to be >= 0"

  double py[2] = { p[0][0] + p[1][0], p[0][1] + p[1][1] };  // marginal probabilities of y
  double px[2] = { p[0][0] + p[0][1], p[1][0] + p[1][1] };  // marginal probabilities of x

  double M = 0.0;              // mutual information M
  for (int j = 0; j < 2; ++j)  // loop over all entries in p
    for (int i = 0; i < 2; ++i)
    {
      // add entropy only for non-zero entries in p
      if (p[i][j] > 0.0) M += p[i][j] * std::log(p[i][j] / px[i] / py[j]);
    }

  constexpr double base = 2.0;
  M /= std::log(base);  // convert to given base

  return M;
}

/*
  Mutual bitwise information of the elements in input arrays A,B.
  A and B have to be of same size and eltype.
*/
static MutualInformation
mutual_information(float *A, float *B, size_t nelements)
{
  constexpr double confidence = 0.99;
  constexpr bool setZeroInsignificant = true;

  BitpairCounters BC = bitpair_count(A, B, nelements);  // nbits x 2 x 2 array of bitpair counters

  MutualInformation MI;
  for (int i = 0; i < NBITS; ++i) MI.M[i] = 0.0;

  double P[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };  // allocate joint probability mass function

  for (int i = 0; i < NBITS; ++i)  // mutual information for every bit position
  {
    for (int j = 0; j < 2; ++j)  // joint prob mass from counter C
      for (int k = 0; k < 2; ++k) P[j][k] = ((double) BC.C[i][j][k]) / nelements;

    MI.M[i] = mutual_information_kernel(P);
  }

  // remove information that is insignificantly different from a random 50/50
  if (setZeroInsignificant) set_zero_insignificant(MI.M, nelements, confidence);

  return MI;
}

/*
  M = bitinformation(A::AbstractArray{T})

  Bitwise real information content of array `A` calculated from the bitwise mutual information
  in adjacent entries in `A`.
*/
MutualInformation
bitinformation(float *A, size_t n)
{
  //  create a BitArray mask if a masked_value is provided, use === to also allow NaN comparison
  //  isnothing(masked_value) || return bitinformation(A,A .=== masked_value;dim,kwargs...)

  //  A = permute_dim_forward(A,dim)  # Permute A to take adjacent entry in dimension dim
  //  n = A.size(); // n elements in dim

  //  create a two views on A for pairs of adjacent entries, dim is always 1st dimension after permutation
  //  A1view = selectdim(A,1,1:n-1)   # no adjacent entries in A array bounds
  //  A2view = selectdim(A,1,2:n)     # for same indices A2 is the adjacent entry to A1

  return mutual_information(A, A + 1, n - 1);
}

static inline uint32_t
signed_exponent_kernel(uint32_t Auint)
{
  constexpr uint32_t float_sign_mask = 0x80000000;
  constexpr uint32_t float_significand_mask = 0x007fffff;
  constexpr uint32_t float_exponent_mask = 0x7f800000;
  constexpr int32_t float_significand_bits = 23;
  constexpr int32_t float_exponent_bias = 127;

  constexpr auto sfmask = (float_sign_mask | float_significand_mask);
  constexpr auto emask = float_exponent_mask;
  constexpr auto esignmask = (float_sign_mask >> 1);  // exponent sign mask (1st exp bit)

  constexpr auto sbits = float_significand_bits;
  constexpr auto bias = float_exponent_bias;

  auto ui = Auint;
  auto sf = ui & sfmask;                                // sign & fraction bits
  auto e = ((int32_t) ((ui & emask) >> sbits)) - bias;  // de-biased exponent
  auto eabs = (uint32_t) abs(e);                        // magnitude of exponent
  auto esign = (e < 0) ? esignmask : 0;                 // determine sign of exponent
  auto esigned = esign | (eabs << sbits);               // concatentate exponent

  return (sf | esigned);  // concatenate everything back together
}

void
signed_exponent(float *A, size_t n)
{
  uint32_t *Auint = (uint32_t *) A;

  for (size_t i = 0; i < n; ++i) Auint[i] = signed_exponent_kernel(Auint[i]);
}
