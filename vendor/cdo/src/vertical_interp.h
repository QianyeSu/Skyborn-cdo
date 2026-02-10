/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef VERTICAL_INTERP_H
#define VERTICAL_INTERP_H

#include <cstddef>
#include "varray.h"

void height_to_pressure(const double *heightLevels, double *pressureLevels, int numLevels);

template <typename T>
void vct_to_hybrid_pressure(T *fullPress, T *halfPress, Varray<double> const &vct, const T *ps, long nhlev, long ngp);

void extrapolate_P(double *slp, const double *halfPress, const double *fullPress, const double *geop, const double *temp, long ngp);

template <typename T>
void vertical_interp_T(const T *geop, const T *gt, T *pt, const T *fullPress, const T *halfPress, const int *vertIndex,
                       const double *plev, long nplev, long ngp, long nhlev, double missval);

template <typename T>
void vertical_interp_Z(const T *geop, const T *gz, T *pz, const T *fullPress, const T *halfPress, const int *vertIndex, const T *gt,
                       const double *plev, long nplev, long ngp, long nhlev, double missval);

template <typename T>
void vertical_interp_X(const T *arrayIn3D, T *arrayOut3D, const T *levels3D, const int *vertIndex3D, const double *levels,
                       long numLevels, long numGP, long nhlev, double missval);

template <typename T>
void gen_vert_index(int *vertIndex, const double *plev, const T *fullPress, long ngp, int nplev, int nhlev, bool lreverse = false);

template <typename T>
void gen_vert_index_mv(int *vertIndex, const double *plev, long ngp, int nplev, const T *psProg, size_t *pnumMissVals,
                       bool lreverse = false);

#endif /* VERTICAL_INTERP_H */
