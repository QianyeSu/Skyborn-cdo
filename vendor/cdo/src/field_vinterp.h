/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_VINTERP_H
#define FIELD_VINTERP_H

#include "field.h"
#include "vertical_interp.h"

void gen_vert_index(std::vector<int> &vertIndex, Varray<double> const &plev, Field3D &full_level, size_t gridsize,
                    bool lreverse = false);

void gen_vert_index_mv(std::vector<int> &vertIndex, Varray<double> const &plev, size_t gridsize, Field &level0,
                       Varray<size_t> &pnumMissVals, bool lreverse = false);

void vertical_interp_T(size_t nlevels, Field3D const &full_level, Field3D const &half_level, Field3D const &field1, Field3D &field2,
                       Field const &sgeopot, std::vector<int> const &vertIndex, Varray<double> const &plev, size_t gridsize);

void vertical_interp_Z(size_t nlevels, Field3D const &full_level, Field3D const &half_level, Field3D const &field1, Field3D &field2,
                       Field3D const &temp, Field const &sgeopot, std::vector<int> const &vertIndex, Varray<double> const &plev,
                       size_t gridsize);

void vertical_interp_X(Field3D const &levels3D, Field3D const &field1, Field3D &field2, std::vector<int> const &vertIndex3D,
                       Varray<double> const &levels2, size_t gridsize);

#endif
