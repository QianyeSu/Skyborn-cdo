/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "field_vinterp.h"

void
gen_vert_index(std::vector<int> &vertIndex, Varray<double> const &plev, Field3D &fullLevel, size_t gridsize, bool lreverse)
{
  auto nplev = plev.size();
  auto nhlevf = fullLevel.nlevels;
  auto func = [&](auto &v) { gen_vert_index(vertIndex.data(), plev.data(), v.data(), gridsize, nplev, nhlevf, lreverse); };
  field_operation(func, fullLevel);
}

void
gen_vert_index_mv(std::vector<int> &vertIndex, Varray<double> const &plev, size_t gridsize, Field &level0,
                  Varray<size_t> &pnumMissVals, bool lreverse)
{
  auto nplev = plev.size();
  auto func = [&](auto &v)
  { gen_vert_index_mv(vertIndex.data(), plev.data(), gridsize, nplev, v.data(), pnumMissVals.data(), lreverse); };
  field_operation(func, level0);
}

void
vertical_interp_T(size_t nlevels, Field3D const &fullLevel, Field3D const &halfLevel, Field3D const &field1, Field3D &field2,
                  Field const &sgeopot, std::vector<int> const &vertIndex, Varray<double> const &plev, size_t gridsize)
{
  auto nplev = plev.size();
  auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_T(sgeopot.vec_f.data(), field1.vec_f.data(), field2.vec_f.data(), fullLevel.vec_f.data(),
                      halfLevel.vec_f.data(), &vertIndex[0], plev.data(), nplev, gridsize, nlevels, missval);
  else
    vertical_interp_T(sgeopot.vec_d.data(), field1.vec_d.data(), field2.vec_d.data(), fullLevel.vec_d.data(),
                      halfLevel.vec_d.data(), &vertIndex[0], plev.data(), nplev, gridsize, nlevels, missval);
}

void
vertical_interp_Z(size_t nlevels, Field3D const &fullLevel, Field3D const &halfLevel, Field3D const &field1, Field3D &field2,
                  Field3D const &temp, Field const &sgeopot, std::vector<int> const &vertIndex, Varray<double> const &plev,
                  size_t gridsize)
{
  auto nplev = plev.size();
  auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_Z(sgeopot.vec_f.data(), field1.vec_f.data(), field2.vec_f.data(), fullLevel.vec_f.data(),
                      halfLevel.vec_f.data(), &vertIndex[0], temp.vec_f.data(), plev.data(), nplev, gridsize, nlevels, missval);
  else
    vertical_interp_Z(sgeopot.vec_d.data(), field1.vec_d.data(), field2.vec_d.data(), fullLevel.vec_d.data(),
                      halfLevel.vec_d.data(), &vertIndex[0], temp.vec_d.data(), plev.data(), nplev, gridsize, nlevels, missval);
}

void
vertical_interp_X(Field3D const &levels3D, Field3D const &field1, Field3D &field2, std::vector<int> const &vertIndex3D,
                  Varray<double> const &levels2, size_t gridsize)
{
  auto numLevels2 = levels2.size();
  auto missval = field1.missval;
  if (field1.memType == MemType::Float)
    vertical_interp_X(field1.vec_f.data(), field2.vec_f.data(), levels3D.vec_f.data(), vertIndex3D.data(), levels2.data(),
                      numLevels2, gridsize, levels3D.nlevels, missval);
  else
    vertical_interp_X(field1.vec_d.data(), field2.vec_d.data(), levels3D.vec_d.data(), vertIndex3D.data(), levels2.data(),
                      numLevels2, gridsize, levels3D.nlevels, missval);
}
