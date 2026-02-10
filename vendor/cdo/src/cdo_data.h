#ifndef CDO_DATA_H
#define CDO_DATA_H

#include "varray.h"

struct PackedData;

namespace cdo
{

extern const PackedData topoData;
extern const PackedData tempData;
extern const PackedData maskData;

Varray<float> unpack_data(const PackedData &packedData);

void fill_random(Varray<float> &varray);
void fill_sincos(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals);
void fill_coshill(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals);
void fill_testfield(Varray<float> &varray, Varray<double> const &xvals, Varray<double> const &yvals);

double std_atm_temperatur(double height);
double std_atm_pressure(double height);

}  // namespace cdo

#endif
