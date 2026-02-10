#ifndef CDO_ZAXIS_H
#define CDO_ZAXIS_H

#include <string>
#include "varray.h"
#include "cdi.h"

int cdo_define_zaxis(std::string const &zaxisfile);
void cdo_set_zaxes(std::string const &zaxisarg);
int zaxis_from_name(std::string const &zaxisname);
int zaxis_from_file(FILE *zfp, const char *filename);
int zaxis_to_ltype(int zaxisID);
double cdo_zaxis_inq_level(int zaxisID, int levelID);
int cdo_zaxis_inq_levels(int zaxisID, double *levels);

void gen_layer_bounds(int nlev, Varray<double> const &levels, Varray<double> &lbounds, Varray<double> &ubounds);
int get_layer_thickness(bool useWeights, bool genBounds, int index, int zaxisID, int nlev, Varray<double> &thickness,
                        Varray<double> &weights);

static inline bool
positive_is_down(int zaxisID)
{
  return (zaxisInqPositive(zaxisID) == 2);
}

#endif
