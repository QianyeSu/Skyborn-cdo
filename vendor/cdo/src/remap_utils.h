/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_UTILS_H
#define REMAP_UTILS_H

#include "remap.h"

enum
{
  REMAP,
  REMAPBIL,
  REMAPBIC,
  REMAPKNN,
  REMAPNN,
  REMAPDIS,
  REMAPLAF,
  REMAPAVG,
  REMAPCON,
  REMAPYCON2,

  GENBIL,
  GENBIC,
  GENKNN,
  GENNN,
  GENDIS,
  GENLAF,
  GENCON,
  GENYCON2
};

struct RemapDefaults
{
  double fracMin{ 0.0 };
  int maxRemaps{ -1 };
  int extrapolate{ -1 };
  bool genMultiWeights{ false };
};

inline bool
remap_func_is_dist(int operfunc)
{
  return (operfunc == REMAPNN || operfunc == GENNN || operfunc == REMAPDIS || operfunc == GENDIS);
}

void remap_print_info(int operfunc, bool remap_genweights, RemapGrid const &srcGrid, RemapGrid const &tgtGrid, size_t numMissVals,
                      KnnParams const &knnParams);
void remap_print_warning(std::string const &remapWeightsFile, int operfunc, RemapGrid const &srcGrid, size_t numMissVals);

RemapDefaults remap_get_params();
void remap_set_params(RemapDefaults const &remapDefaults);

RemapSwitches remap_operfunc_to_maptype(int operfunc);

NormOpt remap_get_normOpt(void);

void remap_gen_weights(RemapMethod mapType, KnnParams const &knnParams, RemapType &remap);

std::vector<bool> remap_set_grids(int vlistID, VarList const &varList, bool findOnlyfirst = false);

int remap_get_max_maps(int vlistID);

#endif  // REMAP_UTILS_H
