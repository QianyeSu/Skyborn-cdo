#ifndef CDO_CMOR_H
#define CDO_CMOR_H

#include <string>
#include <vector>
#include <cdi.h>

struct CmorVar
{
  bool convert = false;
  bool remove = false;
  // missing value
  bool changeMissval = false;
  double missvalOld = 0.0;
  //
  bool applyFactor = false;
  double factor = 0.0;
  //
  bool checkValid = false;
  double valid_min = 0.0;
  double valid_max = 0.0;
  //
  bool check_min_mean_abs = false;
  double ok_min_mean_abs = 0.0;
  //
  bool check_max_mean_abs = false;
  double ok_max_mean_abs = 0.0;
  // units
  bool changeUnits = false;
  char unitsOld[CDI_MAX_NAME] = { 0 };
  char units[CDI_MAX_NAME] = { 0 };
  // varname
  std::string name;
  // converter
  void *ut_converter = nullptr;

  double amean = 0;
  long nvals = 0, n_lower_min = 0, n_greater_max = 0;

  std::string filterSpec;
};

void cmor_check_init(std::vector<CmorVar> &cmorVars);
void cmor_check_eval(int vlistID, std::vector<CmorVar> const &cmorVars);
void cmor_check_prep(CmorVar &var, long gridsize, double missval, const double *const array);

#endif
