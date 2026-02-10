#ifndef GAUSSIAN_LATITUDES_H_
#define GAUSSIAN_LATITUDES_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

  void gaussian_latitudes(size_t nlats, double *latitudes, double *weights);
  bool is_gaussian_latitudes(size_t nlats, const double *latitudes);

#ifdef __cplusplus
}
#endif

#endif /* GAUSSIAN_LATITUDES_H_ */
