#ifdef _OPENMP
#include <omp.h>
#endif

#include "eof_mode.h"
#include "cdo_output.h"
#include "cdo_options.h"

enum T_EIGEN_MODE
get_eigenmode(void)
{
  enum T_EIGEN_MODE eigen_mode = JACOBI;

  char *envstr = getenv("CDO_SVD_MODE");
  if (envstr)
  {
    if (!std::strncmp(envstr, "danielson_lanczos", 17))
      eigen_mode = DANIELSON_LANCZOS;
    else if (!std::strncmp(envstr, "jacobi", 6))
      eigen_mode = JACOBI;
    else
    {
      cdo_warning("Unknown environmental setting %s for CDO_SVD_MODE. Available options are", envstr);
      cdo_warning("  - 'jacobi' for a one-sided parallelized jacobi algorithm");
      cdo_warning("  - 'danielson_lanzcos' for the D/L algorithm");
    }
  }

  if (Options::cdoVerbose)
    cdo_print("Using CDO_SVD_MODE '%s' from %s", eigen_mode == JACOBI ? "jacobi" : "danielson_lanczos",
              envstr ? "Environment" : " default");

#ifdef _OPENMP
  if (omp_get_max_threads() > 1 && eigen_mode == DANIELSON_LANCZOS)
  {
    cdo_warning("Requested parallel computation with %i Threads ", omp_get_max_threads());
    cdo_warning("  but environmental setting CDO_SVD_MODE causes sequential Singular value decomposition");
  }
#endif

  return eigen_mode;
}

enum T_WEIGHT_MODE
get_weightmode(void)
{
  enum T_WEIGHT_MODE weight_mode = WEIGHT_OFF;

  char *envstr = getenv("CDO_WEIGHT_MODE");
  if (envstr)
  {
    if (!std::strncmp(envstr, "off", 3))
      weight_mode = WEIGHT_OFF;
    else if (!std::strncmp(envstr, "on", 2))
      weight_mode = WEIGHT_ON;
    else
      cdo_warning("Unknown environmental setting %s for CDO_WEIGHT_MODE. Available options are: on/off", envstr);
  }

  if (Options::cdoVerbose)
    cdo_print("Using CDO_WEIGHT_MODE '%s' from %s", weight_mode == WEIGHT_OFF ? "off" : "on", envstr ? "Environment" : " default");

  return weight_mode;
}
