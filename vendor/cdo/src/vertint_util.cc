#include <cstdlib>

#include "cdo_output.h"

bool
getenv_extrapolate(void)
{
  auto extrapolate = false;
  const auto envstr = getenv("EXTRAPOLATE");
  if (envstr && std::isdigit((int) envstr[0]))
  {
    if (atoi(envstr) == 1) extrapolate = true;
    if (extrapolate) cdo_print("Extrapolation of missing values enabled!");
  }

  return extrapolate;
}
