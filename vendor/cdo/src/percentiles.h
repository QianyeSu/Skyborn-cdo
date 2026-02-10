/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PERCENTILES_H
#define PERCENTILES_H

#include <string>

void percentile_set_method(std::string const &methodstr);

template <typename T>
double percentile(T *array, size_t len, double pn);

#endif /* PERCENTILES_H */
