/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef MODULES_H
#define MODULES_H

#include <vector>

/***
  type definition for module functions loaded from a custom module
  */
using dyn_oper_t = void (*)(void *arg);

/***
  vector for library handles for loaded custom modules
  */
extern std::vector<void *> custom_modules_lib_handles;

#endif /* MODULES_H */
