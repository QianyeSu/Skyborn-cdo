/*
 * #include "util_string.h"
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/

#include "modules.h"

/* removes '-' from operator string and returns copy of parameter operatorCommand */
/***
  vector for library handles for loaded custom modules
  */
std::vector<void *> custom_modules_lib_handles;

/***
  Contains added modules as values and their names as key
  */
