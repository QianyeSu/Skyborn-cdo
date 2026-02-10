/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_TREND_H
#define FIELD_TREND_H

#include "field.h"

void calc_trend_sum(FieldVector3D &work, Field const &field, double zj, int varID, int levelID);
void sub_trend(double zj, Field &field1, Field const &field2, Field const &field3);
void calc_trend_param(const FieldVector3D &work, Field &field2, Field &field3, int varID, int levelID);

#endif
