/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef EXPR_FUN_H
#define EXPR_FUN_H

#include <cstddef>
#include "field.h"
#include "expr.h"

nodeType *expr_con_con(int oper, const nodeType *p1, const nodeType *p2);
void oper_expr_con_var(int oper, bool hasMV, size_t n, double mv, double *odat, double cval, const double *idat);
void oper_expr_var_con(int oper, bool hasMV, size_t n, double mv, double *odat, const double *idat, double cval);
void oper_expr_var_var(int oper, bool hasMV, size_t n, double mv1, double mv2, double *odat, const double *idat1, double *idat2);

void fld_field_init(Field &field, size_t numMissVals, double missval, size_t ngp, double const *array, double const *w);
void vert_weights(int zaxisID, size_t nlev, Varray<double> &weights);

#endif
