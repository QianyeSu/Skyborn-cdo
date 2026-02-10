#ifndef LIBNCL_H
#define LIBNCL_H

#include "cf_interface.h"
/*
void DCFINDIF(double *,double *,int *,double *,
              double *,int *,int *, double *,
              double *,int *,double *,int *);

void DVRFIDF(double *, double *, double *, double *, int, int, double, int,
double *, int *); void DDVFIDF(double *, double *, double *, double *, int, int,
double, int, double *, int *);
*/

// LIBNCL Fortran routines

#ifdef HAVE_CF_INTERFACE

PROTOCCALLSFSUB10(DDVFIDF, ddvfidf, DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT, INT, DOUBLE, INT, DOUBLEV, PINT)
#define DDVFIDF(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)                                                                           \
  CCALLSFSUB10(DDVFIDF, ddvfidf, DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT, INT, DOUBLE, INT, DOUBLEV, PINT, A1, A2, A3, A4, A5, A6, \
               A7, A8, A9, A10)

PROTOCCALLSFSUB10(DVRFIDF, dvrfidf, DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT, INT, DOUBLE, INT, DOUBLEV, PINT)
#define DVRFIDF(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)                                                                           \
  CCALLSFSUB10(DVRFIDF, dvrfidf, DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT, INT, DOUBLE, INT, DOUBLEV, PINT, A1, A2, A3, A4, A5, A6, \
               A7, A8, A9, A10)

#endif

#endif
