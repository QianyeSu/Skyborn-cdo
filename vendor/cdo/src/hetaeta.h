#ifndef _HETAETA_H
#define _HETAETA_H

#include "varray.h"

template <typename T>
void hetaeta(bool ltq, int ngp, Vmask const &imiss, int nlev1, const double *ah1, const double *bh1,
             Varray<double> const &fis1, Varray<double> const &ps1, Varray<T> const &t1, Varray<T> const &q1, int nlev2,
             const double *ah2, const double *bh2, Varray<double> const &fis2, Varray<double> &ps2, Varray<T> &t2, Varray<T> &q2,
             int nvars, Varray2D<T> const &vars1, Varray2D<T> &vars2, Varray<double> &scor, Varray<double> &pscor,
             Varray<double> &secor);

#endif /* _HETAETA_H */
