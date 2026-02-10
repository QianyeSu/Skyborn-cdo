#ifndef EIGEN_SOLUTION_H
#define EIGEN_SOLUTION_H

#include "varray.h"

void eigen_solution_of_symmetric_matrix(Varray2D<double> &a, Varray<double> &eig_val, size_t n, const char *prompt);

// make parallel eigen solution accessible for eigen value computation in EOF3d.c
void parallel_eigen_solution_of_symmetric_matrix(Varray2D<double> &M, Varray<double> &A, size_t n, const char func[]);

#endif
