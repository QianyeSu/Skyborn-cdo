#ifndef STATISTIC_H
#define STATISTIC_H

namespace cdo
{

double normal_density(double x);
double normal(double x);
double normal_inv(double p);
double student_t_density(double n, double x);
double student_t(double n, double x);
double student_t_inv(double n, double p);
double chi_square_density(double n, double x);
double chi_square(double n, double x);
double chi_square_inv(double n, double p);
void chi_square_constants(double n, double p, double *c1, double *c2);
double beta_distr_density(double a, double b, double x);
double beta_distr(double a, double b, double x);
double beta_distr_inv(double a, double b, double p);
void beta_distr_constants(double a, double b, double p, double *c1, double *c2);
double fisher(double m, double n, double x);

}  // namespace cdo

#endif
