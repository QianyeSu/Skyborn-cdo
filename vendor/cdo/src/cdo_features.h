#ifndef CDO_FEATURES_H
#define CDO_FEATURES_H

#include <string>

namespace cdo
{
namespace features
{

void activate_hdf5_diag();

int print_config(std::string const &option);
void print_features();
void print_libraries();
void print_argument_options();
void print_system_info();
void version();

void print_rusage();
void print_openmp_info();

};  // namespace features
};  // namespace cdo

#endif
