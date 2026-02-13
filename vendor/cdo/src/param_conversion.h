#ifndef PARAM_CONVERSION_H
#define PARAM_CONVERSION_H

#include <string>
#include <vector>

long parameter_to_bytes(std::string const &string);

const char *parameter_to_word(const char *cstring);
double parameter_to_double(const char *cstring);
int parameter_to_int(const char *cstring);
long parameter_to_long(const char *cstring);
size_t parameter_to_size_t(const char *cstring);
int parameter_to_intlist(const char *cstring);

std::string const &parameter_to_word(std::string const &string);
double parameter_to_double(std::string const &string);
bool parameter_to_bool(std::string const &string);
int parameter_to_int(std::string const &string);
long parameter_to_long(std::string const &string);
size_t parameter_to_size_t(std::string const &string);
int parameter_to_intlist(std::string const &string);

double radius_str_to_meter(std::string const &string);
double radius_str_to_deg(std::string const &string);

int string_to_param(std::string const &paramstr);
std::string param_to_string(int param);

/* time/date/season converisons */
/* =================================================================================== */
void season_to_months(std::string const &season, int (&imonths)[13]);
double datestr_to_double(std::string const &datestr, int opt);

/* argv conversions */
std::vector<int> cdo_argv_to_intarr(std::vector<std::string> const &argv);
std::vector<double> cdo_argv_to_fltarr(std::vector<std::string> const &argv);

void split_intstring(std::string const &intstr, int &first, int &last, int &inc);
void split_fltstring(std::string const &fltstr, double &first, double &last, double &inc);

template <typename T>
T convert(std::string const &str_value);

#endif
