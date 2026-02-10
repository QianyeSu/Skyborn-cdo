#ifndef CDO_SETTINGS_H
#define CDO_SETTINGS_H

#include <string>

namespace cdo
{

extern int netcdf_hdr_pad;

int evaluate_except_options(std::string const &arg);
int set_feenableexcept(int excepts);

void set_cdi_options();
void set_external_proj_func();

void signal_handler(int signo);

void set_digits(std::string const &arg);
void set_default_filetype(std::string filetypeString);
void set_default_datatype(std::string const &datatypeString);
void set_filterspec(std::string const &arg);
void set_compression_type(std::string const &arg);
void set_chunktype(std::string const &arg);

void evaluate_color_options(std::string const &arg);

void setup_openMP(int ompNumUserRequestedThreads);

};  // namespace cdo

#endif
