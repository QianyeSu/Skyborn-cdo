#ifndef UTIL_FILES_H
#define UTIL_FILES_H

#include <sys/types.h>
#include <string>
#include <cstdint>

namespace FileUtils
{
bool file_exists(std::string const &fileName);
bool user_file_overwrite(std::string const &fileName);
std::uintmax_t size(std::string const &fileName);
std::string gen_suffix(int filetype, int vlistID, std::string const &refenceName);
}  // namespace FileUtils

#endif
