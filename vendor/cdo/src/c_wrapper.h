#ifndef C_WRAPPER_H
#define C_WRAPPER_H

#include <cstdio>
#include <memory>
#include <string_view>

struct FileDestructor
{
  void
  operator()(std::FILE *f) const
  {
    std::fclose(f);
  }
};

const auto c_fopen = [](std::string_view path, std::string_view mode)
{ return std::unique_ptr<std::FILE, FileDestructor>{ std::fopen(path.data(), mode.data()) }; };

#endif
