#include <vector>
#include <string>
#include "cdo_syntax_error.h"
#include "mpmo_color.h"

InternalCdoSyntaxError::InternalCdoSyntaxError(std::vector<std::string>::const_iterator p_iter, std::string const &p_msg,
                                               std::string const &p_file, std::string const &p_line)
    : CdoException(p_msg, p_file, p_line), iter(p_iter)
{
}

CdoSyntaxError::CdoSyntaxError(InternalCdoSyntaxError &e, std::vector<std::string> &p_argv)

    : InternalCdoSyntaxError(e.iter, e.message, e.file, e.line)
{
  std::string padding = "";
  std::string errLine = "";
  bool addPadding = true;
  for (auto it = p_argv.begin(); it < p_argv.end(); it++)
    {
      if (it != e.iter)
        {
          if (addPadding) padding += std::string((*it).length() + 1, ' ');
          errLine += *it;
        }
      else
        {
          if (addPadding) padding += std::string((*it).length() / 2, ' ');
          addPadding = false;
          errLine += Red("\033[4m" + *it + "\033[0m");
        }
      errLine += " ";
    }
  message = errLine + " \n" + padding + "^ " + e.what();
#ifdef EXCEPTION_EXTRA_INFO
  message += "\n Thrown from " + file + ": " + line + "\n";
#endif
}

const char *
CdoSyntaxError::what() const noexcept
{
  return message.c_str();
};
