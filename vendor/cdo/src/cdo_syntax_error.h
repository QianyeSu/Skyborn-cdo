#ifndef CDO_SYNTAX_ERROR_H
#define CDO_SYNTAX_ERROR_H

#include <vector>
#include <string>
#include "cdo_exception.h"  //TODO MOVE CdoException  from node to cdoException.h
                            //
struct InternalCdoSyntaxError : CdoException
{
  std::vector<std::string>::const_iterator iter;
  std::string message;
  InternalCdoSyntaxError(std::vector<std::string>::const_iterator p_iter, std::string const &p_msg, std::string const &p_file = "?",
                         std::string const &p_line = "?");
};

struct CdoSyntaxError : InternalCdoSyntaxError
{

  CdoSyntaxError(InternalCdoSyntaxError &e, std::vector<std::string> &p_argv);

  CdoSyntaxError(InternalCdoSyntaxError &e);
  const char *what() const noexcept override;
};

#endif
