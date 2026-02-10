/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef UTIL_WILDCARDS_H
#define UTIL_WILDCARDS_H

#include <vector>
#include <string>

char *expand_filename(const char *fileName);

std::vector<std::string> expand_path_names(std::vector<std::string> argv);
std::vector<std::string> expand_wild_cards(std::vector<std::string> argv);

bool wildcard_match(std::string const &text, std::string const &pattern);

#endif
