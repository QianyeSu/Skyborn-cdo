#ifndef PARSE_LITERALS_H
#define PARSE_LITERALS_H

#include <vector>
#include <string>

int literals_find_datatype(int n, std::vector<std::string> const &literals);
int literal_get_datatype(std::string const &literal);
int literal_to_int(std::string const &literal);
double literal_to_double(std::string const &literal);

#endif
