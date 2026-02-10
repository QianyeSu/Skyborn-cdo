#ifndef TABLE_H
#define TABLE_H

#include <string>

namespace cdo
{

int define_table(std::string const &tablearg);
std::string predefined_tables(int p_padding);

}  // namespace cdo

#endif
