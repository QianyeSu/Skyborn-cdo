#ifndef CHUNKSPEC_H
#define CHUNKSPEC_H

#include <string>

struct ChunkSpec
{
  int x{ 0 };
  int y{ 0 };
  int z{ 0 };
  int t{ 0 };
};

namespace cdo
{

ChunkSpec parse_chunkspec_parameter(std::string const &argument);
ChunkSpec get_chunkspec(int vlistID, int varID);

}  // namespace cdo

#endif
