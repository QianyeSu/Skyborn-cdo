#ifndef SELBOXINFO_H
#define SELBOXINFO_H

#include <cstdint>
#include <vector>

struct SelboxInfo
{
  std::vector<int64_t> cellIndices;
  long nvals = 0;
  long lat1 = 0, lat2 = 0, lon11 = 0, lon12 = 0, lon21 = 0, lon22 = 0;
  int gridID1 = -1;
  int gridID2 = -1;
  int gridtype = -1;
};

SelboxInfo gen_lonlat_selbox(int argcOffset, int gridID);
SelboxInfo gen_index_selbox(int argcOffset, int gridID);

#endif
