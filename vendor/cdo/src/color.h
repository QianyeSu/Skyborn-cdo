#ifndef COLOR_H
#define COLOR_H

#include <cstdio>
#include <vector>

struct LUT
{
  double z_low = 0.0, z_high = 0.0, i_dz = 0.0;
  int rgb_low[3] = { 0, 0, 0 }, rgb_high[3] = { 0, 0, 0 }, rgb_diff[3] = { 0, 0, 0 };
  int annot = 0;
  int skip = 0;
};

struct BFN_COLOR
{ /* For back-, fore-, and nan-colors */
  int rgb[3] = { 0, 0, 0 };
  int skip = 1;
};

struct CPT
{
  int ncolors = 0;
  std::vector<LUT> lut;
  BFN_COLOR bfn[3];
};

int cpt_read(std::FILE *fp, CPT *cpt);
int cpt_write(std::FILE *fp, const CPT &cpt);
int cpt_write_c(std::FILE *fp, const CPT &cpt, const char *name);

#endif /* COLOR_H */
