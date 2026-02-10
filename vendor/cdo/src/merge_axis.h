#ifndef MERGEAXIS_H
#define MERGEAXIS_H
#include <vector>

#include "pmlist.h"

struct MergeVarKeys
{
  int vlistID, varID, gridID, projID, zaxisID;
  char datatype = 0;
};

struct MergeVarsOnAxis
{
  std::vector<MergeVarKeys> inputKeys;
  std::vector<double> data;
  KeyValues inputNames;
  MergeVarKeys output;

  void check_axissize_consistency(std::vector<int> axissize);
  std::vector<int> define_new_axes(std::vector<int> axissize);
  void define_var_structure(int vlistID, int ntsteps, std::vector<int> const &axissize);
  void read_cmor_charvar(std::vector<int> const &axissize, int streamID, int oldgridsize);
};
#endif
