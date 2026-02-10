/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_fill.h"
#include "cdo_varlist.h"

void
cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData)
{
  VarList varList(p_vlistID);
  p_varData.resize(varList.numVars());
  for (auto const &var : varList.vars) { p_varData[var.ID].resize(var.nlevels * var.gridsize); }
}

void
cdo_fill_ts(int p_vlistID, Varray2D<double> &p_varData, Varray2D<size_t> &p_varNmiss)
{
  VarList varList(p_vlistID);
  p_varData.resize(varList.numVars());
  p_varNmiss.resize(varList.numVars());
  for (auto const &var : varList.vars)
    {
      p_varData[var.ID].resize(var.nlevels * var.gridsize);
      p_varNmiss[var.ID].resize(var.nlevels);
    }
}
