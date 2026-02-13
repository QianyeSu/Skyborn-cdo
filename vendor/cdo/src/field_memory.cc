/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "field.h"
#include "cdo_options.h"

static void
field2D_init_kernel(FieldVector2D &field2D, VarList const &varList, int ptype, bool doFill, double fillValue)
{
  auto allocateData = (ptype & FIELD_VEC);
  auto numVars = varList.numVars();
  field2D.resize(numVars);

  for (auto const &var : varList.vars)
  {
    auto gridSize = var.gridsize;
    auto size = gridSize * var.nwpv;
    auto numLevels = var.nlevels;
    auto memType = (ptype & FIELD_NAT) ? var.memType : ((ptype & FIELD_FLT) ? MemType::Float : MemType::Double);

    field2D[var.ID].resize(numLevels);

    for (int levelID = 0; levelID < numLevels; ++levelID)
    {
      auto &field = field2D[var.ID][levelID];
      field.nwpv = var.nwpv;
      field.grid = var.gridID;
      field.size = size;
      field.memType = memType;
      field.missval = var.missval;

      if (allocateData)
      {
        if (memType == MemType::Float) { doFill ? field.resizef(size, (float) fillValue) : field.resizef(size); }
        else { doFill ? field.resize(size, fillValue) : field.resize(size); }
      }
    }
  }
}

void
field2D_init(FieldVector2D &field2D, VarList const &varList)
{
  field2D_init_kernel(field2D, varList, 0, false, 0);
}

void
field2D_init(FieldVector2D &field2D, VarList const &varList, int ptype)
{
  field2D_init_kernel(field2D, varList, ptype, false, 0);
}

void
field2D_init(FieldVector2D &field2D, VarList const &varList, int ptype, double fillValue)
{
  field2D_init_kernel(field2D, varList, ptype, true, fillValue);
}

static void
field1Dvars_init_kernel(FieldVector &field1D, VarList const &varList, int ptype, bool lfill, double fillValue)
{
  auto allocateData = (ptype & FIELD_VEC);
  auto numVars = varList.numVars();
  field1D.resize(numVars);

  for (auto const &var : varList.vars)
  {
    auto gridSize = var.gridsize;
    auto size = gridSize * var.nwpv;
    auto dataType = var.dataType;
    auto memType = (ptype & FIELD_FLT) ? MemType::Float : MemType::Double;
    if (ptype & FIELD_NAT)
    {
      if (Options::CDO_Memtype == MemType::Native)
        memType = (dataType == CDI_DATATYPE_FLT32 || dataType == CDI_DATATYPE_CPX32) ? MemType::Float : MemType::Double;
      else
        memType = Options::CDO_Memtype;
    }

    auto &field = field1D[var.ID];

    field.nwpv = var.nwpv;
    field.grid = var.gridID;
    field.size = size;
    field.memType = memType;
    field.missval = var.missval;

    if (allocateData)
    {
      if (memType == MemType::Float)
      {
        if (lfill)
          field.resizef(size, (float) fillValue);
        else
          field.resizef(size);
      }
      else
      {
        if (lfill)
          field.resize(size, fillValue);
        else
          field.resize(size);
      }
    }
  }
}

void
field1Dvars_init(FieldVector &field1D, VarList const &varList)
{
  field1Dvars_init_kernel(field1D, varList, 0, false, 0);
}

void
field1Dvars_init(FieldVector &field1D, VarList const &varList, int ptype)
{
  field1Dvars_init_kernel(field1D, varList, ptype, false, 0);
}

static void
field1Dlevels_init_kernel(FieldVector &field1D, VarList const &varList, int ptype, bool lfill, double fillValue)
{
  auto allocateData = (ptype & FIELD_VEC);

  auto const &var = varList.vars[0];
  auto gridSize = var.gridsize;
  auto size = gridSize * var.nwpv;
  auto numLevels = var.nlevels;
  auto dataType = var.dataType;
  auto memType = (ptype & FIELD_FLT) ? MemType::Float : MemType::Double;
  if (ptype & FIELD_NAT)
  {
    if (Options::CDO_Memtype == MemType::Native)
      memType = (dataType == CDI_DATATYPE_FLT32 || dataType == CDI_DATATYPE_CPX32) ? MemType::Float : MemType::Double;
    else
      memType = Options::CDO_Memtype;
  }

  field1D.resize(numLevels);

  for (int levelID = 0; levelID < numLevels; ++levelID)
  {
    auto &field = field1D[levelID];

    field.nwpv = var.nwpv;
    field.grid = var.gridID;
    field.size = size;
    field.memType = memType;
    field.missval = var.missval;

    if (allocateData)
    {
      if (memType == MemType::Float)
      {
        if (lfill)
          field.resizef(size, (float) fillValue);
        else
          field.resizef(size);
      }
      else
      {
        if (lfill)
          field.resize(size, fillValue);
        else
          field.resize(size);
      }
    }
  }
}

void
field1Dlevels_init(FieldVector &field1D, VarList const &varList)
{
  field1Dlevels_init_kernel(field1D, varList, 0, false, 0);
}

void
field1Dlevels_init(FieldVector &field1D, VarList const &varList, int ptype)
{
  field1Dlevels_init_kernel(field1D, varList, ptype, false, 0);
}
