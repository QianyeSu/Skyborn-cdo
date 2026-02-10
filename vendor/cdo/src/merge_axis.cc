/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

#include <cdi.h>

#include <unistd.h>

#include "process_int.h"
#include "varray.h"

#include <mpim_grid.h>

#include "pmlist.h"
#include "merge_axis.h"

/** All Vars in the input file will be merged.
*** Therefore, they need to have the same structure i.e. axis sizes
*** Axissize contains the sizes of the first variable named in inputNames.
**/
void
MergeVarsOnAxis::check_axissize_consistency(std::vector<int> axissize)
{
  for (int i = 1; i < this->inputNames.nvalues; ++i)
  {
    if (axissize[0] != (int) gridInqXsize(this->inputKeys[i].gridID))
      cdo_abort("ERROR (infile: '%s')! In merging variables to a variable with a character coordinate:\n          Size of x-axis: "
                "'%d' of variable '%s'\n          differ from x-axis size of variable '%s': '%d'.",
                cdo_get_stream_name(0), gridInqXsize(this->inputKeys[i].gridID), this->inputNames.values[i].c_str(),
                this->inputNames.values[0].c_str(), axissize[0]);
    if (axissize[1] != (int) gridInqYsize(this->inputKeys[i].gridID))
      cdo_abort("ERROR (infile: '%s')! In merging variables to a variable with a character coordinate:\n          Size of y-axis: "
                "'%d' of variable '%s'\n          differ from y-axis size of variable '%s': '%d'.",
                cdo_get_stream_name(0), gridInqYsize(this->inputKeys[i].gridID), this->inputNames.values[i].c_str(),
                this->inputNames.values[0].c_str(), axissize[1]);
    if (axissize[2] != zaxisInqSize(this->inputKeys[i].zaxisID))
      cdo_abort("ERROR (infile: '%s')! In merging variables to a variable with a character coordinate:\n          Size of z-axis: "
                "'%d' of variable '%s'\n          differ from z-axis size of variable '%s': '%d'.",
                cdo_get_stream_name(0), zaxisInqSize(this->inputKeys[i].zaxisID), this->inputNames.values[i].c_str(),
                this->inputNames.values[0].c_str(), axissize[2]);
  }
}

/** The next function will define output.gridID, output.zaxisID
*** and reset axissize
***
*** The axis that has size=1 will be used as the merge axis
*/
std::vector<int>
MergeVarsOnAxis::define_new_axes(std::vector<int> axissize)
{
  check_axissize_consistency(axissize);

  int nvertex = 0;
  auto projID = this->inputKeys[0].projID;
  std::vector<double> xcell_bounds, ycell_bounds;
  std::vector<double> pxcell_bounds, pycell_bounds;
  std::vector<double> pxcoord_vals, pycoord_vals;
  std::vector<double> xvals, yvals, zvals;
  std::vector<double> subsvals(this->inputNames.nvalues);
  for (int i = 0; i < this->inputNames.nvalues; ++i) subsvals[i] = i + 1;

  auto gridType = gridInqType(this->inputKeys[0].gridID);
  if (axissize[0] == 1)
  {
    xvals = subsvals;
    yvals.resize(axissize[1]);
    zvals.resize(axissize[2]);
    gridInqYvals(this->inputKeys[0].gridID, yvals.data());
    zaxisInqLevels(this->inputKeys[0].zaxisID, zvals.data());
    axissize[0] = this->inputNames.nvalues;
    this->output.gridID = gridCreate(GRID_GENERIC, axissize[0] * axissize[1]);
    auto zaxisType = zaxisInqType(this->inputKeys[0].zaxisID);
    this->output.zaxisID = zaxisCreate(zaxisType, axissize[2]);
  }
  else if (axissize[1] == 1)
  {
    xvals.resize(axissize[0]);
    yvals = subsvals;
    zvals.resize(axissize[2]);
    gridInqXvals(this->inputKeys[0].gridID, xvals.data());
    zaxisInqLevels(this->inputKeys[0].zaxisID, zvals.data());
    axissize[1] = this->inputNames.nvalues;
    this->output.gridID = gridCreate(GRID_GENERIC, axissize[0] * axissize[1]);
    auto zaxisType = zaxisInqType(this->inputKeys[0].zaxisID);
    this->output.zaxisID = zaxisCreate(zaxisType, axissize[2]);
  }
  else if (axissize[2] == 1)
  {
    zvals = subsvals;
    axissize[2] = this->inputNames.nvalues;
    if (axissize[0] && axissize[1])
    {
      nvertex = gridInqNvertex(this->inputKeys[0].gridID);

      if (gridType == GRID_CURVILINEAR)
      {
        xvals.resize(axissize[0] * axissize[1]);
        yvals.resize(axissize[0] * axissize[1]);
        // maximal 4 gridbounds per gridcell permitted
        xcell_bounds.resize(4 * axissize[0] * axissize[1]);
        ycell_bounds.resize(4 * axissize[0] * axissize[1]);
        gridInqXvals(this->inputKeys[0].gridID, xvals.data());
        gridInqYvals(this->inputKeys[0].gridID, yvals.data());
        gridInqYbounds(this->inputKeys[0].gridID, ycell_bounds.data());
        gridInqXbounds(this->inputKeys[0].gridID, xcell_bounds.data());
      }
      else
      {
        xvals.resize(axissize[0]);
        yvals.resize(axissize[1]);
        gridInqXvals(this->inputKeys[0].gridID, xvals.data());
        gridInqYvals(this->inputKeys[0].gridID, yvals.data());
      }
      if (projID != CDI_UNDEFID)
      {
        auto pylength = gridInqYsize(projID);
        auto pxlength = gridInqXsize(projID);
        pxcoord_vals.resize(pxlength);
        pycoord_vals.resize(pylength);
        pxcell_bounds.resize(nvertex * pxlength);
        pycell_bounds.resize(nvertex * pylength);
        gridInqXvals(projID, pxcoord_vals.data());
        gridInqYvals(projID, pycoord_vals.data());
        gridInqYbounds(projID, pycell_bounds.data());
        gridInqXbounds(projID, pxcell_bounds.data());
      }
      if (gridType == GRID_UNSTRUCTURED)
      {
        this->output.gridID = gridCreate(gridType, axissize[0]);
        xcell_bounds.resize(nvertex * axissize[0]);
        ycell_bounds.resize(nvertex * axissize[0]);
        gridInqYbounds(this->inputKeys[0].gridID, ycell_bounds.data());
        gridInqXbounds(this->inputKeys[0].gridID, xcell_bounds.data());
      }
      else
        // this->output.gridID = gridCreate(gridType, axissize[0]);
        this->output.gridID = gridCreate(gridType, axissize[0] * axissize[1]);
    }
    else
      this->output.gridID = gridCreate(gridType, 1);
    this->output.zaxisID = zaxisCreate(ZAXIS_GENERIC, axissize[2]);
  }

  gridDefXsize(this->output.gridID, axissize[0]);
  gridDefYsize(this->output.gridID, axissize[1]);
  if (axissize[0] == 0)
    axissize[0] = 1;
  else
    gridDefXvals(this->output.gridID, xvals.data());
  if (axissize[1] == 0)
    axissize[1] = 1;
  else
    gridDefYvals(this->output.gridID, yvals.data());
  zaxisDefLevels(this->output.zaxisID, zvals.data());

  if (gridType == GRID_UNSTRUCTURED || gridType == GRID_CURVILINEAR)
  {
    auto xunits = cdo::inq_key_string(this->inputKeys[0].gridID, CDI_XAXIS, CDI_KEY_UNITS);
    gridDefXunits(this->output.gridID, xunits.c_str());
    auto yunits = cdo::inq_key_string(this->inputKeys[0].gridID, CDI_YAXIS, CDI_KEY_UNITS);
    gridDefYunits(this->output.gridID, yunits.c_str());
    gridDefNvertex(this->output.gridID, nvertex);
    gridDefXbounds(this->output.gridID, xcell_bounds.data());
    gridDefYbounds(this->output.gridID, ycell_bounds.data());
    if (gridType == GRID_UNSTRUCTURED) { axissize[1] = 1; }
  }

  if (projID != CDI_UNDEFID)
  {
    auto projID2 = gridCreate(GRID_PROJECTION, axissize[0] * axissize[1]);
    gridDefXsize(projID2, axissize[0]);
    gridDefYsize(projID2, axissize[1]);

    gridDefXvals(projID2, pxcoord_vals.data());
    gridDefYvals(projID2, pycoord_vals.data());

    gridDefNvertex(projID2, nvertex);
    gridDefXbounds(projID2, pxcell_bounds.data());
    gridDefYbounds(projID2, pycell_bounds.data());

    grid_copy_names(projID, projID2);
    grid_copy_mapping(projID, projID2);

    gridDefProj(this->output.gridID, projID);
  }

  std::vector<int> newaxissize = { axissize[0], axissize[1], axissize[2] };
  return newaxissize;
}

/**
*** This function will define
*** output.datatype, output.vlistID, output.varID
*** Therefore, a new var is created in output.vlistID
*** It will allocate output.data for ntsteps*axissizes
**/

void
MergeVarsOnAxis::define_var_structure(int vlistID, int ntsteps, std::vector<int> const &axissize)
{
  this->output.vlistID = vlistID;
  this->output.varID = vlistDefVar(this->output.vlistID, this->output.gridID, this->output.zaxisID, TIME_VARYING);
  auto oldcode = vlistInqVarCode(this->output.vlistID, this->inputKeys[0].varID);
  vlistDefVarCode(this->output.vlistID, this->inputKeys[0].varID, 1);
  vlistDefVarCode(this->output.vlistID, this->output.varID, oldcode);
  cdiDefKeyString(this->output.vlistID, this->inputKeys[0].varID, CDI_KEY_NAME, "ChangedForMap");
  cdiDefKeyString(this->output.vlistID, this->output.varID, CDI_KEY_NAME, this->inputNames.values[0].c_str());
  auto datatype = vlistInqVarDatatype(this->output.vlistID, this->inputKeys[0].varID);
  vlistDefVarDatatype(this->output.vlistID, this->output.varID, datatype);
  auto mv = vlistInqVarMissval(this->output.vlistID, this->inputKeys[0].varID);
  vlistDefVarMissval(this->output.vlistID, this->output.varID, mv);
  this->output.datatype = (datatype == CDI_DATATYPE_FLT64) ? 'd' : 'f';
  this->data.resize(ntsteps * axissize[0] * axissize[1] * axissize[2]);
}
/**
*** This function reads data from the inputKeys.varID from a streamID
*** to data.
*** The index for data is defined for CMOR
**/

void
MergeVarsOnAxis::read_cmor_charvar(std::vector<int> const &axissize, int streamID, int oldgridsize)
{
  Varray<double> buffer_old(oldgridsize);

  auto gridType = gridInqType(this->inputKeys[0].gridID);
  auto ztype = zaxisInqType(this->inputKeys[0].zaxisID);
  int tsID = 0;
  while (true)
  {
    auto numFields = streamInqTimestep(streamID, tsID);
    if (numFields == 0) break;

    while (numFields--)
    {
      int varIDrw, levelIDrw;
      size_t numMissVals;
      streamInqField(streamID, &varIDrw, &levelIDrw);
      for (int i = 0; i < this->inputNames.nvalues; ++i)
        if (varIDrw == this->inputKeys[i].varID)
        {
          streamReadField(streamID, buffer_old.data(), &numMissVals);
          int newIndex;
          for (int j = 0; j < oldgridsize; ++j)
          {
            // (lev x lat, basin )

            /* (lev x lat, basin )
                        newIndex = j * levdim + levelID; */
            /* Use this index because z-axis is registered at the end and rearranged by CMOR */
            if (oldgridsize == axissize[0] * axissize[1])
            {
              if ((gridType == GRID_UNSTRUCTURED || gridType == GRID_CURVILINEAR) && ztype != ZAXIS_HYBRID)
                newIndex = tsID * axissize[2] * axissize[0] * axissize[1] + axissize[0] * axissize[1] * i + j;
              else
                newIndex = tsID * axissize[2] * axissize[0] * axissize[1] + j * axissize[2] + i;
            }
            else if (axissize[0] == this->inputNames.nvalues)
              newIndex
                  = tsID * axissize[2] * axissize[0] * axissize[1] + i * axissize[1] * axissize[2] + j * axissize[2] + levelIDrw;
            else
              newIndex = tsID * axissize[2] * axissize[0] * axissize[1] + levelIDrw * axissize[0] * axissize[1]
                         + i * axissize[0] * axissize[1] / oldgridsize + j;

            this->data[newIndex] = buffer_old[j];
          }
        }
    }
    tsID++;
  }
}
