/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_varlist.h"
#include "compare.h"
#include "cdi_lockedIO.h"
#include "varray.h"

std::vector<bool>
cdo_read_timestepmask(std::string const &maskFileName)
{
  auto streamID = stream_open_read_locked(maskFileName.c_str());
  auto vlistID = streamInqVlist(streamID);
  VarList varList(vlistID);
  auto const &var0 = varList.vars[0];

  if (varList.numVars() > 1) cdo_abort("timestepmask %s contains more than one variable!", maskFileName);
  if (var0.nlevels > 1) cdo_abort("timestepmask %s has more than one level!", maskFileName);
  if (var0.gridsize > 1) cdo_abort("timestepmask %s has more than one gridpoint!", maskFileName);

  auto numSteps = varList.numSteps();
  if (numSteps == -1)
    {
      numSteps = 0;
      while (streamInqTimestep(streamID, numSteps)) numSteps++;

      if (Options::cdoVerbose) cdo_print("%s: counted %i timeSteps in %s", __func__, numSteps, maskFileName);

      streamClose(streamID);
      streamID = stream_open_read_locked(maskFileName.c_str());
    }
  else if (Options::cdoVerbose) { cdo_print("%s: found %i timeSteps in %s", __func__, numSteps, maskFileName); }

  std::vector<bool> imask(numSteps);

  int tsID = 0;
  while (true)
    {
      auto numFields = streamInqTimestep(streamID, tsID);
      if (numFields == 0) break;

      if (numFields != 1) cdo_abort("Internal error; unexprected number of fields!");

      int varID, levelID;
      size_t numMissVals;
      double value;
      streamInqField(streamID, &varID, &levelID);
      streamReadField(streamID, &value, &numMissVals);

      imask[tsID] = !(numMissVals || is_equal(value, 0));

      tsID++;
    }

  streamClose(streamID);

  return imask;
}

static Varray<double>
read_one_field(const char *text, const char *filename)
{
  auto streamID = stream_open_read_locked(filename);
  auto vlistID = streamInqVlist(streamID);
  VarList varList(vlistID);
  auto const &var0 = varList.vars[0];

  if (varList.numVars() > 1) cdo_abort("%s file %s contains more than one variable!", text, filename);
  if (var0.nlevels > 1) cdo_abort("%s file %s has more than one level!", text, filename);
  if (varList.maxFields() != 1) cdo_abort("%s file %s contains more than one field!", text, filename);

  Varray<double> array(var0.gridsize);

  int varID, levelID;
  size_t numMissVals;
  streamInqField(streamID, &varID, &levelID);
  streamReadField(streamID, array.data(), &numMissVals);
  streamClose(streamID);

  return array;
}

std::vector<bool>
cdo_read_mask(const char *maskFileName)
{
  auto array = read_one_field("Mask", maskFileName);
  auto gridsize = array.size();
  std::vector<bool> imask(gridsize);
  for (size_t i = 0; i < gridsize; ++i) imask[i] = is_not_equal(array[i], 0);

  return imask;
}

std::vector<int>
cdo_read_index(const char *indexFileName)
{
  auto array = read_one_field("Index", indexFileName);
  auto gridsize = array.size();
  std::vector<int> index(gridsize);
  for (size_t i = 0; i < gridsize; ++i) index[i] = (int) std::lround(array[i]) - 1;

  return index;
}
