/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "remap_utils.h"
#include "param_conversion.h"
#include "cdo_output.h"
#include "util_string.h"
#include "mpim_grid.h"

void
remap_print_info(int operfunc, bool remap_genweights, const RemapGrid &srcGrid, const RemapGrid &tgtGrid, size_t numMissVals,
                 KnnParams const &knnParams)
{
  std::string outStr;
  // clang-format off
  if      (operfunc == REMAPBIL   || operfunc == GENBIL)   outStr = "Bilinear";
  else if (operfunc == REMAPBIC   || operfunc == GENBIC)   outStr = "Bicubic";
  else if (operfunc == REMAPKNN   || operfunc == GENKNN)   outStr = "K-nearest neighbor";
  else if (operfunc == REMAPNN    || operfunc == GENNN)    outStr = "Nearest neighbor";
  else if (operfunc == REMAPDIS   || operfunc == GENDIS)   outStr = "Distance-weighted averaged";
  else if (operfunc == REMAPLAF   || operfunc == GENLAF)   outStr = "YAC largest area fraction";
  else if (operfunc == REMAPCON   || operfunc == GENCON)   outStr = "YAC first order conservative";
  else if (operfunc == REMAPYCON2 || operfunc == GENYCON2) outStr = "YAC second order conservative";
  else if (operfunc == REMAPAVG)                           outStr = "Average";
  else                                                     outStr = "Unknown";
  // clang-format on

  int numNeighbors = knnParams.k;
  if ((operfunc == REMAPDIS || operfunc == GENDIS) && numNeighbors != 4) outStr += " (k=" + std::to_string(numNeighbors) + ")";

  if (operfunc == REMAPKNN || operfunc == GENKNN)
    outStr += " (k=" + std::to_string(numNeighbors) + " weighted=" + weightingMethod_to_string(knnParams.weighted) + ")";

  outStr += remap_genweights ? " weights from " : " remapping from ";
  outStr += srcGrid.name;
  outStr += " (" + std::to_string(srcGrid.dims[0]);
  if (srcGrid.rank == 2) outStr += "x" + std::to_string(srcGrid.dims[1]);
  outStr += ") to ";
  outStr += tgtGrid.name;
  outStr += " (" + std::to_string(tgtGrid.dims[0]);
  if (tgtGrid.rank == 2) outStr += "x" + std::to_string(tgtGrid.dims[1]);
  outStr += ") grid";

  if (numMissVals) outStr += ", with source mask (" + std::to_string(gridInqSize(srcGrid.gridID) - numMissVals) + ")";

  cdo_print(outStr);
}

void
remap_print_warning(std::string const &remapWeightsFile, int operfunc, const RemapGrid &srcGrid, size_t numMissVals)
{
  (void) operfunc;

  std::string outStr = "Remap weights from ";
  outStr += remapWeightsFile;
  outStr += " not used, ";
  outStr += gridNamePtr(gridInqType(srcGrid.gridID));
  outStr += " (" + std::to_string(srcGrid.dims[0]);
  if (srcGrid.rank == 2) outStr += "x" + std::to_string(srcGrid.dims[1]);
  outStr += ") grid";

  if (numMissVals) outStr += " with mask (" + std::to_string(gridInqSize(srcGrid.gridID) - numMissVals) + ")";

  outStr += " not found!";

  cdo_warning(outStr);
}

RemapDefaults
remap_get_params()
{
  RemapDefaults remapDefaults;

  {
    auto envString = getenv_string("REMAP_EXTRAPOLATE");
    if (envString.size())
    {
      // clang-format off
      if      (envString == "ON"  || envString == "on")  remapDefaults.extrapolate = true;
      else if (envString == "OFF" || envString == "off") remapDefaults.extrapolate = false;
      else cdo_warning("Environment variable REMAP_EXTRAPOLATE has wrong value!");
      // clang-format on

      if (Options::cdoVerbose) cdo_print("Extrapolation %s!", remapDefaults.extrapolate ? "enabled" : "disabled");
    }
  }

  {
    auto envString = getenv_string("CDO_REMAP_GENWEIGHTS");
    if (envString.size())
    {
      // clang-format off
      if      (envString == "ON"  || envString == "on")  Options::RemapGenerateWeights = true;
      else if (envString == "OFF" || envString == "off") Options::RemapGenerateWeights = false;
      else cdo_warning("Environment variable CDO_REMAP_GENWEIGHTS has wrong value!");
      // clang-format on

      if (Options::cdoVerbose) cdo_print("Generation of weights %s!", Options::RemapGenerateWeights ? "enabled" : "disabled");
    }
  }

  {
    auto envString = getenv_string("CDO_REMAP_RADIUS");
    if (envString.size())
    {
      auto fval = radius_str_to_deg(envString);
      if (fval < 0.0 || fval > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", "CDO_REMAP_RADIUS", fval);
      cdo_set_search_radius(fval);
      if (Options::cdoVerbose) cdo_print("Set CDO_REMAP_RADIUS to %g", cdo_get_search_radius());
    }
  }

  {
    auto envString = getenv_string("CDO_GRIDSEARCH_RADIUS");
    if (envString.size())
    {
      auto fval = radius_str_to_deg(envString);
      if (fval < 0.0 || fval > 180.0) cdo_abort("%s=%g out of bounds (0-180 deg)!", "CDO_GRIDSEARCH_RADIUS", fval);
      cdo_set_search_radius(fval);
      if (Options::cdoVerbose) cdo_print("Set CDO_GRIDSEARCH_RADIUS to %g", cdo_get_search_radius());
    }
  }

  if (Options::cdoVerbose) cdo_print("Point search radius = %g deg", cdo_get_search_radius());

  {
    auto envString = getenv_string("REMAP_AREA_MIN");
    if (envString.size())
    {
      auto fval = std::stof(envString);
      if (fval > 0.0)
      {
        remapDefaults.fracMin = fval;
        if (Options::cdoVerbose) cdo_print("Set REMAP_AREA_MIN to %g", remapDefaults.fracMin);
      }
    }
  }

  {
    auto envString = getenv_string("MAX_REMAPS");
    if (envString.size())
    {
      auto ival = std::stoi(envString);
      if (ival > 0)
      {
        remapDefaults.maxRemaps = ival;
        if (Options::cdoVerbose) cdo_print("Set MAX_REMAPS to %d", remapDefaults.maxRemaps);
      }
    }
  }

  {
    auto envString = getenv_string("REMAP_MAP3D");
    if (envString.size())
    {
      auto ival = std::stoi(envString);
      if (ival > 0)
      {
        remapDefaults.genMultiWeights = true;
        if (Options::cdoVerbose) cdo_print("Set REMAP_MAP3D to %d", remapDefaults.genMultiWeights);
      }
    }
  }

  return remapDefaults;
}

RemapSwitches
remap_operfunc_to_maptype(int operfunc)
{
  RemapSwitches remapSwitches;
  remapSwitches.remapOrder = 1;

  switch (operfunc)
  {
    case REMAPBIL:
    case GENBIL: remapSwitches.mapType = RemapMethod::BILINEAR; break;
    case REMAPBIC:
    case GENBIC: remapSwitches.mapType = RemapMethod::BICUBIC; break;
    case REMAPKNN:
    case GENKNN:
      remapSwitches.mapType = RemapMethod::KNN;
      remapSwitches.numNeighbors = -1;
      break;
    case REMAPNN:
    case GENNN:
      remapSwitches.mapType = RemapMethod::KNN;
      remapSwitches.numNeighbors = 1;
      break;
    case REMAPDIS:
    case GENDIS:
      remapSwitches.mapType = RemapMethod::KNN;
      if (remapSwitches.numNeighbors == 0) remapSwitches.numNeighbors = 4;
      break;
    case REMAPCON:
    case GENCON: remapSwitches.mapType = RemapMethod::CONSERV; break;
    case REMAPYCON2:
    case GENYCON2:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.remapOrder = 2;
      break;
    case REMAPLAF:
    case GENLAF:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.submapType = SubmapType::LAF;
      break;
    case REMAPAVG:
      remapSwitches.mapType = RemapMethod::CONSERV;
      remapSwitches.submapType = SubmapType::AVG;
      break;
    default: cdo_abort("Unknown mapping method"); break;
  }

  return remapSwitches;
}

NormOpt
remap_get_normOpt(void)
{
  // clang-format off
  NormOpt normOpt(NormOpt::FRACAREA);

  auto envString = getenv_string("CDO_REMAP_NORM");
  if (envString.size())
  {
    if      (envString == "frac" || envString == "fracarea") normOpt = NormOpt::FRACAREA;
    else if (envString == "dest" || envString == "destarea") normOpt = NormOpt::DESTAREA;
    else if (envString == "none") normOpt = NormOpt::NONE;
    else cdo_warning("CDO_REMAP_NORM=%s unsupported!", envString);
  }

  if (Options::cdoVerbose)
  {
    auto outStr = (normOpt == NormOpt::FRACAREA) ? "frac" : (normOpt == NormOpt::DESTAREA) ? "dest" : "none";
    cdo_print("Normalization option: %s", outStr);
  }
  // clang-format on

  return normOpt;
}

void
remap_gen_weights(RemapMethod mapType, KnnParams const &knnParams, RemapType &remap)
{
  // clang-format off
  if      (mapType == RemapMethod::BILINEAR) remap_bilinear_weights(remap.search, remap.vars);
  else if (mapType == RemapMethod::BICUBIC)  remap_bicubic_weights(remap.search, remap.vars);
  else if (mapType == RemapMethod::KNN)      remap_knn_weights(knnParams, remap.search, remap.vars);
  else if (mapType == RemapMethod::CONSERV)  remap_conserv_weights(remap.search, remap.vars);
  // clang-format on
}

std::vector<bool>
remap_set_grids(int vlistID, VarList const &varList, bool findOnlyfirst)
{
  auto numGrids = vlistNumGrids(vlistID);
  std::vector<bool> remapGrids(numGrids, false);
  for (int index = 0; index < numGrids; ++index)
  {
    auto gridID = vlistGrid(vlistID, index);
    auto gridtype = gridInqType(gridID);
    auto hasProjParams = ((gridtype == GRID_PROJECTION) && grid_has_proj_params(gridID));

    if (gridProjIsSupported(gridID) || hasProjParams || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GME
        || gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || gridtype == GRID_GAUSSIAN_REDUCED
        || gridtype == GRID_HEALPIX)
    {
      remapGrids[index] = true;
      if (findOnlyfirst) break;
    }
    else if (!(gridtype == GRID_GENERIC && gridInqSize(gridID) <= 20))
    {
      for (auto const &var : varList.vars)
        if (gridID == var.gridID)
        {
          cdo_abort("Unsupported %s coordinates (Variable: %s)!", gridNamePtr(gridtype), var.name);
          break;
        }
    }
  }

  return remapGrids;
}

int
remap_get_max_maps(int vlistID)
{
  int maxRemaps = 0;

  auto numZaxis = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxis; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    auto numLevels = zaxisInqSize(zaxisID);
    if (numLevels > maxRemaps) maxRemaps = numLevels;
  }

  auto numVars = vlistNvars(vlistID);
  if (numVars > maxRemaps) maxRemaps = numVars;

  maxRemaps++;

  if (Options::cdoVerbose) cdo_print("Set maxRemaps to %d", maxRemaps);

  return maxRemaps;
}
