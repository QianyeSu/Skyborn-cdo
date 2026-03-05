#include "cdo_query.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "param_conversion.h"
#include "util_string.h"
#include "grid_proj.h"
#include <string_view>

static constexpr std::string_view separator("/to/");

static std::vector<int>
decode_int_parameter(std::string const &param)
{
  std::vector<int> values;
  auto n = param.find(separator);
  if (n == std::string::npos) { values.push_back(parameter_to_int(param)); }
  else
  {
    auto v1 = parameter_to_int(param.substr(0, n));
    auto v2 = parameter_to_int(param.substr(n + separator.size()));
    if (v2 < v1) cdo_abort("Second parameter %d muss be greater than first parameter %d!", v2, v1);
    values.push_back(v1);
    auto numVals = (v2 - v1 + 1);
    if (numVals > 1) values.push_back(numVals);
  }

  return values;
}

static std::vector<size_t>
decode_cell_parameter(std::string const &param)
{
  std::vector<size_t> cells;
  auto n = param.find(separator);
  if (n == std::string::npos) { cells.push_back(parameter_to_size_t(param)); }
  else
  {
    auto v1 = parameter_to_size_t(param.substr(0, n));
    auto v2 = parameter_to_size_t(param.substr(n + separator.size()));
    if (v2 < v1) cdo_abort("Second parameter %zu muss be greater than first parameter %zu!", v2, v1);
    cells.push_back(v1);
    auto numVals = (v2 - v1 + 1);
    if (numVals > 1) cells.push_back(numVals);
  }

  return cells;
}

static std::vector<int>
decode_layer_parameter(std::string const &param)
{
  std::vector<int> layers;
  auto n = param.find(separator);
  if (n == std::string::npos) { layers.push_back(parameter_to_int(param)); }
  else
  {
    auto v1 = parameter_to_int(param.substr(0, n));
    auto v2 = parameter_to_int(param.substr(n + separator.size()));
    if (v2 < v1) cdo_abort("Second parameter %zu muss be greater than first parameter %zu!", v2, v1);
    layers.push_back(v1);
    auto numVals = (v2 - v1 + 1);
    if (numVals > 1) layers.push_back(numVals);
  }

  return layers;
}

std::string
set_query_parameter(const KVList &kvlist, struct CdiQuery *query)
{
  std::string path;

  for (auto const &kv : kvlist)
  {
    auto const &key = kv.key;
    int numValues = kv.nvalues;
    if (numValues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

    if (key == "name")
    {
      std::vector<const char *> queryNames(numValues);
      for (int i = 0; i < numValues; ++i) { queryNames[i] = kv.values[i].c_str(); }
      cdiQuerySetNames(query, queryNames.size(), queryNames.data());
    }
    else if (key == "step")
    {
      if (numValues > 1) cdo_abort("Too many values for key=step (maxvalues=1 or range=start/to/end)");
      auto querySteps = decode_int_parameter(kv.values[0]);
      cdiQuerySetSteps(query, querySteps.size(), querySteps.data());
    }
    else if (key == "layer")
    {
      if (numValues > 1) cdo_abort("Too many values for key=layer (maxvalues=1 or range=start/to/end)");
      auto queryLayers = decode_layer_parameter(kv.values[0]);
      cdiQuerySetLayers(query, queryLayers.size(), queryLayers.data());
    }
    else if (key == "cell")
    {
      if (numValues > 1) cdo_abort("Too many values for key=cell (maxvalues=1 or range=start/to/end)");
      auto queryCells = decode_cell_parameter(kv.values[0]);
      cdiQuerySetCells(query, queryCells.size(), queryCells.data());
    }
    else if (key == "hplonlatbox")
    {
      if (numValues != 5) cdo_abort("hplonlatbox needs 5 parameter (hplonlatbox=zoom,lon1,lon2,lat1,lat2)");
      auto zoom = parameter_to_int(kv.values[0]);
      auto xlon1 = parameter_to_double(kv.values[1]);
      auto xlon2 = parameter_to_double(kv.values[2]);
      auto xlat1 = parameter_to_double(kv.values[3]);
      auto xlat2 = parameter_to_double(kv.values[4]);

      if (xlon1 >= xlon2) std::swap(xlon1, xlon2);
      if (xlat1 >= xlat2) std::swap(xlat1, xlat2);

      auto queryCells = healpix_compute_cell_parameter(zoom, xlon1, xlon2, xlat1, xlat2);
      if (queryCells[1] == 0) cdo_abort("hplonlatbox: No grid points found!");

      cdiQuerySetCells(query, queryCells.size(), queryCells.data());
    }
    else if (key == "path")
    {
      path = kv.values[0];
      for (int i = 1; i < numValues; ++i) { path += "," + kv.values[i]; }
    }
    else { cdo_abort("Invalid parameter key >%s<!", key); }
  }

  return path;
}

std::string
set_query_parameter(std::string const &params, struct CdiQuery *query)
{
  auto paramsArgv = split_string(params, ",");

  KVList kvlist;
  kvlist.name = "QUERY";
  if (kvlist.parse_arguments(paramsArgv) != 0) cdo_abort("Parse error!");
  if (Options::cdoVerbose) kvlist.print();

  return set_query_parameter(kvlist, query);
}
