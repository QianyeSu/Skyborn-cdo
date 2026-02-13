/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <fstream>

#include "cdo_options.h"
#include "cdo_output.h"
#include "util_string.h"
#include "varray.h"
#include "dcw_reader.h"
#include "region.h"

static int
read_coords(size_t segmentNo, Varray<double> &xvals, Varray<double> &yvals, std::string const &polyfile, std::ifstream &file)
{
  auto maxVals = xvals.size();

  size_t number = 0, jumpedlines = 0;
  std::string line;
  while (std::getline(file, line))
  {
    if (line[0] == '#' || line[0] == '\0')
    {
      jumpedlines++;
      continue;
    }

    if (number == 0 && line[0] == '>') continue;  // Dump of DCW-GMT
    if (line[0] == '&' || line[0] == '>') break;

    auto lineNo = number + jumpedlines + 1;

    double xcoord, ycoord;
    auto nread = std::sscanf(line.c_str(), "%lf %lf", &xcoord, &ycoord);
    if (nread != 2)
    {
      if (Options::cdoVerbose) cdo_print("nread=%zu, xcoord=%g, ycoord=%g, line=%s\n", nread, xcoord, ycoord, line);
      cdo_abort("Wrong value format in file %s at segment %zu line %zu", polyfile, segmentNo, lineNo);
    }

    if (number >= maxVals) cdo_abort("Too many polygon points (max=%zu)!", maxVals);
    xvals[number] = xcoord;
    yvals[number] = ycoord;
    number++;
  }

  if ((number != 0) && (!(is_equal(xvals[0], xvals[number - 1]) && is_equal(yvals[0], yvals[number - 1]))))
  {
    xvals[number] = xvals[0];
    yvals[number] = yvals[0];
    number++;
  }

  if (Options::cdoVerbose)
    for (size_t i = 0; i < number; ++i) std::fprintf(stderr, "%zu %g %g\n", i + 1, xvals[i], yvals[i]);

  return number;
}

void
read_regions_from_file(std::string const &filename, Regions &regions)
{
  std::ifstream file(filename);
  if (!file.is_open()) cdo_abort("Open failed on: %s\n", filename);

  constexpr size_t maxVals = 1048576;
  Varray<double> xcoords(maxVals), ycoords(maxVals);

  size_t segmentNo = 0;
  size_t n = 0;
  while (true)
  {
    auto segmentSize = read_coords(segmentNo++, xcoords, ycoords, filename, file);
    if (segmentSize == 0) break;
    if (segmentSize < 3) cdo_abort("Too few point for polygon in file %s (Min=3)!", filename);

    auto &x = regions.x;
    auto &y = regions.y;
    auto offset = x.size();
    regions.segmentSize.push_back(segmentSize);
    regions.segmentOffset.push_back(offset);
    regions.numSegments++;
    n += segmentSize;
    x.resize(n);
    y.resize(n);
    for (int i = 0; i < segmentSize; ++i) x[offset + i] = xcoords[i];
    for (int i = 0; i < segmentSize; ++i) y[offset + i] = ycoords[i];
  }

  file.close();
}

void
read_regions_from_dcw(const char *codeNames, Regions &regions)
{
  if (*codeNames == 0) cdo_abort("DCW country code parameter missing!");

  DCW_Lists dcw_lists;
  if (dcw_load_lists(dcw_lists)) cdo_abort("dcw_load_lists() failed!");

  auto codeList = split_string(codeNames, "\\+");
  dcw_sort_countries(dcw_lists);

  codeList = dcw_expand_code_list(dcw_lists, codeList);

  if (codeList.size() == 0) cdo_abort("Empty country code list!");

  auto &lon = regions.x;
  auto &lat = regions.y;
  if (dcw_get_lonlat(dcw_lists, codeList, lon, lat)) cdo_abort("Reading DCW data failed!");

  auto n = lon.size();
  if (n == 0) cdo_abort("Empty country code list!");

  for (size_t i = 0; i < n; ++i)
  {
    if (is_equal(lon[i], 0.0) && is_equal(lat[i], 0.0))
    {
      regions.segmentOffset.push_back(i + 1);
      regions.numSegments++;
    }
  }

  auto numSegments = regions.numSegments;
  if (numSegments == 0) cdo_abort("Empty polygons!");

  for (size_t i = 0; i < numSegments - 1; ++i)
  {
    auto segmentSize = regions.segmentOffset[i + 1] - regions.segmentOffset[i] - 1;
    regions.segmentSize.push_back(segmentSize);
  }
  regions.segmentSize.push_back(regions.x.size() - regions.segmentOffset[numSegments - 1]);
}
