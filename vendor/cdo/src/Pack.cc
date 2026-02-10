/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Pack    pack         Pack
*/

#include <climits>
#include <fstream>

#include <cdi.h>

#include "process_int.h"
#include "datetime.h"
#include "cdo_default_values.h"
#include "field_functions.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "util_string.h"

namespace
{
struct Parameter
{
  double add_offset{ 0 };
  double scale_factor{ 1 };
  std::string filename{};
  bool printParam{ false };
};
}  // namespace

namespace
{
struct PackEntry
{
  double add_offset{ 0 };
  double scale_factor{ 1 };
  std::string name{};
};
}  // namespace

static std::vector<PackEntry>
read_params_from_file(std::string const &filename)
{
  std::vector<PackEntry> packList;

  if (filename.size())
  {
    std::ifstream file(filename);
    if (!file.is_open()) cdo_abort("Open failed on: %s\n", filename);

    std::string line;
    while (std::getline(file, line))
    {
      auto keyValuesLine = split_string(line, " +");

      KVList kvlist;
      if (kvlist.parse_arguments(keyValuesLine) != 0) cdo_abort("Parse error!");
      if (Options::cdoVerbose) kvlist.print();

      PackEntry packEntry;
      for (auto const &kv : kvlist)
      {
        auto const &key = kv.key;
        if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
        if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
        auto const &value = kv.values[0];

        // clang-format off
        if      (key == "name")         packEntry.name = parameter_to_word(value);
        else if (key == "add_offset")   packEntry.add_offset = parameter_to_double(value);
        else if (key == "scale_factor") packEntry.scale_factor = parameter_to_double(value);
        else cdo_abort("Invalid parameter key >%s<!", key);
        // clang-format on
      }

      if (!packEntry.name.empty()) packList.push_back(packEntry);
    }

    file.close();
  }

  return packList;
}

static int
get_type_values(int datatype, double &tmin, double &tmax, double &tmv)
{
  int status = 0;

  // clang-format off
  switch (datatype)
  {
    case CDI_DATATYPE_INT8:   tmv = -SCHAR_MAX; tmin = -SCHAR_MAX + 1;  tmax = SCHAR_MAX;     break;
    case CDI_DATATYPE_UINT8:  tmv =  UCHAR_MAX; tmin = 0;               tmax = UCHAR_MAX - 1; break;
    case CDI_DATATYPE_INT16:  tmv =  -SHRT_MAX; tmin = -SHRT_MAX + 1;   tmax = SHRT_MAX;      break;
    case CDI_DATATYPE_UINT16: tmv =  USHRT_MAX; tmin = 0;               tmax = USHRT_MAX - 1; break;
    case CDI_DATATYPE_INT32:  tmv =   -INT_MAX; tmin = -INT_MAX + 1;    tmax = INT_MAX;       break;
    case CDI_DATATYPE_UINT32: tmv =   UINT_MAX; tmin = 0;               tmax = UINT_MAX - 1;  break;
    default: status = 1; break;
  }
  // clang-format on

  return status;
}

static int
compute_scale_and_offset(int datatype, double fmin, double fmax, double &scaleFactor, double &addOffset)
{
  scaleFactor = 1.0;
  addOffset = 0.0;

  double tmin, tmax, tmv;
  if (get_type_values(datatype, tmin, tmax, tmv)) return 1;

  if (is_not_equal(fmin, fmax))
  {
    scaleFactor = (fmax - fmin) / (tmax - tmin);
    addOffset = ((fmax + fmin) - scaleFactor * (tmin + tmax)) / 2;
  }

  return 0;
}

template <typename T>
static void
field_change_missval(Varray<T> &v, double mv1, double mv2)
{
  T missval1 = mv1;
  T missval2 = mv2;
  auto n = v.size();
  for (size_t i = 0; i < n; ++i)
    if (fp_is_equal(v[i], missval1)) v[i] = missval2;
}

static void
field_change_missval(Field &field, double missval1, double missval2)
{
  auto func = [&](auto &v) { field_change_missval(v, missval1, missval2); };
  field_operation(func, field);
}

static void
print_parameter(const Parameter &params)
{
  std::stringstream outbuffer;

  outbuffer << "add_offset=" << params.add_offset;
  outbuffer << ", scale_factor=" << params.scale_factor;
  outbuffer << ", printpack=" << params.printParam;
  outbuffer << ", filename=" << params.filename;

  cdo_verbose("%s", outbuffer.str());
}

static void
print_pack_params(std::string const &name, double addOffset, double scaleFactor)
{
  fprintf(stdout, "name=%s  add_offset=%.9g  scale_factor=%.9g\n", name.c_str(), addOffset, scaleFactor);
}

static Parameter
get_parameter()
{
  Parameter params;

  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "add_offset")   params.add_offset = parameter_to_double(value);
      else if (key == "scale_factor") params.scale_factor = parameter_to_double(value);
      else if (key == "printparam")   params.printParam = parameter_to_bool(value);
      else if (key == "filename")     params.filename = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

class Pack : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Pack",
    .operators = { { "pack", PackHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static auto registration = RegisterEntry<Pack>(module);

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  Parameter params;

  int datatype = CDI_DATATYPE_INT16;

  VarList varList1;

private:
  void
  run_method1()
  {
    FieldVector3D varsData;
    DateTimeList dtlist;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      constexpr size_t NALLOC_INC = 1024;
      if ((size_t) tsID >= varsData.size()) varsData.resize(varsData.size() + NALLOC_INC);

      dtlist.taxis_inq_timestep(taxisID1, tsID);

      field2D_init(varsData[tsID], varList1);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto &field = varsData[tsID][varID][levelID];
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);
      }

      tsID++;
    }

    auto nts = tsID;

    constexpr double undefValue = 1.0e300;

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];

      double fmin = undefValue, fmax = -undefValue;
      size_t numMissValspv = 0;

      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        for (int t = 0; t < nts; ++t)
        {
          if (t > 0 && var.isConstant) continue;

          auto const &field = varsData[t][varID][levelID];
          auto numMissVals = field.numMissVals;

          if (numMissVals) numMissValspv += numMissVals;

          if (numMissVals < var.gridsize)
          {
            auto mm = field_min_max(field);
            fmin = std::min(fmin, mm.min);
            fmax = std::max(fmax, mm.max);
          }
        }
      }

      vlistDefVarDatatype(vlistID2, varID, datatype);

      auto hasValidData = (is_not_equal(fmin, undefValue) && is_not_equal(fmax, -undefValue));

      if (numMissValspv)
      {
        double tmin, tmax, missval2;
        if (!get_type_values(datatype, tmin, tmax, missval2))
        {
          vlistDefVarMissval(vlistID2, varID, missval2);

          if (!(missval2 < tmin || missval2 > tmax))
            cdo_warning("new missing value %g is inside data range (%g - %g)!", missval2, tmin, tmax);

          for (int levelID = 0; levelID < var.nlevels; ++levelID)
          {
            for (int t = 0; t < nts; ++t)
            {
              if (t > 0 && var.isConstant) continue;

              auto &field = varsData[t][varID][levelID];
              if (field.numMissVals) field_change_missval(field, var.missval, missval2);
            }
          }
        }
      }

      if (hasValidData)
      {
        double scaleFactor, addOffset;
        if (!compute_scale_and_offset(datatype, fmin, fmax, scaleFactor, addOffset))
        {
          auto memTypeIsFloat = (var.memType == MemType::Float);
          cdiDefKeyFloat(vlistID2, varID, CDI_KEY_ADDOFFSET, memTypeIsFloat ? (float) addOffset : addOffset);
          cdiDefKeyFloat(vlistID2, varID, CDI_KEY_SCALEFACTOR, memTypeIsFloat ? (float) scaleFactor : scaleFactor);

          if (params.printParam) print_pack_params(var.name, addOffset, scaleFactor);
        }
      }
    }

    cdo_def_vlist(streamID2, vlistID2);

    for (tsID = 0; tsID < nts; ++tsID)
    {
      dtlist.taxis_def_timestep(taxisID2, tsID);
      cdo_def_timestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto const &var = varList1.vars[varID];
        if (tsID > 0 && var.isConstant) continue;
        for (int levelID = 0; levelID < var.nlevels; ++levelID)
        {
          auto &field = varsData[tsID][varID][levelID];
          if (field.hasData())
          {
            cdo_def_field(streamID2, varID, levelID);
            cdo_write_field(streamID2, field);
          }
        }
      }
    }
  }

  void
  run_method2()
  {
    const auto packList = read_params_from_file(params.filename);

    auto numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];

      for (auto const &packEntry : packList)
      {
        if (var.name == packEntry.name)
        {
          vlistDefVarDatatype(vlistID2, varID, datatype);

          auto scaleFactor = packEntry.scale_factor;
          auto addOffset = packEntry.add_offset;
          auto memTypeIsFloat = (var.memType == MemType::Float);
          cdiDefKeyFloat(vlistID2, varID, CDI_KEY_ADDOFFSET, memTypeIsFloat ? (float) addOffset : addOffset);
          cdiDefKeyFloat(vlistID2, varID, CDI_KEY_SCALEFACTOR, memTypeIsFloat ? (float) scaleFactor : scaleFactor);

          if (params.printParam) print_pack_params(var.name, addOffset, scaleFactor);

          break;
        }
      }
    }

    cdo_def_vlist(streamID2, vlistID2);

    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        field.init(varList1.vars[varID]);
        cdo_read_field(streamID1, field);

        cdo_def_field(streamID2, varID, levelID);
        cdo_write_field(streamID2, field);
      }

      tsID++;
    }
  }

public:
  void
  init() override
  {
    params = get_parameter();
    if (Options::cdoVerbose) print_parameter(params);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);

    varList1 = VarList(vlistID1);

    if (CdoDefault::DataType != CDI_UNDEFID)
    {
      // 32-bit float rounding error with CDI_DATATYPE_INT32|CDI_DATATYPE_UINT32
      if (CdoDefault::DataType == CDI_DATATYPE_FLT64 || CdoDefault::DataType == CDI_DATATYPE_FLT32
          || CdoDefault::DataType == CDI_DATATYPE_INT32 || CdoDefault::DataType == CDI_DATATYPE_UINT32)
      {
        cdo_warning("Changed default output datatype to int16");
        CdoDefault::DataType = datatype;
      }
      else { datatype = CdoDefault::DataType; }
    }
  }

  void
  run() override
  {
    if (params.filename.empty())
      run_method1();
    else
      run_method2();
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
