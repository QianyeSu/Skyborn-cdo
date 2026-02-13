/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Sinfo      sinfo           Short dataset informationl
*/

#include "cdi.h"
#include "julian_date.h"

#include <cstring>
#include <string>
#include "cdo_options.h"
#include "printinfo.h"
#include "mpmo_color.h"
#include "process_int.h"
#include "util_string.h"
#include "datetime.h"
#include "cdo_default_values.h"
#include "cdo_cdi_wrapper.h"

enum
{
  func_generic,
  func_param,
  func_name,
  func_code
};

static const char *
memtype_to_cstr(MemType memType)
{
  return (memType == MemType::Double) ? "F64" : "F32";
}

static const char *
num_values_to_byte_cstr(size_t numValues)
{
  static char cstring[128];
  cstring[0] = 0;

  size_t memUnitsIdx = 0;
  constexpr std::array<const char *, 6> memUnitsList = { "Bytes", "Kbytes", "Mbytes", "Gbytes", "Tbytes", "Pbytes" };
  while (numValues > 9999 && memUnitsIdx < memUnitsList.size() - 1)
  {
    numValues /= 1024;
    memUnitsIdx++;
  }
  std::snprintf(cstring, sizeof(cstring), "%zu %s", numValues, memUnitsList[memUnitsIdx]);

  return cstring;
}

static size_t
get_num_input_bits(int datatype)
{
  // clang-format off
  if      (datatype == CDI_DATATYPE_INT8  ) return 8;
  else if (datatype == CDI_DATATYPE_UINT8 ) return 8;
  else if (datatype == CDI_DATATYPE_INT16 ) return 16;
  else if (datatype == CDI_DATATYPE_UINT16) return 16;
  else if (datatype == CDI_DATATYPE_INT32 ) return 32;
  else if (datatype == CDI_DATATYPE_UINT32) return 32;
  else if (datatype == CDI_DATATYPE_FLT32 ) return 32;
  else if (datatype == CDI_DATATYPE_FLT64 ) return 64;
  else if (datatype == CDI_DATATYPE_PACK8 ) return 8;
  else if (datatype == CDI_DATATYPE_PACK16) return 16;
  else if (datatype == CDI_DATATYPE_PACK32) return 24;
  else if (datatype == CDI_DATATYPE_PACK  ) return 8;        // unknown
  else if (datatype > 0 && datatype <= 32 ) return datatype; // Number of packed bits in GRIB
  else                                      return 64;
  // clang-format on
}

static size_t
get_num_output_bits(int datatype)
{
  if (CdoDefault::DataType == CDI_UNDEFID)
  {
    if (CdoDefault::FileType != CDI_UNDEFID) {}
  }
  else { datatype = CdoDefault::DataType; }

  return get_num_input_bits(datatype);
}

static void
limit_string_length(char *string, size_t maxlen)
{
  string[maxlen - 1] = 0;
  auto len = std::strlen(string);
  if (len > 10)
  {
    for (size_t i = 3; i < len; ++i)
      if (string[i] == ' ' || string[i] == ',' || (i > 10 && string[i] == '.'))
      {
        string[i] = 0;
        break;
      }
  }
}

static void
print_vars_info(int operfunc, bool ensembleInfo, VarList const &varList, int vlistID, bool xsInfo)
{
  char tmpname[CDI_MAX_NAME];
  char paramstr[32];

  auto numVars = varList.numVars();
  auto nsubtypes = vlistNsubtypes(vlistID);

  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var = varList.vars[varID];

    auto tabnum = tableInqNum(vlistInqVarTable(vlistID, varID));

    std::fprintf(stdout, "%6d", varID + 1);
    std::fprintf(stdout, " : ");

    set_text_color(stdout, BLUE);
    // institute info
    auto instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
    std::strcpy(tmpname, "unknown");
    if (instptr) std::strncpy(tmpname, instptr, CDI_MAX_NAME - 1);
    limit_string_length(tmpname, 32);
    std::fprintf(stdout, "%-8s ", tmpname);

    // source info
    auto modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
    std::strcpy(tmpname, "unknown");
    if (modelptr) std::strncpy(tmpname, modelptr, CDI_MAX_NAME - 1);
    limit_string_length(tmpname, 32);
    std::fprintf(stdout, "%-8s ", tmpname);

    // timetype
    std::fprintf(stdout, "%c ", var.isConstant ? 'c' : 'v');

    // tsteptype
    std::fprintf(stdout, "%-8s ", cdo::get_steptype_name(var.stepType));

    // ensemble information
    if (ensembleInfo)
    {
      int perturbationNumber, numberOfForecastsInEnsemble;
      auto r1 = cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber);
      auto r2 = cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble);
      if (r1 == 0 && r2 == 0)
        std::fprintf(stdout, "%2d/%-2d ", perturbationNumber, numberOfForecastsInEnsemble);
      else
        std::fprintf(stdout, "--/-- ");
    }

    if (nsubtypes > 1)
    {
      auto subtypeID = vlistInqVarSubtype(vlistID, varID);
      auto subtypesize = subtypeInqSize(subtypeID);
      std::fprintf(stdout, " %6d  ", subtypesize);
      std::fprintf(stdout, "%3d ", vlistSubtypeIndex(vlistID, subtypeID) + 1);
    }
    reset_text_color(stdout);

    // layer info
    set_text_color(stdout, GREEN);
    std::fprintf(stdout, "%6d ", var.nlevels);
    reset_text_color(stdout);
    std::fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, var.zaxisID) + 1);

    // grid info
    set_text_color(stdout, GREEN);
    std::fprintf(stdout, "%9zu ", var.gridsize);
    reset_text_color(stdout);
    std::fprintf(stdout, "%3d ", vlistGridIndex(vlistID, var.gridID) + 1);

    // datatype
    set_text_color(stdout, BLUE);
    std::fprintf(stdout, " %-3s", cdo::datatype_to_cstr(var.dataType));

    auto compType = vlistInqVarCompType(vlistID, varID);
    auto isCompressed = (compType != CDI_COMPRESS_NONE);
    std::fprintf(stdout, "%c ", isCompressed ? (int) comptype_to_name(compType)[0] : ' ');

    // memType
    if (xsInfo) std::fprintf(stdout, "   %-3s", memtype_to_cstr(var.memType));

    reset_text_color(stdout);

    std::fprintf(stdout, ": ");

    // parameter info
    cdiParamToString(var.param, paramstr, sizeof(paramstr));

    // set_text_color(stdout, GREEN);
    // clang-format off
    if      (operfunc == func_name) std::fprintf(stdout, "%-14s", var.name.c_str());
    else if (operfunc == func_code) std::fprintf(stdout, " %4d %4d", tabnum, var.code);
    else                            std::fprintf(stdout, "%-14s", paramstr);
    // clang-format on
    if (xsInfo && Options::cdoVerbose && operfunc == func_name && var.units.size()) std::fprintf(stdout, " [%s]", var.units.c_str());
    // reset_text_color(stdout);

    if (Options::cdoVerbose)
    {
      auto chunkSpecString = cdo::get_chunkspec_string(vlistID, varID);
      std::fprintf(stdout, " : %s", chunkSpecString.c_str());
    }

    std::fprintf(stdout, "\n");
  }
}

static void
print_time_info(int ntsteps, int taxisID)
{
  set_text_color(stdout, BRIGHT);
  std::fprintf(stdout, "   Time coordinate");
  reset_text_color(stdout);
  std::fprintf(stdout, " :\n");

  auto taxisName = taxisNamePtr(taxisID);
  auto tname = taxisName ? taxisName : "time";
  std::fprintf(stdout, "%33s : ", tname);

  set_text_color(stdout, GREEN);
  if (ntsteps == CDI_UNDEFID)
    std::fprintf(stdout, "unlimited steps\n");
  else
    std::fprintf(stdout, "%d step%s\n", ntsteps, (ntsteps == 1) ? "" : "s");
  reset_text_color(stdout);

  if (Options::cdoVerbose)
  {
    int datatype;
    cdiInqKeyInt(taxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
    std::fprintf(stdout, "%33s : %s\n", "datatype", cdo::datatype_to_cstr(datatype));
    std::fprintf(stdout, "%33s : %d\n", "taxisID", taxisID);
  }

  if (taxisID != CDI_UNDEFID)
  {
    if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
    {
      auto rDateTime = taxisInqRdatetime(taxisID);
      std::fprintf(stdout, "     RefTime = %s %s", date_to_string(rDateTime.date).c_str(), time_to_string(rDateTime.time).c_str());

      auto tunits = taxisInqTunit(taxisID);
      if (tunits != CDI_UNDEFID) std::fprintf(stdout, "  Units = %s", tunit_to_cstr(tunits));

      auto calendar = taxisInqCalendar(taxisID);
      if (calendar != CDI_UNDEFID) std::fprintf(stdout, "  Calendar = %s", calendar_to_cstr(calendar));

      if (taxisHasBounds(taxisID)) std::fprintf(stdout, "  Bounds = true");

      std::fprintf(stdout, "\n");

      if (taxisInqType(taxisID) == TAXIS_FORECAST)
      {
        auto fDateTime = taxisInqFdatetime(taxisID);
        std::fprintf(stdout, "     ForecastRefTime = %s\n", datetime_to_string(fDateTime).c_str());
      }
    }
  }
}

static int
print_time_info_xs(int ntsteps, int taxisID, CdoStreamID streamID)
{
  int numTimesteps = ntsteps;
  TimeIncrement timeIncrement, timeIncrement0;
  CdiDateTime vDateTimeFirst{};
  CdiDateTime vDateTimeLast{};

  set_text_color(stdout, BRIGHT);
  std::fprintf(stdout, "   Time coordinate");
  reset_text_color(stdout);
  std::fprintf(stdout, " :\n");

  if (taxisID != CDI_UNDEFID)
  {
    auto calendar = taxisInqCalendar(taxisID);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID, tsID);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID);

      if (tsID)
      {
        auto julianDate0 = julianDate_encode(calendar, vDateTimeLast);
        auto julianDate = julianDate_encode(calendar, vDateTime);
        auto jdelta = julianDate_to_seconds(julianDate_sub(julianDate, julianDate0));
        timeIncrement = get_time_increment(jdelta, vDateTimeLast.date, vDateTime.date);
      }
      else { vDateTimeFirst = vDateTime; }

      if (tsID == 1) { timeIncrement0 = timeIncrement; }
      else if (tsID > 1 && timeIncrement0 != timeIncrement) { timeIncrement0.period = 0; }

      vDateTimeLast = vDateTime;

      tsID++;
    }

    numTimesteps = tsID;
  }

  auto taxisName = taxisNamePtr(taxisID);
  auto tname = taxisName ? taxisName : "time";

  std::fprintf(stdout, "%33s : ", "steps");

  set_text_color(stdout, GREEN);
  if (ntsteps == CDI_UNDEFID)
    std::fprintf(stdout, "unlimited (%d currently)\n", numTimesteps);
  else
    std::fprintf(stdout, "%d\n", ntsteps);
  reset_text_color(stdout);

  if (taxisID != CDI_UNDEFID)
  {
    std::fprintf(stdout, "%33s : ", tname);

    std::fprintf(stdout, "%s", datetime_to_string(vDateTimeFirst).c_str());
    if (numTimesteps > 1)
    {
      std::fprintf(stdout, " to %s", datetime_to_string(vDateTimeLast).c_str());
      if (timeIncrement0.period)
        std::fprintf(stdout, " by %d %s%s", (int) timeIncrement0.period, time_units_cstr(timeIncrement0.units),
                (timeIncrement0.period != 1) ? "s" : "");
    }
    std::fprintf(stdout, "\n");

    auto calendar = taxisInqCalendar(taxisID);

    if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
    {
      std::fprintf(stdout, "%33s : ", "units");

      auto tunits = taxisInqTunit(taxisID);
      std::fprintf(stdout, "%s", tunit_to_cstr(tunits));

      auto rDateTime = taxisInqRdatetime(taxisID);
      std::fprintf(stdout, " since %s\n", datetime_to_string(rDateTime).c_str());

      if (calendar != CDI_UNDEFID) std::fprintf(stdout, "%33s : %s\n", "calendar", calendar_to_cstr(calendar));

      if (taxisInqType(taxisID) == TAXIS_FORECAST)
      {
        auto fDateTime = taxisInqFdatetime(taxisID);
        std::fprintf(stdout, "%33s : %s\n", "forecastRefTime", datetime_to_string(fDateTime).c_str());
      }
    }

    if (taxisHasBounds(taxisID)) std::fprintf(stdout, "%33s : %s\n", "available", "bounds");

    if (Options::cdoVerbose)
    {
      int datatype;
      cdiInqKeyInt(taxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
      std::fprintf(stdout, "%33s : %s\n", "datatype", cdo::datatype_to_cstr(datatype));
      std::fprintf(stdout, "%33s : %d\n", "taxisID", taxisID);
    }
  }

  return numTimesteps;
}

class Sinfo : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Sinfo",
    .operators = { { "sinfo", func_generic, 0, SinfoHelp },
                   { "sinfop", func_param, 0, SinfoHelp },
                   { "sinfon", func_name, 0, SinfoHelp },
                   { "sinfoc", func_code, 0, SinfoHelp },
                   { "seinfo", func_generic, 1, SinfoHelp },
                   { "seinfop", func_param, 1, SinfoHelp },
                   { "seinfon", func_name, 1, SinfoHelp },
                   { "seinfoc", func_code, 1, SinfoHelp },
                   { "xsinfo", func_name, 2, XsinfoHelp },
                   { "xsinfop", func_param, 2, XsinfoHelp },
                   { "xsinfon", func_name, 2, XsinfoHelp },
                   { "xsinfoc", func_code, 2, XsinfoHelp } },
    .aliases = { { "infov", "infon" }, { "sinfov", "sinfon" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { -1, 0, NoRestriction },
  };
  inline static RegisterEntry<Sinfo> registration = RegisterEntry<Sinfo>();

private:
  int operfunc{ 0 };
  int ensembleInfo{ 0 };
  bool xsInfo{ false };

public:
  void
  init() override
  {
    cdiDefGlobal("COPY_CHUNKSPEC", true);

    auto operatorID = cdo_operator_id();

    operfunc = cdo_operator_f1(operatorID);
    ensembleInfo = (cdo_operator_f2(operatorID) == 1);
    xsInfo = (cdo_operator_f2(operatorID) == 2);

    operator_check_argc(0);
  }

  void
  run() override
  {
    for (int fileIdx = 0; fileIdx < cdo_stream_cnt(); fileIdx++)
    {
      auto streamID = cdo_open_read(fileIdx);
      auto vlistID = cdo_stream_inq_vlist(streamID);

      VarList varList(vlistID);

      auto nsubtypes = vlistNsubtypes(vlistID);

      set_text_color(stdout, BRIGHT);
      std::fprintf(stdout, "   File format");
      reset_text_color(stdout);
      std::fprintf(stdout, " : ");
      print_filetype(streamID, vlistID);

      set_text_color(stdout, BRIGHT);
      std::fprintf(stdout, "%6d : Institut Source   T Steptype", -(fileIdx + 1));
      if (ensembleInfo) std::fprintf(stdout, " Einfo");
      if (nsubtypes > 1) std::fprintf(stdout, " Subtypes");
      std::fprintf(stdout, " Levels Num    Points Num Dtype");
      if (xsInfo) std::fprintf(stdout, " Mtype");
      std::fprintf(stdout, " : %s",
              (operfunc == func_name) ? "Parameter name" : ((operfunc == func_code) ? "Table Code" : "Parameter ID"));

      if (Options::cdoVerbose) std::fprintf(stdout, " : Chunkspec");
      reset_text_color(stdout);
      std::fprintf(stdout, "\n");

      auto numVars = varList.numVars();
      int numFieldsConst = 0;
      int numFieldsVar = 0;
      size_t numValues = 0;
      // size_t memorySize = 0;
      size_t inputSize = 0;
      size_t outputSize = 0;
      if (Options::test && xsInfo)
      {
        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList.vars[varID];

          if (var.isConstant)
            numFieldsConst += var.nlevels;
          else
            numFieldsVar += var.nlevels;

          auto size = var.nlevels * var.gridsize * var.nwpv;
          numValues += size;

          // auto numBytes = (var.memType == MemType::Double) ? 8 : 4;
          // memorySize += size * numBytes;

          auto numBits = get_num_input_bits(var.dataType);
          inputSize += (size * numBits) / 8;

          numBits = get_num_output_bits(var.dataType);
          outputSize += (size * numBits) / 8;
        }
      }

      print_vars_info(operfunc, ensembleInfo, varList, vlistID, xsInfo);

      set_text_color(stdout, BRIGHT);
      std::fprintf(stdout, "   Grid coordinates");
      reset_text_color(stdout);
      std::fprintf(stdout, " :\n");

      print_grid_info(vlistID);

      set_text_color(stdout, BRIGHT);
      std::fprintf(stdout, "   Vertical coordinates");
      reset_text_color(stdout);
      std::fprintf(stdout, " :\n");

      print_zaxis_info(vlistID);

      if (nsubtypes > 1)
      {
        std::fprintf(stdout, "   Subtypes");
        std::fprintf(stdout, " :\n");

        print_subtype_info(vlistID);
      }

      auto taxisID = vlistInqTaxis(vlistID);
      auto numSteps = varList.numSteps();

      int numTimesteps = numSteps;

      if (numSteps != 0)
      {
        if (xsInfo) { numTimesteps = print_time_info_xs(numSteps, taxisID, streamID); }
        else
        {
          print_time_info(numSteps, taxisID);

          std::fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

          set_text_color(stdout, MAGENTA);
          print_timesteps(streamID, taxisID, Options::cdoVerbose);
          reset_text_color(stdout);
          std::fprintf(stdout, "\n");
        }
      }

      if (Options::test && xsInfo)
      {
        set_text_color(stdout, BRIGHT);
        std::fprintf(stdout, "   Summary");
        reset_text_color(stdout);
        std::fprintf(stdout, " :\n");

        int numFields = numFieldsVar * numTimesteps + numFieldsConst;
        std::fprintf(stdout, "%33s : %d\n", "number of fields", numFields);
        std::fprintf(stdout, "%33s : %d\n", "number of variables", numVars);
        if (numTimesteps) std::fprintf(stdout, "%33s : %d\n", "number of timesteps", numTimesteps);
        if (numTimesteps == 0) numTimesteps = 1;
        std::fprintf(stdout, "%33s : %zu\n", "number of values", numTimesteps * numValues);
        // std::fprintf(stdout, "%33s : %s\n", "required memory", num_values_to_byte_cstr(memorySize));
        std::fprintf(stdout, "%33s : ~%s\n", "input size", num_values_to_byte_cstr(numTimesteps * inputSize));
        std::fprintf(stdout, "%33s : ~%s\n", "output size", num_values_to_byte_cstr(numTimesteps * outputSize));
      }

      cdo_stream_close(streamID);
    }
  }

  void
  close() override
  {
  }
};
