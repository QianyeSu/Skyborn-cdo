/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include "util_wildcards.h"
#include "util_string.h"
#include "cdi_uuid.h"

static void
printf_attname(const char *varName, const char *attname)
{
  (varName) ? std::fprintf(stdout, "%s@%s = ", varName, attname) : std::fprintf(stdout, "%s = ", attname);
}

static void
print_attr_txt(const char *varName, int vlistID, int varOrGlobal, const char *attname, int attlen)
{
  printf_attname(varName, attname);
  std::vector<char> atttxt(attlen + 1);
  cdiInqAttTxt(vlistID, varOrGlobal, attname, attlen, atttxt.data());
  atttxt[attlen] = 0;
  std::fprintf(stdout, "\"");
  for (int i = 0; i < attlen; ++i)
  {
    if (atttxt[i] == '\n')
    {
      printf("\\n");
      /*
      if (atttxt[i + 1] != 0)
      {
        printf("\"\n");
        printf("             \"");
      }
      */
    }
    else if (atttxt[i] == '"') { printf("\\\""); }
    else { printf("%c", atttxt[i]); }
  }
  printf("\"\n");
}

void
print_attr_int(const char *varName, int vlistID, int varOrGlobal, const char *attname, int attlen)
{
  printf_attname(varName, attname);
  std::vector<int> attint(attlen);
  cdiInqAttInt(vlistID, varOrGlobal, attname, attlen, attint.data());
  for (int i = 0; i < attlen; ++i)
  {
    if (i) printf(", ");
    printf("%d", attint[i]);
  }
  printf("\n");
}

static void
print_attr_flt(const char *varName, int vlistID, int varOrGlobal, const char *attname, int attlen, int atttype)
{
  printf_attname(varName, attname);
  char fltstr[128];
  std::vector<double> attflt(attlen);
  cdiInqAttFlt(vlistID, varOrGlobal, attname, attlen, attflt.data());
  for (int i = 0; i < attlen; ++i)
  {
    if (i) printf(", ");
    if (atttype == CDI_DATATYPE_FLT32)
      printf("%sf", double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
    else
      printf("%s", double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
  }
  printf("\n");
}

static void
print_attr_special_global(int vlistID, const char *argument)
{
  auto gridID = vlistInqVarGrid(vlistID, 0);
  if (gridInqType(gridID) == GRID_UNSTRUCTURED)
  {
    if (auto attname = "number_of_grid_used"; argument == nullptr || (argument && wildcard_match(attname, argument)))
    {
      int number = 0;
      cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
      if (number > 0) std::fprintf(stdout, "%s = %d\n", attname, number);
    }
    if (auto attname = "grid_file_uri"; argument == nullptr || (argument && wildcard_match(attname, argument)))
    {
      if (int length = 0; CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
      {
        char referenceLink[8192];
        cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, referenceLink, &length);
        std::fprintf(stdout, "%s = \"%s\"\n", attname, referenceLink);
      }
    }
    if (auto attname = "uuidOfHGrid"; argument == nullptr || (argument && wildcard_match(attname, argument)))
    {
      unsigned char uuid[CDI_UUID_SIZE] = { 0 };
      int length = CDI_UUID_SIZE;
      auto status = cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
      if (status == CDI_NOERR && !cdiUUIDIsNull(uuid))
      {
        char uuidStr[uuidNumHexChars + 1] = { 0 };
        if (cdiUUID2Str(uuid, uuidStr) == uuidNumHexChars) std::fprintf(stdout, "%s = \"%s\"\n", attname, uuidStr);
      }
    }
  }
}

static void
print_attr_special(const char *varName, const CdoVars &cdoVars, int vlistID, int varOrGlobal, const char *argument)
{
  auto const &var = cdoVars[varOrGlobal];
  auto stdname = cdo::inq_key_string(vlistID, varOrGlobal, CDI_KEY_STDNAME);

  double addoffset = 0.0, scalefactor = 1.0;
  auto haveAddoffset = (cdiInqKeyFloat(vlistID, varOrGlobal, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
  auto haveScalefactor = (cdiInqKeyFloat(vlistID, varOrGlobal, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);

  if (argument)
  {
    if (stdname.size() && wildcard_match("standard_name", argument))
      std::fprintf(stdout, "%s@standard_name = \"%s\"\n", varName, stdname.c_str());
    if (var.longname.size() && wildcard_match("long_name", argument))
      std::fprintf(stdout, "%s@long_name = \"%s\"\n", varName, var.longname.c_str());
    if (var.units.size() && wildcard_match("units", argument)) std::fprintf(stdout, "%s@units = \"%s\"\n", varName, var.units.c_str());
    if (wildcard_match("missing_value", argument)) std::fprintf(stdout, "%s@missing_value = %g\n", varName, var.missval);
    if (haveAddoffset && wildcard_match("add_offset", argument)) std::fprintf(stdout, "%s@add_offset = %g\n", varName, addoffset);
    if (haveScalefactor && wildcard_match("scale_factor", argument))
      std::fprintf(stdout, "%s@scale_factor = %g\n", varName, scalefactor);
  }
  else
  {
    if (stdname.size()) std::fprintf(stdout, "%s@standard_name = \"%s\"\n", varName, stdname.c_str());
    if (var.longname.size()) std::fprintf(stdout, "%s@long_name = \"%s\"\n", varName, var.longname.c_str());
    if (var.units.size()) std::fprintf(stdout, "%s@units = \"%s\"\n", varName, var.units.c_str());
    std::fprintf(stdout, "%s@missing_value = %g\n", varName, var.missval);
    if (haveAddoffset) std::fprintf(stdout, "%s@add_offset = %g\n", varName, addoffset);
    if (haveScalefactor) std::fprintf(stdout, "%s@scale_factor = %g\n", varName, scalefactor);
  }
}

static void
print_attributes(const char *varName, const CdoVars &cdoVars, int vlistID, int varOrGlobal, int natts, const char *argument)
{
  if (varOrGlobal != CDI_GLOBAL) print_attr_special(varName, cdoVars, vlistID, varOrGlobal, argument);

  for (int ia = 0; ia < natts; ia++)
  {
    char attname[CDI_MAX_NAME];
    int atttype, attlen;
    cdiInqAtt(vlistID, varOrGlobal, ia, attname, &atttype, &attlen);

    if (argument && !wildcard_match(attname, argument)) continue;

    if (atttype == CDI_DATATYPE_TXT)
      print_attr_txt(varName, vlistID, varOrGlobal, attname, attlen);
    else if (atttype == CDI_DATATYPE_INT32)
      print_attr_int(varName, vlistID, varOrGlobal, attname, attlen);
    else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
      print_attr_flt(varName, vlistID, varOrGlobal, attname, attlen, atttype);
    else
      cdo_warning("Unsupported type %i name %s", atttype, attname);
  }

  if (varOrGlobal == CDI_GLOBAL) print_attr_special_global(vlistID, argument);
}

static void
check_varname_and_print(VarList const &varList, int vlistID, char *checkvarname, char *attname)
{
  auto lfound = false;
  auto numVars = varList.numVars();
  for (int varID = 0; varID < numVars; ++varID)
  {
    auto const &var = varList.vars[varID];
    if (!checkvarname || wildcard_match(var.name, checkvarname))
    {
      lfound = true;
      // std::fprintf(stdout, "%s:\n", var.name.c_str());
      std::fprintf(stdout, "\n");
      int natts;
      cdiInqNatts(vlistID, varID, &natts);
      print_attributes(var.name.c_str(), varList.vars, vlistID, varID, natts, attname);
      if (!checkvarname) break;
    }
  }
  if (!lfound && checkvarname) cdo_abort("Could not find variable %s!", checkvarname);
}

class Showattribute : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Showattribute",
    // clang-format off
    .operators = { { "showattribute", ShowattributeHelp },
                   { "showattsvar", ShowattributeHelp } },
    // clang-format on
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Showattribute> registration = RegisterEntry<Showattribute>();

private:
  int SHOWATTRIBUTE{}, SHOWATTSVAR{};

  CdoStreamID streamID{};
  int operatorID{};
  int vlistID{};

  VarList varList{};

public:
  void
  init() override
  {
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    SHOWATTRIBUTE = module.get_id("showattribute");
    SHOWATTSVAR = module.get_id("showattsvar");

    operatorID = cdo_operator_id();

    streamID = cdo_open_read(0);
    vlistID = cdo_stream_inq_vlist(streamID);

    varList = VarList(vlistID);
  }

  void
  run() override
  {
    auto numVars = varList.numVars();

    auto numArgs = cdo_operator_argc();
    if (numArgs == 0)
    {
      if (operatorID == SHOWATTSVAR) { check_varname_and_print(varList, vlistID, nullptr, nullptr); }
      else
      {
        for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var = varList.vars[varID];
          // std::fprintf(stdout, "%s:\n", var.name.c_str());
          std::fprintf(stdout, "\n");

          int nattsvar;
          cdiInqNatts(vlistID, varID, &nattsvar);
          print_attributes(var.name.c_str(), varList.vars, vlistID, varID, nattsvar, nullptr);
        }

        int natts;
        cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
        // if (natts) std::fprintf(stdout, "Global:\n");
        std::fprintf(stdout, "\n");
        print_attributes(nullptr, varList.vars, vlistID, CDI_GLOBAL, natts, nullptr);
      }
    }
    else
    {
      constexpr int delim = '@';
      auto const &argList = cdo_get_oper_argv();
      char buffer[CDI_MAX_NAME];
      for (int i = 0; i < numArgs; ++i)
      {
        std::strcpy(buffer, argList[i].c_str());
        char *result = strrchr(buffer, delim);
        char *input = buffer;
        if (result == nullptr)
        {
          if (operatorID == SHOWATTRIBUTE)
          {
            int natts;
            cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
            // if (natts) std::fprintf(stdout, "Global:\n");
            std::fprintf(stdout, "\n");
            print_attributes(nullptr, varList.vars, vlistID, CDI_GLOBAL, natts, input);
          }
          else if (operatorID == SHOWATTSVAR) { check_varname_and_print(varList, vlistID, input, nullptr); }
        }
        else
        {
          if (operatorID == SHOWATTRIBUTE)
          {
            input = result + 1;
            if (*input == 0) input = nullptr;
            *result = 0;
            char *varname = buffer;
            if (*varname == 0) cdo_abort("Variable name not specified!");
            check_varname_and_print(varList, vlistID, varname, input);
          }
          else if (operatorID == SHOWATTSVAR) { check_varname_and_print(varList, vlistID, input, nullptr); }
        }
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID);
  }
};
