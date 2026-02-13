/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>
#include "julian_date.h"
#include <cstdlib>
#include <signal.h>
#include <iomanip>

#include "c_wrapper.h"
#include "process_int.h"
#include "param_conversion.h"
#include "util_files.h"
#include "cdi_lockedIO.h"
#include "cdo_options.h"
#include "cdo_cdi_wrapper.h"
#include "varray.h"

#ifndef HAVE_LIBCMOR  //------------------------------------------------------------
class CMOR : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CMOR",
    .operators = { { "cmor", CmorHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<CMOR> registration = RegisterEntry<CMOR>();

  // clang-format off
  void
  init() override
  {
    cdo_abort("CMOR support not compiled in!");
  }
  void
  run() override
  {
    cdo_abort("CMOR support not compiled in!");
  }
  void
  close() override
  {
    cdo_abort("CMOR support not compiled in!");
  }
  // clang-format on
};

#else  //---------------------------------------------------------------------------

#include <cassert>
#include <unistd.h>
#include <netcdf.h>

extern "C"
{
#include "cmor.h"
}

#include "util_string.h"
#include "pmlist.h"
#include "datetime.h"
#include "merge_axis.h"
#include "listbuffer.h"
#include "cdo_cmor.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

static void
get_stringcode(int vlistID, int varID, std::string &varcodestring)
{
  int varcode = vlistInqVarCode(vlistID, varID);
  std::stringstream ss;
  ss << std::setw(3) << std::setfill('0') << varcode;
  varcodestring = ss.str();
}

static int
get_ifilevalue_code(int vlistID, std::string const &value, int nvars)
{
  int code = std::stoi(value);
  if (code > 0 && code < 1000)
  {
    char newcode[4];
    std::snprintf(newcode, sizeof(newcode), "%03d", code);
    for (int varID = 0; varID < nvars; ++varID)
    {
      std::string codestring;
      get_stringcode(vlistID, varID, codestring);
      if (codestring == newcode) return varID;
    }
    return CDI_UNDEFID;
  }
  else
  {
    cdo_print("'%s' could not be transformed into the code format (three digit integer). Code wont be used.", value);
    return CDI_UNDEFID;
  }
}

static int
get_ifilevalue_name(int vlistID, std::string const &value, int nvars)
{
  char ifilevalue[CDI_MAX_NAME];
  for (int varID = 0; varID < nvars; ++varID)
  {
    vlistInqVarName(vlistID, varID, ifilevalue);
    if (value == ifilevalue) return varID;
  }
  return CDI_UNDEFID;
}

static int
getVarIDToMap(int vlistID, int nvars, std::string const &key, std::string const &value)
{
  if (key.compare("code") == 0)
    return get_ifilevalue_code(vlistID, value, nvars);
  else
    return get_ifilevalue_name(vlistID, value, nvars);
}

#if (CMOR_VERSION_MAJOR == 3)
static void
removeDataset()
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  cwd[strlen(cwd)] = '\0';
  int procID = getpid();
  std::string dataset_path = std::string(cwd) + "/dataset" + std::to_string(procID) + ".json";
  remove(dataset_path.c_str());
}
#endif

void
sigfunc(int sig)
{
  if (sig == SIGTERM)
  {
#if (CMOR_VERSION_MAJOR == 3)
    removeDataset();
#endif
    cdo_abort("ERROR (infile: '%s')! Program terminated by CMOR. A temporary ofile can outlive which needs to be deleted manually.",
              cdo_get_stream_name(0));
  }
}

static std::string
kv_get_a_val(KVList *kvl, std::string const &key, std::string const &replacer)
{
  return kvl->get_first_value(key, replacer);
}

KVList
maptab_search_miptab(PMList pmlist, std::string const &cmorname, std::string const &miptab, std::string const &key)
{
  KVList listlatest;
  if (pmlist.size() && !cmorname.empty() && !miptab.empty())
  {
    for (auto &node : pmlist)
    {
      const KeyValues *kvcn = node.search(key);
      if (kvcn && kvcn->nvalues > 0 && kvcn->values[0] == cmorname)
      {
        listlatest = node;
        const KeyValues *kvmt = node.search("project_mip_table");
        if ((kvmt && kvmt->nvalues > 0 && kvmt->values[0] == miptab) || !kvmt) break;
      }
    }
  }

  return listlatest;
}

static void
handleError(std::string filename, int errnum, std::string const &argument)
{
  if (filename.empty())
    cdo_abort("ERROR (infile: '%s')! In parsing the command line:\n          More than 150 values for a key are not supported.",
              cdo_get_stream_name(0));
  switch (errnum)
  {
    case (1):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "Unexpected blank in line:\n          '%s'\n          Check syntax.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (2):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "Unexpected separator sign ',' in a key of line:\n          '%s'\n          Check syntax.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (3):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "More than 150 values for a key are not supported.",
                cdo_get_stream_name(0), filename);
      break;
    case (4):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "No values found for a keyword in line:\n          '%s'.\n          Check syntax.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (5):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "A value for a keyword begins with ',' in line:\n          '%s'.\n          Check syntax.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (6):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "A Value for a keyword has a start quote sign but no end quote sign in line:\n          '%s'.\n"
                "          Check syntax.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (7):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "Unexpected separator sign '=' or ':' is found in values in line:\n          '%s'.",
                cdo_get_stream_name(0), filename, argument);
      break;
    case (8):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "Connected lines for one keyvalue contain more than the allowed 4096 characters.",
                cdo_get_stream_name(0), filename);
      break;
    case (9):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          "
                "A ',' is found at end of a line without information in next line.",
                cdo_get_stream_name(0), filename);
      break;
    case (10):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A ',' is found at end of file.", cdo_get_stream_name(0),
                filename);
      break;
    case (11):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A character is found after the end of a string terminated "
                "with \".",
                cdo_get_stream_name(0), filename);
      break;
    case (12):
      cdo_abort("ERROR (infile: '%s')! In parsing file '%s':\n          A string value is larger than the maximum size allowed by "
                "CMOR (1024 signs).",
                cdo_get_stream_name(0), filename);
      break;
  }
}

static void
free_array(char **tofree)
{
  int i = 0;
  while (tofree[i])
  {
    std::free(tofree[i]);
    i++;
  }
  std::free(tofree);
}

static int
copy_value(std::string const &value, std::vector<std::string> &values)
{
  if (values.size() > 150) return 3;
  values.push_back(value);
  return 0;
}

static char *
to_cmor(std::string const &convert)
{
  if (convert.empty()) return nullptr;
  return ((char *) convert.c_str());
}

static void
quote_replace(std::vector<std::string> &values, int nvalues, int i)
{
  std::string source = values[nvalues];
  source = source.substr(1);
  source.erase(i - 2);
  values[nvalues] = source;
}

static std::vector<std::string>
parse_string_to_values(std::string const &workfile, std::string const &pppline, int *nvalues, std::string const &keyword)
{
  const char *pline = pppline.c_str();
  std::vector<std::string> values;
  size_t start = 0;
  while (std::isspace(static_cast<unsigned char>(pline[start]))) { start++; }
  size_t end = pppline.length();
  while (end > start && std::isspace(static_cast<unsigned char>(pline[end - 1]))) { end--; }

  std::string trimmedLine = pppline.substr(start, end - start);

  if (trimmedLine.empty()) { handleError(workfile, 4, keyword); }

  *nvalues = 0;
  size_t i = 0;
  int errh;
  while (i < trimmedLine.length())
  {
    if (trimmedLine[i] == ',')
    {
      if (i == 0) { handleError(workfile, 5, trimmedLine); }

      errh = copy_value(trimmedLine, values);
      *nvalues = (int) values.size();
      if (errh) { handleError(workfile, errh, trimmedLine); }

      if (values.back().back() == '"' || values.back().back() == '\'') { quote_replace(values, *nvalues - 1, i); }
      else
      {
        values.back().resize(i);  // truncate at the comma
      }

      if (values.back().length() > CMOR_MAX_STRING) { handleError(workfile, 12, trimmedLine); }

      i++;
      trimmedLine = trimmedLine.substr(i);
      i = 0;
    }
    else if (trimmedLine[i] == '"')
    {
      i++;
      while (trimmedLine[i] != '"')
      {
        i++;
        if (i == trimmedLine.length()) { handleError(workfile, 6, trimmedLine); }
      }
      i++;

      if (i != trimmedLine.length() && trimmedLine[i] != ',' && !std::isspace(static_cast<unsigned char>(trimmedLine[i]))
          && trimmedLine[i] != 0 && trimmedLine[i] != '\n')
      {
        handleError(workfile, 11, trimmedLine);
      }
    }
    else if (std::isspace(static_cast<unsigned char>(trimmedLine[i]))) { break; }
    else if (trimmedLine[i] == '=')
    {
      cdo_warning("Value to parse '%s' contains '='. This will not be used as assignment.", trimmedLine);
      i++;
    }
    else
    {
      i++;
    }
  }
  errh = copy_value(trimmedLine, values);
  if (errh) { handleError(workfile, errh, trimmedLine); }

  if (values.back().back() == '"') { quote_replace(values, *nvalues - 1, i); }
  else
  {
    values.back().resize(i);  // truncate at the end
  }

  trimmedLine = trimmedLine.substr(i);

  if (values.back().length() > CMOR_MAX_STRING) { handleError(workfile, 12, trimmedLine); }

  *nvalues = (int) values.size();
  return values;
}

static void
kv_insert_vals(KVList *kvl, std::string const &key, std::string const &stringValue, bool lreplace, bool lparse)
{
  /* For internal model and institution check if string is not null */
  if (!stringValue.empty())
  {
    int nvalues = 0;

    if (lreplace) { kvl->remove(key); }
    if (lparse)
    {
      cdo_print("%s", stringValue);
      const std::vector<std::string> values
          = parse_string_to_values(kv_get_a_val(kvl, "workfile4err", ""), stringValue, &nvalues, key);
      kvl->append(key, values, nvalues);
    }
    else
    {
      std::vector<std::string> values = { stringValue };
      kvl->append(key, values, 1);
    }
  }
}

static std::vector<std::string>
kv_get_vals(KVList *kvl, std::string const &key, int *numvals)
{
  std::vector<std::string> result;
  const KeyValues *kv = kvl->search(key);
  if (kv)
  {
    result = kv->values;
    *numvals = kv->nvalues;
  }
  cdo_print("%s %d", key, result.size());
  return result;
}

PMList
cdo_parse_cmor_file(std::string const &filename, bool lismap)
{
  PMList pml;

  auto fobj = c_fopen(filename, "r");
  if (fobj.get() == nullptr)
  {
    if (lismap)
      std::fprintf(stderr, "In reading the mapping table:\n          Open failed on %s: %s.", filename.c_str(), strerror(errno));
    else
      std::fprintf(stderr, "In reading info table:\n          Open failed on %s: %s.", filename.c_str(), strerror(errno));
    return pml;
  }

  ListBuffer listBuffer;
  auto status = listBuffer.read(fobj.get(), filename.c_str());
  if (status && lismap)
    cdo_abort("Read error on mapping table %s!", filename);
  else if (status)
    cdo_abort("Read error on info table %s!", filename);

  NamelistParser p;
  status = parse_list_buffer(p, listBuffer);
  if (status && lismap)
    cdo_abort("Mapping table not found!");
  else if (status)
    cdo_abort("Info table not found!");

  parse_namelist(pml, p, listBuffer.buffer.data(), true);

  return pml;
}

static const std::string
check_short_key(std::string const &key, bool l_short)
{
  std::vector<std::string> short_keys
      = { "cn", "n",  "c",  "u",   "cm",     "vc",         "p",  "i",  "ca", "za", "gi", "rtu", "mt", "om", "ms", "dr",
          "d",  "lc", "dj", "kaa", "member", "project_id", "vd", "dl", "ta", "ci", "di", "tp",  "sc", "ml", "" };
  std::vector<std::string> long_keys = { "cmor_name",
                                         "name",
                                         "code",
                                         "units",
                                         "cell_methods",
                                         "variable_comment",
                                         "positive",
                                         "info",
                                         "character_axis",
                                         "z_axis",
                                         "grid_info",
                                         "required_time_units",
                                         "mapping_table",
                                         "output_mode",
                                         "max_size",
                                         "drs_root",
                                         "drs",
                                         "last_chunk",
                                         "dataset_json",
                                         "keep_all_attributes",
                                         "member",
                                         "project_id",
                                         "version_date",
                                         "deflate_level",
                                         "t_axis",
                                         "climatology_interval",
                                         "decadal_interval",
                                         "tracking_prefix",
                                         "save_chunk",
                                         "move_longitudes",
                                         "" };

  for (size_t i = 0; i < short_keys.size(); ++i)
  {
    if (key == short_keys[i] || key == long_keys[i]) { return l_short ? short_keys[i] : long_keys[i]; }
  }

  return "";
}

static void
map_it(KVList *kvl, int vlistID, int varID, std::string const &var2map)
{
  for (auto const &kv : *kvl)
  {
    const std::string longkey = !check_short_key(kv.key, false).empty() ? check_short_key(kv.key, false) : "";
    if (longkey.empty())
    {
      if (kv.key != "project_mip_table")
        cdo_warning("In variable mapping:\n           you try to assign '%s' to variable '%s'.\n          This "
                    "mapping table keyword is skipped. Check allowed mapping table keywords.",
                    kv.key.c_str(), var2map);
      continue;
    }
    const std::string value = kv.values[0];
    if (value.empty()) continue;

    bool hasValidMin = false, hasValidMax = false;
    CmorVar kcv;
    kcv.name = "cdocmor";
    KeyValues kvlong;
    kvlong.key = longkey;
    kvlong.nvalues = kv.nvalues;
    kvlong.values.resize(kvlong.nvalues);
    for (int i = 0; i < kvlong.nvalues; ++i) kvlong.values[i] = kv.values[i];
    mapvar(vlistID, varID, kvlong, kcv, hasValidMin, hasValidMax, 0, false);
  }
}

static int
change_name_via_name(int vlistID, std::string const &map_name, std::string const &cmor_name)
{
  char name[CDI_MAX_NAME];
  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
  {
    vlistInqVarName(vlistID, varID, name);
    if (map_name == name)
    {
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word(cmor_name.c_str()));
      return 1;
    }
  }
  return 0;
}

static int
change_name_via_code(int vlistID, std::string const &map_code, std::string const &cmor_name)
{
  int codeproof = std::stol(map_code);
  if (!codeproof || codeproof > 1000)
  {
    cdo_print("'%s' could not be transformed into the code format (three digit integer). Code wont be used.", map_code);
    return 0;
  }

  int code;
  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
  {
    code = vlistInqVarCode(vlistID, varID);
    if (codeproof == code)
    {
      cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, parameter_to_word((const char *) cmor_name.c_str()));
      return 1;
    }
  }
  return 0;
}

struct mapping
{
  int help_var;
  int cdi_varID;
  int cmor_varID;
  int zfactor_id;
  int charvars;
  char datatype;
  void *data;
};

static struct mapping *
construct_var_mapping(int vlistID)
{
  int nvars_max = vlistNvars(vlistID);
  struct mapping *vars = (struct mapping *) std::malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  vars[0].cmor_varID = CMOR_UNDEFID;
  vars[0].data = nullptr;
  vars[0].charvars = 0;
  return vars;
}

static void
destruct_var_mapping(struct mapping vars[])
{
  for (int i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i) std::free(vars[i].data);
  std::free(vars);
}

static struct mapping *
map_var(int cdi_varID, struct mapping vars[])
{
  for (int i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    if (cdi_varID == vars[i].cdi_varID) return &vars[i];
  return nullptr;
}

static struct mapping *
new_var_mapping(struct mapping vars[])
{
  int i;
  for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i);
  vars[i + 1].cdi_varID = CDI_UNDEFID;
  vars[i + 1].cmor_varID = CMOR_UNDEFID;
  vars[i + 1].data = nullptr;
  vars[i + 1].charvars = 0;
  return &vars[i];
}

static int *
new_axis_id(int *axis_ids)
{
  int i;
  for (i = 0; axis_ids[i] != CMOR_UNDEFID; ++i);
  axis_ids[i + 1] = CMOR_UNDEFID;
  return &axis_ids[i];
}

static int
count_axis_ids(int *axis_ids)
{
  int i;
  for (i = 0; axis_ids[i] != CMOR_UNDEFID; ++i);
  return i;
}

static const KVList *
check_for_charvars(KVList *cmorVarLine, std::string const &key)
{
  /***/
  /* If a infile variable selector (name or code) has more than one value, it must be a character coordinate*/
  /* If it is given as a string and the string contains a ',', */
  /* it must be divided into several values and is a variable with character coordinate */
  /***/
  if (Options::cdoVerbose) cdo_print("Start to check for variables to merge as a character coordinate.");
  const KeyValues *kvn = nullptr;
  bool isCode = false;
  if (!key.empty())
    kvn = cmorVarLine->search(key);
  else
  {
    kvn = cmorVarLine->search("name");
    if (!kvn)
    {
      kvn = cmorVarLine->search("code");
      if (kvn) isCode = true;
    }
  }
  if (kvn && kvn->nvalues > 1) return cmorVarLine;

  if (kvn && kvn->values.size() == 1 && kvn->values[0].find(',') != std::string::npos)
  {
    if (Options::cdoVerbose) cdo_print("Start to replace identifier string with its values.");
    std::string unfilteredComma = kvn->values[0];
    if (isCode)
      kv_insert_vals(cmorVarLine, "code", unfilteredComma, true, true);
    else
      kv_insert_vals(cmorVarLine, "name", unfilteredComma, true, true);
    if (Options::cdoVerbose) cdo_print("Successfully replaced identifier string with its values.");
    return cmorVarLine;
  }
  if (Options::cdoVerbose) cdo_print("Successfully checked for variables to merge as a character coordinate.");
  return nullptr;
}

static void
addcharvar(const KeyValues *charvars, int vlistID, std::string const &key, struct mapping vars[])
{
  MergeVarsOnAxis withnewcharaxis;
  withnewcharaxis.inputNames = *charvars;

  if (Options::cdoVerbose) cdo_print("Start to merge variables to one character coordinate.");
  int nvars = vlistNvars(vlistID);
  for (int i = 0; i < charvars->nvalues; ++i)
  {
    MergeVarKeys temp;
    temp.vlistID = vlistID;
    temp.varID = getVarIDToMap(vlistID, nvars, key, withnewcharaxis.inputNames.values[i]);
    if (temp.varID == CDI_UNDEFID)
      cdo_abort("ERROR (infile: '%s')! In merging variables to a variable with a character coordinate:\n          Could not find "
                "'%s' in infile '%s' to build a variable with character coordinate.",
                cdo_get_stream_name(0), withnewcharaxis.inputNames.values[i].c_str(), cdo_get_stream_name(0));
    temp.gridID = vlistInqVarGrid(vlistID, temp.varID);
    temp.projID = gridInqProj(temp.gridID);
    temp.zaxisID = vlistInqVarZaxis(vlistID, temp.varID);
    withnewcharaxis.inputKeys.push_back(temp);
  }

  std::vector<int> axissize(3);
  axissize[0] = gridInqXsize(withnewcharaxis.inputKeys[0].gridID);
  axissize[1] = gridInqYsize(withnewcharaxis.inputKeys[0].gridID);
  axissize[2] = zaxisInqSize(withnewcharaxis.inputKeys[0].zaxisID);

  int oldgridsize = gridInqSize(withnewcharaxis.inputKeys[0].gridID);
  if (axissize[0] * axissize[1] == 0) oldgridsize = 1;

  if (axissize[0] != 1 && axissize[1] != 1 && axissize[2] != 1)
  {
    if (Options::cdoVerbose)
      cdo_print("In merging variables to an axis:\n          Spatial dimensions are already set. A fourth axis is created.");
    char ids[CDI_MAX_NAME];
    std::snprintf(ids, CDI_MAX_NAME, "%d", withnewcharaxis.inputKeys[0].varID);
    for (int i = 1; i < charvars->nvalues; ++i)
    {
      char tempint[sizeof(int)];
      std::snprintf(tempint, sizeof(tempint), "%d", withnewcharaxis.inputKeys[i].varID);
      std::strcat(ids, ",");
      std::strcat(ids, tempint);
    }
    cdiDefAttTxt(vlistID, withnewcharaxis.inputKeys[0].varID, "merge_axis", (int) std::strlen(ids), ids);
    return;
  }
  int ntsteps = vlistNtsteps(vlistID);

  if (cdo_assert_files_only() == false)
    cdo_abort("ERROR (infile: '%s')! Cannot merge variables to a new axis because the input file needs to be opened twice but you "
              "piped operators so the input file is not clear.",
              cdo_get_stream_name(0));

  auto streamID2 = streamOpenRead(cdo_get_stream_name(0));
  if (ntsteps == -1)
  {
    ntsteps = 0;
    int dummy;
    while ((dummy = streamInqTimestep(streamID2, ntsteps++)));
  }

  axissize = withnewcharaxis.define_new_axes(axissize);
  withnewcharaxis.define_var_structure(vlistID, ntsteps, axissize);

  struct mapping *var = new_var_mapping(vars);
  var->cdi_varID = withnewcharaxis.output.varID;
  var->datatype = withnewcharaxis.output.datatype;

  withnewcharaxis.read_cmor_charvar(axissize, streamID2, oldgridsize);

  var->data = std::malloc(ntsteps * axissize[0] * axissize[1] * axissize[2] * sizeof(double));
  for (int i = 0; i < ntsteps * axissize[0] * axissize[1] * axissize[2]; i++)
  {
    if (withnewcharaxis.output.datatype == 'd')
      ((double *) var->data)[i] = withnewcharaxis.data[i];
    else
      ((float *) var->data)[i] = (float) withnewcharaxis.data[i];
  }
  var->charvars = 1;

  streamClose(streamID2);

  if (Options::cdoVerbose)
    cdo_print("Successfully merged variables into one character axis. The final variable is called '%s' and has the ID: '%d'",
              charvars->values[0].c_str(), var->cdi_varID);
}

static int
maptab_via_key(std::string const &tablename, PMList pml, int vlistID, int varID, std::string const &key,
               std::string const &miptabfreq)
{
  char ifilevalue[CDI_MAX_NAME];
  std::string cppifile;
  if ((key.compare("cmor_name") == 0) || (key.compare("name") == 0))
    vlistInqVarName(vlistID, varID, ifilevalue);
  else
  {
    get_stringcode(vlistID, varID, cppifile);
    cppifile.copy(ifilevalue, cppifile.size());
    ifilevalue[cppifile.size()] = '\0';  // Add null terminator
  }

  if (ifilevalue[0])
  {
    KVList kvl = maptab_search_miptab(pml, ifilevalue, miptabfreq, key);
    if (kvl.size())
    {
      if (Options::cdoVerbose) cdo_print("Start to map via '%s'.", key);
      map_it(&kvl, vlistID, varID, ifilevalue);
      return 1;
    }
    if (Options::cdoVerbose)
      cdo_print("In variable mapping with table '%s':\n          "
                "Variable named '%s' with varID '%d' could not be mapped via '%s' "
                "because no corresponding key '%s' was found in mapping table file.",
                tablename, ifilevalue, varID, key, key);
    return 0;
  }
  else
  {
    if (Options::cdoVerbose)
      cdo_print("In variable mapping with table '%s':\n          "
                "Variable with varID '%d' could not be mapped via '%s' because it "
                "does not possess a '%s' in infile.",
                tablename, varID, key, key);
    return 0;
  }
}

static int
maptab_via_cn_and_key(KVList *kvl_oname, int vlistID, int nvars, std::string const &key)
{
  const KeyValues *kv = kvl_oname->search(key);
  const KeyValues *kvcn = kvl_oname->search("cmor_name");
  if (kv)
  {
    int varID = (key.compare("cmor_name") == 0) ? getVarIDToMap(vlistID, nvars, "name", kv->values[0])
                                                : getVarIDToMap(vlistID, nvars, key, kv->values[0]);
    int newvar = getVarIDToMap(vlistID, nvars, "name", kvcn->values[0]);
    if (varID != CDI_UNDEFID)
    {
      if (newvar != CDI_UNDEFID && newvar != varID)
      {
        cdo_warning("In variable mapping : \n          "
                    "You try to define a variable '%s' that is already in the input stream.\n"
                    "The already existing infile variable will be renamed.\n"
                    "This may lead to errors. Choose the respective var with '-selname' before applying cmor.",
                    kvcn->values[0]);
        cdiDefKeyString(vlistID, newvar, CDI_KEY_NAME, "overwritten");
      }
      if (Options::cdoVerbose) cdo_print("Started mapping of variable via '%s'.", key);
      map_it(kvl_oname, vlistID, varID, kv->values[0].c_str());
      return 1;
    }
    cdo_print("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped "
              "via key '%s' because no infile variable '%s' equals '%s'.",
              kvcn->values[0], key, key, kv->values[0]);
  }
  else if (Options::cdoVerbose)
    cdo_print("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped "
              "via key '%s' because it possesses no corresponding key '%s' in mapping file.",
              kvcn->values[0], key, key);
  return 0;
}

static void
maptab_via_cmd(std::string tablename, PMList pml, std::string const &origValue, int vlistID, std::string const &key,
               std::string cmorName, std::string const &miptabfreq, int filetype, std::string const &maptab)
{
  KVList cmorVarLine;

  int nvars = vlistNvars(vlistID);
  int varIDToMap = getVarIDToMap(vlistID, nvars, key, origValue);
  if (varIDToMap == CDI_UNDEFID)
    cdo_abort("ERROR (infile: '%s')! In variable mapping with table '%s':\n          "
              "Variable with '%s': '%s' configured via cmdline could not be "
              "found in infile '%s'.",
              cdo_get_stream_name(0), tablename, key, origValue, cdo_get_stream_name(0));

  cmorVarLine = maptab_search_miptab(pml, cmorName, miptabfreq, "cmor_name");
  const KeyValues *kvcn = cmorVarLine.search("cmor_name");
  int newvar = getVarIDToMap(vlistID, nvars, "name", kvcn->values[0]);
  if (newvar != CDI_UNDEFID && newvar != varIDToMap)
  {
    cdo_warning("In variable mapping : \n          "
                "You try to define a variable '%s' that is already in the input stream.\n"
                "The already existing infile variable will be renamed.\n"
                "This may lead to errors. Choose the respective var with '-selname' before applying cmor.",
                kvcn->values[0]);
    cdiDefKeyString(vlistID, newvar, CDI_KEY_NAME, "overwritten");
  }
  if (!cmorVarLine.size())
  {
    cdo_warning("In variable mapping with table '%s':\n          "
                "The registered cmor_name '%s' via cmdline could not be found in "
                "mapping table.\n          No mapping table is applied.",
                tablename, cmorName);
    cdiDefKeyString(vlistID, varIDToMap, CDI_KEY_NAME, parameter_to_word((const char *) cmorName.c_str()));
  }
  else
  {
    if ((filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2) && key.compare("name") == 0)
    {
      if (Options::cdoVerbose)
        cdo_print("5.1. In applying the mapping table:\n          Note that you use 'name' as selector keyword "
                  "allthough the type of infile is GRB.");
    }
    if (Options::cdoVerbose) cdo_print("Started mapping of variable via '%s'.", key);
    map_it(&cmorVarLine, vlistID, varIDToMap, cmorName);

    if (Options::cdoVerbose) cdo_print("5. Successfully found, read and applied mapping table '%s'.", maptab);
  }
}

static void
maptab_via_cn(std::string const &tablename, PMList pml, std::vector<std::string> request, int vlistID, int numvals,
              std::string const &miptabfreq, int filetype, struct mapping vars[], bool isWarn)
{
  for (int j = 0; j < numvals; ++j)
  {
    KVList kvl_oname = maptab_search_miptab(pml, request[j], miptabfreq, "cmor_name");
    if (kvl_oname.size())
    {
      if (vars && check_for_charvars(&kvl_oname, ""))
      {
        const KeyValues *charcode = kvl_oname.search("name");
        if (!charcode)
        {
          charcode = kvl_oname.search("code");
          addcharvar(charcode, vlistID, "code", vars);
        }
        else
          addcharvar(charcode, vlistID, "name", vars);
      }
      int nvars = vlistNvars(vlistID);
      if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
      {
        if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "code"))
        {
          if (Options::cdoVerbose) cdo_print("Successfully mapped variable via code to cmor_name '%s'.", request[j]);
          continue;
        }
      }
      if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "name"))
      {
        if (Options::cdoVerbose) cdo_print("Successfully mapped variable via name to cmor_name: '%s'.", request[j]);
        continue;
      }
      else
      {
        if (Options::cdoVerbose)
          cdo_print("In variable mapping with table '%s':\n          "
                    "Try to use cmor_name for selecting the infile variable.",
                    tablename, request[j]);
        if (maptab_via_cn_and_key(&kvl_oname, vlistID, nvars, "cmor_name"))
        {
          cdo_print("Successfully mapped variable via cmor_name to cmor_name '%s'.", request[j]);
          continue;
        }
        if (isWarn)
          cdo_warning("In variable mapping with table '%s':\n          "
                      "Mapping table line of cmor_name '%s' could neither be mapped "
                      "via 'name', 'code' nor 'cmor_name'.\n          No mapping for cmor_name: '%s'.",
                      tablename, request[j], request[j]);
        continue;
      }
    }
    else
    {
      if (isWarn)
        cdo_warning("In variable mapping with table '%s':\n          "
                    "Requested cmor_name: '%s' is not found in mapping table.\n       "
                    "   No mapping for cmor_name: '%s'",
                    tablename, request[j], request[j]);
      continue;
    }
  }
}

static bool
file_exist(std::string const &tfilename, bool force, std::string const &fileart, bool print)
{
  assert(!tfilename.empty());
  size_t filesize = FileUtils::size(tfilename);
  if (filesize == 0 && force)
    cdo_abort("ERROR (infile: '%s')!\n          Empty '%s' file: '%s'.", cdo_get_stream_name(0), fileart, tfilename);
  else if (filesize == 0 && !force)
  {
    if (print && Options::cdoVerbose) cdo_print("Empty '%s' file: '%s'.", fileart, tfilename);
    return false;
  }
  if ((tfilename.find(".nc") != std::string::npos) || (tfilename.find(".grb") != std::string::npos)) return 1;
  auto fp = std::fopen(tfilename.c_str(), "r");
  if (fp == nullptr && force)
    cdo_abort("ERROR (infile: '%s')!\n           Open failed on '%s' file: '%s'.", cdo_get_stream_name(0), fileart, tfilename);
  else if (fp == nullptr && !force)
  {
    if (print && Options::cdoVerbose) cdo_print("Open failed on '%s' file: '%s'.", fileart, tfilename);
    return false;
  }

  std::fclose(fp);
  return true;
}

static int
parse_kv_file(KVList *kvl, std::string const &filename)
{
  PMList pmkv = cdo_parse_cmor_file(filename, false);
  if (!pmkv.size()) return 0;

  for (auto &kv : pmkv.front())
  {
    std::string keyname = kv.key;
    const KeyValues *kvfromlist = kvl->search(kv.key.c_str());
    if (kvfromlist) continue;

    std::vector<std::string> values(kv.nvalues + 1);
    int k = 0;
    for (k = 0; k < kv.nvalues; k++) values[k] = kv.values[k];

    kvl->append(keyname, values, kv.nvalues);
  }
  return 1;
}

static void
check_compare_set(std::string &finalset, std::string &attribute, std::string const &attname, std::string const &defaultstr)
{
  if (finalset.empty())
  {
    if (attribute.empty())
    {
      if (!defaultstr.empty())
        finalset = defaultstr;
      else
        cdo_abort("ERROR (infile: '%s')! Function 'attribute check, compare and set':\n          Required value for attribute '%s' "
                  "is neither found in input file nor in the configuration.",
                  cdo_get_stream_name(0), attname);
    }
    else
      finalset = attribute;
  }
  else if (!attribute.empty())
  {
    if (attribute != finalset)
    {
      if (Options::cdoVerbose)
        cdo_print("Function 'attribute check, compare and set':\n          Be aware of differences between infile and "
                  "user specification.\n          Attribute '%s' in input file: '%s' does not agree with user "
                  "specification '%s'.",
                  attname, finalset, attribute);
      cdo_print("Attribute '%s' = '%s'", attname, attribute);
      finalset = attribute;
    }
  }
}

static std::string
get_infile_attvalue(int vlistID, int varID, std::string &name, int type, int len)
{
  std::string infile_attvalue = "";
  if (type == CDI_DATATYPE_INT32)
  {
    std::vector<int> values(len);
    int *c_array = values.data();
    cdiInqAttInt(vlistID, varID, name.c_str(), len, c_array);
    std::ostringstream oss;
    oss << values[0];
    for (int l = 1; l < len; ++l) { oss << values[l]; }
    infile_attvalue = oss.str();
  }
  else if (type == CDI_DATATYPE_FLT32 || type == CDI_DATATYPE_FLT64)
  {
    char fltstr[128];
    std::vector<double> attflt(len);
    cdiInqAttFlt(vlistID, varID, name.c_str(), len, attflt.data());
    for (int i = 0; i < len; ++i)
    {
      if (i > 0) { infile_attvalue += ", "; }
      /*std::strcat(infile_attvalue, i ? "," : " ");*/
      char tempflt[64];
      if (type == CDI_DATATYPE_FLT32)
      {
        std::snprintf(tempflt, sizeof(tempflt), "%sf",
                      double_to_att_str(Options::CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
      }
      else
      {
        std::snprintf(tempflt, sizeof(tempflt), "%s",
                      double_to_att_str(Options::CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
      }
      tempflt[sizeof(tempflt) - 1] = '\0';
      infile_attvalue += tempflt;
    }
  }
  else
  {
    infile_attvalue += cdo::inq_att_string(vlistID, varID, name);
  }

  return infile_attvalue;
}

static std::string
get_infile_attname(int vlistID, int varID, int natt, int *type, int *len)
{
  char infile_attname[CDI_MAX_NAME];
  cdiInqAtt(vlistID, varID, natt, infile_attname, type, len);
  return std::string(infile_attname);
}

static std::string
get_txtatt(int vlistID, int varID, std::string const &key)
{
  int natts;
  std::string returnvalue = "";
  cdiInqNatts(vlistID, varID, &natts);
  for (int i = 0; i < natts; ++i)
  {
    int type, len;
    std::string infile_attname = get_infile_attname(vlistID, varID, i, &type, &len);
    if (infile_attname == key) { returnvalue = get_infile_attvalue(vlistID, varID, infile_attname, type, len); }
  }
  return returnvalue;
}

static void
get_all_atts(KVList *kvl, int vlistID, std::vector<std::string> const &infileAttSpecLong,
             std::vector<std::string> const &infileAttSpecShort)
{
  int natts, type = 0, len = 0;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  for (int i = 0; i < natts; ++i)
  {
    std::string infile_attname = get_infile_attname(vlistID, CDI_GLOBAL, i, &type, &len);

    if (infile_attname == "history") continue;

    size_t j = 0;
    for (j = 0; j < infileAttSpecLong.size(); ++j)
    {
      if (infile_attname.compare(infileAttSpecLong[j]) == 0)
      {
        kv_insert_vals(kvl, infileAttSpecShort[j], get_infile_attvalue(vlistID, CDI_GLOBAL, infile_attname, type, len), false,
                       true);
        break;
      }
    }
    if (j == infileAttSpecLong.size())
    {
      kv_insert_vals(kvl, infile_attname, get_infile_attvalue(vlistID, CDI_GLOBAL, infile_attname, type, len), false, false);
    }
  }

  kv_insert_vals(kvl, "institution", institutInqLongnamePtr(vlistInqInstitut(vlistID)), false, false);
  /*  if ( strcmp(project_id, "CMIP6") != 0 ) */
  kv_insert_vals(kvl, "source", modelInqNamePtr(vlistInqModel(vlistID)), false, false);
}

static void
add_globalhybrids(KVList *kvl, int vlistID)
{
  /* Do not define history, source or institution as one of those: */
  std::vector<std::string> infileAttSpecLong = { "branch_dates", "climatology_interval", "decadal_interval" };
  std::vector<std::string> infileAttSpecShort = { "branch_dates", "ci", "di" };
  if (kv_get_a_val(kvl, "kaa", "n") == "y" || kv_get_a_val(kvl, "keep_all_attributes", "n") == "y")
    get_all_atts(kvl, vlistID, infileAttSpecLong, infileAttSpecShort);
  else
  {
    for (size_t i = 0; i < infileAttSpecLong.size(); ++i)
    {
      std::string att = get_txtatt(vlistID, CDI_GLOBAL, infileAttSpecLong[i]);
      if (!att.empty()) kv_insert_vals(kvl, infileAttSpecShort[i], att, false, false);
    }

    std::vector<std::string> infileAtt = { "member", "variant_label" };
    for (size_t i = 0; i < infileAtt.size(); ++i)
    {
      std::string att = get_txtatt(vlistID, CDI_GLOBAL, infileAtt[i]);
      if (!att.empty()) kv_insert_vals(kvl, infileAtt[i], att, false, false);
    }
  }
  std::vector<std::string> longAtt = {
    "required_time_units",
    "grid_info",
    "mapping_table",
    "keep_all_attributes",
    "drs_root",
    "drs",
    "t_axis",
    "version_date",
    "deflate_level",
    "character_axis",
    "z_axis",
    "output_mode",
    "max_size",
    "last_chunk",
    "save_chunk",
    "move_longitudes",
  };
  std::vector<std::string> shortAtt
      = { "rtu", "gi", "mt", "kaa", "dr", "d", "ta", "vd", "dl", "ca", "za", "om", "ms", "lc", "sc", "ml" };
  for (size_t i = 0; i < longAtt.size(); ++i)
  {
    for (auto &kv_latt : *kvl)
    {
      if (kv_latt.key == longAtt[i])
      {
        kv_insert_vals(kvl, shortAtt[i], kv_latt.values[0], false, false);
        kvl->remove(kv_latt.key);
        break;
      }
    }
  }
}

static int
check_attarray(KVList *kvl, std::vector<std::string> const &reqAtt, int vlistID)
{
  for (size_t i = 0; i < reqAtt.size(); ++i)
  {
    const KeyValues *kv_reqatt = kvl->search(reqAtt[i]);
    if (!kv_reqatt || kv_reqatt->nvalues == 0)
    {
      if (kv_get_a_val(kvl, "kaa", "y") != "y")
      {
        std::string infileatt = get_txtatt(vlistID, CDI_GLOBAL, reqAtt[i]);
        if (!infileatt.empty() && (infileatt != "notSet"))
        {
          if (Options::cdoVerbose) cdo_print("Attribute (from infile) '%s' is '%s'. ", reqAtt[i], infileatt);
          std::vector<std::string> infileatts = { infileatt };
          kvl->append(reqAtt[i], infileatts, 1);
        }
        return i;
      }
      else
        return i;
    }
    else if (kv_reqatt->values[0] == "notSet")
      return i;
    else if (Options::cdoVerbose)
      cdo_print("Attribute (from meta data file) '%s' is '%s'. ", reqAtt[i], kv_reqatt->values[0]);
  }
  return -1;
}

static void
attErr(std::vector<std::string> const &reqAtt, int errnum)
{
  std::string errStr;

  errStr = "ERROR! Attribute '" + reqAtt[errnum]
           + "' is required. Either it is missing, 'notSet', or the value is invalid.\n       "
             "   Make sure that you have configured all following attributes:\n   "
             "       "
           + reqAtt[0];
  for (size_t i = 1; i < reqAtt.size(); ++i)
  {
    errStr += ", ";
    errStr += reqAtt[i];
  }
  cdo_abort(errStr);
}

static void
check_attr(KVList *kvl, std::string const &project_id, int vlistID)
{
  /* Project id moved to main void fct */
  std::vector<std::string> reqAtt = { "institution", "source", "experiment_id", "rtu" };
  std::vector<std::string> reqAttCMIP5 = { "institute_id", "product", "member" };
  std::vector<std::string> reqAttCORDEX
      = { "institute_id", "product", "member", "CORDEX_domain", "driving_model_id", "rcm_version_id" };
  /*  std::vector<std::string>reqAttCMIP6CMOR3[] = {"outpath", "output_path_template", "output_file_template", "tracking_prefix",
  nullptr};
  "further_info_url", "variant_label",*/
  std::vector<std::string> reqAttCMIP6 = { "activity_id", "experiment",         "grid",      "grid_label", "institution_id",
                                           "license",     "nominal_resolution", "source_id", "source_type" };
  std::vector<std::string> expdepAttCMIP6
      = { "parent_experiment_id", "parent_activity_id", "parent_mip_era",    "parent_source_id", "parent_variant_label",
          "parent_time_units",    "sub_experiment",     "sub_experiment_id", "branch_method" };
  /* In all Projects needed Attributes are tested first */

  int errnum = 0;
  if ((errnum = check_attarray(kvl, reqAtt, vlistID)) != -1) attErr(reqAtt, errnum);

#if (CMOR_VERSION_MAJOR == 2)
  /* Set default attributes */
  std::vector<std::string> reqAttCMOR2 = { "contact", "model_id" };
  if ((errnum = check_attarray(kvl, reqAttCMOR2, vlistID)) != -1) attErr(reqAttCMOR2, errnum);

  const KeyValues *kv = kvl->search("references");
  if (!kv || kv->nvalues == 0 || kv->values[0] == "notSet")
  {
    const KeyValues *kv_model_id = kvl->search("model_id");
    std::string references;
    references = "No references available for ";
    std::strcat(references, kv_model_id->values[0].c_str());
    cdo_print("Attribute 'references' is set to '%s' ", references);
    kv_insert_vals(kvl, "references", references, false, false);
  }
#endif
  /* Special check for CMIP or CORDEX projects */
  if (project_id == "CMIP5")
  {
    if (Options::cdoVerbose) cdo_print("Since the project id is CMIP5 further attributes are tested. ");
    if ((errnum = check_attarray(kvl, reqAttCMIP5, vlistID)) != -1) attErr(reqAttCMIP5, errnum);
#if (CMOR_VERSION_MAJOR == 3)
    std::vector<std::string> reqAttCMIP5CMOR3 = { "modeling_realm" };
    /**************/
    /* Add additional attributes for CMIP5 */
    /* allthough using CMOR 3 */
    /**************/
    if ((errnum = check_attarray(kvl, reqAttCMIP5CMOR3, vlistID)) != -1) attErr(reqAttCMIP5CMOR3, errnum);
#endif
  }
  else if (project_id == "CORDEX")
  {
    if (Options::cdoVerbose) cdo_print("Since the project id is CORDEX further attributes are tested.");
    if ((errnum = check_attarray(kvl, reqAttCORDEX, vlistID)) != -1) attErr(reqAttCORDEX, errnum);
  }
  else if (project_id == "CMIP6")
  {
    kv_insert_vals(kvl, "_cmip6_option", "CMIP6", true, false);
    if (Options::cdoVerbose) cdo_print("Since the project_id is CMIP6 further attributes are tested.");
    if ((errnum = check_attarray(kvl, reqAttCMIP6, vlistID)) != -1) attErr(reqAttCMIP6, errnum);
    std::string pei;
    if (kv_get_a_val(kvl, "parent_experiment_id", "no parent") != "no parent")
    {
      pei = kv_get_a_val(kvl, "parent_experiment_id", "");
      cdo_print("Since you set attribute 'parent_experiment_id'='%s', further attributes are checked.", pei);
    }
    for (size_t j = 0; j < expdepAttCMIP6.size(); ++j)
    {
      const KeyValues *kv_reqatt = kvl->search(expdepAttCMIP6[j]);
      if (!kv_reqatt || kv_reqatt->nvalues == 0 || kv_reqatt->values[0] == "notSet")
      {
        if (Options::cdoVerbose || !pei.empty())
          cdo_print("Depending on the experiment, attribute '%s' may be required. Either it is missing or notSet",
                    expdepAttCMIP6[j]);
      }
      else if (Options::cdoVerbose)
      {
        cdo_print("Attribute (from meta data file) '%s' is '%s'. ", expdepAttCMIP6[j], kv_reqatt->values[0]);
      }
    }
  }
#if (CMOR_VERSION_MAJOR == 3)
  else
  {
    std::vector<std::string> customProjectCMOR3
        = { "_controlled_vocabulary_file", "_FORMULA_VAR_FILE",    "_AXIS_ENTRY_FILE",    "_history_template",
            "_further_info_url_tmpl",      "output_path_template", "output_file_template" };
    std::vector<std::string> customProjectCMOR3default = { "CMIP6_CV.json",
                                                           FORMULA_VAR_FILENAME,
                                                           AXIS_ENTRY_FILENAME,
                                                           CMOR_DEFAULT_HISTORY_TEMPLATE,
                                                           CMOR_DEFAULT_FURTHERURL_TEMPLATE,
                                                           CMOR_DEFAULT_PATH_TEMPLATE,
                                                           CMOR_DEFAULT_FILE_TEMPLATE };
    const KeyValues *kv = kvl->search("_controlled_vocabulary_file");
    if (!kv || kv->nvalues == 0 || kv->values[0] == "notSet")
    {
      if (Options::cdoVerbose)
        cdo_print("Since you have specified a project different than CMIP and you use CMOR version 3, we recommend that you "
                  "also specify the following variables: '%s' '%s' '%s'. The default values are: '%s '%s' '%s' and need to be "
                  "available in the MIP-tables directory.",
                  customProjectCMOR3[0], customProjectCMOR3[1], customProjectCMOR3[2], customProjectCMOR3default[0],
                  customProjectCMOR3default[1], customProjectCMOR3default[2]);
    }
  }
#endif
}

static void
check_mem(KVList *kvl, std::string const &project_id)
{
  /*Check if both is registered */
  std::vector<std::string> ensindexCMIP5 = { "realization", "initialization_method", "physics_version" };
  std::vector<std::string> ensindexCMIP6 = { "realization_index", "initialization_index", "physics_index", "forcing_index" };
  /*std::string ripf = "ripf"; */
  /*std::string rip = "rip"; */

  if (project_id == "CMIP5" || project_id == "CORDEX")
  {
    std::string cm = kv_get_a_val(kvl, "cm", " ");

    if (cm[0] == 'n' && ((project_id == "CMIP5") || project_id == "CORDEX"))
    {
      if (project_id == "CORDEX") kv_insert_vals(kvl, "driving_model_ensemble_member", "r0i0p0", true, false);
      kv_insert_vals(kvl, "member", "r0i0p0", true, false);
    }
    else
    {
      std::string vlabel = kv_get_a_val(kvl, "member", "");
      bool lanyindex = false, lallindices = true;
      long int indexvalues[3] = { 0 };
      int indexint = 0;

      for (auto const &index : ensindexCMIP5)
      {
        std::string indexstring = kv_get_a_val(kvl, index, "");
        if (!indexstring.empty())
        {
          indexvalues[indexint] = std::stol(indexstring);
          if (!indexvalues[indexint])
          {
            cdo_warning("Could not parse value '%s' of attribute '%s' to integer.", indexstring, index);
            lallindices = false;
          }
          else
            lanyindex = true;
        }
        else
          lallindices = false;
        indexint++;
      }
      if (lanyindex && !vlabel.empty())
        cdo_warning("You specified both variant_label and all indices attributes.\n"
                    "Note that indices attributes have higher priority.");
      else if (!lallindices && vlabel.empty())
        cdo_abort("ERROR (infile: '%s')!\n"
                  "Could not find all index values required for describing the ensemble member (rip).\n"
                  "Make sure you have specified either 'variant_label' or all 3 'rip' indexes.",
                  cdo_get_stream_name(0));

      if (!lanyindex)
      {
        int scanres = std::sscanf((const char *) vlabel.c_str(), "r%ldi%ldp%ld", &indexvalues[0], &indexvalues[1], &indexvalues[2]);
        if (!scanres)
          cdo_abort("ERROR (infile: '%s')!\n"
                    "Could not scan all integers from attribute 'member'. \n"
                    "Make sure it has the format 'rINTiINTpINT'",
                    cdo_get_stream_name(0));
        indexint = 0;
        for (auto const &index : ensindexCMIP5)
        {
          std::string index2string = std::to_string(indexvalues[indexint]);
          kv_insert_vals(kvl, index, index2string, true, false);
          indexint++;
        }
      }
      else if (!lallindices)
      {
        long int vlvalues[3] = { 0 };
        int scanres = std::sscanf((const char *) vlabel.c_str(), "r%ldi%ldp%ld", &vlvalues[0], &vlvalues[1], &vlvalues[2]);
        if (!scanres)
          cdo_warning("Could not scan all integers from attribute 'member'. \n"
                      "Make sure it has the format 'rINTiINTpINT'");
        indexint = 0;
        while (indexint < 3)
        {
          if (!indexvalues[indexint])
          {
            indexvalues[indexint] = vlvalues[indexint];
            if (!indexvalues[indexint])
              cdo_abort("ERROR (infile: '%s')!\n"
                        "You did not provide a value for attribute '%s'.\n"
                        "Make sure you have specified either 'variant_label' or all 3 'rip' indexes.",
                        cdo_get_stream_name(0), ensindexCMIP5[indexint]);
            std::string index2string = std::to_string(indexvalues[indexint]);
            kv_insert_vals(kvl, ensindexCMIP5[indexint], index2string, true, false);
          }
          indexint++;
        }
      }
      char member[CMOR_MAX_STRING];
      std::snprintf(member, sizeof(member), "r%ldi%ldp%ld", indexvalues[0], indexvalues[1], indexvalues[2]);
      kv_insert_vals(kvl, "member", member, true, false);
      kv_insert_vals(kvl, "variant_label", member, true, false);
      indexint = 0;
    }
  }
  else
  {
    std::string vlabel = kv_get_a_val(kvl, "variant_label", "");
    bool lanyindex = false, lallindices = true;
    long int indexvalues[4] = { 0 };
    int indexint = 0;

    for (auto const &index : ensindexCMIP6)
    {
      std::string indexstring = kv_get_a_val(kvl, index, "");
      if (!indexstring.empty())
      {
        indexvalues[indexint] = std::stol(indexstring);
        if (!indexvalues[indexint])
        {
          cdo_warning("Could not parse value '%s' of attribute '%s' to integer.", indexstring, index);
          lallindices = false;
        }
        else
          lanyindex = true;
      }
      else
        lallindices = false;
      indexint++;
    }
    if (lanyindex && !vlabel.empty())
      cdo_warning("You specified both variant_label and all indices attributes.\n"
                  "Note that indices attributes have higher priority.");
    else if (!lallindices && vlabel.empty())
      cdo_abort("ERROR (infile: '%s')!\n"
                "Could not find all index values required for describing the ensemble member (ripf).\n"
                "Make sure you have specified either 'variant_label' or all 4 'ripf' indexes.",
                cdo_get_stream_name(0));

    if (!lanyindex)
    {
      int scanres = std::sscanf((const char *) vlabel.c_str(), "r%ldi%ldp%ldf%ld", &indexvalues[0], &indexvalues[1],
                                &indexvalues[2], &indexvalues[3]);
      if (!scanres)
        cdo_abort("ERROR (infile: '%s')!\n"
                  "Could not scan all integers from attribute 'variant_label'. \n"
                  "Make sure it has the format 'rINTiINTpINTfINT'",
                  cdo_get_stream_name(0));
      indexint = 0;
      for (auto const &index : ensindexCMIP6)
      {
        std::string index2string = std::to_string(indexvalues[indexint]);
        kv_insert_vals(kvl, index, index2string, true, false);
        indexint++;
      }
    }
    else if (!lallindices)
    {
      long int vlvalues[4] = { 0 };
      int scanres
          = std::sscanf((const char *) vlabel.c_str(), "r%ldi%ldp%ldf%ld", &vlvalues[0], &vlvalues[1], &vlvalues[2], &vlvalues[3]);
      if (!scanres)
        cdo_warning("Could not scan all integers from attribute 'variant_label'. \n"
                    "Make sure it has the format 'rINTiINTpINTfINT'");
      indexint = 0;
      while (indexint < 4)
      {
        if (!indexvalues[indexint])
        {
          indexvalues[indexint] = vlvalues[indexint];
          if (!indexvalues[indexint])
            cdo_abort("ERROR (infile: '%s')!\n"
                      "You did not provide a value for attribute '%s'.\n"
                      "Make sure you have specified either 'variant_label' or all 4 'ripf' indexes.",
                      cdo_get_stream_name(0), ensindexCMIP6[indexint]);
          std::string index2string = std::to_string(indexvalues[indexint]);
          kv_insert_vals(kvl, ensindexCMIP6[indexint], index2string, true, false);
        }
        indexint++;
      }
    }
    char member[CMOR_MAX_STRING];
    std::snprintf(member, sizeof(member), "r%ldi%ldp%ldf%ld", indexvalues[0], indexvalues[1], indexvalues[2], indexvalues[3]);
    kv_insert_vals(kvl, "member", member, true, false);
    kv_insert_vals(kvl, "variant_label", member, true, false);
  }
}

static void
read_config_files(KVList *kvl)
{
  if (Options::cdoVerbose) cdo_print("1. Start to read configuration files.");
  /* Files from info key in command line. */
  const KeyValues *info = kvl->search("i");
  int i = 0;
  if (info)
    while (i < info->nvalues)
    {
      if (Options::cdoVerbose) cdo_print("1.1. Try to parse file: '%s' configured with key 'info'.", info->values[i]);
      if (parse_kv_file(kvl, info->values[i].c_str()) == 0)
        cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), info->values[i]);
      if (Options::cdoVerbose) cdo_print("1.1. Successfully parsed file: '%s' configured with key 'info'.", info->values[i]);
      i++;
    }
  else if (Options::cdoVerbose)
    cdo_print("1.1. No info file was passed to the operator.");
  /* Config file in user's $cwd directory. */
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  cwd[strlen(cwd)] = '\0';
  const std::string dotconfig = ".cdocmorinfo";
  std::string workfile = std::string(cwd) + "/" + dotconfig;
  if (Options::cdoVerbose) cdo_print("1.2. Try to parse default file: '%s'.", workfile);
  if (parse_kv_file(kvl, workfile) == 0 && Options::cdoVerbose)
    cdo_warning("Default file for keyword 'info': '%s' does not exist.", workfile);
  else if (Options::cdoVerbose)
    cdo_print("1.2. Successfully parsed default file: '%s'.", workfile);

  if (i == 0)
  {
    const KeyValues *info2 = kvl->search("i");
    if (info2)
      while (i < info2->nvalues)
      {
        if (Options::cdoVerbose)
          cdo_print("1.3. Try to parse file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
        if (parse_kv_file(kvl, info2->values[i].c_str()) == 0)
          cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), info2->values[i]);
        if (Options::cdoVerbose)
          cdo_print("1.3. Successfully parsed file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
        i++;
      }
  }
  if (Options::cdoVerbose) cdo_print("1. Successfully read configuration files.");
}

static int
in_list(std::vector<std::string> const &list, std::string const &needle, int num)
{
  for (int i = 0; i < num; ++i)
    if (list[i] == needle) return 1;
  return 0;
}

static int
get_netcdf_file_action(KVList *kvl, std::string const &proj)
{
  (void) proj;
  std::string chunk = kv_get_a_val(kvl, "om", "a");
  if (chunk[0] == 'r')
  {
#if (CMOR_VERSION_MAJOR == 3)
    return CMOR_REPLACE;
#else
    return CMOR_REPLACE_4;
#endif
  }
  else if (chunk[0] == 'a')
  {
#if (CMOR_VERSION_MAJOR == 3)
    return CMOR_APPEND;
#else
    return CMOR_APPEND_4;
#endif
  }
  else if (chunk[0] == 'p')
  {
#if (CMOR_VERSION_MAJOR == 3)
    return CMOR_PRESERVE;
#else
    return CMOR_PRESERVE_4;
#endif
  }
  else
  {
    cdo_warning("You set output_mode = '%s', but valid are 'a' for append ,'r' for replace or 'p' for preserve.\n "
                "         CMOR output mode is set to: replace.",
                chunk);
    return CMOR_REPLACE;
  }
}

static int
get_cmor_verbosity(KVList *kvl)
{
  std::string verbos = kv_get_a_val(kvl, "set_verbosity", "");
  if (verbos.empty()) return CMOR_NORMAL;
  if (verbos == "CMOR_QUIET")
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int
get_cmor_exit_control(KVList *kvl)
{
  std::string exit = kv_get_a_val(kvl, "exit_control", "");
  if (exit.empty()) return CMOR_NORMAL;
  if (exit == "cmor_exit_on_major")
    return CMOR_EXIT_ON_MAJOR;
  else if (exit == "cmor_exit_on_warning")
    return CMOR_EXIT_ON_WARNING;
  else
    return CMOR_NORMAL;
}

static std::string
get_calendar_ptr(int calendar)
{
  std::string calendar_ptr;
  switch (calendar)
  {
    case CALENDAR_STANDARD: calendar_ptr = "standard"; break;
    case CALENDAR_GREGORIAN: calendar_ptr = "gregorian"; break;
    case CALENDAR_PROLEPTIC: calendar_ptr = "proleptic_gregorian"; break;
    case CALENDAR_360DAYS: calendar_ptr = "360_day"; break;
    case CALENDAR_365DAYS: calendar_ptr = "noleap"; break;
    case CALENDAR_366DAYS: calendar_ptr = "all_leap"; break;
    default: calendar_ptr = "";
  }
  return calendar_ptr;
}

static int
get_calendar_int(std::string const &calendar)
{
  if (calendar.empty())
    return -1;
  else if (calendar == "standard")
    return CALENDAR_STANDARD;
  else if (calendar == "gregorian")
    return CALENDAR_GREGORIAN;
  else if (calendar == "proleptic_gregorian")
    return CALENDAR_PROLEPTIC;
  else if (calendar == "360_day")
    return CALENDAR_360DAYS;
  else if (calendar == "noleap")
    return CALENDAR_365DAYS;
  else if (calendar == "all_leap")
    return CALENDAR_366DAYS;
  else
  {
    cdo_warning("You set calendar type = '%s' which is not supported by CMOR.", calendar);
    return -1;
  }
}

/***********************************************/
/*Time related functions************************/
/***********************************************/

static std::string
get_time_units(int taxisID)
{
  std::string units;
  units.resize(CMOR_MAX_STRING);
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second, ms;
  auto rDateTime = taxisInqRdatetime(taxisID);
  cdiDate_decode(rDateTime.date, &year, &month, &day);
  cdiTime_decode(rDateTime.time, &hour, &minute, &second, &ms);
  if (timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES) timeunit = TUNIT_MINUTE;
  if (timeunit == TUNIT_3HOURS || timeunit == TUNIT_6HOURS || timeunit == TUNIT_12HOURS) timeunit = TUNIT_HOUR;

  std::snprintf(&units[0], units.size(), "%s since %d-%d-%d %02d:%02d:%02d", tunitNamePtr(timeunit), year, month, day, hour, minute,
                second);
  units.resize(std::strlen(units.c_str()));
  return units;
}

static int
get_time_step_int(std::string const &time_step)
{
  if (time_step == "seconds")
    return TUNIT_SECOND;
  else if (time_step == "minutes")
    return TUNIT_MINUTE;
  else if (time_step == "hours")
    return TUNIT_HOUR;
  else if (time_step == "days")
    return TUNIT_DAY;
  else if (time_step == "months")
    return TUNIT_MONTH;
  else if (time_step == "years")
    return TUNIT_YEAR;
  else
  {
    cdo_warning("You set required_time_units = '%s since...'.\n          This time step is not yet implemented in cmor.",
                time_step);
    return 0;
  }
}

static int
check_time_units(std::string const &time_units)
{
  /* Required attribute in check_att */
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char time_step[CMOR_MAX_STRING];
  if (std::sscanf(time_units.c_str(), "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s", time_step, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond)
      != 7)
  {
    cdo_warning("You set required_time_units = '%s'\n          but it requires the form 'timestep since "
                "year-month-day hour:minute:second.\n          Could not read all 7 required time unit values.",
                time_units);
    return 0;
  }
  std::string cpptime(time_step);
  if (!get_time_step_int(cpptime)) return 0;
  return 1;
}

static void
get_time_method(KVList *kvl, int vlistID, int varID, std::string &cmor_time_name, std::string const &project_id, int miptab_freq,
                int *time_axis)
{
  if ((project_id == "CMIP5" || project_id == "CMIP6") && miptab_freq) switch (miptab_freq)
    {
      case 1:
        cmor_time_name = "time2";
        *time_axis = 2;
        break;
      case 2:
        cmor_time_name = "time";
        *time_axis = 0;
        break;
      case 4:
        cmor_time_name = "time";
        *time_axis = 0;
        break;
      case 5:
        cmor_time_name = "time1";
        *time_axis = 1;
        break;
      case 6:
        cmor_time_name = "time1";
        *time_axis = 1;
        break;
      case 7:
        cmor_time_name = "time3";
        *time_axis = 3;
        break;
    }
  if (cmor_time_name.empty())
  {
    std::string time_method = get_txtatt(vlistID, varID, "cell_methods");
    std::string att_time_method = kv_get_a_val(kvl, "cm", "");
    check_compare_set(time_method, att_time_method, "cell_methods", " ");
    if (time_method[0] == 'm')
    {
      cmor_time_name = "time \0";
      *time_axis = 0;
    }
    else if (time_method[0] == 'p')
    {
      cmor_time_name = "time1\0";
      *time_axis = 1;
    }
    else if (time_method[0] == 'c')
    {
      cmor_time_name = "time2\0";
      *time_axis = 2;
    }
    else if (time_method[0] == 'd')
    {
      cmor_time_name = "time3\0";
      *time_axis = 3;
    }
    else if (time_method[0] == 'n')
    {
      cmor_time_name = "none\0";
      *time_axis = 4;
    }
    else
    {
      if (time_method[0] == ' ')
      {
        if (Options::cdoVerbose)
          cdo_print("No value found for attribute 'cell_methods' of variable with ID '%d'.\n          Depending on "
                    "the aggregation method, you can specifiy one of: \n          'n', 'm', 'p', 'c', 'd'.\n      "
                    "    The default ('m') is used.",
                    varID);
      }
      else
        cdo_warning("The value for attribute cell method = '%s' of variable with ID '%d' is not valid.\n          "
                    "Depending on the aggregation method, you can specifiy one of: \n          'n', 'm', 'p', 'c', "
                    "'d'.\n          The default ('m') is used.",
                    time_method, varID);
      cmor_time_name = "time \0";
    }
  }
  if (Options::cdoVerbose) cdo_print("Successfully determined time_axis = '%d'.", *time_axis);
}

static CdiDateTime
get_taxis(std::string const &required_time_units, int *timeunit)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char atttimeunit[CMOR_MAX_STRING];

  std::sscanf(required_time_units.c_str(), "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s", atttimeunit, &attyear, &attmonth, &attday,
              &atthour, &attminute, &attsecond);
  std::string cppatt(atttimeunit);
  *timeunit = get_time_step_int(cppatt);
  CdiDateTime sDateTime{};
  sDateTime.date = cdiDate_encode(attyear, attmonth, attday);
  sDateTime.time = cdiTime_encode(atthour, attminute, attsecond, 0);
  return sDateTime;
}

static double *
get_branch_times(KVList *kvl, int calendar, std::string const &time_units, std::string const &project_id)
{
  if (Options::cdoVerbose) cdo_print("6.1.2. Start to compute attribute 'branch_time'.");
  std::string btip = kv_get_a_val(kvl, "branch_time_in_parent", "");
  std::string btic = kv_get_a_val(kvl, "branch_time_in_child", "");
  std::string ptu = kv_get_a_val(kvl, "parent_time_units", "");
  double *branch_time = (double *) std::malloc(2 * sizeof(double));
  branch_time[0] = 0.0;
  branch_time[1] = 0.0;

  if (kv_get_a_val(kvl, "parent_experiment_id", "no parent") != "no parent"
      && (btip.empty() || (project_id == "CMIP6" && btic.empty())))
  {

    int numdates = 0;
    std::vector<std::string> branch_dates_p = kv_get_vals(kvl, "branch_dates", &numdates);

    if (numdates == 2 && !ptu.empty())
    {
      int branchdates[2], branchyears[2], branchmonths[2], branchdays[2];
      CdiDateTime branchDateTimes[2]{};
      for (int i = 0; i < 2; ++i)
      {
        branchdates[i] = std::stol(branch_dates_p[i]);
        branchyears[i] = branchdates[i] / 100 / 100;
        branchmonths[i] = (branchdates[i] - branchyears[i] * 100 * 100) / 100;
        branchdays[i] = branchdates[i] - branchyears[i] * 100 * 100 - branchmonths[i] * 100;
        branchDateTimes[i].date = cdiDate_encode(branchyears[i], branchmonths[i], branchdays[i]);
      }

      int parenttimeunit;
      auto parentsDateTime = get_taxis(ptu, &parenttimeunit);
      auto parentstartdate = julianDate_encode(calendar, parentsDateTime);
      auto parentbranchdate = julianDate_encode(calendar, branchDateTimes[0]);

      int childtimeunit;
      auto childsDateTime = get_taxis(time_units, &childtimeunit);
      auto childstartdate = julianDate_encode(calendar, childsDateTime);
      auto childbranchdate = julianDate_encode(calendar, branchDateTimes[1]);

      /* If time unit is always "days since.." */
      branch_time[0] = julianDate_to_seconds(julianDate_sub(parentbranchdate, parentstartdate)) / 86400;
      branch_time[1] = julianDate_to_seconds(julianDate_sub(childbranchdate, childstartdate)) / 86400;
    }
    else
    {
      cdo_warning("Since you specified 'parent_experiment_id' you have to set 'parent_time_units' and either 'branch_dates' or "
                  "both 'branch_time_in_parent' and 'branch_time_in_child'. Could not find a correct configuration.");
    }
  }
  else
  {
    if (!btip.empty()) branch_time[0] = std::stol(btip);
    if (!btic.empty()) branch_time[1] = std::stol(btic);
  }
  if (Options::cdoVerbose)
    cdo_print("6.1.2. Successfully computed 'branch_time_in_parent'='%f' and 'branch_time_in_child'='%f'.", branch_time[0],
              branch_time[1]);
  return branch_time;
}

static std::string
check_required_time_units(KVList *kvl, int taxisID)
{
  if (Options::cdoVerbose) cdo_print("4.1. Start to check attribute 'required_time_units'.");
  std::string time_units = get_time_units(taxisID);
  std::string required_time_units = kv_get_a_val(kvl, "rtu", "");
  if (!required_time_units.empty() && check_time_units(required_time_units))
    check_compare_set(time_units, required_time_units, "time_units", "");
  else
    cdo_warning("Required Attribute 'required_time_units' from configuration is invalid!\n          Continue with infile time "
                "units instead.");
  kv_insert_vals(kvl, "rtu", time_units, true, false);
  if (Options::cdoVerbose) cdo_print("4.1. Successfully checked attribute 'required_time_units'.");
  return time_units;
}

static std::string
check_calendar(KVList *kvl, int taxisID, int *calendar)
{
  if (Options::cdoVerbose) cdo_print("6.1.1. Start to check attribute 'calendar'.");
  std::string attcalendar = kv_get_a_val(kvl, "calendar", "");
  std::string calendarptr = get_calendar_ptr(taxisInqCalendar(taxisID));
  if ((*calendar = get_calendar_int(attcalendar)) > -1)
    check_compare_set(calendarptr, attcalendar, "calendar", "");
  else if (get_calendar_int(calendarptr) < 0)
  {
    if (!attcalendar.empty()) cdo_print("6.1.1. Cannot use calendar '%s' from configuration.", attcalendar);
    if (!calendarptr.empty()) cdo_print("6.1.1. Cannot use calendar '%s' from infile.", calendarptr);
    if (Options::cdoVerbose)
      cdo_print("6.1.1. You did not provide a valid calendar. Default 'standard' is used. Valid calendars are:\n          "
                "'standard', 'gregorian', 'proleptic_gregorian', '360_day', 'noleap' and 'all_leap'.");
    calendarptr = "standard";
  }
  else if (get_calendar_int(calendarptr) > -1 && Options::cdoVerbose)
    cdo_print("6.1.1. Use calendar '%s' from infile.", calendarptr);
  if (Options::cdoVerbose) cdo_print("6.1.1. Successfully retrived calendar: '%s'.", calendarptr);
  return calendarptr;
}

/*********/
/* main: */
/*********/

static bool
keep_this_attribute(KeyValues *kv, std::vector<std::string> const &array)
{
  if (kv->key.size() && kv->nvalues == 1)
  {
    size_t j = 0;
    for (j = 0; j < array.size(); ++j)
    {
      if (kv->key == array[j]) break;
    }
    if (j == array.size() && (kv->values[0] != "notSet")) { return true; }
    else
      return false;
  }
  else
    return false;
}

#if (CMOR_VERSION_MAJOR == 3)

static const std::string
copyCV(std::string directory)
{
  const std::string cvwithout = "CMIP6_CV_without_prefix.json";
  char cvname[CDI_MAX_NAME];
  std::snprintf(cvname, CDI_MAX_NAME, "%s/%s", directory.c_str(), cvwithout.c_str());
  if (Options::cdoVerbose) cdo_print("Check whether a CV without tracking prefix check exists.");
  if (file_exist(std::string(cvname), false, "CV", false))
    return cvwithout;
  else
  {
    if (Options::cdoVerbose) cdo_print("Try to create a CV without tracking prefix check with program 'sed'.");

    char command[1024];
    std::snprintf(command, sizeof(command), "sed 's/\"hdl:21.14100\\/\\.\\*\"/\"\\^\\.*\"/' %s/CMIP6_CV.json >%s",
                  directory.c_str(), cvname);
    int dir_err = system(command);
    if (dir_err != 0)
    {
      if (Options::cdoVerbose)
        cdo_print("Creation of a CV without tracking prefix check failed. "
                  "This can be due to missing 'sed' program or missing write permissions. "
                  "Continue with tracking prefix. "
                  "Consider that the output tracking id needs to be registered as a PID.");
      std::snprintf(cvname, CDI_MAX_NAME, "CMIP6_CV.json");
      return cvname;
    }
    if (Options::cdoVerbose) cdo_print("Successfully created a CV without tracking prefix check.");
    return cvwithout;
  }
}
#endif

static void
setup_dataset(KVList *kvl, CdoStreamID streamID, int *calendar, std::string const &project_id)
{
  if (Options::cdoVerbose) cdo_print("6. Start to process cmor_setup and cmor_dataset.");
  int netcdf_file_action = get_netcdf_file_action(kvl, project_id);
  int set_verbosity = get_cmor_verbosity(kvl);
  int exit_control = get_cmor_exit_control(kvl);
  int creat_subs = 1;
  std::string drs = kv_get_a_val(kvl, "d", "y");
  if (drs[0] == 'n')
    creat_subs = 0;
  else if (drs[0] != 'y')
  {
    cdo_warning("In preparing cmor_setup:\n          You set 'd' = '%s' which is not valid.\n          Allowed are: "
                "'n' or 'y'. d is set to 'y'.",
                drs);
    kv_insert_vals(kvl, "d", "y", true, false);
  }
  if (project_id == "CORDEX" && creat_subs)
  {
    if (Options::cdoVerbose)
      cdo_print("Since you produce CORDEX compliant output, a path is constructed with DRS elements according to the "
                "project defined template.");

    std::string miptabfreq = kv_get_a_val(kvl, "miptab_freq", "");
    std::string freq;
    if (miptabfreq == "6h")
      freq = "6hr";
    else if (miptabfreq == "3h")
      freq = "3hr";
    else if (miptabfreq == "1h")
      freq = "1hr";
    else
      freq = miptabfreq;

    char cordexDir[CDI_MAX_NAME];
    char cordexFileTem[CDI_MAX_NAME];
    std::snprintf(cordexDir, CDI_MAX_NAME, "%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s", kv_get_a_val(kvl, "dr", "./").c_str(),
                  project_id.c_str(), to_cmor(kv_get_a_val(kvl, "product", "")), to_cmor(kv_get_a_val(kvl, "CORDEX_domain", "")),
                  to_cmor(kv_get_a_val(kvl, "institute_id", "")), to_cmor(kv_get_a_val(kvl, "driving_model_id", "")),
                  to_cmor(kv_get_a_val(kvl, "experiment_id", "")), to_cmor(kv_get_a_val(kvl, "member", "")),
                  to_cmor(kv_get_a_val(kvl, "model_id", "")), to_cmor(kv_get_a_val(kvl, "rcm_version_id", "")), freq.c_str());
    std::snprintf(cordexFileTem, CDI_MAX_NAME, "%s_%s_%s_%s_%s_%s_%s", to_cmor(kv_get_a_val(kvl, "CORDEX_domain", "")),
                  to_cmor(kv_get_a_val(kvl, "driving_model_id", "")), to_cmor(kv_get_a_val(kvl, "experiment_id", "")),
                  to_cmor(kv_get_a_val(kvl, "member", "")), to_cmor(kv_get_a_val(kvl, "model_id", "")),
                  to_cmor(kv_get_a_val(kvl, "rcm_version_id", "")), freq.c_str());

    kv_insert_vals(kvl, "dr", "./", true, false);
    kv_insert_vals(kvl, "cordexDir", std::string(cordexDir), true, false);
    kv_insert_vals(kvl, "cordexFileTem", std::string(cordexFileTem), true, false);
    creat_subs = 0;
  }

  int vlistID = cdo_stream_inq_vlist(streamID);

  int cmf = cmor_setup(to_cmor(kv_get_a_val(kvl, "inpath", "/usr/share/cmor/")), &netcdf_file_action, &set_verbosity, &exit_control,
                       to_cmor(kv_get_a_val(kvl, "logfile", "")), &creat_subs);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_setup failed!", cdo_get_stream_name(0));

  signal(SIGTERM, sigfunc);
  int taxisID = vlistInqTaxis(vlistID);

  /*
    char *attcomment = kv_get_a_val(kvl, "comment", nullptr);
    char *comment = get_txtatt(vlistID, CDI_GLOBAL, "comment");
  */

  /* First compare file calendar and config calendar and retrieve pointer and integer
     Then check the required time units from config and retrieve
     Then compute branch_time_in_parent and branch_time_in_child */

  if (Options::cdoVerbose)
    cdo_print("6.1. Start to check model calendar as well as 'required_time_units' and 'branch_times' attributes.");
  std::string calendarptr = check_calendar(kvl, taxisID, calendar);
  std::string time_units = kv_get_a_val(kvl, "rtu", "");
  double *branch_times = get_branch_times(kvl, *calendar, time_units, project_id);

  if (Options::cdoVerbose) cdo_print("6.1. Successfully found valid calendar, 'required_time_units' and 'branch_times'.");

  /* if keep_all_atts is set: */
  std::vector<std::string> kaa_notneeded_general
      = { "cn", "n", "c", "u", "cm", "vc", "p", "i", "ca", "za", "gi", "rtu", "mt", "om", "ms", "dr", "d", "lc", "dj",
          "workfile4err", "kaa", "mtproof", "miptab_freq", "mip_table_dir", "grid_info_dir", "mapping_table_dir", "branch_dates",
          "member", "ta", "firsttimeval", "calendar", "cordexDir", "cordexFileTem", "branch_time", "dl", "vd",
          /*Following attributes are set by CMOR: */
          "tracking_id", "creation_date", "table_id", "di", "ci", "sc", "ml", "grid_file", "cdo_openmp_thread_number",
          "grid_north_pole_latitude", "grid_north_pole_longitude", "switch_xy", "switch_z" };
#if defined(CMOR_VERSION_MAJOR)
#if (CMOR_VERSION_MAJOR == 2)
  std::vector<std::string> kaa_cmor2 = { "parent_experiment", "modeling_realm" };
  std::vector<std::string> datasetvals = {
    "dr",
    "experiment_id",
    "institution",
    "source",
    "realization",
    "contact",
    "history",
    "comment",
    "references",
    "leap_year",
    "leap_month",
    "model_id",
    "forcing",
    "initialization_method",
    "physics_version",
    "institute_id",
    "parent_experiment_id",
    "parent_experiment_rip",
    "product",
  };
  cmf = cmor_dataset((char *) kv_get_a_val(kvl, datasetvals[0], "./").c_str(), to_cmor(kv_get_a_val(kvl, datasetvals[1], "")),
                     to_cmor(kv_get_a_val(kvl, datasetvals[2], "")), to_cmor(kv_get_a_val(kvl, datasetvals[3], "")),
                     calendarptr.c_str(), std::stol(kv_get_a_val(kvl, datasetvals[4], "")),
                     to_cmor(kv_get_a_val(kvl, datasetvals[5], "")), to_cmor(kv_get_a_val(kvl, datasetvals[6], "")),
                     to_cmor(kv_get_a_val(kvl, datasetvals[7], "")), to_cmor(kv_get_a_val(kvl, datasetvals[8], "")),
                     std::stol(kv_get_a_val(kvl, datasetvals[9], "")), std::stol(kv_get_a_val(kvl, datasetvals[10], "")), nullptr,
                     to_cmor(kv_get_a_val(kvl, datasetvals[11], "")), to_cmor(kv_get_a_val(kvl, datasetvals[12], "")),
                     std::stol(kv_get_a_val(kvl, datasetvals[13], "")), std::stol(kv_get_a_val(kvl, datasetvals[14], "")),
                     to_cmor(kv_get_a_val(kvl, datasetvals[15], "")), to_cmor(kv_get_a_val(kvl, datasetvals[16], "")),
                     &(branch_times[0]), to_cmor(kv_get_a_val(kvl, datasetvals[17], "")));
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_dataset failed!", cdo_get_stream_name(0));
  std::vector<std::string> allneeded2 = {
    "CORDEX_domain",           "driving_experiment", "driving_model_id", "driving_model_ensemble_member",
    "driving_experiment_name", "rcm_version_id",
  };
  if (project_id == "CORDEX")
    for (size_t ind = 0; ind < allneeded2.size(); ++ind)
    {
      std::string tmp = kv_get_a_val(kvl, allneeded2[ind], "");
      if (!tmp.empty()) cmf = cmor_set_cur_dataset_attribute((char *) allneeded2[ind].c_str(), tmp.c_str(), 1);
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_cur_dataset_attribute failed!", cdo_get_stream_name(0));
    }
  if ((kv_get_a_val(kvl, "kaa", "n") == "y")
      || ((project_id != "CMIP6")) /* && (project_id != "PalMod2") && (project_id != "CORDEX-CMIP6") && (project_id != "EERIE")) */
  )
  {
    char notincluded[2048];
    std::strcpy(
        notincluded,
        "The following attributes are not included in the global attributes list.\n          Reasons can be: 1. Attribute is "
        "an internal keyword 2. No valaue is available 3. CMOR creates the attribute itself:\n          ");
    size_t inilen = std::strlen(notincluded);
    size_t strlens = inilen;
    for (auto &kv : *kvl)
    {
      if (keep_this_attribute(&kv, kaa_notneeded_general) && keep_this_attribute(&kv, allneeded2)
          && keep_this_attribute(&kv, datasetvals) && keep_this_attribute(&kv, kaa_cmor2))
      {
        cmf = cmor_set_cur_dataset_attribute((char *) kv.key.c_str(), (char *) kv.values[0].c_str(), 1);
        if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_cur_dataset_attribute failed!", cdo_get_stream_name(0));
      }
      else if (!keep_this_attribute(&kv, kaa_notneeded_general) || !keep_this_attribute(&kv, kaa_cmor2))
      {
        strlens += (std::strlen(kv.key.c_str()) + 2);
        if (strlens < 2048)
        {
          std::strcat(notincluded, kv.key.c_str());
          std::strcat(notincluded, ", ");
        }
      }
    }
    if (strlens > inilen)
      if (Options::cdoVerbose) cdo_print("%s", notincluded);
  }
#elif (CMOR_VERSION_MAJOR == 3)
  {
    /***/
    /* Could not give CMOR all attributes separately because some are required to be in a json file (outpath,...). */
    /* Better collect them in this file. */
    /* todo this **/
    /* If a Json file is denoted, read this file and check attributes */
    /***/

    /*
          std::string filename = kv_get_a_val(kvl, "dj", nullptr); */

    size_t cwdsize = 1024;
    char cwd[CMOR_MAX_STRING];
    getcwd(cwd, cwdsize);
    cwd[strlen(cwd)] = '\0';
    int procID = getpid();
    char dataset_path[CMOR_MAX_STRING];

    std::snprintf(dataset_path, CMOR_MAX_STRING, "%s/dataset%d.json", cwd, procID);
    auto dataset_json = std::fopen(dataset_path, "w+");
    if (!dataset_json)
      cdo_abort("ERROR (infile: '%s')! In preparing cmor_dataset:\n          Could not open a dataset file '%s' for cmor_dataset.",
                cdo_get_stream_name(0), dataset_path);
    fputs("{\n", dataset_json);

    int i = 0;
    if ((kv_get_a_val(kvl, "kaa", "n") == "y")
        || ((project_id
             != "CMIP6")) /* && (project_id != "PalMod2") && (project_id != "CORDEX-CMIP6") && (project_id != "EERIE")) */
    )
    {
      for (auto &kv : *kvl)
      {
        if (keep_this_attribute(&kv, kaa_notneeded_general))
        {
          int linelen = std::strlen(kv.key.c_str()) + std::strlen(kv.values[0].c_str()) + 10;
          std::vector<char> line(linelen);
          std::snprintf(line.data(), linelen, "\"%s\" : \"%s\",\n", kv.key.c_str(), kv.values[0].c_str());
          fputs((const char *) line.data(), dataset_json);
          if (Options::cdoVerbose) cdo_print("%s", (const char *) line.data());
        }
      }
    }
    else
    {
      std::vector<std::string> allneeded
          = /*CMIP5*/ { "project_id", "experiment_id", "institution", "source", "realization", "contact", "history", "comment",
                        "references", "leap_year", "leap_month", "source_id", "model_id", "forcing", "initialization_method",
                        "modeling_realm", "physics_version", "institute_id", "parent_experiment_rip",
                        /*CORDEX */
                        "CORDEX_domain", "driving_experiment", "driving_model_id", "driving_model_ensemble_member",
                        "driving_experiment_name", "rcm_version_id",
                        /* CMIP6: */
                        /* Glob Atts */
                        "_cmip6_option", "Conventions", "activity_id", "branch_method", "experiment", "experiment_id",
                        "forcing_index", "further_info_url", "grid", "grid_label", "initialization_index", "institution",
                        "institution_id", "license", "mip_era", "nominal_resolution", "physics_index", "product",
                        "realization_index", "source", "source_id", "source_type", "sub_experiment", "sub_experiment_id",
                        "table_id", "variant_label", "parent_experiment_id", "parent_activity_id", "parent_mip_era",
                        "parent_source_id", "parent_variant_label", "parent_time_units", "variant_info", "title",
                        /* Others */
                        "_controlled_vocabulary_file", "_FORMULA_VAR_FILE", "_AXIS_ENTRY_FILE" };
      for (i = 0; i < (int) allneeded.size(); ++i)
      {
        std::string tmp = kv_get_a_val(kvl, allneeded[i], "notSet");
        if (tmp.substr(0, 6) != "notSet")
        {
          int linelen = std::strlen(allneeded[i].c_str()) + std::strlen(tmp.c_str()) + 10;
          std::vector<char> line(linelen);
          std::snprintf(line.data(), line.size(), "\"%s\" : \"%s\",\n", allneeded[i].c_str(), tmp.c_str());
          fputs((const char *) line.data(), dataset_json);
        }
      }
    }

    char branch_time_in_parent[CMOR_MAX_STRING];
    char branch_time_in_child[CMOR_MAX_STRING];
    std::snprintf(branch_time_in_parent, CMOR_MAX_STRING, "%.12f", branch_times[0]);
    std::snprintf(branch_time_in_child, CMOR_MAX_STRING, "%.12f", branch_times[1]);

    /* CMOR internal */
    fputs("\"outpath\" : \"", dataset_json);
    fputs(kv_get_a_val(kvl, "dr", "./").c_str(), dataset_json);
    fputs("\",\n", dataset_json);
    fputs("\"output_path_template\" : \"", dataset_json);
    if (project_id == "PalMod2")
    {
      fputs(kv_get_a_val(kvl, "output_path_template",
                         "<activity_id><institution_id><source_id><experiment_id><"
                         "member_id><table><variable_id><grid_label>")
                .c_str(),
            dataset_json);
      fputs(kv_get_a_val(kvl, "vd", "<version>").c_str(), dataset_json);
      fputs("\",\n", dataset_json);
      fputs("\"_history_template\" : \"", dataset_json);
      fputs(kv_get_a_val(kvl, "_history_template",
                         "%s ; CMOR rewrote data to be consistent with <activity_id>, <Conventions> and CF standards.")
                .c_str(),
            dataset_json);
    }
    else
    {
      fputs(kv_get_a_val(kvl, "output_path_template",
                         "<mip_era><activity_id><institution_id><source_id><experiment_id><"
                         "member_id><table><variable_id><grid_label>")
                .c_str(),

            dataset_json);
      fputs(kv_get_a_val(kvl, "vd", "<version>").c_str(), dataset_json);
    }
    fputs("\",\n", dataset_json);
    fputs("\"output_file_template\" : \"", dataset_json);
    fputs(
        kv_get_a_val(kvl, "output_file_template", "<variable_id><table><source_id><experiment_id><member_id><grid_label>").c_str(),
        dataset_json);
    fputs("\",\n", dataset_json);
    /*          fputs("\"frequency\" : \"", dataset_json);
              fputs(freq, dataset_json);
              fputs("\",\n", dataset_json);  */

    /* cdo cmor preprocessed: */
    fputs("\"calendar\" : \"", dataset_json);
    fputs(calendarptr.c_str(), dataset_json);
    fputs("\",\n", dataset_json);
    if (kv_get_a_val(kvl, "parent_experiment_id", "no parent") != "no parent")
    {
      fputs("\"branch_time_in_parent\" : \"", dataset_json);
      fputs(branch_time_in_parent, dataset_json);
      fputs("\",\n", dataset_json);
      fputs("\"branch_time_in_child\" : \"", dataset_json);
      fputs(branch_time_in_child, dataset_json);
      fputs("\",\n", dataset_json);
    }

    if (kv_get_a_val(kvl, "tp", "y") == "y")
    {
      if (project_id == "CMIP6")
      {
        fputs("\"tracking_prefix\" : \"hdl:21.14100\"", dataset_json);
        fputs(",\n", dataset_json);
      }
      else if (project_id == "PalMod2")
      {
        fputs("\"tracking_prefix\" : \"hdl:21.14105\"", dataset_json);
        fputs(",\n", dataset_json);
      }
      else if (project_id == "CORDEX-CMIP6")
      {
        fputs("\"tracking_prefix\" : \"hdl:21.14103\"", dataset_json);
        fputs(",\n", dataset_json);
      }
      else if (project_id == "EERIE")
      {
        fputs("\"tracking_prefix\" : \"hdl:21.14102\"", dataset_json);
        fputs(",\n", dataset_json);
      }
      else
        cdo_warning("Don't know the tracking_prefix for project '%s'.\n          "
                    "Continue without setting a tracking prefix.",
                    project_id);
    }
    else if (kv_get_a_val(kvl, "_controlled_vocabulary_file", "notset") == "notset")
    {
      std::string cvname = copyCV(kv_get_a_val(kvl, "mip_table_dir", ""));
      fputs("\"_controlled_vocabulary_file\" : \"", dataset_json);
      fputs(cvname.c_str(), dataset_json);
      fputs("\",\n", dataset_json);
    }
    if (kv_get_a_val(kvl, "cvn", "y") == "y")
    {
      fputs("\"CDO\" : \"", dataset_json);
      fputs(cdo_comment(), dataset_json);
      fputs("\",\n", dataset_json);
    }
    if (kv_get_a_val(kvl, "cdi_grid_type", "n") == "unstructured")
    {
      fputs("\"CDI_grid_type\" : \"unstructured\"", dataset_json);
      fputs(",\n", dataset_json);
      fputs("\"grid_type\" : \"unstructured\"", dataset_json);
      fputs(",\n", dataset_json);
    }

    fputs("}\n", dataset_json);
    std::fclose(dataset_json);
    cmf = cmor_dataset_json(dataset_path);
    if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_dataset_json failed!", cdo_get_stream_name(0));

    removeDataset();
    /*      std::free(freq); */
  }

#else
  cdo_abort("ERROR (infile: '%s')! Cmor version %d not yet enabled!", cdo_get_stream_name(0), (int) CMOR_VERSION_MAJOR);
#endif
#else
  cdo_abort("ERROR (infile: '%s')! It is not clear which CMOR version is installed since\n          Makros CMOR_VERSION_MAJOR and "
            "CMOR_VERSION_MINOR are not available.",
            cdo_get_stream_name(0));
#endif
  std::free(branch_times);
  if (Options::cdoVerbose) cdo_print("6. Successfully finished cmor_setup and cmor_dataset.");
}

static void
gen_bounds(int n, double *vals, double *bounds)
{
  for (int i = 0; i < n - 1; ++i)
  {
    bounds[2 * i + 1] = 0.5 * (vals[i] + vals[i + 1]);
    bounds[2 * (i + 1)] = 0.5 * (vals[i] + vals[i + 1]);
  }

  bounds[0] = 2 * vals[0] - bounds[1];
  bounds[2 * n - 1] = 2 * vals[n - 1] - bounds[2 * (n - 1)];
}

static bool
get_zcell_bounds(int zaxisID, double *zcell_bounds, double *levels, int zsize)
{
  bool selfGenerated = false;
  double *lbounds;
  lbounds = (double *) std::malloc(zsize * sizeof(double));
  zaxisInqLbounds(zaxisID, lbounds);
  double *ubounds;
  ubounds = (double *) std::malloc(zsize * sizeof(double));
  zaxisInqUbounds(zaxisID, ubounds);
  if (!lbounds || !ubounds || std::pow((ubounds[1] - ubounds[0]), 2) < 0.001 || std::pow((lbounds[1] - lbounds[0]), 2) < 0.001)
  {
    gen_bounds(zsize, levels, zcell_bounds);
    selfGenerated = true;
  }
  else
  {
    if (lbounds)
      zcell_bounds[0] = lbounds[0];
    else
      zcell_bounds[0] = 0;
    for (int i = 0; i < zsize - 1; ++i)
    {
      zcell_bounds[2 * i + 1] = ubounds[i];
      zcell_bounds[2 * (i + 1)] = lbounds[i + 1];
    }
    if (ubounds)
      zcell_bounds[2 * zsize - 1] = ubounds[zsize - 1];
    else
      zcell_bounds[2 * zsize - 1] = levels[zsize - 1] + (levels[zsize - 1] - zcell_bounds[2 * zsize - 2]);
    selfGenerated = false;
  }
  std::free(lbounds);
  std::free(ubounds);
  return selfGenerated;
}

static void
get_zhybrid_half(int zaxisID, double *p0, double *alev_val, double *b_val, double *ap_val)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  if (vctsize < 1)
    cdo_abort("ERROR (infile: '%s')! Missing z-axis description. Please provide all parameters"
              " to calculate the formula of a hybrid sigma pressure axis:"
              "\n          p = ap + b *ps",
              cdo_get_stream_name(0));
  double *vct = (double *) std::malloc(vctsize * sizeof(double));
  zaxisInqVct(zaxisID, vct);
  for (int i = 0; i < zsize; ++i)
  {
    ap_val[i] = vct[i];
    b_val[i] = vct[zsize + i];
  }
  for (int i = 0; i < zsize; ++i) { alev_val[i] = ap_val[i] / p0[0] + b_val[i]; }
  std::free(vct);
}

static void
get_zhybrid(int zaxisID, double *p0, double *alev_val, double *alev_bnds, double *b_val, double *b_bnds, double *ap_val,
            double *ap_bnds)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  if (vctsize < 1)
    cdo_abort("ERROR (infile: '%s')! Missing z-axis description. Please provide all parameters"
              " to calculate the formula of a hybrid sigma pressure axis or a hybrid height axis:"
              "\n          p = ap + b * ps||orog",
              cdo_get_stream_name(0));
  double *vct = (double *) std::malloc(vctsize * sizeof(double));
  zaxisInqVct(zaxisID, vct);
  for (int i = 0; i < (zsize + 1); ++i)
  {
    ap_bnds[i] = vct[i];
    b_bnds[i] = vct[zsize + 1 + i];
  }
  for (int i = 0; i < zsize; ++i)
  {
    ap_val[i] = (ap_bnds[i] + ap_bnds[i + 1]) / 2.0;
    b_val[i] = (b_bnds[i] + b_bnds[i + 1]) / 2.0;
    alev_val[i] = ap_val[i] / p0[0] + b_val[i];
    alev_bnds[i] = ap_bnds[i] / p0[0] + b_bnds[i];
  }
  alev_bnds[zsize] = ap_bnds[zsize] / p0[0] + b_bnds[zsize];
  std::free(vct);
}

static size_t
get_strmaxlen(std::vector<std::string> const &array, size_t len)
{
  size_t result = 0;
  for (size_t i = 0; i < len; ++i)
    if (result < array[i].length()) result = array[i].length();
  return result;
}

static void
get_charvals_and_bnds(KVList *kvl, std::string const &chardim, std::vector<std::string> &fvalss, std::vector<std::string> &fbndss,
                      std::string &funits, int *nofvals, int *nofbnds, std::string cmor_name)
{
  bool fivedim = true;
  char charvalstring[CMOR_MAX_STRING];
  std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_%s", chardim.c_str(), cmor_name.c_str());
  fvalss = kv_get_vals(kvl, charvalstring, nofvals);
  if (!fvalss.size())
  {
    if (Options::cdoVerbose) cdo_print("Start to register char_axis_%s", chardim);
    std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s", chardim.c_str());
    fvalss = kv_get_vals(kvl, charvalstring, nofvals);
    fivedim = false;
  }
  else
  {
    if (Options::cdoVerbose) cdo_print("Start to register char_axis_%s_%s", chardim, cmor_name);
  }
  if (!fvalss.size())
    cdo_warning("You specify variables to merge to an axis and the axis name but its values are missing!"
                "\n          Specify its values via attribute 'char_axis_$name' in infofile.");

  if (fivedim)
    std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_%s_bounds", chardim.c_str(), cmor_name.c_str());
  else
    std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_bounds", chardim.c_str());
  fbndss = kv_get_vals(kvl, charvalstring, nofbnds);
  if (!fbndss.size())
    if (Options::cdoVerbose)
      cdo_print("You specified variables to merge to an axis, the axis name and its values."
                "\n          Note that sometimes bounds are required for these axes."
                "\n          Specify its bounds via attribute 'char_axis_$name_bounds' in infofile.");

  if (fivedim)
    std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_%s_units", chardim.c_str(), cmor_name.c_str());
  else
    std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_units", chardim.c_str());
  funits = kv_get_a_val(kvl, std::string(charvalstring), "");
  if (funits.empty())
    if (Options::cdoVerbose)
      cdo_print("You specified variables to merge to an axis, the axis name and its values."
                "\n          Note that units are required if the axis has digital values."
                "\n          Specify its units via attribute 'char_axis_$name_units' in infofile.");
}

static void
register_char_axis(int numchar, std::vector<std::string> const &charvals, int *axis_ids, std::string chardim)
{
  if (Options::cdoVerbose) cdo_print("Start to register a character axis.");
  size_t maxlen = get_strmaxlen(charvals, numchar);
  size_t len = numchar * maxlen + 1;
  void *charcmor = (void *) std::malloc(len);
  std::snprintf((char *) charcmor, len, "%.*s", (int) charvals[0].size(), charvals[0].c_str());
  std::vector<char> blanks(maxlen);
  for (size_t i = 0; i < maxlen; ++i) blanks[i] = ' ';
  char tempb[CMOR_MAX_STRING];
  std::snprintf(tempb, CMOR_MAX_STRING, "%.*s", (int) (maxlen - charvals[0].size()), blanks.data());
  std::strcat((char *) charcmor, tempb);
  for (int i = 1; i < numchar; ++i)
  {
    std::strcat((char *) charcmor, charvals[i].c_str());
    char tempblanks[CMOR_MAX_STRING];
    std::snprintf(tempblanks, CMOR_MAX_STRING, "%.*s", (int) (maxlen - charvals[i].size()), blanks.data());
    std::strcat((char *) charcmor, tempblanks);
  }
  int cmf = cmor_axis(new_axis_id(axis_ids), (char *) chardim.c_str(), (char *) "", numchar, (void *) charcmor, 'c', NULL, maxlen,
                      NULL);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
  std::free(charcmor);
  if (Options::cdoVerbose) cdo_print("Successfully registered a character axis.");
}

static void
register_fourth_axis(KVList *kvl, int vlistID, int varID, std::string const &varname, int *axis_ids, std::string const &project_id,
                     int miptab_freq, int *mergeIDs)
{
  (void) project_id;
  (void) miptab_freq;
  std::string fa = get_txtatt(vlistID, varID, "merge_axis");
  std::string chardimatt = kv_get_a_val(kvl, "ca", "");
  std::string chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(chardim, chardimatt, "character_axis", "notSet");
  if (chardim == "notSet" && !fa.empty())
    cdo_abort(
        "ERROR (infile: '%s')! You specify variables to merge to an axis but the axis name and its values are missing!"
        "\n          Specify its name via variable attribute 'character_axis' and its values via 'char_axis_$name' in infofile.",
        cdo_get_stream_name(0));
  if ((chardim != "notSet") && kv_get_a_val(kvl, "capr", "n") == "n")
  {
    int nofvals = 0, nofbnds = 0;
    std::vector<std::string> fvalss, fbndss;
    std::string funits;
    get_charvals_and_bnds(kvl, chardim, fvalss, fbndss, funits, &nofvals, &nofbnds, varname);
    if (nofvals)
    {
      if (!fa.empty())
      {
        if (Options::cdoVerbose) cdo_print("Try to merge variables to a fourth axis.");
        int nvalues = 0;
        std::vector<std::string> idss = parse_string_to_values(kv_get_a_val(kvl, "workfile4err", ""), fa, &nvalues, "fourth");

        if (nvalues > 150)
          cdo_abort("ERROR (infile: '%s')! Only supported for a maximum of 150 variables.", cdo_get_stream_name(0));
        for (int i = 0; i < nvalues; ++i)
        {
          if (std::sscanf(idss[i].c_str(), "%i", &mergeIDs[i]) == EOF)
            cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
        }
        mergeIDs[nvalues] = CMOR_UNDEFID;
        if (nofvals != nvalues)
          cdo_abort("ERROR (infile: '%s')! Number of variables '%d' to merge to a fourth axis '%s' is not equal to specified "
                    "values for this axis: '%d'.",
                    cdo_get_stream_name(0), nvalues, chardim, nofvals);
      }
      else if (nofvals == zaxisInqSize(vlistInqVarZaxis(vlistID, varID)))
      {
        if (Options::cdoVerbose) cdo_print("Try to replace the vertical axis with a character_axis.");
      }
      else if (nofvals == 1)
      {
        if (Options::cdoVerbose) cdo_print("Try to register a forth axis with one value.");
      }
      else
      {
        cdo_abort("ERROR (infile: '%s')! No data to substitute found in infile for specified number of values '%d' for axis '%s'.",
                  cdo_get_stream_name(0), nofvals, chardim);
      }
      double *fbnds = NULL;
      bool daxis = false;
      int i = 0;
      if (fbndss.size())
      {
        fbnds = (double *) std::malloc(nofbnds * sizeof(double));
        for (i = 0; i < nofbnds; ++i)
        {
          int scanreturn = std::sscanf(fbndss[i].c_str(), "%lf", &fbnds[i]);
          if (scanreturn == EOF || scanreturn == 0)
          {
            cdo_warning("Could not convert the '%d'th bnds value of the fourth spatial axis which is to merge into a "
                        "double format.",
                        i);
            break;
          }
        }
      }
      if (i == nofbnds && i > 0) daxis = true;
      double *fvals = NULL;
      i = 0;
      if (nofvals > 0)
      {
        fvals = (double *) std::malloc(nofvals * sizeof(double));
        for (i = 0; i < nofvals; ++i)
        {
          int scanreturn = std::sscanf(fvalss[i].c_str(), "%lf", &fvals[i]);
          if (scanreturn == EOF || scanreturn == 0)
          {
            if (Options::cdoVerbose)
              cdo_print("Could not convert the '%d'th value of the fourth spatial axis which is to merge into a double format.", i);
            break;
          }
        }
      }
      if (i != nofvals && daxis)
        cdo_abort("ERROR (infile: '%s')! Could convert bnds of '%s' to double values but failed to convert values of '%s'."
                  "Check axis attribute 'char_axis_%s'.",
                  cdo_get_stream_name(0), chardim, chardim, chardim);
      if (i == nofvals)
      {
        if (cmor_axis(new_axis_id(axis_ids), (char *) chardim.c_str(), (char *) funits.c_str(), nofvals, (void *) fvals, 'd', fbnds,
                      2, NULL)
            != 0)
          cdo_abort("ERROR (infile: '%s')! Could not register axis '%s' with double values.", cdo_get_stream_name(0), chardim);
      }
      else
      {
        register_char_axis(nofvals, fvalss, axis_ids, chardim);
      }
    }
  }
  else if (Options::cdoVerbose)
    cdo_print("Fourth axis is '%s'.", chardim);
}

static int
getRegisteredPsid(struct mapping vars[], int ps_index)
{
  int psID = CDI_UNDEFID;
  int i = 0;
  for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
  {
    if (vars[i].cdi_varID == ps_index) psID = i;
  }
  return psID;
}

static int
registerPsid(struct mapping vars[], int psindex, int vlistID)
{
  struct mapping *var = new_var_mapping(vars);
  size_t gridsize = cdo_vlist_gridsizemax(vlistID);
  var->cdi_varID = psindex;
  var->cmor_varID = -5;
  cdiDefKeyString(vlistID, psindex, CDI_KEY_NAME, parameter_to_word("ps"));
  if (vlistInqVarDatatype(vlistID, psindex) == CDI_DATATYPE_FLT64)
  {
    var->datatype = 'd';
    var->data = std::malloc(gridsize * sizeof(double));
  }
  else
  {
    var->datatype = 'f';
    var->data = std::malloc(gridsize * sizeof(float));
  }
  int psID = getRegisteredPsid(vars, psindex);
  if (Options::cdoVerbose) cdo_print("9. Successfully registered surface pressure '%d'.", psID);

  return psID;
}

static void
register_z_axis(KVList *kvl, int vlistID, int varID, int zaxisID, std::string const &varname, std::string &zaxis, int *axis_ids,
                int *zfactor_id, std::string const &project_id, int miptab_freq, int *mergeIDs, struct mapping vars[],
                CdoStreamID streamID)
{
  char zaxisunits[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zaxisunits, &length);
  *zfactor_id = 0;
  int cmf = 0;
  int zsize = zaxisInqSize(zaxisID);
  double *levels;

  std::string chardimatt = kv_get_a_val(kvl, "ca", "");
  std::string chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(chardim, chardimatt, "character_axis", "notSet");
  std::string capr = kv_get_a_val(kvl, "capr", "n");
  if (capr != "n") chardim = "notSet";

  /*
    if (strcmp(chardim, "vegtype") == 0 || strcmp(chardim, "landUse") == 0 || strcmp(chardim, "soilpools") == 0 ) */

  if (zsize > 1)
  {
    levels = (double *) std::malloc(zsize * sizeof(double));
    zaxisInqLevels(zaxisID, levels);
    double *zcell_bounds;
    zcell_bounds = (double *) std::malloc(2 * zsize * sizeof(double));
    bool selfGenerated = get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
    char zaxisname[CDI_MAX_NAME];
    length = CDI_MAX_NAME;
    cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zaxisname, &length);
    if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
    {
      if (zaxisunits[0] == ' ') std::strcpy(zaxisunits, "Pa");
      if (zaxis == "alevel")
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                        nullptr);
      else if (zaxis != "notSet")
      {
        if (zaxis.substr(0, 5) == "plev7")
          cmf = cmor_axis(new_axis_id(axis_ids), (char *) zaxis.c_str(), (char *) zaxisunits, zsize, (void *) levels, 'd',
                          zcell_bounds, 2, nullptr);
        else
          cmf = cmor_axis(new_axis_id(axis_ids), (char *) zaxis.c_str(), (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr,
                          0, nullptr);
      }
      else if (project_id != "CMIP5" && project_id != "CMIP6")
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                        nullptr);
      else
      {
        if (project_id == "CMIP6")
        {
          switch (miptab_freq)
          {
            case 4:
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev8", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
            default:
              cdo_warning("CMIP6 requires a zaxis name for zaxis type PRESSURE.\n          The operator "
                          "tries to use a zaxis name matching the number of levels found in infile.");
              char zaxisname2[CDI_MAX_NAME];
              std::snprintf(zaxisname2, CDI_MAX_NAME, "plev%d", zsize);
              cmf = cmor_axis(new_axis_id(axis_ids), zaxisname2, (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
          }
        }
        else
        {
          switch (miptab_freq)
          {
            case 3:
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev7", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
            case 4:
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev8", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
            case 5:
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plev3", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
            default:
              cmf = cmor_axis(new_axis_id(axis_ids), (char *) "plevs", (char *) zaxisunits, zsize, (void *) levels, 'd', nullptr, 0,
                              nullptr);
              break;
          }
        }
      }
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF || (zaxis == "hybrid_height"))
    {
      std::string time_method = get_txtatt(vlistID, varID, "cell_methods");
      std::string att_time_method = kv_get_a_val(kvl, "cm", "");
      check_compare_set(time_method, att_time_method, "cell_methods", " ");
      int zaxistype = zaxisInqType(zaxisID);
      int vctsize = zaxisInqVctSize(zaxisID);
      if (2 * zsize == vctsize) zaxistype = ZAXIS_HYBRID_HALF;
      double *alev_val, *alev_bnds = nullptr, *ap_val, *ap_bnds = nullptr, *b_val, *b_bnds = nullptr;
      double *p0 = (double *) std::malloc(sizeof(double));
      p0[0] = 101325.0;
      if (zaxis == "hybrid_height") p0[0] = 1;

      if (zaxistype == ZAXIS_HYBRID)
      {
        alev_val = (double *) std::malloc(zsize * sizeof(double));
        alev_bnds = (double *) std::malloc((zsize + 1) * sizeof(double));
        ap_val = (double *) std::malloc(zsize * sizeof(double));
        ap_bnds = (double *) std::malloc((zsize + 1) * sizeof(double));
        b_val = (double *) std::malloc(zsize * sizeof(double));
        b_bnds = (double *) std::malloc((zsize + 1) * sizeof(double));
      }
      else
      {
        alev_val = (double *) std::malloc(zsize * sizeof(double));
        ap_val = (double *) std::malloc(zsize * sizeof(double));
        b_val = (double *) std::malloc(zsize * sizeof(double));
      }

      std::string mtproof = kv_get_a_val(kvl, "mtproof", "");
      if (!mtproof.empty())
      {
        if (Options::cdoVerbose) cdo_print("Mapping table: '%s' is applied for ps or orog. ", mtproof);
        /*                  kv_insert_vals(kvl, "capr", (std::string )"y", 1); */
        int filetype = cdo_inq_filetype(streamID);
        kv_insert_vals(kvl, "workfile4err", mtproof, true, false);
        PMList pml = cdo_parse_cmor_file(mtproof, true);
        if (zaxis != "hybrid_height")
        {
          std::vector<std::string> tempo = { "ps" };
          maptab_via_cn(mtproof, pml, tempo, vlistID, 1, kv_get_a_val(kvl, "miptab_freq", ""), filetype, nullptr, false);
        }
      }
      int psindex = CDI_UNDEFID;
      if (zaxis != "hybrid_height")
      {
        psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
        if (psindex == CDI_UNDEFID) psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "code", "134");
        if (psindex == CDI_UNDEFID)
          cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Could not find a surface pressure "
                    "variable in infile. Cannot register a hybrid zaxis without surface pressure.",
                    cdo_get_stream_name(0));
      }
      int psID = getRegisteredPsid(vars, psindex);
      if (zaxis != "hybrid_height")
      {
        if (psID == CDI_UNDEFID)
        {
          {
            if (Options::cdoVerbose) cdo_print("Start to register auxiliary ps. ");
            psID = registerPsid(vars, psindex, vlistID);
            if (Options::cdoVerbose) cdo_print("Successfully registered auxiliary ps. ");
          }
        }
      }
      void *orogdata = nullptr;
      char orogtype = 'f';

      if (zaxis == "hybrid_height")
      {
        get_zhybrid(zaxisID, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "hybrid_height", (char *) "m", zsize, (void *) alev_val, 'd', alev_bnds, 1,
                        nullptr);

        int tsID = 0;

        auto streamID2 = streamOpenRead(cdo_get_stream_name(0));
        int vlistID2 = streamInqVlist(streamID2);
        psindex = getVarIDToMap(vlistID2, vlistNvars(vlistID2), "name", "orog");
        if (vlistInqVarDatatype(vlistID2, psindex) == CDI_DATATYPE_FLT64) orogtype = 'd';

        int gridID = vlistInqVarGrid(vlistID2, psindex);
        auto gridsize = gridInqSize(gridID);
        Varray<double> buffer(gridsize);

        if (vlistInqVarDatatype(vlistID2, psindex) == CDI_DATATYPE_FLT64)
          orogdata = std::malloc(gridsize * sizeof(double));
        else
          orogdata = std::malloc(gridsize * sizeof(float));

        while (true)
        {
          auto numFields = streamInqTimestep(streamID2, tsID);
          if (numFields == 0) break;

          while (numFields--)
          {
            int varIDrw, levelIDrw;
            size_t nmiss;
            streamInqField(streamID2, &varIDrw, &levelIDrw);
            if (varIDrw == psindex)
            {
              cdo_print("read the record");
              streamReadField(streamID2, buffer.data(), &nmiss);
              for (size_t i = 0; i < gridsize; ++i)
              {
                if (vlistInqVarDatatype(vlistID2, psindex) == CDI_DATATYPE_FLT64)
                  ((double *) orogdata)[i] = (double) buffer.data()[i];
                else
                  ((float *) orogdata)[i] = (float) buffer.data()[i];
              }
            }
          }
          tsID++;
        }
        streamClose(streamID2);
      }
      else if (zaxistype == ZAXIS_HYBRID)
      {
        get_zhybrid(zaxisID, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alternate_hybrid_sigma", (char *) "", zsize, (void *) alev_val, 'd',
                        alev_bnds, 1, nullptr);
      }
      else
      {
        get_zhybrid_half(zaxisID, p0, alev_val, b_val, ap_val);
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alternate_hybrid_sigma_half", (char *) "", zsize, (void *) alev_val, 'd',
                        nullptr, 1, nullptr);
      }
      /*cmor_zfactor (int *zfactor_id,int zaxis_id, char *zfactor_name, char *units, int ndims, int axis_ids[],
       * char type, void *zfactor_values, void *zfactor_bounds)*/

      int lev_id = axis_ids[count_axis_ids(axis_ids) - 1];
      int lev_id_array[2];
      lev_id_array[0] = lev_id;
      std::vector<int> hharray(count_axis_ids(axis_ids) - 2);
      for (int i = 0; i < count_axis_ids(axis_ids) - 2; i++) hharray[i] = axis_ids[i + 1];

      if (zaxistype == ZAXIS_HYBRID)
      {
        if (zaxis == "hybrid_height")
        {
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "a", (char *) "", 1, &lev_id_array[0], 'd', (void *) ap_val,
                             (void *) ap_bnds);
        }
        else
        {
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "p0", (char *) "Pa", 0, 0, 'd', (void *) p0, nullptr);
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ap", (char *) "Pa", 1, &lev_id_array[0], 'd', (void *) ap_val,
                             (void *) ap_bnds);
        }
        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "b", (char *) "", 1, &lev_id_array[0], 'd', (void *) b_val,
                           (void *) b_bnds);
        if (zaxis == "hybrid_height")
        {
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "orog", (char *) "m", count_axis_ids(axis_ids) - 2, hharray.data(),
                             orogtype, (void *) orogdata, nullptr);
        }
        else if (time_method[0] == 'p')
        {
          if (count_axis_ids(axis_ids) == 3)
            cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps3", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                               vars[psID].datatype, nullptr, nullptr);
          else
            cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps1", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                               vars[psID].datatype, nullptr, nullptr);
        }
        else if (time_method[0] == 'c')
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps2", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                             vars[psID].datatype, nullptr, nullptr);
        else
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                             vars[psID].datatype, nullptr, nullptr);
      }
      else
      {
        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "b_half", (char *) "", 1, &lev_id_array[0], 'd', (void *) b_val, nullptr);
        cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ap_half", (char *) "Pa", 1, &lev_id_array[0], 'd', (void *) ap_val,
                           nullptr);
        if (time_method[0] == 'p')
        {
          if (count_axis_ids(axis_ids) == 3)
            cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps3", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                               vars[psID].datatype, nullptr, nullptr);
          else
            cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps1", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                               vars[psID].datatype, nullptr, nullptr);
        }
        else if (time_method[0] == 'c')
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps2", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                             vars[psID].datatype, nullptr, nullptr);
        else
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *) "ps", (char *) "Pa", count_axis_ids(axis_ids) - 1, axis_ids,
                             vars[psID].datatype, nullptr, nullptr);
      }
      std::free(alev_val);
      std::free(ap_val);
      std::free(b_val);
      if (zaxistype == ZAXIS_HYBRID)
      {
        if (zaxis == "hybrid_height") std::free(orogdata);
        std::free(ap_bnds);
        std::free(alev_bnds);
        std::free(b_bnds);
      }
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA || zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_LAND)
    {
      if (!zaxisunits[0] || zaxisunits[0] == ' ' || zaxisunits[0] == '\0')
      {
        if (zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA)
          std::strcpy(zaxisunits, "m");
        else
          std::strcpy(zaxisunits, "cm");
      }
      if (selfGenerated)
      {
        zcell_bounds[0] = (double) 0;
        zcell_bounds[2 * zsize - 1] = levels[zsize - 1];
      }
      char longname[CDI_MAX_NAME];
      length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname, &length);
      if (std::string(longname) == "depth_below_sea" || zaxis == "olevel" || std::string(longname) == "ocean depth coordinate")
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "depth_coord", (char *) zaxisunits, zsize, (void *) levels, 'd',
                        zcell_bounds, 2, nullptr);
      else
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "sdepth", (char *) zaxisunits, zsize, (void *) levels, 'd', zcell_bounds, 2,
                        nullptr);
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_ALTITUDE)
    {
      if (zaxis == "alevel")
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "alt", (char *) zaxisunits, zsize, (void *) levels, 'd', zcell_bounds, 2,
                        nullptr);
      else if (zaxis.substr(0, 3) == "alt")
      {
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) zaxis.c_str(), (char *) zaxisunits, zsize, (void *) levels, 'd',
                        zcell_bounds, 2, nullptr);
      }
      else
      {
        char zaxisname2[CDI_MAX_NAME];
        std::snprintf(zaxisname2, CDI_MAX_NAME, "alt%d", zsize);
        cmf = cmor_axis(new_axis_id(axis_ids), zaxisname2, (char *) zaxisunits, zsize, (void *) levels, 'd', zcell_bounds, 2,
                        nullptr);
      }
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_GENERIC && (chardim != "notSet"))
    {
      register_fourth_axis(kvl, vlistID, varID, varname, axis_ids, project_id, miptab_freq, mergeIDs);
      kv_insert_vals(kvl, "capr", "y", true, false);
    }
    else if (zaxisInqType(zaxisID) == ZAXIS_GENERIC || zaxisInqType(zaxisID) == ZAXIS_HEIGHT)
    {
      char longname[CDI_MAX_NAME];
      length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname, &length);

      if (std::string(zaxisname) == "rho")
      {
        if (std::string(zaxisunits) != "kg m-3")
        {
          cdo_abort("ERROR (infile: '%s')! For zaxis with name 'rho' the units must be kg m-3 but are: '%s'",
                    cdo_get_stream_name(0), zaxisunits);
        }
        else
        {
          cmf = cmor_axis(new_axis_id(axis_ids), (char *) "rho", (char *) "kg m-3", zsize, (void *) levels, 'd', zcell_bounds, 2,
                          nullptr);
        }
      }
      else if (std::strcmp(longname, "Z-coordinate in Cartesian system") == 0)
      {
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "height_coord", (char *) zaxisunits, zsize, (void *) levels, 'd',
                        zcell_bounds, 2, nullptr);
      }

      else
        cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Z-axis type %d with name '%s' not yet "
                  "enabled.",
                  cdo_get_stream_name(0), zaxisInqType(zaxisID), zaxisname);
    }
    else
      cdo_abort("ERROR (infile: '%s')! In registration of a vertical axis:\n          Invalid Z-axis type %d . ",
                cdo_get_stream_name(0), zaxisInqType(zaxisID));
    std::free(zcell_bounds);
    std::free(levels);
  }

  else if (zsize == 1 && (zaxis != "notSet"))
  {
    /*      if ( strcmp(zaxis, "notSet") == 0 )
            {
              if (
            }
          else
            {*/
    std::string szc_value = kv_get_a_val(kvl, zaxis, "");
    if (!szc_value.empty())
    {
      levels = (double *) std::malloc(sizeof(double));
      levels[0] = (double) atof(szc_value.c_str());

      std::string szc_key;
      szc_key = zaxis + "_bounds";

      int numchar = 0;
      std::vector<std::string> szc_bndss = kv_get_vals(kvl, szc_key, &numchar);
      if (numchar != 0 && numchar != 2)
        cdo_warning("Scalar z coordinate bounds need to be exactly two values! You configured '%d'.", numchar);

      double szc_bnds[2];
      if (!szc_bndss.empty())
      {
        if (std::sscanf(szc_bndss[0].c_str(), "%lf", &szc_bnds[0]) == EOF)
          cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
        if (std::sscanf(szc_bndss[1].c_str(), "%lf", &szc_bnds[1]) == EOF)
          cdo_abort("ERROR (infile: '%s')! Internal error.", cdo_get_stream_name(0));
      }

      szc_key = zaxis + "_units";
      std::string szcunits = kv_get_a_val(kvl, szc_key, "m");

      if (Options::cdoVerbose)
        cdo_print("Scalar z coordinate name is: '%s'\n          Scalar z coordinate value is: '%f' meter", zaxis, levels[0]);
      if (!szc_bndss.empty())
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) zaxis.c_str(), (char *) szcunits.c_str(), zsize, (void *) levels, 'd',
                        szc_bnds, 2, nullptr);
      else
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) zaxis.c_str(), (char *) szcunits.c_str(), zsize, (void *) levels, 'd',
                        nullptr, 0, nullptr);
      std::free(levels);
    }
    else
      cdo_print("You specified z_axis='%s'.\n          No value has been specified, axis will be created with the "
                "default value.",
                zaxis);
  }
  else
    cdo_print("Vertical axis is either default and scalar or not available.");

  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
}

/*
static void register_character_dimension(int *axis_ids, std::string filename)
{
  printf("The grid type is generic and a dimension 'basin' is found.\nTherefore, it is tried to read the character
dimension.\n");
  int nc_file_id, nfiledims, nvars, ngatts, unlimdimid;
  nc_type xtypep;
  int varndims, varnattsp;
  int *vardimids;

  std::string varname = std::malloc(36);
  std::string dimname = std::malloc(36);

  size_t dimlength, dimstrlength;

  nc_open(filename, NC_NOWRITE, &nc_file_id);
  nc_inq(nc_file_id, &nfiledims, &nvars, &ngatts, &unlimdimid);
  vardimids = std::malloc(nfiledims * sizeof(int));
  void *final_chardim;
  for ( int i = 0; i < nvars; i++ )
    {
      nc_inq_var(nc_file_id, i, varname, &xtypep, &varndims, vardimids, &varnattsp);
      if ( strcmp(varname, "region") == 0 )
        {
          nc_inq_dim(nc_file_id, vardimids[1], dimname, &dimstrlength);
          nc_inq_dim(nc_file_id, vardimids[0], dimname, &dimlength);

          final_chardim = (void *)std::malloc(dimstrlength * dimlength *sizeof(char));
          nc_get_var(nc_file_id, i, final_chardim);
        }
    }
  nc_close(nc_file_id);
  cmor_axis(new_axis_id(axis_ids), dimname, "", dimlength, final_chardim, 'c',  nullptr, dimstrlength, nullptr);
  std::free(varname);
  std::free(dimname);
  std::free(vardimids);
}
*/
static void
change_zaxis(KVList *kvl, int vlistID, int zaxisID, int zaxisID2, std::string const &grid_file)
{
  int a, b;
  a = zaxisInqSize(zaxisID);
  b = zaxisInqSize(zaxisID2);

  bool noZaxis = false;
  if (zaxisInqType(zaxisID2) == ZAXIS_SURFACE && b == 1)
  {
    double level[1];
    zaxisInqLevels(zaxisID2, level);
    if (level[0] < 1) noZaxis = true;
  }

  if (!noZaxis)
  {
    if (a != b)
    {
      cdo_warning("Could not use zaxis from file '%s' configured via attribute 'ginfo'\n          because total "
                  "size of infile: '%d' is not identical to total size of ginfo file: '%d'.",
                  grid_file, a, b);
    }
    else if (zaxisID != zaxisID2)
    {
      vlistChangeZaxis(vlistID, zaxisID, zaxisID2);
      char zaxisname[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zaxisname, &length);
      kv_insert_vals(kvl, "za", zaxisname, false, false);
      cdo_print("Successfully substituted zaxis.");
    }
    else
      cdo_print("Zaxis from grid info file %s is equal to infile Zaxis.", grid_file);
  }
}

static void
change_grid(int vlistID, int gridID, int gridID2, std::string grid_file)
{
  if (!gridID2)
    cdo_abort(
        "ERROR (infile: '%s')! Could not use grid from file '%s' configured via attribute 'ginfo'\n          because of internal "
        "problems.",
        cdo_get_stream_name(0), grid_file);

  int a, b;
  bool lswitch = true;
  a = gridInqSize(gridID);
  b = gridInqSize(gridID2);

  if (b == 1)
    cdo_print("Grid is not substituted because the size of the grid found in the grid info file is one.");
  else
  {
    if (a != b)
    {
      lswitch = false;
      cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because total "
                  "size of infile: '%d' is not identical to total size of ginfo file: '%d'.",
                  grid_file, a, b);
    }

    a = gridInqYsize(gridID);
    b = gridInqYsize(gridID2);
    if (a != b)
    {
      if (gridInqType(gridID) != GRID_UNSTRUCTURED) lswitch = false;
      cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because ysize of grid info "
                  "file '%d' does not match ysize of infile '%d'.",
                  grid_file, b, a);
    }

    a = gridInqXsize(gridID);
    b = gridInqXsize(gridID2);
    if (a != b)
    {
      if (gridInqType(gridID) != GRID_UNSTRUCTURED) lswitch = false;
      lswitch = false;
      cdo_warning("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because xsize of grid info "
                  "file '%d' does not match xsize of infile '%d'.",
                  grid_file, b, a);
    }

    if (lswitch && gridID != gridID2)
    {
      vlistChangeGrid(vlistID, gridID, gridID2);
      cdo_print("Successfully substituted grid.");
    }
    else if (lswitch)
      cdo_print("Grid from grid info file %s is equal to infile grid.", grid_file);
  }
}

static void
move_lons(double *xcoord_vals, double *xcell_bounds, int xsize, int xboundsize, int xnbounds)
{
  int testbool = 0;
  for (int i = 0; i < xsize; ++i)
    if (xcoord_vals[i] < 0.0)
    {
      testbool = 1;
      break;
    }
  if (testbool > 0)
    for (int i = 0; i < xsize; ++i)
      if (xcoord_vals[i] < 0) xcoord_vals[i] += 360.0;
  if (xnbounds > 1 && testbool > 0)
    for (int j = 0; j < xboundsize; ++j)
      if (xcell_bounds[j] < 0) xcell_bounds[j] += 360.0;
}

static void
inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals, double *xcell_bounds,
                        double *ycell_bounds)
{
  char unitstring[CDI_MAX_NAME];
  gridInqYvals(gridID, ycoord_vals);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, unitstring, &length);
  if (std::strcmp(unitstring, "radian") == 0)
  {
    for (size_t i = 0; i < gridInqXsize(gridID); i++) ycoord_vals[i] = (180. / M_PI) * ycoord_vals[i];
    for (size_t i = 0; i < gridInqNvertex(gridID) * gridInqSize(gridID); i++) ycell_bounds[i] = (180. / M_PI) * ycell_bounds[i];
  }
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, unitstring, &length);
  if (std::strcmp(unitstring, "radian") == 0)
  {
    for (size_t i = 0; i < gridInqXsize(gridID); i++) xcoord_vals[i] = (180. / M_PI) * xcoord_vals[i];
    for (size_t i = 0; i < gridInqNvertex(gridID) * gridInqSize(gridID); i++) xcell_bounds[i] = (180. / M_PI) * xcell_bounds[i];
  }
}

static void
get_cmor_table(KVList *kvl, std::string const &project_id)
{
  int gridtable_id;
  int cmf = 0;
  char gridtable[CMOR_MAX_STRING];
  std::string mip_table_dir = kv_get_a_val(kvl, "mip_table_dir", "");
  if (!mip_table_dir.empty())
#if (CMOR_VERSION_MAJOR == 2)
  {
    if (mip_table_dir[strlen(mip_table_dir) - 1] == '/')
      std::snprintf(gridtable, CMOR_MAX_STRING, "%s%s_grids", mip_table_dir.c_str(), project_id.c_str());
    else
      std::snprintf(gridtable, CMOR_MAX_STRING, "%s/%s_grids", mip_table_dir.c_str(), project_id.c_str());
  }
#elif (CMOR_VERSION_MAJOR == 3)
  {
    std::snprintf(gridtable, CMOR_MAX_STRING, "%s/%s_grids.json", mip_table_dir.c_str(), project_id.c_str());
  }
#endif
  if (!mip_table_dir.empty())
  {
    if (file_exist(std::string(gridtable), false, "Cmor-grid_table", false))
    {
      cmf = cmor_load_table(gridtable, &gridtable_id);
      cmf = cmor_set_table(gridtable_id);
    }
    else
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          File '%s' not found!\n          "
                "A project grid table is required for this type of grid but not "
                "found in the mip table directory '%s'.",
                cdo_get_stream_name(0), gridtable, mip_table_dir);
  }
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_load_table or cmor_set_table failed!", cdo_get_stream_name(0));
}

static void
check_and_gen_bounds(int gridID, int nbounds, int length, double *coord_vals, double *cell_bounds, int x)
{
  if (nbounds != 2 * length)
  {
    gen_bounds(length, coord_vals, cell_bounds);
    if (x)
    {
      gridDefNvertex(gridID, 2);
      gridDefXbounds(gridID, cell_bounds);
    }
    else
      gridDefYbounds(gridID, cell_bounds);
  }
}

static double
lonbnds_mids_trans_check(double value1, double value2)
{
  if (std::fabs(value1 - value2) < 180.0)
    return (value1 + value2) * 0.5;
  else
  {
    if (value1 + value2 < 360.0)
      return (value1 + value2 + 360.0) * 0.5;
    else
      return (value1 + value2 + 360.0) * 0.5 - 360.0;
  }
}

static double
lonbnds_bnds_trans_check(double value1, double value2)
{
  if (std::fabs(value1 - value2) < 180)
  {
    if (2 * value1 < value2)
      return (2 * value1 - value2 + 360.0);
    else if (2 * value1 > value2 + 360.0)
      return (2 * value1 - value2 - 360.0);
    else
      return (2 * value1 - value2);
  }
  else if (value1 - value2 > 180)
    return (2 * value1 - value2 - 360.0);
  else
    return (2 * value1 - value2 + 360.0);
}

static void
check_and_gen_bounds_curv(int gridID, int totalsize, int xnbounds, int xlength, double *xcoord_vals, double *xcell_bounds,
                          int ynbounds, int ylength, double *ycoord_vals, double *ycell_bounds)
{
  if (xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (is_equal(xcell_bounds[1], 0.0) && is_equal(xcell_bounds[2], 0.0))
      || (is_equal(ycell_bounds[1], 0.0) && is_equal(ycell_bounds[2], 0.0)))
  {
    Varray2D<double> halflons(xlength + 1, Varray<double>(ylength));
    Varray2D<double> halflats(xlength, Varray<double>(ylength + 1));
    Varray2D<double> halflonsOnhalflats(xlength + 1, Varray<double>(ylength + 1));
    Varray2D<double> halflatsOnhalflons(xlength + 1, Varray<double>(ylength + 1));

    /**/
    /*************Half-lons with 360-0 transmission check**************/
    /**/
    for (int j = 0; j < ylength; ++j)
    {
      for (int i = 1; i < xlength; ++i)
        halflons[i][j] = lonbnds_mids_trans_check(xcoord_vals[i - 1 + j * xlength], xcoord_vals[i + j * xlength]);
      /*left and right boundary: */
      halflons[0][j] = lonbnds_bnds_trans_check(xcoord_vals[j * xlength], halflons[1][j]);
      halflons[xlength][j] = lonbnds_bnds_trans_check(xcoord_vals[j * xlength - 1], halflons[xlength - 1][j]);
    }
    /**/
    /*************Half-lats **************/
    /**/
    for (int i = 0; i < xlength; ++i)
    {
      for (int j = 1; j < ylength; ++j) halflats[i][j] = (ycoord_vals[i + (j - 1) * xlength] + ycoord_vals[i + j * xlength]) * 0.5;
      /*upper and lower boundary: */
      halflats[i][0] = 2 * ycoord_vals[i] - halflats[i][1];
      halflats[i][ylength] = 2 * ycoord_vals[i + (ylength - 1) * xlength] - halflats[i][ylength - 1];
    }
    /**/
    /****************Half-lons-on-half-lats with 0-360 transmission check**********/
    /****************Half-lats-on-half-lons                              **********/
    /**/

    for (int i = 1; i < xlength; ++i)
    {
      for (int j = 1; j < ylength; ++j)
      {
        halflonsOnhalflats[i][j] = lonbnds_mids_trans_check(halflons[i][j - 1], halflons[i][j]);
        halflatsOnhalflons[i][j] = (halflats[i - 1][j] + halflats[i][j]) * 0.5;
      }
      /*upper and lower boundary: */
      halflonsOnhalflats[i][0] = lonbnds_bnds_trans_check(halflons[i][0], halflonsOnhalflats[i][1]);
      halflonsOnhalflats[i][ylength] = lonbnds_bnds_trans_check(halflons[i][ylength - 1], halflonsOnhalflats[i][ylength - 1]);
      halflatsOnhalflons[i][0] = (halflats[i - 1][0] + halflats[i][0]) * 0.5;
      halflatsOnhalflons[i][ylength] = (halflats[i - 1][ylength] + halflats[i][ylength]) * 0.5;
    }

    /*left and right boundary: */
    for (int j = 1; j < ylength; ++j)
    {
      halflonsOnhalflats[0][j] = lonbnds_mids_trans_check(halflons[0][j - 1], halflons[0][j]);
      halflonsOnhalflats[xlength][j] = lonbnds_mids_trans_check(halflons[xlength][j - 1], halflons[xlength][j]);

      halflatsOnhalflons[0][j] = 2 * halflats[0][j] - halflatsOnhalflons[1][j];
      halflatsOnhalflons[xlength][j] = 2 * halflats[xlength - 1][j] - halflatsOnhalflons[xlength - 1][j];
    }
    halflatsOnhalflons[0][0] = 2 * halflats[0][0] - halflatsOnhalflons[1][0];
    halflatsOnhalflons[0][ylength] = 2 * halflats[0][ylength] - halflatsOnhalflons[1][ylength];
    halflatsOnhalflons[xlength][0] = 2 * halflats[xlength - 1][0] - halflatsOnhalflons[xlength - 1][0];
    halflatsOnhalflons[xlength][ylength] = 2 * halflats[xlength - 1][ylength] - halflatsOnhalflons[xlength - 1][ylength];

    halflonsOnhalflats[0][0] = lonbnds_bnds_trans_check(halflons[0][0], halflonsOnhalflats[0][1]);
    halflonsOnhalflats[0][ylength] = lonbnds_bnds_trans_check(halflons[0][ylength - 1], halflonsOnhalflats[0][ylength - 1]);
    halflonsOnhalflats[xlength][0] = lonbnds_bnds_trans_check(halflons[xlength][0], halflonsOnhalflats[xlength][1]);
    halflonsOnhalflats[xlength][ylength]
        = lonbnds_bnds_trans_check(halflons[xlength][ylength - 1], halflonsOnhalflats[xlength - 1][ylength]);

    for (int i = 0; i < xlength; ++i)
      for (int j = 0; j < ylength; ++j)
      {
        xcell_bounds[4 * (j * xlength + i)] = halflonsOnhalflats[i][j + 1];
        xcell_bounds[4 * (j * xlength + i) + 1] = halflonsOnhalflats[i][j];
        xcell_bounds[4 * (j * xlength + i) + 2] = halflonsOnhalflats[i + 1][j];
        xcell_bounds[4 * (j * xlength + i) + 3] = halflonsOnhalflats[i + 1][j + 1];
        ycell_bounds[4 * (j * xlength + i)] = halflatsOnhalflons[i][j + 1];
        ycell_bounds[4 * (j * xlength + i) + 1] = halflatsOnhalflons[i][j];
        ycell_bounds[4 * (j * xlength + i) + 2] = halflatsOnhalflons[i + 1][j];
        ycell_bounds[4 * (j * xlength + i) + 3] = halflatsOnhalflons[i + 1][j + 1];
      }
    gridDefNvertex(gridID, 4);
    gridDefXbounds(gridID, xcell_bounds);
    gridDefYbounds(gridID, ycell_bounds);
  }
}
static void
register_lon_axis(int gridID, int xlength, int *axis_ids)
{
  Varray<double> xcoord_vals(xlength);
  if (gridInqXvals(gridID, xcoord_vals.data()))
  {
    Varray<double> xcell_bounds(2 * xlength);
    int xnbounds = gridInqXbounds(gridID, xcell_bounds.data());
    check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals.data(), xcell_bounds.data(), 1);
    int cmf = cmor_axis(new_axis_id(axis_ids), (char *) "longitude", (char *) "degrees_east", xlength, (void *) xcoord_vals.data(),
                        'd', (void *) xcell_bounds.data(), 2, nullptr);
    if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
  }
}

static void
register_lat_axis(int gridID, int ylength, int *axis_ids)
{
  Varray<double> ycoord_vals(ylength);
  if (gridInqYvals(gridID, ycoord_vals.data()))
  {
    Varray<double> ycell_bounds(2 * ylength);
    int ynbounds = gridInqYbounds(gridID, ycell_bounds.data());
    check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals.data(), ycell_bounds.data(), 0);
    int cmf = cmor_axis(new_axis_id(axis_ids), (char *) "latitude", (char *) "degrees_north", ylength, (void *) ycoord_vals.data(),
                        'd', (void *) ycell_bounds.data(), 2, nullptr);
    if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
  }
}

static void
register_projection(int *grid_ids, int projID, double *ycoord_vals, double *xcoord_vals, double *ycell_bounds, double *xcell_bounds,
                    int xlength, int ylength, int projtype)
{
  int cmf = 0;
  int pxnbounds;
  int pynbounds;
  int pylength = gridInqYsize(projID);
  int pxlength = gridInqXsize(projID);
  double *pxcoord_vals = (double *) std::malloc(pxlength * sizeof(double));
  double *pycoord_vals = (double *) std::malloc(pylength * sizeof(double));
  double *pxcell_bounds = (double *) std::malloc(2 * pxlength * sizeof(double));
  double *pycell_bounds = (double *) std::malloc(2 * pylength * sizeof(double));
  inquire_vals_and_bounds(projID, &pxnbounds, &pynbounds, pxcoord_vals, pycoord_vals, pxcell_bounds, pycell_bounds);
  check_and_gen_bounds(projID, pxnbounds, pxlength, pxcoord_vals, pxcell_bounds, 1);
  check_and_gen_bounds(projID, pynbounds, pylength, pycoord_vals, pycell_bounds, 0);

  char p_rll_cmor[CMOR_MAX_STRING];
  int l_p_rll = std::strlen("grid_north_pole_longitude") + 1;
  std::memcpy(p_rll_cmor,
              "grid_north_pole_latitude\0 "
              "grid_north_pole_longitude\0north_pole_grid_longitude\0",
              3 * l_p_rll);

  char u_rll_cmor[CMOR_MAX_STRING];
  int l_u_rll = std::strlen("degrees_north") + 1;
  std::memcpy(u_rll_cmor, "degrees_north\0degrees_east\0 degrees_east\0 ", 3 * l_u_rll);

  char p_lcc_cmor[CMOR_MAX_STRING];
  int l_p_lcc = std::strlen("longitude_of_central_meridian") + 1;
  std::memcpy(p_lcc_cmor,
              "standard_parallel1\0           "
              "longitude_of_central_meridian\0latitude_of_projection_"
              "origin\0standard_parallel2\0           ",
              4 * l_p_lcc);

  char u_cmor[CMOR_MAX_STRING];
  int l_u_cmor = 6;

  std::vector<std::string> p_rll = { "grid_north_pole_latitude", "grid_north_pole_longitude", "north_pole_grid_longitude" };

  std::vector<std::string> p_lcc
      = { "standard_parallel1", "longitude_of_central_meridian", "latitude_of_projection_origin", "standard_parallel2" };
  std::vector<std::string> p_ps_falselsp = { "straight_vertical_longitude_from_pole", "latitude_of_projection_origin",
                                             "scale_factor_at_projection_origin", "false_easting", "false_northing" };
  std::vector<std::string> p_ps_falsestandard = { "straight_vertical_longitude_from_pole", "latitude_of_projection_origin",
                                                  "standard_parallel", "false_easting", "false_northing" };
  std::vector<std::string> p_ps_truestandard
      = { "straight_vertical_longitude_from_pole", "latitude_of_projection_origin", "standard_parallel" };
  std::vector<std::string> p_ps;
  /*char *p_ps = nullptr;*/

  int atttype, attlen;
  char attname[CDI_MAX_NAME];

  int natts;
  cdiInqNatts(projID, CDI_GLOBAL, &natts);

  char p_ps_cmor[CMOR_MAX_STRING];
  int l_p_ps = std::strlen("straight_vertical_longitude_from_pole") + 1;

  if (projtype == CDI_PROJ_STERE)
  {
    bool lsp = false;
    bool lfalse = false;
    for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
      if (std::strcmp(attname, "standard_parallel") == 0)
        lsp = true;
      else if (std::strcmp(attname, "false_easting") == 0)
        lfalse = true;
    }

    if (lfalse == true)
    {
      std::memcpy(u_cmor,
                  "degrees_east\0 "
                  "degrees_north\0"
                  "degrees_north\0"
                  "degrees_east\0 "
                  "degrees_north\0",
                  3 * l_u_rll);
      if (lsp == false)
      {
        p_ps = p_ps_falselsp;
        std::memcpy(p_ps_cmor,
                    "straight_vertical_longitude_from_pole\0"
                    "latitude_of_projection_origin\0        "
                    "scale_factor_at_projection_origin\0    "
                    "false_easting\0                        "
                    "false_northing\0                       ",

                    5 * l_p_ps);
      }
      else
      {
        p_ps = p_ps_falsestandard;
        std::memcpy(p_ps_cmor,
                    "straight_vertical_longitude_from_pole\0"
                    "latitude_of_projection_origin\0        "
                    "standard_parallel\0                    "
                    "false_easting\0                        "
                    "false_northing\0                       ",

                    5 * l_p_ps);
      }
    }
    else
    {
      p_ps = p_ps_truestandard;
      std::memcpy(p_ps_cmor,
                  "straight_vertical_longitude_from_pole\0"
                  "latitude_of_projection_origin\0        "
                  "standard_parallel\0                    ",
                  3 * l_p_ps);
      std::memcpy(u_cmor,
                  "degrees_east\0 "
                  "degrees_north\0"
                  "degrees_north\0",
                  3 * l_u_rll);
    }
  }
  double *parameter_values = nullptr;

  char mapping[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(projID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, mapping, &length);

  int p_len = 0;
  switch (projtype)
  {
    case CDI_PROJ_RLL: p_len = p_rll.size(); break;
    case CDI_PROJ_LAEA:
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          This grid projection is not yet enabled.",
                cdo_get_stream_name(0));
      break;
    case CDI_PROJ_LCC: p_len = p_lcc.size(); break;
    case CDI_PROJ_SINU:
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          This grid projection is not yet enabled.",
                cdo_get_stream_name(0));
      break;
    case CDI_PROJ_STERE: p_len = p_ps.size(); break;
  }
  if (natts < p_len)
    cdo_warning("In grid registration:\n          Number of required grid mapping attributes '%d' is larger than the "
                "number of given grid mapping attributes '%d'.\n          Note that all required mapping attributes are "
                "set to 0.0 by default in case they are not given.",
                p_len, natts);

  parameter_values = (double *) std::malloc(p_len * sizeof(double));
  for (int i = 0; i < p_len; ++i) parameter_values[i] = 0.0;

  for (int iatt = 0; iatt < natts; ++iatt)
  {
    cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
    if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
    {
      if (attlen > 1)
        cdo_abort("ERROR (infile: '%s')! In grid registration:\n          Dont know what to do with grid mapping attribute '%s'.",
                  cdo_get_stream_name(0), attname);
      Varray<double> attflt(attlen);
      cdiInqAttFlt(projID, CDI_GLOBAL, attname, attlen, attflt.data());
      int i = 0;
      for (i = 0; i < p_len; ++i)
      {
        if (projtype == CDI_PROJ_RLL)
        {
          if (p_rll[i] == attname)
          {
            parameter_values[i] = attflt[0];
            break;
          }
        }
        else if (projtype == CDI_PROJ_LCC)
        {
          if (p_lcc[i] == attname)
          {
            parameter_values[i] = attflt[0];
            break;
          }
        }
        else if (projtype == CDI_PROJ_STERE)
        {
          if (p_ps[i] == attname)
          {
            parameter_values[i] = attflt[0];
            break;
          }
        }
      }
      if (i == p_len) cdo_warning("In grid registration:\n          grid mapping attribute '%s' is neglected.", attname);
    }
    else if (atttype == CDI_DATATYPE_TXT)
    {
      std::vector<char> atttxt(attlen + 1);
      cdiInqAttTxt(projID, CDI_GLOBAL, attname, attlen, atttxt.data());
      atttxt[attlen] = 0;
    }
  }

  int grid_axis[2];
  if (projtype == CDI_PROJ_RLL)
  {
    cmf = cmor_axis(&grid_axis[0], (char *) "grid_latitude", (char *) "degrees_north", pylength, (void *) pycoord_vals, 'd',
                    (void *) pycell_bounds, 2, nullptr);
    cmf = cmor_axis(&grid_axis[1], (char *) "grid_longitude", (char *) "degrees_east", pxlength, (void *) pxcoord_vals, 'd',
                    (void *) pxcell_bounds, 2, nullptr);
    cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4, (void *) ycell_bounds,
                    (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
    cmf = cmor_set_grid_mapping(grid_ids[0], "rotated_latitude_longitude", p_len, (char **) p_rll_cmor, l_p_rll, parameter_values,
                                (char **) u_rll_cmor, l_u_rll);
#elif (CMOR_VERSION_MAJOR == 3)
    cmf = cmor_set_grid_mapping(grid_ids[0], (char *) "rotated_latitude_longitude", p_len, p_rll_cmor, l_p_rll, parameter_values,
                                u_rll_cmor, l_u_rll);
#endif
  }
  else if (projtype == CDI_PROJ_LCC)
  {
    std::memcpy(u_cmor, "      \0      \0      \0      \0", 4 * l_u_cmor);
    double *xii = (double *) std::malloc(xlength * sizeof(double));
    double *yii = (double *) std::malloc(ylength * sizeof(double));
    for (int i = 0; i < xlength; ++i) xii[i] = (double) i;
    for (int i = 0; i < ylength; ++i) yii[i] = (double) i;
    cmf = cmor_axis(&grid_axis[0], (char *) "x", (char *) "m", ylength, (void *) yii, 'd', 0, 0, nullptr);
    cmf = cmor_axis(&grid_axis[1], (char *) "y", (char *) "m", xlength, (void *) xii, 'd', 0, 0, nullptr);
    cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4, (void *) ycell_bounds,
                    (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
    cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len, (char **) p_lcc_cmor, l_p_lcc, parameter_values, (char **) u_cmor,
                                l_u_cmor);
#elif (CMOR_VERSION_MAJOR == 3)
    cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len, p_lcc_cmor, l_p_lcc, parameter_values, u_cmor, l_u_cmor);
#endif
    std::free(xii);
    std::free(yii);
  }
  else if (projtype == CDI_PROJ_STERE)
  {
    cmf = cmor_axis(&grid_axis[0], (char *) "y", (char *) "m", ylength, (void *) pycoord_vals, 'd', (void *) pycell_bounds, 2,
                    nullptr);
    cmf = cmor_axis(&grid_axis[1], (char *) "x", (char *) "m", xlength, (void *) pxcoord_vals, 'd', (void *) pxcell_bounds, 2,
                    nullptr);
    cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4, (void *) ycell_bounds,
                    (void *) xcell_bounds);
#if (CMOR_VERSION_MAJOR == 2)
    cmf = cmor_set_grid_mapping(grid_ids[0], "polar_stereographic", p_len, (char **) p_ps_cmor, l_p_ps, parameter_values,
                                (char **) u_ps_cmor, l_u_ps);
#elif (CMOR_VERSION_MAJOR == 3)
    cmf = cmor_set_grid_mapping(grid_ids[0], (char *) "polar_stereographic", p_len, p_ps_cmor, l_p_ps, parameter_values, u_cmor,
                                l_u_rll);
#endif
  }

  std::free(parameter_values);
  std::free(pxcell_bounds);
  std::free(pycell_bounds);
  std::free(pxcoord_vals);
  std::free(pycoord_vals);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis or cmor_set_grid_mapping failed!", cdo_get_stream_name(0));
}

static void
register_grid(KVList *kvl, int vlistID, int varID, int *axis_ids, int *grid_ids, std::string const &project_id,
              std::string const &cmor_name)
{
  int cmf = 0;
  int gridID = vlistInqVarGrid(vlistID, varID);

  std::string chardimatt = kv_get_a_val(kvl, "ca", "");
  std::string movelons = kv_get_a_val(kvl, "ml", "y");
  std::string chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(chardim, chardimatt, "character_axis", "notSet");

  char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, xname, &length);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, yname, &length);

  char xdimname[CDI_MAX_NAME], ydimname[CDI_MAX_NAME];
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_DIMNAME, xdimname, &length);
  length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_DIMNAME, ydimname, &length);

  auto totalsize = gridInqSize(gridID);

  if (totalsize > 1)
  {
    int projID = gridInqProj(gridID);
    int projtype = CDI_UNDEFID;
    if (projID != CDI_UNDEFID) projtype = gridInqProjType(projID);
    int type = gridInqType(gridID);
    int ylength = gridInqYsize(gridID);
    int xlength = gridInqXsize(gridID);
    double *xcoord_vals = nullptr;
    double *ycoord_vals = nullptr;
    double *xcell_bounds = nullptr;
    double *ycell_bounds = nullptr;
    int xnbounds;
    int ynbounds;

    if (type == GRID_GAUSSIAN || type == GRID_LONLAT)
    {
      grid_ids[0] = 0;
      xcoord_vals = (double *) std::malloc(xlength * sizeof(double));
      ycoord_vals = (double *) std::malloc(ylength * sizeof(double));
      xcell_bounds = (double *) std::malloc(2 * xlength * sizeof(double));
      ycell_bounds = (double *) std::malloc(2 * ylength * sizeof(double));
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);

      check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
      check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
      if (ylength > 1)
      {
        if (ycell_bounds[0] < -90) ycell_bounds[0] = -90;
        if (ycell_bounds[2 * ylength - 1] > 90) ycell_bounds[2 * ylength - 1] = 90;
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "latitude", (char *) "degrees_north", ylength, (void *) ycoord_vals, 'd',
                        (void *) ycell_bounds, 2, nullptr);
      }
      if (xlength > 1)
        cmf = cmor_axis(new_axis_id(axis_ids), (char *) "longitude", (char *) "degrees_east", xlength, (void *) xcoord_vals, 'd',
                        (void *) xcell_bounds, 2, nullptr);

      std::free(xcell_bounds);
      std::free(ycell_bounds);
      std::free(xcoord_vals);
      std::free(ycoord_vals);
    }
    else if (type == GRID_UNSTRUCTURED)
    {
      int nvertex = gridInqNvertex(gridID);
      xcoord_vals = (double *) std::malloc(totalsize * sizeof(double));
      ycoord_vals = (double *) std::malloc(totalsize * sizeof(double));
      /* maximal 4 gridbounds per gridcell permitted */
      if (nvertex)
      {
        xcell_bounds = (double *) std::malloc(nvertex * totalsize * sizeof(double));
        ycell_bounds = (double *) std::malloc(nvertex * totalsize * sizeof(double));
      }
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
      /* In a projection, this is done by setting mapping parameter */
      if (movelons == "y") move_lons(xcoord_vals, xcell_bounds, totalsize, nvertex * totalsize, xnbounds);
      int grid_axis[2];
      double *coord_vals;
      coord_vals = (double *) std::malloc(xlength * sizeof(double));
      for (int j = 0; j < xlength; ++j) coord_vals[j] = (double) j;
      if (chardim == "site")
      {
        cmf = cmor_axis(&grid_axis[0], (char *) "site", (char *) "1", xlength, (void *) coord_vals, 'd', 0, 0, nullptr);
        get_cmor_table(kvl, project_id);
        cmf = cmor_grid(&grid_ids[0], 1, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 0, nullptr, nullptr);
        kv_insert_vals(kvl, "capr", (char *) "y", true, false);
      }
      else
      {
        get_cmor_table(kvl, project_id);
        cmf = cmor_axis(&grid_axis[0], (char *) "i_index", (char *) "1", xlength, (void *) coord_vals, 'd', 0, 0, nullptr);
        cmf = cmor_grid(&grid_ids[0], 1, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, nvertex, (void *) ycell_bounds,
                        (void *) xcell_bounds);
      }
      std::free(coord_vals);
      std::free(xcoord_vals);
      std::free(ycoord_vals);
      if (xcell_bounds) std::free(xcell_bounds);
      if (ycell_bounds) std::free(ycell_bounds);
    }
    else if (type == GRID_CURVILINEAR)
    {
      xcoord_vals = (double *) std::malloc(totalsize * sizeof(double));
      ycoord_vals = (double *) std::malloc(totalsize * sizeof(double));
      /* maximal 4 gridbounds per gridcell permitted */
      xcell_bounds = (double *) std::malloc(4 * totalsize * sizeof(double));
      ycell_bounds = (double *) std::malloc(4 * totalsize * sizeof(double));
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
      /* In a projection, this is done by setting mapping parameter */
      if (movelons == "y") move_lons(xcoord_vals, xcell_bounds, totalsize, 4 * totalsize, xnbounds);
      get_cmor_table(kvl, project_id);
      int grid_axis[2];
      check_and_gen_bounds_curv(gridID, totalsize, xnbounds, xlength, xcoord_vals, xcell_bounds, ynbounds, ylength, ycoord_vals,
                                ycell_bounds);
      if (projID == CDI_UNDEFID || projtype == CDI_UNDEFID)
      {
        double *xncoord_vals;
        double *yncoord_vals;
        xncoord_vals = (double *) std::malloc(xlength * sizeof(double));
        yncoord_vals = (double *) std::malloc(ylength * sizeof(double));
        for (int j = 0; j < ylength; ++j) yncoord_vals[j] = (double) j;
        for (int j = 0; j < xlength; ++j) xncoord_vals[j] = (double) j;
        cmf = cmor_axis(&grid_axis[0], (char *) "j_index", (char *) "1", ylength, (void *) yncoord_vals, 'd', 0, 0, nullptr);
        cmf = cmor_axis(&grid_axis[1], (char *) "i_index", (char *) "1", xlength, (void *) xncoord_vals, 'd', 0, 0, nullptr);
        cmf = cmor_grid(&grid_ids[0], 2, grid_axis, 'd', (void *) ycoord_vals, (void *) xcoord_vals, 4, (void *) ycell_bounds,
                        (void *) xcell_bounds);
        std::free(xncoord_vals);
        std::free(yncoord_vals);
        std::free(xcoord_vals);
        std::free(ycoord_vals);
        std::free(xcell_bounds);
        std::free(ycell_bounds);
      }
      /*else
        {
          cmf = cmor_axis(&grid_axis[0],    "grid_longitude",   "degrees",    xlength,    (void *)xcoord_vals,
        'd', 0, 0, nullptr);
          cmf = cmor_axis(&grid_axis[1],    "grid_latitude",    "degrees",    ylength,    (void *)ycoord_vals,
        'd', 0, 0, nullptr);
          cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,
        2,     (void *)ycell_bounds,    (void *)xcell_bounds);
        }*/
    }
    else if (type == GRID_GENERIC && !chardim.empty()
             && ((!strstr(xname, "lon") && !strstr(xdimname, "lon")) || (!strstr(yname, "lat") && !strstr(ydimname, "lat")))

                 ) /* && (strcmp(chardim, "oline") == 0 || strcmp(chardim, "basin") == 0 || strcmp(chardim, "siline") == 0)) */
    {
      if (Options::cdoVerbose) cdo_print("Unknown grid type.");
      if (Options::cdoVerbose) cdo_print("Start to define a character axis '%s' instead of a grid axis'.", chardim);
      grid_ids[0] = 0;
      int numchar = 0;
      char charvalstring[CMOR_MAX_STRING];
      std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s_%s", chardim.c_str(), cmor_name.c_str());
      std::vector<std::string> charvals = kv_get_vals(kvl, charvalstring, &numchar);
      if (numchar == 0)
      {
        std::snprintf(charvalstring, CMOR_MAX_STRING, "char_axis_%s", chardim.c_str());
        charvals = kv_get_vals(kvl, charvalstring, &numchar);
      }

      if ((xlength > 0 && xlength != numchar) && (ylength > 0 && ylength != numchar))
        cdo_abort("ERROR (infile: '%s')! In registration of a character coordinate as substitution for a horizontal axis:\n        "
                  "  You configured a character coordinate '%s' with '%d' string values but you also registered a "
                  "grid with '%d' numerical values on X axis and '%d' numerical values on Y axis. One axis must "
                  "match the number of string values.",
                  cdo_get_stream_name(0), chardim, numchar, xlength, ylength);
      if (!charvals.size())
        cdo_abort("ERROR (infile: '%s')! In registration of a character coordinate as substitution for a horizontal axis:\n        "
                  "  You configured a character coordinate '%s' but no values are found! Configure values via "
                  "attribute 'char_axis_vals'!",
                  cdo_get_stream_name(0), chardim);
      if (charvals.size() && (xlength == numchar || xlength == 0))
      {
        register_char_axis(numchar, charvals, axis_ids, chardim);
        if (ylength > 1) register_lat_axis(gridID, ylength, axis_ids);
      }
      else
      {
        if (xlength > 1) register_lon_axis(gridID, xlength, axis_ids);
        register_char_axis(numchar, charvals, axis_ids, chardim);
      }
      if (Options::cdoVerbose) cdo_print("Successfully defined a character axis '%s' instead of a grid axis.", chardim);
      kv_insert_vals(kvl, "capr", (char *) "y", true, false);
    }
    else if (type == GRID_CHARXY)
    {
      grid_ids[0] = 0;
      if (std::strcmp(xdimname, "line") == 0) std::strcpy(xdimname, "oline");
      int dimstrlen;
      if ((dimstrlen = gridInqXIsc(gridID)))
      {
        char **xchars = (char **) std::malloc((xlength + 1) * sizeof(char *));
        for (int i = 0; i < xlength; ++i) xchars[i] = (char *) std::malloc((dimstrlen + 1) * sizeof(char));
        gridInqXCvals(gridID, xchars);
        std::vector<std::string> xcharspp(xlength + 1);
        // Assign each C-style string to the corresponding element in the vector
        for (int i = 0; i < xlength; ++i) { xcharspp[i] = xchars[i]; }
        free_array(xchars);
        register_char_axis(xlength, xcharspp, axis_ids, xdimname);
      }
      else if (xlength)
        register_lon_axis(gridID, xlength, axis_ids);

      if ((dimstrlen = gridInqYIsc(gridID)))
      {
        char **ychars = (char **) std::malloc((ylength + 1) * sizeof(char *));
        for (int i = 0; i < ylength; ++i) ychars[i] = (char *) std::malloc((dimstrlen + 1) * sizeof(char));
        gridInqYCvals(gridID, ychars);
        std::vector<std::string> ycharspp(ylength + 1);
        // Assign each C-style string to the corresponding element in the vector
        for (int i = 0; i < ylength; ++i) { ycharspp[i] = ychars[i]; }
        free_array(ychars);
        register_char_axis(ylength, ycharspp, axis_ids, ydimname);
      }
      else if (ylength)
        register_lat_axis(gridID, ylength, axis_ids);
    }
    else if (type == GRID_PROJECTION)
    {
      cdo_abort("ERROR (infile: '%s')! In grid registration:\n          For a 'rotated_lat_lon' projection, both grids, the "
                "unprojected lat/lon and the projected rlat/rlon are required.",
                cdo_get_stream_name(0));
    }
    else
    {
      grid_ids[0] = 0;
      cdo_warning("Registration of a grid is skipped. Either the grid type is unknown or a registration is not necessary.");
    }

    if (projtype != CDI_UNDEFID)
    {
      register_projection(grid_ids, projID, ycoord_vals, xcoord_vals, ycell_bounds, xcell_bounds, xlength, ylength, projtype);
      std::free(xcoord_vals);
      std::free(ycoord_vals);
      std::free(xcell_bounds);
      std::free(ycell_bounds);
    }
  }
  else
    grid_ids[0] = 0;
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
}

static void
register_variable(KVList *kvl, int vlistID, int varID, int *axis_ids, struct mapping *var, int *grid_ids, std::string name)
{
  int cmf = 0;
  if (Options::cdoVerbose) cdo_print("8.5.1. Start to retrieve 'positive' and 'units'.");
  std::string positive = get_txtatt(vlistID, varID, "positive");
  std::string origname = get_txtatt(vlistID, varID, "original_name");
  std::string history = get_txtatt(vlistID, varID, "history");
  std::string varcom = get_txtatt(vlistID, varID, "variable_comment");
  char unitsc[CDI_MAX_NAME];
  vlistInqVarUnits(vlistID, varID, unitsc);
  std::string units = std::string(unitsc);
  std::string attunits = kv_get_a_val(kvl, "u", "");
  std::string attp = kv_get_a_val(kvl, "p", "");
  std::string attorigname = kv_get_a_val(kvl, "original_name", "");
  std::string attvarcom = kv_get_a_val(kvl, "vc", "");
  check_compare_set(positive, attp, "positive", " ");
  if (positive == " ") positive = "";
  check_compare_set(units, attunits, "units", "");
  check_compare_set(origname, attorigname, "original_name", " ");
  char *corigname = nullptr;
  char *cvarcom = nullptr;
  if (origname.at(0) != ' ') { corigname = (char *) origname.c_str(); }
  check_compare_set(varcom, attvarcom, "variable_comment", " ");

  if (varcom.at(0) != ' ') { cvarcom = (char *) varcom.c_str(); }
  if (Options::cdoVerbose)
    cdo_print("8.5.1. Successfully retrieved 'positive': '%s' and 'units' : '%s' and 'variable_comment': '%s'", positive, units,
              varcom);
  if (units.empty()) cdo_abort("ERROR (infile: '%s')! No units found for CMOR variable '%s'.", cdo_get_stream_name(0), name);
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = cdo_vlist_gridsizemax(vlistID);
  int zsize = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->help_var = 0;
  if (vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_FLT64)
  {
    if (!var->data)
    {
      var->charvars = 0;
      var->datatype = 'd';
      var->data = std::malloc(gridsize * zsize * sizeof(double));
    }
    *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
  }
  else
  {
    if (!var->data)
    {
      var->charvars = 0;
      var->datatype = 'f';
      var->data = std::malloc(gridsize * zsize * sizeof(float));
    }
    *(float *) missing_value = vlistInqVarMissval(vlistID, varID);
  }
  if (Options::cdoVerbose) cdo_print("8.5.2. Start to call cmor_variable.");
  if ((zaxisInqType(vlistInqVarZaxis(vlistID, varID)) != ZAXIS_HYBRID
       && zaxisInqType(vlistInqVarZaxis(vlistID, varID)) != ZAXIS_HYBRID_HALF)
      && grid_ids[0] != 0)
  {
    int *tmp_id = new_axis_id(axis_ids);
    *tmp_id = grid_ids[0];
    cmf = cmor_variable(&var->cmor_varID, (char *) name.c_str(), (char *) units.c_str(), (count_axis_ids(axis_ids)), axis_ids,
                        var->datatype, (void *) missing_value, &tolerance, (char *) positive.c_str(), corigname,
                        (char *) history.c_str(), cvarcom);
  }
  else
  {
    cmf = cmor_variable(&var->cmor_varID, (char *) name.c_str(), (char *) units.c_str(), count_axis_ids(axis_ids), axis_ids,
                        var->datatype, (void *) missing_value, &tolerance, (char *) positive.c_str(), corigname,
                        (char *) history.c_str(), cvarcom);
  }
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_variable failed!", cdo_get_stream_name(0));
  if (Options::cdoVerbose) cdo_print("8.5.2. Successfully called cmor_variable.");
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR >= 3)
  std::string deflates = kv_get_a_val(kvl, "dl", "");
  long int deflate = 0;
  if (!deflates.empty())
  {
    if (Options::cdoVerbose) cdo_print("8.5.3. Start to set deflate for variable '%s'.", name);
    if ((deflate = std::stol(deflates)))
    {
      if (deflate == -1)
        cmf = cmor_set_deflate(var->cmor_varID, 0, 0, 0);
      else
        cmf = cmor_set_deflate(var->cmor_varID, 1, 1, deflate);
    }
    else
      cmf = cmor_set_deflate(var->cmor_varID, 1, 1, 1);
    if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_variable failed!", cdo_get_stream_name(0));
    if (Options::cdoVerbose) cdo_print("8.5.3. Successfully set deflate for variable '%s'.", name);
  }
#endif
}

static void
switch_grid_info(KVList *kvl, CdoStreamID streamID, std::string grid_file, int varID)
{
  if (Options::cdoVerbose)
    cdo_print("You configured a grid_info file: '%s'. It is tested for a valid use as substitution.\n", grid_file);
  int vlistID = cdo_stream_inq_vlist(streamID);
  int nvars = vlistNvars(vlistID);
  if (nvars > 1) cdo_print("Note that the grids of the variables found first in both files are switched.");

  int byteorder = 0;
  int filetype = cdiGetFiletype(grid_file.c_str(), &byteorder);
  if ((filetype == CDI_FILETYPE_NC) || (filetype == CDI_FILETYPE_NC2) || (filetype == CDI_FILETYPE_NC4)
      || (filetype == CDI_FILETYPE_NC4C))
  {
    int streamID2 = stream_open_read_locked(grid_file.c_str());
    int vlistID2 = streamInqVlist(streamID2);
    int gridID2 = vlistInqVarGrid(vlistID2, 0);
    int zaxisID2 = vlistInqVarZaxis(vlistID2, 0);

    if (kv_get_a_val(kvl, "switch_xy", "y") == "y")
    {
      int gridID = vlistInqVarGrid(vlistID, varID);
      change_grid(vlistID, gridID, gridID2, grid_file);
    }

    if (kv_get_a_val(kvl, "switch_z", "y") == "y")
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      change_zaxis(kvl, vlistID, zaxisID, zaxisID2, grid_file);
    }

    streamClose(streamID2);
  }
  else
  {
    if (parse_kv_file(kvl, grid_file) == 0)
      cdo_abort("ERROR (infile: '%s')! File '%s' does not exist.", cdo_get_stream_name(0), grid_file);
  }
}

static void
register_all_dimensions(KVList *kvl, CdoStreamID streamID, struct mapping vars[], int table_id, std::string const &project_id,
                        int miptab_freq, int *time_axis, int *mergeIDs)
{
  int cmf = 0;
  int vlistID = cdo_stream_inq_vlist(streamID);

  std::string time_units = kv_get_a_val(kvl, "rtu", "");
  std::string attzaxis = kv_get_a_val(kvl, "za", "");

  if (Options::cdoVerbose) cdo_print("7. Start to retrieve requested variables.");

  int numvals = 0;
  std::vector<std::string> cmor_names = kv_get_vals(kvl, "cn", &numvals);

  /* Cmdlinemapping: */
  std::string mapname, mapcode;
  if (kv_get_a_val(kvl, "mt", "").empty() && numvals)
  {
    if (Options::cdoVerbose) cdo_print("7.1. Start to search for a command line mapping");
    mapname = kv_get_a_val(kvl, "n", "");
    mapcode = kv_get_a_val(kvl, "c", "");
    if (!mapname.empty())
    {
      if (change_name_via_name(vlistID, mapname, cmor_names[0].c_str()))
      {
        if (Options::cdoVerbose)
          cdo_print("7.1. Successfully mapped '%s' on '%s' via command line attribute name.", mapname, cmor_names[0]);
      }
    }
    else if (!mapcode.empty())
    {
      if (change_name_via_code(vlistID, mapcode, cmor_names[0].c_str()))
      {
        if (Options::cdoVerbose) cdo_print("7.1. Successfully mapped via command line attribute name.");
      }
    }
    else if (Options::cdoVerbose)
      cdo_print("7.1. No command line mapping given.");
  }

  if (!cmor_names.size() && vlistNvars(vlistID) > 1)
    cdo_print("Function 'all axes registration':\n          You have not requested a particular variable via "
              "'cmor_name'.\n          There are several in infile and all will be processed.\n          Notice that "
              "attributes specified in the cmdline will be used for all infile variables.");
  if (Options::cdoVerbose) cdo_print("7. Successfully retrieved requested variables");
  int foundName = 0, psRequired = 0;
  std::map<int, int> timeAxes;

  for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
  {
    std::string zaxis = get_txtatt(vlistID, varID, "z_axis");
    check_compare_set(zaxis, attzaxis, "z_axis", "notSet");
    if (zaxis == "hybrid_height") kv_insert_vals(kvl, "za", zaxis, true, false);

    char name[CDI_MAX_NAME];
    vlistInqVarName(vlistID, varID, name);
    if (!cmor_names.size() || in_list(cmor_names, std::string(name), numvals))
    {
      struct mapping *var = map_var(varID, vars);
      if (!var)
        var = new_var_mapping(vars);
      else if (Options::cdoVerbose)
        cdo_print("Already mapped '%d'", varID);
      var->cdi_varID = varID;
      std::string grid_file = kv_get_a_val(kvl, "gi", "");
      if (!grid_file.empty()) switch_grid_info(kvl, streamID, grid_file, varID);
      int axis_ids[CMOR_MAX_AXES];
      axis_ids[0] = CMOR_UNDEFID;
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int zsize = zaxisInqSize(zaxisID);
      if (Options::cdoVerbose) cdo_print("8. Start to define variable with ID: '%d' and name: '%s'", varID, name);
      if ((zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF || zaxis == "hybrid_height")
          && zsize > 1)
      {
        if (Options::cdoVerbose)
          cdo_print("Since the zaxis of variable '%s' is of type HYBRID, surface pressure or orog need to be available in "
                    "the input file to fully describe the axis.",
                    name);
        if (zaxis != "hybrid_height") psRequired++;
      }
      /* Time-Axis */
      if (Options::cdoVerbose) cdo_print("8.1. Start to register time axis of '%s'", name);
      std::string cmor_time_name;
      /*char cmor_time_name[CMOR_MAX_STRING]; */
      /*cmor_time_name[0] = '\0'; */
      get_time_method(kvl, vlistID, varID, cmor_time_name, project_id, miptab_freq, time_axis);
      if (cmor_time_name != "none")
      {
        auto search = timeAxes.find(*time_axis);
        int *cmorTimeAxisID = new_axis_id(axis_ids);
        if (foundName)
        {
          if (search != timeAxes.end())
          {
            *cmorTimeAxisID = search->second;
            if (Options::cdoVerbose) cdo_print("8.1. Use already defined time axis '%d'", *cmorTimeAxisID);
          }
        }
        if (*cmorTimeAxisID == CMOR_UNDEFID)
        {
          cmf = cmor_axis(cmorTimeAxisID, (char *) cmor_time_name.c_str(), (char *) time_units.c_str(), 0, nullptr, 0, nullptr, 0,
                          nullptr);
          if (((project_id == "CMIP5") || (project_id == "CORDEX"))
              && ((kv_get_a_val(kvl, "realization", " ").at(0) == '0')
                  || (kv_get_a_val(kvl, "initialization_method", " ").at(0) == '0')
                  || (kv_get_a_val(kvl, "physics_version", " ").at(0) == '0')))
            cdo_warning("At least one ensemble index is set to '0' while cell_methods is not 'none'!\n"

                        "          '0' is usually reserved for fixed fields!");
        }
        if (search == timeAxes.end()) timeAxes.insert(std::pair<int, int>(*time_axis, axis_ids[count_axis_ids(axis_ids) - 1]));
      }
      if (Options::cdoVerbose && cmf == 0)
        cdo_print("8.1. Successfully handled time axis registration.");
      else if (cmf != 0)
        cdo_abort("ERROR (infile: '%s')! Function cmor_axis failed!", cdo_get_stream_name(0));
      /* Grid: */
      if (Options::cdoVerbose) cdo_print("8.2. Start to register grid of '%s'", name);
      int grid_ids[CMOR_MAX_GRIDS];
      register_grid(kvl, vlistID, varID, axis_ids, grid_ids, project_id, name);
      if ((zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) && grid_ids[0] != 0)
      {
        int *tmp_id = new_axis_id(axis_ids);
        *tmp_id = grid_ids[0];
      }
      cmf = cmor_set_table(table_id);
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_table failed!", cdo_get_stream_name(0));
      if (Options::cdoVerbose) cdo_print("8.2. Successfully handled grid registration.");
      if (miptab_freq == 19)
      {
        /* Possible 4th spatial dimension */
        /* Z-Axis */
        if (Options::cdoVerbose) cdo_print("8.4. Start to register 4th spatial dimension of '%s'", name);
        register_fourth_axis(kvl, vlistID, varID, std::string(name), axis_ids, project_id, miptab_freq, mergeIDs);
        if (Options::cdoVerbose) cdo_print("8.4. Successfully handled 4th spatial dimension registration.");
        /* Z-Axis */
        if (Options::cdoVerbose) cdo_print("8.3. Start to register zaxis of '%s'", name);
        register_z_axis(kvl, vlistID, varID, zaxisID, std::string(name), zaxis, axis_ids, &var->zfactor_id, project_id, miptab_freq,
                        mergeIDs, vars, streamID);
        if (Options::cdoVerbose) cdo_print("8.3. Successfully handled zaxis registration.");
      }
      else
      {
        /* Z-Axis */
        if (Options::cdoVerbose) cdo_print("8.3. Start to register zaxis of '%s'", name);
        register_z_axis(kvl, vlistID, varID, zaxisID, std::string(name), zaxis, axis_ids, &var->zfactor_id, project_id, miptab_freq,
                        mergeIDs, vars, streamID);
        if (Options::cdoVerbose) cdo_print("8.3. Successfully handled zaxis registration.");
        /* Possible 4th spatial dimension */
        /* Z-Axis */
        if (Options::cdoVerbose) cdo_print("8.4. Start to register 4th spatial dimension of '%s'", name);
        register_fourth_axis(kvl, vlistID, varID, std::string(name), axis_ids, project_id, miptab_freq, mergeIDs);
        if (Options::cdoVerbose) cdo_print("8.4. Successfully handled 4th spatial dimension registration.");
      }
      /* Variable */
      if (Options::cdoVerbose) cdo_print("8.5. Start to register variable '%s'", name);
      register_variable(kvl, vlistID, varID, axis_ids, var, grid_ids, name);
      if (Options::cdoVerbose) cdo_print("8.5. Successfully handled variable registration.");
      if (Options::cdoVerbose) cdo_print("8. Successfully defined variable with ID: '%d' and name: '%s'.", varID, name);
      kv_insert_vals(kvl, "capr", "n", true, false);
      foundName++;
    }
  }
  if (mergeIDs[0] != CMOR_UNDEFID)
  {
    int refdatatype = vlistInqVarDatatype(vlistID, mergeIDs[0]);
    int mergeID = 0;
    while (mergeIDs[mergeID] != CMOR_UNDEFID)
    {
      struct mapping *var = map_var(mergeIDs[mergeID], vars);
      if (!var)
      {
        var = new_var_mapping(vars);
        var->charvars = 0;
        var->cdi_varID = mergeIDs[mergeID];
        if (vlistInqVarDatatype(vlistID, mergeIDs[mergeID]) != refdatatype)
          cdo_abort("ERROR (infile: '%s')! Variable with ID '%d' has datatype '%d' but"
                    " variable with id '%d' has datatype '%d'."
                    " All variables that should be merged needs to be"
                    " of the same datatype.",
                    cdo_get_stream_name(0), mergeIDs[0], refdatatype, mergeIDs[mergeID],
                    vlistInqVarDatatype(vlistID, mergeIDs[mergeID]));
        if (vlistInqVarDatatype(vlistID, mergeIDs[mergeID]) == CDI_DATATYPE_FLT64)
        {
          var->datatype = 'd';
          var->data = std::malloc(gridInqSize(vlistInqVarGrid(vlistID, mergeIDs[mergeID]))
                                  * zaxisInqSize(vlistInqVarZaxis(vlistID, mergeIDs[mergeID])) * sizeof(double));
        }
        else
        {
          var->datatype = 'f';
          var->data = std::malloc(gridInqSize(vlistInqVarGrid(vlistID, mergeIDs[mergeID]))
                                  * zaxisInqSize(vlistInqVarZaxis(vlistID, mergeIDs[mergeID])) * sizeof(float));
        }
      }
      mergeID++;
    }
  }
  if (psRequired && cmor_names.size() && !in_list(cmor_names, "ps", numvals))
  {
    int psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
    if (psindex == CDI_UNDEFID) psindex = getVarIDToMap(vlistID, vlistNvars(vlistID), "code", "134");
    int psID = getRegisteredPsid(vars, psindex);
    if (psID == CDI_UNDEFID) cdo_abort("ERROR (infile: '%s')! Could not find ps or orog for hybrid axis.");
    vars[psID].help_var = 1;
    if (Options::cdoVerbose) cdo_print("9. Set ps as a auxiliary variable.");
  }
  if (!foundName && cmor_names.size())
    cdo_abort(
        "ERROR (infile: '%s')! After registration of all dimensions for all variables:\n          None of the given variables to "
        "process by attribute 'cmor_name' is found in infile.",
        cdo_get_stream_name(0));
  if (Options::cdoVerbose) cdo_print("Successfully registered all dimensions for %d variables successfully.", foundName);
}

static std::string
get_frequency(/*KVList *kvl,*/ int vlistID, int miptab_freq)
{
  char frequency[CDI_MAX_NAME];
  std::strcpy(frequency, "no");
  int ntsteps = vlistNtsteps(vlistID);
  int reccounter = 0;
  int recdummy = 0;

  switch (miptab_freq)
  {
    case 11: std::strcpy(frequency, "yr"); break;
    case 2: std::strcpy(frequency, "yr"); break;
    case 12: std::strcpy(frequency, "mon"); break;
    case 3: std::strcpy(frequency, "mon"); break;
    case 13: std::strcpy(frequency, "day"); break;
    case 4: std::strcpy(frequency, "day"); break;
    case 14: std::strcpy(frequency, "6hr"); break;
    case 5: std::strcpy(frequency, "6hr"); break;
    case 6: std::strcpy(frequency, "6hr"); break;
    case 15: std::strcpy(frequency, "3hr"); break;
    case 7: std::strcpy(frequency, "1hr"); break;
    case 8: std::strcpy(frequency, "3hr"); break;
    case 16: std::strcpy(frequency, "1hr"); break;
    case 17: std::strcpy(frequency, "sem"); break;
    case 18: std::strcpy(frequency, "dec"); break;
    case 19: std::strcpy(frequency, "subhr"); break;
    default:
    {
      if (cdo_assert_files_only() == false)
      {
        cdo_abort("ERROR (infile: '%s')! No frequency could be determined from MIP-table and, additionally,"
                  " cdo cmor cannot check frequency of "
                  "Ifile recs since you piped several cdo operators.",
                  cdo_get_stream_name(0));
        /*          char *dummy;
                  cdo_warning("Cdo cmor cannot check frequency of Ifile recs since you piped several cdo
           operators.\nIt is tried to use a configuration attribute frequency.");
                  if ( !(dummy = kv_get_a_val(kvl, "frequency", nullptr)) )
                    cdo_abort("ERROR (infile: '%s')! No attribute frequency is found.");
                  else
                    {
                      std::strcpy(frequency, dummy);
                      return frequency;
                    }
        */
      }

      CdiStreamID streamID2 = streamOpenRead(cdo_get_stream_name(0));
      int vlistID2 = streamInqVlist(streamID2);
      int taxisID2 = vlistInqTaxis(vlistID2);
      if (ntsteps < 0)
      {
        while ((recdummy = streamInqTimestep(streamID2, reccounter++)));
        ntsteps = reccounter;
      }
      ntsteps -= 1;
      int fyear, lyear, fmonth, lmonth, dummytwo;

      if (ntsteps > 2)
      {
        streamInqTimestep(streamID2, 0);
        cdiDate_decode(taxisInqVdatetime(taxisID2).date, &fyear, &fmonth, &dummytwo);
        streamInqTimestep(streamID2, ntsteps);
        cdiDate_decode(taxisInqVdatetime(taxisID2).date, &lyear, &lmonth, &dummytwo);

        double covered_years = lyear - fyear + 1.0;
        double ntperyr = (double) ((ntsteps + 1) / covered_years);
        if (fp_is_equal(ntperyr, (double) 1))
          std::strcpy(frequency, "yr");
        else if (fp_is_equal(ntperyr, (double) 12))
          std::strcpy(frequency, "mon");
        else if (fp_is_equal(ntperyr, (double) 365) || fp_is_equal(ntperyr, (double) 365.25) || fp_is_equal(ntperyr, (double) 366))
          std::strcpy(frequency, "day");
        else if (fp_is_equal(ntperyr, (double) 365 * 4) || fp_is_equal(ntperyr, (double) 365.25 * 4)
                 || fp_is_equal(ntperyr, (double) 366 * 4))
          std::strcpy(frequency, "6hr");
        else if (fp_is_equal(ntperyr, (double) 365 * 8) || fp_is_equal(ntperyr, (double) 365.25 * 8)
                 || fp_is_equal(ntperyr, (double) 366 * 8))
          std::strcpy(frequency, "3hr");
        else
        {
          int step_per_year = 0;
          reccounter = 0;
          if (Options::cdoVerbose)
            cdo_print("Frequency could not be determined by comparing all time steps (%d) divided by covered "
                      "years (%f).\n          It is now calculated by counting all timesteps in year %d\n         "
                      " in order to calculate time bounds in case they are not given.",
                      ntsteps, covered_years, fyear);
          while ((recdummy = streamInqTimestep(streamID2, reccounter++)))
          {
            int reqyear;
            cdiDate_decode(taxisInqVdatetime(taxisID2).date, &reqyear, &lmonth, &dummytwo);
            if (reqyear == (fyear + 1)) break;
            step_per_year++;
          }
          int covered_months = lmonth - fmonth + 1;
          if (step_per_year > 366 * 8)
            cdo_abort("ERROR (infile: '%s')! In estimating frequency:\n          Frequency is sub-3hourly! Not yet enabled.",
                      cdo_get_stream_name(0));
          else
          {
            if ((double) step_per_year / (double) covered_months > 31 * 8)
              cdo_abort("ERROR (infile: '%s')! Frequency is sub-3hourly! Not yet enabled.", cdo_get_stream_name(0));
            else if ((double) step_per_year / (double) covered_months > 31 * 4)
              std::strcpy(frequency, "3hr");
            else if ((double) step_per_year / (double) covered_months > 31)
              std::strcpy(frequency, "6hr");
            else if ((double) step_per_year / (double) covered_months > 1)
              std::strcpy(frequency, "day");
            else
              std::strcpy(frequency, "mon");
          }
          if (Options::cdoVerbose)
            cdo_print("Found %d time steps in year %d.\n          Therefore, the frequency is %s.", step_per_year, fyear,
                      frequency);
        }
      }
      else
      {
        if (!taxisHasBounds(taxisID2) && ntsteps > 0)
          cdo_abort("ERROR (infile: '%s')! In estimating frequency:\n          No time bounds are found in infile and for %d found "
                    "timesteps no frequency can be computed - at least 3 timesteps are required.\n          Define "
                    "time bounds before cdo cmor.",
                    cdo_get_stream_name(0), ntsteps);
        else
          cdo_warning("In frequency estimation:\n          For %d found timesteps no frequency can be computed - at "
                      "least 3 timesteps are required.\n          Time bounds of the rec are used.",
                      ntsteps);
      }
      streamClose(streamID2);
    }
  }
  std::string cppfrequency(frequency);
  return cppfrequency;
}

static int
get_tunitsec(int tunit)
{
  switch (tunit)
  {
    case TUNIT_MINUTE: return 60;
    case TUNIT_HOUR: return 3600;
    case TUNIT_DAY: return 86400;
    default: return 3600;
  }
}

static JulianDate
get_cmor_time_val(KVList *kvl, int taxisID, JulianDate ref_date, int /*tunitsec*/, int calendar, std::string const &frequency,
                  int ts_id, int time_axis)
{
  auto vDateTime = taxisInqVdatetime(taxisID);
  int year, month, day, hour, min, sec, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  auto juldate = julianDate_encode(calendar, vDateTime);

  if (month == 0 || day == 0)
  {
    int timeoffset;
    if ((timeoffset = std::stoi(kv_get_a_val(kvl, "firsttimeval", "-99"))) < 0)
      cdo_abort("ERROR (infile: '%s')! Time axis is broken (month or day = 0).\n          Provide 'timeoffset' and the operator "
                "tries to calculate time values with frequency.",
                cdo_get_stream_name(0));
    else
    {
      int ryear, rmonth, rday, addseconds = 0;
      auto rDateTime = julianDate_decode(calendar, ref_date);
      cdiDate_decode(rDateTime.date, &ryear, &rmonth, &rday);
      /* Only print this for the first time step */
      if (ts_id < 2)
        cdo_warning("In writing the data:\n          Time axis is broken (month or day = 0). It is tried to "
                    "calculate time values with frequency and timeoffset ignoring the time stamp year.\n          "
                    "Note: These can only be valid if\n           - cm=m \n           - a equally spaced "
                    "monotonical time axis exist according to the frequency \n           - a correct calendar "
                    "exist!");
      /**/
      /* First record is valid for correfdate = refdate + timeoffset. */
      /**/
      while (timeoffset != 0)
      {
        if (timeoffset > 11)
        {
          ryear += 1;
          timeoffset -= 12;
        }
        else if (timeoffset != 0)
        {
          rmonth += timeoffset;
          if (rmonth > 12)
          {
            ryear += 1;
            rmonth -= 12;
          }
          timeoffset = 0;
        }
      }
      rDateTime.date = cdiDate_encode(ryear, rmonth, rday);
      ref_date = julianDate_encode(calendar, rDateTime);
      /**/
      /* Add time index * frequency on correfdate */
      /**/
      if (frequency == "yr")
      {
        year = ryear + ts_id;
        month = 6; /* Is set to mid point by CMOR */
        day = 14;  /* Is set to mid point by CMOR */
      }
      else if (frequency == "mon")
      {
        year = ryear + std::floor(((double) (ts_id - 1)) / 12);
        month = (ts_id % 12);
        if (month == 0) month = 12;
        day = 14; /* Is set to mid point by CMOR */
      }
      else if (frequency == "day")
      {
        addseconds = ts_id * 24 * 60 * 60 + 60 * 60 * 12;
        juldate = julianDate_add_seconds(ref_date, addseconds);
      }
      else if (frequency == "6hr")
      {
        addseconds = ts_id * 6 * 60 * 60;
        juldate = julianDate_add_seconds(ref_date, addseconds);
      }
      else if (frequency == "3hr")
      {
        addseconds = ts_id * 3 * 60 * 60;
        juldate = julianDate_add_seconds(ref_date, addseconds);
      }
      else if (frequency == "1hr")
      {
        addseconds = ts_id * 60 * 60;
        juldate = julianDate_add_seconds(ref_date, addseconds);
      }
      if (addseconds == 0)
      {
        CdiDateTime dt{};
        dt.date = cdiDate_encode(year, month, 1);
        juldate = julianDate_encode(calendar, dt);
      }
    }
  }
  if (time_axis == 1 && (kv_get_a_val(kvl, "ta", "n") == "cmip"))
  {
    auto julDateTime = julianDate_decode(calendar, juldate);
    cdiTime_decode(julDateTime.time, &hour, &min, &sec, &ms);
    if (frequency == "6hr")
    {
      int iv[] = { 0, 6, 12, 18, 24 };
      size_t ivsize = sizeof(iv) / sizeof(iv[0]);
      int minid = 0;
      julDateTime.time = cdiTime_encode(0, 0, 0, 0);
      double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

      for (size_t loopid = 1; loopid < ivsize; loopid++)
      {
        julDateTime.time = cdiTime_encode(iv[loopid], 0, 0, 0);
        if (std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff)
        {
          minid = loopid;
          julDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
          diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
        }
      }

      vDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
      juldate = julianDate_encode(calendar, vDateTime);
    }
    else if (frequency == "3hr")
    {
      int iv[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24 };
      size_t ivsize = sizeof(iv) / sizeof(iv[0]);
      int minid = 0;
      julDateTime.time = cdiTime_encode(0, 0, 0, 0);
      double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

      for (size_t loopid = 1; loopid < ivsize; loopid++)
      {
        julDateTime.time = cdiTime_encode(iv[loopid], 0, 0, 0);
        if (std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff)
        {
          minid = loopid;
          julDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
          diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
        }
      }

      vDateTime.time = cdiTime_encode(iv[minid], 0, 0, 0);
      juldate = julianDate_encode(calendar, vDateTime);
    }
    else if (frequency == "1hr")
    {
      size_t ivsize = 25;
      int minid = 0;
      julDateTime.time = cdiTime_encode(0, 0, 0, 0);
      double diff = julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)));

      for (size_t loopid = 1; loopid < ivsize; loopid++)
      {
        julDateTime.time = cdiTime_encode(loopid, 0, 0, 0);
        if (std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime)))) < diff)
        {
          minid = loopid;
          julDateTime.time = cdiTime_encode(minid, 0, 0, 0);
          diff = std::fabs(julianDate_to_seconds(julianDate_sub(juldate, julianDate_encode(calendar, julDateTime))));
        }
      }

      julDateTime.time = cdiTime_encode(minid, 0, 0, 0);
      juldate = julianDate_encode(calendar, julDateTime);
    }
  }
  return juldate;
}

static double *
get_time_bounds(KVList *kvl, int taxisID, int ifreq, JulianDate ref_date, JulianDate jtime_val, int calendar, int tunitsec,
                double *time_bnds, int time_axis, int /*vlistID*/)
{
  double time_val = julianDate_to_seconds(julianDate_sub(jtime_val, ref_date)) / tunitsec;
  auto vDateTime = taxisInqVdatetime(taxisID);
  auto vDateTimeCorr = vDateTime;
  CdiDateTime vDateTime0b{};
  CdiDateTime vDateTime1b{};
  int year, month, day;
  int hour, min, sec, ms;
  cdiDate_decode(vDateTime.date, &year, &month, &day);
  if (month == 0 || day == 0)
  {
    vDateTimeCorr = julianDate_decode(calendar, jtime_val);
    cdiDate_decode(vDateTimeCorr.date, &year, &month, &day);
  }
  /***/
  /* If file time axis has bounds use them, otherwise use cmor time axis deduced from miptable frequency and
   * cell_methods or frequency itself*/
  /***/

  if (!taxisHasBounds(taxisID) || (kv_get_a_val(kvl, "ta", "n") == "cmip") || time_axis == 2 || ifreq == 8)
  {
    /***/
    /* Climatologies */
    /***/

    if (time_axis == 2)
    {
      if (Options::cdoVerbose) cdo_print("10.4.1. Get climatology interval.");
      int numdates = 0;
      std::vector<std::string> climyears = kv_get_vals(kvl, "ci", &numdates);
      if (!climyears.size())
        cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for climatology time "
                  "axis because attribute 'climatology_interval' is not available.",
                  cdo_get_stream_name(0));
      if (numdates != 2)
        cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for climatology time "
                  "axis because attribute 'climatology_interval' has not two values.",
                  cdo_get_stream_name(0));
      int expstartyear = std::stoi(climyears[0]);
      int expendyear = std::stoi(climyears[1]);

      vDateTime0b.date = cdiDate_encode(expstartyear, month, 1);
      month++;
      if (month != 12)
        vDateTime1b.date = cdiDate_encode(expendyear, month, 1);
      else
        vDateTime1b.date = cdiDate_encode(expendyear + 1, 1, 1);
    }
    /***/
    /* Diurnal cycle */
    /***/
    else if (time_axis == 3)
    {
      vDateTime0b.date = cdiDate_encode(year, month, 1);
      cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
      vDateTime0b.time = cdiTime_encode(hour, 0, 0, 0);

      hour++;
      month++;
      vDateTime1b.time = cdiTime_encode(hour, 0, 0, 0);
      vDateTime1b.date = (hour > 23) ? cdiDate_encode(year, month, 2) : cdiDate_encode(year, month, 1);
    }
    else
    {
      /***/
      /* Frequency dependent: */
      /***/
      if (ifreq == 1)
      {
        vDateTime0b.date = cdiDate_encode(year, 1, 1);
        vDateTime1b.date = cdiDate_encode(year + 1, 1, 1);
      }
      else if (ifreq == 7)
      {
        if (month == 12 || month == 1 || month == 2)
        {
          if (month == 12)
          {
            vDateTime0b.date = cdiDate_encode(year, 12, 1);
            vDateTime1b.date = cdiDate_encode(year + 1, 3, 1);
          }
          else
          {
            vDateTime0b.date = cdiDate_encode(year - 1, 12, 1);
            vDateTime1b.date = cdiDate_encode(year, 3, 1);
          }
        }
        else if (month == 3 || month == 4 || month == 5)
        {
          vDateTime0b.date = cdiDate_encode(year, 3, 1);
          vDateTime1b.date = cdiDate_encode(year, 6, 1);
        }
        else if (month == 6 || month == 7 || month == 8)
        {
          vDateTime0b.date = cdiDate_encode(year, 6, 1);
          vDateTime1b.date = cdiDate_encode(year, 9, 1);
        }
        else
        {
          vDateTime0b.date = cdiDate_encode(year, 9, 1);
          vDateTime1b.date = cdiDate_encode(year, 12, 1);
        }
      }
      else if (ifreq == 2)
      {
        vDateTime0b.date = cdiDate_encode(year, month, 1);
        month++;
        if (month > 12)
        {
          month = 1;
          year++;
        }
        vDateTime1b.date = cdiDate_encode(year, month, 1);
      }
      else if (ifreq == 3)
      {
        /*Assuming that the requested axis is always "days since.. 00:00:00" */
        time_bnds[0] = std::floor(time_val);
        time_bnds[1] = std::ceil(time_val);
        /*Assuming that the averaged value is written for and at the end of the interval */
        if (fp_is_equal(time_bnds[0], time_bnds[1])) time_bnds[0] -= 1.;
        return time_bnds;
      }
      /***/
      /* Note that time_val must be correct in Infile for subdaily frequencies */
      /***/
      else if (ifreq == 4)
      {
        cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
        int iv[] = { 0, 6, 12, 18, 24 };
        size_t ivsize = sizeof(iv) / sizeof(iv[0]);
        int minid = 0, newhour[2];
        vDateTimeCorr.time = cdiTime_encode(0, 0, 0, 0);
        double diff = julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)));
        for (size_t loopid = 1; loopid < ivsize; loopid++)
        {
          vDateTimeCorr.time = cdiTime_encode(iv[loopid], 0, 0, 0);
          if (std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)))) < diff)
          {
            minid = loopid;
            vDateTimeCorr.time = cdiTime_encode(iv[minid], 0, 0, 0);
            diff = std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr))));
          }
        }
        newhour[0] = (hour - iv[minid] < 0) ? iv[minid - 1] : iv[minid];
        newhour[1] = (hour - iv[minid] < 0) ? iv[minid] : iv[minid + 1];

        vDateTime0b.time = cdiTime_encode(newhour[0], 0, 0, 0);
        vDateTime1b.time = cdiTime_encode(newhour[1], 0, 0, 0);

        vDateTime0b.date = cdiDate_encode(year, month, day);
        vDateTime1b.date = vDateTime0b.date;
      }
      else if (ifreq == 5)
      {
        cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
        int iv[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24 };
        size_t ivsize = sizeof(iv) / sizeof(iv[0]);
        int minid = 0, newhour[2];
        vDateTimeCorr.time = cdiTime_encode(0, 0, 0, 0);
        double diff = julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)));
        for (size_t loopid = 1; loopid < ivsize; loopid++)
        {
          vDateTimeCorr.time = cdiTime_encode(iv[loopid], 0, 0, 0);
          if (std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr)))) < diff)
          {
            minid = loopid;
            vDateTimeCorr.time = cdiTime_encode(iv[minid], 0, 0, 0);
            diff = std::fabs(julianDate_to_seconds(julianDate_sub(jtime_val, julianDate_encode(calendar, vDateTimeCorr))));
          }
        }
        newhour[0] = (hour - iv[minid] < 0) ? iv[minid - 1] : iv[minid];
        newhour[1] = (hour - iv[minid] < 0) ? iv[minid] : iv[minid + 1];

        vDateTime0b.time = cdiTime_encode(newhour[0], 0, 0, 0);
        vDateTime1b.time = cdiTime_encode(newhour[1], 0, 0, 0);

        vDateTime0b.date = cdiDate_encode(year, month, day);
        vDateTime1b.date = vDateTime0b.date;
      }
      else if (ifreq == 6)
      {
        cdiTime_decode(vDateTime.time, &hour, &min, &sec, &ms);
        vDateTime0b.time = cdiTime_encode(hour, 0, 0, 0);
        vDateTime1b.time = cdiTime_encode(hour + 1, 0, 0, 0);

        vDateTime0b.date = cdiDate_encode(year, month, day);
        vDateTime1b.date = vDateTime0b.date;
      }
      else if (ifreq == 8)
      {
        if (Options::cdoVerbose) cdo_print("10.4.1. Get decadal interval.");
        int numdates = 0;
        std::vector<std::string> decyears = kv_get_vals(kvl, "di", &numdates);
        if (!decyears.size())
          cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for decadal time "
                    "axis because attribute 'decadal_interval' is not available.",
                    cdo_get_stream_name(0));
        if (numdates != 2)
          cdo_abort("ERROR (infile: '%s')! In writing model output:\n          Could not calculate time bounds for decadal time "
                    "axis because attribute 'decadal_interval' has not two values.",
                    cdo_get_stream_name(0));
        int expstartyear = std::stoi(decyears[0]);
        int expendyear = std::stoi(decyears[1]);

        vDateTime0b.date = cdiDate_encode(expstartyear, 1, 1);
        vDateTime1b.date = cdiDate_encode(expendyear + 1, 1, 1);
      }
    }
  }
  else
  {
    taxisInqVdatetimeBounds(taxisID, &vDateTime0b, &vDateTime1b);
  }

  auto juldate = julianDate_encode(calendar, vDateTime0b);
  time_bnds[0] = julianDate_to_seconds(julianDate_sub(juldate, ref_date)) / tunitsec;

  juldate = julianDate_encode(calendar, vDateTime1b);
  time_bnds[1] = julianDate_to_seconds(julianDate_sub(juldate, ref_date)) / tunitsec;
  if (time_axis == 3) time_bnds[1] -= 1;

  return time_bnds;
}

static void
read_record(CdoStreamID streamID, struct mapping vars[], int vlistID)
{
  auto [varID, levelID] = cdo_inq_field(streamID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int type = gridInqType(gridID);
  auto gridsize = gridInqSize(gridID);
  double *buffer = (double *) std::malloc(gridsize * sizeof(double));

  struct mapping *var = map_var(varID, vars);
  if (var && var->charvars != 1)
  {
    int zaxisID = vlistInqVarZaxis(vlistID, varID);
    int ztype = zaxisInqType(zaxisID);
    /*      int latdim = gridInqYsize(gridID); */
    int levdim = zaxisInqSize(zaxisID);
    size_t nmiss;
    cdo_read_field(streamID, buffer, &nmiss);
    for (size_t i = 0; i < gridsize; ++i)
    {
      // Wrong:  (lat x basin, lev ) gridsize * levelID + i
      // Wrong:  (basin x lat, lev) gridsize * levelID + i * chardim - ( int ) std::floor(i / latdim) * gridsize + ( int
      // ) std::floor(i/latdim)
      // Wrong:  (basin x lev, lat ) gridsize/latdim * levdim * ( i - ( int ) std::floor(i/latdim) * latdim ) + ( int )
      // std::floor(i/latdim) + gridsize/latdim * levelID;
      // Wrong:  (lat x lev, basin ) latdim * levdim * ( int ) std::floor(i/latdim) + ( i - ( int ) std::floor(i/latdim) *
      // latdim ) + levelID * latdim
      // (lev x lat, basin )
      int newIndex;
      if (levdim > 1)
      {
        if ((type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR) && ztype != ZAXIS_HYBRID) newIndex = i + gridsize * levelID;
        //              else if ( type == GRID_CURVILINEAR )
        //                newIndex = i + gridsize * levelID;
        else
          newIndex = i * levdim + levelID;
      }
      else
        newIndex = i;
      if (var->datatype == 'f') { ((float *) var->data)[newIndex] = (float) buffer[i]; }
      else
      {
        ((double *) var->data)[newIndex] = (double) buffer[i];
      }
    }
  }
  std::free(buffer);
}

static int
check_append_and_size(KVList *kvl, char *testIn, int ifreq, int calendar)
{
  char *test = testIn;
  size_t filesize = FileUtils::size((const char *) testIn);
  char old_start_date[CMOR_MAX_STRING];
  char old_end_date[CMOR_MAX_STRING];
  int i = 0, j = 0;
  /* Get dates from chunk string */
  if (Options::cdoVerbose) cdo_print("Start to retrieve dates from chunk string.");
  while (*(test + i) != 0)
  {
    if (*(test + i) == '_')
    {
      test += (i + 1);
      i = 0;
    }
    if (*(test + i) == '-') j = i;
    i++;
  }
  if (!i || !j || *(test + j + 1) == 0 || *(test + 2 * j) == 0)
  {
    if (Options::cdoVerbose)
      cdo_print("In checking the last chunk:\n          Date from filename of the chunk cannot be read.\n          "
                "Switched to replace mode for this variable.");
    return 0;
  }

  strncpy(old_start_date, test, j);
  old_start_date[j] = 0;
  test += (j + 1);
  strncpy(old_end_date, test, j);
  old_end_date[j] = 0;

  if (Options::cdoVerbose)
    cdo_print("Successfully retrieved start date: '%s' and end date: '%s' chunk string.", old_start_date, old_end_date);
  /* Check frequency of chunk with frequency of file */

  if ((j == 12 && ifreq < 4) || (j == 8 && ifreq != 3) || (j == 6 && ifreq != 2) || (j == 4 && ifreq != 1 && ifreq != 8))
  {
    if (Options::cdoVerbose)
      cdo_print("In checking last chunk:\n          Frequency of chunk file does not agree with frequency of the "
                "working file.\n         Switched to replace mode for this variable.");
    return 0;
  }

  /* Encode in julseconds depending on frequency */
  if (Options::cdoVerbose) cdo_print("Start to encode dates with frequencies to julseconds.");

  int old_start_year = 0, old_start_month = 1, old_start_day = 1, old_start_hr = 1, old_start_min = 1;
  int old_end_year = 0, old_end_month = 1, old_end_day = 1, old_end_hr = 1, old_end_min = 1;

  switch (j)
  {
    case (12):
      std::sscanf(old_start_date, "%04d%02d%02d%02d%02d", &old_start_year, &old_start_month, &old_start_day, &old_start_hr,
                  &old_start_min);
      std::sscanf(old_end_date, "%04d%02d%02d%02d%02d", &old_end_year, &old_end_month, &old_end_day, &old_end_hr, &old_end_min);
      break;
    case (8):
      std::sscanf(old_start_date, "%04d%02d%02d", &old_start_year, &old_start_month, &old_start_day);
      std::sscanf(old_end_date, "%04d%02d%02d", &old_end_year, &old_end_month, &old_end_day);
      break;
    case (6):
      std::sscanf(old_start_date, "%04d%02d", &old_start_year, &old_start_month);
      std::sscanf(old_end_date, "%04d%02d", &old_end_year, &old_end_month);
      break;
    case (4):
      old_start_year = atol(old_start_date);
      old_end_year = atol(old_end_date);
      break;
    default:
    {
      if (kv_get_a_val(kvl, "om", "r") == "a")
      {
        cdo_warning("Last chunk's frequency cannot yet be tested for being suitable.");
        return 1;
      }
      else if (Options::cdoVerbose)
      {
        cdo_print("In checking last chunk:\n          Last chunk has unknown frequency "
                  "which is yet not enabled to be appended by "
                  "cdo cmor.\n          Switched to replace mode for this variable.");
        return 0;
      }
    }
  }

  CdiDateTime startDateTime{};
  CdiDateTime endDateTime{};
  startDateTime.date = cdiDate_encode(old_start_year, old_start_month, old_start_day);
  endDateTime.date = cdiDate_encode(old_end_year, old_end_month, old_end_day);

  startDateTime.time = cdiTime_encode(old_start_hr, old_start_min, 0, 0);
  endDateTime.time = cdiTime_encode(old_end_hr, old_end_min, 0, 0);
  auto julostart = julianDate_encode(calendar, startDateTime);
  auto juloend = julianDate_encode(calendar, endDateTime);

  if (Options::cdoVerbose) cdo_print("Successfully calculated juldates.");
  /* Read in first vdate in case not piped */
  if (Options::cdoVerbose) cdo_print("Start to calculate temporal gap between chunk and working file.");
  if (cdo_assert_files_only() == false)
  {
    if (Options::cdoVerbose)
      cdo_print("Cdo cmor cannot enable append mode since you piped several cdo operators.\n          Switched to "
                "replace mode for this variable.");
    if (kv_get_a_val(kvl, "om", "r") == "a")
    {
      cdo_warning("Could not check whether chunk is suitable to be appended since you piped operators.\n"
                  "          Note that the operator did not check 1. for time gaps and 2. for the max size.");
      cdo_print("Output mode: (A)pend.");
      return 1;
    }
    return 0;
  }

  CdiStreamID streamID2 = streamOpenRead(cdo_get_stream_name(0));
  int vlistID2 = streamInqVlist(streamID2);
  int taxisID2 = vlistInqTaxis(vlistID2);
  const auto vDateTime2 = taxisInqVdatetime(taxisID2);
  auto firstdate = julianDate_encode(calendar, vDateTime2);

  int fyear, fmonth, dummy;
  cdiDate_decode(vDateTime2.date, &fyear, &fmonth, &dummy);
  if (ifreq == 1 && fyear == old_end_year)
  {
    if (Options::cdoVerbose)
      cdo_print("In checking the last chunk:\n          The years of the end date of the chunk file "
                "and the first date of the working file are the same: '%d'."
                "   Switched to replace mode for this variable.",
                fyear);
    return 0;
  }

  /* Check temporal distance between last chunk date and first file date */
  double append_distance = julianDate_to_seconds(julianDate_sub(firstdate, juloend)) / 3600.0;
  if ((ifreq == 6 && (append_distance > 2.0 || append_distance < 0))
      || (ifreq == 5 && (append_distance > 6.0 || append_distance < 0))
      || (ifreq == 4 && (append_distance > 12.0 || append_distance < 0))
      || (ifreq == 3 && (append_distance > 48.0 || append_distance < 0))
      || (ifreq == 2 && (append_distance / 24.0 > 62.0 || append_distance < 0))
      || (ifreq == 1 && (append_distance / 24.0 / 30.5 > 24.0 || append_distance < 0))
      || (ifreq == 1 && (append_distance / 24.0 / 30.5 < 1.0))
      || (ifreq == 8 && (append_distance / 24.0 / 30.5 / 12.0 < 1.0 || append_distance < 0))
      || (ifreq == 8 && (append_distance / 24.0 / 30.5 / 12.0 > 20.0)))
  {
    if (Options::cdoVerbose)
      cdo_print("In checking the last chunk:\n          A temporal gap is diagnosed between end date of chunk file "
                "and first date of working file of: '%f' hours ( '%f' days, '%f' months, '%f' years)"
                ". Maximal valid gaps are:\n"
                "          2 hours for 1-hourly frequency\n"
                "          6 hours for 3-hourly frequency\n"
                "          12 hours for 6-hourly frequency\n"
                "          48 hours for daily frequency\n"
                "          62 days for monthly frequency\n"
                "          24 months for yearly frequency\n"
                "          20 years for decadal frequency\n"
                "Minimal valid gaps are:\n"
                "          1 month for yearly frequency\n"
                "          1 year for decadal frequency\n"
                "          Switched to replace mode for this variable.",
                append_distance, append_distance / 24.0, append_distance / 24.0 / 30.5, append_distance / 24.0 / 30.5 / 12.0);
    streamClose(streamID2);
    return 0;
  }
  else if (Options::cdoVerbose)
    cdo_print("The temporal gap between end date of chunk file and first date of working file is '%f' and therefore valid.",
              append_distance);

  if (Options::cdoVerbose) cdo_print("Successfully checked temporal gap.");
  /* Check file size */
  if (Options::cdoVerbose) cdo_print("Start to check file size of chunk + working file.");
  double old_interval_sec = julianDate_to_seconds(julianDate_sub(juloend, julostart));
  double size_per_sec = (double) filesize / old_interval_sec;

  int maxsizegb = std::stoi(kv_get_a_val(kvl, "ms", "2"));
  int maxsizeb = maxsizegb * 1024 * 1024 * 1024;

  int ntsteps = vlistNtsteps(vlistID2);
  if (ntsteps < 0)
  {
    ntsteps = 0;
    while (streamInqTimestep(streamID2, ntsteps++));
    if (ntsteps == 0)
    {
      if (Options::cdoVerbose)
        cdo_print("In checking whether append mode is possible:\n          No time steps found in infile.\n         "
                  " Switched to replace mode for this variable.");
      streamClose(streamID2);
      return 0;
    }
  }

  double estimated_size;
  switch (ifreq)
  {
    case (6): estimated_size = ntsteps * 60 * 60 * 1 * size_per_sec + (double) filesize; break;
    case (5): estimated_size = ntsteps * 60 * 60 * 3 * size_per_sec + (double) filesize; break;
    case (4): estimated_size = ntsteps * 60 * 60 * 6 * size_per_sec + (double) filesize; break;
    case (3): estimated_size = ntsteps * 60 * 60 * 24 * size_per_sec + (double) filesize; break;
    case (2): estimated_size = ntsteps * 60 * 60 * 24 * 30.5 * size_per_sec + (double) filesize; break;
    case (1): estimated_size = ntsteps * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize; break;
    case (8): estimated_size = ntsteps * 10 * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize; break;
    default:
    {
      if (Options::cdoVerbose)
        cdo_print("In checking whether append mode is valid:\n          Selected chunk to append data has subdaily frequency "
                  "which is yet not enabled by cdo cmor.\n          Switched to replace mode for this variable.");
      streamClose(streamID2);
      return 0;
    }
  }

  if (maxsizeb != 0 && (unsigned int) estimated_size > (unsigned int) maxsizeb)
  {
    if (Options::cdoVerbose)
      cdo_print("In checking whether append mode is valid:\n          Estimated file size of appended file is : '%f'gb and "
                "exceeds maximal allowed file size: '%d'gb.\n          Switched to replace mode for this variable.",
                estimated_size / 1024.0 / 1024.0 / 1024.0, maxsizegb);
    streamClose(streamID2);
    return 0;
  }
  streamClose(streamID2);
  if (Options::cdoVerbose) cdo_print("Successfully checked file size of chunk + working file.");
  return 1;
}

static char *
use_chunk_des_files(KVList *kvl, int vlistID, int /*var_id*/, char *chunk_des_file, int ifreq, int calendar)
{
  (void) vlistID;
  char *chunk_file = (char *) std::malloc(4096 * sizeof(char));
  if (file_exist(std::string(chunk_des_file), false, "chunk_description", false))
  {
    auto fp = std::fopen(chunk_des_file, "r");
    if (fp == nullptr) cdo_abort("Could not open chunk description file '%s'", chunk_des_file);
    ListBuffer listBuffer;
    auto status = listBuffer.read(fp, chunk_des_file);
    if (status) cdo_abort("Read error on chunk_description %s!", chunk_des_file);

    char *cbuffer = (char *) (listBuffer.buffer.data());
    size_t pos = 0;
    for (pos = 0; pos < listBuffer.buffer.size(); pos++)
      if (cbuffer[pos] == '\n')
      {
        cbuffer[pos] = 0;
        break;
      }
    if (pos == listBuffer.buffer.size())
    {
      cbuffer[pos - 1] = 0;
      pos--;
    }
    std::snprintf(chunk_file, pos + 1, "%s", cbuffer);

    if (file_exist(chunk_file, false, "chunk_description", false))
    {
      if (check_append_and_size(kvl, chunk_file, ifreq, calendar))
        return chunk_file;
      else
      {
        if (Options::cdoVerbose)
          cdo_print("In checking the last chunk:\n          Chunk '%s' configured via chunk description file is not suitable "
                    "to be appended.",
                    chunk_file);
        cdo_print("Output mode: (R)eplace.");
      }
    }
    else
    {
      if (chunk_file[0])
      {
        if (Options::cdoVerbose)
          cdo_print("In checking the last chunk:\n          Chunk '%s' configured via chunk description file does not exist.",
                    chunk_file);
        cdo_print("Output mode: (R)eplace.");
      }
      else
      {
        if (Options::cdoVerbose) cdo_print("In checking the last chunk:\n          No name found in chunk description file.");
        cdo_print("Output mode: (R)eplace.");
      }
    }
  }
  else
  {
    if (Options::cdoVerbose) cdo_print("Chunk description file '%s' could not be opened.", chunk_des_file);
  }
  strcpy(chunk_file, " \0");
  return chunk_file;
}
static char **
empty_array(struct mapping const vars[], char ***chunk_files)
{
  for (int i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i) (*chunk_files)[i] = nullptr;
  return *chunk_files;
}

static char **
get_chunk_des_files(KVList *kvl, struct mapping vars[], std::string const &miptab_freqptr, int nreq, int vlistID, char *charname,
                    std::string const &project_id)
{
  char **chunk_des_files = (char **) std::malloc((nreq + 1) * sizeof(char *));
  chunk_des_files[nreq] = nullptr;

  char trunk[CMOR_MAX_STRING];
  if (project_id == "CMIP6")
    std::snprintf(trunk, CMOR_MAX_STRING, "%s_", kv_get_a_val(kvl, "source_id", "").c_str());
  else
    std::snprintf(trunk, CMOR_MAX_STRING, "%s_", kv_get_a_val(kvl, "model_id", "").c_str());
  const char *description_atts[] = { "experiment_id", "member", "sub_experiment_id", nullptr };
  strcpy(trunk, miptab_freqptr.c_str());
  for (int i = 0; description_atts[i]; ++i)
  {
    std::strcat(trunk, "_");
    std::strcat(trunk, kv_get_a_val(kvl, description_atts[i], "").c_str());
  }

  for (int j = 0; vars[j].cmor_varID != CMOR_UNDEFID; ++j)
  {
    char *name = (char *) std::malloc(CDI_MAX_NAME * sizeof(char));
    if (charname)
    {
      if (std::strcmp(charname, " ") != 0)
        strcpy(name, charname);
      else
        vlistInqVarName(vlistID, vars[j].cdi_varID, name);
    }
    else
      vlistInqVarName(vlistID, vars[j].cdi_varID, name);

    chunk_des_files[j] = (char *) std::malloc(CMOR_MAX_STRING);
    std::snprintf(chunk_des_files[j], CMOR_MAX_STRING, ".CHUNK_FILE_%s_%s.txt", name, trunk);
    std::free(name);
  }
  return chunk_des_files;
}

static char **
get_chunk_files(KVList *kvl, struct mapping vars[], int vlistID, int ifreq, int time_axis, int calendar,
                std::string const &miptab_freqptr, std::string const &project_id, int *mergeIDs, int psID)
{
  int i = 0;
  for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i);
  if (mergeIDs[0] != CMOR_UNDEFID) i = 1;
  char **chunk_files = (char **) std::malloc((i + 1) * sizeof(char *));
  chunk_files[i] = nullptr;

  if (Options::cdoVerbose) cdo_print("10.2.1. Start to validate append mode.");
  if (kv_get_a_val(kvl, "om", "a") != "a")
    return empty_array(vars, &chunk_files);
  else if (time_axis == 4)
  {
    if (Options::cdoVerbose)
      cdo_print("In validating append mode:\n          CMOR APPEND mode not possible for time independent "
                "variables.");
    cdo_print("Output mode: (R)eplace.");
    return empty_array(vars, &chunk_files);
  }
  if (Options::cdoVerbose) cdo_print("10.2.1. Successfully validated append mode.");

  if (Options::cdoVerbose) cdo_print("10.2.2. Start to get chunk names.");
  int num_aaf = 0;
  std::vector<std::string> chunk_att_files = kv_get_vals(kvl, "lc", &num_aaf);
  char **chunk_des_files = nullptr;
  if (num_aaf != i && num_aaf > 0)
  {
    if (Options::cdoVerbose)
      cdo_print("Number of chunk files '%d' disagree with number of requested variables '%d'.\n Switched to replace mode.\n",
                num_aaf, i);
    cdo_print("Output mode: (R)eplace.");
    return empty_array(vars, &chunk_files);
  }
  else if (num_aaf == 0)
  {
    /* For chunk description file : */
    if (kv_get_a_val(kvl, "d", "y") == "y")
      chunk_des_files = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, nullptr, project_id);
    else
    {
      if (Options::cdoVerbose)
        cdo_print("In getting chunk names:\n          Automatic chunk configuration via file not possible if DRS is "
                  "not created.");
      cdo_print("Output mode: (R)eplace.");
      return empty_array(vars, &chunk_files);
    }
  }
  if (Options::cdoVerbose) cdo_print("10.2.2. Successfully retrieved chunk names.");

  for (int j = 0; vars[j].cmor_varID != CMOR_UNDEFID; ++j)
  {
    if (vars[j].cmor_varID == psID and vars[j].help_var)
    {
      if (Options::cdoVerbose) cdo_print("Chunkfile for ps which was registered for hybrid axis is skipped.");
      continue;
    }
    if (num_aaf != 0)
    {
      if (file_exist(chunk_att_files[j], false, "chunk file", false)
          && check_append_and_size(kvl, (char *) chunk_att_files[j].c_str(), ifreq, calendar))
        chunk_files[j] = strdup(chunk_att_files[j].c_str());
      else
      {
        if (Options::cdoVerbose) cdo_print("Chunk '%s' could not be used.", chunk_att_files[j]);
        cdo_print("Output mode: (R)eplace.");
        chunk_files[j] = strdup(" ");
      }
    }
    else
    {
      if (Options::cdoVerbose)
        cdo_print("It is tried to open a chunk description file for varID: '%d': '%s'.", vars[j].cdi_varID, chunk_des_files[j]);
      chunk_files[j] = use_chunk_des_files(kvl, vlistID, vars[j].cdi_varID, chunk_des_files[j], ifreq, calendar);
      printf("%d %s\n", j, chunk_files[j]);
    }
    if (std::strcmp(chunk_files[j], " ") != 0)
      cdo_print("Output mode: (A)ppend.\n          (Chunk file for var ID %d is: '%s')", vars[j].cdi_varID, chunk_files[j]);
  }
  if (chunk_des_files) free_array(chunk_des_files);
  if (Options::cdoVerbose) cdo_print("Successfully processed chunk file retrieval.");
  return chunk_files;
}

static void
write_variables(KVList *kvl, CdoStreamID streamID, struct mapping vars[], int miptab_freq, int time_axis, int calendar,
                std::string const &miptab_freqptr, std::string const &project_id, int *mergeIDs)
{
  int cmf = 0;
  int vlistID = cdo_stream_inq_vlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tsID = 0;
  int numFields;
  size_t gridsize = cdo_vlist_gridsizemax(vlistID);

  if (Options::cdoVerbose) cdo_print("10. Start to write variables via cmor_write.");
  if (Options::cdoVerbose) cdo_print("10.1. Start to get frequency.");
  int time_unit;
  auto sDateTime = get_taxis(kv_get_a_val(kvl, "rtu", ""), &time_unit);
  int tunitsec = get_tunitsec(time_unit);
  auto ref_date = julianDate_encode(calendar, sDateTime);
  std::string frequency = "";
  if (time_axis != 4)
  {
    frequency = get_frequency(/*kvl,*/ vlistID, miptab_freq);
    if (Options::cdoVerbose && (frequency != "no")) cdo_print("10.1. Successfully retrieved frequency %s.", frequency);
  }
  else if (Options::cdoVerbose)
    cdo_print("10.1. Successfully retrieved fixed time axis.");

  int ifreq = 0;
  if (!frequency.empty())
  {
    if (frequency == "yr") ifreq = 1;
    if (frequency == "mon") ifreq = 2;
    if (frequency == "day") ifreq = 3;
    if (frequency.find("6hr") != std::string::npos) ifreq = 4;
    if (frequency.find("3hr") != std::string::npos) ifreq = 5;
    if (frequency.find("1hr") != std::string::npos) ifreq = 6;
    if (frequency.find("sem") != std::string::npos) ifreq = 7;
    if (frequency.find("dec") != std::string::npos) ifreq = 8;
    if (frequency.find("clim") != std::string::npos) ifreq = 9;
    if (frequency.find("subhr") != std::string::npos) ifreq = 10;
  }

  std::string zaxis = kv_get_a_val(kvl, "za", "notSet");
  bool l_hh = (zaxis == "hybrid_height");
  int ps_index = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
  int psID = getRegisteredPsid(vars, ps_index);

  if (Options::cdoVerbose) cdo_print("10.2. Start to get chunk files.");
  char **chunk_files = nullptr;
  if (ifreq > 0 && ifreq != 7)
    chunk_files = get_chunk_files(kvl, vars, vlistID, ifreq, time_axis, calendar, miptab_freqptr, project_id, mergeIDs, psID);
  else
  {
    int number = 0;
    while (vars[number].cmor_varID != CMOR_UNDEFID) number++;
    chunk_files = (char **) std::malloc((number + 1) * sizeof(char *));
    empty_array(vars, &chunk_files);
  }
  if (ifreq == 7) cdo_print("10.2. Append mode not possible for frequency '%s'. Switch to replace mode.", frequency);
  if (Options::cdoVerbose) cdo_print("10.2. Successfully retrieved chunk files.");
  int i = 0;
  if (chunk_files[0] && chunk_files[0][0] != ' ' && kv_get_a_val(kvl, "sc", "n") == "y")
  {
    while (chunk_files[i])
    {
      char command[CDI_MAX_NAME];
      std::snprintf(command, CDI_MAX_NAME, "cp %s %s.save", chunk_files[i], chunk_files[i]);
      int dir_err = system(command);
      if (dir_err != 0) cdo_warning("Could not create a .save file out of the previous chunk '%s'.", chunk_files[i]);
      i++;
    }
  }
  i = 0;

  int zaxisID, zsize = 0, pscheck = 1;
  char charname[CDI_MAX_NAME] = " ";
  CdoStreamID newstreamID = nullptr;
  for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
    if (vars[i].charvars)
    {
      if (Options::cdoVerbose) cdo_print("10.3. Start to get auxiliary variables.");
      zaxisID = vlistInqVarZaxis(vlistID, vars[i].cdi_varID);
      zsize = zaxisInqSize(zaxisID);
      vlistInqVarName(vlistID, vars[i].cdi_varID, charname);

      cdo_stream_close(streamID);
      newstreamID = cdo_open_read(0);

      vlistID = cdo_stream_inq_vlist(newstreamID);
      taxisID = vlistInqTaxis(vlistID);

      pscheck = 0;
      if (Options::cdoVerbose) cdo_print("10.3. Successfully retrieved auxiliary variables.");
      break;
    }
  if (pscheck == 0)
    if (Options::cdoVerbose)
      cdo_print("Since you defined a variable with character coordinate axis you cannot write another variable with zaxis "
                "of type ZAXIS_HYBRID.");
  if (!newstreamID) newstreamID = streamID;

  int fsize = 0;
  int *mergeIdx = nullptr;
  if (mergeIDs[0] != CMOR_UNDEFID)
  {
    while (mergeIDs[fsize] != CMOR_UNDEFID) fsize++;
    ;
    if (Options::cdoVerbose) cdo_print("10.3. '%d' Variables will be merged.", fsize);
    zaxisID = vlistInqVarZaxis(vlistID, mergeIDs[0]);
    zsize = zaxisInqSize(zaxisID);

    mergeIdx = (int *) std::malloc(fsize * sizeof(int));
    mergeIdx[0] = -1;
    for (int j = 0; j < fsize; j++)
    {
      for (i = 0; vars[i].cdi_varID != CDI_UNDEFID; ++i)
      {
        if (vars[i].cdi_varID == mergeIDs[j]) mergeIdx[j] = i;
      }
    }
    if (mergeIdx[0] == -1) cdo_abort("Could not find registered CMOR varID.");
  }

  if (Options::cdoVerbose) cdo_print("10.4. Start to loop over time steps.");
  while ((numFields = cdo_stream_inq_timestep(newstreamID, tsID++)))
  {
    double time_bnds[2];
    double *time_bndsp = nullptr;
    JulianDate jtime_val;
    double time_val;
    if (time_axis != 4)
    {
      jtime_val = get_cmor_time_val(kvl, taxisID, ref_date, tunitsec, calendar, frequency, tsID, time_axis);
      time_val = julianDate_to_seconds(julianDate_sub(jtime_val, ref_date)) / tunitsec;
      time_bndsp = (time_axis != 1) ? get_time_bounds(kvl, taxisID, ifreq, ref_date, jtime_val, calendar, tunitsec, time_bnds,
                                                      time_axis, vlistID)
                                    : 0;
    }
    while (numFields--) read_record(newstreamID, vars, vlistID);
    if (mergeIDs[0] != CMOR_UNDEFID)
    {
      void *dataslice;
      if (vars[mergeIdx[0]].datatype == 'd')
        dataslice = (void *) std::malloc(gridsize * zsize * fsize * sizeof(double));
      else
        dataslice = (void *) std::malloc(gridsize * zsize * fsize * sizeof(float));

      for (i = 0; i < fsize; i++)
      {
        for (int j = 0; j < (int) gridsize * zsize; ++j)
        {
          if (miptab_freq == 8)
          {
            if (vars[mergeIdx[0]].datatype == 'd')
              ((double *) dataslice)[j * fsize + i] = ((double *) vars[mergeIdx[i]].data)[j];
            else
              ((float *) dataslice)[j * fsize + i] = ((float *) vars[mergeIdx[i]].data)[j];
          }
          else
          {
            if (vars[mergeIdx[0]].datatype == 'd')
              ((double *) dataslice)[j + i * gridsize * zsize] = ((double *) vars[mergeIdx[i]].data)[j];
            else
              ((float *) dataslice)[j + i * gridsize * zsize] = ((float *) vars[mergeIdx[i]].data)[j];
          }
        }
      }
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
      cmf = cmor_write(vars[mergeIdx[0]].cmor_varID, dataslice, vars[mergeIdx[0]].datatype, 1, &time_val, time_bndsp, NULL);
#else
      cmf = cmor_write(vars[mergeIdx[0]].cmor_varID, dataslice, vars[mergeIdx[0]].datatype, chunk_files[0], 1, &time_val,
                       time_bndsp, NULL);
#endif
      std::free(dataslice);
    }
    for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i)
    {
      /*          char name[CDI_MAX_NAME];
                vlistInqVarName(vlistID, vars[i].cdi_varID, name); */
      if (!vars[i].help_var)
      {
        if (time_axis != 4)
        {
          if (vars[i].charvars)
          {
            void *dataslice;
            if (vars[i].datatype == 'd')
            {
              dataslice = (void *) std::malloc(gridsize * zsize * sizeof(double));
              for (int j = 0; j < (int) gridsize * zsize; ++j)
                ((double *) dataslice)[j] = ((double *) vars[i].data)[(tsID - 1) * gridsize * zsize + j];
            }
            else
            {
              dataslice = (void *) std::malloc(gridsize * zsize * sizeof(float));
              for (int j = 0; j < (int) gridsize * zsize; ++j)
                ((float *) dataslice)[j] = ((float *) vars[i].data)[(tsID - 1) * gridsize * zsize + j];
            }
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
            cmf = cmor_write(vars[i].cmor_varID, dataslice, vars[i].datatype, 1, &time_val, time_bndsp, nullptr);
#else
            cmf = cmor_write(vars[i].cmor_varID, dataslice, vars[i].datatype, chunk_files[i], 1, &time_val, time_bndsp, nullptr);
#endif
            std::free(dataslice);
          }
          else if (vars[i].cdi_varID != mergeIDs[0])
          {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)

            cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, 1, &time_val, time_bndsp, nullptr);
#else
            cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, chunk_files[i], 1, &time_val, time_bndsp, nullptr);
#endif
          }
          if (vars[i].zfactor_id > 0 && !l_hh)
          {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
            cmf = cmor_write(vars[i].zfactor_id, vars[psID].data, vars[psID].datatype, 1, &time_val, time_bndsp,
                             &vars[i].cmor_varID);
#else
            cmf = cmor_write(vars[i].zfactor_id, vars[psID].data, vars[psID].datatype, chunk_files[i], 1, &time_val, time_bndsp,
                             &vars[i].cmor_varID);
#endif
          }
        }
        else
        {
#if (CMOR_VERSION_MAJOR == 3 && CMOR_VERSION_MINOR <= 2 && CMOR_VERSION_PATCH <= 7)
          cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, 0, 0, 0, nullptr);
#else
          cmf = cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype, chunk_files[i], 0, 0, 0, nullptr);
#endif
        }
      }
      if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_write failed!", cdo_get_stream_name(0));
    }
  }

  if (Options::cdoVerbose) cdo_print("10.4. Successfully looped over time steps.");
  if (Options::cdoVerbose) cdo_print("10. Successfully written variables via cmor_write.");
  if (Options::cdoVerbose) cdo_print("11. Start to close files, free allocated memory and, if necessary, write chunk files.");
  char **chunkdf = nullptr;
  if (mergeIDs[0] != CMOR_UNDEFID) i = 1;
  if (kv_get_a_val(kvl, "d", "y") == "y")
    chunkdf = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, charname, project_id);

  char file_name[CMOR_MAX_STRING];
  for (i = 0; vars[i].cmor_varID != CMOR_UNDEFID; ++i)
  {
    if (!vars[i].help_var)
    {
      cmf = cmor_close_variable(vars[i].cmor_varID, file_name, nullptr);
      if (std::strcmp(file_name, "") == 0) cdo_abort("ERROR (infile: '%s')! Function cmor_write failed!", cdo_get_stream_name(0));
#if (CMOR_VERSION_MAJOR == 2)
      if (project_id == "CORDEX")
      {
        if (kv_get_a_val(kvl, "tp", "y") == "y")
        {
          if (Options::cdoVerbose) cdo_print("11.1. Start to set a prefix for tracking_id.");
          int ncid, status;
          status = nc_open((const char *) file_name, NC_WRITE, &ncid);
          status = nc_redef(ncid);
          char prelim[CMOR_MAX_STRING];
          status = nc_get_att_text(ncid, NC_GLOBAL, "tracking_id", prelim);
          std::string prefixCordex = "hdl:21.14103/";
          int lengthCombi = std::strlen(prelim) + std::strlen(prefixCordex);
          std::string track;
          std::snprintf(track, lengthCombi, "%s%s", prefixCordex, prelim);
          status = nc_put_att_text(ncid, NC_GLOBAL, "tracking_id", (size_t) lengthCombi, (const char *) track);
          status = nc_enddef(ncid);
          status = nc_close(ncid);
          if (status != NC_NOERR) cdo_abort("ERROR (infile: '%s')! Could not set a prefix for tracking_id", cdo_get_stream_name(0));
          if (Options::cdoVerbose) cdo_print("11.1. Successfully set a prefix for tracking_id.");
        }
      }
#endif
      /*
                    if ( strcmp(kv_get_a_val(kvl, "tracking_prefix", "n"), "y") == 0 )
                      {
                        if (Options::cdoVerbose) cdo_print("11.1. Start to set a prefix for tracking_id.");
                        CdiStreamID streamIDF = streamOpenRead(file_name);
                        int vlistIDF = streamInqVlist(streamIDF);
                        std::string prelim = get_txtatt(vlistIDF, CDI_GLOBAL, "tracking_id");
                        std::string prefixCordex = strdup("21.14103/");
                        size_t lengthCombi = (size_t) (std::strlen(prelim) + std::strlen(prefixCordex));
                        std::string track = (char *) std::malloc( lengthCombi *sizeof(char));
                        sprintf(track, "%s%s", prefixCordex, prelim);
                        cdiDefAttTxt(vlistIDF, CDI_GLOBAL, "tracking_id", lengthCombi, (const char *)track);
                        streamClose(streamIDF);
                        if (Options::cdoVerbose) cdo_print("11.1. Successfully set a prefix for tracking_id.");
                      }
      */
      bool isCordexName = false;
      char cordex_file_name[CMOR_MAX_STRING];
      std::string cordexdir = kv_get_a_val(kvl, "cordexDir", "");
      if ((project_id == "CORDEX") && !cordexdir.empty() && !kv_get_a_val(kvl, "cordexFileTem", "").empty())
      {
        char varname[CMOR_MAX_STRING], timename[CMOR_MAX_STRING];
        char *dummy = file_name;
        int count = 0, firsts = 0, lasts = 0;
        while (file_name[count])
        {
          if (file_name[count] == '_')
          {
            if (firsts == 0) firsts = count;

            lasts = count;
          }
          count++;
        }
        strncpy(varname, file_name, firsts);
        varname[firsts] = '\0';
        dummy += lasts - 1;
        /* Check for CMOR-bug */
        while (*dummy != '_')
        {
          if (*dummy == '-') break;
          dummy--;
        }
        if (*dummy != '_')
        {
          while (*dummy != '_') dummy--;
          std::strcpy(timename, dummy);
          lasts = 1;
          dummy++;
          while (*dummy != '_')
          {
            lasts++;
            dummy++;
          }
          timename[lasts] = '.';
          timename[lasts + 1] = 'n';
          timename[lasts + 2] = 'c';
          timename[lasts + 3] = '\0';
        }
        /* end check for CMOR-bug */
        else
        {
          dummy = file_name;
          dummy += lasts;
          std::strcpy(timename, dummy);
        }

        if (ifreq == 7)
        {
          char smon1[12], smon2[12];
          std::memcpy(smon1, &timename[5], 2);
          std::memcpy(smon2, &timename[12], 2);
          smon1[2] = '\0';
          smon2[2] = '\0';
          if (atol(smon1) != 1)
            std::snprintf(smon1, sizeof(smon1), "%02d", atoi(smon1) - 1);
          else
          {
            char syr[12];
            std::memcpy(syr, &timename[1], 4);
            syr[4] = '\0';
            std::snprintf(syr, sizeof(syr), "%04d", atoi(syr) - 1);
            std::memcpy(&timename[1], syr, 4);

            std::snprintf(smon1, sizeof(smon1), "12");
          }
          if (atol(smon2) != 1)
            std::snprintf(smon2, sizeof(smon2), "%02d", atoi(smon2) + 1);
          else
          {
            char syr[12];
            std::memcpy(syr, &timename[8], 4);
            syr[4] = '\0';
            std::snprintf(syr, sizeof(syr), "%04d", atoi(syr) - 1);
            std::memcpy(&timename[8], syr, 4);

            std::snprintf(smon2, sizeof(smon2), "12");
          }
          std::memcpy(&timename[5], smon1, 2);
          std::memcpy(&timename[12], smon2, 2);
        }

        int cmdlen = 11 + kv_get_a_val(kvl, "cordexDir", "").length() + std::strlen(varname);
        std::vector<char> command1(cmdlen);
        std::snprintf(command1.data(), cmdlen, "mkdir -p %s/%s", kv_get_a_val(kvl, "cordexDir", "").c_str(), varname);

        int dir_err = system(command1.data());
        if (dir_err != 0)
        {
          cdo_warning("Could not create CORDEX compliant path for output files of cdo cmor. Files are created "
                      "in current working directory.");
        }

        std::snprintf(cordex_file_name, CMOR_MAX_STRING, "%s/%s/%s_%s%s", kv_get_a_val(kvl, "cordexDir", "").c_str(), varname,
                      varname, kv_get_a_val(kvl, "cordexFileTem", "").c_str(), timename);

        cmdlen = 5 + std::strlen(file_name) + std::strlen(cordex_file_name);
        std::vector<char> command2(cmdlen);

        std::snprintf(command2.data(), cmdlen, "mv %s %s", file_name, cordex_file_name);
        dir_err = system(command2.data());
        if (dir_err != 0)
        {
          cdo_warning("Could not move cdo cmor output file to CORDEX compliant path.");
          cdo_print("     File stored in:  '%s' with cmor!", file_name);
          if (Options::silentMode) cdo_warning("     File stored in:  '%s' with cmor!", file_name);
        }
        else
        {
          isCordexName = true;
          cdo_print("     File stored in:  '%s' with cmor!", cordex_file_name);
          if (Options::silentMode) cdo_warning("     File stored in:  '%s' with cmor!", cordex_file_name);
        }
      }
      else
      {
        std::string realization = kv_get_a_val(kvl, "realization", "");
        if (realization.empty()) realization = kv_get_a_val(kvl, "realization_index", "");

        if (realization[0] == '0' && realization[1])
        {
          char newname[CDI_MAX_NAME], oldmember[CDI_MAX_NAME], newmember[CDI_MAX_NAME], chunkpath[CDI_MAX_NAME],
              oldchunkpath[CDI_MAX_NAME];
          std::snprintf(oldmember, CDI_MAX_NAME, "r%di", std::stoi(realization));
          std::snprintf(newmember, CDI_MAX_NAME, "r%si", realization.c_str());

          std::string startcmp = std::string(file_name);
          int startpattern = 0, lastSlash = 0;
          std::strcpy(chunkpath, file_name);
          int patternlength = std::strlen(oldmember);
          bool oldchunkcopied = false;
          /* member is once in the path, once in the file name */
          while (file_name[startpattern])
          {
            if (startcmp.substr(0, patternlength) == oldmember)
            {
              if (!oldchunkcopied)
                chunkpath[startpattern] = '\0';
              else
                chunkpath[startpattern + std::strlen(newmember) - patternlength] = '\0';
              startcmp.erase(0, patternlength);
              startpattern += patternlength;
              std::snprintf(newname, CDI_MAX_NAME, "%s%s%s", chunkpath, newmember, startcmp.c_str());
              if (!oldchunkcopied)
              {
                std::snprintf(oldchunkpath, CDI_MAX_NAME, "%s%s", chunkpath, oldmember);
                oldchunkcopied = true;
              }
              std::strcpy(chunkpath, newname);
            }
            startcmp.erase(0, 1);
            startpattern++;
          }
          startpattern = 0;
          while (newname[startpattern])
          {
            if (newname[startpattern] == '/') lastSlash = startpattern;
            startpattern++;
          }

          std::strcpy(chunkpath, newname);
          chunkpath[lastSlash] = '\0';
          char command[CDI_MAX_NAME];
          std::snprintf(command, CDI_MAX_NAME, "mkdir -p %s; mv %s %s;", chunkpath, file_name, newname);

          int dir_err = system(command);
          if (dir_err != 0)
          {
            cdo_warning("Could not move cdo cmor output file to a path with realization=0*.");
            cdo_print("     File stored in:  '%s' with cmor!", file_name);
            if (Options::silentMode) cdo_warning("     File stored in:  '%s' with cmor!", file_name);
          }
          else
          {
            cdo_print("     File stored in:  '%s' with cmor!", newname);
            if (Options::silentMode) cdo_warning("     File stored in:  '%s' with cmor!", newname);
          }

          std::snprintf(command, CDI_MAX_NAME, "rmdir %s*;", oldchunkpath);
          if (Options::cdoVerbose) cdo_print("Start to remove wrong data path (r1 instead of r0*).");
          dir_err = system(command);
          if (dir_err != 0)
            if (Options::cdoVerbose) cdo_print("Failed to remove '%s'", oldchunkpath);
        }
        else
        {
          cdo_print("     File stored in:  '%s' with cmor!", file_name);
          if (Options::silentMode) cdo_warning("     File stored in:  '%s' with cmor!", file_name);
        }
      }
      if (chunkdf)
      {
        if (Options::cdoVerbose) cdo_print("11.2. Start to write a chunk description file.");
        auto fp = std::fopen(chunkdf[i], "w+");
        if (fp)
        {
          if (isCordexName)
            std::fprintf(fp, "%s\n", cordex_file_name);
          else
            std::fprintf(fp, "%s\n", file_name);
        }
        else
        {
          cdo_print("Could not open a chunk description file '%s'.", chunkdf[i]);
          continue;
        }
        std::fclose(fp);
        if (Options::cdoVerbose) cdo_print("11.2. Successfully written a chunk description file '%s'.", chunkdf[i]);
      }
    }
  }
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_close_variable failed!", cdo_get_stream_name(0));
  if (mergeIDs[0] != CMOR_UNDEFID and mergeIdx) std::free(mergeIdx);
  cdo_stream_close(newstreamID);
  if (Options::cdoVerbose) cdo_print("11. Successfully closed files and freed allocated memory.");
}
/*
static int find_cmorvar(list_t *tester, std::string cn, std::string  miptabfreq)
{
  KeyValues *kvcharcn = nullptr, *kvmip = nullptr;
  kvcharcn = tester->search("cmor_name");

  if ( kvcharcn )
    {
      if ( strcmp(kvcharcn->values[0], cn) == 0 )
        {
          kvmip = tester->search("project_mip_table");
          if ( kvmip )
            {
              if ( strcmp(kvmip->values[0], miptabfreq) == 0 )
                return 1;
              else
                return 0;
            }
          return 2;
        }
    }
  return 0;
}
*/
static void
read_maptab(KVList *kvl, CdoStreamID streamID, std::string const &miptabfreq, struct mapping vars[])
{
  /***/
  /* Build mapping table from a combination of two attributes if mt does not begin with / and a directory path is given
   */
  /***/
  if (Options::cdoVerbose) cdo_print("5. Start to find, read and apply mapping table.");
  std::string maptab = kv_get_a_val(kvl, "mt", "");
  std::string maptabdir = kv_get_a_val(kvl, "mapping_table_dir", "");
  std::string maptabbuild;
  const KeyValues *kvn = kvl->search("n");
  const KeyValues *kvc = kvl->search("c");
  const KeyValues *kvcn = kvl->search("cn");
  int filetype = cdo_inq_filetype(streamID);

  if (!maptab.empty() && !maptabdir.empty())
    if (maptab[0] != '/') { maptabbuild = maptabdir + "/" + maptab; }
  if (!maptab.empty())
  {
    if (!maptabbuild.empty()) maptab = maptabbuild;
    int vlistID = cdo_stream_inq_vlist(streamID);

    /***/
    /* Parse the table as a fortran namelist wich contains lists (=lines) of keyvalues */
    /***/
    if (Options::cdoVerbose) cdo_print("5.1 Try to read mapping table: '%s'", maptab);
    kv_insert_vals(kvl, "workfile4err", maptab, true, false);
    PMList pml = cdo_parse_cmor_file(maptab, true);
    if (!pml.size())
    {
      cdo_warning("5.1. In parsing the mapping table '%s':\n          Mapping table could not be parsed. Operator "
                  "continues.",
                  maptab);
      return;
    }
    /***/
    /* If a variable selector name or code is given in cmdline, the corresponding variable is picked from Infile and
     * mapped. */
    /* Only the first value of name/code given in the cmdline is processed */
    /* If no variable selector is given, process all variables and map via name and code */
    /***/
    /* However, if the mapping table contains a keyvalue pair for name or code with more than one value, */
    /* the corresponding variable has a character coordinate and requires special treatment */
    /* This is tested once before mapping. If the special variable equals the variable which is to map, */
    /* the special treatment begins with fct addcharvar */
    /***/
    /* Different CMOR variables are built with one model variable. */
    /* Consequently, for one model variable more than one mapping table entry can exist */
    /* As a second identification argument, the mapping table name (miptabfreq) is used */
    /***/
    /* If no variable selector is given in the mapping table, it is assumed that the infile variable is already named
     * like cmor_name */
    /***/

    if (kvn)
    {
      if (Options::cdoVerbose)
        cdo_print("5. No character axis is built from several variables"
                  "(only possible if values for 'name' are provided in the mapping table).");
      if (kvn->nvalues > 1)
        cdo_warning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline "
                    "variable selection key 'name' is processed.",
                    maptab);
      maptab_via_cmd(maptab, pml, kvn->values[0].c_str(), vlistID, "name", kvcn->values[0], miptabfreq, filetype, maptab);
    }
    else if (kvc)
    {
      if (Options::cdoVerbose)
        cdo_print("5. No character axis is built from several variables"
                  "(only possible if values for 'code' are provided in the mapping table).");
      if (kvc->nvalues > 1)
        cdo_warning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline "
                    "variable selection key 'code' is processed.",
                    maptab);
      maptab_via_cmd(maptab, pml, kvc->values[0].c_str(), vlistID, "code", kvcn->values[0], miptabfreq, filetype, maptab);
    }
    else if (kvcn) { maptab_via_cn(maptab, pml, kvcn->values, vlistID, kvcn->nvalues, miptabfreq, filetype, vars, true); }
    else
    {
      if (Options::cdoVerbose)
        cdo_print("5. No character axis is built from several variables"
                  "(only possible if cmor_name is provided in command line).");
      for (int varID = 0; varID < vlistNvars(vlistID); ++varID)
      {
        /***/
        /* Begin with Code in case infile is of type GRB */
        /***/
        if (filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_GRB2)
          if (maptab_via_key(maptab, pml, vlistID, varID, "code", miptabfreq))
          {
            if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via code.", varID);
            continue;
          }
        if (maptab_via_key(maptab, pml, vlistID, varID, "name", miptabfreq))
        {
          if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via name.", varID);
          continue;
        }
        if (maptab_via_key(maptab, pml, vlistID, varID, "code", miptabfreq))
        {
          if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via code.", varID);
          continue;
        }
        /***/
        /* In case corresponding mapping table entry does not contain a variable selector attribute */
        /***/
        if (maptab_via_key(maptab, pml, vlistID, varID, "cmor_name", miptabfreq))
        {
          if (Options::cdoVerbose) cdo_print("5.1. Successfully mapped varID '%d' via cmor_name.", varID);
          continue;
        }
        cdo_warning("5.1. In applying the mapping table '%s':\n          Could not map variable with id '%d'.", maptab, varID);
      }
    }
    /***/
    /* In case a requested variable needs an auxilliary variable, the latter may be mapped later. */
    /* If a mapping table exists is saved here */
    /***/
    kv_insert_vals(kvl, "mtproof", maptab, true, false);
    cdo_print("Mapping Table = '%s'.", maptab);
  }
  else if (Options::cdoVerbose)
    cdo_print("5. No mapping table found.");
}

static void
replace_key(KVList *kvl, const KeyValues &kv, std::string const &newkey)
{
  std::vector<std::string> values(kv.nvalues);
  int k = 0;
  for (k = 0; k < kv.nvalues; k++) values[k] = kv.values[k];
  kvl->remove(kv.key);
  kvl->append(newkey, values, k);
}

static void
parse_cmdline(KVList *kvl, std::vector<std::string> &params)
{
  /* Already set params++ in main function */
  if (kvl->parse_arguments(params) != 0) cdo_abort("ERROR (infile: '%s')! Could not parse command line.", cdo_get_stream_name(0));

  std::vector<KeyValues> keystorm, keystosubs;
  for (auto const &kv : *kvl)
  {
    const std::string short_key = check_short_key(kv.key, true);
    if (!short_key.empty())
    {
      if (kv.key != short_key) keystosubs.push_back(kv);
    }
    else
    {
      cdo_warning("Unknown commandline keyword: '%s'\n", kv.key);
      keystorm.push_back(kv);
    }
  }
  for (size_t i = 0; i < keystosubs.size(); ++i) { replace_key(kvl, keystosubs[i], check_short_key(keystosubs[i].key, true)); }
  for (size_t i = 0; i < keystorm.size(); ++i) { kvl->remove(keystorm[i].key); }
}

static std::string
get_mip_table(std::string const &params, KVList *kvl, std::string const &project_id, bool print)
{
  std::string miptab;
  if (print && Options::cdoVerbose) cdo_print("2.2. Start to find a MIP table file.");
  if (params.empty())
    cdo_abort("ERROR (infile: '%s')! First parameter not passed. A MIP table file is required.", cdo_get_stream_name(0));
  if (file_exist(params, false, "MIP table", print))
  {
    miptab = params;
    int j = 0;
    for (size_t i = 0; i < params.length(); ++i)
    {
      if (params.at(i) == '/') j = i;
    }
    char miptabdir[1024];
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    cwd[strlen(cwd)] = '\0';
    if (params.at(0) == '/')
    {
      strncpy(miptabdir, params.c_str(), j + 1);
      miptabdir[j + 1] = '\0';
    }
    else if (j == 0)
      std::strcpy(miptabdir, cwd);
    else
    {
      std::strcpy(miptabdir, cwd);
      std::strcat(miptabdir, "/");
      strncat(miptabdir, params.c_str(), j + 1);
      miptabdir[strlen(cwd) + j + 1] = '\0';
    }
    kv_insert_vals(kvl, "mip_table_dir", std::string(miptabdir), true, false);
    if (print) cdo_print("MIP table file = '%s'.", miptab);
    return miptab;
  }
  else
  {
    if (print && Options::cdoVerbose)
      cdo_print("Try to build a path with additional configuration attributes:\n          'mip_table_dir' and "
                "'project_id'\n          in order to use '%s' as MIP-table.",
                params);
    std::string miptabdir = kv_get_a_val(kvl, "mip_table_dir", "");
    if (!miptabdir.empty() && !project_id.empty())
    {
#if (CMOR_VERSION_MAJOR == 2)
      {
        miptab = miptabdir + "/" + project_id + "_" + params;
      }
#elif (CMOR_VERSION_MAJOR == 3)
      {
        miptab = miptabdir + "/" + project_id + "_" + params + ".json";
      }
#endif
      file_exist(miptab, true, "MIP table", print);
      if (print) cdo_print("MIP table file = '%s'", miptab);
      return miptab;
    }
    else
      cdo_abort("ERROR (infile: '%s')! In finding the MIP table:\n          Could not find attribute 'mip_table_dir'.",
                cdo_get_stream_name(0));
  }

  return miptab;
}

static std::string
freq_from_path(std::string const &mip_table)
{
  std::string freq;

  // Find the position of the last '/' character
  size_t lastSlashPos = mip_table.find_last_of('/');

  if (lastSlashPos != std::string::npos) freq = mip_table.substr(lastSlashPos + 1);

  size_t lastUnderscorePos = mip_table.find_last_of('_');
  if ((lastUnderscorePos != std::string::npos) && lastUnderscorePos > lastSlashPos) freq = mip_table.substr(lastUnderscorePos + 1);

  return freq;
}

static int
get_miptab_freq(std::string const &mip_table, std::string const &project_id)
{
  int miptab_freq = 0;
  std::string freq = freq_from_path(mip_table);
  if (!freq.empty())
  {
    if (freq.find("yr") != std::string::npos || freq.find("Yr") != std::string::npos)
      miptab_freq = 11;
    else if (freq.find("mon") != std::string::npos || freq.find("Mon") != std::string::npos)
      miptab_freq = 12;
    else if (freq.find("day") != std::string::npos || freq.find("Day") != std::string::npos)
      miptab_freq = 13;
    else if (freq.find("6h") != std::string::npos)
      miptab_freq = 14;
    else if (freq.find("3h") != std::string::npos)
      miptab_freq = 15;
    else if (freq.find("1h") != std::string::npos || freq.find("AERhr") != std::string::npos)
      miptab_freq = 16;
    else if (freq.find("sem") != std::string::npos)
      miptab_freq = 17;
    else if (freq.find("dec") != std::string::npos)
      miptab_freq = 18;
    else if (freq.find("subhr") != std::string::npos)
      miptab_freq = 19;

    if (freq == "Oclim")
      miptab_freq = 1;
    else if (freq == "Oyr")
      miptab_freq = 2;
    else if (freq == "cfMon")
      miptab_freq = 3;
    else if (freq == "day")
      miptab_freq = 4;
    else if ((freq == "6hrPlev") && (project_id == "CMIP5"))
      miptab_freq = 5;
    else if (freq == "6hrPlevPt")
      miptab_freq = 5;
    else if (freq == "6hrLev")
      miptab_freq = 6;
    else if (freq == "E1hrClimMon")
      miptab_freq = 7;
    else if (freq == "E3hrPt")
      miptab_freq = 8;
  }
  return miptab_freq;
}

static void
check_cmdline_mapping(KVList *kvl)
{
  std::string name = kv_get_a_val(kvl, "n", "");
  std::string code = kv_get_a_val(kvl, "c", "");
  std::string cn = kv_get_a_val(kvl, "cn", "");
  if (!name.empty() && !code.empty())
    cdo_abort("ERROR (infile: '%s')! Mapping via command line failed. Only one variable selector of 'name' and 'code' is allowed.",
              cdo_get_stream_name(0));
  if ((!name.empty() && cn.empty()) || (!code.empty() && cn.empty()))
    cdo_abort("ERROR (infile: '%s')! Mapping via command line failed. A corresponding 'cmor_name' is needed.",
              cdo_get_stream_name(0));
}

static std::string
get_project_id(KVList *kvl, std::string const &params)
{
  if (Options::cdoVerbose) cdo_print("2.1. Start to check whether 'project_id' or 'mip_era' is denoted.");
  std::string project_id = "", dummy, dummy2;
  dummy = kv_get_a_val(kvl, "project_id", "");
  dummy2 = kv_get_a_val(kvl, "mip_era", "");

  char tester[CDI_MAX_NAME];
  std::strcpy(tester, params.c_str());
  char *testerpointer = tester;
  int underscore = 0, slash = 0, testint = 0;
  while (tester[testint])
  {
    if (tester[testint] == '_') underscore = testint;
    if (tester[testint] == '/') slash = testint;
    testint++;
  }
  if (underscore > slash)
  {
    if (slash) testerpointer += slash + 1;
    testerpointer[underscore - slash - 1] = '\0';
  }
#if defined(CMOR_VERSION_MAJOR)
#if (CMOR_VERSION_MAJOR == 2)
  {
    if (dummy.empty() && dummy2.empty())
    {
      if (testerpointer != tester)
      {
        if (Options::cdoVerbose)
          cdo_print("Could not find attribute 'project_id'.\n          "
                    "Try to use substring from MIP-table input '%s' as project_id.",
                    testerpointer);
        project_id = std::string(testerpointer);
      }
      else
        cdo_abort("ERROR (infile: '%s')! Attribute 'project_id' is required.", cdo_get_stream_name(0));
    }
    else if (dummy.empty())
      cdo_abort("ERROR (infile: '%s')! Cannot produce CMIP6 standard with CMOR2.\n          "
                "Value for attribute 'project_id' is required.",
                cdo_get_stream_name(0));
    else
      project_id = dummy;
  }
#elif (CMOR_VERSION_MAJOR == 3)
  {
    if (dummy.empty() && dummy2.empty())
    {
      if (std::strcmp(testerpointer, tester) != 0)
      {
        if (Options::cdoVerbose)
          cdo_print("Could not find attribute 'project_id'.\n          "
                    "Try to use substring from MIP-table input '%s' as project_id.",
                    testerpointer);
        project_id = std::string(testerpointer);
      }
      else
        cdo_abort("ERROR (infile: '%s')! Attribute 'mip_era' or 'project_id' is required.", cdo_get_stream_name(0));
    }
    else if (dummy2.empty())
    {
      if (Options::cdoVerbose)
        cdo_print("You have not provided 'mip_era' but only 'project_id'."
                  " If you try to produce CMIP5 standard,\n          It is recommended to use CMOR2 for this job instead.");
      project_id = dummy;
    }
    else
      project_id = dummy2;
  }
#endif
#else
  cdo_abort("ERROR (infile: '%s')! Cannot check CMOR version: Missing makro CMOR_VERSION_MAJOR", cdo_get_stream_name(0));
#endif

  if (Options::cdoVerbose) cdo_print("2.1. Successfully found project_id / mip_era: '%s'.", project_id);
  return project_id;
}

static int
cmor_load_and_set_table(KVList *kvl, std::string const &param0, std::string const &project_id, std::string &mip_table)
{
  int table_id = 0, cmf = 0;
#if (CMOR_VERSION_MAJOR == 3)
  mip_table = get_mip_table(param0, kvl, project_id, false);
#endif
  cmf = cmor_load_table((char *) mip_table.c_str(), &table_id);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_load_table failed!", cdo_get_stream_name(0));
  cmf = cmor_set_table(table_id);
  if (cmf != 0) cdo_abort("ERROR (infile: '%s')! Function cmor_set_table failed!", cdo_get_stream_name(0));
  return table_id;
}

class CMOR : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "CMOR",
    .operators = { { "cmor", CmorHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<CMOR> registration = RegisterEntry<CMOR>();

  CdoStreamID streamID;
  KVList kvl;
  int vlistID;

  int miptab_freq;

  std::string project_id;
  std::string miptab_freqptr;
  std::string miptableInput;
  std::string mip_table;

  struct mapping *vars;

public:
  void
  init() override
  {
    int nparams = cdo_operator_argc();
    if (nparams < 1) cdo_abort("ERROR (infile: '%s')! No parameter found. Need at least a MIP-table.", cdo_get_stream_name(0));

    auto params = cdo_get_oper_argv();
    miptableInput = params[0];
    // copy everything from params except the first
    if (nparams > 1)
    {
      /* Define cmdline list and read cmdline */
      params = std::vector<std::string>(params.begin() + 1, params.end());
      parse_cmdline(&kvl, params);
    }

    /* Check whether a command line mapping is active */
    check_cmdline_mapping(&kvl);

    /* Config files are read with descending priority. */
    read_config_files(&kvl);

    /* Get project_id, mip_table and mip_table frequency*/
    if (Options::cdoVerbose) cdo_print("2. Start to find a MIP table and to deduce a frequency from MIP table file.");
    project_id = get_project_id(&kvl, miptableInput);
    mip_table = get_mip_table(miptableInput, &kvl, project_id, true);
#if (CMOR_VERSION_MAJOR == 3)
    mip_table.erase(mip_table.length() - 5);
#endif
    miptab_freq = get_miptab_freq(mip_table, project_id);

    miptab_freqptr = freq_from_path(mip_table);
    kv_insert_vals(&kvl, "miptab_freq", miptab_freqptr, true, false);

    if (Options::cdoVerbose)
      cdo_print("2. Successfully found a MIP table '%s' and deduced a MIP table frequency '%s'.", mip_table, miptab_freqptr);

    if (Options::cdoVerbose) cdo_print("3. Start to open infile '%s'.", cdo_get_stream_name(0));
    streamID = cdo_open_read(0);
    vlistID = cdo_stream_inq_vlist(streamID);

    if (Options::cdoVerbose) cdo_print("3. Successfully opened infile '%s'.", cdo_get_stream_name(0));
  }

  void
  run() override
  {
    if (Options::cdoVerbose) cdo_print("4. Start to check attributes.");
    /* Short keys from rtu, mt, gi must be included similar to global atts */
    add_globalhybrids(&kvl, vlistID);
    /* Allow time units from infile */
    check_required_time_units(&kvl, vlistInqTaxis(vlistID));

    /* Check for attributes and member name */
    check_attr(&kvl, project_id, vlistID);
    check_mem(&kvl, project_id);
    if (Options::cdoVerbose) cdo_print("4. Successfully checked global attributes.");

    /* dump_global_attributes(pml, streamID); */

    vars = construct_var_mapping(vlistID);

    /* read mapping table */
    read_maptab(&kvl, streamID, miptab_freqptr, vars);

    int time_axis = 0, calendar = 0;

    setup_dataset(&kvl, streamID, &calendar, project_id);

    int table_id = cmor_load_and_set_table(&kvl, miptableInput, project_id, mip_table);

    int mergeIDs[150];
    mergeIDs[0] = CMOR_UNDEFID;
    register_all_dimensions(&kvl, streamID, vars, table_id, project_id, miptab_freq, &time_axis, mergeIDs);
    write_variables(&kvl, streamID, vars, miptab_freq, time_axis, calendar, miptab_freqptr, project_id, mergeIDs);
  }

  void
  close() override
  {
    destruct_var_mapping(vars);
    /* std::free(miptableInput); */
  }
};
#endif  // HAVE_LIBCMOR
