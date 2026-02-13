/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstring>
#include <cstdlib>
#include <vector>

#include "compare.h"
#include "pmlist.h"
#include "namelist.h"
#include "cdo_output.h"
#include "util_string.h"

static void
keyValuesPrint(std::FILE *fp, const KeyValues &keyValues)
{
  std::fprintf(fp, "  %s =", keyValues.key.c_str());
  for (int i = 0; i < keyValues.nvalues; ++i) std::fprintf(fp, " '%s'", keyValues.values[i].c_str());
  std::fprintf(fp, "\n");
}

void
KVList::print(std::FILE *fp) const
{
  for (auto const &keyValues : *this) keyValuesPrint(fp, keyValues);
}

int
KVList::parse_arguments(std::vector<std::string> const &argv)
{
  /*  this function assumes input in the form of
   *  value_name1=10,2,3,4,value_name2=201,23
   **/

  size_t equalPos = argv[0].find('=');
  if (equalPos == std::string::npos)
  {
    std::fprintf(stderr, "missing '=' in key/value string: >%s<\n", argv[0].c_str());
    return -1;
  }

  for (auto &it : argv)
  {
    equalPos = it.find('=');
    if (equalPos != std::string::npos)
    {
      auto key = it.substr(0, equalPos);

      std::string value = Util::String::trim(it.substr(equalPos + 1));
      {  // value vector gets moved at the end of scope
        std::vector<std::string> valueVector;
        value.empty() ? valueVector = {} : valueVector = { value };

        push_back(KeyValues());
        KeyValues &newValue = back();
        newValue.nvalues = 0;
        newValue.values = std::move(valueVector);
        newValue.key = key;
      }
      back().nvalues = back().values.size();
    }
    else
    {
      auto &last = back();
      last.values.push_back(it);
      last.nvalues = last.values.size();
    }
  }

  return 0;
}

const KeyValues *
KVList::search(std::string const &key) const
{
  for (auto const &kv : *this)
  {
    if (kv.key == key) return &kv;
  }

  return nullptr;
}

const std::string
KVList::get_first_value(std::string const &key, std::string const &replacer)
{
  auto kv = this->search(key);
  if (kv && kv->nvalues > 0) return kv->values[0];
  return replacer;
}

void
KVList::append(std::string const &key, std::vector<std::string> const &values, int nvalues)
{
  KeyValues kv;
  kv.key = key;
  kv.nvalues = nvalues;
  kv.values.resize(nvalues);
  for (int i = 0; i < nvalues; ++i) kv.values[i] = values[i];
  this->push_back(kv);
}

// Remove only one list item
void
KVList::remove(std::string const &inkey)
{
  std::list<KeyValues>::iterator i;
  for (i = this->begin(); i != this->end(); ++i)
    if (i->key == inkey) break;
  if (i->key == inkey) this->erase(i);
}

const KVList *
PMList::searchKVListVentry(std::string const &key, std::string const &value, std::vector<std::string> const &entry)
{
  for (auto const &kvlist : *this)
  {
    for (auto const &s : entry)
      if (kvlist.name == s)
      {
        auto kv = kvlist.search(key);
        if (kv && kv->nvalues > 0 && kv->values[0] == value) return &kvlist;
      }
  }

  return nullptr;
}

const KVList *
PMList::getKVListVentry(std::vector<std::string> const &entry)
{
  for (auto const &kvlist : *this)
  {
    for (auto const &s : entry)
      if (kvlist.name == s) return &kvlist;
  }

  return nullptr;
}

static void
KVListAppendNamelist(KVList &kvlist, const char *key, const char *buffer, NamelistToken *t, int nvalues, bool cdocmor)
{
  std::vector<char> value;
  KeyValues kv;
  kv.key = key;
  kv.nvalues = nvalues;
  if (nvalues > 0) kv.values.resize(nvalues);

  for (int i = 0; i < nvalues; ++i)
  {
    size_t len = t[i].end - t[i].start;
    /** CMOR_MAX_STRING cannot be used **/
    if (cdocmor && len > 1024) cdo_abort("A string value is larger than the maximum size allowed by CMOR (1024 signs).");

    if (value.size() < (len + 1)) value.resize(len + 1);
    auto pval = buffer + t[i].start;

    if (cdocmor && kv.key == "code" && !std::strstr(pval, ","))
    {
      value.resize(4);
      auto code = atol(pval);
      if (code > 0 && code < 1000)
        std::snprintf(value.data(), 4, "%03ld", code);
      else
        cdo_warning("In parsing a line of a file:\n          "
                    "Codes could not be transformed into the code format (three digit integer). Codes wont be used.");
    }
    else
    {
      std::memcpy(value.data(), pval, len);
      value[len] = 0;
    }

    kv.values[i] = value.data();
  }

  kvlist.push_back(std::move(kv));
}

static unsigned long
get_number_of_values(const unsigned long ntok, const NamelistToken *tokens)
{
  unsigned long it = 0;
  for (; it < ntok; ++it)
  {
    auto type = tokens[it].type;
    if (type != NamelistType::WORD && type != NamelistType::STRING) break;
  }

  return it;
}

static void
replace_name(char *name)
{
  for (size_t pos = 0; pos < std::strlen(name); pos++) name[pos] = std::tolower(name[pos]);

  if (cdo_cmpstr(name, "conventions")) std::strcpy(name, "Conventions");
  if (cdo_cmpstr(name, "cn")) std::strcpy(name, "cmor_name");
  if (cdo_cmpstr(name, "c")) std::strcpy(name, "code");
  if (cdo_cmpstr(name, "n")) std::strcpy(name, "name");
  if (cdo_cmpstr(name, "pmt")) std::strcpy(name, "project_mip_table");
  if (cdo_cmpstr(name, "cordex_domain")) std::strcpy(name, "CORDEX_domain");
  if (cdo_cmpstr(name, "char_axis_landuse")) std::strcpy(name, "char_axis_landUse");
  if (cdo_cmpstr(name, "_formula_var_file")) std::strcpy(name, "_FORMULA_VAR_FILE");
  if (cdo_cmpstr(name, "_axis_entry_file")) std::strcpy(name, "_AXIS_ENTRY_FILE");
}

int
parse_namelist(PMList &pmlist, NamelistParser &parser, char *buf, bool cdocmor)
{
  char name[4096];
  KVList kvlist;
  auto &tokens = parser.tokens;
  auto ntok = parser.toknext;

  for (unsigned long it = 0; it < ntok; ++it)
  {
    auto const &t = tokens[it];
    // printf("Token %u", it+1);
    if (t.type == NamelistType::OBJECT)
    {
      name[0] = 0;
      if (it + 1 < ntok && tokens[it + 1].type == NamelistType::WORD)
      {
        it++;
        auto const &t2 = tokens[it];
        std::snprintf(name, sizeof(name), "%.*s", (int) (t2.end - t2.start), buf + t2.start);
        name[sizeof(name) - 1] = 0;
      }

      if (kvlist.size())
      {
        pmlist.push_back(kvlist);
        kvlist.clear();
      }

      kvlist.name = name;
    }
    else if (t.type == NamelistType::KEY)
    {
      // printf(" key >%.*s<\n", t.end - t.start, buf+t.start);
      std::snprintf(name, sizeof(name), "%.*s", (int) (t.end - t.start), buf + t.start);
      name[sizeof(name) - 1] = 0;
      auto nvalues = get_number_of_values(ntok - it - 1, &tokens[it + 1]);

      if (cdocmor)
      {
        if (nvalues == 0)
        {
          cdo_warning("Could not find values for key '%s'.", name);
          continue;
        }
        replace_name(name);
      }

      KVListAppendNamelist(kvlist, name, buf, &tokens[it + 1], nvalues, cdocmor);
      it += nvalues;
    }
    else
    {
      // printf(" token >%.*s<\n", (int)(t.end - t.start), buf+t.start);
      break;
    }
  }

  if (kvlist.size()) pmlist.push_back(kvlist);

  return 0;
}

int
parse_list_buffer(NamelistParser &p, ListBuffer const &listBuffer)
{
  const char *errMsg = "Namelist error";
  const auto name = listBuffer.name.c_str();

  auto status = p.parse(listBuffer.buffer.data(), listBuffer.buffer.size());
  if (status != NamelistError::UNDEFINED)
  {
    switch (status)
    {
      case NamelistError::INVAL:
        std::fprintf(stderr, "%s: Invalid character in %s (line=%lu character='%c' dec=%u)!\n", errMsg, name, p.lineno,
                     listBuffer.buffer[p.pos], (unsigned char) listBuffer.buffer[p.pos]);
        break;
      case NamelistError::PART:
        std::fprintf(stderr, "%s: End of string not found in %s (line=%lu)!\n", errMsg, name, p.lineno);
        break;
      case NamelistError::INKEY: std::fprintf(stderr, "%s: Invalid keyword in %s (line=%lu)!\n", errMsg, name, p.lineno); break;
      case NamelistError::INTYP:
        std::fprintf(stderr, "%s: Invalid keyword type in %s (line=%lu)!\n", errMsg, name, p.lineno);
        break;
      case NamelistError::INOBJ: std::fprintf(stderr, "%s: Invalid object in %s (line=%lu)!\n", errMsg, name, p.lineno); break;
      case NamelistError::EMKEY: std::fprintf(stderr, "%s: Empty key name in %s (line=%lu)!\n", errMsg, name, p.lineno); break;
      default: std::fprintf(stderr, "%s in %s (line=%lu)!\n", errMsg, name, p.lineno); break;
    }
    cdo_abort("%s!", errMsg);
  }

  // p.dump(listBuffer.buffer.data());
  if (p.verify())
  {
    std::fprintf(stderr, "%s: Invalid contents in %s!\n", errMsg, name);
    cdo_abort("Namelist error!");
  }

  return 0;
}

void
PMList::read_namelist(std::FILE *fp, std::string_view name)
{
  ListBuffer listBuffer;
  if (listBuffer.read(fp, name.data())) cdo_abort("Read error on namelist %s!", name);

  NamelistParser p;
  auto status = parse_list_buffer(p, listBuffer);
  if (status) cdo_abort("Namelist not found!");

  parse_namelist(*this, p, listBuffer.buffer.data(), false);
}

void
PMList::print(std::FILE *fp)
{
  for (auto const &kvlist : *this)
  {
    std::fprintf(fp, "\nFound %s list with %zu key/values: \n", kvlist.name.c_str(), kvlist.size());
    kvlist.print(fp);
  }
}
