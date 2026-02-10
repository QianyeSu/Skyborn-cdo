/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef PMLIST_H
#define PMLIST_H

#include <list>
#include <vector>
#include <string>

#include "listbuffer.h"
#include "namelist.h"
#include "cdo_cmor.h"

struct KeyValues
{
  int nvalues{ 0 };
  std::string key;
  std::vector<std::string> values;
};

// clang-format off
class  // KVList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
KVList : public std::list<KeyValues>
// clang-format on
{
public:
  std::string name;

  void print(FILE *fp = stderr) const;
  int parse_arguments(std::vector<std::string> const &argv);
  const KeyValues *search(std::string const &key) const;
  void remove(std::string const &inkey);
  void append(std::string const &key, std::vector<std::string> const &values, int nvalues);
  const std::string get_first_value(std::string const &key, std::string const &replacer);
};

// clang-format off
class  // PMList
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
PMList : public std::list<KVList>
// clang-format on
{
public:
  const KVList *searchKVListVentry(std::string const &key, std::string const &value, std::vector<std::string> const &entry);
  const KVList *getKVListVentry(std::vector<std::string> const &entry);
  void print(FILE *fp = stderr);
  void read_namelist(FILE *fp, const char *name);
  void read_cmor_table(FILE *fp, const char *name);
};

void mapvar(int vlistID, int varID, const KeyValues &kv, CmorVar &cmorVar, bool &hasValidMin, bool &hasValidMax, int ptab,
            bool isnPtmodeName);
int parse_namelist(PMList &pmlist, NamelistParser &parser, char *buf, bool cdocmor);
int parse_list_buffer(NamelistParser &p, ListBuffer &listBuffer);

#endif
