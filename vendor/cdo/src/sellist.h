/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef SELLIST_H
#define SELLIST_H

#include "pmlist.h"

enum struct SelType
{
  UNDEF,
  INT,
  FLT,
  WORD
};

struct SelectEntry
{
  int nvalues = 0;
  std::string key;
  std::vector<std::string> values;
  std::vector<bool> flag;
  std::string description;
  SelType type = { SelType::UNDEF };
  std::vector<int> ivalues;
  std::vector<double> dvalues;
  std::vector<const char *> cvalues;
};

class SelectInfo
{
private:
  void init(const KVList &kvlist);

public:
  std::vector<SelectEntry> selList;

  explicit SelectInfo(const KVList &kvlist) { init(kvlist); }

  void verify() const;
  void print() const;
  int nvalues(int listIdx) const;
  bool
  isValidListIdx(int listIdx) const noexcept
  {
    return (listIdx >= 0 && listIdx < (int) selList.size());
  }
};

#define SELINFO_ADD_INT(name, description) \
  int name = 0;                            \
  int listIdx_##name = selinfo_add(selInfo, description, #name, SelType::INT)
#define SELINFO_ADD_FLT(name, description) \
  double name = 0;                         \
  int listIdx_##name = selinfo_add(selInfo, description, #name, SelType::FLT)
#define SELINFO_ADD_WORD(name, description) \
  const char *name = nullptr;               \
  int listIdx_##name = selinfo_add(selInfo, description, #name, SelType::WORD)
#define SELINFO_NVAL(name) selInfo.nvalues(listIdx_##name)
#define SELINFO_CHECK_FLAG(name) selinfo_check_flag(selInfo, listIdx_##name)
#define SELINFO_CHECK_RANGE_FLAG(name) selinfo_check_range_flag(selInfo, listIdx_##name)
#define SELINFO_CHECK(name) selinfo_check(selInfo, listIdx_##name, &name)
#define SELINFO_CHECK_INDEX(name, numValues) selinfo_check_index(selInfo, listIdx_##name, name, numValues)
#define SELINFO_CHECK_DATE(name) selinfo_check_date(selInfo, listIdx_##name, name)
#define SELINFO_CHECK_SEASON(name, month) selinfo_check_season(selInfo, listIdx_##name, month)
#define SELINFO_CHECK_RANGE(name, value) selinfo_check_range(selInfo, listIdx_##name, value)
#define SELINFO_DEF_FLAG(name, valIdx, flag) selinfo_def_flag(selInfo, listIdx_##name, valIdx, flag)
#define SELINFO_GET_VAL(name, valIdx, val) selinfo_get_val(selInfo, listIdx_##name, valIdx, val)
#define SELINFO_DEF_VAL(name, valIdx, val) selinfo_def_val(selInfo, listIdx_##name, valIdx, val)

int selinfo_add(SelectInfo &selInfo, std::string const &description, std::string const &name, SelType type);
void selinfo_check_flag(SelectInfo const &selInfo, int listIdx);
void selinfo_check_range_flag(SelectInfo const &selInfo, int listIdx);
bool selinfo_check(SelectInfo &selInfo, int listIdx, void *par);
bool selinfo_check_index(SelectInfo &selInfo, int listIdx, int par, int maxValues);
bool selinfo_check_date(SelectInfo &selInfo, int listIdx, const char *par);
bool selinfo_check_season(SelectInfo &selInfo, int listIdx, int month);
bool selinfo_check_range(SelectInfo &selInfo, int listIdx, double value);
void selinfo_def_flag(SelectInfo &selInfo, int listIdx, int valIdx, bool flag);
void selinfo_get_val(SelectInfo const &selInfo, int listIdx, int valIdx, void *val);
void selinfo_def_val(SelectInfo &selInfo, int listIdx, int valIdx, void *val);

#endif
