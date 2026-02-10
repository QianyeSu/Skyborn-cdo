/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cassert>
#include <cmath>

#include "param_conversion.h"
#include "sellist.h"
#include "util_wildcards.h"
#include "util_string.h"
#include "cdo_output.h"

// #define SELDEBUG 1

void
SelectInfo::init(const KVList &kvlist)
{
  selList.resize(kvlist.size());

  int i = 0;
  for (auto const &kv : kvlist)
  {
    auto &e = selList[i];
    e.type = SelType::UNDEF;
    e.key = kv.key;
    e.nvalues = kv.nvalues;
    e.values.resize(kv.nvalues);
    for (int k = 0; k < kv.nvalues; ++k) e.values[k] = kv.values[k];
#ifdef SELDEBUG
    printf("%s =", e.key.c_str());
    for (auto const &value : e.values) printf(" '%s'", value.c_str());
    printf("\n");
#endif
    ++i;
  }

#ifdef SELDEBUG
  for (auto const &e : selList)
  {
    printf("%s =", e.key.c_str());
    for (auto const &value : e.values) printf(" '%s'", value.c_str());
    printf("\n");
  }
#endif
}

void
SelectInfo::verify() const
{
  for (auto const &e : selList)
    if (e.type == SelType::UNDEF) cdo_abort("Unsupported selection keyword: '%s'!", e.key);
}

int
selinfo_add(SelectInfo &selInfo, std::string const &description, std::string const &name, SelType type)
{
  auto &selList = selInfo.selList;
  int listIdx = -1;

  for (int i = 0, n = selList.size(); i < n; ++i)
  {
    if (selList[i].key == name)
    {
      listIdx = i;
      break;
    }
  }

  if (selInfo.isValidListIdx(listIdx))
  {
    auto &e = selList[listIdx];
    e.type = type;
    e.description = description;
    if (e.nvalues)
    {
      switch (type)
      {
        case SelType::INT: e.ivalues.resize(e.nvalues); break;
        case SelType::FLT: e.dvalues.resize(e.nvalues); break;
        case SelType::WORD: e.cvalues.resize(e.nvalues); break;
        case SelType::UNDEF: break;
      }
    }

    int j = 0;
    auto nvalues = e.nvalues;
    for (int i = 0; i < nvalues; ++i)
    {
      switch (type)
      {
        case SelType::INT:
        {
          int first, last, inc;
          split_intstring(e.values[i], first, last, inc);

          if (first == last) { e.ivalues[j++] = first; }
          else
          {
            int k = 0;
            if (inc >= 0)
              for (int ival = first; ival <= last; ival += inc) k++;
            else
              for (int ival = first; ival >= last; ival += inc) k++;

            e.nvalues += k - 1;
            if (e.nvalues)
            {
              e.ivalues.resize(e.nvalues);

              if (inc >= 0)
                for (int ival = first; ival <= last; ival += inc) e.ivalues[j++] = ival;
              else
                for (int ival = first; ival >= last; ival += inc) e.ivalues[j++] = ival;
            }
          }

          break;
        }
        case SelType::FLT: e.dvalues[i] = parameter_to_double(e.values[i]); break;
        case SelType::WORD: e.cvalues[i] = parameter_to_word(e.values[i].c_str()); break;
        case SelType::UNDEF: break;
      }
    }

    if (e.nvalues) e.flag.resize(e.nvalues, false);
#ifdef SELDEBUG
    printf("%s =", e.key.c_str());
    for (int i = 0; i < e.nvalues; ++i)
    {
      switch (type)
      {
        case SelType::INT: printf(" %d", e.ivalues[i]); break;
        case SelType::FLT: printf(" %g", e.dvalues[i]); break;
        case SelType::WORD: printf(" %s", e.cvalues[i]); break;
        case SelType::UNDEF: break;
      }
    }
    printf("\n");
#endif
  }

  return listIdx;
}

int
SelectInfo::nvalues(int listIdx) const
{
  return (isValidListIdx(listIdx) ? selList[listIdx].nvalues : 0);
}

void
selinfo_check_flag(SelectInfo const &selInfo, int listIdx)
{
  auto const &selList = selInfo.selList;
  if (!selInfo.isValidListIdx(listIdx)) return;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    auto const &e = selList[listIdx];
    for (int i = 0; i < nvalues; ++i)
      if (!e.flag[i])
      {
        switch (e.type)
        {
          case SelType::INT: cdo_warning("%s >%d< not found!", e.description, e.ivalues[i]); break;
          case SelType::FLT: cdo_warning("%s >%g< not found!", e.description, e.dvalues[i]); break;
          case SelType::WORD: cdo_warning("%s >%s< not found!", e.description, e.cvalues[i]); break;
          case SelType::UNDEF: break;
        }
      }
  }
}

void
selinfo_check_range_flag(SelectInfo const &selInfo, int listIdx)
{
  auto const &selList = selInfo.selList;
  if (!selInfo.isValidListIdx(listIdx)) return;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues == 2)
  {
    auto const &e = selList[listIdx];
    if (!e.flag[0] && e.type == SelType::FLT) cdo_warning("%s %g to %g not found!", e.description, e.dvalues[0], e.dvalues[1]);
  }
}

bool
selinfo_check(SelectInfo &selInfo, int listIdx, void *par)
{
  auto &selList = selInfo.selList;
  auto found = false;

  if (!selInfo.isValidListIdx(listIdx)) return found;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    auto &e = selList[listIdx];
    switch (e.type)
    {
      case SelType::INT:
      {
        int ival = *static_cast<int *>(par);
        for (int i = 0; i < nvalues; ++i)
        {
          if (ival == e.ivalues[i]) { e.flag[i] = found = true; }
        }
        break;
      }
      case SelType::FLT:
      {
        const double dval = *static_cast<double *>(par);
        for (int i = 0; i < nvalues; ++i)
        {
          if (std::fabs(dval - e.dvalues[i]) < 1.e-4) { e.flag[i] = found = true; }
        }
        break;
      }
      case SelType::WORD:
      {
        const char *cval = *static_cast<char **>(par);
        for (int i = 0; i < nvalues; ++i)
        {
          if (wildcard_match(cval, e.cvalues[i])) { e.flag[i] = found = true; }
        }
        break;
      }
      case SelType::UNDEF: break;
    }
  }

  return found;
}

bool
selinfo_check_index(SelectInfo &selInfo, int listIdx, int ival, int maxValues)
{
  auto &selList = selInfo.selList;
  auto found = false;

  if (!selInfo.isValidListIdx(listIdx)) return found;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    auto &e = selList[listIdx];
    if (e.type == SelType::INT)
    {
      for (int i = 0; i < nvalues; ++i)
      {
        int eval = e.ivalues[i];
        if (eval < 0) eval = maxValues + eval + 1;
        if (ival == eval) { e.flag[i] = found = true; }
      }
    }
  }

  return found;
}

bool
selinfo_check_date(SelectInfo &selInfo, int listIdx, const char *par)
{
  auto &selList = selInfo.selList;
  auto found = false;

  if (!selInfo.isValidListIdx(listIdx)) return found;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    auto &e = selList[listIdx];

    if (*par == ' ') ++par;

    for (int i = 0; i < nvalues; ++i)
    {
      auto wcdate = string_to_upper(e.values[i]) + '*';
      if (wildcard_match(par, wcdate)) { e.flag[i] = found = true; }
    }
  }

  return found;
}

bool
selinfo_check_season(SelectInfo &selInfo, int listIdx, int month)
{
  auto &selList = selInfo.selList;
  assert(month >= 1 && month <= 12);
  auto found = false;

  if (!selInfo.isValidListIdx(listIdx)) return found;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    int imon[13];  // 1-12 !
    auto &e = selList[listIdx];

    for (int i = 0; i < nvalues; ++i)
    {
      for (int m = 0; m < 13; ++m) imon[m] = 0;
      season_to_months(e.values[i], imon);
      if (imon[month]) { e.flag[i] = found = true; }
    }
  }

  return found;
}

bool
selinfo_check_range(SelectInfo &selInfo, int listIdx, double value)
{
  auto &selList = selInfo.selList;
  auto found = false;

  if (!selInfo.isValidListIdx(listIdx)) return found;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues == 2)
  {
    auto &e = selList[listIdx];
    if (e.type == SelType::FLT)
    {
      auto rmin = e.dvalues[0];
      auto rmax = e.dvalues[1];
      if (value >= rmin && value <= rmax) { e.flag[0] = e.flag[1] = found = true; }
    }
  }

  return found;
}

void
selinfo_def_flag(SelectInfo &selInfo, int listIdx, int valIdx, bool flag)
{
  auto &selList = selInfo.selList;
  if (!selInfo.isValidListIdx(listIdx)) return;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues)
  {
    if (valIdx >= 0 && valIdx < nvalues)
    {
      auto &e = selList[listIdx];
      e.flag[valIdx] = flag;
    }
  }
}

void
selinfo_get_val(SelectInfo const &selInfo, int listIdx, int valIdx, void *val)
{
  auto const &selList = selInfo.selList;
  if (!selInfo.isValidListIdx(listIdx)) return;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues && valIdx >= 0 && valIdx < nvalues)
  {
    auto const &e = selList[listIdx];
    switch (e.type)
    {
      case SelType::INT: *static_cast<int *>(val) = e.ivalues[valIdx]; break;
      case SelType::FLT: *static_cast<double *>(val) = e.dvalues[valIdx]; break;
      case SelType::WORD: *static_cast<const char **>(val) = e.cvalues[valIdx]; break;
      case SelType::UNDEF: break;
    }
  }
}

void
selinfo_def_val(SelectInfo &selInfo, int listIdx, int valIdx, void *val)
{
  auto &selList = selInfo.selList;
  if (!selInfo.isValidListIdx(listIdx)) return;

  auto nvalues = selInfo.nvalues(listIdx);
  if (nvalues && valIdx >= 0 && valIdx < nvalues)
  {
    auto &e = selList[listIdx];
    switch (e.type)
    {
      case SelType::INT: e.ivalues[valIdx] = *static_cast<int *>(val); break;
      case SelType::FLT: e.dvalues[valIdx] = *static_cast<double *>(val); break;
      case SelType::WORD: e.cvalues[valIdx] = *static_cast<char **>(val); break;
      case SelType::UNDEF: break;
    }
  }
}

static void
selinfo_print_val(const SelectEntry &e, int valIdx)
{
  switch (e.type)
  {
    case SelType::INT: printf(" %d", e.ivalues[valIdx]); break;
    case SelType::FLT: printf(" %g", e.dvalues[valIdx]); break;
    case SelType::WORD: printf(" %s", e.cvalues[valIdx]); break;
    case SelType::UNDEF: break;
  }
}

void
SelectInfo::print() const
{
  if (selList.size() > 0)
  {
    printf("Num  Name             Type  Size  Entries\n");
    for (int listIdx = 0, n = selList.size(); listIdx < n; ++listIdx)
    {
      auto const &e = selList[listIdx];
      printf("%3d  %-16s %4d  %4d ", listIdx + 1, e.key.c_str(), (int) e.type, e.nvalues);
      auto numValues = e.nvalues;
      if (numValues > 12) numValues = 11;
      for (int valIdx = 0; valIdx < numValues; ++valIdx) selinfo_print_val(e, valIdx);
      if (numValues < e.nvalues)
      {
        printf(" ...");
        selinfo_print_val(e, e.nvalues - 1);
      }
      printf("\n");
    }
  }
}
