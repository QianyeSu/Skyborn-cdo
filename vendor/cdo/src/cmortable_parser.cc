/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_output.h"
#include "pmlist.h"
#include "json/jsmn.h"

static char *
readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  size_t ipos = 0;

  while (*buffersize)
    {
      int ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if (ichar == '\r')
        {
          if (*buffersize)
            {
              ichar = *buffer;
              if (ichar == '\n')
                {
                  (*buffersize)--;
                  buffer++;
                }
            }
          break;
        }
      if (ichar == '\n') break;
      line[ipos++] = ichar;
      if (ipos >= len)
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %zu)!\n", len);
          break;
        }
    }

  if (ipos < len) line[ipos] = 0;

  if (*buffersize == 0 && ipos == 0) buffer = nullptr;

  return buffer;
}

static char *
skipSeparator(char *pline)
{
  while (std::isspace((int) *pline)) pline++;
  if (*pline == '=' || *pline == ':') pline++;
  while (std::isspace((int) *pline)) pline++;

  return pline;
}

static char *
getElementName(char *pline, char *name)
{
  while (std::isspace((int) *pline)) pline++;
  auto len = std::strlen(pline);
  size_t pos = 0;
  while (pos < len && !std::isspace((int) *(pline + pos)) && *(pline + pos) != '=' && *(pline + pos) != ':') pos++;

  strncpy(name, pline, pos);
  name[pos] = 0;

  pline += pos;
  return pline;
}

static char *
getElementValue(char *pline)
{
  while (std::isspace((int) *pline)) pline++;
  auto len = std::strlen(pline);
  if (*pline != '"' && *pline != '\'')
    for (size_t i = 1; i < len; ++i)
      if (pline[i] == '!')
        {
          pline[i] = 0;
          len = i;
          break;
        }
  while (std::isspace((int) *(pline + len - 1)) && len)
    {
      *(pline + len - 1) = 0;
      len--;
    }

  return pline;
}

static void
parse_cmortable_buffer(PMList &pmlist, size_t buffersize, char *buffer)
{
  char line[4096], name[256];
  const char *listentry[] = { "axis_entry", "variable_entry" };
  int nentry = sizeof(listentry) / sizeof(listentry[0]);
  // int linenumber = 0;
  KVList kvlist;

  while ((buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))))
    {
      // linenumber++;
      auto pline = line;
      while (std::isspace((int) *pline)) pline++;
      if (*pline == '#' || *pline == '!' || *pline == '\0') continue;
      //  len = (int) std::strlen(pline);

      int ientry = -1;
      for (ientry = 0; ientry < nentry; ++ientry)
        if (std::strncmp(pline, listentry[ientry], std::strlen(listentry[ientry])) == 0) break;

      if (ientry < nentry)
        {
          pline += std::strlen(listentry[ientry]);

          if (kvlist.size())
            {
              pmlist.push_back(kvlist);
              kvlist.clear();
            }

          kvlist.name = listentry[ientry];

          pline = skipSeparator(pline);
          pline = getElementValue(pline);

          std::vector<std::string> pppline = { std::string(pline) };
          if (*pline) kvlist.append("name", pppline, 1);
        }
      else
        {
          pline = getElementName(pline, name);
          pline = skipSeparator(pline);
          pline = getElementValue(pline);

          if (kvlist.size() == 0) kvlist.name = "Header";

          std::vector<std::string> pppline = { std::string(pline) };
          if (*pline) kvlist.append(name, pppline, 1);
        }
    }

  if (kvlist.size()) pmlist.push_back(kvlist);
}

// not used
int
dump_json(const char *js, jsmntok_t *t, size_t count, int level)
{
  if (count == 0) return 0;

  if (t->type == JSMN_PRIMITIVE)
    {
      printf("%.*s", t->end - t->start, js + t->start);
      return 1;
    }
  else if (t->type == JSMN_STRING)
    {
      printf("'%.*s'", t->end - t->start, js + t->start);
      return 1;
    }
  else if (t->type == JSMN_OBJECT)
    {
      printf("\n");
      //  printf("Object: size %d\n", t->size);
      printf("Object: size %d count %d level %d\n", t->size, (int) count, level);
      int j = 0;
      for (int i = 0; i < t->size; ++i)
        {
          for (int k = 0; k < level; ++k) printf("  ");
          j += dump_json(js, t + 1 + j, count - j, level + 1);
          printf(": ");
          j += dump_json(js, t + 1 + j, count - j, level + 1);
          printf("\n");
        }
      return j + 1;
    }
  else if (t->type == JSMN_ARRAY)
    {
      int j = 0;
      printf("\n");
      for (int i = 0; i < t->size; ++i)
        {
          for (int k = 0; k < level - 1; ++k) printf("  ");
          printf("   - ");
          j += dump_json(js, t + 1 + j, count - j, level + 1);
          printf("\n");
        }
      return j + 1;
    }

  return 0;
}

static void
KVList_append_json(KVList &kvlist, const char *key, const char *js, jsmntok_t *t, int nvalues)
{
  KeyValues kv;
  kv.key = strdup(key);
  kv.nvalues = nvalues;
  kv.values.resize(nvalues);
  for (int i = 0; i < nvalues; ++i)
    {
      auto len = t[i].end - t[i].start;
      std::vector<char> value(len + 1);
      std::snprintf(value.data(), len + 1, "%.*s", (int) len, js + t[i].start);
      value[len] = 0;
      // printf("set %s: '%s'\n", key, value);
      kv.values[i] = value.data();
    }
  kvlist.push_back(kv);
}

static int
add_tokens_json(PMList &pmlist, const char *js, jsmntok_t *t, int count)
{
  bool debug = false;
  char name[4096];
  int i = 0;
  int nobj = t[0].size;

  if (t[0].type == JSMN_OBJECT)
    {
      KVList kvlist;
      while (nobj--)
        {
          ++i;
          auto pmlname = i;
          if (debug) printf("  object: %.*s\n", t[i].end - t[i].start, js + t[i].start);
          ++i;
          if (t[i].type == JSMN_OBJECT)
            {
              int ic = 0;
            NEXT:
              std::snprintf(name, sizeof(name), "%.*s", t[pmlname].end - t[pmlname].start, js + t[pmlname].start);
              name[sizeof(name) - 1] = 0;
              // printf("new object: %s\n", name);
              if (kvlist.size())
                {
                  pmlist.push_back(kvlist);
                  kvlist.clear();
                }

              kvlist.name = name;

              if (t[i + 2].type == JSMN_OBJECT)
                {
                  if (ic == 0)
                    ic = t[i].size;
                  else
                    ic--;

                  ++i;
                  KVList_append_json(kvlist, "name", js, &t[i], 1);
                  if (debug) printf("    name: '%.*s'\n", t[i].end - t[i].start, js + t[i].start);
                  ++i;
                }
              int n = t[i].size;
              while (n--)
                {
                  ++i;
                  std::snprintf(name, sizeof(name), "%.*s", t[i].end - t[i].start, js + t[i].start);
                  name[sizeof(name) - 1] = 0;
                  if (debug) printf("    %.*s:", t[i].end - t[i].start, js + t[i].start);
                  ++i;
                  if (t[i].type == JSMN_ARRAY)
                    {
                      int nae = t[i].size;
                      KVList_append_json(kvlist, name, js, &t[i + 1], nae);
                      while (nae--)
                        {
                          ++i;
                          if (debug) printf(" '%.*s'", t[i].end - t[i].start, js + t[i].start);
                        }
                    }
                  else
                    {
                      KVList_append_json(kvlist, name, js, &t[i], 1);
                      if (debug) printf(" '%.*s'", t[i].end - t[i].start, js + t[i].start);
                    }
                  if (debug) printf("\n");
                }
              if (ic > 1) goto NEXT;
            }
        }

      if (kvlist.size()) pmlist.push_back(kvlist);
    }

  if (debug) printf("Processed %d of %d tokens!\n", i, count - 1);

  return 0;
}

static void
parse_cmortable_buffer_json(PMList &pmlist, const std::vector<char> &buffer, const char *filename)
{
  // Prepare parser
  auto p = jsmn_new();

  auto status = jsmn_parse(p, buffer.data(), buffer.size());
  if (status != 0)
    {
      switch (status)
        {
        case JSMN_ERROR_INVAL:
          fprintf(stderr, "JSON error: Invalid character in %s (line=%u character='%c')!\n", filename, p->lineno, buffer[p->pos]);
          break;
        case JSMN_ERROR_PART: fprintf(stderr, "JSON error: End of string not found in %s (line=%u)!\n", filename, p->lineno); break;
        default: fprintf(stderr, "JSON error in %s (line=%u)\n", filename, p->lineno); break;
        }
    }

  add_tokens_json(pmlist, buffer.data(), p->tokens, (int) p->toknext);
  jsmn_destroy(p);
}

void
PMList::read_cmor_table(FILE *fp, const char *name)
{
  ListBuffer listBuffer;
  if (listBuffer.read(fp, name)) cdo_abort("Read error on CMOR table %s!", name);

  int buffer0 = listBuffer.buffer[0];

  if (buffer0 == '{') { parse_cmortable_buffer_json(*this, listBuffer.buffer, name); }
  else if (std::strncmp(listBuffer.buffer.data(), "table_id:", 9) == 0)
    {
      parse_cmortable_buffer(*this, listBuffer.buffer.size(), listBuffer.buffer.data());
    }
  else if (buffer0 == '&' || buffer0 == '#')
    {
      NamelistParser p;
      auto status = parse_list_buffer(p, listBuffer);
      if (status) cdo_abort("Namelist not found!");

      parse_namelist(*this, p, listBuffer.buffer.data(), false);
    }
  else
    cdo_abort("Invalid CMOR table (file: %s)!", name);
}
