/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <vector>
#include <string>

#if __has_include(<wordexp.h>)
#define HAVE_WORDEXP_H 1
#include <wordexp.h>
#endif
#if __has_include(<glob.h>)
#define HAVE_GLOB_H 1
#include <glob.h>
#endif

#include "util_wildcards.h"

#ifdef HAVE_GLOB_H
static int
get_glob_flags(void)
{
  int glob_flags = 0;

#ifdef GLOB_NOCHECK
  glob_flags |= GLOB_NOCHECK;
#endif
#ifdef GLOB_TILDE
  glob_flags |= GLOB_TILDE;
#endif

  return glob_flags;
}
#endif

static int
find_wildcard(const char *text)
{
  auto len = std::strlen(text);
  if (len > 0)
  {
    if (text[0] == '~') { return true; }
    for (size_t i = 0; i < len; ++i)
      if (text[i] == '?' || text[i] == '*' || text[i] == '[') { return true; }
  }

  return false;
}

// used in griddes.cc
char *
expand_filename(const char *fileName)
{
  char *fileNameOut = nullptr;

  if (find_wildcard(fileName))
  {
#ifdef HAVE_GLOB_H
    auto glob_flags = get_glob_flags();
    glob_t glob_results;
    glob(fileName, glob_flags, 0, &glob_results);
    if (glob_results.gl_pathc == 1) fileNameOut = strdup(glob_results.gl_pathv[0]);
    globfree(&glob_results);
#endif
  }

  return fileNameOut;
}

#ifdef HAVE_WORDEXP_H
void
wordexp_error_status(int status, const char *argument)
{
  if (status == 0) return;

  if (status == WRDE_BADCHAR)
  {
    fprintf(stderr,
            "Argument '%s' contains one of the following unsupported unquoted characters: <newline>, `|', "
            "`&', `;', `<', `>', `(', `)', `{', `}'.\n",
            argument);
  }
  else if (status == WRDE_NOSPACE) { fprintf(stderr, "Not enough memory to store the result.\n"); }
  else if (status == WRDE_SYNTAX) { fprintf(stderr, "Shell syntax error in '%s'\n", argument); }
  else if (status == WRDE_BADVAL) { fprintf(stderr, "Undefined shell variable in '%s'\n", argument); }
  else { fprintf(stderr, "wordexp() returns an error.\n"); }

  exit(EXIT_FAILURE);
}

// Expands all input file wildcards and removes the wildcard while inserting all expanded files into argv
std::vector<std::string>
expand_path_names(std::vector<std::string> argv)
{
  for (size_t idx = 0; idx < argv.size(); idx++)
  {
    // if argv[idx] contains wildcard (* or [?]+), multiple ** are ignored
    if (argv[idx].find_first_of("*?[ ") != std::string::npos)
    {
      constexpr int flags = WRDE_UNDEF;
      wordexp_t glob_results;
      auto status = wordexp(argv[idx].c_str(), &glob_results, flags);
      if (status != 0) wordexp_error_status(status, argv[idx].c_str());

      // range based insert (glob_results.we_wordv is inserted before wildcard
      if (std::string(glob_results.we_wordv[0]).find_first_of("*?[ ") == std::string::npos)
      {
        argv.insert(argv.begin() + idx + 1, glob_results.we_wordv, glob_results.we_wordv + glob_results.we_wordc);
        argv.erase(argv.begin() + idx);
      }
      // delete wildcard
      wordfree(&glob_results);
    }
  }

  return argv;
}
#else
std::vector<std::string>
expand_path_names(std::vector<std::string> argv)
{
  return argv;
}
#endif

#ifdef HAVE_WORDEXP_H
// Expands all input file wildcards and removes the wildcard while inserting all expanded files into argv
std::vector<std::string>
expand_wild_cards(std::vector<std::string> argv)
{
  bool applyActive = false;
  int bracketsOpen = 0;
  for (size_t idx = 1; idx < argv.size(); idx++)
  {
    // if argv[idx] contains wildcard (* or [?]+), multiple ** are ignored
    if (argv[idx].compare(0, 6, "-apply") == 0)
    {
      applyActive = true;
      continue;
    }
    if (argv[idx].size() == 1 && argv[idx][0] == '[')
    {
      bracketsOpen++;
      continue;
    }
    if (argv[idx].size() == 1 && argv[idx][0] == ']')
    {
      bracketsOpen--;
      if (bracketsOpen == 0) { applyActive = false; }
      continue;
    }
    if (argv[idx][0] != '-' && argv[idx].find_first_of("*?[ ") != std::string::npos)
    {
      constexpr int flags = WRDE_UNDEF;
      wordexp_t glob_results;
      auto status = wordexp(argv[idx].c_str(), &glob_results, flags);
      if (status != 0) wordexp_error_status(status, argv[idx].c_str());

      // range based insert (glob_results.we_wordv is inserted before wildcard
      if (std::string(glob_results.we_wordv[0]).find_first_of("*?[ ") == std::string::npos)
      {
        auto insertAt = idx + 1;
        if (applyActive == false)
        {
          argv.insert(argv.begin() + insertAt, "]");
          argv.insert(argv.begin() + insertAt, "[");
          insertAt += 1;
        }
        argv.insert(argv.begin() + insertAt, glob_results.we_wordv, glob_results.we_wordv + glob_results.we_wordc);
        argv.erase(argv.begin() + idx);
      }
      // delete wildcard
      wordfree(&glob_results);
    }
  }

  return argv;
}
#else
std::vector<std::string>
expand_wild_cards(std::vector<std::string> argv)
{
  return argv;
}
#endif

// Wild card matching using single traversal from https://www.geeksforgeeks.org/dsa/wildcard-pattern-matching
bool
wildcard_match(std::string const &text, std::string const &pattern)
{
  int n = text.length();
  int m = pattern.length();
  int i = 0, j = 0, startIndex = -1, match = 0;
  while (i < n)
  {
    // Characters match or '?' in pattern matches any character.
    if (j < m && (pattern[j] == '?' || pattern[j] == text[i])) { i++, j++; }
    // Wildcard character '*', mark the current position in the pattern and the text as a proper match.
    else if (j < m && pattern[j] == '*') { startIndex = j++, match = i; }
    // No match, but a previous wildcard was found.
    // Backtrack to the last '*' character position and try for a different match.
    else if (startIndex != -1) { j = startIndex + 1, i = ++match; }
    // If none of the above cases comply, the pattern does not match.
    else { return false; }
  }
  // Consume any remaining '*' characters in the given pattern.
  while (j < m && pattern[j] == '*') { j++; }
  // If we have reached the end of both the pattern and the text, the pattern matches the text.
  return j == m;
}
