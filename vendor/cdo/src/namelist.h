/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef NAMELIST_H_
#define NAMELIST_H_

#include <vector>

enum struct NamelistType
{
  UNDEFINED = 0,
  OBJECT = 1,
  KEY = 2,
  STRING = 3,
  WORD = 4
};

enum struct NamelistError
{
  UNDEFINED = 0,
  INVAL = -1,  // Invalid character inside NAMELIST string/word
  PART = -2,   // The string is not a full NAMELIST packet, more bytes expected
  INKEY = -3,  // Invalid character inside NAMELIST key
  INTYP = -4,  // Invalid NAMELIST key type
  INOBJ = -5,  // Invalid NAMELIST object
  EMKEY = -6   // Empty key name
};

// NAMELIST token description.
class NamelistToken
{
public:
  NamelistType type{};  // type (object, key, string word)
  long start{};         // start position in NAMELIST buffer
  long end{};           // end position in NAMELIST buffer

  void fill(NamelistType type, long start, long end);
};

class NamelistParser
{
public:
  std::vector<NamelistToken> tokens;
  unsigned long num_tokens{};
  unsigned long toknext{};
  unsigned long pos{};
  unsigned long lineno{};

  NamelistError parse(const char *buf, size_t len);
  void dump(const char *buf);
  int verify();

private:
  NamelistError checkKeyname(const char *buf, NamelistToken *t);
  NamelistError parseString(const char *buf, size_t len, char quote);
  NamelistError parseWord(const char *buf, size_t len);
  void newObject();
  NamelistToken *allocToken();
};

#endif  // NAMELIST_H_
