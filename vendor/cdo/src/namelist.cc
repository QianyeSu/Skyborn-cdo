/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdio>
#include <cctype>

#include "namelist.h"

// Allocates a fresh unused token from the token pull.
NamelistToken *
NamelistParser::allocToken()
{
  if (this->toknext >= this->num_tokens)
  {
    constexpr unsigned long TOK_MEM_INCR = 1024;
    this->num_tokens += TOK_MEM_INCR;
    this->tokens.resize(this->num_tokens);
  }

  auto tok = &this->tokens[this->toknext++];
  tok->start = tok->end = -1;
  return tok;
}

// Fills token type and boundaries.
void
NamelistToken::fill(NamelistType _type, long _start, long _end)
{
  this->type = _type;
  this->start = _start;
  this->end = _end;
}

void
NamelistParser::newObject()
{
  auto token = this->allocToken();
  token->type = NamelistType::OBJECT;
  token->start = this->pos;
}

// Fills next available token with NAMELIST word.
NamelistError
NamelistParser::parseWord(const char *buf, size_t len)
{
  auto start = this->pos;

  for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
  {
    switch (buf[this->pos])
    {
      case ':':
      case '=':
      case ',':
      case '&':
      case '/':
      case '\r':
      case '\n':
      case '\t':
      case ' ': goto found;
    }

    if (buf[this->pos] < 32 || buf[this->pos] >= 127)
    {
      this->pos = start;
      return NamelistError::INVAL;
    }
  }

found:

  auto token = this->allocToken();
  token->fill(NamelistType::WORD, start, this->pos);
  this->pos--;

  return NamelistError::UNDEFINED;
}

// Fills next token with NAMELIST string.
NamelistError
NamelistParser::parseString(const char *buf, size_t len, char quote)
{
  auto start = this->pos;

  this->pos++;

  // Skip starting quote
  for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
  {
    auto c = buf[this->pos];

    // Quote: end of string
    if (c == quote)
    {
      auto token = this->allocToken();
      token->fill(NamelistType::STRING, start + 1, this->pos);
      return NamelistError::UNDEFINED;
    }

    // Backslash: Quoted symbol expected
    if (c == '\\' && this->pos + 1 < len)
    {
      this->pos++;
      switch (buf[this->pos])
      {
        // Allowed escaped symbols
        case '\"':
        case '\\':
        case 'b':
        case 'f':
        case 'r':
        case 'n':
        case 't': break;
        // Allows escaped symbol \uXXXX
        case 'u':
          this->pos++;
          for (long i = 0; i < 4 && this->pos < len && buf[this->pos] != '\0'; ++i)
          {
            // If it isn't a hex character we have an error
            if (!((buf[this->pos] >= 48 && buf[this->pos] <= 57) ||  // 0-9
                  (buf[this->pos] >= 65 && buf[this->pos] <= 70) ||  // A-F
                  (buf[this->pos] >= 97 && buf[this->pos] <= 102)))  // a-f
            {
              return NamelistError::INVAL;
            }
            this->pos++;
          }
          this->pos--;
          break;
        // Unexpected symbol
        default: return NamelistError::INVAL;
      }
    }
  }

  this->pos = start;

  return NamelistError::PART;
}

NamelistError
NamelistParser::checkKeyname(const char *buf, NamelistToken *t)
{
  switch (t->type)
  {
    case NamelistType::STRING:
      while (std::isspace((int) buf[t->start]) && t->start < t->end) t->start++;
      while (std::isspace((int) buf[t->end - 1]) && t->start < t->end) t->end--;
      if ((t->end - t->start) < 1) return NamelistError::EMKEY;
      for (long i = t->start; i < t->end; ++i)
        if (std::isspace((int) buf[i])) return NamelistError::INKEY;
      t->type = NamelistType::KEY;
      break;
    case NamelistType::WORD: t->type = NamelistType::KEY; break;
    default: return NamelistError::INTYP;
  }

  return NamelistError::UNDEFINED;
}

NamelistError
NamelistParser::parse(const char *buf, size_t len)
{
  auto status = NamelistError::UNDEFINED;

  this->lineno = 1;

  for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
  {
    auto c = buf[this->pos];
    switch (c)
    {
      case '&': this->newObject(); break;
      case '/':
        for (long i = this->toknext - 1; i >= 0; i--)
        {
          auto token = &this->tokens[i];
          if (token->start != -1 && token->end == -1)
          {
            if (token->type != NamelistType::OBJECT) return NamelistError::INOBJ;
            token->end = this->pos + 1;
            break;
          }
        }
        break;
      case '\t':
      case ' ': break;
      case '\r':
        if (this->pos + 1 < len && buf[this->pos + 1] == '\n') this->pos++;
        this->lineno++;
        break;
      case '\n': this->lineno++; break;
      case ',': break;
      case '#':
      case '!':  // Skip to end of line
        for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
          if (buf[this->pos] == '\r' || buf[this->pos] == '\n')
          {
            this->pos--;
            break;
          }
        break;
      case ':':
      case '=': status = this->checkKeyname(buf, &this->tokens[this->toknext - 1]); break;
      case '\"':
      case '\'': status = this->parseString(buf, len, c); break;
      default: status = this->parseWord(buf, len); break;
    }

    if (status != NamelistError::UNDEFINED) return status;
  }

  return status;
}

void
NamelistParser::dump(const char *buf)
{
  auto ntok = this->toknext;
  printf("Number of tokens %lu\n", ntok);

  for (unsigned long it = 0; it < ntok; ++it)
  {
    auto t = &this->tokens[it];
    int length = t->end - t->start;
    auto start = buf + t->start;
    printf("Token %lu", it + 1);
    if (t->type == NamelistType::OBJECT)
    {
      printf(" NAMELIST=");
      if (length > 80) length = 80;
      printf("'%.*s'", length, start);
    }
    else if (t->type == NamelistType::KEY)
    {
      printf(" KEY=");
      printf("'%.*s'", length, start);
    }
    else if (t->type == NamelistType::WORD)
    {
      printf(" WORD=");
      printf("'%.*s'", length, start);
    }
    else if (t->type == NamelistType::STRING)
    {
      printf(" STRING=");
      printf("'%.*s'", length, start);
    }
    printf("\n");
  }
}

int
NamelistParser::verify()
{
  auto ntok = this->toknext;

  if (ntok)
  {
    auto t = &this->tokens[0];
    if (t->type != NamelistType::OBJECT && t->type != NamelistType::KEY) return -1;
  }

  return 0;
}
