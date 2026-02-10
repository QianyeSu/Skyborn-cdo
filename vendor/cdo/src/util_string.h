/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef UTIL_STRING_H
#define UTIL_STRING_H

#include <string>
#include <tuple>
#include <vector>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <concepts>

#define ADD_PLURAL(n) ((n) != 1 ? "s" : "")

std::vector<std::string> get_operator_argv(std::string operatorArguments);
std::vector<std::string> split_args(std::string operatorArguments);

std::string getenv_string(std::string const &envVar);

std::vector<std::string> split_string(std::string const &str, std::string const &regex_str);

std::string string_to_upper(std::string str);
std::string string_to_lower(std::string str);

void cstr_to_lower(char *cstr);
void cstr_to_upper(char *cstr);
char *double_to_att_str(int digits, char *str, size_t len, double value);

const char *tunit_to_cstr(int tunits);
const char *calendar_to_cstr(int calendar);

std::string get_scientific(double p_float_string);

bool string_is_float(std::string const &str);
bool string_is_int(std::string const &str);
void cstr_replace_char(char *str_in, char orig_char, char rep_char);

std::tuple<bool, std::vector<std::string>> tokenize_comma_seperated_int_list(std::string const &args);

template <typename... Args>
std::string
string_format(std::string const &format, Args... args)
{
  auto size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;  // Extra space for '\0'
  if (size_s <= 0) { throw std::runtime_error("Error during formatting."); }
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(), static_cast<char *>(buf.get()) + size - 1);  // We don't want the '\0' inside
}

inline bool
string_contains(std::string const &s, unsigned char ch)
{
  return (s.find(ch) != std::string::npos);
}

template <typename Precission>
Precission
string_to_number(std::string const &str)
{
  static_assert(std::is_arithmetic<Precission>::value, "Only arithmetic types allowed");
  if (str.empty()) { throw std::invalid_argument("Error, conversion of " + str + " not possible, string empty"); }

  std::stringstream ss(str);
  Precission number = 0.0;
  ss >> number;
  if (!ss.eof()) { throw std::invalid_argument("Error, conversion of " + str + " not possible"); }
  return number;
}

template <typename Precission>
std::tuple<bool, Precission>
string_to_integral(std::string const &str)
{
  static_assert(std::is_integral<Precission>::value, "Only integral number types allowed");
  Precission number;
  /* clang-format off */
  try
    { number = string_to_number<Precission>(str); }
  catch (std::invalid_argument &e)
    { return { false, std::numeric_limits<Precission>::max() }; }
  /* clang-format on */

  return { true, number };
}

template <typename Precission>
std::tuple<bool, Precission>
string_to_floating(std::string const &str)
{
  static_assert(std::is_floating_point<Precission>::value, "Only floating point number types allowed");
  Precission number;
  /* clang-format off */
  try
    { number = string_to_number<Precission>(str); }
  catch (std::invalid_argument &e)
    { return { false, std::numeric_limits<Precission>::quiet_NaN() }; }
  /* clang-format on */

  return { true, number };
}
namespace Util
{
namespace String
{
static inline std::string
ltrim(std::string s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  return s;
}

static inline std::string
rtrim(std::string s)
{
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
  return s;
}

static inline std::string
trim(std::string s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
  return s;
}
}  // namespace String
}  // namespace Util

#endif
