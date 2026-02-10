#include "mpmo.h"

namespace MpMO
{

bool silentMode = false;
bool warningsEnabled = true;
bool verbose = false;
bool pedantic = false;
bool exitOnError = true;
unsigned DebugLevel = 0;
int padding_width = 40;
int context_padding = 18;  // determined by currently longest operator name

void
Debug_(const char *p_file, const char *p_func, int p_line, const char *context, int p_debugScope, std::function<void()> p_function)
{
  (void) p_func;
  (void) p_file;
  (void) p_line;
  (void) context;
  if (p_debugScope) p_function();
}

/** Function for padding debug information
 * If the string is larger than the padding the given variable for the padding
 * is permanently increased!
 */
std::string
get_padding(std::string const &debug_scope_string, int &p_padding_width)
{
  int len = debug_scope_string.size();
  while (p_padding_width - len <= 0) { p_padding_width += 5; }
  return std::string(p_padding_width - len, ' ');
}

void
Debug_(const char *p_file, const char *p_func, int p_line, const char *p_context, std::function<void()> p_function)
{

  (void) p_context;
  (void) p_func;
  (void) p_file;
  (void) p_line;
  if (DebugLevel > 0) p_function();
}

std::string
debug_scope_string(const char *p_file, const char *p_func, int p_line, const char *context)
{
  auto file = std::string(p_file);
  file = file.substr(file.find_last_of("/\\") + 1);
  std::string context_string = Cyan(context);
  context_string = context_string + get_padding(context_string, context_padding);
  std::string scope_string = context_string + std::string(p_func) + ": " + std::string(file) + ":" + std::to_string(p_line);

  scope_string = (scope_string + get_padding(scope_string, padding_width));

  return scope_string;
}

void
Verbose_(bool p_verbose, std::function<void()> p_function) noexcept
{
  if (p_verbose) p_function();
}

void
enable_silent_mode(bool enable)
{
  silentMode = enable;
}

void
enable_warnings(bool enable)
{
  warningsEnabled = enable;
}

void
enable_pedantic(bool enable)
{
  pedantic = enable;
}

void
enable_verbose(bool enable)
{
  verbose = enable;
}
}  // namespace MpMO
