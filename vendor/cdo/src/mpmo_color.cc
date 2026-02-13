/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <functional>
#include "cdo_options.h"
#include "mpmo_color.h"

int colorBehaviour = Auto;

bool
color_enabled()
{
  if (colorBehaviour == Auto)
  {
    if (cdo::stderrIsTerminal && cdo::stdoutIsTerminal) return true;
  }
  if (colorBehaviour == All) return true;
  if (colorBehaviour == No) return false;
  return false;
}

void
mpmo_color_set(int p_colorBehaviour)
{
  colorBehaviour = p_colorBehaviour;
}

int
mpmo_get_color_mode()
{
  return colorBehaviour;
}

int CDO_Color = 1;

void
set_text_color(std::FILE *fp, TextMode attr, TextColor fg)
{
  if (!color_enabled()) return;
  if (!color_enabled()) return;

  // std::fprintf(fp, "%c[%d", 0x1B, attr);
  // if (fg >= 0)
  //  {
  //    std::fprintf(fp, ";%d", fg + 30);
  //    if (bg >= 0) std::fprintf(fp, ";%d", bg + 40);
  //  }
  //
  std::string color_str = text_color(fg, attr, NO_COLOR);
  std::fprintf(fp, "%s", color_str.c_str());
}

void
set_text_color(std::FILE *fp, TextMode attr)
{
  set_text_color(fp, attr, NO_COLOR);
}

void
set_text_color(std::FILE *fp, TextColor fg)
{
  set_text_color(fp, MODELESS, fg);
}

void
reset_text_color(std::FILE *fp)
{
  if (fp == stdout && !color_enabled()) return;
  if (fp == stderr && !color_enabled()) return;

  std::fprintf(fp, "%c[%dm", 0x1B, MODELESS);
}

void
colorize(std::FILE *fp, TextColor fg, TextMode mode, std::function<void()> func)
{
  set_text_color(fp, mode, fg);
  func();
  reset_text_color(fp);
}

void
colorize(std::FILE *fp, TextColor fg, std::function<void()> func)
{
  colorize(fp, fg, MODELESS, func);
}
