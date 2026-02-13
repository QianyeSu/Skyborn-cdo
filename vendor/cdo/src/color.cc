/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "color.h"

#define RGB 0
#define HSV 1
#define CMYK 2

#define irint(x) ((int) rint(x))

int
check_rgb(int rgb[])
{
  return (((rgb[0] < 0 || rgb[0] > 255) || (rgb[1] < 0 || rgb[1] > 255) || (rgb[2] < 0 || rgb[2] > 255)));
}

int
check_hsv(double h, double s, double v)
{
  return (h < 0.0 || h > 360.0 || s < 0.0 || s > 1.0 || v > 1.0);
}

static bool
check_cmyk(double cmyk[])
{
  for (int i = 0; i < 4; ++i)
    if (cmyk[i] < 0.0 || cmyk[i] > 100.0) return true;
  return false;
}

void
hsv_to_rgb(int rgb[], double h, double s, double v)
{
  int i;
  double f, p, q, t, rr = 0, gg = 0, bb = 0;

  if (std::fabs(s) <= 0)
    rgb[0] = rgb[1] = rgb[2] = (int) std::floor(255.999 * v);
  else
  {
    while (h >= 360.0) h -= 360.0;
    h /= 60.0;
    i = (int) h;
    f = h - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * f));
    t = v * (1.0 - (s * (1.0 - f)));
    switch (i)
    {
      case 0:
        rr = v;
        gg = t;
        bb = p;
        break;
      case 1:
        rr = q;
        gg = v;
        bb = p;
        break;
      case 2:
        rr = p;
        gg = v;
        bb = t;
        break;
      case 3:
        rr = p;
        gg = q;
        bb = v;
        break;
      case 4:
        rr = t;
        gg = p;
        bb = v;
        break;
      case 5:
        rr = v;
        gg = p;
        bb = q;
        break;
    }

    rgb[0] = (rr < 0.0) ? 0 : (int) std::floor(rr * 255.999);
    rgb[1] = (gg < 0.0) ? 0 : (int) std::floor(gg * 255.999);
    rgb[2] = (bb < 0.0) ? 0 : (int) std::floor(bb * 255.999);
  }
}

void
cmyk_to_rgb(int rgb[], double cmyk[])
{
  // Plain conversion; no undercolor removal or blackgeneration

  // CMYK is in 0-100, RGB will be in 0-255 range
  for (int i = 0; i < 3; ++i) rgb[i] = (int) std::floor((100.0 - cmyk[i] - cmyk[3]) * 2.55999);
}

int
slash_count(char *txt)
{
  int i = 0, n = 0;
  while (txt[i])
    if (txt[i++] == '/') n++;
  return n;
}

static bool
getrgb(char *line, int rgb[], int color_model)
{
  int n;

  int count = slash_count(line);

  if (count == 3)
  { /* c/m/y/k */
    double cmyk[4];
    n = std::sscanf(line, "%lf/%lf/%lf/%lf", &cmyk[0], &cmyk[1], &cmyk[2], &cmyk[3]);
    if (n != 4 || check_cmyk(cmyk)) return true;
    cmyk_to_rgb(rgb, cmyk);
    return false;
  }

  if (count == 2)
  { /* r/g/b or h/s/v */
    if (color_model == RGB)
    { /* r/g/b */
      double r, g, b;
      n = std::sscanf(line, "%lf/%lf/%lf", &r, &g, &b);
      rgb[0] = r;
      rgb[1] = g;
      rgb[2] = b;
      if (n != 3 || check_rgb(rgb)) return true;
    }
    else
    { /* h/s/v */
      double h, s, v;
      n = std::sscanf(line, "%lf/%lf/%lf", &h, &s, &v);
      if (n != 3 || check_hsv(h, s, v)) return true;
      hsv_to_rgb(rgb, h, s, v);
    }
    return false;
  }

  if (count == 0)
  { /* gray */
    n = std::sscanf(line, "%d", &rgb[0]);
    rgb[1] = rgb[2] = rgb[0];
    if (n != 1 || check_rgb(rgb)) return true;
    return false;
  }

  /* Get here if there is a problem */

  return true;
}

#define READERR -1

int
cpt_read(std::FILE *fp, CPT *cpt)
{
  int ncolors;
  int status = 0;
  /* Opens and reads a color palette file in RGB, HSV, or CMYK of arbitrary length */
  int small_chunk = 64;
  int n = 0, i, nread, annot, n_alloc = small_chunk, color_model, id;
  double dz;
  bool gap;
  int error = 0;
  char T0[64], T1[64], T2[64], T3[64], T4[64], T5[64], T6[64], T7[64], T8[64], T9[64];
  char line[BUFSIZ], option[260], c;

  if (fp == nullptr) return (READERR);

  cpt->lut.resize(n_alloc);

  /* Save the original setting since it may be modified by settings in the CPT file */
  color_model = RGB;

  while (!error && fgets(line, BUFSIZ, fp))
  {
    if (std::strstr(line, "COLOR_MODEL"))
    { /* cpt file overrides default color model */
      if (std::strstr(line, "RGB") || std::strstr(line, "rgb"))
        color_model = RGB;
      else if (std::strstr(line, "HSV") || std::strstr(line, "hsv"))
        color_model = HSV;
      else if (std::strstr(line, "CMYK") || std::strstr(line, "cmyk"))
        color_model = CMYK;
      else
      {
        std::fprintf(stderr, "%s: unrecognized COLOR_MODEL\n", __func__);
        return (READERR);
      }
    }

    c = line[0];
    if (c == '#' || c == '\n') continue; /* Comment or blank */

    T1[0] = T2[0] = T3[0] = T4[0] = T5[0] = T6[0] = T7[0] = T8[0] = T9[0] = 0;
    switch (c)
    {
      case 'B': id = 0; break;
      case 'F': id = 1; break;
      case 'N': id = 2; break;
      default: id = 3; break;
    }

    if (id < 3)
    { /* Foreground, background, or nan color */
      cpt->bfn[id].skip = false;
      if ((nread = std::sscanf(&line[2], "%s %s %s %s", T1, T2, T3, T4)) < 1) error++;
      if (T1[0] == 'p' || T1[0] == 'P')
      { /* Gave a pattern */
        std::fprintf(stderr, "%s: CPT Pattern fill (%s) unsupported!\n", __func__, T1);
        return (READERR);
      }
      else
      {                   /* Shades, RGB, HSV, or CMYK */
        if (T1[0] == '-') /* Skip this slice */
          cpt->bfn[id].skip = true;
        else if (nread == 1)
        { /* Gray shade */
          std::snprintf(option, sizeof(option), "%s", T1);
          if (getrgb(option, cpt->bfn[id].rgb, color_model)) error++;
        }
        else if (color_model == CMYK)
        {
          std::snprintf(option, sizeof(option), "%s/%s/%s/%s", T1, T2, T3, T4);
          if (getrgb(option, cpt->bfn[id].rgb, color_model)) error++;
        }
        else
        {
          std::snprintf(option, sizeof(option), "%s/%s/%s", T1, T2, T3);
          if (getrgb(option, cpt->bfn[id].rgb, color_model)) error++;
        }
      }
      continue;
    }

    /* Here we have regular z-slices.  Allowable formats are
     *
     * z0 - z1 - [LUB]
     * z0 pattern z1 - [LUB]
     * z0 r0 z1 r1 [LUB]
     * z0 r0 g0 b0 z1 r1 g1 b1 [LUB]
     * z0 h0 s0 v0 z1 h1 s1 v1 [LUB]
     * z0 c0 m0 y0 k0 z1 c1 m1 y1 k1 [LUB]
     */

    /* Determine if psscale need to label these steps by examining for the
     * optional L|U|B character at the end */

    c = line[strlen(line) - 2];
    if (c == 'L')
      cpt->lut[n].annot = 1;
    else if (c == 'U')
      cpt->lut[n].annot = 2;
    else if (c == 'B')
      cpt->lut[n].annot = 3;
    // Chop off this information so it does not affect our column count below
    if (cpt->lut[n].annot) line[strlen(line) - 2] = '\0';

    // Chop off this information so it does not affect our column count below
    nread = std::sscanf(line, "%s %s %s %s %s %s %s %s %s %s", T0, T1, T2, T3, T4, T5, T6, T7, T8, T9);

    if (nread <= 0) continue;                          /* Probably a line with spaces - skip */
    if (color_model == CMYK && nread != 10) error = 1; /* CMYK should results in 10 fields */
    if (color_model != CMYK && !(nread == 4 || nread == 8))
      error = 1; /* HSV or RGB should result in 8 fields, gray, patterns, or skips in 4 */

    cpt->lut[n].z_low = atof(T0);
    cpt->lut[n].skip = false;
    if (T1[0] == '-')
    { /* Skip this slice */
      if (nread != 4)
      {
        std::fprintf(stderr, "%s: z-slice to skip not in [z0 - z1 -] format!\n", __func__);
        return (READERR);
      }
      cpt->lut[n].z_high = atof(T2);
      cpt->lut[n].skip = true;                                                        /* Don't paint this slice if possible*/
      for (i = 0; i < 3; ++i) cpt->lut[n].rgb_low[i] = cpt->lut[n].rgb_high[i] = 255; /* If you must, use page color */
    }
    else if (T1[0] == 'p' || T1[0] == 'P')
    { /* Gave pattern fill */
      std::fprintf(stderr, "%s: CPT Pattern fill (%s) unsupported!\n", __func__, T1);
      return (READERR);
    }
    else
    { /* Shades, RGB, HSV, or CMYK */
      if (nread == 4)
      { /* gray shades */
        cpt->lut[n].z_high = atof(T2);
        cpt->lut[n].rgb_low[0] = cpt->lut[n].rgb_low[1] = cpt->lut[n].rgb_low[2] = irint(atof(T1));
        cpt->lut[n].rgb_high[0] = cpt->lut[n].rgb_high[1] = cpt->lut[n].rgb_high[2] = irint(atof(T3));
        if (cpt->lut[n].rgb_low[0] < 0 || cpt->lut[n].rgb_high[0] < 0) error++;
      }
      else if (color_model == CMYK)
      {
        cpt->lut[n].z_high = atof(T5);
        std::snprintf(option, sizeof(option), "%s/%s/%s/%s", T1, T2, T3, T4);
        if (getrgb(option, cpt->lut[n].rgb_low, color_model)) error++;
        std::snprintf(option, sizeof(option), "%s/%s/%s/%s", T6, T7, T8, T9);
        if (getrgb(option, cpt->lut[n].rgb_high, color_model)) error++;
      }
      else
      { /* RGB or HSV */
        cpt->lut[n].z_high = atof(T4);
        std::snprintf(option, sizeof(option), "%s/%s/%s", T1, T2, T3);
        if (getrgb(option, cpt->lut[n].rgb_low, color_model)) error++;
        std::snprintf(option, sizeof(option), "%s/%s/%s", T5, T6, T7);
        if (getrgb(option, cpt->lut[n].rgb_high, color_model)) error++;
      }

      dz = cpt->lut[n].z_high - cpt->lut[n].z_low;
      if (std::fabs(dz) <= 0)
      {
        std::fprintf(stderr, "%s: Z-slice with dz = 0\n", __func__);
        return (READERR);
      }
      cpt->lut[n].i_dz = 1.0 / dz;

      for (i = 0; i < 3; ++i) cpt->lut[n].rgb_diff[i] = cpt->lut[n].rgb_high[i] - cpt->lut[n].rgb_low[i]; /* Used in get_rgb24 */
    }

    n++;
    if (n == n_alloc)
    {
      i = n_alloc;
      n_alloc += small_chunk;
      cpt->lut.resize(n_alloc);
    }
  }

  std::fclose(fp);

  if (error)
  {
    std::fprintf(stderr, "%s: Decoding error\n", __func__);
    return (READERR);
  }

  if (n == 0)
  {
    std::fprintf(stderr, "%s: CPT file has no z-slices!\n", __func__);
    return (READERR);
  }

  cpt->lut.resize(n);
  ncolors = n;
  for (i = annot = 0, gap = false; i < ncolors - 1; ++i)
  {
    if (std::fabs(cpt->lut[i].z_high - cpt->lut[i + 1].z_low) > 0) gap = true;
    annot += cpt->lut[i].annot;
  }

  annot += cpt->lut[i].annot;
  if (gap)
  {
    std::fprintf(stderr, "%s: Color palette table has gaps - aborts!\n", __func__);
    return (READERR);
  }

  if (!annot)
  { /* Must set default annotation flags */
    for (i = 0; i < ncolors; ++i) cpt->lut[i].annot = 1;
    cpt->lut[i - 1].annot = 3;
  }

  cpt->ncolors = ncolors;

  return (status);
}

int
cpt_write(std::FILE *fp, const CPT &cpt)
{
  char code[3] = { 'B', 'F', 'N' };
  int n, k;
  int status = 0;

  for (n = 0; n < cpt.ncolors; ++n)
  {
    std::fprintf(fp, "%g\t%d\t%d\t%d\t%g\t%d\t%d\t%d\n", cpt.lut[n].z_low, cpt.lut[n].rgb_low[0], cpt.lut[n].rgb_low[1],
            cpt.lut[n].rgb_low[2], cpt.lut[n].z_high, cpt.lut[n].rgb_high[0], cpt.lut[n].rgb_high[1], cpt.lut[n].rgb_high[2]);
  }

  for (k = 0; k < 3; ++k)
  {
    if (cpt.bfn[k].skip)
      std::fprintf(fp, "%c -\n", code[k]);
    else
      std::fprintf(fp, "%c\t%d\t%d\t%d\n", code[k], cpt.bfn[k].rgb[0], cpt.bfn[k].rgb[1], cpt.bfn[k].rgb[2]);
  }

  return status;
}

int
cpt_write_c(std::FILE *fp, const CPT &cpt, const char *name)
{
  char lut_name[4096];
  char cpt_name[4096];
  int n, k;
  int status = 0;

  std::strcpy(lut_name, name);
  std::strcat(lut_name, "_lut");
  std::strcpy(cpt_name, name);
  std::strcat(cpt_name, "_cpt");

  std::fprintf(fp, "\nstatic LUT %s[] = {\n", lut_name);
  for (n = 0; n < cpt.ncolors; ++n)
  {
    std::fprintf(fp, "  { %7g, %7g, %7g, {%3d, %3d, %3d}, {%3d, %3d, %3d}, {%3d, %3d, %3d}, %d, %d},\n", cpt.lut[n].z_low,
            cpt.lut[n].z_high, cpt.lut[n].i_dz, cpt.lut[n].rgb_low[0], cpt.lut[n].rgb_low[1], cpt.lut[n].rgb_low[2],
            cpt.lut[n].rgb_high[0], cpt.lut[n].rgb_high[1], cpt.lut[n].rgb_high[2], cpt.lut[n].rgb_diff[0], cpt.lut[n].rgb_diff[1],
            cpt.lut[n].rgb_diff[2], cpt.lut[n].annot, cpt.lut[n].skip);
  }
  std::fprintf(fp, "};\n");

  std::fprintf(fp, "\nstatic const CPT %s = {\n", cpt_name);
  std::fprintf(fp, "  %d,\n", cpt.ncolors);
  std::fprintf(fp, "  %s,\n", lut_name);
  std::fprintf(fp, "  {\n");
  for (k = 0; k < 3; ++k)
  {
    std::fprintf(fp, "    {{%3d, %3d, %3d}, %d},\n", cpt.bfn[k].rgb[0], cpt.bfn[k].rgb[1], cpt.bfn[k].rgb[2], cpt.bfn[k].skip);
  }
  std::fprintf(fp, "  }\n");
  std::fprintf(fp, "};\n");

  return status;
}
