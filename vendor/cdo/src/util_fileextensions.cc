
#include <cdi.h>
#include <cstring>

/*
 * Return the filetype extension (const char) for a given filetype (int)
 * TODO: handle lists of extensions i.e. grb and grb2 for GRIB2-format
 */
const char *
filetypeext(int filetype)
{
  switch (filetype)
  {
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2: return ".grb";
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NCZARR:
    case CDI_FILETYPE_NC5:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C: return ".nc";
    case CDI_FILETYPE_SRV: return ".srv";
    case CDI_FILETYPE_EXT: return ".ext";
    case CDI_FILETYPE_IEG: return ".ieg";
    default: return "";
  }
}

/*
 * Remove file extension:
 * -------------------------------------------------
 * Remove file extension if it is the expected one
 * Do nothing otherwise
 */
static void
rm_filetypeext(char *file, const char *ext)
{
  // length of filename
  auto namelen = std::strlen(file);
  // length of the original file extension
  auto extlen = std::strlen(ext);

  // delete original extension if it is the expected one
  if (std::strcmp(&file[namelen - extlen], ext) == 0) file[namelen - extlen] = 0;
}

/*
 * Replace or just add file extension:
 * -------------------------------------------------
 * Replace file extension with new one or just add the new file extension
 * if the original extension is not the expected one
 */
void
repl_filetypeext(char *file, const char *oldext, const char *newext)
{
  // delete original extension if it is the expected one
  rm_filetypeext(file, oldext);

  // add new file extension
  std::strcat(file, newext);
}
