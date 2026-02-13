#include <cdi.h>

#include <cstdio>
#include <cstdint>
#include <filesystem>

#include "util_files.h"
#include "cdo_options.h"
#include "cdo_vlist.h"

#include "cdo_default_values.h"

bool
FileUtils::file_exists(std::string const &fileName)
{
  namespace fs = std::filesystem;
  /*
   auto isZarr = (fileName.starts_with("file://") && fileName.find("zarr", 6) != std::string::npos);
   if (isZarr)
     {
       cdo_abort("Enlargement of %s not possible!", fileName);
       int start = 7;
       auto pos = fileName.find("#mode", start);
       if (pos == std::string::npos) return false;
       auto zarrName = fileName.substr(start, pos - start);
       struct stat buf;
       auto status = stat(zarrName.c_str(), &buf);
       return (status == 0) && (S_ISDIR(buf.st_mode) && buf.st_size > 0);
     }
   else
   */
  {
    return (fs::exists(fileName) && fs::is_regular_file(fileName) && fs::file_size(fileName) > 0);
  }
}

bool
FileUtils::user_file_overwrite(std::string const &fileName)
{
  auto status = false;

  if (!Options::silentMode && cdo::stdinIsTerminal && cdo::stderrIsTerminal)
  {
    std::fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", fileName.c_str());
    std::string line;
    std::getline(std::cin, line);
    while (std::isspace((int) line[0])) line.erase(0, 1);
    if (line.size() == 3)
    {
      if (line.starts_with("yes") || line.starts_with("YES")) status = true;
    }
    else if (line.size() == 1)
    {
      if (line[0] == 'y' || line[0] == 'Y') status = true;
    }
  }

  return status;
}

std::uintmax_t
FileUtils::size(std::string const &fileName)
{
  namespace fs = std::filesystem;
  if (fileName[0] != '(' /* && filename[1] != 'p' */ && fs::exists(fileName)) { return fs::file_size(fileName); }
  return 0;
}

static std::string
gen_filesuffix(int filetype, std::string const &referenceName, int vlistID)
{
  std::string suffix;
  auto foundSuffix = false;
  auto isCompSZ = false;

  if (filetype == CdoDefault::FileType && CdoDefault::DataType == -1 && CdoDefault::Byteorder == -1)
  {
    size_t len = 0;
    if (referenceName.size() > 0 && referenceName[0] != '-' && referenceName[0] != '.') len = referenceName.size();

    if (len > 2)
    {
      auto pos = referenceName.find_last_of('.');
      if (pos > 1 && pos < (referenceName.size() - 1))
      {
        auto result = referenceName.substr(pos + 1);
        auto firstchar = std::tolower(result[1]);
        switch (firstchar)
        {
          case 'g':
            if (CdoDefault::FileType == CDI_FILETYPE_GRB || CdoDefault::FileType == CDI_FILETYPE_GRB2) foundSuffix = true;
            break;
          case 'n':
            if (CdoDefault::FileType == CDI_FILETYPE_NC || CdoDefault::FileType == CDI_FILETYPE_NC2
                || CdoDefault::FileType == CDI_FILETYPE_NC4 || CdoDefault::FileType == CDI_FILETYPE_NC4C
                || CdoDefault::FileType == CDI_FILETYPE_NC5)
              foundSuffix = true;
            break;
          case 's':
            if (CdoDefault::FileType == CDI_FILETYPE_SRV) foundSuffix = true;
            break;
          case 'e':
            if (CdoDefault::FileType == CDI_FILETYPE_EXT) foundSuffix = true;
            break;
          case 'i':
            if (CdoDefault::FileType == CDI_FILETYPE_IEG) foundSuffix = true;
            break;
        }

        if (foundSuffix)
        {
          for (int i = 0, n = result.size(); i < n; ++i)
          {
            if (result[i] == '.' || std::isalnum(result[i])) suffix += result[i];
          }
        }
      }
    }
  }

  if (!foundSuffix)
  {
    suffix += streamFilesuffix(CdoDefault::FileType);
    if (CdoDefault::FileType == CDI_FILETYPE_GRB && vlist_is_szipped(vlistID)) isCompSZ = true;
  }

  if (CdoDefault::FileType == CDI_FILETYPE_GRB && Options::cdoCompType == CDI_COMPRESS_SZIP) isCompSZ = true;
  if (isCompSZ) suffix += ".sz";

  return suffix;
}

std::string
FileUtils::gen_suffix(int filetype, int vlistID, std::string const &referenceName)
{
  std::string suffix;
  if (cdo::FileSuffix != "NULL")
  {
    if (cdo::FileSuffix.size()) { suffix = cdo::FileSuffix; }
    else
    {
      suffix = gen_filesuffix(filetype, referenceName, vlistID);
    }
  }
  return suffix;
}
