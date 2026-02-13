/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#include <cstring>

#include <cdi.h>

#include "cdi_uuid.h"
#include "cdi_lockedIO.h"
#include "grid_options.h"
#include "gridreference.h"
#include "cdo_output.h"
#include <mpim_grid.h>

// callback function for curl for writing the network retrieved grid file
#ifdef HAVE_LIBCURL
static size_t
write_data(void *ptr, size_t size, size_t nmemb, std::FILE *stream)
{
  return std::fwrite(ptr, size, nmemb, stream);
}
#endif

// code from grid_tools.2
static int
download_gridfile(const char *uri, const char *basename)
{
  int rval = 1;
#ifdef HAVE_LIBCURL
  // As curl_easy_init calls non-thread safe curl_global_init the libcurl
  // developer advice to call curl_global_init first and before potential thread spawning.
  int curlflags = CURL_GLOBAL_DEFAULT;

#ifdef CURL_GLOBAL_ACK_EINTR
  curlflags |= CURL_GLOBAL_ACK_EINTR;
#endif

  auto ret = curl_global_init(curlflags);
  if (ret != 0)
  {
    std::fprintf(stderr, "ERROR: %s!\n", curl_easy_strerror(ret));
    return -1;
  }

  auto curl = curl_easy_init();
  if (curl == nullptr)
  {
    std::fprintf(stderr, "ERROR: could not get curl handler.\n");
    return -1;
  }
  else
  {
    auto fp = std::fopen(basename, "w");
    if (fp == nullptr)
    {
      std::fprintf(stderr, "ERROR: could not open local output file %s. %s.\n", basename, strerror(errno));
      return -1;
    }

    // if (gridVerbose) curl_easy_setopt(curl, CURLOPT_VERBOSE, 1);
    // curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1);
    curl_easy_setopt(curl, CURLOPT_URL, uri);
    // curl_easy_setopt(curl, CURLOPT_BUFFERSIZE, 131072L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    ret = curl_easy_perform(curl);
    std::fclose(fp);
    if (ret == 0)
    {
      /*
      int ihead;
      curl_easy_getinfo(curl, CURLINFO_HEADER_SIZE, &ihead);
      printf("ihead %d\n", ihead);
      */
      char *ctype;
      curl_easy_getinfo(curl, CURLINFO_CONTENT_TYPE, &ctype);

      if (std::strstr(ctype, "html") == nullptr)  // no html content
      {
        if (gridVerbose)
        {
          curl_off_t total = 0, length = 0;
          curl_easy_getinfo(curl, CURLINFO_TOTAL_TIME_T, &total);
          curl_easy_getinfo(curl, CURLINFO_SIZE_DOWNLOAD_T, &length);
          cdo_print("File %s downloaded: %dMB/s", basename, (int) ((length / (1024. * 1024.)) / (total / 1000000.)));
        }
        rval = 0;
      }
      else
      {
        int status = remove(basename);
        if (status == -1) std::perror(basename);
        cdo_warning("The requested URL was not found (%s)!", uri);
      }
    }
    else
    {
      int status = remove(basename);
      if (status == -1) std::perror(basename);
      std::fprintf(stderr, "ERROR: %s. Download %s failed.\n\n", curl_easy_strerror(ret), basename);
    }

    curl_easy_cleanup(curl);
  }
#else
  (void) uri;
  (void) basename;

  cdo_warning("CURL support not compiled in!");
#endif

  return rval;
}

// Search for directory
static int
search_directory(const char *directory)
{
#ifdef HAVE_SYS_STAT_H
  struct stat buf;
  if (stat(directory, &buf) == 0) return 0;
#endif

  return 1;
}

// Search for filename
static int
search_file(const char *directory, const char *filename)
{
#ifdef HAVE_SYS_STAT_H
  struct stat buf;
  if (stat(directory, &buf) == 0)
  {
    if (stat(filename, &buf) == 0)
    {
      if (buf.st_size != 0 && !(buf.st_mode & S_IFDIR)) return 0;
    }
  }
  else
  {
    std::perror(directory);
  }
#endif

  return 1;
}

static int
grid_from_URI(char *griduri, char *gridfilepath)
{
  int status = -1;

  if (griduri[0])
  {
    cdo_print("Download horizontal grid file to %s", gridfilepath);
    if (gridVerbose) cdo_print("Download horizontal grid file from %s", griduri);
    status = download_gridfile(griduri, gridfilepath);
  }

  return status;
}

static int
grid_from_file(int gridID1, char *gridfilepath)
{
  int gridID2 = -1;

  if (gridVerbose) cdo_print("Horizontal grid file used: %s", gridfilepath);

  auto gridsize = gridInqSize(gridID1);

  int position = 0;
  cdiInqKeyInt(gridID1, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);

  open_lock();
  auto streamID = streamOpenRead(gridfilepath);
  if (streamID < 0) cdi_open_error(streamID, "Open failed on horizontal grid file >%s<", gridfilepath);
  open_unlock();

  auto vlistID = streamInqVlist(streamID);
  auto ngrids = vlistNumGrids(vlistID);
  if (position > 0 && position <= ngrids)
  {
    auto gridID = vlistGrid(vlistID, position - 1);
    if (gridInqSize(gridID) == gridsize)
      gridID2 = gridDuplicate(gridID);
    else
      cdo_warning("Grid size %zu on position %d do not match! Reference=%s", gridsize, position, gridfilepath);
  }
  else if (position == 0)
  {
    for (int grididx = 0; grididx < ngrids; ++grididx)
    {
      auto gridID = vlistGrid(vlistID, grididx);
      if (gridInqSize(gridID) == gridsize)
      {
        gridID2 = gridDuplicate(gridID);
        break;
      }
    }
  }
  else
    cdo_warning("Number of grid in reference %d not available! Reference=%s", position, gridfilepath);

  if (gridID2 != -1)
  {
    unsigned char uuidOfHGrid1[CDI_UUID_SIZE] = { 0 };
    unsigned char uuidOfHGrid2[CDI_UUID_SIZE] = { 0 };

    int length = CDI_UUID_SIZE;
    cdiInqKeyBytes(gridID1, CDI_GLOBAL, CDI_KEY_UUID, uuidOfHGrid1, &length);
    length = CDI_UUID_SIZE;
    cdiInqKeyBytes(gridID2, CDI_GLOBAL, CDI_KEY_UUID, uuidOfHGrid2, &length);

    if (!cdiUUIDIsNull(uuidOfHGrid1) && !cdiUUIDIsNull(uuidOfHGrid2) && memcmp(uuidOfHGrid1, uuidOfHGrid2, CDI_UUID_SIZE))
      cdo_warning("UUID of horizontal grids differ!");

    int number1 = 0, number2 = 0;
    cdiInqKeyInt(gridID1, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number1);
    cdiInqKeyInt(gridID2, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number2);
    if (number1 > 0 && number2 > 0 && number1 != number2) cdo_warning("Number of grid used of horizontal grids differ!");
  }

  streamClose(streamID);

  return gridID2;
}

static int
referenceToGrid(int gridID1)
{
  int gridID2 = -1;

  char griduri[8192] = { 0 };
  int length = 0;
  if (CDI_NOERR == cdiInqKeyLen(gridID1, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
    cdiInqKeyString(gridID1, CDI_GLOBAL, CDI_KEY_REFERENCEURI, griduri, &length);

  if (griduri[0] == 0) { cdo_warning("Reference to horizontal grid not available!"); }
  else
  {
    auto filename = strrchr(griduri, '/');
    if (filename == nullptr)
    {
      filename = griduri;
      griduri[0] = 0;
    }
    else
    {
      filename++;
    }

    int status = 0;
    char griddir[8192] = { 0 };
    char gridfilepath[8192] = { 0 };

    if (!IconGrids.empty() && griduri[0])
    {
      if (search_directory(IconGrids.c_str())) cdo_abort("CDO_ICON_GRIDS not found: %s!", IconGrids);

      auto wpath = griduri;
      if (strncmp(griduri, "http://", 7) == 0) { wpath += 7; }
      else if (strncmp(griduri, "https://", 8) == 0) { wpath += 8; }

      if (wpath != griduri)
      {
        auto const *gridpath = strchr(wpath, '/');
        if (gridpath)
        {
          strcpy(griddir, IconGrids.c_str());
          strcpy(gridfilepath, griddir);
          strcat(gridfilepath, gridpath);
          status = search_file(griddir, gridfilepath);
          if (status != 0) gridfilepath[0] = 0;
        }
      }
    }

    if (gridfilepath[0] == 0)
    {
      if (!DownloadPath.empty())
      {
        strcpy(griddir, DownloadPath.c_str());
        strcat(griddir, "/");

        status = search_directory(griddir);
        if (status != 0) cdo_abort("Download path not found: %s!", griddir);

        strcpy(gridfilepath, griddir);
        strcat(gridfilepath, filename);
      }
      else
      {
        strcpy(griddir, "./");

        strcpy(gridfilepath, griddir);
        strcat(gridfilepath, filename);

#ifdef HAVE_SYS_STAT_H
        status = search_file(griddir, gridfilepath);
        if (status != 0)
#endif
        {
          cdo_print("CDO_DOWNLOAD_PATH not set!");
          cdo_print("Set the environment variable CDO_DOWNLOAD_PATH to download gridfile %s.", filename);
          return gridID2;
        }
      }
    }

    if (gridVerbose) cdo_print("Search for horizontal grid file: %s", gridfilepath);

    // scan local directory for file
    status = search_file(griddir, gridfilepath);
    if (status != 0) status = grid_from_URI(griduri, gridfilepath);
    if (status == 0) gridID2 = grid_from_file(gridID1, gridfilepath);
  }

  return gridID2;
}

RefGrid
dereferenceGrid(int gridID)
{
  RefGrid reference;

  int number_of_grid_used = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number_of_grid_used);
  if (number_of_grid_used > 0)
  {
    reference.exists = true;
    reference.gridID = referenceToGrid(gridID);
    reference.isValid = (reference.gridID != -1);
    reference.notFound = (reference.gridID == -1);
  }

  return reference;
}
