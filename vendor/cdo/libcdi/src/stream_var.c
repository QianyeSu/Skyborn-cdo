/* DO NOT REMOVE the config.h include file under any circumstances,
 * it's very much needed on some platforms */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
/* DO NOT REMOVE the above config.h include file under any
 * circumstances as long as it's the autoconf configuration header
 * used to build this package. When it's missing on some platforms,
 * some poor person has to do long, tedious debugging sessions, where
 * struct offsets almost imperceptibly change from one file to the
 * next to find out what happened */

#include "dmemory.h"

#include "cdi.h"
#include "cdi_int.h"

static void
record_init(sleveltable_t *record)
{
  record->nlevs = 0;
  record->recordID = NULL;
  record->lindex = NULL;
}

static void
var_init(svarinfo_t *var)
{
#ifdef HAVE_LIBNETCDF
  var->cdfCache = NULL;
#endif
  var->recordTable = NULL;

  var->ncvarid = CDI_UNDEFID;
  var->defmiss = false;

  var->subtypeSize = 0;

  var->gridID = CDI_UNDEFID;
  var->zaxisID = CDI_UNDEFID;
  var->tsteptype = CDI_UNDEFID;
  var->subtypeID = CDI_UNDEFID;
}

static int
streamvar_new_entry(stream_t *streamptr)
{
  int varID = 0;
  int streamvarSize = streamptr->varsAllocated;
  svarinfo_t *streamvar = streamptr->vars;
  /*
    Look for a free slot in streamvar.
    (Create the table the first time through).
  */
  if (!streamvarSize)
  {
    streamvarSize = 2;
    streamvar = (svarinfo_t *) Malloc((size_t) streamvarSize * sizeof(svarinfo_t));
    if (streamvar == NULL)
    {
      Message("streamvarSize = %d", streamvarSize);
      SysError("Allocation of svarinfo_t failed");
    }

    for (int i = 0; i < streamvarSize; i++) streamvar[i].isUsed = false;
  }
  else
  {
    while (varID < streamvarSize)
    {
      if (!streamvar[varID].isUsed) break;
      varID++;
    }
  }
  /*
    If the table overflows, double its size.
  */
  if (varID == streamvarSize)
  {
    streamvarSize = 2 * streamvarSize;
    streamvar = (svarinfo_t *) Realloc(streamvar, (size_t) streamvarSize * sizeof(svarinfo_t));
    if (streamvar == NULL)
    {
      Message("streamvarSize = %d", streamvarSize);
      SysError("Reallocation of svarinfo_t failed");
    }
    varID = streamvarSize / 2;

    for (int i = varID; i < streamvarSize; i++) streamvar[i].isUsed = false;
  }

  streamptr->varsAllocated = streamvarSize;
  streamptr->vars = streamvar;

  var_init(&streamptr->vars[varID]);

  streamptr->vars[varID].isUsed = true;

  return varID;
}

static void
allocate_record_table_entry(stream_t *streamptr, int varID, int subID, int nlevs)
{
  int *level = (int *) Malloc((size_t) nlevs * sizeof(int));
  int *lindex = (int *) Malloc((size_t) nlevs * sizeof(int));

  for (int levID = 0; levID < nlevs; levID++)
  {
    level[levID] = CDI_UNDEFID;
    lindex[levID] = levID;
  }

  sleveltable_t *record = &streamptr->vars[varID].recordTable[subID];
  record->nlevs = nlevs;
  record->recordID = level;
  record->lindex = lindex;
}

int
stream_new_var(stream_t *streamptr, int gridID, int zaxisID, int tilesetID)
{
  if (CDI_Debug) Message("gridID = %d  zaxisID = %d", gridID, zaxisID);

  int varID = streamvar_new_entry(streamptr);
  int nlevs = zaxisInqSize(zaxisID);

  streamptr->nvars++;

  svarinfo_t *var = &streamptr->vars[varID];
  var->gridID = gridID;
  var->zaxisID = zaxisID;

  int nsub = 1;
  if (tilesetID != CDI_UNDEFID) nsub = subtypeInqSize(tilesetID); /* e.g. no of tiles */
  if (CDI_Debug) Message("varID %d: create %d tiles with %d level(s), zaxisID=%d", varID, nsub, nlevs, zaxisID);
  var->recordTable = (sleveltable_t *) Malloc((size_t) nsub * sizeof(sleveltable_t));
  if (var->recordTable == NULL) SysError("Allocation of leveltable failed!");
  var->subtypeSize = nsub;

  for (int isub = 0; isub < nsub; isub++)
  {
    record_init(&streamptr->vars[varID].recordTable[isub]);
    allocate_record_table_entry(streamptr, varID, isub, nlevs);
    if (CDI_Debug) Message("streamptr->vars[varID].recordTable[isub].recordID[0]=%d", var->recordTable[isub].recordID[0]);
  }

  var->subtypeID = tilesetID;

  return varID;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
