#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"

void
recinfoInitEntry(recinfo_t *recinfo)
{
  recinfo->varID = CDI_UNDEFID;
  recinfo->levelID = CDI_UNDEFID;
  recinfo->used = false;
}

static void
recordInitEntry(record_t *record)
{
  record->position = CDI_UNDEFID;
  record->size = 0;
  record->gridsize = 0;
  record->param = 0;
  record->ilevel = CDI_UNDEFID;
  record->tsteptype = CDI_UNDEFID;
#ifdef HAVE_LIBGRIB
  varScanKeysInit(&record->scanKeys);
  memset(&record->tiles, 0, sizeof(record->tiles));
#ifdef HAVE_LIBGRIB_API
  memset(record->varname, 0, sizeof(record->varname));
#endif
#endif
#ifdef HAVE_LIBFDB5
  record->fdbItemIndex = -1;
#endif
}

int
recordNewEntry(stream_t *streamptr, int tsID)
{
  int recordSize = streamptr->tsteps[tsID].recordSize;
  recinfo_t *recinfo = streamptr->tsteps[tsID].recinfo;
  record_t *records = streamptr->tsteps[tsID].records;

  // Look for a free slot in record.
  int recordID = 0;
  if (recordSize)
  {
    while (recordID < recordSize && recinfo[recordID].used != CDI_UNDEFID) ++recordID;
  }
  else  // Create the table the first time through.
  {
    recordSize = 1;  //  <<<<----
    recinfo = (recinfo_t *) Malloc((size_t) recordSize * sizeof(recinfo_t));
    records = (record_t *) Malloc((size_t) recordSize * sizeof(record_t));
    for (int i = recordID; i < recordSize; i++) recinfo[i].used = CDI_UNDEFID;
  }

  // If the table overflows, double its size.
  if (recordID == recordSize)
  {
    // clang-format off
      if      (recordSize <= INT_MAX / 2) recordSize *= 2;
      else if (recordSize < INT_MAX)      recordSize = INT_MAX;
      else Error("Cannot handle this many records!\n");
    // clang-format on
    recinfo = (recinfo_t *) Realloc(recinfo, (size_t) recordSize * sizeof(recinfo_t));
    records = (record_t *) Realloc(records, (size_t) recordSize * sizeof(record_t));

    for (int i = recordID; i < recordSize; i++) recinfo[i].used = CDI_UNDEFID;
  }

  recinfoInitEntry(&recinfo[recordID]);
  recordInitEntry(&records[recordID]);

  recinfo[recordID].used = true;

  streamptr->tsteps[tsID].recordSize = recordSize;
  streamptr->tsteps[tsID].recinfo = recinfo;
  streamptr->tsteps[tsID].records = records;

  return recordID;
}

static void
cdiInitRecord(stream_t *streamptr)
{
  Record *record = (Record *) Malloc(sizeof(Record));
  streamptr->record = record;

  record->param = 0;
  record->ilevel = 0;
  record->vdate = 0;
  record->vtime = 0;
  record->gridID = 0;
  record->buffer = NULL;
  record->buffersize = 0;
  record->position = 0;
  record->varID = 0;
  record->levelID = CDI_UNDEFID;
}

void
stream_inq_field(stream_t *streamPtr, int *varID, int *levelID)
{
  check_parg(varID);
  check_parg(levelID);

  stream_def_accesstype(streamPtr, TYPE_REC);

  if (!streamPtr->record) cdiInitRecord(streamPtr);

  const int tsID = streamPtr->curTsID;
  tsteps_t *tstep = &(streamPtr->tsteps[tsID]);

  int rindex = tstep->curRecID + 1;
  if (rindex >= tstep->nrecs) Error("record %d not available at timestep %d", rindex + 1, tsID + 1);

  int recID = tstep->recIDs[rindex];
  if (recID == -1 || recID >= tstep->nallrecs) Error("Internal problem! tsID = %d recID = %d", tsID, recID);

  *varID = tstep->recinfo[recID].varID;
  if (*varID == -1) Error("Internal problem! varID = %d recID = %d", *varID, recID);

  int lindex = tstep->recinfo[recID].levelID;
  int isub = subtypeInqActiveIndex(streamPtr->vars[*varID].subtypeID);
  *levelID = streamPtr->vars[*varID].recordTable[isub].lindex[lindex];

  if (CDI_Debug)
    Message("streamID = %d tsID = %d, recID = %d, varID = %d, levelID = %d", streamPtr->self, tsID, recID, *varID, *levelID);

  streamPtr->curTsID = tsID;
  tstep->curRecID = rindex;
}

void
pstreamInqField(void *streamPtr, int *varID, int *levelID)
{
  stream_inq_field((stream_t *) streamPtr, varID, levelID);
}

void
streamInqField(int streamID, int *varID, int *levelID)
{
  stream_inq_field(stream_to_pointer(streamID), varID, levelID);
}

/*
@Function  streamDefRecord
@Title     Define the next record

@Prototype void streamDefRecord(int streamID, int varID, int levelID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.

@Description
The function streamDefRecord defines the meta-data of the next record.
@EndFunction
*/
void
streamDefField(int streamID, int varID, int levelID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  int tsID = streamptr->curTsID;
  if (tsID == CDI_UNDEFID)
  {
    tsID++;
    streamDefTimestep(streamID, tsID);
  }

  if (!streamptr->record) cdiInitRecord(streamptr);

  int vlistID = streamptr->vlistID;
  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int param = vlistInqVarParam(vlistID, varID);
  int ilevel = (int) lround(zaxisInqLevel(zaxisID, levelID));

  Record *record = streamptr->record;
  record->varID = varID;
  record->levelID = levelID;
  record->param = param;
  record->ilevel = ilevel;
  record->vdate = (int) cdiDate_get(streamptr->tsteps[tsID].taxis.vDateTime.date);
  record->vtime = cdiTime_get(streamptr->tsteps[tsID].taxis.vDateTime.time);
  record->gridID = gridID;
  record->prec = vlistInqVarDatatype(vlistID, varID);

  switch (cdiBaseFiletype(streamptr->filetype))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grbDefField(streamptr); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvDefField(streamptr); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extDefField(streamptr); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegDefField(streamptr); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF:
      if (streamptr->accessmode == 0) cdfEndDef(streamptr);
      cdfDefField(streamptr);
      break;
#endif
    default: Error("%s support not compiled in!", strfiletype(streamptr->filetype)); break;
  }
}

void
streamCopyField(int streamID2, int streamID1)
{
  stream_t *streamptr1 = stream_to_pointer(streamID1);
  stream_t *streamptr2 = stream_to_pointer(streamID2);
  int filetype1 = streamptr1->filetype;
  int filetype2 = streamptr2->filetype;
  int filetype = CDI_FILETYPE_UNDEF;

  if (cdiBaseFiletype(filetype1) == cdiBaseFiletype(filetype2)) filetype = filetype2;

  if (filetype == CDI_FILETYPE_UNDEF)
    Error("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));

  switch (cdiBaseFiletype(filetype))
  {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grbCopyField(streamptr2, streamptr1); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvCopyField(streamptr2, streamptr1); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extCopyField(streamptr2, streamptr1); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegCopyField(streamptr2, streamptr1); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdfCopyField(streamptr2, streamptr1); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype));
  }
}

void
cdi_create_records(stream_t *streamptr, int tsID, bool allocRecords)
{
  unsigned nrecords, maxrecords;

  tsteps_t *sourceTstep = streamptr->tsteps;
  tsteps_t *destTstep = sourceTstep + tsID;

  if (destTstep->recinfo) return;

  int vlistID = streamptr->vlistID;

  if (tsID == 0)
  {
    maxrecords = 0;
    int nvars = streamptr->nvars;
    for (int varID = 0; varID < nvars; varID++)
      for (int isub = 0; isub < streamptr->vars[varID].subtypeSize; isub++)
        maxrecords += (unsigned) streamptr->vars[varID].recordTable[isub].nlevs;
  }
  else { maxrecords = (unsigned) sourceTstep->recordSize; }

  if (tsID == 0) { nrecords = maxrecords; }
  else if (tsID == 1)
  {
    nrecords = 0;
    maxrecords = (unsigned) sourceTstep->recordSize;
    if (sourceTstep->records)
    {
      for (size_t recID = 0; recID < maxrecords; recID++)
      {
        int varID = sourceTstep->recinfo[recID].varID;
        nrecords += (varID == CDI_UNDEFID /* varID = CDI_UNDEFID for write mode !!! */
                     || vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT);
        //    printf("varID nrecords %d %d %d \n", varID, nrecords, vlistInqVarTsteptype(vlistID, varID));
      }
    }
    else { nrecords = maxrecords; }
  }
  else { nrecords = (unsigned) streamptr->tsteps[1].nallrecs; }
  //  printf("tsID, nrecords %d %d\n", tsID, nrecords);

  recinfo_t *recinfo = (maxrecords > 0) ? (recinfo_t *) Malloc(maxrecords * sizeof(recinfo_t)) : (recinfo_t *) NULL;
  record_t *records = (allocRecords && maxrecords > 0) ? (record_t *) Malloc(maxrecords * sizeof(record_t)) : (record_t *) NULL;

  destTstep->recinfo = recinfo;
  destTstep->records = records;
  destTstep->recordSize = (int) maxrecords;
  destTstep->nallrecs = (int) nrecords;
#ifdef HAVE_LIBFDB5
  if (destTstep->records) destTstep->records->fdbItemIndex = -1;
#endif

  if (tsID == 0)
  {
    for (unsigned recID = 0; recID < maxrecords; recID++) recinfoInitEntry(&destTstep->recinfo[recID]);
    if (allocRecords)
      for (unsigned recID = 0; recID < maxrecords; recID++) recordInitEntry(&destTstep->records[recID]);
  }
  else if (sourceTstep->recinfo)
  {
    memcpy(destTstep->recinfo, sourceTstep->recinfo, (size_t) maxrecords * sizeof(recinfo_t));
    if (allocRecords) memcpy(destTstep->records, sourceTstep->records, (size_t) maxrecords * sizeof(record_t));

    for (size_t recID = 0; recID < maxrecords; recID++)
    {
      recinfo_t *curRecord = &sourceTstep->recinfo[recID];
      destTstep->recinfo[recID].used = curRecord->used;
      if (curRecord->used != CDI_UNDEFID && curRecord->varID != -1)  // curRecord->varID = -1 for write mode !!!
      {
        if (vlistInqVarTimetype(vlistID, curRecord->varID) != TIME_CONSTANT)
        {
          if (allocRecords) destTstep->records[recID].position = CDI_UNDEFID;
          if (allocRecords) destTstep->records[recID].size = 0;
          destTstep->recinfo[recID].used = false;
        }
      }
    }
  }
}

#include "file.h"

void
streamFCopyRecord(stream_t *streamptr2, stream_t *streamptr1, const char *container_name)
{
  int fileID1 = streamptr1->fileID;
  int fileID2 = streamptr2->fileID;

  int tsID = streamptr1->curTsID;
  int vrecID = streamptr1->tsteps[tsID].curRecID;
  int recID = streamptr1->tsteps[tsID].recIDs[vrecID];
  off_t recpos = streamptr1->tsteps[tsID].records[recID].position;
  size_t recsize = streamptr1->tsteps[tsID].records[recID].size;

  if (fileSetPos(fileID1, recpos, SEEK_SET) != 0) Error("Cannot seek input file for %s record copy!", container_name);

  char *buffer = (char *) Malloc(recsize);

  if (fileRead(fileID1, buffer, recsize) != recsize) Error("Failed to read record from %s file for copying!", container_name);

  if (fileWrite(fileID2, buffer, recsize) != recsize) Error("Failed to write record to %s file when copying!", container_name);

  Free(buffer);
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
