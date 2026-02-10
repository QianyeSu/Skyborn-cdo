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
#include "cdi_int.h"

static void
cdf_init_timestep(tsteps_t *timeStep, int numRecs, int numRecsAvail)
{
  timeStep->recinfo = (recinfo_t *) Malloc((size_t) numRecs * sizeof(recinfo_t));
  timeStep->nrecs = numRecsAvail;
  timeStep->nallrecs = numRecs;
  timeStep->recordSize = numRecs;
  timeStep->curRecID = CDI_UNDEFID;
}

static int
cdf_get_numRecsAvail(int vlistID)
{
  int numRecsAvail = 0;
  int numVars = vlistNvars(vlistID);
  for (int varID = 0; varID < numVars; varID++)
  {
    if (vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT) { numRecsAvail += zaxisInqSize(vlistInqVarZaxis(vlistID, varID)); }
  }
  return numRecsAvail;
}

static void
cdf_init_records_step0(int numRecs, int *recIDs, recinfo_t *recinfo, int vlistID)
{
  for (int recID = 0; recID < numRecs; recID++) recIDs[recID] = recID;

  int numVars = vlistNvars(vlistID);
  for (int varID = 0, recID = 0; varID < numVars; varID++)
  {
    int zaxisID = vlistInqVarZaxis(vlistID, varID);
    int nlevels = zaxisInqSize(zaxisID);
    for (int levelID = 0; levelID < nlevels; levelID++)
    {
      recinfoInitEntry(&recinfo[recID]);
      recinfo[recID].varID = (short) varID;
      recinfo[recID].levelID = levelID;
      recID++;
    }
  }
}

static void
cdf_init_records_step1(int numRecs, int *recIDs, recinfo_t *recinfo, int vlistID)
{
  for (int recID = 0, vrecID = 0; recID < numRecs; recID++)
  {
    if (vlistInqVarTimetype(vlistID, recinfo[recID].varID) != TIME_CONSTANT) { recIDs[vrecID++] = recID; }
  }
}

void
cdf_create_records(stream_t *streamptr, size_t tsID)
{
  if ((streamptr->ntsteps < 0 || tsID >= (size_t) streamptr->ntsteps) && tsID > 0) return;

  if (streamptr->tsteps[tsID].nallrecs > 0) return;

  int vlistID = streamptr->vlistID;

  tsteps_t *sourceTstep = streamptr->tsteps;
  tsteps_t *destTstep = sourceTstep + tsID;

  int numFields = vlistNumFields(vlistID);
  if (numFields <= 0) return;

  if (tsID == 0)
  {
    int numRecsAvail = numFields;  // use all records at first timestep

    streamptr->nrecs += numFields;

    cdf_init_timestep(destTstep, numFields, numRecsAvail);

    destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
    cdf_init_records_step0(numFields, destTstep->recIDs, destTstep->recinfo, vlistID);
  }
  else if (tsID == 1)
  {
    int numRecsAvail = cdf_get_numRecsAvail(vlistID);

    streamptr->nrecs += numRecsAvail;

    cdf_init_timestep(destTstep, numFields, numRecsAvail);

    memcpy(destTstep->recinfo, sourceTstep->recinfo, (size_t) numFields * sizeof(recinfo_t));

    if (numRecsAvail)
    {
      destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
      cdf_init_records_step1(numFields, destTstep->recIDs, destTstep->recinfo, vlistID);
    }
  }
  else
  {
    if (streamptr->tsteps[1].recinfo == 0) cdf_create_records(streamptr, 1);

    int numRecsAvail = streamptr->tsteps[1].nrecs;

    streamptr->nrecs += numRecsAvail;

    cdf_init_timestep(destTstep, numFields, numRecsAvail);

    memcpy(destTstep->recinfo, sourceTstep->recinfo, (size_t) numFields * sizeof(recinfo_t));

    if (numRecsAvail)
    {
      destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
      memcpy(destTstep->recIDs, streamptr->tsteps[1].recIDs, (size_t) numRecsAvail * sizeof(int));
    }
  }
}
