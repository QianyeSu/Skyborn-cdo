#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "cdi_key.h"
#include "dmemory.h"
#include "resource_handle.h"
#include "varscan.h"
#include "vlist.h"
#include "zaxis.h"
#include "subtype.h"

static size_t Vctsize = 0;
static double *Vct = NULL;

static int numberOfVerticalLevels = 0;
static int numberOfVerticalGrid = 0;
static unsigned char uuidVGrid[CDI_UUID_SIZE];

typedef struct
{
  int level1;
  int level2;
  int recID;
  int lindex;
} leveltable_t;

typedef struct
{
  int subtypeIndex;  //  corresponding tile in subtype_t structure (subtype->self)
  int nlevels;
  int levelTableSize;
  leveltable_t *levelTable;
} subtypetable_t;

typedef struct
{
  int varID;
  int param;
  int prec;
  int tsteptype;
  VarScanKeys scanKeys;
  int gridID;
  int zaxistype;
  int ltype1;  // GRIB first level type
  int ltype2;  // GRIB second level type
  int hasBounds;
  int level_sf;
  int level_unit;
  int zaxisID;

  int nsubtypes_alloc;
  int nsubtypes;
  subtypetable_t *recordTable;  // ~ two-dimensional record list [nsubtypes_alloc][levelTableSize]

  int instID;
  int modelID;
  int tableID;
  int comptype;   // compression type
  int complevel;  // compression level
  bool lmissval;
  double missval;
  char *name;

  // meta-data for specification of tiles (currently only GRIB-API:
  subtype_t *tiles;

  cdi_keys_t keys;

  int opt_grib_nentries;                // current no. key-value pairs
  int opt_grib_kvpair_size;             // current allocated size
  opt_key_val_pair_t *opt_grib_kvpair;  // (optional) list of keyword/value pairs
} VarInfo;

static VarInfo *varInfoList = NULL;
static int varInfoListSize = 0;
static int varInfoListUsed = 0;

static void
paramInitEntry(int varID, int param)
{
  VarInfo *var = &varInfoList[varID];
  var->varID = varID;
  var->param = param;
  var->prec = 0;
  var->tsteptype = TSTEP_INSTANT;
  varScanKeysInit(&(var->scanKeys));
  var->gridID = CDI_UNDEFID;
  var->zaxistype = 0;
  var->ltype1 = 0;
  var->ltype2 = -1;
  var->hasBounds = 0;
  var->level_sf = 0;
  var->level_unit = 0;
  var->recordTable = NULL;
  var->nsubtypes_alloc = 0;
  var->nsubtypes = 0;
  var->instID = CDI_UNDEFID;
  var->modelID = CDI_UNDEFID;
  var->tableID = CDI_UNDEFID;
  cdiInitKeys(&(var->keys));
  var->comptype = CDI_COMPRESS_NONE;
  var->complevel = 1;
  var->lmissval = false;
  var->missval = 0;
  var->name = NULL;
  var->tiles = NULL;
}

// Test if a variable specified by the given meta-data has already been registered in "vartable".
static int
varGetEntry(int param, int gridID, int zaxistype, int ltype1, int tsteptype, const char *name, const VarScanKeys *scanKeys,
            const var_tile_t *tiles)
{
  for (int varID = 0; varID < varInfoListSize; ++varID)
  {
    VarInfo *var = &varInfoList[varID];
    // testing for "param" implicitly checks if we are beyond the current vartable size:
    if (var->param == param)
    {
      int no_of_tiles = tiles ? tiles->numberOfTiles : -1;
      int vt_no_of_tiles = var->tiles ? subtypeGetGlobalDataP(var->tiles, SUBTYPE_ATT_NUMBER_OF_TILES) : -1;
      if ((var->zaxistype == zaxistype) && (var->ltype1 == ltype1) && (var->tsteptype == tsteptype)
          && (scanKeys == NULL || varScanKeysIsEqual(&var->scanKeys, scanKeys)) && (var->gridID == gridID)
          && (vt_no_of_tiles == no_of_tiles))
      {
        if (name && name[0] && var->name && var->name[0])
        {
          if (str_is_equal(name, var->name)) return varID;
        }
        else { return varID; }
      }
    }
  }

  return -1;
}

static void
varFree(void)
{
  if (CDI_Debug) Message("call to varFree");

  for (int varID = 0; varID < varInfoListUsed; ++varID)
  {
    VarInfo *var = &varInfoList[varID];
    if (var->recordTable)
    {
      for (int isub = 0; isub < var->nsubtypes_alloc; isub++) Free(var->recordTable[isub].levelTable);
      Free(var->recordTable);
    }

    if (var->name) Free(var->name);
    if (var->tiles) subtypeDestroyPtr(var->tiles);

    cdi_keys_t *keysp = &(var->keys);
    cdiDeleteVarKeys(keysp);

    if (var->opt_grib_kvpair)
    {
      for (int i = 0; i < var->opt_grib_nentries; i++)
      {
        if (var->opt_grib_kvpair[i].keyword) Free(var->opt_grib_kvpair[i].keyword);
      }
      Free(var->opt_grib_kvpair);
    }
    var->opt_grib_nentries = 0;
    var->opt_grib_kvpair_size = 0;
    var->opt_grib_kvpair = NULL;
  }

  if (varInfoList) Free(varInfoList);
  varInfoList = NULL;
  varInfoListSize = 0;
  varInfoListUsed = 0;

  if (Vct) Free(Vct);
  Vct = NULL;
  Vctsize = 0;
}

// Search for a tile subtype with subtypeIndex == tile_index.
static int
tileGetEntry(int varID, int tile_index)
{
  VarInfo *var = &varInfoList[varID];
  for (int isub = 0; isub < var->nsubtypes; isub++)
    if (var->recordTable[isub].subtypeIndex == tile_index) return isub;
  return CDI_UNDEFID;
}

/* Resizes vartable:recordTable data structure, if necessary. */
static int
tileNewEntry(int varID)
{
  VarInfo *var = &varInfoList[varID];
  int tileID = 0;
  if (var->nsubtypes_alloc == 0)
  {
    /* create table for the first time. */
    var->nsubtypes_alloc = 2;
    var->nsubtypes = 0;
    var->recordTable = (subtypetable_t *) Malloc((size_t) var->nsubtypes_alloc * sizeof(subtypetable_t));
    if (var->recordTable == NULL) SysError("Allocation of leveltable failed!");

    for (int isub = 0; isub < var->nsubtypes_alloc; isub++)
    {
      var->recordTable[isub].levelTable = NULL;
      var->recordTable[isub].levelTableSize = 0;
      var->recordTable[isub].nlevels = 0;
      var->recordTable[isub].subtypeIndex = CDI_UNDEFID;
    }
  }
  else
  {
    /* data structure large enough; find a free entry. */
    while (tileID < var->nsubtypes_alloc)
    {
      if (var->recordTable[tileID].levelTable == NULL) break;
      tileID++;
    }
  }

  /* If the table overflows, double its size. */
  if (tileID == var->nsubtypes_alloc)
  {
    tileID = var->nsubtypes_alloc;
    var->nsubtypes_alloc *= 2;
    var->recordTable = (subtypetable_t *) Realloc(var->recordTable, (size_t) var->nsubtypes_alloc * sizeof(subtypetable_t));
    if (var->recordTable == NULL) SysError("Reallocation of leveltable failed");
    for (int isub = tileID; isub < var->nsubtypes_alloc; isub++)
    {
      var->recordTable[isub].levelTable = NULL;
      var->recordTable[isub].levelTableSize = 0;
      var->recordTable[isub].nlevels = 0;
      var->recordTable[isub].subtypeIndex = CDI_UNDEFID;
    }
  }

  return tileID;
}

static int
levelNewEntry(int varID, int level1, int level2, int tileID)
{
  VarInfo *var = &varInfoList[varID];
  int levelID = 0;
  int levelTableSize = var->recordTable[tileID].levelTableSize;
  leveltable_t *levelTable = var->recordTable[tileID].levelTable;

  // Look for a free slot in levelTable. (Create the table the first time through).
  if (!levelTableSize)
  {
    levelTableSize = 2;
    levelTable = (leveltable_t *) Malloc((size_t) levelTableSize * sizeof(leveltable_t));
    for (int i = 0; i < levelTableSize; i++) levelTable[i].recID = CDI_UNDEFID;
  }
  else
  {
    while (levelID < levelTableSize && levelTable[levelID].recID != CDI_UNDEFID) ++levelID;
  }

  // If the table overflows, double its size.
  if (levelID == levelTableSize)
  {
    levelTable = (leveltable_t *) Realloc(levelTable, (size_t) (levelTableSize *= 2) * sizeof(leveltable_t));
    for (int i = levelID; i < levelTableSize; i++) levelTable[i].recID = CDI_UNDEFID;
  }

  levelTable[levelID].level1 = level1;
  levelTable[levelID].level2 = level2;
  levelTable[levelID].lindex = levelID;

  var->recordTable[tileID].nlevels = levelID + 1;
  var->recordTable[tileID].levelTableSize = levelTableSize;
  var->recordTable[tileID].levelTable = levelTable;

  return levelID;
}

#define UNDEF_PARAM -4711

static int
paramNewEntry(int param)
{
  int varID = 0;

  // Look for a free slot in vartable. (Create the table the first time through).
  if (!varInfoListSize)
  {
    varInfoListSize = 2;
    varInfoList = (VarInfo *) Malloc((size_t) varInfoListSize * sizeof(VarInfo));
    if (varInfoList == NULL)
    {
      Message("varTableSize = %d", varInfoListSize);
      SysError("Allocation of vartable failed");
    }

    for (int i = 0; i < varInfoListSize; i++)
    {
      varInfoList[i].param = UNDEF_PARAM;
      varInfoList[i].opt_grib_kvpair = NULL;
      varInfoList[i].opt_grib_kvpair_size = 0;
      varInfoList[i].opt_grib_nentries = 0;
    }
  }
  else
  {
    while (varID < varInfoListSize)
    {
      if (varInfoList[varID].param == UNDEF_PARAM) break;
      varID++;
    }
  }

  // If the table overflows, double its size.
  if (varID == varInfoListSize)
  {
    varInfoList = (VarInfo *) Realloc(varInfoList, (size_t) (varInfoListSize *= 2) * sizeof(VarInfo));
    for (int i = varID; i < varInfoListSize; i++)
    {
      varInfoList[i].param = UNDEF_PARAM;
      varInfoList[i].opt_grib_kvpair = NULL;
      varInfoList[i].opt_grib_kvpair_size = 0;
      varInfoList[i].opt_grib_nentries = 0;
    }
  }

  paramInitEntry(varID, param);

  return varID;
}

// Append tile set to a subtype. Return index of the new tile (i.e. the "entry->self" value).
static int
varInsertTileSubtype(VarInfo *vptr, const var_tile_t *tiles)
{
  if (tiles == NULL) return 0;

  // first, generate a subtype based on the info in "tiles".
  subtype_t *subtype_ptr;
  subtypeAllocate(&subtype_ptr, SUBTYPE_TILES);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS, tiles->totalno_of_tileattr_pairs);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TILE_CLASSIFICATION, tiles->tileClassification);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_NUMBER_OF_TILES, tiles->numberOfTiles);

  // Here, we create a tile set for comparison that contains only one tile/attribute pair (based on "tiles").
  struct subtype_entry_t *entry = subtypeEntryInsert(subtype_ptr);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_NUMBER_OF_ATTR, tiles->numberOfAttributes);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEINDEX, tiles->tileindex);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEATTRIBUTE, tiles->attribute);

  if (vptr->tiles == NULL)
  {
    vptr->tiles = subtype_ptr;
    return 0;
  }
  else
  {
    tilesetInsertP(vptr->tiles, subtype_ptr);
    subtypeDestroyPtr(subtype_ptr);
    return vptr->tiles->nentries - 1;
  }
}

void
varAddRecord(int recID, int param, int gridID, int zaxistype, int hasBounds, int level1, int level2, int level_sf, int level_unit,
             int prec, int *pvarID, int *plevelID, int tsteptype, int ltype1, int ltype2, const char *name,
             const VarScanKeys *scanKeys, const var_tile_t *tiles, int *tile_index)
{
  int varID = (CDI_Split_Ltype105 != 1 || zaxistype != ZAXIS_HEIGHT)
                  ? varGetEntry(param, gridID, zaxistype, ltype1, tsteptype, name, scanKeys, tiles)
                  : CDI_UNDEFID;

  if (varID == CDI_UNDEFID)
  {
    varInfoListUsed++;
    varID = paramNewEntry(param);
    VarInfo *var = &varInfoList[varID];
    var->gridID = gridID;
    var->zaxistype = zaxistype;
    var->ltype1 = ltype1;
    var->ltype2 = ltype2;
    var->hasBounds = hasBounds;
    var->level_sf = level_sf;
    var->level_unit = level_unit;
    var->tsteptype = tsteptype;
    if (scanKeys) var->scanKeys = *scanKeys;

    if (name && name[0]) var->name = strdup(name);
  }
  else
  {
    char paramstr[32];
    cdiParamToString(param, paramstr, sizeof(paramstr));

    VarInfo *var = &varInfoList[varID];
    if (var->gridID != gridID)
    {
      Message("param = %s gridID = %d", paramstr, gridID);
      Error("horizontal grid must not change for same parameter!");
    }
    if (var->zaxistype != zaxistype)
    {
      Message("param = %s zaxistype = %d", paramstr, zaxistype);
      Error("zaxistype must not change for same parameter!");
    }
  }

  VarInfo *var = &varInfoList[varID];
  if (prec > var->prec) var->prec = prec;

  // append current tile to tile subtype info.
  int this_tile = varInsertTileSubtype(&varInfoList[varID], tiles);
  int tileID = tileGetEntry(varID, this_tile);
  if (tile_index) (*tile_index) = this_tile;
  if (tileID == CDI_UNDEFID)
  {
    tileID = tileNewEntry((int) varID);
    var->recordTable[tileID].subtypeIndex = this_tile;
    var->nsubtypes++;
  }

  // append current level to level table info
  int levelID = levelNewEntry(varID, level1, level2, tileID);
  if (CDI_Debug)
    Message("vartable[%d].recordTable[%d].levelTable[%d].recID = %d; level1,2=%d,%d", varID, tileID, levelID, recID, level1,
            level2);
  var->recordTable[tileID].levelTable[levelID].recID = recID;

  *pvarID = (int) varID;
  *plevelID = levelID;
}

/*
static
int dblcmp(const void *s1, const void *s2)
{
  int cmp = 0;

  if      ( *((double *) s1) < *((double *) s2) ) cmp = -1;
  else if ( *((double *) s1) > *((double *) s2) ) cmp =  1;

  return cmp;
}
*/
static int
cmpLevelTable(const void *s1, const void *s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t *) s1;
  const leveltable_t *y = (const leveltable_t *) s2;
  // printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  if (x->level1 < y->level1)
    cmp = -1;
  else if (x->level1 > y->level1)
    cmp = 1;

  return cmp;
}

static int
cmpLevelTableInv(const void *s1, const void *s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t *) s1;
  const leveltable_t *y = (const leveltable_t *) s2;
  // printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  if (x->level1 < y->level1)
    cmp = 1;
  else if (x->level1 > y->level1)
    cmp = -1;

  return cmp;
}

void
varCopyKeys(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  cdiInitKeys(&vlistptr->vars[varID].keys);
  cdiCopyVarKeys(&(varInfoList[varID].keys), &(vlistptr->vars[varID].keys));
}
/*
struct cdi_generate_varinfo
{
  int varid;
  const char *name;
};

static int
cdi_generate_cmp_varname(const void *s1, const void *s2)
{
  const struct cdi_generate_varinfo *x = (const struct cdi_generate_varinfo *) s1, *y = (const struct cdi_generate_varinfo *) s2;
  return strcmp(x->name, y->name);
}
*/
void
cdi_generate_vars(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;

  int *varids = (int *) Malloc((size_t) varInfoListUsed * sizeof(int));
  for (int varID = 0; varID < varInfoListUsed; varID++) varids[varID] = (int) varID;
  /*
    if (streamptr->sortname)
      {
        bool hasName = true;
        for (int varID = 0; varID < varTableUsed; varID++)
          if (!vartable[varID].name) hasName = false;

        if (hasName)
          {
            struct cdi_generate_varinfo *varInfo
                = (struct cdi_generate_varinfo *) Malloc((size_t) varTableUsed * sizeof(struct cdi_generate_varinfo));

            for (int varID = 0; varID < varTableUsed; varID++)
              {
                varInfo[varID].varid = varids[varID];
                varInfo[varID].name = vartable[varids[varID]].name;
              }
            qsort(varInfo, varTableUsed, sizeof(varInfo[0]), cdi_generate_cmp_varname);
            for (int varID = 0; varID < varTableUsed; varID++)
              {
                varids[varID] = varInfo[varID].varid;
              }
            Free(varInfo);
          }
      }
  */
  for (int index = 0; index < varInfoListUsed; index++)
  {
    int varid = varids[index];
    VarInfo *var = &varInfoList[varid];

    int gridID = var->gridID;
    int param = var->param;
    int ltype1 = var->ltype1;
    int ltype2 = var->ltype2;
    int zaxistype = var->zaxistype;
    if (ltype1 == 0 && zaxistype == ZAXIS_GENERIC && cdiDefaultLeveltype != -1) zaxistype = cdiDefaultLeveltype;
    int hasBounds = var->hasBounds;
    int prec = var->prec;
    int instID = var->instID;
    int modelID = var->modelID;
    int tableID = var->tableID;
    int tsteptype = var->tsteptype;
    int comptype = var->comptype;

    double level_sf = (var->level_sf != 0) ? (1.0 / var->level_sf) : 1;

    /* consistency check: test if all subtypes have the same levels: */
    int nlevels = var->recordTable[0].nlevels;
    for (int isub = 1; isub < var->nsubtypes; isub++)
    {
      if (var->recordTable[isub].nlevels != nlevels)
      {
        fprintf(stderr,
                "var \"%s\": isub = %d / %d :: "
                "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                var->name, isub, var->nsubtypes, nlevels, var->recordTable[isub].nlevels);
        Error("zaxis size must not change for same parameter!");
      }

      const leveltable_t *t1 = var->recordTable[isub - 1].levelTable;
      const leveltable_t *t2 = var->recordTable[isub].levelTable;
      for (int ilev = 0; ilev < nlevels; ilev++)
        if ((t1[ilev].level1 != t2[ilev].level1) || (t1[ilev].level2 != t2[ilev].level2) || (t1[ilev].lindex != t2[ilev].lindex))
        {
          fprintf(stderr,
                  "var \"%s\", varID=%d: isub = %d / %d :: "
                  "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                  var->name, varid, isub, var->nsubtypes, nlevels, var->recordTable[isub].nlevels);
          Message("t1[ilev].level1=%d / t2[ilev].level1=%d", t1[ilev].level1, t2[ilev].level1);
          Message("t1[ilev].level2=%d / t2[ilev].level2=%d", t1[ilev].level2, t2[ilev].level2);
          Message("t1[ilev].lindex=%d / t2[ilev].lindex=%d", t1[ilev].lindex, t2[ilev].lindex);
          Error("zaxis type must not change for same parameter!");
        }
    }
    leveltable_t *levelTable = var->recordTable[0].levelTable;

    if (ltype1 == 0 && zaxistype == ZAXIS_GENERIC && nlevels == 1 && levelTable[0].level1 == 0) zaxistype = ZAXIS_SURFACE;

    double *dlevels = (double *) Malloc((size_t) nlevels * sizeof(double));

    /*
    if ( hasBounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
      for (int levelID = 0; levelID < nlevels; levelID++)
        dlevels[levelID] = (level_sf*levelTable[levelID].level1 + level_sf*levelTable[levelID].level2) / 2.0;
    else
    */
    for (int levelID = 0; levelID < nlevels; levelID++) dlevels[levelID] = level_sf * levelTable[levelID].level1;

    if (nlevels > 1)
    {
      bool linc = true, ldec = true, lsort = false;
      for (int levelID = 1; levelID < nlevels; levelID++)
      {
        // check increasing of levels
        linc &= (dlevels[levelID] > dlevels[levelID - 1]);
        // check decreasing of levels
        ldec &= (dlevels[levelID] < dlevels[levelID - 1]);
      }
      /*
       * always sort pressure z-axis to ensure
       * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 > levelID2
       * unless already sorted in decreasing order
       */
      if ((!linc && !ldec) && zaxistype == ZAXIS_PRESSURE)
      {
        qsort(levelTable, (size_t) nlevels, sizeof(leveltable_t), cmpLevelTableInv);
        lsort = true;
      }
      /*
       * always sort hybrid and depth-below-land z-axis to ensure
       * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 < levelID2
       * unless already sorted in increasing order
       */
      else if ((!linc && !ldec) || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_DEPTH_BELOW_LAND)
      {
        qsort(levelTable, (size_t) nlevels, sizeof(leveltable_t), cmpLevelTable);
        lsort = true;
      }

      if (lsort)
      {
        /*
        if ( hasBounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
          for (int levelID = 0; levelID < nlevels; levelID++)
            dlevels[levelID] = (level_sf*levelTable[levelID].level1 + level_sf*levelTable[levelID].level2) / 2.0;
        else
        */
        for (int levelID = 0; levelID < nlevels; levelID++) dlevels[levelID] = level_sf * levelTable[levelID].level1;
      }
    }

    double *dlevels1 = NULL, *dlevels2 = NULL;
    if (hasBounds)
    {
      dlevels1 = (double *) Malloc(2 * (size_t) nlevels * sizeof(double));
      dlevels2 = dlevels1 + nlevels;
      for (int levelID = 0; levelID < nlevels; levelID++) dlevels1[levelID] = level_sf * levelTable[levelID].level1;
      for (int levelID = 0; levelID < nlevels; levelID++) dlevels2[levelID] = level_sf * levelTable[levelID].level2;
    }

    const char **cvals = NULL;
    const char *unitptr = cdiUnitNamePtr(var->level_unit);
    int zaxisID = varDefZaxis(vlistID, zaxistype, nlevels, dlevels, cvals, 0, dlevels1, dlevels2, (int) Vctsize, Vct, NULL, NULL,
                              unitptr, 0, 0, ltype1, ltype2);

    if (CDI_CMOR_Mode && nlevels == 1 && zaxistype != ZAXIS_HYBRID) zaxisDefScalar(zaxisID);

    if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
    {
      if (numberOfVerticalLevels > 0) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NLEV, numberOfVerticalLevels);
      if (numberOfVerticalGrid > 0) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, numberOfVerticalGrid);
      if (!cdiUUIDIsNull(uuidVGrid)) cdiDefKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuidVGrid, CDI_UUID_SIZE);
    }

    if (hasBounds) Free(dlevels1);
    Free(dlevels);

    // define new subtype for tile set
    int tilesetID = CDI_UNDEFID;
    if (var->tiles) tilesetID = vlistDefTileSubtype(vlistID, var->tiles);

    // generate new variable
    (void) stream_new_var(streamptr, gridID, zaxisID, tilesetID);
    int varID = vlistDefVarTiles(vlistID, gridID, zaxisID, TIME_VARYING, tilesetID);

    vlistDefVarTsteptype(vlistID, varID, tsteptype);
    vlistDefVarParam(vlistID, varID, param);
    vlistDefVarDatatype(vlistID, varID, prec);
    vlistDefVarCompType(vlistID, varID, comptype);

    varCopyKeys(vlistID, varID);

    if (var->lmissval) vlistDefVarMissval(vlistID, varID, var->missval);
    if (var->name) cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, var->name);

    vlist_t *vlistptr = vlist_to_pointer(vlistID);
    for (int i = 0; i < var->opt_grib_nentries; i++)
    {
      resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
      vlistptr->vars[varID].opt_grib_nentries += 1;
      int idx = vlistptr->vars[varID].opt_grib_nentries - 1;

      vlistptr->vars[varID].opt_grib_kvpair[idx] = var->opt_grib_kvpair[i];
      vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = NULL;
      if (var->opt_grib_kvpair[i].keyword)
        vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdup(var->opt_grib_kvpair[i].keyword);
      vlistptr->vars[varID].opt_grib_kvpair[i].update = true;
    }
    // note: if the key is not defined, we do not throw an error!

    if (CDI_Default_TableID != CDI_UNDEFID)
    {
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      char name[CDI_MAX_NAME];
      name[0] = 0;
      char longname[CDI_MAX_NAME];
      longname[0] = 0;
      char units[CDI_MAX_NAME];
      units[0] = 0;
      tableInqEntry(CDI_Default_TableID, pnum, -1, name, longname, units);
      if (name[0])
      {
        if (tableID != CDI_UNDEFID)
        {
          cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name);
          if (longname[0]) cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname);
          if (units[0]) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units);
        }
        else
          tableID = CDI_Default_TableID;
      }
      if (CDI_Default_ModelID != CDI_UNDEFID) modelID = CDI_Default_ModelID;
      if (CDI_Default_InstID != CDI_UNDEFID) instID = CDI_Default_InstID;
    }

    if (instID != CDI_UNDEFID) vlistDefVarInstitut(vlistID, varID, instID);
    if (modelID != CDI_UNDEFID) vlistDefVarModel(vlistID, varID, modelID);
    if (tableID != CDI_UNDEFID) vlistDefVarTable(vlistID, varID, tableID);
  }

  for (int index = 0; index < varInfoListUsed; index++)
  {
    int varid = varids[index];
    VarInfo *var = &varInfoList[varid];
    int nlevels = var->recordTable[0].nlevels;

    int nsub = (var->nsubtypes >= 0) ? var->nsubtypes : 0;
    for (int isub = 0; isub < nsub; isub++)
    {
      sleveltable_t *streamRecordTable = streamptr->vars[index].recordTable + isub;
      leveltable_t *vartableLevelTable = var->recordTable[isub].levelTable;
      for (int levelID = 0; levelID < nlevels; levelID++)
      {
        streamRecordTable->recordID[levelID] = vartableLevelTable[levelID].recID;
        int lindex;
        for (lindex = 0; lindex < nlevels; lindex++)
          if (levelID == vartableLevelTable[lindex].lindex) break;
        if (lindex == nlevels) Error("Internal problem! lindex not found.");
        streamRecordTable->lindex[levelID] = (int) lindex;
      }
    }
  }

  Free(varids);

  varFree();
}

void
varDefVCT(size_t vctsize, double *vctptr)
{
  if (Vct == NULL && vctptr != NULL && vctsize > 0)
  {
    Vctsize = vctsize;
    Vct = (double *) Malloc(vctsize * sizeof(double));
    memcpy(Vct, vctptr, vctsize * sizeof(double));
  }
}

void
varDefZAxisReference(int nhlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE])
{
  numberOfVerticalLevels = nhlev;
  numberOfVerticalGrid = nvgrid;
  memcpy(uuidVGrid, uuid, CDI_UUID_SIZE);
}

bool
zaxis_compare(int zaxisID, int zaxistype, int nlevels, const double *levels, const double *lbounds, const double *ubounds,
              const char *longname, const char *units, int ltype1, int ltype2)
{
  bool differ = true;

  int ltype1_0 = 0, ltype2_0 = -1;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype1_0);
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, &ltype2_0);
  bool ltype1IsEqual = (ltype1 == ltype1_0);
  bool ltype2IsEqual = (ltype2 == ltype2_0);
  bool hasBounds = (lbounds && ubounds);

  if (ltype1IsEqual && ltype2IsEqual && (zaxistype == zaxisInqType(zaxisID) || zaxistype == ZAXIS_GENERIC))
  {
    bool hasBoundsZ = (zaxisInqLbounds(zaxisID, NULL) > 0 && zaxisInqUbounds(zaxisID, NULL) > 0);
    if (nlevels == zaxisInqSize(zaxisID) && hasBoundsZ == hasBounds)
    {
      const double *dlevels = zaxisInqLevelsPtr(zaxisID);
      if (dlevels && levels)
      {
        int levelID;
        for (levelID = 0; levelID < nlevels; levelID++)
          if (fabs(dlevels[levelID] - levels[levelID]) > 1.e-9) break;
        if (levelID == nlevels) differ = false;
      }

      if (!differ && hasBounds)
      {
        double *bounds = (double *) Malloc(2 * (size_t) nlevels * sizeof(double));
        zaxisInqLbounds(zaxisID, bounds);
        zaxisInqUbounds(zaxisID, bounds + nlevels);
        for (int levelID = 0; levelID < nlevels; levelID++)
        {
          if (fabs(lbounds[levelID] - bounds[levelID]) > 1.e-9 || fabs(ubounds[levelID] - bounds[levelID + nlevels]) > 1.e-9)
          {
            differ = true;
            break;
          }
        }
        Free(bounds);
      }

      if (!differ)
      {
        if (longname && longname[0])
        {
          char zlongname[CDI_MAX_NAME];
          int length = CDI_MAX_NAME;
          cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, zlongname, &length);
          if (zlongname[0] && !str_is_equal(longname, zlongname)) differ = true;
        }
        if (units && units[0])
        {
          char zunits[CDI_MAX_NAME];
          int length = CDI_MAX_NAME;
          cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zunits, &length);
          if (zunits[0] && !str_is_equal(units, zunits)) differ = true;
        }
      }
    }
  }

  return differ;
}

struct varDefZAxisSearchState
{
  int resIDValue;
  int zaxistype;
  int nlevels;
  const double *levels;
  const double *lbounds;
  const double *ubounds;
  const char *longname;
  const char *units;
  int ltype1;
  int ltype2;
};

static enum cdiApplyRet
varDefZAxisSearch(int id, void *res, void *data)
{
  struct varDefZAxisSearchState *state = (struct varDefZAxisSearchState *) data;
  (void) res;
  if (zaxis_compare(id, state->zaxistype, state->nlevels, state->levels, state->lbounds, state->ubounds, state->longname,
                    state->units, state->ltype1, state->ltype2)
      == false)
  {
    state->resIDValue = id;
    return CDI_APPLY_STOP;
  }
  else
    return CDI_APPLY_GO_ON;
}

int
varDefZaxis(int vlistID, int zaxistype, int nlevels, const double *levels, const char **cvals, size_t clength,
            const double *levels1, const double *levels2, int vctsize, const double *vct, char *name, const char *longname,
            const char *units, int prec, int mode, int ltype1, int ltype2)
{
  /*
    mode: 0 search in vlist and zaxis table
          1 search in zaxis table
   */
  int zaxisID = CDI_UNDEFID;
  bool zaxisdefined = false;
  bool zaxisglobdefined = false;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int nzaxis = vlistptr->nzaxis;

  if (ltype2 == 255) ltype2 = -1;

  if (mode == 0)
    for (int index = 0; index < nzaxis; index++)
    {
      zaxisID = vlistptr->zaxisIDs[index];

      if (!zaxis_compare(zaxisID, zaxistype, nlevels, levels, levels1, levels2, longname, units, ltype1, ltype2))
      {
        zaxisdefined = true;
        break;
      }
    }

  if (!zaxisdefined)
  {
    struct varDefZAxisSearchState query;
    query.zaxistype = zaxistype;
    query.nlevels = nlevels;
    query.levels = levels;
    query.lbounds = levels1;
    query.ubounds = levels2;
    query.longname = longname;
    query.units = units;
    query.ltype1 = ltype1;
    query.ltype2 = ltype2;

    if ((zaxisglobdefined = (cdiResHFilterApply(getZaxisOps(), varDefZAxisSearch, &query) == CDI_APPLY_STOP)))
      zaxisID = query.resIDValue;

    if (mode == 1 && zaxisglobdefined)
      for (int index = 0; index < nzaxis; index++)
        if (vlistptr->zaxisIDs[index] == zaxisID)
        {
          zaxisglobdefined = false;
          break;
        }
  }

  if (!zaxisdefined)
  {
    if (!zaxisglobdefined)
    {
      zaxisID = zaxisCreate(zaxistype, nlevels);
      if (levels) zaxisDefLevels(zaxisID, levels);
      if (levels1 && levels2)
      {
        zaxisDefLbounds(zaxisID, levels1);
        zaxisDefUbounds(zaxisID, levels2);
      }

      if (cvals != NULL && nlevels != 0 && clength != 0) zaxisDefCvals(zaxisID, cvals, (int) clength);

      if ((zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && vctsize > 0) zaxisDefVct(zaxisID, vctsize, vct);

      if (name && name[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name);
      if (longname && longname[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname);
      if (units && units[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units);
      zaxisDefDatatype(zaxisID, prec);
      cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, ltype1);
      if (ltype2 != -1) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, ltype2);
    }

    vlistptr->zaxisIDs[nzaxis] = zaxisID;
    vlistptr->nzaxis++;
  }

  return zaxisID;
}

void
varDefMissval(int varID, double missval)
{
  varInfoList[varID].lmissval = true;
  varInfoList[varID].missval = missval;
}

void
varDefCompType(int varID, int comptype)
{
  if (varInfoList[varID].comptype == CDI_COMPRESS_NONE) varInfoList[varID].comptype = comptype;
}

void
varDefCompLevel(int varID, int complevel)
{
  varInfoList[varID].complevel = complevel;
}

int
varInqInst(int varID)
{
  return varInfoList[varID].instID;
}

void
varDefInst(int varID, int instID)
{
  varInfoList[varID].instID = instID;
}

int
varInqModel(int varID)
{
  return varInfoList[varID].modelID;
}

void
varDefModel(int varID, int modelID)
{
  varInfoList[varID].modelID = modelID;
}

int
varInqTable(int varID)
{
  return varInfoList[varID].tableID;
}

void
varDefTable(int varID, int tableID)
{
  varInfoList[varID].tableID = tableID;
}

void
varDefKeyInt(int varID, int key, int value)
{
  cdi_keys_t *keysp = &(varInfoList[varID].keys);
  cdiDefVarKeyInt(keysp, key, value);
}

void
varDefKeyBytes(int varID, int key, const unsigned char *bytes, int length)
{
  cdi_keys_t *keysp = &(varInfoList[varID].keys);
  cdiDefVarKeyBytes(keysp, key, bytes, length);
}

void
varDefKeyString(int varID, int key, const char *string)
{
  int length = (int) strlen(string) + 1;
  cdi_keys_t *keysp = &(varInfoList[varID].keys);
  cdiDefVarKeyBytes(keysp, key, (const unsigned char *) string, length);
}

#ifdef HAVE_LIBGRIB_API
// Resizes and initializes opt_grib_kvpair data structure.
static void
resize_vartable_opt_grib_entries(VarInfo *var, int nentries)
{
  if (var->opt_grib_kvpair_size < nentries)
  {
    if (CDI_Debug) Message("resize data structure, %d -> %d", var->opt_grib_kvpair_size, nentries);

    int new_size = ((2 * var->opt_grib_kvpair_size) > nentries) ? (2 * var->opt_grib_kvpair_size) : nentries;
    if (CDI_Debug) Message("resize vartable opt_grib_entries array to size %d", new_size);
    opt_key_val_pair_t *tmp = (opt_key_val_pair_t *) Malloc((size_t) new_size * sizeof(opt_key_val_pair_t));
    for (int i = 0; i < var->opt_grib_kvpair_size; i++) { tmp[i] = var->opt_grib_kvpair[i]; }
    for (int i = var->opt_grib_kvpair_size; i < new_size; i++)
    {
      tmp[i].int_val = 0;
      tmp[i].dbl_val = 0;
      tmp[i].int_arr = 0;
      tmp[i].dbl_arr = 0;
      tmp[i].update = false;
      tmp[i].keyword = NULL;
    }  // for
    var->opt_grib_kvpair_size = new_size;
    Free(var->opt_grib_kvpair);
    var->opt_grib_kvpair = tmp;
  }
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribInt(int varID, int tile_index, long lval, const char *keyword)
{
  VarInfo *var = &varInfoList[varID];

  int idx = -1;
  for (int i = 0; i < var->opt_grib_nentries; i++)
  {
    if (str_is_equal(keyword, var->opt_grib_kvpair[i].keyword) && (var->opt_grib_kvpair[i].data_type == t_int)
        && (var->opt_grib_kvpair[i].subtype_index == tile_index))
      idx = i;
  }

  if (idx == -1)
  {
    resize_vartable_opt_grib_entries(&varInfoList[varID], var->opt_grib_nentries + 1);
    var->opt_grib_nentries += 1;
    idx = var->opt_grib_nentries - 1;
  }
  else
  {
    if (var->opt_grib_kvpair[idx].keyword) Free(var->opt_grib_kvpair[idx].keyword);
  }
  var->opt_grib_kvpair[idx].data_type = t_int;
  var->opt_grib_kvpair[idx].int_val = (int) lval;
  var->opt_grib_kvpair[idx].keyword = strdup(keyword);
  var->opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribDbl(int varID, int tile_index, double dval, const char *keyword)
{
  VarInfo *var = &varInfoList[varID];
  int idx = -1;
  for (int i = 0; i < var->opt_grib_nentries; i++)
  {
    if (str_is_equal(keyword, var->opt_grib_kvpair[i].keyword) && (var->opt_grib_kvpair[i].data_type == t_double)
        && (var->opt_grib_kvpair[i].subtype_index == tile_index))
      idx = i;
  }

  if (idx == -1)
  {
    resize_vartable_opt_grib_entries(&varInfoList[varID], var->opt_grib_nentries + 1);
    var->opt_grib_nentries += 1;
    idx = var->opt_grib_nentries - 1;
  }
  else
  {
    if (var->opt_grib_kvpair[idx].keyword) Free(var->opt_grib_kvpair[idx].keyword);
  }
  var->opt_grib_kvpair[idx].data_type = t_double;
  var->opt_grib_kvpair[idx].dbl_val = dval;
  var->opt_grib_kvpair[idx].keyword = strdup(keyword);
  var->opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribIntArr(int varID, int tile_index, const long *larr, size_t arr_len, const char *keyword)
{
  VarInfo *var = &varInfoList[varID];
  int idx = -1;

  // Search for an existing entry with the same keyword and tile_index
  for (int i = 0; i < var->opt_grib_nentries; i++)
  {
    if (str_is_equal(keyword, var->opt_grib_kvpair[i].keyword) && (var->opt_grib_kvpair[i].data_type == t_intarr)
        && (var->opt_grib_kvpair[i].subtype_index == tile_index))
    {
      idx = i;
      break;
    }
  }

  // If entry not found, create a new one
  if (idx == -1)
  {
    resize_vartable_opt_grib_entries(&varInfoList[varID], var->opt_grib_nentries + 1);
    var->opt_grib_nentries += 1;
    idx = var->opt_grib_nentries - 1;
  }
  else
  {
    // Free existing keyword and array
    if (var->opt_grib_kvpair[idx].keyword) Free(var->opt_grib_kvpair[idx].keyword);
    if (var->opt_grib_kvpair[idx].int_arr) Free(var->opt_grib_kvpair[idx].int_arr);
  }

  // Assign keyword, data type, and tile index
  var->opt_grib_kvpair[idx].data_type = t_intarr;
  var->opt_grib_kvpair[idx].keyword = strdup(keyword);
  var->opt_grib_kvpair[idx].subtype_index = tile_index;

  // Allocate memory for integer array
  var->opt_grib_kvpair[idx].arr_len = arr_len;
  if (arr_len > 0)
  {
    var->opt_grib_kvpair[idx].int_arr = (int *) Malloc(arr_len * sizeof(int));
    for (size_t i = 0; i < arr_len; i++)
    {
      var->opt_grib_kvpair[idx].int_arr[i] = (int) larr[i];  // convert long to int
    }
  }
  else { var->opt_grib_kvpair[idx].int_arr = NULL; }
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribDblArr(int varID, int tile_index, const double *darr, size_t arr_len, const char *keyword)
{
  VarInfo *var = &varInfoList[varID];
  int idx = -1;

  // Search for an existing entry with the same keyword and tile_index
  for (int i = 0; i < var->opt_grib_nentries; i++)
  {
    if (str_is_equal(keyword, var->opt_grib_kvpair[i].keyword) && (var->opt_grib_kvpair[i].data_type == t_doublearr)
        && (var->opt_grib_kvpair[i].subtype_index == tile_index))
    {
      idx = i;
      break;
    }
  }

  // If entry not found, create a new one
  if (idx == -1)
  {
    resize_vartable_opt_grib_entries(&varInfoList[varID], var->opt_grib_nentries + 1);
    var->opt_grib_nentries += 1;
    idx = var->opt_grib_nentries - 1;
  }
  else
  {
    // Free existing keyword and array
    if (var->opt_grib_kvpair[idx].keyword) Free(var->opt_grib_kvpair[idx].keyword);
    if (var->opt_grib_kvpair[idx].dbl_arr) Free(var->opt_grib_kvpair[idx].dbl_arr);
  }

  // Assign keyword, data type, and tile index
  var->opt_grib_kvpair[idx].data_type = t_doublearr;
  var->opt_grib_kvpair[idx].keyword = strdup(keyword);
  var->opt_grib_kvpair[idx].subtype_index = tile_index;

  // Allocate memory for double array
  var->opt_grib_kvpair[idx].arr_len = arr_len;
  if (arr_len > 0)
  {
    var->opt_grib_kvpair[idx].dbl_arr = (double *) Malloc(arr_len * sizeof(double));
    for (size_t i = 0; i < arr_len; i++) { var->opt_grib_kvpair[idx].dbl_arr[i] = darr[i]; }
  }
  else { var->opt_grib_kvpair[idx].dbl_arr = NULL; }
}

#endif

#ifdef HAVE_LIBGRIB_API
int
varOptGribNentries(int varID)
{
  return varInfoList[varID].opt_grib_nentries;
}
#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
