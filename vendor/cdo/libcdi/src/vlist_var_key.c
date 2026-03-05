#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"

/* vlistDefVarIntKey: Set an arbitrary keyword/integer value pair for GRIB API */
void
vlistDefVarIntKey(int vlistID, int varID, const char *name, int value)
{
#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL) Error("Internal error!");
  int idx;

  if (vlistptr->immutable)
    Error("vlistDefVarIntKey() was called on an immutable vlist object (vlistID = %d)\n"
          "Either call vlistDefVarIntKey() before passing the vlist object to streamDefVlist(),\n"
          "or use the stream-internal vlist by calling streamInqVlist().",
          vlistID);

  for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int))
      break;

  if (idx < vlistptr->vars[varID].opt_grib_nentries)
  {
    vlistptr->vars[varID].opt_grib_kvpair[idx].int_val = value;
    vlistptr->vars[varID].opt_grib_kvpair[idx].update = true;
  }
  else
  {
    resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
    vlistptr->vars[varID].opt_grib_nentries += 1;
    idx = vlistptr->vars[varID].opt_grib_nentries - 1;
    vlistptr->vars[varID].opt_grib_kvpair[idx].data_type = t_int;
    vlistptr->vars[varID].opt_grib_kvpair[idx].int_val = value;
    vlistptr->vars[varID].opt_grib_kvpair[idx].update = true;
    if (name)
      vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdup(name);
    else
      Error("Internal error, name undefined!");
  }

  if (CDI_Debug)
  {
    Message("define additional GRIB2 key \"%s\" (integer): %d", name, value);
    Message("total list of registered, additional GRIB2 keys (total: %d):", vlistptr->vars[varID].opt_grib_nentries);
    for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
      if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int)
        Message("%s -> integer %d", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                vlistptr->vars[varID].opt_grib_kvpair[idx].int_val);
      else if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double)
        Message("%s -> double %d", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val);
      else
        Message("%s -> unknown", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword);
  }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void) vlistID;
  (void) varID;
  (void) name;
  (void) value;
#endif
}

/* vlistDefVarDblKey: Set an arbitrary keyword/double value pair for GRIB API */
void
vlistDefVarDblKey(int vlistID, int varID, const char *name, double value)
{
#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL) Error("Internal error!");
  int idx;

  if (vlistptr->immutable)
    Error("vlistDefVarDblKey() was called on an immutable vlist object (vlistID = %d)\n"
          "Either call vlistDefVarIntKey() before passing the vlist object to streamDefVlist(),\n"
          "or use the stream-internal vlist by calling streamInqVlist().",
          vlistID);

  for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double))
      break;

  if (idx < vlistptr->vars[varID].opt_grib_nentries)
  {
    vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val = value;
    vlistptr->vars[varID].opt_grib_kvpair[idx].update = true;
  }
  else
  {
    resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
    vlistptr->vars[varID].opt_grib_nentries += 1;
    idx = vlistptr->vars[varID].opt_grib_nentries - 1;
    vlistptr->vars[varID].opt_grib_kvpair[idx].data_type = t_double;
    vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val = value;
    vlistptr->vars[varID].opt_grib_kvpair[idx].update = true;
    if (name)
      vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdup(name);
    else
      Error("Internal error, name undefined!");
  }

  if (CDI_Debug)
  {
    Message("define additional GRIB2 key \"%s\" (double): %d", name, value);
    Message("total list of registered, additional GRIB2 keys (total: %d):", vlistptr->vars[varID].opt_grib_nentries);
    for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
      if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int)
        Message("%s -> integer %d", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                vlistptr->vars[varID].opt_grib_kvpair[idx].int_val);
      else if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double)
        Message("%s -> double %d", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val);
      else
        Message("%s -> unknown", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword);
  }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void) vlistID;
  (void) varID;
  (void) name;
  (void) value;
#endif
}

void
vlistDefVarIntArrKey(int vlistID, int varID, const char *name, int *values, int nvalues)
{
#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL) Error("Internal error!");
  if (name == NULL) Error("Internal error, name undefined!");

  if (vlistptr->immutable)
    Error("vlistDefVarIntArrKey() was called on an immutable vlist object (vlistID = %d)\n"
          "Either call this before passing the vlist object to streamDefVlist(),\n"
          "or use the stream-internal vlist by calling streamInqVlist().",
          vlistID);

  int idx;
  for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_intarr))
      break;

  opt_key_val_pair_t *entry = NULL;

  if (idx < vlistptr->vars[varID].opt_grib_nentries)
  {
    entry = &vlistptr->vars[varID].opt_grib_kvpair[idx];
    free(entry->int_arr);
  }
  else
  {
    resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
    vlistptr->vars[varID].opt_grib_nentries += 1;
    idx = vlistptr->vars[varID].opt_grib_nentries - 1;
    entry = &vlistptr->vars[varID].opt_grib_kvpair[idx];

    entry->keyword = strdup(name);
    entry->data_type = t_intarr;
  }

  entry->update = true;
  entry->arr_len = nvalues;
  entry->int_arr = (int *) malloc(nvalues * sizeof(int));
  if (!entry->int_arr) Error("Memory allocation failed for int_arr");
  memcpy(entry->int_arr, values, nvalues * sizeof(int));

  entry->dbl_arr = NULL;  // safety

  if (CDI_Debug)
  {
    Message("define additional GRIB2 key \"%s\" (integer array, len=%zu)", name, nvalues);
    for (int i = 0; i < nvalues; ++i) Message("  [%zu] = %d", i, entry->int_arr[i]);
  }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void) vlistID;
  (void) varID;
  (void) name;
  (void) values;
  (void) nvalues;
#endif
}

void
vlistDefVarDblArrKey(int vlistID, int varID, const char *name, double *values, int nvalues)
{
#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL) Error("Internal error!");
  if (name == NULL) Error("Internal error, name undefined!");

  if (vlistptr->immutable)
    Error("vlistDefVarDblArrKey() was called on an immutable vlist object (vlistID = %d)\n"
          "Either call this before passing the vlist object to streamDefVlist(),\n"
          "or use the stream-internal vlist by calling streamInqVlist().",
          vlistID);

  int idx;
  for (idx = 0; idx < vlistptr->vars[varID].opt_grib_nentries; idx++)
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_doublearr))
      break;

  opt_key_val_pair_t *entry;

  if (idx < vlistptr->vars[varID].opt_grib_nentries)
  {
    entry = &vlistptr->vars[varID].opt_grib_kvpair[idx];
    free(entry->dbl_arr);
  }
  else
  {
    resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
    vlistptr->vars[varID].opt_grib_nentries += 1;
    idx = vlistptr->vars[varID].opt_grib_nentries - 1;
    entry = &vlistptr->vars[varID].opt_grib_kvpair[idx];

    entry->keyword = strdup(name);
    entry->data_type = t_doublearr;
  }

  entry->update = true;
  entry->arr_len = nvalues;
  entry->dbl_arr = (double *) malloc(nvalues * sizeof(double));
  if (!entry->dbl_arr) Error("Memory allocation failed for dbl_arr");
  memcpy(entry->dbl_arr, values, nvalues * sizeof(double));

  entry->int_arr = NULL;  // safety

  if (CDI_Debug)
  {
    Message("define additional GRIB2 key \"%s\" (double array, len=%zu)", name, nvalues);
    for (int i = 0; i < nvalues; ++i) Message("  [%zu] = %g", i, entry->dbl_arr[i]);
  }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void) vlistID;
  (void) varID;
  (void) name;
  (void) values;
  (void) nvalues;
#endif
}

/* cdiClearAdditionalKeys: Clears the list of additional GRIB keys. */
void
cdiClearAdditionalKeys(void)
{
#ifdef HAVE_LIBGRIB_API
  for (int i = 0; i < cdiNAdditionalGRIBKeys; ++i) free(cdiAdditionalGRIBKeys[i]);
  cdiNAdditionalGRIBKeys = 0;
#endif
}

/* cdiDefAdditionalKey: Register an additional GRIB key which is read when file is opened. */
void
cdiDefAdditionalKey(const char *name)
{
#ifdef HAVE_LIBGRIB_API
  int idx = cdiNAdditionalGRIBKeys;
  cdiNAdditionalGRIBKeys++;
  if (idx >= MAX_OPT_GRIB_ENTRIES) Error("Too many additional keywords!");
  if (name)
    cdiAdditionalGRIBKeys[idx] = strdup(name);
  else
    Error("Internal error!");
#else
  (void) name;
#endif
}

/* vlistHasVarKey: returns 1 if meta-data key was read, 0 otherwise. */
int
vlistHasVarKey(int vlistID, int varID, const char *name)
{
#ifdef HAVE_LIBGRIB_API
  /* check if the GRIB key was previously read and is stored */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; ++i)
  {
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword)) return 1;
  }
#else
  (void) vlistID;
  (void) varID;
  (void) name;
#endif
  return 0;
}

/* vlistInqVarDblKey: raw access to GRIB meta-data */
double
vlistInqVarDblKey(int vlistID, int varID, const char *name)
{
  double value = 0;
#ifdef HAVE_LIBGRIB_API
  /* check if the GRIB key was previously read and is stored in "opt_grib_dbl_val" */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; ++i)
  {
    int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_double)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub))
      return vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val;
  }
#else
  (void) vlistID;
  (void) varID;
  (void) name;
#endif
  return value;
}

/* vlistInqVarIntKey: raw access to GRIB meta-data */
int
vlistInqVarIntKey(int vlistID, int varID, const char *name)
{
  long value = 0;
#ifdef HAVE_LIBGRIB_API
  /* check if the GRIB key was previously read and is stored in "opt_grib_int_val" */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; ++i)
  {
    int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_int)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub))
      return vlistptr->vars[varID].opt_grib_kvpair[i].int_val;
  }

#else
  (void) vlistID;
  (void) varID;
  (void) name;
#endif
  return (int) value;
}

/* vlistInqVarIntArrKey: raw access to GRIB integer array meta-data */
int *
vlistInqVarIntArrKey(int vlistID, int varID, const char *name)
{
  int *int_arr = NULL;

#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; ++i)
  {
    int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_intarr)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub))
    {
      int_arr = vlistptr->vars[varID].opt_grib_kvpair[i].int_arr;
      break;
    }
  }
#else
  (void) vlistID;
  (void) varID;
  (void) name;
#endif

  return int_arr;
}

/* vlistInqVarDblArrKey: raw access to GRIB double array meta-data */
double *
vlistInqVarDblArrKey(int vlistID, int varID, const char *name)
{
  double *dbl_arr = NULL;

#ifdef HAVE_LIBGRIB_API
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; ++i)
  {
    int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
    if (str_is_equal(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_doublearr)
        && (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub))
    {
      dbl_arr = vlistptr->vars[varID].opt_grib_kvpair[i].dbl_arr;
      break;
    }
  }
#else
  (void) vlistID;
  (void) varID;
  (void) name;
#endif

  return dbl_arr;
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
