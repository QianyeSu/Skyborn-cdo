#ifndef _EXTRA_H
#define _EXTRA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

enum
{
  EXT_REAL = 1,
  EXT_COMP = 2,
};

typedef struct
{
  int checked;
  int byteswap;
  int header[4];
  int prec;    // single or double precison
  int number;  // real or complex
  size_t datasize;
  size_t buffersize;
  void *buffer;
} extrec_t;

const char *extLibraryVersion(void);

void extDebug(int debug);

int extCheckFiletype(int fileID, int *swap);

void *extNew(void);
void extDelete(void *ext);

int extRead(int fileID, void *ext);
int extWrite(int fileID, void *ext);

int extInqHeader(void *ext, int *header);
int extInqDataFP32(void *ext, float *data);
int extInqDataFP64(void *ext, double *data);

int extDefHeader(void *ext, const int *header);
int extDefDataFP32(void *ext, const float *data);
int extDefDataFP64(void *ext, const double *data);

#ifdef HAVE__FLOAT16
int extInqDataFP16(void *ext, _Float16 *data);
int extDefDataFP16(void *ext, const _Float16 *data);
#endif

#endif /* _EXTRA_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
