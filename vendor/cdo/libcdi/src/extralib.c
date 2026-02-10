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

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "extra.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "exse.h"
#include "swap.h"

enum
{
  EXT_HEADER_LEN = 4,
};

union EXT_HEADER
{
  int32_t i32[EXT_HEADER_LEN];
  int64_t i64[EXT_HEADER_LEN];
};

static int initExtLib = 0;
static int extDefaultPrec = 0;
static int extDefaultNumber = EXT_REAL;

// A version string.
#undef LIBVERSION
#define LIBVERSION 2.0.0
#define XSTRING(x) #x
#define STRING(x) XSTRING(x)
static const char ext_libvers[] = STRING(LIBVERSION);

const char *
extLibraryVersion(void)
{
  return ext_libvers;
}

static int EXT_Debug = 0;  // If set to 1, debugging

void
extDebug(int debug)
{
  if (debug) Message("debug level %d", debug);
  EXT_Debug = debug;
}

static void
extLibInit(void)
{
  const char *envName = "EXT_PRECISION";

  char *envString = getenv(envName);
  if (envString)
  {
    if (strlen(envString) == 2)
    {
      switch (tolower((int) envString[0]))
      {
        case 'r':
        {
          extDefaultNumber = EXT_REAL;
          switch ((int) envString[1])
          {
            case '2': extDefaultPrec = EXSE_PREC_FP16; break;
            case '4': extDefaultPrec = EXSE_PREC_FP32; break;
            case '8': extDefaultPrec = EXSE_PREC_FP64; break;
            default: Warning("Invalid digit in %s: %s", envName, envString);
          }
          break;
        }
        case 'c':
        {
          extDefaultNumber = EXT_COMP;
          switch ((int) envString[1])
          {
            case '4': extDefaultPrec = EXSE_PREC_FP32; break;
            case '8': extDefaultPrec = EXSE_PREC_FP64; break;
            default: Warning("Invalid digit in %s: %s", envName, envString);
          }
          break;
        }
        default:
        {
          Warning("Invalid character in %s: %s", envName, envString);
          break;
        }
      }
    }
  }

  initExtLib = 1;
}

static void
extInit(extrec_t *extp)
{
  extp->checked = 0;
  extp->byteswap = 0;
  extp->prec = 0;
  extp->number = extDefaultNumber;
  extp->datasize = 0;
  extp->buffersize = 0;
  extp->buffer = NULL;
}

void *
extNew(void)
{
  if (!initExtLib) extLibInit();

  extrec_t *extp = (extrec_t *) Malloc(sizeof(extrec_t));

  extInit(extp);

  return (void *) extp;
}

void
extDelete(void *ext)
{
  extrec_t *extp = (extrec_t *) ext;

  if (extp)
  {
    if (extp->buffer) Free(extp->buffer);
    Free(extp);
  }
}

int
extCheckFiletype(int fileID, int *swap)
{
  size_t fact = 0;
  size_t data = 0;
  size_t dimxy = 0;
  unsigned char buffer[40], *pbuf;

  if (fileRead(fileID, buffer, 4) != 4) return 0;

  size_t blocklen = (size_t) get_uint32(buffer);
  size_t sblocklen = (size_t) get_swap_uint32(buffer);

  if (EXT_Debug) Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  // clang-format off
  if (blocklen == 16)
    {
     *swap = 0;
      fact = blocklen / 4;
      if (fileRead(fileID, buffer, blocklen+8) != blocklen+8) return 0;
      pbuf = buffer+3*fact;      dimxy = (size_t) get_uint32(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_uint32(pbuf);
    }
  else if (blocklen == 32)
    {
     *swap = 0;
      fact = blocklen / 4;
      if (fileRead(fileID, buffer, blocklen+8) != blocklen+8) return 0;
      pbuf = buffer+3*fact;      dimxy = (size_t) get_uint64(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_uint32(pbuf);
    }
  else if (sblocklen == 16)
    {
     *swap = 1;
      fact = sblocklen / 4;
      if (fileRead(fileID, buffer, sblocklen+8) != sblocklen+8) return 0;
      pbuf = buffer+3*fact;       dimxy = (size_t) get_swap_uint32(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_swap_uint32(pbuf);
    }
  else if (sblocklen == 32)
    {
     *swap = 1;
      fact = sblocklen / 4;
      if (fileRead(fileID, buffer, sblocklen+8) != sblocklen+8) return 0;
      pbuf = buffer+3*fact;       dimxy = (size_t) get_swap_uint64(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_swap_uint32(pbuf);
    }
  // clang-format on

  fileRewind(fileID);

  if (EXT_Debug) Message("swap = %d fact = %d", *swap, fact);
  if (EXT_Debug) Message("dimxy = %lu data = %lu", dimxy, data);

  int found = data && (dimxy * fact == data || dimxy * fact * 2 == data || dimxy * fact / 2 == data);
  return found;
}

int
extInqHeader(void *ext, int *header)
{
  extrec_t *extp = (extrec_t *) ext;

  for (int i = 0; i < EXT_HEADER_LEN; i++) header[i] = extp->header[i];

  if (EXT_Debug) Message("datasize = %zu", extp->datasize);

  return 0;
}

int
extDefHeader(void *ext, const int *header)
{
  extrec_t *extp = (extrec_t *) ext;

  for (int i = 0; i < EXT_HEADER_LEN; i++) extp->header[i] = header[i];

  extp->datasize = (size_t) header[3];
  if (extp->number == EXT_COMP) extp->datasize *= 2;

  if (EXT_Debug) Message("datasize = %zu", extp->datasize);

  return 0;
}

static int
extInqData(extrec_t *extp, int prec, void *data)
{
  int ierr = 0;
  int byteswap = extp->byteswap;
  size_t datasize = extp->datasize, buffer_size = datasize * (size_t) prec;
  void *buffer = extp->buffer;
  int rprec = extp->prec;

  switch (rprec)
  {
    case EXSE_PREC_FP32:
    {
      if (byteswap) swap4byte(buffer, datasize);

      if (EXSE_PREC_FP32 == prec)
        memcpy(data, buffer, buffer_size);
      else if (EXSE_PREC_FP64 == prec)
        for (size_t i = 0; i < datasize; ++i) ((double *) data)[i] = (double) ((float *) buffer)[i];
#ifdef HAVE__FLOAT16
      else if (EXSE_PREC_FP16 == prec)
        for (size_t i = 0; i < datasize; ++i) ((_Float16 *) data)[i] = (_Float16) ((float *) buffer)[i];
#endif
      break;
    }
    case EXSE_PREC_FP64:
    {
      if (byteswap) swap8byte(buffer, datasize);

      if (EXSE_PREC_FP64 == prec)
        memcpy(data, buffer, buffer_size);
      else if (EXSE_PREC_FP32 == prec)
        for (size_t i = 0; i < datasize; ++i) ((float *) data)[i] = (float) ((double *) buffer)[i];
#ifdef HAVE__FLOAT16
      else if (EXSE_PREC_FP16 == prec)
        for (size_t i = 0; i < datasize; ++i) ((_Float16 *) data)[i] = (_Float16) ((double *) buffer)[i];
#endif
      break;
    }
#ifdef HAVE__FLOAT16
    case EXSE_PREC_FP16:
    {
      if (EXSE_PREC_FP16 == prec)
        memcpy(data, buffer, buffer_size);
      else if (EXSE_PREC_FP64 == prec)
        for (size_t i = 0; i < datasize; ++i) ((double *) data)[i] = (double) ((_Float16 *) buffer)[i];
      else if (EXSE_PREC_FP32 == prec)
        for (size_t i = 0; i < datasize; ++i) ((float *) data)[i] = (float) ((_Float16 *) buffer)[i];

      break;
    }
#endif
    default:
    {
      Error("unexpected data precision %d", rprec);
      break;
    }
  }

  return ierr;
}

int
extInqDataFP32(void *ext, float *data)
{
  return extInqData((extrec_t *) ext, EXSE_PREC_FP32, (void *) data);
}

int
extInqDataFP64(void *ext, double *data)
{
  return extInqData((extrec_t *) ext, EXSE_PREC_FP64, (void *) data);
}
#ifdef HAVE__FLOAT16
int
extInqDataFP16(void *ext, _Float16 *data)
{
  return extInqData((extrec_t *) ext, EXSE_PREC_FP16, (void *) data);
}
#endif

static int
extDefData(void *ext, int prec, const void *data)
{
  extrec_t *extp = (extrec_t *) ext;

  int rprec = extDefaultPrec ? extDefaultPrec : extp->prec;
  extp->prec = rprec ? rprec : prec;

  int *header = extp->header;

  size_t datasize = (size_t) header[3];
  if (extp->number == EXT_COMP) datasize *= 2;
  size_t blocklen = datasize * (size_t) rprec;

  extp->datasize = datasize;

  if (extp->buffersize != blocklen)
  {
    extp->buffersize = blocklen;
    extp->buffer = Realloc(extp->buffer, extp->buffersize);
  }

  switch (rprec)
  {
    case EXSE_PREC_FP32:
    {
      if (EXSE_PREC_FP32 == prec)
        memcpy(extp->buffer, data, blocklen);
      else if (EXSE_PREC_FP64 == prec)
        for (size_t i = 0; i < datasize; i++) ((float *) extp->buffer)[i] = (float) ((double *) data)[i];
#ifdef HAVE__FLOAT16
      else if (EXSE_PREC_FP16 == prec)
        for (size_t i = 0; i < datasize; i++) ((float *) extp->buffer)[i] = (float) ((_Float16 *) data)[i];
#endif
      break;
    }
    case EXSE_PREC_FP64:
    {
      if (EXSE_PREC_FP64 == prec)
        memcpy(extp->buffer, data, blocklen);
      else if (EXSE_PREC_FP32 == prec)
        for (size_t i = 0; i < datasize; i++) ((double *) extp->buffer)[i] = (double) ((float *) data)[i];
#ifdef HAVE__FLOAT16
      else if (EXSE_PREC_FP16 == prec)
        for (size_t i = 0; i < datasize; i++) ((double *) extp->buffer)[i] = (double) ((_Float16 *) data)[i];
#endif
      break;
    }
#ifdef HAVE__FLOAT16
    case EXSE_PREC_FP16:
    {
      if (EXSE_PREC_FP16 == prec)
        memcpy(extp->buffer, data, blocklen);
      else if (EXSE_PREC_FP32 == prec)
        for (size_t i = 0; i < datasize; i++) ((_Float16 *) extp->buffer)[i] = (_Float16) ((float *) data)[i];
      else if (EXSE_PREC_FP64 == prec)
        for (size_t i = 0; i < datasize; i++) ((_Float16 *) extp->buffer)[i] = (_Float16) ((double *) data)[i];

      break;
    }
#endif
    default:
    {
      Error("unexpected data precision %d", rprec);
      break;
    }
  }

  return 0;
}

int
extDefDataFP32(void *ext, const float *data)
{
  return extDefData(ext, EXSE_PREC_FP32, (void *) data);
}

int
extDefDataFP64(void *ext, const double *data)
{
  return extDefData(ext, EXSE_PREC_FP64, (void *) data);
}
#ifdef HAVE__FLOAT16
int
extDefDataFP16(void *ext, const _Float16 *data)
{
  return extDefData(ext, EXSE_PREC_FP16, (void *) data);
}
#endif

int
extRead(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;

  if (!extp->checked)
  {
    int status = extCheckFiletype(fileID, &extp->byteswap);
    if (status == 0) Error("Not a EXTRA file!");
    extp->checked = 1;
  }

  int byteswap = extp->byteswap;

  // read header record
  size_t blocklen = binReadF77Block(fileID, byteswap);

  if (fileEOF(fileID)) return -1;

  if (EXT_Debug) Message("blocklen = %lu", blocklen);

  size_t hprec = blocklen / EXT_HEADER_LEN;
  // extp->prec = (int) hprec;

  union EXT_HEADER tempheader;
  switch (hprec)
  {
    case EXSE_PREC_FP32:
    case EXSE_PREC_FP16:
    {
      binReadInt32(fileID, byteswap, EXT_HEADER_LEN, tempheader.i32);
      for (int i = 0; i < EXT_HEADER_LEN; i++) extp->header[i] = (int) tempheader.i32[i];
      break;
    }
    case EXSE_PREC_FP64:
    {
      binReadInt64(fileID, byteswap, EXT_HEADER_LEN, tempheader.i64);
      for (int i = 0; i < EXT_HEADER_LEN; i++) extp->header[i] = (int) tempheader.i64[i];
      break;
    }
    default:
    {
      Error("Unexpected header precision %d", hprec);
      break;
    }
  }

  size_t blocklen2 = binReadF77Block(fileID, byteswap);

  if (blocklen2 != blocklen)
  {
    Warning("Header blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
    if (blocklen2 != 0) return -1;
  }

  extp->datasize = (size_t) extp->header[3];

  if (EXT_Debug) Message("datasize = %zu", extp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  if (extp->buffersize < blocklen)
  {
    extp->buffersize = blocklen;
    extp->buffer = Realloc(extp->buffer, extp->buffersize);
  }

  size_t dprec = blocklen / extp->datasize;
  extp->prec = (int) dprec;

  if (dprec == hprec || dprec == hprec / 2) { extp->number = EXT_REAL; }
  else if (dprec == hprec * 2)
  {
    dprec /= 2;
    extp->datasize *= 2;
    extp->number = EXT_COMP;
  }

  if (dprec != EXSE_PREC_FP32 && dprec != EXSE_PREC_FP64 && dprec != EXSE_PREC_FP16)
  {
    Warning("Unexpected data precision %d", dprec);
    return -1;
  }

  fileRead(fileID, extp->buffer, blocklen);

  blocklen2 = binReadF77Block(fileID, byteswap);

  if (blocklen2 != blocklen)
  {
    Warning("Data blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
    if (blocklen2 != 0) return -1;
  }

  return 0;
}

int
extWrite(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;
  union EXT_HEADER tempheader;
  int byteswap = extp->byteswap;
  int rprec = extp->prec;
  int number = extp->number;
  int *header = extp->header;

  // write header record
  size_t blocklen = EXT_HEADER_LEN * (size_t) ((rprec == EXSE_PREC_FP16) ? EXSE_PREC_FP32 : rprec);

  binWriteF77Block(fileID, byteswap, blocklen);

  switch (rprec)
  {
    case EXSE_PREC_FP16:
    case EXSE_PREC_FP32:
    {
      for (int i = 0; i < EXT_HEADER_LEN; i++) tempheader.i32[i] = (int32_t) header[i];
      binWriteInt32(fileID, byteswap, EXT_HEADER_LEN, tempheader.i32);
      break;
    }
    case EXSE_PREC_FP64:
    {
      for (int i = 0; i < EXT_HEADER_LEN; i++) tempheader.i64[i] = (int64_t) header[i];
      binWriteInt64(fileID, byteswap, EXT_HEADER_LEN, tempheader.i64);
      break;
    }
    default:
    {
      Error("unexpected header precision %d", rprec);
      break;
    }
  }

  binWriteF77Block(fileID, byteswap, blocklen);

  extp->datasize = (size_t) header[3];
  if (number == EXT_COMP) extp->datasize *= 2;
  blocklen = extp->datasize * (size_t) rprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch (rprec)
  {
#ifdef HAVE__FLOAT16
    case EXSE_PREC_FP16:
    {
      binWriteFlt16(fileID, extp->datasize, (_Float16 *) extp->buffer);
      break;
    }
#endif
    case EXSE_PREC_FP32:
    {
      binWriteFlt32(fileID, byteswap, extp->datasize, (float *) extp->buffer);
      break;
    }
    case EXSE_PREC_FP64:
    {
      binWriteFlt64(fileID, byteswap, extp->datasize, (double *) extp->buffer);
      break;
    }
    default:
    {
      Error("unexpected data precision %d", rprec);
      break;
    }
  }

  binWriteF77Block(fileID, byteswap, blocklen);

  return 0;
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
