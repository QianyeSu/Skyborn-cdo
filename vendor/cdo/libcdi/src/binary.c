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

#define CDI_BIGENDIAN 0     // Byte order BIGENDIAN
#define CDI_LITTLEENDIAN 1  // Byte order LITTLEENDIAN
#include "error.h"
#include "file.h"
#include "swap.h"
#include "binary.h"

uint32_t
get_uint32(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((uint32_t)x[0])<<24) + (((uint32_t)x[1])<<16) + (((uint32_t)x[2])<< 8) + (uint32_t)x[3];
    case CDI_LITTLEENDIAN:
      return (((uint32_t)x[3])<<24) + (((uint32_t)x[2])<<16) + (((uint32_t)x[1])<< 8) + (uint32_t)x[0];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
  // clang-format on
}

uint32_t
get_swap_uint32(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((uint32_t)x[3])<<24) + (((uint32_t)x[2])<<16) + (((uint32_t)x[1])<< 8) + (uint32_t)x[0];
    case CDI_LITTLEENDIAN:
      return (((uint32_t)x[0])<<24) + (((uint32_t)x[1])<<16) + (((uint32_t)x[2])<< 8) + (uint32_t)x[3];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
  // clang-format on
}

uint64_t
get_uint64(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((uint64_t)x[0])<<56) + (((uint64_t)x[1])<<48) + (((uint64_t)x[2])<<40) + (((uint64_t)x[3])<<32) +
             (((uint64_t)x[4])<<24) + (((uint64_t)x[5])<<16) + (((uint64_t)x[6])<< 8) +   (uint64_t)x[7];
    case CDI_LITTLEENDIAN:
      return (((uint64_t)x[7])<<56) + (((uint64_t)x[6])<<48) + (((uint64_t)x[5])<<40) + (((uint64_t)x[4])<<32) +
             (((uint64_t)x[3])<<24) + (((uint64_t)x[2])<<16) + (((uint64_t)x[1])<< 8) +   (uint64_t)x[0];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
  // clang-format on
}

uint64_t
get_swap_uint64(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((uint64_t)x[7])<<56) + (((uint64_t)x[6])<<48) + (((uint64_t)x[5])<<40) + (((uint64_t)x[4])<<32) +
             (((uint64_t)x[3])<<24) + (((uint64_t)x[2])<<16) + (((uint64_t)x[1])<< 8) +   (uint64_t)x[0];
    case CDI_LITTLEENDIAN:
      return (((uint64_t)x[0])<<56) + (((uint64_t)x[1])<<48) + (((uint64_t)x[2])<<40) + (((uint64_t)x[3])<<32) +
             (((uint64_t)x[4])<<24) + (((uint64_t)x[5])<<16) + (((uint64_t)x[6])<< 8) +   (uint64_t)x[7];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
  // clang-format on
}

size_t
binReadF77Block(int fileID, int byteswap)
{
  unsigned char f77block[4];
  size_t blocklen = 0;

  if (fileRead(fileID, f77block, 4) == 4) { blocklen = byteswap ? get_swap_uint32(f77block) : get_uint32(f77block); }

  return blocklen;
}

void
binWriteF77Block(int fileID, int byteswap, size_t blocksize)
{
  static const unsigned int s[4] = { 0, 8, 16, 24 };
  const unsigned long ublocksize = (unsigned long) blocksize;
  unsigned char f77block[4];

  switch (HOST_ENDIANNESS)
  {
    case CDI_BIGENDIAN:
      if (byteswap)
      {
        for (int i = 0; i <= 3; ++i) f77block[i] = (unsigned char) (ublocksize >> s[i]);
      }
      else
      {
        for (int i = 0; i <= 3; ++i) f77block[3 - i] = (unsigned char) (ublocksize >> s[i]);
      }
      break;
    case CDI_LITTLEENDIAN:
      if (byteswap)
      {
        for (int i = 0; i <= 3; ++i) f77block[3 - i] = (unsigned char) (ublocksize >> s[i]);
      }
      else
      {
        for (int i = 0; i <= 3; ++i) f77block[i] = (unsigned char) (ublocksize >> s[i]);
      }
      break;
    default: Error("Unhandled endianness %d", HOST_ENDIANNESS);
  }

  if (fileWrite(fileID, f77block, 4) != 4) Error("Write failed on %s", fileInqName(fileID));
}

int
binReadInt32(int fileID, int byteswap, size_t size, int32_t *ptr)
{
  fileRead(fileID, (void *) ptr, 4 * size);
  if (byteswap) swap4byte(ptr, size);
  return 0;
}

int
binReadInt64(int fileID, int byteswap, size_t size, int64_t *ptr)
{
  fileRead(fileID, (void *) ptr, 8 * size);
  if (byteswap) swap8byte(ptr, size);
  return 0;
}

int
binWriteInt32(int fileID, int byteswap, size_t size, int32_t *ptr)
{
  if (byteswap) swap4byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 4 * size);
  return 0;
}

int
binWriteInt64(int fileID, int byteswap, size_t size, int64_t *ptr)
{
  if (byteswap) swap8byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 8 * size);
  return 0;
}

int
binWriteFlt32(int fileID, int byteswap, size_t size, float *ptr)
{
  if (byteswap) swap4byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 4 * size);
  return 0;
}

int
binWriteFlt64(int fileID, int byteswap, size_t size, double *ptr)
{
  if (byteswap) swap8byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 8 * size);
  return 0;
}

#ifdef HAVE__FLOAT16
int
binWriteFlt16(int fileID, size_t size, _Float16 *ptr)
{
  fileWrite(fileID, (void *) ptr, 2 * size);
  return 0;
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
