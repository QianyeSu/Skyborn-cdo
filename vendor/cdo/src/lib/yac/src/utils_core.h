// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef UTILS_CORE_H
#define UTILS_CORE_H

#include <string.h>
#include <stdint.h>

#include "utils_common.h"

void yac_quicksort_index ( int * a, size_t n, int * idx);
void yac_quicksort_index_yac_int_size_t ( yac_int * a, size_t n, size_t * idx);
void yac_quicksort_index_yac_int_yac_int ( yac_int * a, size_t n, yac_int * idx);
void yac_quicksort_index_size_t_yac_int ( size_t * a, size_t n, yac_int * idx);
void yac_quicksort_index_yac_int_uint64_t ( yac_int * a, size_t n, uint64_t * idx);
void yac_quicksort_index_yac_int_int ( yac_int * a, size_t n, int * idx);
void yac_quicksort_index_size_t_size_t ( size_t * a, size_t n, size_t * idx);
void yac_quicksort_index_uint64_t_size_t ( uint64_t * a, size_t n, size_t * idx);
void yac_quicksort_index_int_size_t ( int * a, size_t n, size_t * idx);
void yac_quicksort_index_size_t_int ( size_t * a, size_t n, int * idx);
void yac_quicksort_index_size_t_void_p ( size_t * a, size_t n, void * * idx);
void yac_quicksort_index_int_yac_int ( int * a, size_t n, yac_int * idx);
void yac_quicksort_index_int_double ( int * a, size_t n, double * idx);
void yac_quicksort_index_size_t_size_t_double (
  size_t * a, size_t n, size_t * b, double * c );
void yac_quicksort_index_yac_int_yac_int_double (
  yac_int * a, size_t n, yac_int * b, double * c );
void yac_quicksort_index_yac_int_yac_int_size_t (
  yac_int * a, size_t n, yac_int * b, size_t * c );
void yac_quicksort_index_int_size_t_size_t (
  int * a, size_t n, size_t * b, size_t * c );
void yac_quicksort_index_int_size_t_yac_int (
  int * a, size_t n, size_t * b, yac_int * c );

void yac_mergesort(void* base, size_t num, size_t size,
                   int (*compar)(const void*,const void*));

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_uint(unsigned * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   unsigned prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of doubles
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_double(double * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   double prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      // use memcmp to handle NaNs correctly
      if (!memcmp(&array[i], &prev, sizeof(prev))) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of size_t
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t(size_t * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of size_t[2]
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t_2(
  size_t (*array)[2], size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev[2] = {array[0][0],
                     array[0][1]};

   for (size_t i = 1; i < N; ++i) {

      if ((array[i][0] == prev[0]) &&
          (array[i][1] == prev[1]))continue;

      prev[0] = array[i][0];
      prev[1] = array[i][1];
      ++pos;

      if (pos != i) {
        array[pos][0] = array[i][0];
        array[pos][1] = array[i][1];
      }
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of size_t[3]
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t_3(
  size_t (*array)[3], size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev[3] = {array[0][0],
                     array[0][1],
                     array[0][2]};

   for (size_t i = 1; i < N; ++i) {

      if ((array[i][0] == prev[0]) &&
          (array[i][1] == prev[1]) &&
          (array[i][2] == prev[2]))continue;

      prev[0] = array[i][0];
      prev[1] = array[i][1];
      prev[2] = array[i][2];
      ++pos;

      if (pos != i) {
        array[pos][0] = array[i][0];
        array[pos][1] = array[i][1];
        array[pos][2] = array[i][2];
      }
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of yac_int
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_yac_int(
   yac_int * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   yac_int prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

// Sorts the provided arrays based on a flag-array (containing the
// values "0" and "!= 0"). After the sort, all array elements whose associated
// flag value is "0" are the front of the array.
//
// This sort is:
//   * not stable
//   * has a time complexity of O(n)
static inline void yac_flag_sort_size_t(
  size_t * array_size_t, int * flag, size_t false_count) {

  // The number of "true" elements in the 0...false_count-1 range of the
  // array is identical to the number of "false" elements in the
  // false_count...false_count+true_count-1. We just have to find matching
  // pairs and swap them.
  for (size_t i = 0, j = false_count; i < false_count; ++i) {
    // if there is a wrongfully placed "true" element
    if (flag[i]) {
      // find a wrongfully place "false" element
      for (; flag[j]; ++j);
      // swap elements
      size_t temp_size_t = array_size_t[i];
      array_size_t[i] = array_size_t[j];
      array_size_t[j] = temp_size_t;
      // set to next element in "true" list
      ++j;
    }
  }
}

#define ASSERT(c) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s in %s:%d\n",\
           #c, __FILE__, __LINE__);\
   abort ();\
}

#define COPY_DATA(data, count) \
  (memcpy( \
    xmalloc((size_t)(count) * sizeof(*(data))), \
    (data), (size_t)(count) * sizeof(*(data))))

#endif // UTILS_CORE_H

