// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_ASSERT_H
#define YAC_ASSERT_H

// YAC PUBLIC HEADER START

void yac_abort_message(char const*, char const*, int);

#ifdef YAC_FOR_CDO
#define die(msg) abort()
#else
#define die(msg) yac_abort_message((msg), __FILE__, __LINE__)
#endif

// if YAC is built for code coverage test or for CDO
#if defined(YAC_CODE_COVERAGE_TEST) || defined(YAC_FOR_CDO)

// to avoid reduction of code coverage due to error handling code not being
// executed, redefine respective macros
#define YAC_ASSERT(exp, msg) {(void)(exp);}
#define YAC_ASSERT_F(exp, format, ...) {(void)(exp);}
#define YAC_UNREACHABLE_DEFAULT(msg) {(void)(exp);}
#define YAC_UNREACHABLE_DEFAULT_F(format, ...) {(void)(exp);}

#else // defined(YAC_CODE_COVERAGE_TEST) || defined(YAC_FOR_CDO)

#define YAC_ASSERT(exp, msg) \
  {if(!((exp))) die(msg);}

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      char msg_buffer[1024]; \
      int ret = snprintf( \
        msg_buffer, sizeof(msg_buffer), ((format)), __VA_ARGS__); \
      if ((ret >= 0) && ((size_t)ret < sizeof(msg_buffer))) \
        die(msg_buffer); \
      else \
        die("an error occured, but error message could not be generated"); \
    } \
  }

#define YAC_UNREACHABLE_DEFAULT(msg) \
    default: \
      { \
        die(msg); \
        __attribute__((fallthrough)); \
      }

#define YAC_UNREACHABLE_DEFAULT_F(format, ...) \
    default: \
      { \
        char msg_buffer[1024]; \
        int ret = snprintf( \
          msg_buffer, sizeof(msg_buffer), ((format)), __VA_ARGS__); \
        if ((ret >= 0) && ((size_t)ret < sizeof(msg_buffer))) \
          die(msg_buffer); \
        else \
          die("an error occured, but error message could not be generated"); \
        __attribute__((fallthrough)); \
      }

#endif // defined(YAC_CODE_COVERAGE_TEST)

#define YAC_UNREACHABLE(msg) \
  YAC_ASSERT(0, (msg));

#define YAC_UNREACHABLE_F(format, ...) \
  YAC_ASSERT_F(0, (format), __VA_ARGS__)

// YAC PUBLIC HEADER STOP

#endif // YAC_ASSERT_H

