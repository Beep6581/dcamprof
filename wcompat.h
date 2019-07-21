/*
 * (c) Copyright 2017 - 2018 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef WCOMPAT_H_
#define WCOMPAT_H_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#ifndef bit_swap16
#define bit_swap16(value)                                               \
    (((((uint16_t)(value)) & 0xFF00) >> 8) |                            \
     ((((uint16_t)(value)) & 0x00FF) << 8))
#endif
#ifndef bit_swap32
#define bit_swap32(value)                                               \
    (((((uint32_t)(value)) & 0xFF000000) >> 24) |                       \
     ((((uint32_t)(value)) & 0x00FF0000) >>  8) |                       \
     ((((uint32_t)(value)) & 0x0000FF00) <<  8) |                       \
     ((((uint32_t)(value)) & 0x000000FF) << 24))
#endif

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
//#include <winsock2.h> requires lib, so we defined inline functions instead
#ifndef htons
static inline uint16_t htons(uint16_t x) {
#ifndef __BYTE_ORDER__
# error "__BYTE_ORDER__ not defined by the compiler, check how your compiler defines byte order and update the code accordingly"
#endif
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return x;
#elif __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return bit_swap16(x);
#else
# error "unknown byte order, or __ORDER_BIG_ENDIAN__/__ORDER_LITTLE_ENDIAN__ not defined. Check how your complier defines byte order and update the code accordingly"
#endif
}
#define ntohs htons
#endif

#ifndef htonl
static inline uint32_t htonl(uint32_t x) {
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  return x;
#elif __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  return bit_swap32(x);
#else
# error "unknown byte order, or __ORDER_BIG_ENDIAN__/__ORDER_LITTLE_ENDIAN__ not defined. Check how your complier defines byte order and update the code accordingly"
#endif
}
#define ntohl htonl
#endif

#else
#include <arpa/inet.h>
#endif

bool
file_exists(const char filename[]);

FILE *
utf8_fopen(const char filename[],
           const char mode[]);

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
wchar_t *
utf8to16(const char utf8[]);
#endif

#endif
