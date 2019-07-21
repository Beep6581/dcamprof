/*
 * (c) Copyright 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef ELOG_H_
#define ELOG_H_

void
elog(const char *format,
     ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 1, 2)));
#else
    ;
#endif

#endif
