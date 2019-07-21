/*
 * (c) Copyright 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdio.h>
#include <stdarg.h>

#include <elog.h>

void
elog(const char *format,
     ...)
{
    va_list arg;
    va_start(arg, format);
    vfprintf(stderr, format, arg);
    va_end(arg);
}
