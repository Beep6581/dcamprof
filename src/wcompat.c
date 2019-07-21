/*
 * (c) Copyright 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <sys/types.h>
#include <sys/stat.h>

#include <wcompat.h>

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
#define _WINSOCKAPI_ // hack to avoid winsock.h include on msys2
#include <windows.h>
#endif

bool
file_exists(const char filename[])
{
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    struct _stat buf;
    wchar_t *wfname = utf8to16(filename);
    int res = _wstat(wfname, &buf);
    free(wfname);
    return res == 0;
#else
    struct stat buf;
    return (stat(filename, &buf) == 0);
#endif
}

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
wchar_t *
utf8to16(const char utf8[])
{
    if (utf8 == NULL) {
        return NULL;
    }
    int len = MultiByteToWideChar(CP_UTF8, 0, utf8, -1, NULL, 0);
    wchar_t *utf16 = (wchar_t *)malloc(len * sizeof(utf16[0]));
    MultiByteToWideChar(CP_UTF8, 0, utf8, -1, utf16, len);
    return utf16;
}
#endif

FILE *
utf8_fopen(const char filename[],
           const char mode[])
{
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    if (filename == NULL || mode == NULL) {
        return fopen(filename, mode);
    }
    int len = MultiByteToWideChar(CP_UTF8, 0, filename, -1, NULL, 0);
    wchar_t *fname = (wchar_t *)malloc(len * sizeof(fname[0]));
    MultiByteToWideChar(CP_UTF8, 0, filename, -1, fname, len);
    len = MultiByteToWideChar(CP_UTF8, 0, mode, -1, NULL, 0);
    wchar_t wmode[len];
    MultiByteToWideChar(CP_UTF8, 0, mode, -1, wmode, len);
    FILE *stream = _wfopen(fname, wmode);
    free(fname);
    return stream;
#else
    return fopen(filename, mode);
#endif
}
