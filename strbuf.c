/*
 * (c) Copyright 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdarg.h>
#include <stdio.h>

#include <strbuf.h>

static void
increase_sb_capacity(strbuf_t *sb)
{
    /* we double the capacity (up to a certain limit) to avoid frequent
       reallocs() in length loops with lots of asprintf() calls */
    if (sb->capacity_ < 16 * 1024 * 1024) {
        if (sb->capacity_ == 0) {
            sb->capacity_ = 2048;
        } else {
            sb->capacity_ <<= 1;
        }
    } else {
        sb->capacity_ += 1024 * 1024;
    }
    sb->buf_ = realloc(sb->buf_, sb->capacity_);
}

int
strbuf_vasprintf(strbuf_t *sb,
                 const char *format,
                 va_list arg)
{
    if (sb->capacity_ - sb->offset_ < 1024) {
        increase_sb_capacity(sb);
    }
    va_list arg_copy;
    va_copy(arg_copy, arg);
    int space = sb->capacity_ - sb->offset_;
#ifdef WIN32
    /* later Windows has a C99 compliant vsnprintf(), older return -1 on truncation,
       we force-use that here to make most compatible code */
    int fmtlen = _vsnprintf(&sb->buf_[sb->offset_], space, format, arg);
#else
    int fmtlen = vsnprintf(&sb->buf_[sb->offset_], space, format, arg);
#endif

    if (fmtlen < 0) {
#ifdef WIN32
        /* capacity error */
        for (int i = 0; i < 10000; i++) {
            increase_sb_capacity(sb);
            va_list arg_copy1;
            va_copy(arg_copy1, arg_copy);
            space = sb->capacity_ - sb->offset_;
            fmtlen = _vsnprintf(&sb->buf_[sb->offset_], space, format, arg_copy1);
            if (fmtlen != -1) {
                strbuf_inc(sb, fmtlen);
                break;
            }
        }
#endif
        /* should not really happen */
        return fmtlen;
    }
    if (fmtlen >= space) {
        strbuf_reserve(sb, fmtlen + 4);
        vsnprintf(&sb->buf_[sb->offset_], fmtlen+1, format, arg_copy);
    }
    strbuf_inc(sb, fmtlen);
    return fmtlen;
}

int
strbuf_asprintf(strbuf_t *sb,
                const char *format,
                ...)
{
    int ret;
    va_list arg;
    va_start(arg, format);
    ret = strbuf_vasprintf(sb, format, arg);
    va_end(arg);
    return ret;
}
