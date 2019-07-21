/*
 * (c) Copyright 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef STRBUF_H_
#define STRBUF_H_

#include <inttypes.h>
#include <string.h>
#include <strings.h>
#ifndef _WIN32
#include <alloca.h>
#endif
#include <stdlib.h>
#include <stdarg.h>

#ifdef strdupa
#undef strdupa
#endif

#define calloca(nmemb, size)                                                   \
({                                                                             \
    size_t size_ = (nmemb) * (size);                                           \
    void *buf_ = alloca(size_);                                                \
    memset(buf_, 0, size_);                                                    \
})

#define strdupa(str)                                                           \
({                                                                             \
    char *s_;                                                                  \
    int n_;                                                                    \
    n_ = strlen(str) + 1;                                                      \
    s_ = alloca(n_);                                                           \
    (char *)memcpy(s_, str, n_);                                               \
})

#ifdef strndupa
#undef strndupa
#endif

#define strndupa(str, maxlen)                                                  \
({                                                                             \
    char *s_;                                                                  \
    int n_;                                                                    \
    n_ = strnlen(str, maxlen);                                                 \
    s_ = alloca(n_ + 1);                                                       \
    s_[n_] = '\0';                                                             \
    (char *)memcpy(s_, str, n_);                                               \
})

// string buffer api

struct strbuf_t_ {
    char *buf_;
    size_t capacity_;
    size_t offset_;
};
typedef struct strbuf_t_ strbuf_t;

#define STRBUF_INITIALIZER ((strbuf_t){ .buf_ = NULL, .capacity_ = 0, .offset_ = 0 })

/*
  strbuf_reserve: reserve capacity from current offset. Len may exclude
  trailing '\0', that is actually there will be space for len + 1 bytes.
*/
static inline void
strbuf_reserve(strbuf_t *sb,
               size_t len)
{
    if (sb->capacity_ - sb->offset_ > len) {
        return;
    }
    sb->capacity_ += len + 256;
    sb->buf_ = realloc(sb->buf_, sb->capacity_);
}

#define strbuf_init(sb) do { memset((sb), 0, sizeof(*(sb))); } while (0)
#define strbuf_freebuf(sb) do { free((sb)->buf_); } while (0)
#define strbuf_reset(sb) do { (sb)->offset_ = 0; } while (0)

#define strbuf_unused(sb) ((sb)->capacity_ - ((sb)->offset_) - 1)
#define strbuf_len(sb) ((sb)->offset_)
#define strbuf_base(sb) ((sb)->buf_)
#define strbuf_cur(sb) (&(sb)->buf_[(sb)->offset_])
#define strbuf_inc(sb, len)                  \
do {                                         \
    (sb)->offset_ += (len);                  \
    if ((sb)->offset_ >= (sb)->capacity_) {  \
        abort();                             \
    }                                        \
} while (0)
#define strbuf_dec(sb, len)                  \
do {                                         \
    size_t len_ = (len);                     \
    if (len_ > (sb)->offset_) {              \
        abort();                             \
    }                                        \
    (sb)->offset_ -= (len);                  \
} while (0)

static inline void
strbuf_cat(strbuf_t *sb,
           const char str[])
{
    size_t len = strlen(str);
    strbuf_reserve(sb, len);
    memcpy(strbuf_cur(sb), str, len + 1);
    strbuf_inc(sb, len);
}

static inline void
strbuf_catl(strbuf_t *sb,
            const char str[],
            size_t len)
{
    strbuf_reserve(sb, len);
    memcpy(strbuf_cur(sb), str, len + 1);
    strbuf_inc(sb, len);
}

static inline void
strbuf_catb(strbuf_t *sb,
            const char chars[],
            size_t len)
{
    strbuf_reserve(sb, len);
    memcpy(strbuf_cur(sb), chars, len);
    strbuf_inc(sb, len);
    strbuf_cur(sb)[0] = '\0';
}

static inline void
strbuf_catc(strbuf_t *sb,
            const char character)
{
    strbuf_reserve(sb, 16);
    strbuf_cur(sb)[0] = character;
    strbuf_inc(sb, 1);
    strbuf_cur(sb)[0] = '\0';
}

static inline void
strbuf_catr(strbuf_t *sb,
            const char *first,
            const char *last)
{
    size_t len = (uintptr_t)last - (uintptr_t)first + 1;
    strbuf_reserve(sb, len);
    memcpy(strbuf_cur(sb), first, len);
    strbuf_inc(sb, len);
    strbuf_cur(sb)[0] = '\0';
}

// space does not need to be pre-reserved
int
strbuf_asprintf(strbuf_t *sb,
                const char *format,
                ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 2, 3)));
#else
    ;
#endif

int
strbuf_vasprintf(strbuf_t *sb,
                 const char *format,
                 va_list arg);


// WARNING: space must be pre-reserved
#define strbuf_sprintf(sb, ...)                                         \
({                                                                      \
    strbuf_t *sb_ = (sb);                                               \
    int sprintf_ret_ = sprintf(&sb_->buf_[sb_->offset_], __VA_ARGS__);  \
    strbuf_inc(sb_, sprintf_ret_);                                      \
    sprintf_ret_;                                                       \
})

// WARNING: space must be pre-reserved
#define strbuf_snprintf(sb, size, ...)                                         \
({                                                                             \
    strbuf_t *sb_ = (sb);                                                      \
    int sprintf_ret_ = snprintf(&sb_->buf_[sb_->offset_], (size), __VA_ARGS__);\
    strbuf_inc(sb_, sprintf_ret_);                                             \
    sprintf_ret_;                                                              \
})

#endif
