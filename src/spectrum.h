/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>

struct spectrum_band_ {
    double f;
    double e;
};

struct spectrum_t_ {
    int band_count;
    struct spectrum_band_ v[];
};

struct spectrum_t_;
typedef struct spectrum_t_ spectrum_t;

spectrum_t *
spec_multiply(const spectrum_t *a,
              const spectrum_t *b);

spectrum_t *
spec_add(const spectrum_t *a,
         const spectrum_t *b);

void
spec_scalar_multiply(spectrum_t *a,
                     double b);

double
spec_integrate(const spectrum_t *a,
               double start,
               double end);

static inline double
spec_at_idx(const spectrum_t *a,
            int idx)
{
    return a->v[idx].e;
}

static inline int
spec_band_count(const spectrum_t *a)
{
    return a->band_count;
}

static inline double
spec_band_start(const spectrum_t *a)
{
    return a->v[0].f;
}

static inline double
spec_band_end(const spectrum_t *a)
{
    return a->v[a->band_count-1].f;
}

double
spec_at(const spectrum_t *a,
        double f);

double
spec_max(const spectrum_t *a);

double
spec_mulint(const spectrum_t *a,
            const spectrum_t *b,
            double start,
            double end);

spectrum_t *
spec_invert(const spectrum_t *a);

spectrum_t *
spec_alloc(const double amplitudes[],
           int start,
           int end,
           int spacing);

spectrum_t *
spec_alloc1(const double bands[],
            const double amplitudes[],
            int band_count);

spectrum_t *
spec_copy(const spectrum_t *s);

bool
spec_equal(const spectrum_t *a,
           const spectrum_t *b);

void
spec_print(const spectrum_t *a,
           FILE *stream);

void
spec_print_gpcolor(const spectrum_t *a,
                   uint32_t gp_color,
                   FILE *stream);

void
spec_getraw(const spectrum_t *s,
            int *band_count,
            double **bands,
            double **amp);

spectrum_t *
spec_guess_lows(const spectrum_t *s,
                double lim);

#endif
