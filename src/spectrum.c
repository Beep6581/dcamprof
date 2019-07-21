/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <elog.h>
#include <spectrum.h>

static double
traparea(double f1,
         double f2,
         double e1,
         double e2)
{
    assert(f2 > f1);
    if (e2 > e1) {
        return (e2 - e1) * (f2 - f1) * 0.5 + e1 * (f2 - f1);
    } else {
        return (e1 - e2) * (f2 - f1) * 0.5 + e2 * (f2 - f1);
    }
}

static double
midpoint_e(double f1,
           double mid_f,
           double f2,
           double e1,
           double e2)
{
    assert(f2 > f1);
    if (mid_f == f1) return e1;
    if (mid_f == f2) return e2;
    assert(mid_f < f2 && mid_f > f1);

    double b = e2 - e1;
    double a = f2 - f1;
    double a2 = mid_f - f1;
    return e1 + a2*(b/a);
}

spectrum_t *
spec_multiply(const spectrum_t *a,
              const spectrum_t *b)
{
    spectrum_t *c = malloc(sizeof(*c)+(a->band_count+b->band_count)*sizeof(c->v[0]));
    int i = 0, j = 0, k = 0;
    while (i < a->band_count && j < b->band_count) {
        if (a->v[i].f == b->v[j].f) {
            c->v[k].f = a->v[i].f;
            c->v[k].e = a->v[i].e * b->v[j].e;
            i++;
            j++;
            k++;
        } else if (a->v[i].f > b->v[j].f) {
            if (i == 0) {
                j++;
            } else {
                double f1 = a->v[i-1].f;
                double f2 = a->v[i].f;
                double e1 = a->v[i-1].e;
                double e2 = a->v[i].e;
                c->v[k].f = b->v[j].f;
                c->v[k].e = b->v[j].e * midpoint_e(f1, b->v[j].f, f2, e1, e2);
                k++;
                j++;
            }
        } else { // b->v[j].f > a->v[i].f
            if (j == 0) {
                i++;
            } else {
                double f1 = b->v[j-1].f;
                double f2 = b->v[j].f;
                double e1 = b->v[j-1].e;
                double e2 = b->v[j].e;
                c->v[k].f = a->v[i].f;
                c->v[k].e = a->v[i].e * midpoint_e(f1, a->v[i].f, f2, e1, e2);
                k++;
                i++;
            }
        }
    }
    c->band_count = k;
    if (k < 2) {
        c->band_count = 2;
        c->v[0].f = a->v[0].f < b->v[0].f ? a->v[0].f : b->v[0].f;
        c->v[1].f = a->v[a->band_count-1].f > b->v[b->band_count-1].f ? a->v[a->band_count-1].f : b->v[b->band_count-1].f;
        c->v[0].e = c->v[1].e = 0;
    }
    c = realloc(c, sizeof(*c)+c->band_count*sizeof(c->v[0]));
    return c;
}

spectrum_t *
spec_add(const spectrum_t *a,
         const spectrum_t *b)
{
    // FIXME: copy of spec_multiply(), merge
    spectrum_t *c = malloc(sizeof(*c)+(a->band_count+b->band_count)*sizeof(c->v[0]));
    int i = 0, j = 0, k = 0;
    while (i < a->band_count && j < b->band_count) {
        if (a->v[i].f == b->v[j].f) {
            c->v[k].f = a->v[i].f;
            c->v[k].e = a->v[i].e + b->v[j].e;
            i++;
            j++;
            k++;
        } else if (a->v[i].f > b->v[j].f) {
            if (i == 0) {
                j++;
            } else {
                double f1 = a->v[i-1].f;
                double f2 = a->v[i].f;
                double e1 = a->v[i-1].e;
                double e2 = a->v[i].e;
                c->v[k].f = b->v[j].f;
                c->v[k].e = b->v[j].e + midpoint_e(f1, b->v[j].f, f2, e1, e2);
                k++;
                j++;
            }
        } else { // b->v[j].f > a->v[i].f
            if (j == 0) {
                i++;
            } else {
                double f1 = b->v[j-1].f;
                double f2 = b->v[j].f;
                double e1 = b->v[j-1].e;
                double e2 = b->v[j].e;
                c->v[k].f = a->v[i].f;
                c->v[k].e = a->v[i].e + midpoint_e(f1, a->v[i].f, f2, e1, e2);
                k++;
                i++;
            }
        }
    }
    c->band_count = k;
    if (k < 2) {
        c->band_count = 2;
        c->v[0].f = a->v[0].f < b->v[0].f ? a->v[0].f : b->v[0].f;
        c->v[1].f = a->v[a->band_count-1].f > b->v[b->band_count-1].f ? a->v[a->band_count-1].f : b->v[b->band_count-1].f;
        c->v[0].e = c->v[1].e = 0;
    }
    c = realloc(c, sizeof(*c)+c->band_count*sizeof(c->v[0]));
    return c;
}


void
spec_scalar_multiply(spectrum_t *a,
                     double b)
{
    for (int i = 0; i < a->band_count; i++) {
        a->v[i].e *= b;
    }
}

double
spec_at(const spectrum_t *a,
        double f)
{
    if (f > a->v[a->band_count-1].f) {
        return 0;
    }
    for (int i = 0; i < a->band_count; i++) {
        if (f <= a->v[i].f) {
            if (f == a->v[i].f) {
                return a->v[i].e;
            }
            if (i == 0) {
                return 0;
            }
            return midpoint_e(a->v[i-1].f, f, a->v[i].f, a->v[i-1].e, a->v[i].e);
        }
    }
    return 0;
}

double
spec_max(const spectrum_t *a)
{
    double maxe = a->v[0].e;
    for (int i = 1; i < a->band_count; i++) {
        if (a->v[i].e > maxe) {
            maxe = a->v[i].e;
        }
    }
    return maxe;
}

double
spec_integrate(const spectrum_t *a,
               double start,
               double end)
{
    if (start < 0) {
        start = a->v[0].f;
    }
    if (end < 0) {
        end = a->v[a->band_count-1].f;
    }
    if (start >= end) {
//        elog("a %f %f %d\n", start, end, a->band_count);
        return 0.0;
    }
    // trapezoid integration
    double sum = 0;
    int i = 0;
    // skip bands before start
    while (a->v[i].f < start && i < a->band_count) {
        i++;
    }
    if (i == a->band_count) {
        return 0;
    }
    if (i > 0 && a->v[i].f > start) {
        // special case first trapeziod
        double f1 = a->v[i-1].f;
        double f2 = a->v[i].f;
        double e1 = a->v[i-1].e;
        double e2 = a->v[1].e;
        if (a->v[i].f > end) {
            e2 = midpoint_e(f1, end, f2, e1, e2);
            f2 = end;
        }
        sum += traparea(start, f2, midpoint_e(f1, start, f2, e1, e2), e2);
        i++;
    }
    if (i == 0) {
        //      elog("c\n");
        i++;
    }
    while (i < a->band_count && a->v[i].f <= end) {
        sum += traparea(a->v[i-1].f, a->v[i].f, a->v[i-1].e, a->v[i].e);
//        elog("%f %f\n", a->v[i-1].f, traparea(a->v[i-1].f, a->v[i].f, a->v[i-1].e, a->v[i].e));
        i++;
    }
    if (i == a->band_count || a->v[i-1].f == end) {
        return sum;
    }
    assert(a->v[i].f > end);

    {
        // special case last trapeziod
        double f1 = a->v[i-1].f;
        double f2 = a->v[i].f;
        double e1 = a->v[i-1].e;
        double e2 = a->v[1].e;
        sum += traparea(f1, end, e1, midpoint_e(f1, end, f2, e1, e2));
    }

    return sum;
}

double
spec_mulint(const spectrum_t *a,
            const spectrum_t *b,
            double start,
            double end)
{
    spectrum_t *c = spec_multiply(a, b);
    double res = spec_integrate(c, start, end);
    free(c);
    return res;
}

spectrum_t *
spec_alloc1(const double bands[],
            const double amplitudes[],
            int band_count)
{
    if (band_count <= 0) {
        return NULL;
    }
    spectrum_t *c = malloc(sizeof(*c)+band_count*sizeof(c->v[0]));
    c->band_count = band_count;
    for (int i = 0; i < band_count; i++) {
        if (i > 0 && bands[i] <= bands[i-1]) {
            elog("%f <= %f\n", bands[i], bands[i-1]);
            abort();
        }
        if (bands[i] < 0) {
            abort();
        }
        c->v[i].f = bands[i];
        c->v[i].e = amplitudes[i];
    }
    return c;
}

spectrum_t *
spec_invert(const spectrum_t *a)
{
    spectrum_t *c = malloc(sizeof(*c) + a->band_count * sizeof(c->v[0]));
    memcpy(c, a, sizeof(*a));
    for (int i = 0; i < a->band_count; i++) {
        c->v[i].f = a->v[i].f;
        c->v[i].e = 1.0 / a->v[i].e;
        if (!isfinite(c->v[i].e)) c->v[i].e = 0.0;
    }
    return c;
}

spectrum_t *
spec_alloc(const double amplitudes[],
           int start,
           int end,
           int spacing)
{
    int band_count = (end - start) / spacing + 1;
    assert(band_count > 0);
    assert(end > start);
    spectrum_t *c = malloc(sizeof(*c) + band_count * sizeof(c->v[0]));
    c->band_count = band_count;
    int f, i;
    for (f = start, i = 0; i < band_count; f += spacing, i++) {
        c->v[i].f = (double)f;
        c->v[i].e = amplitudes[i];
    }
    assert(f == end + spacing);
    return c;
}

spectrum_t *
spec_copy(const spectrum_t *s)
{
    if (s == NULL) {
        return NULL;
    }
    spectrum_t *c = malloc(sizeof(*c) + s->band_count * sizeof(c->v[0]));
    memcpy(c, s, sizeof(*c) + s->band_count * sizeof(c->v[0]));
    return c;
}

bool
spec_equal(const spectrum_t *a,
           const spectrum_t *b)
{
    if (a == b) return true;
    if (a->band_count != b->band_count) return false;
    for (int i = 0; i < a->band_count; i++) {
        if (a->v[i].f != b->v[i].f || a->v[i].e != b->v[i].e) return false;
    }
    return true;
}

void
spec_print(const spectrum_t *a,
           FILE *stream)
{
    for (int i = 0; i < a->band_count; i++) {
        fprintf(stream, "%.1f %f\n", a->v[i].f, a->v[i].e);
    }
}

void
spec_print_gpcolor(const spectrum_t *a,
                   uint32_t gp_color,
                   FILE *stream)
{
    for (int i = 0; i < a->band_count; i++) {
        fprintf(stream, "%.1f %f 0x%06x\n", a->v[i].f, a->v[i].e, gp_color);
    }
}

void
spec_getraw(const spectrum_t *s,
            int *band_count,
            double **bands,
            double **amp)
{
    if (band_count) *band_count = s->band_count;
    if (bands) *bands = malloc(s->band_count * sizeof(double));
    if (amp) *amp = malloc(s->band_count * sizeof(double));
    for (int i = 0; i < s->band_count; i++) {
        if (bands) (*bands)[i] = s->v[i].f;
        if (amp) (*amp)[i] = s->v[i].e;
    }
}

spectrum_t *
spec_guess_lows(const spectrum_t *s_,
                double lim)
{
    // FIXME: just stupid temp test... should do some spline or something
    spectrum_t *s = spec_copy(s_);
    double lim_a = spec_at(s, lim);
    const double base = 380;
    for (int i = 0; i < s->band_count; i++) {
        if (s->v[i].f < base) {
            s->v[i].e = 0.0;
        } else if (s->v[i].f < lim) {
            s->v[i].e = lim_a * (s->v[i].f - base) / (lim - base);
        } else {
            break;
        }
    }
    return s;
}
