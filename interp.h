/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef INTERP_H_
#define INTERP_H_

#include <stdbool.h>
#include <inttypes.h>

static inline double
curve_y(double x,
        const double tc[],
        const int tc_len)
{
    double xm = x * (tc_len-1);
    int idx = (int)xm;
    if (idx >= tc_len-1) {
        return tc[tc_len-1];
    }
    if (idx < 0) {
        return tc[0];
    }
    double d = xm - (double)idx; // [0 .. 1]
    return (1.0 - d) * tc[idx] + d * tc[idx+1];
}

static inline double
curve_y_rng(double x,
            double range[2],
            const double tc[],
            const int tc_len)
{
    double xm = (x - range[0]) / (range[1] - range[0]) * (tc_len - 1);
    if (xm <= 0) {
        return tc[0];
    }
    int idx = (int)xm;
    if (idx >= tc_len-1) {
        return tc[tc_len-1];
    }
    double d = xm - (double)idx; // [0 .. 1]
    return (1.0 - d) * tc[idx] + d * tc[idx+1];
}

double
inverse_curve(double y,
              const double tc[],
              int tc_len);

void
cubic_spline(const double x[],
             const double y[],
             int len,
             const double out_x[],
             double out_y[],
             int out_len);

void
rounded_linear_interpolation(const double x[],
                             const double y[],
                             const int n,
                             const double out_x[],
                             double out_y[],
                             const int out_n);

void
linear_interpolation(const double x[],
                     const double y[],
                     const int n,
                     const double out_x[],
                     double out_y[],
                     const int out_n);

void
reconstruct_tonecurve(double tc[],
                      int tc_len);

void
boxblur_lut(double *clut,
            int cs,
            int radius);

void
boxblur_1d(double array[],
           int len,
           int radius);

typedef struct p2d_ {
    double x, y;
} p2d;

bool
point_inside_polygon(const double vertx[],
                     const double verty[],
                     int nvert,
                     double testx,
                     double testy);

double
point_to_line_distance(double x1,
                       double y1,
                       double x2,
                       double y2,
                       double testx,
                       double testy);

p2d *
convex_hull(const p2d pnts_[],
            int n,
            int *nvert);

static inline double
scurve(const double x) // x must be in 0..1 range
{
    if (x < 0.5) {
        return 2*x*x;
    } else {
        return 1 - 2*(1-x)*(1-x);
    }
}

static inline double
fade(double x,
     double lolim,
     double hilim,
     double out_lo,
     double out_hi,
     bool apply_scurve)
{
    if (x < lolim) {
        return out_lo;
    } else if (x < hilim) {
        double d = (x - lolim) / (hilim - lolim);
        if (apply_scurve) {
            d = scurve(d);
        }
        return out_lo * (1-d) + out_hi * d;
    } else {
        return out_hi;
    }
}

/*
  This is a maximally equidistributed combined Tausworthe generator by L'Ecuyer.

  Generates pseudorandom numbers between 0x0 - 0xFFFFFFFF.

  Don't forget to initialise the state. The seed is 1 or higher, if set to 0,
  the default seed of 1 is used.
*/
void
random_tausworthe_init(uint32_t state[3],
                       uint32_t seed);

static inline uint32_t
random_tausworthe(uint32_t state[3])
{
#define RANDOM_TAUSWORTHE(s,a,b,c,d) ((s & c) << d) ^ (((s << a) ^ s) >> b)
    state[0] = RANDOM_TAUSWORTHE(state[0], 13, 19, ((uint32_t)0xFFFFFFFEUL), 12);
    state[1] = RANDOM_TAUSWORTHE(state[1], 2, 25, ((uint32_t)0xFFFFFFF8UL), 4);
    state[2] = RANDOM_TAUSWORTHE(state[2], 3, 11, ((uint32_t)0xFFFFFFF0UL), 17);
#undef RANDOM_TAUSWORTHE
    return (state[0] ^ state[1] ^ state[2]);
}

#endif
