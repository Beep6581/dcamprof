/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

#include <bisection.h>
#include <interp.h>

void
cubic_spline(const double x[],
             const double y[],
             const int len,
             const double out_x[],
             double out_y[],
             const int out_len)
{
    int i, j;

    double **A = malloc(2 * len * sizeof(*A));
    double *As = calloc(1, 2 * len * 2 * len * sizeof(*As));
    double *b = calloc(1, 2*len*sizeof(*b));
    double *c = calloc(1, 2*len*sizeof(*c));
    double *d = calloc(1, 2*len*sizeof(*d));

    for (i = 0; i < 2*len; i++) {
        A[i] = &As[2*len*i];
    }

    for (i = len-1; i > 0; i--) {
        b[i] = (y[i] - y[i-1]) / (x[i] - x[i-1]);
        d[i-1] = x[i] - x[i-1];
    }

    for (i = 1; i < len-1; i++) {
        A[i][i] = 2 * (d[i-1] + d[i]);
        if (i > 1) {
            A[i][i-1] = d[i-1];
            A[i-1][i] = d[i-1];
        }
        A[i][len-1] = 6 * (b[i+1] - b[i]);
    }

    for(i = 1; i < len-2; i++) {
        double v = A[i+1][i] / A[i][i];
        for(j = 1; j <= len-1; j++) {
            A[i+1][j] -= v * A[i][j];
        }
    }

    for(i = len-2; i > 0; i--) {
        double acc = 0;
        for(j = i; j <= len-2; j++) {
            acc += A[i][j]*c[j];
        }
        c[i] = (A[i][len-1] - acc) / A[i][i];
    }

    for (i = 0; i < out_len; i++) {
        double x_out = out_x == NULL ? (double)i / (out_len-1) : out_x[i];
        double y_out = 0;
        for (j = 0; j < len-1; j++) {
            if (x[j] <= x_out && x_out <= x[j+1]) {
                double v = x_out - x[j];
                y_out = y[j] +
                    ((y[j+1] - y[j]) / d[j] - (2 * d[j] * c[j] + c[j+1] * d[j]) / 6) * v +
                    (c[j] * 0.5) * v*v +
                    ((c[j+1] - c[j]) / (6 * d[j])) * v*v*v;
            }
        }
        out_y[i] = y_out;
    }
    free(A);
    free(As);
    free(b);
    free(c);
    free(d);
}

void
rounded_linear_interpolation(const double x[],
                             const double y[],
                             const int n,
                             const double out_x[],
                             double out_y[],
                             const int out_n)
{
    for (int k = 0; k < out_n; k++) {
        double p = out_x == NULL ? ((double)k / (out_n - 1)) : out_x[k];
        for (int i = 0; i < n-1; i++) {
            if (x[i] <= p && p <= x[i+1]) {
                double d = (p - x[i]) / (x[i+1] - x[i]);
                d = scurve(d);
                out_y[k] = y[i]*(1-d) + y[i+1]*d;
                break;
            }
        }
    }
}

void
linear_interpolation(const double x[],
                     const double y[],
                     const int n,
                     const double out_x[],
                     double out_y[],
                     const int out_n)
{
    for (int k = 0; k < out_n; k++) {
        double p = out_x == NULL ? ((double)k / (out_n - 1)) : out_x[k];
        for (int i = 0; i < n-1; i++) {
            if (x[i] <= p && p <= x[i+1]) {
                double d = (p - x[i]) / (x[i+1] - x[i]);
                out_y[k] = y[i]*(1-d) + y[i+1]*d;
                break;
            }
        }
    }
}

struct inverse_curve_bisection_fun_arg {
    const double *tc;
    int tc_len;
    double y;
};

static double
inverse_curve_bisection_fun(double x,
                            void *arg)
{
    struct inverse_curve_bisection_fun_arg *a = (struct inverse_curve_bisection_fun_arg *)arg;
    return curve_y(x, a->tc, a->tc_len) - a->y;
}

double
inverse_curve(double y,
              const double tc[],
              int tc_len)
{
    // handle diagonal curve edge cases (had problems with bisection() getting to 0, not checked why, anyway this special case fixes it)
    if (y == 0 && tc[0] == 0) {
        return 0;
    }
    if (y == 1 && tc[tc_len-1] == 1) {
        return 1;
    }

    struct inverse_curve_bisection_fun_arg arg;
    arg.y = y;
    arg.tc = tc;
    arg.tc_len = tc_len;
    return bisection(inverse_curve_bisection_fun, &arg, 0, 1, 1.0e-7, 10000);
}

void
reconstruct_tonecurve(double tc[],
                      int tc_len)
{
    double x0 = 0;
    double y0 = tc[0];
    int j = 0;
    for (int i = 1; i < tc_len; i++) {
        double x1 = (double)i / (tc_len - 1);
        double y1 = tc[i];
        if (y1 != y0) {
            for (int k = j + 1; k < i; k++) {
                double x = (double)k / (tc_len - 1);
                double d = (x - x0) / (x1 - x0);
                tc[k] = y0 + (y1 - y0) * d;
            }
            y0 = y1;
            x0 = x1;
            j = i;
        }
    }
}

void
boxblur_lut(double *clut,
            int cs,
            int dist)
{
    double *clut1 = malloc(sizeof(clut[0]) * cs * cs * cs);
#pragma omp parallel for
    for (int r_ = 0; r_ < cs; r_++) {
        for (int g_ = 0; g_ < cs; g_++) {
            for (int b_ = 0; b_ < cs; b_++) {
                double sum = 0;
                int acc = 0;
                int rs = r_-dist;
                if (rs < 0) rs = 0;
                int re = r_+dist;
                if (re > cs-1) re = cs-1;
                for (int r = rs; r <= re; r++) {
                    int gs = g_-dist;
                    if (gs < 0) gs = 0;
                    int ge = g_+dist;
                    if (ge > cs-1) ge = cs-1;
                    for (int g = gs; g <= ge; g++) {
                        int bs = b_-dist;
                        if (bs < 0) bs = 0;
                        int be = b_+dist;
                        if (be > cs-1) be = cs-1;
                        for (int b = bs; b <= be; b++) {
                            sum += clut[cs*cs*r+cs*g+b];
                            acc++;
                        }
                    }
                }
                clut1[cs*cs*r_+cs*g_+b_] = sum / acc;
            }
        }
    }
    memcpy(clut, clut1, sizeof(clut[0]) * cs * cs * cs);
    free(clut1);
}

void
boxblur_1d(double in[],
           int length,
           int radius)
{
    int w = length;
    int len = radius + 1;
    double *buf = malloc(w * sizeof(buf[0]));
    buf[0] = in[0] / len;
    for (int j = 1; j <= radius; j++) {
        buf[0] += in[j] / len;
    }
    for (int x = 1; x <= radius; x++) {
        buf[x] = (buf[x-1]*len + in[x+radius]) / (len+1);
        len++;
    }
    for (int x = radius + 1; x < w - radius; x++) {
        buf[x] = buf[x-1] +(in[x+radius] - in[x-radius-1]) / len;
    }
    for (int x = w - radius; x < w; x++) {
        buf[x] = (buf[x-1]*len - in[x-radius-1]) / (len-1);
        len--;
    }
    memcpy(in, buf, length * sizeof(in[0]));
    free(buf);
}

// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
bool
point_inside_polygon(const double vertx[],
                     const double verty[],
                     int nvert,
                     double testx,
                     double testy)
{
    int i, j;
    bool c = false;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
             (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
            c = !c;
    }
    return c;
}

static int
p2d_cmp(const void *a_,
        const void *b_)
{
    const p2d *a = (const p2d *)a_;
    const p2d *b = (const p2d *)b_;
    if (a->x < b->x) return -1;
    if (b->x < a->x) return 1;
    if (a->y < b->y) return -1;
    if (b->y < a->y) return 1;
    return 0;
}

static void
sort_p2d(p2d p[],
         int count)
{
    qsort(p, count, sizeof(p[0]), p2d_cmp);
}

/* Three points are a counter-clockwise turn if ccw > 0, clockwise if
 * ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 * gives the signed area of the triangle formed by p1, p2 and p3.
 */
static inline double
hull_ccw(p2d p1, p2d p2, p2d p3)
{
  return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
p2d *
convex_hull(const p2d points_[],
            int npoints,
            int *out_hullsize)
{
    p2d *points = malloc(npoints * sizeof(points[0]));
    memcpy(points, points_, npoints * sizeof(points[0]));
    sort_p2d(points, npoints); // lexical order
    int j = 0;
    // remove duplicates
    for (int i = 1; i < npoints; i++) {
        if (points[i].x == points[j].x && points[i].y == points[j].y) {
            continue;
        }
        j++;
        points[j] = points[i];
    }
    npoints = j+1;
    if (npoints < 3) {
        free(points);
        *out_hullsize = 0;
        return NULL;
    }

    p2d *hull = malloc((npoints + 1) * sizeof(*hull));
    ssize_t i, t, k = 0;

    // lower hull
    for (i = 0; i < npoints; ++i) {
        while (k >= 2 && hull_ccw(hull[k-2], hull[k-1], points[i]) <= 0) --k;
        hull[k++] = points[i];
    }

    // upper hull
    for (i = npoints-2, t = k+1; i >= 0; --i) {
        while (k >= t && hull_ccw(hull[k-2], hull[k-1], points[i]) <= 0) --k;
        hull[k++] = points[i];
    }
    *out_hullsize = k-1; // last point is same as first, so we remove that

    hull = realloc(hull, k * sizeof(hull[0]));
    free(points);
    return hull;
}

static double
distance_2d(const double xy1[2],
            const double xy2[2])
{
    return sqrt((xy1[0]-xy2[0])*(xy1[0]-xy2[0]) + (xy1[1]-xy2[1])*(xy1[1]-xy2[1]));
}

static double
dot_product_2d(const double xy1[2],
               const double xy2[2])
{
    return xy1[0]*xy2[0] + xy1[1]*xy2[1];
}

double
point_to_line_distance(double x1,
                       double y1,
                       double x2,
                       double y2,
                       double testx,
                       double testy)
{
    double p[2] = { testx, testy };
    double v[2] = { x1, y1 };
    double w[2] = { x2, y2 };
    // Return minimum distance between line segment vw and point p
    //const double l2 = length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
    const double l2 = (w[0]-v[0])*(w[0]-v[0]) + (w[1]-v[1])*(w[1]-v[1]);
    if (l2 == 0.0) {
        // v == w case
        //return distance(p, v);
        return distance_2d(p, v);
    }
    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    //const double t = dot_product_2d(p - v, w - v) / l2;
    const double p_minus_v[2] = { p[0] - v[0], p[1] - v[1] };
    const double w_minus_v[2] = { w[0] - v[0], w[1] - v[1] };
    const double t = dot_product_2d(p_minus_v, w_minus_v) / l2;
    if (t < 0.0) { // Beyond the 'v' end of the segment
        return distance_2d(p, v);
    } else if (t > 1.0) { // Beyond the 'w' end of the segment
        return distance_2d(p, w);
    }
    //const vec2 projection = v + t * (w - v);  // Projection falls on the segment
    double projection[2] = { v[0] + t * (w[0] - v[0]), v[1] + t * (w[1] - v[1]) };
    return distance_2d(p, projection);
}

void
random_tausworthe_init(uint32_t state[3],
                       uint32_t seed)
{
    /* default seed is 1 */
    if (seed == 0) {
        seed = 1;
    }

#define LCG(n) ((69069 * n) & 0xFFFFFFFFL)
    state[0] = LCG(seed);
    state[1] = LCG(state[0]);
    state[2] = LCG(state[1]);
#undef LCG
    /* "warm it up" */
    random_tausworthe(state);
    random_tausworthe(state);
    random_tausworthe(state);
    random_tausworthe(state);
    random_tausworthe(state);
    random_tausworthe(state);
}
