/*
 * (c) Copyright 2015, 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef MATMATH_H_
#define MATMATH_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <elog.h>

typedef struct v3_ {
    double v[3];
} v3;

typedef struct m3x3_ {
    double v[3][3];
} m3x3;

typedef struct matrix_ {
    int rows;
    int cols;
    double v[];
} matrix;

static inline void
v3_print(const char prefix[], const v3 v)
{
    elog("%s %f %f %f\n", prefix, v.v[0],v.v[1],v.v[2]);
}

static inline void
m3x3_print(const char prefix[], const m3x3 m)
{
    elog("%s %+.4f %+.4f %+.4f\n"
         "%s %+.4f %+.4f %+.4f\n"
         "%s %+.4f %+.4f %+.4f\n\n",
         prefix, m.v[0][0], m.v[0][1], m.v[0][2],
         prefix, m.v[1][0], m.v[1][1], m.v[1][2],
         prefix, m.v[2][0], m.v[2][1], m.v[2][2]);
}

static inline double
v3_max(const v3 v)
{
    double m = v.v[0];
    if (v.v[1] > m) m = v.v[1];
    if (v.v[2] > m) m = v.v[2];
    return m;
}

static inline double
v3_min(const v3 v)
{
    double m = v.v[0];
    if (v.v[1] < m) m = v.v[1];
    if (v.v[2] < m) m = v.v[2];
    return m;
}

static inline double
v3_mid(const v3 v)
{
    double mx = v3_max(v);
    double mn = v3_min(v);
    if (v.v[0] > mn && v.v[0] < mx) return v.v[0];
    if (v.v[1] > mn && v.v[1] < mx) return v.v[2];
    return v.v[2];
}

static inline void
v3_minmidmax(const v3 v, int *mini, int *midi, int *maxi)
{
    *mini = 0;
    *maxi = 2;
    if (v.v[1] < v.v[0]) *mini = 1;
    if (v.v[2] < v.v[*mini]) *mini = 2;
    if (v.v[0] > v.v[2]) *maxi = 0;
    if (v.v[1] > v.v[*maxi]) *maxi = 1;
    if (*mini == 0) {
        if (*maxi == 2) *midi = 1;
        else *midi = 2;
    } else if (*mini == 1) {
        if (*maxi == 2) *midi = 0;
        else *midi = 2;
    } else {
        if (*maxi == 1) *midi = 0;
        else *midi = 1;
    }
}

static inline v3
v3_clip01(v3 v)
{
    if (v.v[0] < 0) v.v[0] = 0; else if (v.v[0] > 1) v.v[0] = 1;
    if (v.v[1] < 0) v.v[1] = 0; else if (v.v[1] > 1) v.v[1] = 1;
    if (v.v[2] < 0) v.v[2] = 0; else if (v.v[2] > 1) v.v[2] = 1;
    return v;
}

static inline v3
v3_clip0(v3 v)
{
    if (v.v[0] < 0) v.v[0] = 0;
    if (v.v[1] < 0) v.v[1] = 0;
    if (v.v[2] < 0) v.v[2] = 0;
    return v;
}

static inline v3
v3_clip1(v3 v)
{
    if (v.v[0] > 1) v.v[0] = 1;
    if (v.v[1] > 1) v.v[1] = 1;
    if (v.v[2] > 1) v.v[2] = 1;
    return v;
}

static inline double
v3_avg(const v3 v)
{
    return (v.v[0] + v.v[1] + v.v[2]) / 3.0;
}

static inline bool
v3_isequal(const v3 a,
           const v3 b)
{
    return (a.v[0] == b.v[0] && a.v[1] == b.v[1] && a.v[2] == b.v[2]);
}

static inline bool
v3_isfinite(const v3 v)
{
    return (isfinite(v.v[0]) && isfinite(v.v[1]) && isfinite(v.v[2]));
}

static inline bool
v3_iszero(const v3 v)
{
    return (v.v[0] == 0 && v.v[1] == 0 && v.v[2] == 0);
}

static inline bool
v3_isclip01(const v3 v)
{
    return (v.v[0] < 0 || v.v[0] > 1 || v.v[1] < 0 || v.v[1] > 1 || v.v[2] < 0 || v.v[2] > 1);
}

static inline bool
m3x3_isequal(const m3x3 A,
              const m3x3 B)
{
    return memcmp(&A, &B, sizeof(A)) == 0;
}

static inline m3x3
m3x3_mul_m3x3(const m3x3 A,
              const m3x3 B)
{
    m3x3 M;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M.v[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                M.v[i][j] += A.v[i][k] * B.v[k][j];
            }
        }
    }
    return M;
}

static inline m3x3
m3x3_add_m3x3(const m3x3 A,
              const m3x3 B)
{
    m3x3 M;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M.v[i][j] = A.v[i][j] + B.v[i][j];
        }
    }
    return M;
}

static inline v3
m3x3_mul_v3(const m3x3 A,
            const v3 B)
{

    v3 M = { .v = { 0, 0, 0 } };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M.v[i] += A.v[i][j] * B.v[j];
        }
    }
    return M;
}

static inline m3x3
k_mul_m3x3(const double k,
           const m3x3 A)
{

    m3x3 M = A;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M.v[i][j] *= k;
        }
    }
    return M;
}

static inline v3
k_mul_v3(const double k,
         const v3 v)
{
    return (v3){{ k*v.v[0], k*v.v[1], k*v.v[2] }};
}

static inline v3
v3_norm(const v3 v)
{
    return k_mul_v3(1.0 / v3_max(v), v);
}

static inline v3
v3_add(v3 v, const v3 v2)
{
    v.v[0] += v2.v[0];
    v.v[1] += v2.v[1];
    v.v[2] += v2.v[2];
    return v;
}

static inline v3
v3_sub(v3 v, const v3 v2)
{
    v.v[0] -= v2.v[0];
    v.v[1] -= v2.v[1];
    v.v[2] -= v2.v[2];
    return v;
}

static inline v3
v3_abs(v3 v)
{
    return (v3){{ fabs(v.v[0]), fabs(v.v[1]), fabs(v.v[2]) }};
}

static inline v3
v3_invert(const v3 v)
{
    return (v3){{ 1.0/v.v[0], 1.0/v.v[1], 1.0/v.v[2] }};
}

static inline m3x3
v3_asdiagonal(const v3 v)
{
    return (m3x3){ .v = { { v.v[0], 0, 0 }, { 0, v.v[1], 0 }, { 0, 0, v.v[2] } } };
}

static inline m3x3
m3x3_invert(const m3x3 A)
{
    const double temp[3][3] = {
        { A.v[1][1] * A.v[2][2] - A.v[2][1] * A.v[1][2],
          A.v[2][1] * A.v[0][2] - A.v[0][1] * A.v[2][2],
          A.v[0][1] * A.v[1][2] - A.v[1][1] * A.v[0][2] },
        { A.v[2][0] * A.v[1][2] - A.v[1][0] * A.v[2][2],
          A.v[0][0] * A.v[2][2] - A.v[2][0] * A.v[0][2],
          A.v[1][0] * A.v[0][2] - A.v[0][0] * A.v[1][2] },
        { A.v[1][0] * A.v[2][1] - A.v[2][0] * A.v[1][1],
          A.v[2][0] * A.v[0][1] - A.v[0][0] * A.v[2][1],
          A.v[0][0] * A.v[1][1] - A.v[1][0] * A.v[0][1] }
    };
    const double det = A.v[0][0] * temp[0][0] + A.v[0][1] * temp[1][0] + A.v[0][2] * temp[2][0];
    return (m3x3){ .v = {
            { temp[0][0] / det, temp[0][1] / det, temp[0][2] / det },
            { temp[1][0] / det, temp[1][1] / det, temp[1][2] / det },
            { temp[2][0] / det, temp[2][1] / det, temp[2][2] / det }
        }};
}

static inline double *
rc(matrix *m,
   const int row,
   const int col)
{
    return &m->v[row*m->cols+col];
}

static inline matrix *
alloc_matrix(const int rows,
             const int cols)
{
    matrix *m = (matrix *)calloc(1, sizeof(*m) + rows * cols * sizeof(m->v[0]));
    m->rows = rows;
    m->cols = cols;
    return m;
}

static inline matrix *
matrix_mul(const matrix *a,
           const matrix *b)
{
    matrix *c = alloc_matrix(a->rows, b->cols);
    for (int i = 0; i < c->rows; i++) {
        for (int j = 0; j < c->cols; j++) {
            double *e = rc(c, i, j);
            *e = 0.0;
            for (int k = 0; k < a->cols; k++) {
                *e += a->v[i*a->cols+k] * b->v[k*b->cols+j];
            }
        }
    }
    return c;
}

static inline double *
copy_mcol(const matrix *m,
          int col)
{
    double *v = (double *)malloc(m->rows * sizeof(*v));
    for (int i = 0; i < m->rows; i++) {
        v[i] = m->v[i*m->cols+col];
    }
    return v;
}

static inline double *
copy_mrow(const matrix *m,
          int row)
{
    double *v = (double *)malloc(m->cols * sizeof(*v));
    memcpy(v, &m->v[m->cols*row], m->cols * sizeof(*v));
    return v;
}

static inline void
swap_mrow(matrix *m,
          int r1,
          int r2)
{
    double *tmp = copy_mrow(m, r1);
    memcpy(&m->v[r1*m->cols], &m->v[r2*m->cols], m->cols*sizeof(m->v[0]));
    memcpy(&m->v[r2*m->cols], tmp, m->cols*sizeof(m->v[0]));
    free(tmp);
}

#endif
