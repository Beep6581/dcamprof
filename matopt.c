/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdbool.h>

#include <matopt.h>
#include <elog.h>
#include <nmsimplex.h>

#define REPORT_PROGRESS(val) if (progress) { (progress(val, progress_arg)); }
#define debug(...)
//#define debug(...) elog(__VA_ARGS__)

struct simplex_base_arg {
    const v3 wp;
    const struct patch_set *ps;
    test_matrix_callback_t test_fun;
    void *test_arg;
};

struct simplex_fun_arg {
    const m3x3 base_matrix;
    int row;
    struct simplex_base_arg a;
};

static double
simplex_fun(double m[],
            void *arg)
{
    struct simplex_fun_arg *a = (struct simplex_fun_arg *)arg;
    m3x3 matrix = { .v = { { m[0], m[1], m[2]}, { m[3], m[4], m[5]}, { m[6], m[7], m[8]} } };
    return a->a.test_fun(&matrix, &a->a.wp, a->a.ps, a->a.test_arg);
}

static double
simplex_fun_diagonal(double diag[],
                     void *arg)
{
    struct simplex_fun_arg *a = (struct simplex_fun_arg *)arg;
    m3x3 matrix = a->base_matrix;
    matrix.v[0][0] = diag[0];
    matrix.v[1][1] = diag[1];
    matrix.v[2][2] = diag[2];
    return a->a.test_fun(&matrix, &a->a.wp, a->a.ps, a->a.test_arg);
}

static double
simplex_fun_row(double row[],
                void *arg)
{
    struct simplex_fun_arg *a = (struct simplex_fun_arg *)arg;
    m3x3 matrix = a->base_matrix;
    matrix.v[a->row][0] = row[0];
    matrix.v[a->row][1] = row[1];
    matrix.v[a->row][2] = row[2];
    return a->a.test_fun(&matrix, &a->a.wp, a->a.ps, a->a.test_arg);
}

static double
simplex_find_diag(m3x3 *m_,
                  struct simplex_base_arg base_arg)
{
    double m[3] = { m_->v[0][0],m_->v[1][1],m_->v[2][2] };
    struct simplex_fun_arg arg = { .base_matrix = *m_, .a = base_arg };
    double min = simplex(simplex_fun_diagonal, &arg, m, 3, 1.0e-7, 1, NULL);
    m3x3 matrix = arg.base_matrix;
    matrix.v[0][0] = m[0];
    matrix.v[1][1] = m[1];
    matrix.v[2][2] = m[2];
    *m_ = matrix;
    return min;
}

static double
simplex_find_row(m3x3 *m_,
                 int row,
                 struct simplex_base_arg base_arg)
{
    double m[3] = { m_->v[row][0],m_->v[row][1],m_->v[row][2] };
    struct simplex_fun_arg arg = { .base_matrix = *m_, .row = row, .a = base_arg };
    double min = simplex(simplex_fun_row, &arg, m, 3, 1.0e-7, 1, NULL);
    m3x3 matrix = arg.base_matrix;
    matrix.v[row][0] = m[0];
    matrix.v[row][1] = m[1];
    matrix.v[row][2] = m[2];
    *m_ = matrix;
    return min;
}

m3x3
matopt_find_matrix_brute_force2(const v3 wp,
                                const struct patch_set *ps,
                                test_matrix_callback_t test_matrix,
                                void *arg)
{
#define S 3
    m3x3 best_matrix;
    memset(&best_matrix, 0, sizeof(best_matrix));
    double old_min_de = test_matrix(&best_matrix, &wp, ps, arg);
    double min_de = old_min_de;
    double adj_range = 8.0;
    for (;;) {
        m3x3 base_matrix = best_matrix;
        double ref = adj_range;
#pragma omp parallel for
        for (int i = 0; i < S*S*S*S*S*S*S*S*S; i++) {
            const double rs[S] = { -1.0, 0, +1.0 };
            int q[9];
            q[0] = i % S;
            q[1] = (i / (S)) % S;
            q[2] = (i / (S*S)) % S;
            q[3] = (i / (S*S*S)) % S;
            q[4] = (i / (S*S*S*S)) % S;
            q[5] = (i / (S*S*S*S*S)) % S;
            q[6] = (i / (S*S*S*S*S*S)) % S;
            q[7] = (i / (S*S*S*S*S*S*S)) % S;
            q[8] = i / (S*S*S*S*S*S*S*S);
            m3x3 matrix = base_matrix;
            matrix.v[0][0] += ref * rs[q[0]];
            matrix.v[0][1] += ref * rs[q[1]];
            matrix.v[0][2] += ref * rs[q[2]];
            matrix.v[1][0] += ref * rs[q[3]];
            matrix.v[1][1] += ref * rs[q[4]];
            matrix.v[1][2] += ref * rs[q[5]];
            matrix.v[2][0] += ref * rs[q[6]];
            matrix.v[2][1] += ref * rs[q[7]];
            matrix.v[2][2] += ref * rs[q[8]];
            double de = test_matrix(&matrix, &wp, ps, arg);
#pragma omp critical
            {
                if (de < min_de) {
                    printf("residual error: %f\n", de);
                    printf("%+.4f %+.4f %+.4f\n", matrix.v[0][0], matrix.v[0][1], matrix.v[0][2]);
                    printf("%+.4f %+.4f %+.4f\n", matrix.v[1][0], matrix.v[1][1], matrix.v[1][2]);
                    printf("%+.4f %+.4f %+.4f\n\n", matrix.v[2][0], matrix.v[2][1], matrix.v[2][2]);
                    min_de = de;
                    best_matrix = matrix;
                }
            }
        }
        elog("loop iteration complete!\n");
        if (min_de == old_min_de) {
            adj_range *= 0.5;
            elog("adjust range: -%f - +%f\n", adj_range, adj_range);
        } else {
            bool differs = false;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    if (round(100 * best_matrix.v[row][col]) !=
                        round(100 * base_matrix.v[row][col]))
                    {
                        differs = true;
                    }
                }
            }
            if (!differs) {
                break;
            }
        }
        old_min_de = min_de;
    }

#undef S
    return best_matrix;
}

m3x3
matopt_find_matrix_brute_force(const v3 wp,
                               const struct patch_set *ps,
                               test_matrix_callback_t test_matrix,
                               void *arg)
{
#define S 5
    m3x3 best_matrix;
    memset(&best_matrix, 0, sizeof(best_matrix));
    double old_min_de = test_matrix(&best_matrix, &wp, ps, arg);
    double min_de = old_min_de;
    double adj_range = 8.0;
    for (;;) {
        m3x3 base_matrix = best_matrix;
        double ref = adj_range;
#pragma omp parallel for
        for (int i = 0; i < S*S*S*S*S*S*S*S*S; i++) {
            const double rs[S] = { -1.0, -0.5, 0, +0.5, +1.0 };
            int q[9];
            q[0] = i % S;
            q[1] = (i / (S)) % S;
            q[2] = (i / (S*S)) % S;
            q[3] = (i / (S*S*S)) % S;
            q[4] = (i / (S*S*S*S)) % S;
            q[5] = (i / (S*S*S*S*S)) % S;
            q[6] = (i / (S*S*S*S*S*S)) % S;
            q[7] = (i / (S*S*S*S*S*S*S)) % S;
            q[8] = i / (S*S*S*S*S*S*S*S);
            m3x3 matrix = base_matrix;
            matrix.v[0][0] += ref * rs[q[0]];
            matrix.v[0][1] += ref * rs[q[1]];
            matrix.v[0][2] += ref * rs[q[2]];
            matrix.v[1][0] += ref * rs[q[3]];
            matrix.v[1][1] += ref * rs[q[4]];
            matrix.v[1][2] += ref * rs[q[5]];
            matrix.v[2][0] += ref * rs[q[6]];
            matrix.v[2][1] += ref * rs[q[7]];
            matrix.v[2][2] += ref * rs[q[8]];
            double de = test_matrix(&matrix, &wp, ps, arg);
#pragma omp critical
            {
                if (de < min_de) {
                    printf("residual error: %f\n", de);
                    printf("%+.4f %+.4f %+.4f\n", matrix.v[0][0], matrix.v[0][1], matrix.v[0][2]);
                    printf("%+.4f %+.4f %+.4f\n", matrix.v[1][0], matrix.v[1][1], matrix.v[1][2]);
                    printf("%+.4f %+.4f %+.4f\n\n", matrix.v[2][0], matrix.v[2][1], matrix.v[2][2]);
                    min_de = de;
                    best_matrix = matrix;
                }
            }
        }
        elog("loop iteration complete!\n");
        if (min_de == old_min_de) {
            adj_range *= 0.5;
            elog("adjust range: -%f - +%f\n", adj_range, adj_range);
        } else {
            bool differs = false;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    if (round(100 * best_matrix.v[row][col]) !=
                        round(100 * base_matrix.v[row][col]))
                    {
                        differs = true;
                    }
                }
            }
            if (!differs) {
                break;
            }
        }
        old_min_de = min_de;
    }

#undef S
    return best_matrix;
}

double
matopt_refine_matrix(m3x3 *m_,
                     const v3 wp,
                     const struct patch_set *ps,
                     test_matrix_callback_t test_matrix,
                     void *arg_)
{
    double m[9] = { m_->v[0][0],m_->v[0][1],m_->v[0][2],
                    m_->v[1][0],m_->v[1][1],m_->v[1][2],
                    m_->v[2][0],m_->v[2][1],m_->v[2][2] };
    struct simplex_base_arg base_arg = { .wp = wp, .ps = ps, .test_fun = test_matrix, .test_arg = arg_ };
    struct simplex_fun_arg arg = { .a = base_arg };
    double old_min = 100000000000;
    double min;
    double old_m[9];
    memcpy(old_m, m, sizeof(old_m));
    for (int i = 0; i < 10; i++) {
        min = simplex(simplex_fun, &arg, m, 9, 1.0e-7, 0.01, NULL);
        if (min + 1e-5 > old_min) {
            if (min > old_min) {
                memcpy(m, old_m, sizeof(old_m));
            }
            break;
        }
        old_min = min;
        memcpy(old_m, m, sizeof(old_m));
    }
    m3x3 matrix = { .v = { { m[0], m[1], m[2]}, { m[3], m[4], m[5]}, { m[6], m[7], m[8]} } };
    *m_ = matrix;
    return min;
}

m3x3
matopt_find_matrix(const v3 wp,
                   const struct patch_set *ps,
                   const double smallest_per_row[],
                   test_matrix_callback_t test_matrix,
                   void *arg,
                   void (*progress)(double prog, void *arg),
                   void *progress_arg)
{
    REPORT_PROGRESS(0);
    m3x3 m = { .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

    if (smallest_per_row != NULL) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (smallest_per_row[i] > m.v[i][j]) m.v[i][j] = smallest_per_row[i];
            }
        }
    }
    struct simplex_base_arg base_arg = { .wp = wp, .ps = ps, .test_fun = test_matrix, .test_arg = arg };
    simplex_find_diag(&m, base_arg);
    REPORT_PROGRESS(0.2);
    simplex_find_row(&m, 1, base_arg);
    REPORT_PROGRESS(0.4);
    simplex_find_row(&m, 0, base_arg);
    REPORT_PROGRESS(0.6);
    simplex_find_row(&m, 2, base_arg);
    REPORT_PROGRESS(0.8);
    matopt_refine_matrix(&m, wp, ps, test_matrix, arg);
    REPORT_PROGRESS(1.0);

    return m;
}
