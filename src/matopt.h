/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef MATOPT_H_
#define MATOPT_H_

#include <matmath.h>
#include <target.h>

typedef double (*test_matrix_callback_t)(const m3x3 *cam2xyz,
                                         const v3 *wp,
                                         const struct patch_set *ps,
                                         void *arg);

double
matopt_refine_matrix(m3x3 *m_,
                     const v3 wp,
                     const struct patch_set *ps,
                     test_matrix_callback_t test_matrix,
                     void *arg);

m3x3
matopt_find_matrix(const v3 wp,
                   const struct patch_set *ps,
                   const double smallest_per_row[],
                   test_matrix_callback_t test_matrix,
                   void *arg,
                   void (*progress)(double prog, void *arg),
                   void *progress_arg);

m3x3
matopt_find_matrix_brute_force2(const v3 wp,
                                const struct patch_set *ps,
                                test_matrix_callback_t test_matrix,
                                void *arg);

m3x3
matopt_find_matrix_brute_force(const v3 wp,
                               const struct patch_set *ps,
                               test_matrix_callback_t test_matrix,
                               void *arg);

#endif
