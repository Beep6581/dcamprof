/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef LOOK_H_
#define LOOK_H_

#include <dcamprof.h>

struct look_tone_rep_op_t_;
typedef struct look_tone_rep_op_t_ look_tone_rep_op_t;

void
look_neutral_tone_rep_op_tc_analysis(const double tc[],
                                     int tc_len,
                                     double *cv,
                                     double *cs);

double
look_lookop_curve_y(const struct lookop_curve *c,
                    double x);

look_tone_rep_op_t *
look_tone_rep_op_new(const m3x3 rgb2xyz,
                     const double tc[],
                     int tc_len,
                     enum gc_type gc_type,
                     const struct tone_rep_op_config *conf);

void *
look_tone_rep_op_get_lcmsh02(look_tone_rep_op_t *state);

void
look_tone_rep_op_delete(look_tone_rep_op_t *state);

void
look_tone_rep_op(v3 *out,
                 v3 xyz, // D50!
                 look_tone_rep_op_t *state);

v3
look_wp(void);

void
look_apply_target_adjustment(struct patch_set *ps,
                             const struct target_adjustment_config *tadj,
                             v3 wp);

#endif
