/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef ICCLUT_H_
#define ICCLUT_H_

#include <dcamprof.h>
#include <profio.h>

struct icc_lut *
icclut_new(const m3x3 cam2xyz,
           const v3 prof_wp,
           const v3 pcs_wp,
           const v3 wb,
           double *trc[3],
           int trc_count,
           const double tc[],
           int tc_len,
           enum gc_type gc_type,
           enum tc_type tc_type,
           bool skip_tc_apply,
           chromalut_t *prof_lut,
           bool pcs_is_lab,
           int clut_side,
           const struct tone_rep_op_config *ntro_conf);

void
icclut_init(void);

#endif
