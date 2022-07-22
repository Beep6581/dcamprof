 /*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <dcamprof.h>
#include <lut.h>

void
dnglut_test_discontinuity(v3 *hsm,
                          int hcount,
                          int scount,
                          int vcount,
                          const char tablename[],
                          bool allow_discontinuity_hue_shifts);

v3 *
dnglut_looktable_new(uint32_t hsmdims[3],
                     int hcount,
                     int scount,
                     int vcount,
                     const double tc[],
                     int tc_len,
                     enum tc_type tc_type,
                     bool skip_tc_apply,
                     bool srgb_gamma,
                     enum gc_type gc_type,
                     const struct tone_rep_op_config *ntro_conf);

v3 *
dnglut_huesatmap_new(uint32_t hsmdims[3],
                     int hcount,
                     int scount,
                     int vcount,
                     bool srgb_gamma,
                     chromalut_t *lut,
                     const v3 lut_d50,
                     const m3x3 cam2xyz_fm,
                     const m3x3 cam2xyz_lm);
