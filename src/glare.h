/*
 * (c) Copyright 2015, 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef GLARE_H_
#define GLARE_H_

#include <dcamprof.h>

struct glare_data {
    double avg_black_g;
    double avg_black_y;
    double avg_white_g;
    double avg_white_y;
};

void
glare_test(const struct patch_set *ps,
           struct glare_data *gd);

void
glare_add(struct patch_set *ps,
          v3 wp,
          double mul);

bool
glare_match(const struct test_chart_spec *tcs, // should contain neutrals
            struct patch_set *ps, // should be flatfield corrected
            const struct observer *obs,
            const spectrum_t *illuminant,
            const v3 *illuminant_wp,
            const char report_dir[],
            bool verbose,
            bool exit_on_failure);

double
flatfield_correct_targets(const struct test_chart_spec *tcs,
                          struct patch_set *ps,
                          struct patch_set *ps_ex,
                          bool verbose,
                          bool exit_on_failure,
                          bool *bad_correction);

void
glare_linearize_target(const struct patch_set *ti1, // should contain neutrals
                       struct patch_set *ps, // should be flatfield corrected
                       const struct test_chart_layout *layout,
                       const char report_dir[]);

#endif
