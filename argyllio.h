/*
 * (c) Copyright 2015, 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef ARGYLLIO_H_
#define ARGYLLIO_H_

#include <stdbool.h>

#include <target.h>

struct argyll_patch {
    double xyz[3];
    double rgb[3];
    double stdev[3];
};

struct argyll_ti3 {
    int patch_count;
    struct argyll_patch patches[];
};

struct patch_set *
argyllio_ti3_parse(const char filename[],
                   bool error_is_fatal,
                   bool verbose);

bool
argyllio_ti3_write(const char filename[],
                   bool print_spectral_data,
                   const struct patch_set *ps,
                   bool error_is_fatal);

bool
argyllio_ti1_write(const char filename[],
                   struct patch_set *ps,
                   bool fill_in_xyz,
                   bool error_is_fatal);

int
textio_process_exclude_list(struct patch_set *ps,
                            const char filename[],
                            int patch_remapped_indexes[]);

int
textio_process_keep_list(struct patch_set *ps,
                         const char filename[]);

#endif
