/*
 * (c) Copyright 2015, 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef TARGET_H_
#define TARGET_H_

#include <matmath.h>
#include <colmath.h>
#include <jsonio.h>
#include <dcamprof.h>

bool
target_isvalid_target_name(const char name[],
                           bool *is_measured_data);

struct grid_entry {
    double u;
    double v;
    struct cam2xyz_patch *p;
};

struct patch_grid {
    int rows;
    int cols;
    struct grid_entry p[];
};

struct patch_set *
target_make_virtual(const char name[],
                    const struct dcam_ssf *ssf,
                    const struct observer *obs,
                    const spectrum_t *illuminant,
                    const spectrum_t *ref_illuminant,
                    bool prefer_cat,
                    double grid_spacing,
                    bool fill_grid,
                    bool skip_border,
                    struct patch_grid **patch_grid, // FIXME: obsolete, clean up code
                    const char report_dir[]);

void
target_delete_patch_set(struct patch_set *ps);

struct patch_set *
target_merge(struct patch_set *inps[], // will be deleted
             int ps_count,
             double min_chroma_distance,
             bool error_is_fatal);

struct patch_set *
target_copy(const struct patch_set *ps);

void
target_set_name(struct patch_set *ps,
                const char class_name[]);

int
target_class_or_patch_to_idx(const struct patch_set *ps,
                             const char name[],
                             bool *is_class);

int
target_name2idx(const struct patch_set *ps,
                const char name[],
                const char class_name[]);

const char *
target_patch_name(char buf[PATCH_STRID_SIZE],
                  const struct patch_set *ps,
                  int idx);

const char *
target_measured_list(void);

const char *
target_generated_list(void);

void
target_print(const char prefix[],
             const struct patch_set *ps,
             const v3 wp,
             const struct observer *obs,
             const char report_dir[]);

struct patch_set *
target_parse_text(const char filename[],
                  const char class_name[],
                  const struct observer *obs,
                  const spectrum_t *illuminant,
                  double scale,
                  double clip_level,
                  int low,
                  int high,
                  int spacing,
                  bool column_layout,
                  bool error_is_fatal);

#endif
