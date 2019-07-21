/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef LUT_H_
#define LUT_H_

#include <stdio.h>

#include <dcamprof.h>

v3
xyz2lutspace(v3 xyz);

v3
lutspace2xyz(v3 lutspace);

bool
chromalut_isnoop(chromalut_t *lut);

chromalut_t *
chromalut_new_from_data(v3 *cp[3],
                        int cp_count,
                        double regularization[3],
                        const m3x3 fm,
                        double compress_factor);

chromalut_t *
chromalut_new(const m3x3 fm,
              const v3 wb,
              const v3 d50,
              const struct patch_set *normalized_ps,
              const double min_chroma_distance,
              const double compress_factor);

void
chromalut_delete(chromalut_t *lut);

v3
chromalut_lookup(chromalut_t *lut,
                 v3 xyz);

v3
chromalut_lookup_noclip(chromalut_t *lut,
                        v3 xyz);

v3
chromalut_lookup_diff(chromalut_t *lut,
                      double u,
                      double v);

/*
v3
chromalut_reverse_lookup(chromalut_t *lut,
                         v3 xyz);
*/

void
chromalut_json_print(FILE *stream,
                     const char indent[],
                     const char name[],
                     chromalut_t *lut);

#endif
