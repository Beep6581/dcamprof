/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef GAMUT_H_
#define GAMUT_H_

#include <matmath.h>
#include <dcamprof.h>

void
gamut_compress_state_alloc(struct gamut_compression *gc,
                           m3x3 rgb2xyz,
                           v3 wp);

void
gamut_compress_state_free(struct gamut_compression *gc);

v3
gamut_compress(const struct gamut_compression *gc,
               m3x3 rgb2xyz,
               m3x3 xyz2rgb,
               v3 xyz,
               v3 wp,
               void *lcms2_h02);

m3x3
gamut_remap_xyz2rgb(m3x3 xyz2rgb,
                    v3 wp);

bool
gamut_isvalid(v3 xyz);

double
gamut_max_valid_chroma(v3 jch,
                       v3 wp,
                       double clip_hi_limit,
                       void *lcms2_h02);

v3
gamut_clip_to_valid(v3 xyz,
                    v3 wp);

v3
gamut_acr_rgb_clip(v3 rgb);

v3
gamut_saturated_rgb_clip(v3 rgb);

void
gamut_init(const struct observer *obs,
           const spectrum_t *illuminant);

void
gamut_exit(void);

#endif
