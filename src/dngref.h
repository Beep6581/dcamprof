/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef DNGREF_H_
#define DNGREF_H_

#include <profio.h>

// XYZ D50 to Prophoto RGB, according to DNG SDK reference code
/*
  0.7976749  0.1351917  0.0313534
  0.2880402  0.7118741  0.0000857
  0.0000000  0.0000000  0.8252100
*/
static const m3x3 dcp_xyz_D50_to_prophoto_rgb =
{{ { +1.3457530096, -0.2556031200, -0.0510246682 },
   { -0.5444259963, +1.5080962194, +0.0204718467 },
   { +0.0000000000, +0.0000000000, +1.2119675456 } }};

void
dcp_xy2temp(const double xy[2],
            double *temp,
            double *tint);

void
dcp_temp2xy(double temp,
            double tint,
            double xy[]);

m3x3
dcp_map_white(const v3 white1_xyz,
              const v3 white2_xyz);

v3
dcp_rgb2hsv(v3 rgb);

v3
dcp_hsv2rgb(const v3 hsv);

v3
dcp_lookup(const v3 *hsm,
           const uint32_t hsmdims[3],
           bool srgb_gamma,
           v3 cxyz,
           double bse_scale,
           const double tc[],
           int tc_len,
           bool *xyz_clip,
           bool *hsv_clip);

v3
dcp_lookup_rgb(const v3 *hsm,
               const uint32_t hsmdims[3],
               const bool srgb_gamma,
               const v3 rgb,
               bool *hsv_clip);

m3x3
dcp_find_xyz2cam(const struct profio_dcp *dcp,
                 const double white_xy[2],
                 m3x3 *fm,
                 m3x3 *rm,
                 const m3x3 *cc1,
                 const m3x3 *cc2,
                 const char cc_sig[],
                 m3x3 *cc);

void
dcp_neutral2xy(const struct profio_dcp *dcp,
               const m3x3 *cc1,
               const m3x3 *cc2,
               const char cc_sig[],
               const v3 *ab,
               v3 neutral,
               double xy[2]);

m3x3
dcp_make_cam2xyz_d50(const struct profio_dcp *dcp,
                     const m3x3 *cc1,
                     const m3x3 *cc2,
                     const char cc_sig[],
                     const v3 *ab,
                     const v3 wb);

v3 *
dcp_make_lut(const struct profio_dcp *dcp,
             const m3x3 *cc1,
             const m3x3 *cc2,
             const char cc_sig[],
             const v3 *ab,
             const v3 wb);

v3
dcp_d50(void);

double *
dcp_acr_tonecurve(int *len);

v3
dcp_apply_tonecurve(v3 rgb,
                    const double tc[],
                    int tc_len);

#endif
