/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdbool.h>
#include <assert.h>

#include <bisection.h>
#include <dngref.h>
#include <gamut.h>
#include <interp.h>
#include <lut3d.h>
#include <elog.h>
#include <look.h>

static struct {
    int locus_poly_count;
    double *locus_poly_x;
    double *locus_poly_y;
    double prophoto_xy[2][3];
} glob;

struct gamut_compress_state_t_ {
    int mixlut_cs[4];
    double *mixlut[4];
    double mixlim[3][2];
    double *hue2max_curve;
    double *hue2min_curve;
    int hue_curve_len;
};

struct find_gamut_edge_arg {
    double j, h, limit;
    v3 wp;
    void *h02;
    const m3x3 *xyz2rgb;
    int poly_count;
    const double *poly_x;
    const double *poly_y;
    double last_valid_chroma;
    double hilim;
    double lolim;
    bool alt_space;
    m3x3 outer_rgb2xyz;
};

static v3
lsh2xyz(v3 lsh,
        const m3x3 rgb2xyz)
{
    v3 hsl = (v3){{ lsh.v[2], lsh.v[1], lsh.v[0] }};
    v3 rgb = hsv2rgb(hsl);
    return m3x3_mul_v3(rgb2xyz, rgb);
}

static v3
xyz2lsh(v3 xyz,
        const m3x3 xyz2rgb)
{
    v3 rgb = m3x3_mul_v3(xyz2rgb, xyz);
    v3 hsl = rgb2hsv(rgb);
    return (v3){{ hsl.v[2], hsl.v[1], hsl.v[0] }};
}

static double
find_gamut_edge_func(double x,
                     void *arg)
{
    struct find_gamut_edge_arg *a = (struct find_gamut_edge_arg *)arg;
    v3 xyz;
    if (a->alt_space) {
        xyz = lsh2xyz((v3){{ a->j, x, a->h }}, a->outer_rgb2xyz);
    } else {
        xyz = jch2xyz_s((v3){{ a->j, x, a->h }}, a->h02);
    }
    if (!v3_isfinite(xyz) || xyz.v[0] < 0 || xyz.v[1] < 0 || xyz.v[2] < 0) {
        return a->limit + x;
    }
    if (a->poly_count > 0) {
        v3 xyY = xyz2xyY(xyz);
        if (!point_inside_polygon(a->poly_x, a->poly_y, a->poly_count, xyY.v[0], xyY.v[1])) {
            return a->limit + x;
        }
    }
    if (a->xyz2rgb != NULL) {
        v3 rgb = m3x3_mul_v3(*a->xyz2rgb, xyz);
        for (int i = 0; i < 3; i++) {
            if (rgb.v[i] > a->hilim || rgb.v[i] < a->lolim) {
                return a->limit + x;
            }
        }
    }
    a->last_valid_chroma = x;
    return a->limit - x;
}

static double
find_gamut_edge(v3 jch,
                const m3x3 outer_rgb2xyz,
                const m3x3 *xyz2rgb,
                int poly_count,
                const double poly_x[],
                const double poly_y[],
                double hilim,
                double lolim,
                void *lcms2_h02,
                v3 wp,
                bool alt_space)
{
    struct find_gamut_edge_arg arg = {
        .j = jch.v[0],
        .h = jch.v[2],
        .wp = wp,
        .h02 = lcms2_h02,
        .xyz2rgb = xyz2rgb,
        .poly_count = poly_count,
        .poly_x = poly_x,
        .poly_y = poly_y,
        .hilim = hilim,
        .lolim = lolim,
        .alt_space = alt_space,
        .outer_rgb2xyz = outer_rgb2xyz,
        .limit = 300 // JCh: chroma 300 covers covers both locus and prophotoRGB
    };
    if (alt_space) {
        arg.limit = 1.0;
    }
    interval_halving_min(find_gamut_edge_func, &arg, 0, arg.limit, 0.0001, 1000);
    return arg.last_valid_chroma;
}

static double *
make_hsv_mix_lut(const m3x3 outer_rgb2xyz,
                 const m3x3 inner_xyz2rgb,
                 const int cs,
                 const int blur_len,
                 const double lo_sat_blend_lim,
                 const double hi_sat_blend_lim)
{
    double *clut = calloc(1, sizeof(clut[0]) * cs * cs * cs);

#pragma omp parallel for
    for (int r_ = 0; r_ < cs; r_++) {
        for (int g_ = 0; g_ < cs; g_++) {
            for (int b_ = 0; b_ < cs; b_++) {
                v3 rgb = {{ (double)r_ / (cs-1), (double)g_ / (cs-1), (double)b_ / (cs-1) }};

                for (int i = 0; i < 3; i++) {
                    rgb.v[i] = srgb_gamma_inverse(rgb.v[i]);
                }
                v3 xyz = m3x3_mul_v3(outer_rgb2xyz, rgb);
                v3 irgb = m3x3_mul_v3(inner_xyz2rgb, xyz);
                double alt_mix;
                if (v3_min(irgb) < -0.001 || v3_max(irgb) > 1.001) {
                    alt_mix = 1;
                } else {
                    irgb = v3_clip01(irgb);
                    irgb.v[0] = srgb_gamma_forward(irgb.v[0]);
                    irgb.v[1] = srgb_gamma_forward(irgb.v[1]);
                    irgb.v[2] = srgb_gamma_forward(irgb.v[2]);
                    v3 hsv = rgb2hsv(irgb);
                    const double lolim = lo_sat_blend_lim;
                    const double hilim = hi_sat_blend_lim;
                    if (hsv.v[1] < lolim) {
                        alt_mix = 0;
                    } else if (hsv.v[1] > hilim) {
                        alt_mix = 1;
                    } else {
                        alt_mix = scurve((hsv.v[1] - lolim) / (hilim - lolim));
                    }
                }

                clut[cs*cs*r_+cs*g_+b_] = alt_mix;
            }
        }
    }

    // approx gaussian blur
    for (int k = 0; k < 3; k++) {
        boxblur_lut(clut, cs, blur_len);
    }
    return clut;
}

// 0 .. 1 value
static double
get_hsv_mix_value(v3 rgb, // gamma version
                  const struct gamut_compression *gc)
{
    gamut_compress_state_t *state = gc->state;
    double mix = lut3d_lookup_1d(state->mixlut[0], state->mixlut_cs[0], rgb);
    v3 hsl = rgb2hsl(rgb);
    const double test_val = hsl.v[2];
    for (size_t i = 0; i < sizeof(state->mixlim) / sizeof(state->mixlim[0]); i++) {
        double alt_mix;
        const double lolim = state->mixlim[i][0];
        const double hilim = state->mixlim[i][1];
        if (test_val < lolim) {
            alt_mix = 0;
        } else if (test_val > hilim) {
            alt_mix = 1;
        } else {
            alt_mix = scurve((test_val - lolim) / (hilim - lolim));
        }
        double mix1 = lut3d_lookup_1d(state->mixlut[1+i], state->mixlut_cs[1+i], rgb);
        mix = mix * alt_mix + mix1 * (1 - alt_mix);
    }
    return mix;
}

static int
get_hsv_compress_curve(double x[5],
                       double y[5],
                       double fact)
{
    if (fact <= 1) {
        x[0] = 0;
        y[0] = 0;
        x[1] = 1;
        y[1] = 1;
        return 2;
    }
    x[4] = 1;
    y[4] = 1;
    x[3] = 1 - 0.6 * (1 - (1.0 / fact));
    y[3] = 1.0 / fact;

    double k = 0.8;
    double A = x[3];
    double B = y[3];
    double c = -(B-k*A)/(k-1);
    x[2] = c;
    y[2] = c;
    double d = c - 0.01;
    if (d < 0) d = c * 0.5;
    x[1] = d;
    y[1] = d;
    x[0] = 0;
    y[0] = 0;

    return 5;
}

static v3
gamut_compress_rgb(const struct gamut_compression *gc,
                   const m3x3 rgb2xyz,
                   const m3x3 xyz2rgb,
                   v3 rgb,
                   v3 wp, // d50
                   void *lcms2_h02,
                   bool use_alt_space)
{
    if (gc->hsv.enabled) { // HSV value compression step
        rgb.v[0] = srgb_gamma_forward(rgb.v[0]);
        rgb.v[1] = srgb_gamma_forward(rgb.v[1]);
        rgb.v[2] = srgb_gamma_forward(rgb.v[2]);
        const v3 hsv = rgb2hsv(rgb);

        v3 alt_rgb = rgb;
        {
            double curve_factor;
            {
                double extension = curve_y(hsv.v[0], gc->state->hue2max_curve, gc->state->hue_curve_len);
                if (extension < 1.0) extension = 1.0;

                const double cutoff = gc->hsv.top.cutoff_limit;
                const double linear = gc->hsv.top.linear_limit;
                const double output = gc->hsv.top.output_limit;
                if (extension > cutoff) {
                    extension = output;
                } else if (extension > linear) {
                    const double cf = (output - linear) / (cutoff - linear);
                    double x = (extension - linear) / (cutoff - linear);
                    x = (x - (1-cf) * pow(x, 1/(1-cf)));
                    x *= cutoff - linear;
                    extension = linear + x;
                }
                curve_factor = extension;
            }
            double handle_x[5], handle_y[5];
            int handle_count = get_hsv_compress_curve(handle_x, handle_y, curve_factor);

            int mini, midi, maxi;
            v3_minmidmax(alt_rgb, &mini, &midi, &maxi);
            double new_max;
            double test_val = alt_rgb.v[maxi];
            cubic_spline(handle_x, handle_y, handle_count, &test_val, &new_max, 1);
            double diff = alt_rgb.v[maxi] - new_max;
            if (diff < 0) {
                new_max = alt_rgb.v[maxi];
                diff = 0;
            }
            double scale = new_max / alt_rgb.v[maxi];
            double new_min = alt_rgb.v[mini] * scale + gc->hsv.top.min_lift * diff;
            if (new_min >= new_max) {
                alt_rgb = (v3){{ new_max, new_max, new_max }};
            } else {
                double mid_mul = (alt_rgb.v[midi] - alt_rgb.v[mini]) / (alt_rgb.v[maxi] - alt_rgb.v[mini]);
                double new_mid = new_min + (new_max - new_min) * mid_mul;
                alt_rgb.v[mini] = new_min;
                alt_rgb.v[midi] = new_mid;
                alt_rgb.v[maxi] = new_max;
            }
        }

        {
            double curve_factor;
            {
                double extension = -curve_y(hsv.v[0], gc->state->hue2min_curve, gc->state->hue_curve_len);
                if (extension < 0) extension = 0;

                const double cutoff = -gc->hsv.bottom.cutoff_limit;
                const double linear = -gc->hsv.bottom.linear_limit;
                const double output = -gc->hsv.bottom.output_limit;
                if (extension > cutoff) {
                    extension = output;
                } else if (extension > linear) {
                    const double cf = (output - linear) / (cutoff - linear);
                    double x = (extension - linear) / (cutoff - linear);
                    x = (x - (1-cf) * pow(x, 1/(1-cf)));
                    x *= cutoff - linear;
                    extension = linear + x;
                }
                curve_factor = 1 + extension;
            }
            double handle_x[5], handle_y[5];
            int handle_count = get_hsv_compress_curve(handle_x, handle_y, curve_factor);

            int mini, midi, maxi;
            v3_minmidmax(alt_rgb, &mini, &midi, &maxi);
            double new_min;
            double test_val = 1 - alt_rgb.v[mini];
            cubic_spline(handle_x, handle_y, handle_count, &test_val, &new_min, 1);
            new_min = 1 - new_min;
            double diff = new_min - alt_rgb.v[mini];
            if (diff < 0) {
                new_min = alt_rgb.v[mini];
                diff = 0;
            }
            double new_max = alt_rgb.v[maxi] + gc->hsv.bottom.max_lift * diff;
            if (new_max > 1) new_max = 1;
            if (new_min >= new_max) {
                alt_rgb = (v3){{ new_max, new_max, new_max }};
            } else {
                double mid_mul = (alt_rgb.v[midi] - alt_rgb.v[mini]) / (alt_rgb.v[maxi] - alt_rgb.v[mini]);
                double new_mid = new_min + (new_max - new_min) * mid_mul;
                alt_rgb.v[mini] = new_min;
                alt_rgb.v[midi] = new_mid;
                alt_rgb.v[maxi] = new_max;
            }
        }

        double alt_mix = get_hsv_mix_value(rgb, gc);
        for (int i = 0; i < 3; i++) {
            rgb.v[i] = rgb.v[i] * (1 - alt_mix) + alt_rgb.v[i] * alt_mix;
            rgb.v[i] = srgb_gamma_inverse(rgb.v[i]);
        }
    }

    v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
    v3 jch;
    if (use_alt_space) {
        jch = xyz2lsh(xyz, rgb2xyz);
    } else {
        jch = xyz2jch_s(xyz, lcms2_h02);
        if (!v3_isfinite(jch)) {
            elog("Bug: out of gamut values fed to CIECAM02.\n");
            v3_print("rgb", rgb);
            v3_print("xyz", xyz);
            abort();
        }
    }

    // for the current J and h, find out chroma for inner, clip and destination gamuts
    double iC = find_gamut_edge(jch, rgb2xyz, gc->inner.xyz2rgb_enabled ? &gc->inner.xyz2rgb : NULL, gc->inner.locus_enabled ? glob.locus_poly_count : 0, glob.locus_poly_x, glob.locus_poly_y, gc->hi_rgb_lim, gc->lo_rgb_lim, lcms2_h02, wp, use_alt_space);
    double dC = find_gamut_edge(jch, rgb2xyz, gc->dest.xyz2rgb_enabled ? &gc->dest.xyz2rgb : NULL, gc->dest.locus_enabled ? glob.locus_poly_count : 0, glob.locus_poly_x, glob.locus_poly_y, gc->hi_rgb_lim, gc->lo_rgb_lim, lcms2_h02, wp, use_alt_space);
    double cC = find_gamut_edge(jch, rgb2xyz, gc->clip.xyz2rgb_enabled ? &gc->clip.xyz2rgb : NULL, gc->clip.locus_enabled ? glob.locus_poly_count : 0, glob.locus_poly_x, glob.locus_poly_y, gc->hi_rgb_lim, gc->lo_rgb_lim, lcms2_h02, wp, use_alt_space);

    iC *= gc->inner.chroma_scale;
    dC *= gc->dest.chroma_scale;
    cC *= gc->clip.chroma_scale;

    // make sure that iC <= dC <= cC
    if (dC > cC) dC = cC;
    if (iC > dC) iC = dC;

    if (jch.v[1] <= iC) {
        // already inside inner gamut, do nothing
        return rgb;
    } else if (jch.v[1] >= cC) {
        jch.v[1] = dC; // clip to destination gamut
    } else {
        // handle unexpected gamut edge order
        if (cC == iC) {
            jch.v[1] = dC; // clip to destination gamut
        } else {
            // compress chroma
            double x = (jch.v[1] - iC) / (cC - iC); // x = [0..1], 0 at inner chroma, 1 at clip chroma
            double d = (dC - iC) / (cC - iC); // target range = [0..d]
            // fit x to target range

            /*
              rolloff function, idea is f'(0) == 1 (starting point) and f'(q) == 0 (end point)
              f(x) = x-x^k
              f'(x) = 1-k*x^(k-1)
              f'(q) = 0 => q = 1/k^(1/(k-1))

              f(q) / q = d =>
                ((1/k^(1/(k-1)))-(1/k^(1/(k-1)))^k) / (1/k^(1/(k-1))) = d =>
                k = 1/(1-d)
            */
            double k = 1.0/(1.0 - d);
            double q = 1.0/pow(k, 1.0/(k-1.0)); // q is f'(q) == 0
            x *= q; // x range is 0..1, make it 0..q
            x = x - pow(x, k);
            x /= q;
            jch.v[1] = iC + x * (cC - iC);
        }
    }
    if (use_alt_space) {
        xyz = lsh2xyz(jch, xyz2rgb);
    } else {
        xyz = jch2xyz_s(jch, lcms2_h02);
    }
    return m3x3_mul_v3(xyz2rgb, xyz);
}

v3
gamut_compress(const struct gamut_compression *gc,
               const m3x3 rgb2xyz,
               const m3x3 xyz2rgb,
               v3 xyz,
               v3 wp, // d50
               void *lcms2_h02)
{
    v3 rgb = m3x3_mul_v3(xyz2rgb, xyz);
    rgb = gamut_saturated_rgb_clip(rgb);
    const double hue = rgb2hsl(rgb).v[0] * 360;

    // An alternate color space is used in the cyan to magenta range, as CIECAM02 is problematic there for high saturations
    double alt_mix;
    {
        // alt_mix = 1 in the Cyan to Magenta HSL-hue range, else 0 with a crossover zone
        double hsl_lolim = 180; // Cyan
        double hsl_hilim = 300; // Magenta
        double xover = 60;
        if (hue > hsl_hilim+xover || hue < hsl_lolim-xover) {
            alt_mix = 0;
        } else if (hue < hsl_lolim) {
            double x = (hue - (hsl_lolim-xover)) / xover;
            x = scurve(x);
            alt_mix = x;
        } else if (hue > hsl_hilim) {
            double x = (hue - hsl_hilim) / xover;
            x = scurve(x);
            alt_mix = 1-x;
        } else {
            alt_mix = 1;
        }
    }
    bool no_change = true;
    v3 rgb_norm;
    if (alt_mix != 1) {
        rgb_norm = gamut_compress_rgb(gc, rgb2xyz, xyz2rgb, rgb, wp, lcms2_h02, false);
        no_change = v3_isequal(rgb_norm, rgb);
    } else {
        rgb_norm = (v3){{0,0,0}};
    }
    bool alt_no_change = true;
    v3 rgb_alt;
    if (alt_mix != 0) {
        rgb_alt = gamut_compress_rgb(gc, rgb2xyz, xyz2rgb, rgb, wp, lcms2_h02, true);
        alt_no_change = v3_isequal(rgb_alt, rgb);
    } else {
        rgb_alt = (v3){{0,0,0}};
    }
    if (no_change && alt_no_change) {
        return xyz;
    }
    // make s-curve transition with sRGB gamma (rather than linear) for better perceptual mapping of the curve
    rgb.v[0] = srgb_gamma_forward(rgb_norm.v[0]) * (1 - alt_mix) + srgb_gamma_forward(rgb_alt.v[0]) * alt_mix;
    rgb.v[1] = srgb_gamma_forward(rgb_norm.v[1]) * (1 - alt_mix) + srgb_gamma_forward(rgb_alt.v[1]) * alt_mix;
    rgb.v[2] = srgb_gamma_forward(rgb_norm.v[2]) * (1 - alt_mix) + srgb_gamma_forward(rgb_alt.v[2]) * alt_mix;
    for (int i = 0; i < 3; i++) {
        rgb.v[i] = srgb_gamma_inverse(rgb.v[i]);
    }
    rgb = gamut_saturated_rgb_clip(rgb);
    xyz = m3x3_mul_v3(rgb2xyz, rgb);

    //v3 ltu = xyz2lutspace(xyz); printf("%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
    return xyz;
}

void
gamut_compress_state_alloc(struct gamut_compression *gc,
                           const m3x3 rgb2xyz,
                           const v3 wp)
{
    if (!gc->hsv.enabled) {
        return;
    }
    gamut_compress_state_t *state = calloc(1, sizeof(*state));

    { // mix LUTs
        double mixlim[3][2] = { { 0.4, 0.7 }, { 0.2, 0.6 }, { 0.1, 0.4 } };
        memcpy(state->mixlim, mixlim, sizeof(mixlim));
        double bs[4] = { 1, 2, 2, 1 };
        state->mixlut_cs[0] = 17;
        state->mixlut_cs[1] = 33;
        state->mixlut_cs[2] = 65;
        state->mixlut_cs[3] = 65;
        for (int i = 0; i < 4; i++) {
            state->mixlut[i] = make_hsv_mix_lut(rgb2xyz, gc->inner.xyz2rgb, state->mixlut_cs[i], bs[i],
                                                gc->hsv.lo_sat_blend_lim, gc->hsv.hi_sat_blend_lim);
        }
    }

    { // hue2max curves
        state->hue_curve_len = 65536;
        state->hue2max_curve = malloc(sizeof(state->hue2max_curve[0]) * state->hue_curve_len);
        state->hue2min_curve = malloc(sizeof(state->hue2min_curve[0]) * state->hue_curve_len);
        for (int h_ = 0; h_ < state->hue_curve_len; h_++) {
            v3 max_hsv = {{ (double)h_ / state->hue_curve_len, 1, 1 }};
            v3 rgb = hsv2rgb(max_hsv);
            v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
            xyz = gamut_clip_to_valid(xyz, wp);
            v3 irgb = m3x3_mul_v3(gc->inner.xyz2rgb, xyz);
            double max = v3_max(irgb);
            double min = v3_min(irgb);

            if (min < gc->hsv.bottom.cutoff_limit) min = gc->hsv.bottom.cutoff_limit;
            if (max > gc->hsv.top.cutoff_limit) max = gc->hsv.top.cutoff_limit;

            state->hue2max_curve[h_] = max;
            state->hue2min_curve[h_] = min;
        }

        for (int k = 0; k < 2; k++) {
            double *curve = k == 0 ? state->hue2max_curve : state->hue2min_curve;
            const int radius = state->hue_curve_len / 72;
            double arr[state->hue_curve_len+2*radius];
            memcpy(&arr[radius], curve, state->hue_curve_len * sizeof(curve[0]));
            memcpy(&arr[0], &curve[state->hue_curve_len-radius], radius * sizeof(curve[0]));
            memcpy(&arr[radius+state->hue_curve_len], curve, radius * sizeof(curve[0]));
            for (int k = 0; k < 3; k++) {
                boxblur_1d(arr, state->hue_curve_len+2*radius, radius);
            }
            memcpy(curve, &arr[radius], state->hue_curve_len * sizeof(curve[0]));
        }
    }

    gc->state = state;
}

void
gamut_compress_state_free(struct gamut_compression *gc)
{
    if (gc->state == NULL) {
        return;
    }
    for (size_t i = 0; i < sizeof(gc->state->mixlut) / sizeof(gc->state->mixlut[0]); i++) {
        free(gc->state->mixlut[i]);
    };
    free(gc->state->hue2max_curve);
    free(gc->state->hue2min_curve);
    free(gc->state);
    gc->state = NULL;
}

m3x3
gamut_remap_xyz2rgb(m3x3 xyz2rgb,
                    v3 wp)
{
    m3x3 rgb2xyz = m3x3_invert(xyz2rgb);
    v3 iwp = m3x3_mul_v3(rgb2xyz, (v3){{1,1,1}});
    m3x3 map_white = bradford_map_white(iwp, wp);
    rgb2xyz = m3x3_mul_m3x3(map_white, rgb2xyz);
    m3x3 out = m3x3_invert(rgb2xyz);
    return out;
}

bool
gamut_isvalid(v3 xyz)
{
    if (xyz.v[0] == 0 && xyz.v[1] == 0 && xyz.v[2] == 0) {
        return true;
    }
    if (xyz.v[0] < 0 || xyz.v[1] < 0 || xyz.v[2] < 0) {
        return false;
    }
    if (!v3_isfinite(xyz)) {
        return false;
    }
    v3 xyY = xyz2xyY(xyz);
    assert(glob.locus_poly_count > 0);
    if (!point_inside_polygon(glob.locus_poly_x, glob.locus_poly_y, glob.locus_poly_count, xyY.v[0], xyY.v[1])) {
        return false;
    }
    if (!point_inside_polygon(glob.prophoto_xy[0], glob.prophoto_xy[1], 3, xyY.v[0], xyY.v[1])) {
        return false;
    }
    return true;
}

double
gamut_max_valid_chroma(v3 jch,
                       v3 wp,
                       double clip_hi_limit,
                       void *lcms2_h02)
{
    const m3x3 xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, wp);
    return find_gamut_edge(jch, m3x3_invert(xyz2rgb), &xyz2rgb, glob.locus_poly_count, glob.locus_poly_x, glob.locus_poly_y, clip_hi_limit, 0, lcms2_h02, wp, false);
}

struct extreme_compress_hsv_ge_arg {
    v3 hsv;
    double scale;
    const m3x3 *cam2xyz;
    v3 last_valid_hsv;
};

static double
extreme_compress_hsv_ge(double x, void *arg)
{
    struct extreme_compress_hsv_ge_arg *a = (struct extreme_compress_hsv_ge_arg *)arg;
    v3 hsv = a->hsv;
    hsv.v[1] = x;
    v3 xyz = m3x3_mul_v3(*a->cam2xyz, k_mul_v3(a->scale, hsv2rgb(hsv)));
    if (!gamut_isvalid(xyz)) {
        return 1 + x;
    }
    a->last_valid_hsv = hsv;
    return 1 - x;
}

v3
gamut_clip_to_valid(v3 xyz,
                    v3 wp)
{
    const m3x3 xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, wp);
    v3 rgb = m3x3_mul_v3(xyz2rgb, xyz);
    rgb = v3_clip0(rgb);//gamut_saturated_rgb_clip(rgb);
    xyz = m3x3_mul_v3(m3x3_invert(xyz2rgb), rgb);
    if (!gamut_isvalid(xyz)) {
        // clipping is in blue range, hsv suits best there (based on visual tests)
        v3 hsv = rgb2hsv(v3_norm(rgb));
        m3x3 rgb2xyz = m3x3_invert(xyz2rgb);
        struct extreme_compress_hsv_ge_arg arg = {
            .hsv = hsv,
            .scale = v3_max(rgb),
            .cam2xyz = &rgb2xyz
        };
        interval_halving_min(extreme_compress_hsv_ge, &arg, 0, 1, 1e-6, 100);
        rgb = k_mul_v3(arg.scale, hsv2rgb(arg.last_valid_hsv));
        xyz = m3x3_mul_v3(rgb2xyz, rgb);
    }
    return xyz;
}

void
gamut_init(const struct observer *obs,
           const spectrum_t *illuminant) // should be D50
{
    assert(obs != NULL && illuminant != NULL);

    const double peak[3] = { 0, 1, 0 };
    int start = spec_band_start(obs->cmf[0]);
    int end = spec_band_end(obs->cmf[0]);
    int i;
    for (i = start; i <= end-10; i += 5) {
        spectrum_t *s = spec_alloc(peak, i, i+10, 5);
        v3 xyz = spec2xyz_ill(s, obs, illuminant, 1);
        free(s);
        if (xyz.v[1] >= 0.003) {
            break;
        }
    }
    start = i;
    for (i = end-10; i >= start; i -= 5) {
        spectrum_t *s = spec_alloc(peak, i, i+10, 5);
        v3 xyz = spec2xyz_ill(s, obs, illuminant, 1);
        free(s);
        if (xyz.v[1] >= 0.003) {
            break;
        }
    }
    end = i;
    const int lop_count = 20;
    int poly_count = (end - start) / 5 + 1 + lop_count;
    double *locus_poly_x = malloc(poly_count * sizeof(*locus_poly_x));
    double *locus_poly_y = malloc(poly_count * sizeof(*locus_poly_y));
    for (int i = 0; i < poly_count; i++) {
        spectrum_t *s;
        if (i < poly_count - lop_count) {
            s = spec_alloc(peak, start+i*5, start+i*5+10, 5);
        } else {
            // line of purples
            double freq[6] = { start, start+5, start+10,  end, end+5, end+10 };
            int k = i - (poly_count - lop_count);
            double a = (double)(k+1) / (lop_count+1);
            a = pow(a, 2);
            double dual_peak[6] = { 0,a,0, 0,1-a,0 };
            s = spec_alloc1(freq, dual_peak, 6);
        }
        v3 xyz = spec2xyz_ill(s, obs, illuminant, 1);
        v3 xyY = xyz2xyY(xyz);
        locus_poly_x[i] = xyY.v[0];
        locus_poly_y[i] = xyY.v[1];
        //v3 ltu = xyz2lutspace(xyz);
        //printf("%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
        free(s);
    }
    assert(poly_count > 0);
    glob.locus_poly_x = locus_poly_x;
    glob.locus_poly_y = locus_poly_y;
    glob.locus_poly_count = poly_count;

    {
        const m3x3 rgb2xyz = m3x3_invert(dcp_xyz_D50_to_prophoto_rgb);
        v3 xyY[3];
        xyY[0] = xyz2xyY(m3x3_mul_v3(rgb2xyz, (v3){{1,0,0}}));
        xyY[1] = xyz2xyY(m3x3_mul_v3(rgb2xyz, (v3){{0,1,0}}));
        xyY[2] = xyz2xyY(m3x3_mul_v3(rgb2xyz, (v3){{0,0,1}}));
        for (int i = 0; i < 3; i++) {
            glob.prophoto_xy[0][i] = xyY[i].v[0];
            glob.prophoto_xy[1][i] = xyY[i].v[1];
        }
    }
}

void
gamut_exit(void)
{
    free(glob.locus_poly_x);
    free(glob.locus_poly_y);
}

static inline void
acr_clip_rgb_tone(double *r, double *g, double *b, const double L)
{
	double r_ = *r > L ? L : *r;
	double b_ = *b > L ? L : *b;
	double g_ = b_ + ((r_ - b_) * (*g - *b) / (*r - *b));
	*r = r_;
	*g = g_;
	*b = b_;
}

v3
gamut_acr_rgb_clip(v3 rgb)
{
    double *r = &rgb.v[0];
    double *g = &rgb.v[1];
    double *b = &rgb.v[2];
    if (*r < 0) *r = 0;
    if (*g < 0) *g = 0;
    if (*b < 0) *b = 0;

    // This is Adobe's hue-stable film-like curve with a diagonal, ie only used for clipping. Can probably be further optimized.
    const double L = 1.0;
    if (*r >= *g) {
        if (*g > *b) {         // Case 1: r >= g >  b
            acr_clip_rgb_tone(r, g, b, L);
        } else if (*b > *r) {  // Case 2: b >  r >= g
            acr_clip_rgb_tone(b, r, g, L);
        } else if (*b > *g) {  // Case 3: r >= b >  g
            acr_clip_rgb_tone(r, b, g, L);
        } else {               // Case 4: r >= g == b
            *r = *r > L ? L : *r;
            *g = *g > L ? L : *g;
            *b = *g;
        }
    } else {
        if (*r >= *b) {        // Case 5: g >  r >= b
            acr_clip_rgb_tone(g, r, b, L);
        } else if (*b > *g) {  // Case 6: b >  g >  r
            acr_clip_rgb_tone(b, g, r, L);
        } else {               // Case 7: g >= b >  r
            acr_clip_rgb_tone(g, b, r, L);
        }
    }
    return rgb;
}

v3
gamut_saturated_rgb_clip(v3 rgb)
{
    int one_count = 0;
    int zero_count = 0;
    for (int i = 0; i < 3; i++) {
        if (rgb.v[i] < 0) { rgb.v[i] = 0; zero_count++; }
        if (rgb.v[i] > 1) one_count++;
    }
    if (one_count == 3) {
        return (v3){{ 1,1,1 }};
    }
    if (zero_count == 3) {
        return (v3){{ 0,0,0 }};
    }
    if (one_count == 0) {
        return rgb;
    }

    const double hue = rgb2hsl(v3_norm(rgb)).v[0] * 360;

    v3 argb = gamut_acr_rgb_clip(rgb);

    if (hue >= 180 && hue <= 270) {
        // In the blue range normal ACR clip has hue shift issues (blue becomes purple), so we do a color correction

        const m3x3 rgb2xyz = m3x3_invert(xyz2rgb_srgb);
        const v3 wp = wp_from_rgb2xyz(rgb2xyz);
        v3 xyz = m3x3_mul_v3(rgb2xyz, v3_norm(rgb));
        v3 jch_o = xyz2jch(xyz, wp);
        xyz = m3x3_mul_v3(rgb2xyz, argb);
        v3 jch_a = xyz2jch(xyz, wp);
        jch_o.v[0] = jch_a.v[0];
        jch_o.v[1] = jch_a.v[1];
        xyz = jch2xyz(jch_o, wp);
        v3 nrgb = m3x3_mul_v3(xyz2rgb_srgb, xyz);
        nrgb = v3_clip0(nrgb); // should be minimal clipping here, so we can ignore the hue shift this clip causes
        nrgb = v3_norm(nrgb);

        // blend to normal ACR clip in edges
        if (hue < 200) {
            double x = (hue - 180) / 20;
            argb = (v3){{ nrgb.v[0] * (1-x) + argb.v[0] * x, nrgb.v[1] * (1-x) + argb.v[1] * x, nrgb.v[2] * (1-x) + argb.v[2] * x }};
        } else if (hue > 250) {
            double x = (hue - 250) / 20;
            argb = (v3){{ nrgb.v[0] * (1-x) + argb.v[0] * x, nrgb.v[1] * (1-x) + argb.v[1] * x, nrgb.v[2] * (1-x) + argb.v[2] * x }};
        } else {
            argb = nrgb;
        }
    }

    if (hue >= 240 && hue <= 360) {
        // in this range blue and red turns magenta (strong effect), then we rather use desaturated acr clip
        // however, acr clip (hsl-style) still suffers from a slight magenta cast int the blue-to-white
        // transition, so we go via CIECAM02
        return argb;
    }

    v3 crgb = v3_clip01(rgb);
    if (hue >= 180 && hue <= 240) {
        // in this range blue turns cyan (strong effect), but as blues are really flat we keep some of the cyan effect for better tonality
        double x = 0.25;
        if (hue > 180) {
            x = 0.25 + 0.75 * (hue - 180) / 60;
        }
        return (v3){{ crgb.v[0] * (1-x) + argb.v[0] * x, crgb.v[1] * (1-x) + argb.v[1] * x, crgb.v[2] * (1-x) + argb.v[2] * x }};
    }
    // in this range red and green turns yellow (weak effect), we accept that as it's more important to keep saturation, this is clearly seen in sunset type of scenes
    return (v3){{ crgb.v[0] * 0.9 + argb.v[0] * 0.1, crgb.v[1] * 0.9 + argb.v[1] * 0.1, crgb.v[2] * 0.9 + argb.v[2] * 0.1 }};
}
