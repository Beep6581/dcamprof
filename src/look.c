/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <lcms2.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <bisection.h>
#include <nmsimplex.h>
#include <dngref.h>
#include <interp.h>
#include <gamut.h>
#include <elog.h>
#include <look.h>

struct look_tone_rep_op_t_ {
    cmsHANDLE **h02;
    int h02_count;
    v3 wp;
    const double *tc;
    int tc_len;
    m3x3 rgb2xyz;
    m3x3 xyz2rgb;
    struct tone_rep_op_config conf;
    enum tc_type tc_type;
};

struct find_tc_slope_fun_arg {
    const double *tc;
    int tc_len;
};

static double
find_tc_slope_fun(double k,
                  void *arg)
{
    struct find_tc_slope_fun_arg *a = (struct find_tc_slope_fun_arg *)arg;
    double areasum = 0;
    const int steps = 256;
    for (int i = 0; i < steps; i++) {
        double x = 0.1 + ((double)i / (steps-1)) * 0.5; // testing (sRGB) range [0.1 - 0.6], ie ignore highligths and dark shadows
        double y = srgb_gamma_forward(curve_y(srgb_gamma_inverse(x), a->tc, a->tc_len));
        double y1 = k * x;
        if (y1 > 1) y1 = 1;
        areasum += (y - y1) * (y - y1); // square is a rough approx of (twice) the area, but it's fine for our purposes
    }
    return areasum;
}

// calculate a single value that represents the contrast of the tone curve
static double
calculate_tonecurve_contrast_value(const double tc[],
                                   int tc_len)
{

    // find linear y = k*x the best approximates the curve, which is the linear scaling/exposure component that does not contribute any contrast

    // Note: the analysis is made on the gamma encoded curve, as the LUT is linear we make backwards gamma to
    struct find_tc_slope_fun_arg arg = { tc, tc_len };
    double k = interval_halving_min(find_tc_slope_fun, &arg, 0.1, 5.0, 0.001, 10000);
    //elog("average slope: %f\n", k);

    double maxslope = 0;
    {
        // look at midtone slope
        const double xd = 0.07;
        const double tx[] = { 0.30, 0.35, 0.40, 0.45 }; // we only look in the midtone range
        for (int i = 0; i < (int)(sizeof(tx)/sizeof(tx[0])); i++) {
            double x0 = tx[i] - xd;
            double y0 = srgb_gamma_forward(curve_y(srgb_gamma_inverse(x0), tc, tc_len)) - k * x0;
            double x1 = tx[i] + xd;
            double y1 = srgb_gamma_forward(curve_y(srgb_gamma_inverse(x1), tc, tc_len)) - k * x1;
            double slope = 1.0 + (y1 - y0) / (x1 - x0);
            if (slope > maxslope) {
                maxslope = slope;
            }
        }

        // look at slope at (light) shadows and (dark) highlights
        double e_maxslope = 0;
        {
            const double tx[] = { 0.20, 0.25, 0.50, 0.55 }; // we only look in the midtone range
            for (int i = 0; i < (int)(sizeof(tx)/sizeof(tx[0])); i++) {
                double x0 = tx[i] - xd;
                double y0 = srgb_gamma_forward(curve_y(srgb_gamma_inverse(x0), tc, tc_len)) - k * x0;
                double x1 = tx[i] + xd;
                double y1 = srgb_gamma_forward(curve_y(srgb_gamma_inverse(x1), tc, tc_len)) - k * x1;
                double slope = 1.0 + (y1 - y0) / (x1 - x0);
                if (slope > e_maxslope) {
                    e_maxslope = slope;
                }
            }
        }
        //elog("%.3f %.3f\n", maxslope, e_maxslope);
        // midtone slope is more important for contrast, but weigh in some slope from brights and darks too.
        maxslope = maxslope * 0.7 + e_maxslope * 0.3;
    }
    return maxslope;
}

void
look_neutral_tone_rep_op_tc_analysis(const double tc[],
                                     int tc_len,
                                     double *cv,
                                     double *cs)
{
    const double p[] = {
        0.60, 0.70, // lowest contrast
        0.70, 0.80,
        0.90, 0.94,
        0.99, 1.00,
        1.00, 1.00, // 1.0 (linear curve) to 1.0, no scaling
        1.07, 1.00,
        1.08, 1.00,
        1.11, 1.02,
        1.20, 1.08,
        1.30, 1.12,
        1.80, 1.20,
        2.00, 1.22  // highest contrast
    };

    double cf[1000];
    double cf_range[2];
    const size_t in_len = sizeof(p)/sizeof(p[0])/2;
    double in_x[in_len];
    double in_y[in_len];
    for (size_t i = 0; i < in_len; i++) {
        in_x[i]= p[2*i+0];
        in_y[i]= p[2*i+1];
    }
    const size_t out_len = sizeof(cf)/sizeof(cf[0]);
    double out_x[out_len];
    for (size_t i = 0; i < out_len; i++) {
        out_x[i] = in_x[0] + (in_x[in_len-1] - in_x[0]) * (double)i / (out_len-1);
    }
    cubic_spline(in_x, in_y, in_len, out_x, cf, out_len);
    cf_range[0] = in_x[0];
    cf_range[1] = in_x[in_len-1];

    *cv = calculate_tonecurve_contrast_value(tc, tc_len);
    *cs = curve_y_rng(*cv, cf_range, cf, sizeof(cf)/sizeof(cf[0]));
}

v3
look_wp(void)
{
    cmsCIEXYZ XYZ = *cmsD50_XYZ();
    return (v3){{ XYZ.X, XYZ.Y, XYZ.Z }};
}

look_tone_rep_op_t *
look_tone_rep_op_new(const m3x3 rgb2xyz,
                     const double tc[],
                     int tc_len,
                     enum tc_type tc_type,
                     enum gc_type gc_type,
                     const struct tone_rep_op_config *conf)
{
    look_tone_rep_op_t *state = calloc(1, sizeof(*state));
    cmsViewingConditions vc;
    state->wp = look_wp();
    vc.whitePoint.X = state->wp.v[0] * 100;
    vc.whitePoint.Y = state->wp.v[1] * 100;
    vc.whitePoint.Z = state->wp.v[2] * 100;
    vc.Yb = 20;
    vc.La = 20;
    vc.surround = AVG_SURROUND;
    vc.D_value = 1.0;

#ifdef _OPENMP
    state->h02_count = omp_get_max_threads();
    if (state->h02_count == 0) {
        state->h02_count = 1;
    }
#else
    state->h02_count = 1;
#endif
    state->h02 = malloc(state->h02_count * sizeof(state->h02[0]));
    for (int i = 0; i < state->h02_count; i++) {
        state->h02[i] = cmsCIECAM02Init(0, &vc);
    }

    if (conf != NULL) {
        state->conf = *conf;
    } else {
        state->conf.chroma_scaling = -1;

        state->conf.curve.keep_factor = 0.0;
        state->conf.curve.low_chroma = 25;
        state->conf.curve.high_chroma = 65;

        state->conf.saturated.adjust_factor = 0.92;
        state->conf.saturated.low_chroma = 35;
        state->conf.saturated.high_chroma = 60;

        state->conf.shadows.adjust_factor = 1.20;
        state->conf.shadows.low_lightness = 0.15;
        state->conf.shadows.high_lightness = 0.50;
        state->conf.shadows.adjust_factor_high_chroma = 1.08;
        state->conf.shadows.low_chroma = 25;
        state->conf.shadows.high_chroma = 40;

        state->conf.rolloff.keep_factor = calloc(1, sizeof(struct lookop_curve) + sizeof(double) * 4);
        state->conf.rolloff.keep_factor->type = LOOKOP_CURVE_LINEAR;
        state->conf.rolloff.keep_factor->handle_count = 2;
        state->conf.rolloff.keep_factor->handles[0] = 0;
        state->conf.rolloff.keep_factor->handles[1] = 0.20;
        state->conf.rolloff.keep_factor->handles[2] = 360;
        state->conf.rolloff.keep_factor->handles[3] = 0.20;
        state->conf.rolloff.low_satscale = calloc(1, sizeof(struct lookop_curve) + sizeof(double) * 4);
        state->conf.rolloff.low_satscale->type = LOOKOP_CURVE_LINEAR;
        state->conf.rolloff.low_satscale->handle_count = 2;
        state->conf.rolloff.low_satscale->handles[0] = 0;
        state->conf.rolloff.low_satscale->handles[1] = 1.00;
        state->conf.rolloff.low_satscale->handles[2] = 360;
        state->conf.rolloff.low_satscale->handles[3] = 1.00;
        state->conf.rolloff.high_satscale = calloc(1, sizeof(struct lookop_curve) + sizeof(double) * 4);
        state->conf.rolloff.high_satscale->type = LOOKOP_CURVE_LINEAR;
        state->conf.rolloff.high_satscale->handle_count = 2;
        state->conf.rolloff.high_satscale->handles[0] = 0;
        state->conf.rolloff.high_satscale->handles[1] = 1.20;
        state->conf.rolloff.high_satscale->handles[2] = 360;
        state->conf.rolloff.high_satscale->handles[3] = 1.20;
    }

    switch (gc_type) {
    case GC_FROM_FILE:
        break;
    case GC_NONE:
        state->conf.gamut_enabled = false;
        break;
    case GC_SRGB:
    case GC_SRGB_STRONG:
        memset(&state->conf.gamut, 0, sizeof(state->conf.gamut));
        state->conf.gamut_enabled = true;

        state->conf.gamut.inner.xyz2rgb_enabled = true;
        state->conf.gamut.inner.xyz2rgb = gamut_remap_xyz2rgb(xyz2rgb_srgb, state->wp);
        state->conf.gamut.inner.chroma_scale = 0.6;

        state->conf.gamut.dest.xyz2rgb_enabled = true;
        state->conf.gamut.dest.xyz2rgb = gamut_remap_xyz2rgb(xyz2rgb_srgb, state->wp);
        state->conf.gamut.dest.chroma_scale = 1.0;

        state->conf.gamut.clip.locus_enabled = true;
        state->conf.gamut.clip.xyz2rgb_enabled = true;
        state->conf.gamut.clip.xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, state->wp);
        state->conf.gamut.clip.chroma_scale = 1.0;

        state->conf.gamut.lo_rgb_lim = -0.005;
        state->conf.gamut.hi_rgb_lim = 1000000000;

        if (gc_type == GC_SRGB_STRONG) {
            state->conf.gamut.hsv.enabled = true;

            state->conf.gamut.hsv.lo_sat_blend_lim = 0.6;
            state->conf.gamut.hsv.hi_sat_blend_lim = 0.9;

            state->conf.gamut.hsv.top.linear_limit = 1.2;
            state->conf.gamut.hsv.top.cutoff_limit = 2.0;
            state->conf.gamut.hsv.top.output_limit = 1.4;
            state->conf.gamut.hsv.top.min_lift = 0.0;

            state->conf.gamut.hsv.bottom.linear_limit = -0.05;
            state->conf.gamut.hsv.bottom.cutoff_limit = -0.2;
            state->conf.gamut.hsv.bottom.output_limit = -0.1;
            state->conf.gamut.hsv.bottom.max_lift = 0.3;
        }
        break;

    case GC_ADOBERGB:
    case GC_ADOBERGB_STRONG:
        memset(&state->conf.gamut, 0, sizeof(state->conf.gamut));
        state->conf.gamut_enabled = true;

        state->conf.gamut.inner.xyz2rgb_enabled = true;
        state->conf.gamut.inner.xyz2rgb = gamut_remap_xyz2rgb(xyz2rgb_adobergb, state->wp);
        state->conf.gamut.inner.chroma_scale = 0.6;

        state->conf.gamut.dest.xyz2rgb_enabled = true;
        state->conf.gamut.dest.xyz2rgb = gamut_remap_xyz2rgb(xyz2rgb_adobergb, state->wp);
        state->conf.gamut.dest.chroma_scale = 1.0;

        state->conf.gamut.clip.locus_enabled = true;
        state->conf.gamut.clip.xyz2rgb_enabled = true;
        state->conf.gamut.clip.xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, state->wp);
        state->conf.gamut.clip.chroma_scale = 1.0;

        state->conf.gamut.lo_rgb_lim = -0.005;
        state->conf.gamut.hi_rgb_lim = 1000000000;

        if (gc_type == GC_ADOBERGB_STRONG) {
            state->conf.gamut.hsv.enabled = true;

            state->conf.gamut.hsv.lo_sat_blend_lim = 0.6;
            state->conf.gamut.hsv.hi_sat_blend_lim = 0.9;

            state->conf.gamut.hsv.top.linear_limit = 1000;
            state->conf.gamut.hsv.top.cutoff_limit = 1000;
            state->conf.gamut.hsv.top.output_limit = 1000;
            state->conf.gamut.hsv.top.min_lift = 0.0;

            state->conf.gamut.hsv.bottom.linear_limit = -0.05;
            state->conf.gamut.hsv.bottom.cutoff_limit = -0.2;
            state->conf.gamut.hsv.bottom.output_limit = -0.1;
            state->conf.gamut.hsv.bottom.max_lift = 0.3;
        }
        break;
    }

    if (state->conf.chroma_scaling < 0) {
        double cv, cs;
        look_neutral_tone_rep_op_tc_analysis(tc, tc_len, &cv, &cs);
        state->conf.chroma_scaling = cs;
    }

    state->tc = tc;
    state->tc_len = tc_len;
    state->tc_type = tc_type;
    state->rgb2xyz = rgb2xyz;
    state->xyz2rgb = m3x3_invert(rgb2xyz);

    gamut_compress_state_alloc(&state->conf.gamut, rgb2xyz, state->wp);
    return state;
}

void
look_tone_rep_op_delete(look_tone_rep_op_t *state)
{
    if (state == NULL) {
        return;
    }
    for (int i = 0; i < state->h02_count; i++) {
        cmsCIECAM02Done(state->h02[i]);
    }
    gamut_compress_state_free(&state->conf.gamut);
    free(state->h02);
    free(state);
}

static double
lookop_curve(const struct lookop_curve *c,
             double x)
{
    double in_x[c->handle_count];
    double in_y[c->handle_count];
    for (int i = 0; i < c->handle_count; i++) {
        in_x[i]= c->handles[2*i+0];
        in_y[i]= c->handles[2*i+1];
    }
    double y;
    switch (c->type) {
    case LOOKOP_CURVE_SPLINE:
        cubic_spline(in_x, in_y, c->handle_count, &x, &y, 1);
        break;
    case LOOKOP_CURVE_ROUNDEDSTEP:
        rounded_linear_interpolation(in_x, in_y, c->handle_count, &x, &y, 1);
        break;
    case LOOKOP_CURVE_LINEAR:
        linear_interpolation(in_x, in_y, c->handle_count, &x, &y, 1);
        break;
    }
    return y;
}

double
look_lookop_curve_y(const struct lookop_curve *c,
                    double x)
{
    return lookop_curve(c, x);
}

static double
get_vs(v3 jch,
       v3 rgb,
       int vs_count,
       const struct lookop_vs *vs_array)
{
    double s = 1.0;

    for (int i = 0; i < vs_count; i++) {
        const struct lookop_vs *vs = &vs_array[i];
        double x = 0;
        double diff = vs->xrange[1] - vs->xrange[0];
        double start = vs->xrange[0];
        switch (vs->x_type) {
        case LOOKOP_VS_X_LIGHTNESS: x = jch.v[0]; break;
        case LOOKOP_VS_X_CHROMA: x = jch.v[1]; break;
        case LOOKOP_VS_X_HUE:
            x = jch.v[2];
            if (start < 0) start += 360;
            if (start + diff > 360 && x < vs->xrange[1]) x += 360;
            break;
        case LOOKOP_VS_X_HSL_LIGHTNESS:
        case LOOKOP_VS_X_HSL_SATURATION:
        case LOOKOP_VS_X_HSL_HUE:
        case LOOKOP_VS_X_HSV_VALUE:
        case LOOKOP_VS_X_HSV_SATURATION:
        case LOOKOP_VS_X_HSV_HUE: {
            v3 hsl = rgb2hsl(rgb);
            v3 hsv = rgb2hsv(rgb);
            switch (vs->x_type) {
            default: abort();
            case LOOKOP_VS_X_HSL_LIGHTNESS: x = 100 * hsl.v[2]; break;
            case LOOKOP_VS_X_HSL_SATURATION: x = 100 * hsl.v[1]; break;
            case LOOKOP_VS_X_HSL_HUE:
                x = hsl.v[0] * 360.0;
                if (start < 0) start += 360;
                if (start + diff > 360 && x < vs->xrange[1]) x += 360;
                break;
            case LOOKOP_VS_X_HSV_VALUE: x = 100 * hsv.v[2]; break;
            case LOOKOP_VS_X_HSV_SATURATION: x = 100 * hsv.v[1]; break;
            case LOOKOP_VS_X_HSV_HUE:
                x = hsv.v[0] * 360.0;
                if (start < 0) start += 360;
                if (start + diff > 360 && x < vs->xrange[1]) x += 360;
                break;
            }
            break;
        }
        }
        x -= start;
        x /= diff;
        if (x < 0) x = 0;
        if (x > 1) x = 1;
        if (!isfinite(x)) x = 0;
        double y = lookop_curve(vs->curve, x);
        y *= vs->yscale;
        s *= y;
    }

    return s;
}

static void
apply_look_operators(v3 jch,
                     v3 *out,
                     look_tone_rep_op_t *state,
                     cmsHANDLE *h02)
{
    if (state->conf.lookop_count <= 0) {
        *out = jch;
        return;
    }

    for (int i = 0; i < state->conf.lookop_count; i++) {
        const struct lookop *lo = &state->conf.lookops[i];

        v3 rgb, xyz;
        {
            xyz = jch2xyz_s(jch, h02);
            rgb = m3x3_mul_v3(state->xyz2rgb, xyz);
            bool changed = false;
            if (v3_isclip01(rgb)) {
                rgb = gamut_saturated_rgb_clip(rgb);
                xyz = m3x3_mul_v3(state->rgb2xyz, rgb);
                changed = true;
            }
            if (!gamut_isvalid(xyz)) {
                xyz = gamut_clip_to_valid(xyz, state->wp);
                rgb = gamut_saturated_rgb_clip(rgb);
                changed = true;
            }
            if (changed) {
                jch = xyz2jch_s(xyz, h02);
            }
        }
        double vs = get_vs(jch, rgb, lo->vs_count, (const struct lookop_vs *)lo->vs);
        if (lo->blend_invert) vs = 1.0 - vs;

        v3 new_jch = jch;
        switch (lo->type) {
        case LOOKOP_STRETCH:
            for (int j = 0; j < lo->p.stretch.dim_count; j++) {
                double x = 0;
                switch (lo->p.stretch.dims[j]) {
                case LOOKOP_VS_X_LIGHTNESS: x = jch.v[0]; break;
                case LOOKOP_VS_X_CHROMA: x = jch.v[1]; break;
                case LOOKOP_VS_X_HUE: x = jch.v[2]; break;
                default: abort();
                }
                x -= lo->p.stretch.ranges[j][0];
                x /= lo->p.stretch.ranges[j][1] - lo->p.stretch.ranges[j][0];
                if (x < 0 || x > 1 || !isfinite(x)) {
                    continue;
                }
                double y = lookop_curve(lo->p.stretch.curves[j], x);
                y *= lo->p.stretch.ranges[j][1] - lo->p.stretch.ranges[j][0];
                y += lo->p.stretch.ranges[j][0];
                switch (lo->p.stretch.dims[j]) {
                case LOOKOP_VS_X_LIGHTNESS: new_jch.v[0] = y; break;
                case LOOKOP_VS_X_CHROMA: new_jch.v[1] = y; break;
                case LOOKOP_VS_X_HUE: new_jch.v[2] = y; break;
                default: abort();
                }
            }
            break;
        case LOOKOP_ADDHUE:
            new_jch.v[2] += lo->p.gen.value;
            if (new_jch.v[2] > 360) new_jch.v[2] -= 360;
            if (new_jch.v[2] < 0)   new_jch.v[2] += 360;
            break;
        case LOOKOP_ADDCHROMA:
            new_jch.v[1] += lo->p.gen.value;
            break;
        case LOOKOP_ADDLIGHTNESS:
            new_jch.v[0] += lo->p.gen.value;
            break;
        case LOOKOP_SCALECHROMA:
            new_jch.v[1] *= lo->p.gen.value;
            break;
        case LOOKOP_SCALELIGHTNESS:
            new_jch.v[0] *= lo->p.gen.value;
            break;
        case LOOKOP_SETTEMPERATURE: {
            v3 wp = dcp_d50();
            wp = xyY2xyz((v3){{ lo->p.settemp.xy[0], lo->p.settemp.xy[1], wp.v[1] }});
            v3 new_xyz = chromatic_adaptation_transform(CAT_CAT02, xyz, wp, dcp_d50());
            if (!gamut_isvalid(new_xyz)) {
                new_xyz = gamut_clip_to_valid(xyz, state->wp);
            }
            new_jch = xyz2jch_s(new_xyz, h02);
            break;
        }
        case LOOKOP_CURVES: {
            v3 new_rgb;
            for (int i = 0; i < 3; i++) {
                double x = rgb.v[i];
                if (lo->p.curves.is_srgb_gamma) {
                    x = srgb_gamma_forward(x);
                } else if (lo->p.curves.gamma != 1.0) {
                    x = pow(x, 1.0 / lo->p.curves.gamma);
                }
                double y = lookop_curve(lo->p.curves.curves[i], x);
                if (lo->p.curves.is_srgb_gamma) {
                    y = srgb_gamma_inverse(y);
                } else if (lo->p.curves.gamma != 1.0) {
                    y = pow(y, lo->p.curves.gamma);
                }
                new_rgb.v[i] = y;
            }
            v3 new_xyz = m3x3_mul_v3(state->rgb2xyz, new_rgb);
            if (!gamut_isvalid(new_xyz)) {
                new_xyz = gamut_clip_to_valid(xyz, state->wp);
            }
            new_jch = xyz2jch_s(new_xyz, h02);
            if (lo->p.curves.keep_lightness) {
                new_jch.v[0] = jch.v[0];
            }
            break;
        }
        }

        if (lo->blend_rgb) {
            v3 new_rgb = gamut_saturated_rgb_clip(jch2rgb_s(new_jch, state->xyz2rgb, h02));

            rgb.v[0] = srgb_gamma_forward(rgb.v[0]) * (1 - vs) + srgb_gamma_forward(new_rgb.v[0]) * vs;
            rgb.v[1] = srgb_gamma_forward(rgb.v[1]) * (1 - vs) + srgb_gamma_forward(new_rgb.v[1]) * vs;
            rgb.v[2] = srgb_gamma_forward(rgb.v[2]) * (1 - vs) + srgb_gamma_forward(new_rgb.v[2]) * vs;
            for (int i = 0; i < 3; i++) {
                rgb.v[i] = srgb_gamma_inverse(rgb.v[i]);
            }
            xyz = m3x3_mul_v3(state->rgb2xyz, rgb);
            if (!gamut_isvalid(xyz)) {
                xyz = gamut_clip_to_valid(xyz, state->wp);
            }
            jch = xyz2jch_s(xyz, h02);
        } else {
            jch = jch_blend(new_jch, jch, vs);
            xyz = jch2xyz_s(jch, h02);
            if (!gamut_isvalid(xyz)) {
                xyz = gamut_clip_to_valid(xyz, state->wp);
                jch = xyz2jch_s(xyz, h02);
            }
        }

    }
    *out = jch;
}

void
look_apply_target_adjustment(struct patch_set *ps,
                             const struct target_adjustment_config *conf,
                             v3 wp)
{
    for (int i = 0; i < ps->patch_count; i++) {
        v3 jch = xyz2jch(ps->patch[i].xyz, wp);
        v3 orig_jch = jch;

        // global adjustments
        jch.v[0] *= conf->global.scale_lightness;
        jch.v[1] *= conf->global.scale_chroma;
        if (conf->global.scale_lightness_per_hue != NULL) jch.v[0] *= lookop_curve(conf->global.scale_lightness_per_hue, jch.v[2]);
        if (conf->global.scale_chroma_per_hue != NULL)    jch.v[1] *= lookop_curve(conf->global.scale_chroma_per_hue, jch.v[2]);
        if (conf->global.add_hue_per_hue != NULL)         jch.v[2] += lookop_curve(conf->global.add_hue_per_hue, jch.v[2]);

        // per patch adjustments
        for (int k = 0; k < conf->patch_count; k++) {
            if (strcmp(conf->patch[k].name, ps->patch[i].str_id) == 0 && (conf->patch[k].class[0] == '\0' || strcmp(conf->patch[k].class, ps->class_names[ps->patch[i].class]) == 0 || ps->class_names[ps->patch[i].class][0] == '\0')) {
                jch.v[0] *= conf->patch[k].adjust_jch[0];
                jch.v[1] *= conf->patch[k].adjust_jch[1];
                jch.v[2] += conf->patch[k].adjust_jch[2];
                if (conf->patch[k].scale_rgb[0] != 1 || conf->patch[k].scale_rgb[1] != 1 || conf->patch[k].scale_rgb[2] != 1) {
                    v3 xyz = jch2xyz(jch, wp);
                    m3x3 xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, wp);
                    v3 rgb = v3_clip01(m3x3_mul_v3(xyz2rgb, xyz));
                    rgb.v[0] *= conf->patch[k].scale_rgb[0];
                    rgb.v[1] *= conf->patch[k].scale_rgb[1];
                    rgb.v[2] *= conf->patch[k].scale_rgb[2];
                    xyz = m3x3_mul_v3(m3x3_invert(xyz2rgb), rgb);
                    v3 jch1 = xyz2jch(xyz, wp);
                    jch.v[1] = jch1.v[1];
                    jch.v[2] = jch1.v[2];
                }
            }
        }

        char cname[TARGET_CLASS_NAME_SIZE + 16];
        if (ps->patch[i].class == -1 || ps->class_names[ps->patch[i].class][0] == '\0') {
            strcpy(cname, "");
        } else {
            sprintf(cname, " (class \"%s\")", ps->class_names[ps->patch[i].class]);
        }
        if (!v3_isequal(jch, orig_jch)) {
            elog("Patch \"%s\"%s reference value was adjusted from JCh %5.1f %5.1f %5.1f to %5.1f %5.1f %5.1f\n",
                 ps->patch[i].str_id, cname, orig_jch.v[0], orig_jch.v[1], orig_jch.v[2], jch.v[0], jch.v[1], jch.v[2]);
            ps->patch[i].xyz = jch2xyz(jch, wp);
        } else {
            elog("Patch \"%s\"%s didn't match any adjustment, kept JCh %5.1f %5.1f %5.1f\n",
                 ps->patch[i].str_id, cname, orig_jch.v[0], orig_jch.v[1], orig_jch.v[2]);
        }
    }
}

void *
look_tone_rep_op_get_lcmsh02(look_tone_rep_op_t *state)
{
    cmsHANDLE *h02;
#ifdef _OPENMP
    {
        int tidx = omp_get_thread_num();
        assert(tidx >= 0 && tidx < state->h02_count);
        h02 = state->h02[tidx];
    }
#else
    h02 = state->h02[0];
#endif
    return (void *)h02;
}

void
look_tone_rep_op(v3 *out,
                 v3 xyz, // D50!
                 look_tone_rep_op_t *state)
{
    const struct tone_rep_op_config *conf = &state->conf;
    const v3 orig_xyz = xyz;

    *out = (v3){{0,0,0}};
    if (xyz.v[0] == 0.0 && xyz.v[1] == 0.0 && xyz.v[2] == 0.0) {
        return;
    }

    cmsHANDLE *h02 = (cmsHANDLE *)look_tone_rep_op_get_lcmsh02(state);

    if (conf->gamut_enabled) {
        // Gamut compression. This is done on the linear "colorimetric" color, meaning that the following
        // tone reproduction can push it outside the gamut again, so it's not "exact". However, as the curve
        // affects saturation in non-linear ways (rolloff towards whitepoint etc) it still makes more sense
        // to apply gamut compression at this stage.
        v3 compressed_xyz = gamut_compress(&conf->gamut, state->rgb2xyz, state->xyz2rgb, xyz, state->wp, h02);
        if (!v3_isequal(compressed_xyz, xyz)) {
            if (!v3_isfinite(compressed_xyz)) {
                elog("Bug: gamut compression returned non-finite value.\n");
                v3_print("xyz", xyz);
                v3_print("compressed_xyz", compressed_xyz);
                abort();
            }
            xyz = compressed_xyz;
            if (xyz.v[0] == 0.0 && xyz.v[1] == 0.0 && xyz.v[2] == 0.0) {
                return;
            }
        }
    }

    v3 orig_rgb = m3x3_mul_v3(state->xyz2rgb, xyz);
    orig_rgb = gamut_saturated_rgb_clip(orig_rgb); // small/no clipping expected due to prior gamut compression in the coloriemtric profile

    bool is_linear = false;
    if (state->tc == NULL || state->tc_len == 2) {
        // special case, linear curve
        is_linear = true;
    }

    const v3 prgb = {{ curve_y(orig_rgb.v[0], state->tc, state->tc_len),
                       curve_y(orig_rgb.v[1], state->tc, state->tc_len),
                       curve_y(orig_rgb.v[2], state->tc, state->tc_len) }};
    const v3 pxyz = m3x3_mul_v3(state->rgb2xyz, prgb);

    const v3 argb = dcp_apply_tonecurve(orig_rgb, state->tc, state->tc_len);
    const v3 axyz = m3x3_mul_v3(state->rgb2xyz, argb);

    if (state->tc_type == TC_RGB
        || state->tc_type == TC_SIMPLE
        || state->tc_type == TC_SIMPLE_ACR_HUE
        || state->tc_type == TC_SIMPLE_RGB_HUE) {

        const v3 hsv = dcp_rgb2hsv(orig_rgb);
        const v3 prgb_hsv = dcp_rgb2hsv(prgb);
        const v3 argb_hsv = dcp_rgb2hsv(argb);
        // orig H == argb H
        assert(fabs(hsv.v[0] - argb_hsv.v[0]) < 1e-5 || fabs(hsv.v[0] - argb_hsv.v[0] - 6) < 1e-5);
        // prgb V == argb V (max channel)
        assert(fabs(argb_hsv.v[2] / prgb_hsv.v[2] - 1) < 1e-5);
        // prgb S == argb S (max - min)
        assert(fabs(prgb_hsv.v[1] / argb_hsv.v[1] - 1) < 1e-5);

        v3 ref_rgb;

        if (state->tc_type == TC_RGB) {

            ref_rgb = prgb;

        } else if (state->tc_type == TC_SIMPLE
                || state->tc_type == TC_SIMPLE_ACR_HUE
                || state->tc_type == TC_SIMPLE_RGB_HUE) {

            // pure rgb by default
            v3 rgb1 = prgb;

            if (state->tc_type == TC_SIMPLE) {

                if (hsv.v[0] > 0 && hsv.v[0] < 6 * 60 / 360) {
                    // hue 0...60 - red-yellow

                    // skintone CCSG data (-1EV, ProPhoto linear):
                    // H 21...37 deg
                    // S 0.42...0.84 (0.42...0.72 excluding dark skintone < -3EV)
                    // V -4.1...-1.7EV
                    // This is reference data, keep in mind highlights and shadows.

                    // keep argb hue for value from 0 to -2EV with roll-off to 0EV
                    // (for reference: dng curve max correction at x = 0.178 (-2.5EV)
                    //  by 1.1EV, from 50L to 68L)
                    const double v_low = pow(2.0, -2.0);  // dng_curve(2**-2) == 0.52 (77L)
                    const double v_high = pow(2.0, 0); // dng_curve(2**-0.5) == 0.92 (97L)
                    // weight of pure rgb
                    double wv = (hsv.v[2] - v_low) / (v_high - v_low);
                    wv = (wv < 0) ? 0 : wv;
                    wv = (wv > 1) ? 1 : wv;

                    // pure rgb for high saturated AND bright colors
                    // tested on sunny images with hair highlights
                    const double sun_s_low = 0.65;
                    const double sun_s_high = 0.80;
                    // weight of pure rgb
                    double w_sun_s = (hsv.v[1] - sun_s_low) / (sun_s_high - sun_s_low);
                    w_sun_s = (w_sun_s < 0) ? 0 : w_sun_s;
                    w_sun_s = (w_sun_s > 1) ? 1 : w_sun_s;

                    const double sun_v_low = pow(2.0, -3);
                    const double sun_v_high = pow(2.0, -1);
                    double w_sun_v = (hsv.v[2] - sun_v_low) / (sun_v_high - sun_v_low);
                    w_sun_v = (w_sun_v < 0) ? 0 : w_sun_v;
                    w_sun_v = (w_sun_v > 1) ? 1 : w_sun_v;

                    // sun is high sat AND high value
                    const double w_sun = w_sun_s * w_sun_v;

                    // weight for pure RGB
                    // double w = wv;
                    // pure rgb for high value OR sun
                    double w = 1 - (1 - wv) * (1 - w_sun);

                    rgb1.v[0] = (1 - w) * argb.v[0] + w * prgb.v[0];
                    rgb1.v[1] = (1 - w) * argb.v[1] + w * prgb.v[1];
                    rgb1.v[2] = (1 - w) * argb.v[2] + w * prgb.v[2];
                }
                else if (hsv.v[0] > 6 * 240 / 360 && hsv.v[0] < 6 * 300 / 360) {
                    // prevent blue shift to purple
                    rgb1 = argb;
                } else {
                    // // acr rgb hue to pure rgb hue roll-off from 0 to -3ev, above -3ev pure rgb hue
                    // const double v_low = 0;
                    // const double v_high = pow(2.0, -3.0);

                    // acr rgb hue to pure rgb hue roll-off from -2ev to 0ev
                    // tested on bright green leafs, that must not become yellow
                    const double v_low = pow(2.0, -2.0);
                    const double v_high = pow(2.0, 0.0);
                    double w = (hsv.v[2] - v_low) / (v_high - v_low);
                    w = (w < 0) ? 0 : w;
                    w = (w > 1) ? 1 : w;
                    rgb1.v[0] = (1 - w) * argb.v[0] + w * prgb.v[0];
                    rgb1.v[1] = (1 - w) * argb.v[1] + w * prgb.v[1];
                    rgb1.v[2] = (1 - w) * argb.v[2] + w * prgb.v[2];
                }
            }
            else if (state->tc_type == TC_SIMPLE_ACR_HUE) {

                rgb1 = argb;

            } else if (state->tc_type == TC_SIMPLE_RGB_HUE) {

                rgb1 = prgb;
                if (hsv.v[0] > 6 * 240 / 360 && hsv.v[0] < 6 * 300 / 360) {
                    // prevent blue shift to purple
                    rgb1 = argb;
                }
            }
            else
                assert(false);

            const v3 hsv1 = dcp_rgb2hsv(rgb1);
            const v3 xyz1 = m3x3_mul_v3(state->rgb2xyz, rgb1);

            // limit saturation

            // desaturate, keep hue & Y
            // bisection between rgb1 and gray of rgb1
            v3 rgb2 = rgb1;
            v3 hsv2 = hsv1;
            v3 xyz2 = xyz1;
            if ((    state->tc_type == TC_SIMPLE
                  || state->tc_type == TC_SIMPLE_ACR_HUE
                  || state->tc_type == TC_SIMPLE_RGB_HUE)
                && hsv2.v[1] > hsv.v[1]) {
                // bisection
                double y1 = xyz1.v[1];
                double k1 = 0, k2 = 1;
                while (k2 - k1 > 1e-9) {
                    double k = 0.5 * (k1 + k2);                    
                    rgb2.v[0] = k * rgb1.v[0] + (1 - k) * y1;
                    rgb2.v[1] = k * rgb1.v[1] + (1 - k) * y1;
                    rgb2.v[2] = k * rgb1.v[2] + (1 - k) * y1;
                    hsv2 = dcp_rgb2hsv(rgb2);
                    if (hsv2.v[1] > hsv.v[1])
                        k2 = k;
                    else
                        k1 = k;
                }
                xyz2 = m3x3_mul_v3(state->rgb2xyz, rgb2);
            }
            assert(fabs(xyz1.v[1] / xyz2.v[1] - 1) < 1e-8);

            v3 ref_xyz = xyz2;
            ref_rgb = m3x3_mul_v3(state->xyz2rgb, ref_xyz);

            // limit result luminance Y
            // darken red-yellow range
            if (state->tc_type == TC_SIMPLE
                || state->tc_type == TC_SIMPLE_ACR_HUE
                || state->tc_type == TC_SIMPLE_RGB_HUE) {
                double y2 = xyz2.v[1];
                double y_tc = curve_y(xyz.v[1], state->tc, state->tc_len);
                // double y_tc = axyz.v[1];
                // variants: to acr or pure rgb lum
                double ky = y_tc / y2;
                if ((y2 > 1e-8) && (ky  < (1 - 1e-6))) {
                    ref_xyz = k_mul_v3(ky, xyz2);
                    ref_rgb = m3x3_mul_v3(state->xyz2rgb, ref_xyz);
                }
            }
            if (v3_isclip01(ref_rgb))
                exit(1);
            v3 ref_hsv = dcp_rgb2hsv(ref_rgb);
            assert(!v3_isclip01(ref_rgb));
            ref_rgb = v3_clip01(ref_rgb);

        } else {
            assert(false);
        }

        if (!v3_isfinite(ref_rgb)) {
            elog("Bug: non-finite RGB value in tone reproduction.\n");
            abort();
        }

        *out = m3x3_mul_v3(state->rgb2xyz, ref_rgb);
        return;
    }

    double Yr = state->rgb2xyz.v[1][0];
    double Yg = state->rgb2xyz.v[1][1];
    double Yb = state->rgb2xyz.v[1][2];
    assert(Yr > 0);
    assert(Yg > 0);
    assert(Yb > 0);
    { // make channel luminance more even in extreme blue range, which protects deep blues from becoming too dark

        double blue_factor = (srgb_gamma_forward(orig_rgb.v[2]) - srgb_gamma_forward(orig_rgb.v[0])) / srgb_gamma_forward(v3_max(orig_rgb));
        if (blue_factor < 0) blue_factor = 0;
        // this shift towards 0.60 on blue may look extreme, but as blue is natural weak in luminance the effect is smaller than it seems
        const double eYr = blue_factor * 0.20 + (1-blue_factor) * Yr;
        const double eYb = blue_factor * 0.60 + (1-blue_factor) * Yb;

        // if low green (luminance/center band) and high saturation then even out more
        double x = 1 - srgb_gamma_forward(orig_rgb.v[1]);
        x *= rgb2hsv(orig_rgb).v[1];
        Yb = (1-x)*Yb + x*eYb;
        Yr = (1-x)*Yr + x*eYr;
        Yg = 1 - Yr - Yb;
    }

    double oldLuminance = orig_rgb.v[0]*Yr + orig_rgb.v[1]*Yg + orig_rgb.v[2]*Yb;
    double acurve_newLuminance = argb.v[0]*Yr + argb.v[1]*Yg + argb.v[2]*Yb;
    double newLuminance = acurve_newLuminance;

    if (conf->curve.keep_factor != 1.0) {
        // here we use LCh chroma instead of JCh due to better stability in extreme blue range
        v3 lch = xyz2lch(orig_xyz, state->wp);
        const double origC = lch.v[1];
        const double lolim = conf->curve.low_chroma;
        const double hilim = conf->curve.high_chroma;
        double useL = 1.0;
        if (origC < lolim) {
            useL = 0.0;
        } else if (origC < hilim) {
            double x = (origC - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
            x = scurve(x);
            useL *= x;
        }
        useL *= (1.0 - conf->curve.keep_factor);
        double lumcurve_newLuminance = curve_y(oldLuminance, state->tc, state->tc_len);
        newLuminance = useL * lumcurve_newLuminance + (1.0 - useL) * acurve_newLuminance;
    }
    double Lcoef = newLuminance/oldLuminance;

    // This scaling can push rgb values >1, say to 2.1 max on blue.
    v3 rgb = {{ orig_rgb.v[0] * Lcoef, orig_rgb.v[1] * Lcoef, orig_rgb.v[2] * Lcoef }};

    { // mix in pure RGB curve for the problematic red-orange-yellow highlight range.
        // This provides better gradients in sunsets (avoid "fried egg" effect) thanks to less
        // aggressive clipping, at the cost of some hue accuracy, meaning that reds become orange
        // near clipping
        //
        // This is very good for sunsets, but hurts performance a little bit on say high saturation red flowers.
        // Everyone does it though, so it's an expected effect.

        double x = (v3_max(rgb) - 0.65) * 1 / (1-0.65);
        //double x = (v3_max(rgb) - 0.50) * 1 / (1-0.50);
        if (x < 0) {
            x = 0;
        } else if (x > 1) {
            x = 1;
        }
        if (x > 0) {
            double hue = rgb2hsv(prgb).v[0] * 360;
            v3 alt_rgb = prgb;
            if (hue < 300) {
                x = fade(hue, 60, 75, x, 0, true);
            } else {
                // in this fadeout zone we use argb (stable RGB hue)
                alt_rgb = argb;
                x = fade(hue, 300, 360, 0, x, true);
            }

            { // small lightness compensation for better matching as luminance curve is a bit brighter
                double lum1 = rgb.v[0]*Yr + rgb.v[1]*Yg + rgb.v[2]*Yb;
                double lum2 = prgb.v[0]*Yr + prgb.v[1]*Yg + prgb.v[2]*Yb;
                double mul = lum2 == 0 ? 1 : lum1/lum2;
                if (mul > 1) {
                    mul = 1 + (mul-1) * 0.5;
                    alt_rgb = k_mul_v3(mul, alt_rgb);
                    alt_rgb = gamut_saturated_rgb_clip(alt_rgb);
                }
            }
            if (x > 0) {
                for (int i = 0; i < 3; i++) {
                    rgb.v[i] = srgb_gamma_forward(alt_rgb.v[i]) * x + srgb_gamma_forward(rgb.v[i]) * (1 - x);
                    rgb.v[i] = srgb_gamma_inverse(rgb.v[i]);
                }
            }
        }
    }

    xyz = m3x3_mul_v3(state->rgb2xyz, rgb);
    xyz = gamut_clip_to_valid(xyz, state->wp); // clip to valid CIECAM02 gamut, rgb values can still be >1
    rgb = m3x3_mul_v3(state->xyz2rgb, xyz);

    v3 jch = xyz2jch_s(xyz, h02);
    if (!v3_isfinite(jch)) {
        elog("Bug: CIECAM02 failed despite pre-clipping to valid gamut.\n");
        abort();
    }

    double cmul = conf->chroma_scaling;

    // depending on color, the chroma scaling factor can be fine-tuned below

    if (cmul > 1 && conf->saturated.adjust_factor != 1.0 && conf->saturated.adjust_factor >= 0.0) {
        // decrease chroma scaling sligthly of extremely saturated colors

        double saturated_scale_factor = conf->saturated.adjust_factor;
        double chroma = jch.v[1];

        { // blue correction: decrease faster and more
            double hue = rgb2hsv(orig_rgb).v[0] * 360;
            if (hue > 200 && hue < 280) {
                double x;
                if (hue < 220) {
                    x = (hue - 200) / 20;
                } else if (hue > 260) {
                    x = (280 - hue) / 20;
                } else {
                    x = 1;
                }
                chroma *= 1.2 * x;
                saturated_scale_factor *= 1.0 - 0.05 * x;
            }
        }

        // limit decrease to canceling out the chroma scaling
        if (1 / conf->chroma_scaling > saturated_scale_factor) {
            saturated_scale_factor = 1 / conf->chroma_scaling;
        }
        const double lolim = conf->saturated.low_chroma;
        const double hilim = conf->saturated.high_chroma;

        if (chroma < lolim) {
            // chroma is low enough, don't scale
            saturated_scale_factor = 1.0;
        } else if (chroma < hilim) {
            // S-curve transition between low and high limit
            double x = (chroma - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
            x = scurve(x);
            saturated_scale_factor = 1.0*(1.0-x) + saturated_scale_factor*x;
        } else {
            // do nothing, high saturation color, keep scale factor
        }
        cmul *= saturated_scale_factor;
    }

    if (!is_linear && conf->shadows.adjust_factor != 1.0 && conf->shadows.adjust_factor >= 0.0) {
        // increase chroma scaling slightly of shadows
        // the purpose of this is to avoid the "gray dullness" of the luminance curve in shadowy parts of a face for example

        double nL = srgb_gamma_forward(newLuminance); // apply gamma so we make comparison and transition with a more perceptual lightness scale
        double dark_scale_factor = conf->shadows.adjust_factor;

        { // transition into high chroma adjust factor
            double sat_sf = conf->shadows.adjust_factor_high_chroma;
            const double lolim = conf->shadows.low_chroma;
            const double hilim = conf->shadows.high_chroma;
            if (jch.v[1] < lolim) {
                // chroma is low enough, don't change
            } else if (jch.v[1] < hilim) {
                double x = (jch.v[1] - lolim) / (hilim - lolim);
                x = scurve(x);
                dark_scale_factor = dark_scale_factor*(1.0-x) + sat_sf*x;
            } else {
                dark_scale_factor = sat_sf;
            }
        }
        const double lolim = conf->shadows.low_lightness;
        const double hilim = conf->shadows.high_lightness;
        if (nL < lolim) {
            // do nothing, keep scale factor
        } else if (nL < hilim) {
            // S-curve transition
            double x = (nL - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
            x = scurve(x);
            dark_scale_factor = dark_scale_factor*(1.0-x) + 1.0*x;
        } else {
            dark_scale_factor = 1.0;
        }
        cmul *= dark_scale_factor;
    }

    if (cmul > 1) { // to avoid strange CIECAM02 chroma errors on close-to-shadow-clipping colors we reduce chroma scaling towards 1.0 for very-close-to-black
        const double lolim = 4;
        const double hilim = 7;
        if (jch.v[0] < lolim) {
            cmul = 1.0; // very dark color: don't scale chroma
        } else if (jch.v[0] < hilim) {
            double x = (jch.v[0] - lolim) / (hilim - lolim);
            x = scurve(x);
            cmul = (1.0 - x) + cmul * x;
        }
    }

    if (cmul > 1) { // decrease scaling towards gamut limit
        double max_chroma = gamut_max_valid_chroma(jch, state->wp, 1000, h02);
        const double lolim = max_chroma * 0.5;
        const double hilim = max_chroma;
        //char name[64]; elog("%f %f %f %s\n", lolim, hilim, jch.v[2], xyz2name(xyz, state->wp, name));
        if (jch.v[1] < lolim) {
            // do nothing, keep scale factor
        } else if (jch.v[1] < hilim) {
            double x = (jch.v[1] - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
            x = scurve(x);
            cmul = 1 + (cmul - 1) * (1 - x);
        } else {
            cmul = 1;
        }
    }

    if (cmul != 1) { // apply new chroma, with safety limits
        double new_chroma = jch.v[1] * cmul;
        const double old_chroma = jch.v[1];

        jch.v[1] = new_chroma;
        xyz = jch2xyz_s(jch, h02);
        while (!v3_isfinite(xyz)) {
            // Not sure this can happen (any longer) when we decrease scaling towards gamut limit
            new_chroma *= 0.99;
            if (new_chroma < old_chroma) new_chroma = old_chroma;
            jch.v[1] = new_chroma;
            xyz = jch2xyz_s(jch, h02);

            if (v3_isfinite(xyz)) {
                new_chroma *= 0.98; // add some distance to the limit
                if (new_chroma < old_chroma) new_chroma = old_chroma;
                jch.v[1] = new_chroma;
                xyz = jch2xyz_s(jch, h02);
                break;
            }
        }
        if (!gamut_isvalid(xyz)) {
            xyz = gamut_clip_to_valid(xyz, state->wp);
            jch = xyz2jch_s(xyz, h02);
        }
    }

    { // rolloff desaturation matched with an RGB curve, as it gives best perceptual result (based on tests)
        rgb = m3x3_mul_v3(state->xyz2rgb, xyz);
        const v3 urgb = rgb;
        rgb = gamut_saturated_rgb_clip(rgb);
        const v3 hsv = rgb2hsv(rgb);
        const v3 ahsv = rgb2hsv(argb);
        double keep = lookop_curve(conf->rolloff.keep_factor, hsv.v[0] * 360);
        const double lolim = lookop_curve(conf->rolloff.low_satscale, hsv.v[0] * 360);
        const double hilim = lookop_curve(conf->rolloff.high_satscale, hsv.v[0] * 360);
        double sat_scale = ahsv.v[1] <= 0.0 ? 1.0 : hsv.v[1] / ahsv.v[1]; // saturation scale compared to Adobe curve
        if (sat_scale < lolim) {
            keep = 1.0;
        } else if (sat_scale < hilim) {
            double x = (sat_scale - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
            x = scurve(x);
            keep = 1.0*(1.0-x) + keep*x;
        } else {
            // do nothing, very high increase, keep minimum amount
        }
        if (newLuminance < oldLuminance) { // only do this for colors that are brightened
            double ls = srgb_gamma_forward(newLuminance) / srgb_gamma_forward(oldLuminance);
            const double lolim = 0.75;
            const double hilim = 1.0;
            if (ls < lolim) {
                keep = 1.0;
            } if (ls < hilim) {
                double x = (ls - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
                x = scurve(x);
                keep = 1.0*(1.0-x) + keep*x;
            }
        }
        if (keep < 1.0) {
            // mix in some of the Adobe chroma (if it's lower)

            // Note: HSV-hue matching is probably overkill, mostly (only?) clipping situations
            // and then furhter desaturation takes place anyway. But we have it here for safety.
            v3 argb1 = hsv2rgb((v3){{ rgb2hsv(v3_norm(urgb)).v[0], ahsv.v[1], ahsv.v[2] }});
            v3 acr_jch = rgb2jch_s(argb1, state->rgb2xyz, h02);
            if (!v3_isfinite(acr_jch)) {
                // rare, but can happen for extreme blue
                v3 axyz = m3x3_mul_v3(state->rgb2xyz, argb1);
                axyz = gamut_clip_to_valid(xyz, state->wp);
                acr_jch = xyz2jch_s(axyz, h02);
            }
            double acr_chroma = acr_jch.v[1];
            if (acr_chroma < jch.v[1]) {
                // new color is made from the unclipped jch value
                jch = (v3){{ jch.v[0], jch.v[1] * keep + acr_chroma * (1.0 - keep), jch.v[2] }};
                rgb = jch2rgb_s(jch, state->xyz2rgb, h02);
                rgb = gamut_saturated_rgb_clip(rgb); // this may reduce chroma further
            }
        }
    }

    if (state->conf.lookop_count > 0) { // apply look operators, if any
        jch = rgb2jch_s(rgb, state->rgb2xyz, h02);
        if (!v3_isfinite(jch)) {
            elog("Bug: CIECAM02 failed despite pre-clipping to valid gamut.\n");
            abort();
        }
        apply_look_operators(jch, &jch, state, h02);
        rgb = jch2rgb_s(jch, state->xyz2rgb, h02);
        rgb = gamut_saturated_rgb_clip(rgb);
    }
    if (!v3_isfinite(rgb)) {
        elog("Bug: non-finite RGB value in tone reproduction.\n");
        abort();
    }
    *out = m3x3_mul_v3(state->rgb2xyz, rgb);
}
