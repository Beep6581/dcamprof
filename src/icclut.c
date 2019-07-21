 /*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include <spectraldb.h>
#include <observers.h>
#include <interp.h>
#include <look.h>
#include <elog.h>
#include <lut.h>
#include <gamut.h>
#include <dngref.h>
#include <icclut.h>

/*
static uint32_t
fletcher32(uint16_t const *data,
           size_t words)
{
    uint32_t sum1 = 0xffff, sum2 = 0xffff;

    while (words) {
        unsigned tlen = words > 359 ? 359 : (unsigned)words;
        words -= tlen;
        do {
            sum2 += sum1 += *data++;
        } while (--tlen);
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    return sum2 << 16 | sum1;
}
*/

static v3
icclut_iccv2clip(v3 lab)
{
    if (!v3_isfinite(lab)) {
        v3_print("bad lab", lab);
        abort();
    }
    // limit to the maximum representable values in ICCv2
    if (lab.v[0] < 0) lab.v[0] = 0;
    if (lab.v[0] > 100 + 25500/65280.0) lab.v[0] = 100 + 25500/65280.0;
    v3 lch = lab2lch(lab);
    while (lab.v[1] < -128.0  || lab.v[1] > 127 + 255/256.0 || lab.v[2] < -128.0  || lab.v[2] > 127 + 255/256.0) {
        double h = lch.v[2] / 180.0 * M_PI;
        if (lab.v[1] < -128.0) {
            if (lab.v[1] > -129.0) {
                lab.v[1] = -128.0;
            } else {
                lch.v[1] = -128.0 / cos(h);
                lab = lch2lab(lch);
            }
        } else if (lab.v[1] > 127 + 255/256.0) {
            if (lab.v[1] < 129.0) {
                lab.v[1] = 127 + 255/256.0;
            } else {
                lch.v[1] = (127 + 255/256.0) / cos(h);
                lab = lch2lab(lch);
            }
        } else if (lab.v[2] < -128.0) {
            if (lab.v[2] > -129.0) {
                lab.v[2] = -128.0;
            } else {
                lch.v[1] = -128.0 / sin(h);
                lab = lch2lab(lch);
            }
        } else if (lab.v[2] > 127 + 255/256.0) {
            if (lab.v[2] < 129.0) {
                lab.v[2] = 127 + 255/256.0;
            } else {
                lch.v[1] = (127 + 255/256.0) / sin(h);
                lab = lch2lab(lch);
            }
        }
    }
    /*
    if (!v3_isequal(lab, orig_lab)) {
        v3 xyz = lab2xyz(orig_lab, icc_d50());
        v3 ltu = xyz2lutspace(xyz);
        printf("%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
        v3 rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, xyz);
        //v3_print("rgb", rgb);
    }
    */
    return lab;
}

static double
cam_similarity(v3 cam1,
               v3 cam2)
{
    double max1 = v3_max(cam1);
    double max2 = v3_max(cam2);
    assert(max1 != 0.0);
    assert(max2 != 0.0);
    cam1 = k_mul_v3(1.0 / max1, cam1);
    cam2 = k_mul_v3(1.0 / max2, cam2);

    double d = euclidean_dist(cam1, cam2);

    return d;
}

static inline double
icclut_incurve(double x,
               bool invert)
{
    if (!invert) {
        const double p = pow(3.0, 9/2.0);
        const double t = 1.0/3.0;
        return pow(sqrt(6075*x*x+1)/p+5*x/9, t) - 1/(27*pow(sqrt(6075*x*x+1)/p+5*x/9, t));
    } else {
        return 0.1*x+0.9*pow(x,3);
    }
}

static int
ig_offsets_cmp(const void *a_, const void *b_)
{
    const uint8_t *a = (const uint8_t *)a_;
    const uint8_t *b = (const uint8_t *)b_;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    if (a[2] < b[2]) return -1;
    if (a[2] > b[2]) return 1;
    return 0;
}

static void
smooth_lut_curve(double curve[],
                 const bool ingamut[],
                 const double chroma[],
                 const int cs)
{
    double smoothed_curve[cs];
    for (int i = 1; i < cs-1; i++) {
        smoothed_curve[i] = curve[i-1] + curve[i] + curve[i+1];
        smoothed_curve[i] /= 3;
    }
    // always keep original endpoints
    smoothed_curve[0] = curve[0];
    smoothed_curve[cs-1] = curve[cs-1];

    double smix[cs];
    {
        double smix_iny[cs];
        double smix_inx[cs];
        double smix_outx[cs];
        int smix_len = 0;
        bool gaps = false;
        for (int i = 0; i < cs; i++) {
            smix_outx[i] = (double)i;
            smix[i] = 0.0;
            if (!ingamut[i]) {
                smix_iny[smix_len] = 1.0;
                smix[i] = 1.0;
            } else if (chroma[i] < 90) {
                // only smooth at high saturation
                smix_iny[smix_len] = 0.0;
            } else {
                // don't smooth too far in to the defined gamut
                const int maxdist = 3;
                bool has_og = false;
                for (int j = i-maxdist < 0 ? 0 : i-maxdist; j <= i+maxdist && j < cs; j++) {
                    if (!ingamut[j]) {
                        has_og = true;
                        break;
                    }
                }
                if (has_og) {
                    gaps = true;
                    continue;
                }
                smix_iny[smix_len] = 0.0;
            }
            smix_inx[smix_len] = (double)smix_outx[i];
            smix_len++;
        }
        if (smix_len > 1 && gaps) {
            rounded_linear_interpolation(smix_inx, smix_iny, smix_len, smix_outx, smix, cs);
        }
    }
    for (int i = 0; i < cs; i++) {
        curve[i] = smix[i] * smoothed_curve[i] + (1.0 - smix[i]) * curve[i];
    }
}

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
           const struct tone_rep_op_config *ntro_conf)
{
    const int cs = clut_side;
    struct icc_lut *lut = calloc(1, sizeof(*lut) + 3*(256 + 2 + cs*cs*cs) * sizeof(lut->data[0]));
    double *pd = lut->data;

    if (clut_side > 255) {
        elog("Error: 3D LUT side may not be larger than 255.\n");
        exit(EXIT_FAILURE);
    }
    if (clut_side < 5) {
        elog("Error: 3D LUT side may not be smaller than 5.\n");
        exit(EXIT_FAILURE);
    }

    lut->out_len = 2;
    lut->in_len = 256;
    lut->clut_side = cs;
    for (int i = 0; i < 3; i++) {
        lut->out[i] = pd;
        pd += lut->out_len;
        lut->in[i] = pd;
        pd += lut->in_len;
    }
    lut->clut = pd;

    for (int i = 0; i < 3; i++) {
        lut->out[i][0] = 0;
        lut->out[i][1] = 1;
        for (int j = 0; j < lut->in_len; j++) {
            double x = (double)j / (lut->in_len-1);
            lut->in[i][j] = icclut_incurve(x, false);
        }
    }
    const m3x3 xyz2rgb = dcp_xyz_D50_to_prophoto_rgb;
    const m3x3 rgb2xyz = m3x3_invert(xyz2rgb);
    look_tone_rep_op_t *ntro = NULL;
    const int itc_len = 65535;
    double *itc = NULL;
    if (tc != NULL) {
        if (tc_type == TC_NEUTRAL) {
            ntro = look_tone_rep_op_new(rgb2xyz, tc, tc_len, gc_type, ntro_conf);
            itc = malloc(itc_len * sizeof(itc[0]));
            for (int i = 0; i < itc_len; i++) {
                double x = (double)i / (itc_len - 1);
                itc[i] = inverse_curve(x, tc, tc_len);
            }
        }
    }

    if (pcs_is_lab) {
        lut->cm.v[0][0] = lut->cm.v[1][1] = lut->cm.v[2][2] = 1.0;
    } else {
        lut->cm = cam2xyz;
        v3 wb_mul = v3_norm(v3_invert(wb));
        lut->cm = m3x3_mul_m3x3(lut->cm, v3_asdiagonal(wb_mul));
    }
    // normally it should be tiny differences between pcs_wp and prof_wp
    const m3x3 map_forward = bradford_map_white(pcs_wp, prof_wp);
    const m3x3 map_reverse = bradford_map_white(prof_wp, pcs_wp);
    const v3 wb_mul = v3_norm(v3_invert(wb));
    elog("Generating ICC RGB to %s LUT with %dx%dx%d = %d entries...\n  ", pcs_is_lab ? "Lab" : "XYZ", cs, cs, cs, cs*cs*cs);
    int complete_count, last_percent;

    elog("Calculating valid in-gamut entries...");
    v3 *cammap = malloc(cs * cs * cs * sizeof(cammap[0]));
    bool *ingamut_map = malloc(cs * cs * cs * sizeof(ingamut_map[0]));
    double *chroma_map = calloc(1, cs * cs * cs * sizeof(chroma_map[0]));
    uint8_t *offsets = malloc(3 * cs * cs * cs * sizeof(offsets[0]));
    int offsets_len = 0;
    uint8_t *ig_offsets = malloc(3 * cs * cs * cs * sizeof(ig_offsets[0]));
    int ig_offsets_len = 0;
    complete_count = last_percent = 0;
#pragma omp parallel for
    for (int r_ = 0; r_ < cs; r_++) {
#pragma omp critical
        {
            elog(".");
            int percent = (100 * complete_count) / cs;
            if (percent < 100 && percent >= last_percent + 10) {
                elog("%d%%", percent);
                last_percent = percent;
            }
            complete_count++;
        }
        for (int g_ = 0; g_ < cs; g_++) {
            for (int b_ = 0; b_ < cs; b_++) {
                double r = (double)r_ / (cs-1);
                double g = (double)g_ / (cs-1);
                double b = (double)b_ / (cs-1);
                v3 cam = {{ r, g, b}};
                for (int i = 0; i < 3; i++) {
                    cam.v[i] = icclut_incurve(cam.v[i], true);
                    if (pcs_is_lab) {
                        cam.v[i] *= wb_mul.v[i];
                    }
                }
                if (trc_count > 0) {
                    for (int i = 0; i < 3; i++) {
                        double v = cam.v[i];
                        if (v < 0) v = 0;
                        if (v > 1) v = 1;
                        if (trc_count == 1) {
                            cam.v[i] = pow(cam.v[i], trc[i][0]);
                        } else {
                            int idx = floor(v * (trc_count-1));
                            double d = v * (trc_count-1) - idx;
                            if (d > 1e-7) {
                                cam.v[i] = (1.0 - d) * trc[i][idx] + d * trc[i][idx+1];
                            } else {
                                cam.v[i] = trc[i][idx];
                            }
                        }
                    }
                }
                v3 cxyz = pcs_is_lab ? m3x3_mul_v3(cam2xyz, cam) : cam;
                bool in_gamut = true;
                if (prof_lut != NULL) {
                    cxyz = m3x3_mul_v3(map_forward, cxyz);
                    cxyz = chromalut_lookup(prof_lut, cxyz);
                    cxyz = m3x3_mul_v3(map_reverse, cxyz);
                    v3 rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, cxyz);
                    rgb = gamut_saturated_rgb_clip(rgb);
                    cxyz = m3x3_mul_v3(m3x3_invert(dcp_xyz_D50_to_prophoto_rgb), rgb);
                } else {
                    in_gamut = gamut_isvalid(cxyz);
                }
                /*
                if (in_gamut) {
                    if (ntro != NULL) {
                        v3 tcxyz;
                        if (look_tone_rep_op(&tcxyz, rgb2xyz, xyz2rgb, cxyz, ntro)) {
                            cxyz = tcxyz;
                        } else {
                            in_gamut = false;
                        }
                        if (skip_tc_apply) {
                            v3 rgb = v3_clip01(m3x3_mul_v3(xyz2rgb, cxyz));
                            rgb.v[0] = curve_y(rgb.v[0], itc, itc_len);
                            rgb.v[1] = curve_y(rgb.v[1], itc, itc_len);
                            rgb.v[2] = curve_y(rgb.v[2], itc, itc_len);
                            cxyz = m3x3_mul_v3(rgb2xyz, rgb);
                        }
                    } else if (tc != NULL && !skip_tc_apply) {
                        v3 rgb = v3_clip01(m3x3_mul_v3(xyz2rgb, cxyz));
                        rgb.v[0] = curve_y(rgb.v[0], tc, tc_len);
                        rgb.v[1] = curve_y(rgb.v[1], tc, tc_len);
                        rgb.v[2] = curve_y(rgb.v[2], tc, tc_len);
                        cxyz = m3x3_mul_v3(rgb2xyz, rgb);
                    }
                    if (in_gamut) {
                        in_gamut = gamut_isvalid(cxyz);
                    }
                }
                */
                if (!in_gamut) {
#pragma omp critical
                    {
                        offsets[3*offsets_len+0] = r_;
                        offsets[3*offsets_len+1] = g_;
                        offsets[3*offsets_len+2] = b_;
                        offsets_len++;
                    }
                } else {
                    if (cam.v[0] > 0 || cam.v[1] > 0 || cam.v[2] > 0) {
#pragma omp critical
                        {
                            ig_offsets[3*ig_offsets_len+0] = r_;
                            ig_offsets[3*ig_offsets_len+1] = g_;
                            ig_offsets[3*ig_offsets_len+2] = b_;
                            ig_offsets_len++;
                        }
                    }
                }
                cammap[cs*cs*r_+cs*g_+b_] = cam;
                ingamut_map[cs*cs*r_+cs*g_+b_] = in_gamut;
                chroma_map[cs*cs*r_+cs*g_+b_] = lab2lch(xyz2lab(cxyz, pcs_wp)).v[1];
                for (int i = 0; i < 3; i++) {
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+i] = cxyz.v[i];
                }
            }
        }
    }

    // the only reason for sorting this array is to get same result each time in the followup interpolation.
    qsort(ig_offsets, ig_offsets_len, 3 * sizeof(ig_offsets[0]), ig_offsets_cmp);
    elog("..100%%\n");

    if (offsets_len > 0) {
        elog("  Interpolating out of gamut entries (%.2f%% of total)...", (100.0 * offsets_len) / (cs * cs * cs));
        complete_count = last_percent = 0;
        int last_percent2 = 0;
#pragma omp parallel for
        for (int i = 0; i < offsets_len; i++) {
#pragma omp critical
            {
                int percent = (100 * complete_count) / offsets_len;
                if (percent < 100 && percent >= last_percent2 + 1) {
                    elog(".");
                    last_percent2 = percent;
                }
                if (percent < 100 && percent >= last_percent + 5) {
                    elog(".%d%%", percent);
                    last_percent = percent;
                }
                complete_count++;
            }
            const int r_ = offsets[3*i+0];
            const int g_ = offsets[3*i+1];
            const int b_ = offsets[3*i+2];
            double min_dist = 10000000;
            v3 cam = cammap[cs*cs*r_+cs*g_+b_];
            double *da = malloc(ig_offsets_len * sizeof(da[0]));
            for (int j = 0; j < ig_offsets_len; j++) {
                const int r1 = ig_offsets[3*j+0];
                const int g1 = ig_offsets[3*j+1];
                const int b1 = ig_offsets[3*j+2];
                double d = cam_similarity(cam, cammap[cs*cs*r1+cs*g1+b1]);
                da[j] = d;
                if (d < min_dist) {
                    min_dist = d;
                }
            }
            int m_count = 0;
            v3 cam_sum={{0,0,0}};
            v3 xyz_sum={{0,0,0}};
            // by multiplying the min_distance we will match few if the min_dist is small,
            // and many if it's large, which makes sense: if we found a good match we don't
            // want/need to average that much, if we didn't find a good match averaging over
            // more entries is desirable.
            min_dist *= 1.2;
            const float gamma = 2.2;
            for (int j = 0; j < ig_offsets_len; j++) {
                if (da[j] <= min_dist) {
                    const int r1 = ig_offsets[3*j+0];
                    const int g1 = ig_offsets[3*j+1];
                    const int b1 = ig_offsets[3*j+2];
                    v3 cam1 = cammap[cs*cs*r1+cs*g1+b1];
                    v3 xyz1;
                    for (int k = 0; k < 3; k++) {
                        cam1.v[k] = pow(cam1.v[k], 1/gamma);
                        xyz1.v[k] = pow(lut->clut[3*(cs*cs*r1+cs*g1+b1)+k], 1/gamma);
                    }
                    double mxc = v3_max(cam1);
                    double mxy = v3_max(xyz1);
                    double mx = mxc > mxy ? mxc : mxy;
                    cam1 = k_mul_v3(1.0 / mx, cam1);
                    xyz1 = k_mul_v3(1.0 / mx, xyz1);
                    for (int k = 0; k < 3; k++) {
                        cam_sum.v[k] += cam1.v[k];
                        xyz_sum.v[k] += xyz1.v[k];
                    }
                    m_count++;
                }
            }
            free(da);
            cam_sum = k_mul_v3(1.0 / m_count, cam_sum);
            xyz_sum = k_mul_v3(1.0 / m_count, xyz_sum);
            // Scale based on raw RGB samples (cam) rather than calculated XYZ value. Doing it on
            // XYZ may seem more "true", but since we're in out of gamut area the XYZ value doesn't
            // make much sense anyway (can often be zero or negative), and the profile becomes more
            // robust if out of gamut inputs gives some sort of sane output even if it's not
            // "correct".

            cam.v[0] = pow(cam.v[0], 1/gamma);
            cam.v[1] = pow(cam.v[1], 1/gamma);
            cam.v[2] = pow(cam.v[2], 1/gamma);
            double lum1 = v3_max(cam);
            double lum2 = v3_max(cam_sum);
            //elog("%f %f\n", lum1, lum2);
            xyz_sum = k_mul_v3(lum1 / lum2, xyz_sum);
            double mx = v3_max(xyz_sum);
            if (mx > 1.0) {
                xyz_sum = k_mul_v3(1.0 / mx, xyz_sum);
            }
            for (int i = 0; i < 3; i++) {
                lut->clut[3*(cs*cs*r_+cs*g_+b_)+i] = pow(xyz_sum.v[i], gamma);
            }
        }
        elog("..100%%\n");
    }
    free(ig_offsets);
    free(cammap);
    free(offsets);

    if (offsets_len > 0) {
        const int perm[6][3] = {
            { 1, 2, 0 },
            { 2, 1, 0 },
            { 0, 2, 1 },
            { 2, 0, 1 },
            { 1, 0, 2 },
            { 0, 1, 2 },
        };
        // here we smooth all curves, and since the all affect each-other we run it a number of times to even out
        elog("  Smoothing out of gamut entries...");

        // convert to Lab to smooth in perceptually uniform space
#pragma omp parallel for
        for (int r_ = 0; r_ < cs; r_++) {
            for (int g_ = 0; g_ < cs; g_++) {
                for (int b_ = 0; b_ < cs; b_++) {
                    v3 cxyz = {{ lut->clut[3*(cs*cs*r_+cs*g_+b_)+0], lut->clut[3*(cs*cs*r_+cs*g_+b_)+1], lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] }};
                    v3 lab = xyz2lab(cxyz, pcs_wp);
                    lab = icclut_iccv2clip(lab);
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+0] = lab.v[0];
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+1] = lab.v[1];
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] = lab.v[2];
                }
            }
        }
        {
            // map 0,0,0 to exactly 0,0,0
            lut->clut[0] = lut->clut[1] = lut->clut[2] = 0;

            // if chromaticity of 1,1,1 is close to zero, set to zero
            int offset = 3*(cs*cs*(cs-1)+cs*(cs-1)+cs-1);
            for (int i = 0; i < 2; i++) {
                double ov = lut->clut[offset+1+i];
                if (fabs(ov) < 0.01) {
                    ov = 0;
                }
                lut->clut[offset+1+i] = ov;
            }
            // if lightness of 1,1,1 is close to 100 set to 100
            if (fabs(lut->clut[offset+0] - 100.0) < 0.03) {
                lut->clut[offset+0] = 100.0;
            }
        }

        const int passes = 12;
        for (int j = 0; j < passes; j++) {
            for (int k = 0; k < 6; k++) {
#pragma omp parallel for
                for (int i1 = 0; i1 < cs; i1++) {
                    for (int i2 = 0; i2 < cs; i2++) {
                        double curve[3][cs];
                        double chroma[cs];
                        bool ingamut[cs];
                        for (int i3 = 0; i3 < cs; i3++) {
                            const int idx[3] = { i1, i2, i3 };
                            int r_ = idx[perm[k][0]], g_ = idx[perm[k][1]], b_ = idx[perm[k][2]];
                            for (int i = 0; i < 3; i++) curve[i][i3] = lut->clut[3*(cs*cs*r_+cs*g_+b_)+i];
                            ingamut[i3] = ingamut_map[cs*cs*r_+cs*g_+b_];
                            chroma[i3] = chroma_map[cs*cs*r_+cs*g_+b_];
                        }
                        for (int i = 0; i < 3; i++)  {
                            smooth_lut_curve(curve[i], ingamut, chroma, cs);
                        }
                        for (int i3 = 0; i3 < cs; i3++) {
                            const int idx[3] = { i1, i2, i3 };
                            int r_ = idx[perm[k][0]], g_ = idx[perm[k][1]], b_ = idx[perm[k][2]];
                            for (int i = 0; i < 3; i++) lut->clut[3*(cs*cs*r_+cs*g_+b_)+i] = curve[i][i3];
                        }
                    }
                }
            }
        }

        // convert back to xyz
#pragma omp parallel for
        for (int r_ = 0; r_ < cs; r_++) {
            for (int g_ = 0; g_ < cs; g_++) {
                for (int b_ = 0; b_ < cs; b_++) {
                    v3 lab = {{ lut->clut[3*(cs*cs*r_+cs*g_+b_)+0], lut->clut[3*(cs*cs*r_+cs*g_+b_)+1], lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] }};
                    v3 xyz = lab2xyz(lab, pcs_wp);
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+0] = xyz.v[0];
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+1] = xyz.v[1];
                    lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] = xyz.v[2];
                }
            }
        }
        elog("done!\n");
    }
    free(ingamut_map);
    free(chroma_map);

    if (ntro != NULL || (tc != NULL && !skip_tc_apply)) {
        // finally we apply tone reproduction / tone curve, can't do it earlier as we up till now need a linear table
        if (ntro != NULL) {
            elog("  Applying tone reproduction operator...");
        } else {
            elog("  Applying tone curve...");
        }
#pragma omp parallel for
        for (int r_ = 0; r_ < cs; r_++) {
            for (int g_ = 0; g_ < cs; g_++) {
                for (int b_ = 0; b_ < cs; b_++) {
                    v3 cxyz = {{ lut->clut[3*(cs*cs*r_+cs*g_+b_)+0], lut->clut[3*(cs*cs*r_+cs*g_+b_)+1], lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] }};
                    if (ntro != NULL) {
                        look_tone_rep_op(&cxyz, cxyz, ntro);
                        if (skip_tc_apply) {
                            v3 rgb = gamut_saturated_rgb_clip(m3x3_mul_v3(xyz2rgb, cxyz));
                            rgb.v[0] = curve_y(rgb.v[0], itc, itc_len);
                            rgb.v[1] = curve_y(rgb.v[1], itc, itc_len);
                            rgb.v[2] = curve_y(rgb.v[2], itc, itc_len);
                            cxyz = m3x3_mul_v3(rgb2xyz, rgb);
                        }
                    } else if (tc != NULL && !skip_tc_apply) {
                        v3 rgb = gamut_saturated_rgb_clip(m3x3_mul_v3(xyz2rgb, cxyz));
                        rgb.v[0] = curve_y(rgb.v[0], tc, tc_len);
                        rgb.v[1] = curve_y(rgb.v[1], tc, tc_len);
                        rgb.v[2] = curve_y(rgb.v[2], tc, tc_len);
                        cxyz = m3x3_mul_v3(rgb2xyz, rgb);
                    }
                    for (int i = 0; i < 3; i++) {
                        lut->clut[3*(cs*cs*r_+cs*g_+b_)+i] = cxyz.v[i];
                    }
                }
            }
        }
        elog("done!\n");
    }

    if (pcs_is_lab) {
        elog("  Converting from XYZ to Lab...");

        { // convert to lab
#pragma omp parallel for
            for (int r_ = 0; r_ < cs; r_++) {
                for (int g_ = 0; g_ < cs; g_++) {
                    for (int b_ = 0; b_ < cs; b_++) {
                        v3 cxyz = {{ lut->clut[3*(cs*cs*r_+cs*g_+b_)+0], lut->clut[3*(cs*cs*r_+cs*g_+b_)+1], lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] }};
                        v3 lab = xyz2lab(cxyz, pcs_wp);
                        lab = icclut_iccv2clip(lab);
                        lut->clut[3*(cs*cs*r_+cs*g_+b_)+0] = lab.v[0];
                        lut->clut[3*(cs*cs*r_+cs*g_+b_)+1] = lab.v[1];
                        lut->clut[3*(cs*cs*r_+cs*g_+b_)+2] = lab.v[2];
                    }
                }
            }
        }

        {
            // map 0,0,0 to exactly 0,0,0
            lut->clut[0] = lut->clut[1] = lut->clut[2] = 0;

            // if chromaticity of 1,1,1 is close to zero, set to zero
            int offset = 3*(cs*cs*(cs-1)+cs*(cs-1)+cs-1);
            for (int i = 0; i < 2; i++) {
                double ov = lut->clut[offset+1+i];
                if (fabs(ov) < 0.01) {
                    ov = 0;
                }
                lut->clut[offset+1+i] = ov;
            }
            // if lightness of 1,1,1 is close to 100 set to 100
            if (fabs(lut->clut[offset+0] - 100.0) < 0.03) {
                lut->clut[offset+0] = 100.0;
            }
        }

        elog("done!\n");
    }

    /*
    {
        uint32_t checksum = fletcher32((const uint16_t *)lut->data, 3*(256 + 2 + cs*cs*cs) * sizeof(lut->data[0]) / 2);
        elog("LUT data checksum: 0x%08X\n", checksum);
    }
    */

    look_tone_rep_op_delete(ntro);
    free(itc);
    elog("done!\n");
    return lut;
}
