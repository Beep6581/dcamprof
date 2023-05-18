 /*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <assert.h>

#include <look.h>
#include <interp.h>
#include <dngref.h>
#include <dnglut.h>
#include <gamut.h>
#include <bisection.h>
#include <elog.h>

static void
fill_undef_hsm_entries(v3 *hsm,
                       int hcount,
                       int scount,
                       int vcount)
{
    v3 *hsm1 = calloc(1, hcount*vcount*scount * sizeof(hsm1[0]));
    bool has_undef;
    do {
        has_undef = false;
#pragma omp parallel for
        for (int v = 0; v < vcount; v++) {
            for (int h = 0; h < hcount; h++) {
                for (int s = 0; s < scount; s++) {
                    v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
                    if (v3_isfinite(*entry)) {
                        hsm1[(hcount*scount)*v + scount*h + s] = *entry;
                        continue;
                    }

                    // make average of surrounding
                    int count = 0;
                    v3 val = {{ 0, 0, 0 }};
                    for (int h1 = h-2; h1 <= h+2; h1++) {
                        if (h1 < 0 || h1 >= hcount) continue;
                        for (int s1 = s-2; s1 <= s+2; s1++) {
                            if (s1 < 0 || s1 >= scount) continue;
                            for (int v1 = v-2; v1 <= v+2; v1++) {
                                if (v1 < 0 || v1 >= vcount) continue;
                                v3 *e1 = &hsm[(hcount*scount)*v1 + scount*h1 + s1];
                                //elog("  %d %d %d | %f %f %f\n", h1, s1, v1, e1->v[0], e1->v[1], e1->v[2]);
                                if (!v3_isfinite(*e1)) continue;
                                // we don't handle hue angle wrapping here, but as the DNG/Adobe reference code doesn't either we don't need to
                                count++;
                                val.v[0] += e1->v[0];
                                val.v[1] += e1->v[1];
                                val.v[2] += e1->v[2];
                            }
                        }
                    }
                    if (count == 0) {
                        has_undef = true;
                        hsm1[(hcount*scount)*v + scount*h + s] = *entry;
                    } else {
                        hsm1[(hcount*scount)*v + scount*h + s] = (v3){{ val.v[0] / count, val.v[1] / count, val.v[2] / count }};
                    }
                    //elog("%d %d %d (%d) | %f %f %f\n", h, s, v, count, entry->v[0], entry->v[1], entry->v[2]);
                }
            }
        }
        memcpy(hsm, hsm1, hcount*vcount*scount * sizeof(hsm1[0]));
    } while (has_undef);
    free (hsm1);
}

static void
clip_and_smooth_hsm(v3 *hsm,
                    int hcount,
                    int scount,
                    int vcount)
{
#pragma omp parallel for
    for (int v = 0; v < vcount; v++) {
        for (int h = 0; h < hcount; h++) {
            for (int s = 0; s < scount; s++) {
                v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
                if ((v == 0 && vcount > 1) || s == 0) {
                    entry->v[0] = 0.0;
                    entry->v[1] = 1.0;
                    entry->v[2] = 1.0;
                    continue;
                }
                v3 hsv = {{ h / (double)hcount * 6.0, s / (double)(scount-1), v / (double)(vcount-1) }};
                if (hsv.v[1] * entry->v[1] > 1.0) {
                    entry->v[1] = 1.0 / hsv.v[1];
                }
                if (hsv.v[2] * entry->v[2] > 1.0) {
                    entry->v[2] = 1.0 / hsv.v[2];
                }
                // round meaninglessy small deviations
                if (fabs(entry->v[0]) < 0.000200) {
                    entry->v[0] = 0;
                }
                if (entry->v[1] > 0.999990 && entry->v[1] < 1.000010) {
                    entry->v[1] = 1.0;
                }
                if (entry->v[2] > 0.999990 && entry->v[2] < 1.000010) {
                    entry->v[2] = 1.0;
                }
            }
        }
    }
}

void
dnglut_test_discontinuity(v3 *hsm,
                          int hcount,
                          int scount,
                          int vcount,
                          const char tablename[],
                          bool allow_discontinuity_hue_shifts)
{
    int discontinuity_count = 0;
#pragma omp parallel for
    for (int v = 0; v < vcount; v++) {
        for (int h = 0; h < hcount; h++) {
            for (int s = 0; s < scount; s++) {
                v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
                double hs = entry->v[0];
                bool has_discontinuity = false;
                for (int h1 = h-1; h1 <= h+1 && !has_discontinuity; h1++) {
                    for (int s1 = s-1; s1 <= s+1 && !has_discontinuity; s1++) {
                        for (int v1 = v-1; v1 <= v+1 && !has_discontinuity; v1++) {
                            if (h1 < 0 || h1 >= hcount) continue;
                            if (s1 < 0 || s1 >= scount) continue;
                            if (v1 < 0 || v1 >= vcount) continue;
                            v3 *entry1 = &hsm[(hcount*scount)*v1 + scount*h1 + s1];
                            double hs1 = entry1->v[0];
                            double avg1 = (hs1 + hs) * 0.5;
                            double avg2 = 0;
                            {
                                double h1 = fmod(hs1, 360);
                                double h2 = fmod(hs, 360);
                                double d = h2 - h1;
                                if (h1 > h2) {
                                    double tmp = h1;
                                    h1 = h2;
                                    h2 = tmp;
                                    d = -d;
                                }
                                if (d > 180) {
                                    h1 = h1 + 360;
                                    avg2 = fmod((h1 + 0.5 * (h2 - h1)), 360);
                                }
                                if (d <= 180) {
                                    avg2 = h1 + 0.5 * d;
                                }
                            }
                            if (fabs(avg1 - avg2) > 0.1) {
#pragma omp atomic
                                discontinuity_count++;
                                //elog("%d %d %d %f %f (%f %f)\n", h, s, v, hs, hs1, avg1, avg2);
                                has_discontinuity = true;
                            }
                        }
                    }
                }
            }
        }
    }
    if (discontinuity_count > 0) {
        elog("%d (%.3f%%) of the %s entries have hue shift discontinuity with neighbors.\n",
                discontinuity_count, 100.0 * (double)discontinuity_count / (hcount * scount * vcount), tablename);
        if (!allow_discontinuity_hue_shifts) {
            elog("Hue shift discontinuity between LUT entry neighbors are not allowed. Aborting.\n");
            exit(EXIT_FAILURE);
        } else {
            elog("Warning: hue shift discontinuity between LUT entry neighbors are not supported by most DNG pipelines.\n");
        }
    }

}

v3 *
dnglut_looktable_new(uint32_t hsmdims[3],
                     int hcount,
                     int scount,
                     int vcount,
                     const double tc[],
                     int tc_len,
                     bool skip_tc_apply,
                     bool srgb_gamma,
                     enum gc_type gc_type,
                     const struct tone_rep_op_config *ntro_conf)
{
    elog("Generating 3D LookTable with %dx%dx%d = %d entries for the tone reproduction operator...\n  0%%", hcount, scount, vcount, hcount*scount*vcount);
    const m3x3 xyz2rgb = dcp_xyz_D50_to_prophoto_rgb;
    const m3x3 rgb2xyz = m3x3_invert(xyz2rgb);
    v3 *hsm = malloc(hcount * scount * vcount * sizeof(*hsm));
    look_tone_rep_op_t *ntro = look_tone_rep_op_new(rgb2xyz, tc, tc_len, gc_type, ntro_conf);

    const int itc_len = 65535;
    double *itc = malloc(itc_len * sizeof(itc[0]));
    for (int i = 0; i < itc_len; i++) {
        double x = (double)i / (itc_len - 1);
        itc[i] = inverse_curve(x, tc, tc_len);
    }
    int complete_count = 0;
    int last_percent = 0;

#pragma omp parallel for
    for (int v = 0; v < vcount; v++) {
        for (int h = 0; h < hcount; h++) {
            for (int s = 0; s < scount; s++) {
                int percent = (100 * complete_count) / (vcount * hcount * scount);
#pragma omp critical
                {
                    if (percent < 100 && percent >= last_percent + 10) {
                        elog("..%d%%", percent);
                        last_percent = percent;
                    }
                    complete_count++;
                }
                v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
                if (v == 0 || s == 0) {
                    /* DNG profile format requires these values for s or v = 0. That saturation can't
                       be scaled from zero is logical, but that we cannot change value of a nonzero
                       gray is unfortunate. */
                    entry->v[0] = 0.0;
                    entry->v[1] = 1.0;
                    entry->v[2] = 1.0;
                    continue;
                }
                v3 hsv = {{ h / (double)hcount * 6.0, s / (double)(scount-1), v / (double)(vcount-1) }};
                if (srgb_gamma) {
                    hsv.v[2] = srgb_gamma_inverse(hsv.v[2]);
                }
                v3 ref_rgb;
                { // calculate corrected HSV value, and HSV after tone curve application
                    v3 rgb = dcp_hsv2rgb(hsv);
                    v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
                    v3 ref_xyz;
                    look_tone_rep_op(&ref_xyz, xyz, ntro);
                    ref_rgb = m3x3_mul_v3(xyz2rgb, ref_xyz);
                    bool clip = false;
                    const double tolerance = 1e-4;
                    for (int i = 0; i < 3; i++) {
                        if (ref_rgb.v[i] < 0.0) {
                            if (ref_rgb.v[i] < -tolerance) {
                                clip = true;
                            } else {
                                ref_rgb.v[i] = 0.0;
                            }
                        }
                        if (ref_rgb.v[i] > 1.0) {
                            if (ref_rgb.v[i] > 1.0 + tolerance) {
                                clip = true;
                            } else {
                                ref_rgb.v[i] = 1.0;
                            }
                        }
                    }
                    if (clip || !v3_isfinite(ref_rgb)) {
#pragma omp critical
                        {
                            v3_print("xyz", ref_xyz);
                            v3_print("rgb", ref_rgb);
                            elog("Bug: LookTable: RGB out of range.\n");
                            abort();
                        }
                    }
                    if (skip_tc_apply) {
                        ref_rgb = dcp_apply_tonecurve(ref_rgb, itc, itc_len);
                    }
                }
                v3 ref_hsv = dcp_rgb2hsv(ref_rgb);

                entry->v[0] = (ref_hsv.v[0] - hsv.v[0]) / 6.0 * 360.0;
                if (entry->v[0] < -180) {
                    entry->v[0] += 360.0;
                } else if (entry->v[0] > 180) {
                    entry->v[0] -= 360.0;
                }
                if (fabs(entry->v[0]) < 1e-7) {
                    entry->v[0] = 0.0;
                }
                entry->v[1] = ref_hsv.v[1] / hsv.v[1];
                entry->v[2] = ref_hsv.v[2] / hsv.v[2];

                if (srgb_gamma) {
                    entry->v[2] = srgb_gamma_forward(hsv.v[2] * entry->v[2]) / srgb_gamma_forward(hsv.v[2]);
                }
                assert(isfinite(entry->v[0]));
                assert(isfinite(entry->v[1]));
                assert(isfinite(entry->v[2]));
            }
        }
    }
    hsmdims[0] = hcount;
    hsmdims[1] = scount;
    hsmdims[2] = vcount;
    look_tone_rep_op_delete(ntro);
    free(itc);

    clip_and_smooth_hsm(hsm, hcount, scount, vcount);
    elog("..100%%\n");

    return hsm;
}

struct find_max_cam_arg {
    v3 hsv;
    m3x3 rgb2xyz;
    m3x3 xyz2cam;
    v3 last_valid_cam;
};

static double
find_max_cam(double x, void *arg)
{
    struct find_max_cam_arg *a = (struct find_max_cam_arg *)arg;
    v3 hsv = a->hsv;
    hsv.v[2] = x;
    v3 rgb = dcp_hsv2rgb(hsv);
    v3 xyz = m3x3_mul_v3(a->rgb2xyz, rgb);
    v3 cam = m3x3_mul_v3(a->xyz2cam, xyz);
    if (v3_max(cam) > 1.0) {
        return 1 + x;
    }
    a->last_valid_cam = cam;
    return 1 - x;
}

v3 *
dnglut_huesatmap_new(uint32_t hsmdims[3],
                     int hcount,
                     int scount,
                     int vcount,
                     bool srgb_gamma,
                     chromalut_t *lut,
                     const v3 lut_d50,
                     const m3x3 cam2xyz_fm,
                     const m3x3 cam2xyz_lm)
{
    const double inf = 1.0 / 0.0;
    m3x3 rgb2xyz = m3x3_invert(dcp_xyz_D50_to_prophoto_rgb);
    v3 *hsm = malloc(hcount * scount * vcount * sizeof(*hsm));

    if (vcount == 1) {
        // srgb gamma only makes sense if we have a value axis
        srgb_gamma = false;
    }
#pragma omp parallel for
    for (int vh = 0; vh < vcount * hcount; vh++) {
        int v = vh / hcount;
        int h = vh % hcount;
        for (int s = 0; s < scount; s++) {
            v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
            if ((v == 0 && vcount > 1) || s == 0) {
                entry->v[0] = 0.0;
                entry->v[1] = 1.0;
                entry->v[2] = 1.0;
                continue;
            }
            // for 2.5D tables (no value axis) we sample at value 0.5
            v3 hsv = {{ h / (double)hcount * 6.0, s / (double)(scount-1), vcount == 1 ? 0.5 : v / (double)(vcount - 1) }};
            if (srgb_gamma) {
                hsv.v[2] = srgb_gamma_inverse(hsv.v[2]);
            }
            v3 ref_hsv;
            { // calculate corrected HSV value
                v3 rgb = dcp_hsv2rgb(hsv);
                v3 xyz_lm = m3x3_mul_v3(rgb2xyz, rgb);
                v3 cam = m3x3_mul_v3(m3x3_invert(cam2xyz_lm), xyz_lm);
                double mx = v3_max(cam);
                if (mx > 1.0) {
                    // this won't normally happen for vcount == 1 (value set to 0.5)
                    struct find_max_cam_arg arg = {
                        .hsv = hsv,
                        .rgb2xyz = rgb2xyz,
                        .xyz2cam = m3x3_invert(cam2xyz_lm),
                    };
                    interval_halving_min(find_max_cam, &arg, 0, 1, 1e-9, 1000);
                    cam = arg.last_valid_cam;
                }
                if (v3_isclip01(cam)) {
                    // negative clipping, chromaticity outside the camera matrix coverage
                    // this entry in the table cannot really be reached as this means the matrix doesn't cover ProPhotoRGB
                    //v3_print("cam1", cam);

                    // set to undefined to fill in later
                    entry->v[0] = inf;
                    entry->v[1] = inf;
                    entry->v[2] = inf;
                    continue;
                }
                v3 xyz = m3x3_mul_v3(cam2xyz_fm, cam);
                //printf("%f %f %f\n", xyz.v[0], xyz.v[1], xyz.v[2]);
                xyz = m3x3_mul_v3(bradford_map_white(dcp_d50(), lut_d50), xyz);
                xyz = chromalut_lookup(lut, xyz);
                xyz = m3x3_mul_v3(bradford_map_white(lut_d50, dcp_d50()), xyz);
                rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, xyz);
                rgb = gamut_saturated_rgb_clip(rgb);
                ref_hsv = dcp_rgb2hsv(rgb);
            }
            //elog("%f %f %f\n", hsv.v[0], hsv.v[1], hsv.v[2]);
            //elog("%f %f %f\n", ref_hsv.v[0], ref_hsv.v[1], ref_hsv.v[2]);

            entry->v[0] = (ref_hsv.v[0] - hsv.v[0]) / 6.0 * 360.0;
            if (entry->v[0] < -180) {
                entry->v[0] += 360.0;
            } else if (entry->v[0] > 180) {
                entry->v[0] -= 360.0;
            }
            entry->v[1] = ref_hsv.v[1] / hsv.v[1];
            entry->v[2] = ref_hsv.v[2] / hsv.v[2];
            //elog("%f %f %f\n", entry->v[0], entry->v[1], entry->v[2]);
            if (srgb_gamma) {
                entry->v[2] = srgb_gamma_forward(hsv.v[2] * entry->v[2]) / srgb_gamma_forward(hsv.v[2]);
            }
        }
        //elog("\n");
    }
    fill_undef_hsm_entries(hsm, hcount, scount, vcount);
    if (vcount > 1) {
        clip_and_smooth_hsm(hsm, hcount, scount, vcount);
    } else {
        // clipping doesn't make sense for vcount == 1, as we then would clip for the worst case and reduce gamut for normal colors
    }

    hsmdims[0] = hcount;
    hsmdims[1] = scount;
    hsmdims[2] = vcount;

    return hsm;
}
