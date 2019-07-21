/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <assert.h>

#include <interp.h>
#include <dngref.h>
#include <gamut.h>
#include <tps.h>
#include <elog.h>
#include <observers.h>
#include <spectraldb.h>

#include <glare.h>

void
glare_test(const struct patch_set *ps,
           struct glare_data *gd)
{
    double min_y = ps->patch[0].xyz.v[1];
    double max_y = ps->patch[0].xyz.v[1];
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].xyz.v[1] < min_y) min_y = ps->patch[i].xyz.v[1];
        if (ps->patch[i].xyz.v[1] > max_y) max_y = ps->patch[i].xyz.v[1];
    }
    min_y *= 1.05;
    max_y /= 1.05;
    double avg_white_g = 0;
    double avg_white_y = 0;
    double avg_black_g = 0;
    double avg_black_y = 0;
    int black_count = 0;
    int white_count = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].xyz.v[1] < min_y) {
            avg_black_g += ps->patch[i].cam.v[1];
            avg_black_y += ps->patch[i].xyz.v[1];
            black_count++;
        }
        if (ps->patch[i].xyz.v[1] > max_y) {
            avg_white_g += ps->patch[i].cam.v[1];
            avg_white_y += ps->patch[i].xyz.v[1];
            white_count++;
        }
    }
    avg_white_g /= white_count;
    avg_white_y /= white_count;
    avg_black_g /= black_count;
    avg_black_y /= black_count;
    double yd = avg_black_y / avg_white_y;
    double ab = avg_white_g * yd;
    if (fabs(log2(avg_white_y / avg_black_y) - log2(avg_white_g / avg_black_g)) > 0.5) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: large dynamic range difference detected. Likely glare issue.\n");
    }
    elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
        "Camera G on darkest patch(es) is %.1f%% lighter compared to observer Y.\n"
        "  Y dynamic range is %.2f stops, G dynamic range is %.2f stops, difference\n"
        "  %.2f stops. A small difference is normal, while a large indicates that there\n"
        "  is glare.\n",
        100 * (avg_black_g / ab - 1.0), log2(avg_white_y / avg_black_y),
        log2(avg_white_g / avg_black_g),
        log2(avg_white_y / avg_black_y) - log2(avg_white_g / avg_black_g));
    gd->avg_black_g = avg_black_g;
    gd->avg_black_y = avg_black_y;
    gd->avg_white_g = avg_white_g;
    gd->avg_white_y = avg_white_y;
}

void
glare_add(struct patch_set *ps,
          v3 wp,
          double mul)
{
    double max_y = ps->patch[0].xyz.v[1];
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].xyz.v[1] > max_y) max_y = ps->patch[i].xyz.v[1];
    }
    m3x3 xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, wp);
    for (int i = 0; i < ps->patch_count; i++) {
        v3 rgb = m3x3_mul_v3(xyz2rgb, ps->patch[i].xyz);
        for (int i = 0; i < 3; i++) {
            rgb.v[i] += mul * max_y;
        }
        ps->patch[i].xyz = m3x3_mul_v3(m3x3_invert(xyz2rgb), rgb);
    }
}

double
flatfield_correct_targets(const struct test_chart_spec *tcs,
                          struct patch_set *ps,
                          struct patch_set *ps_ex,
                          bool verbose,
                          bool exit_on_failure,
                          bool *bad_correction)
{
    v3 avg_white = {{ 0,0,0 }};
    v3 avg_white_xyz = {{ 0,0,0 }};
    spectrum_t *avg_white_spec = NULL;
    int spec_count = 0;

    if (bad_correction != NULL) *bad_correction = false;

    for (int i = 0; i < tcs->white_count; i++) {
        int idx = tcs->white[i].idx;
        avg_white = v3_add(avg_white, ps->patch[idx].cam);
        avg_white_xyz = v3_add(avg_white_xyz, ps->patch[idx].xyz);
        if (ps->patch[idx].spectrum != NULL) {
            spectrum_t *s;
            if (avg_white_spec == NULL) {
                s = spec_copy(ps->patch[idx].spectrum);
            } else {
                s = spec_add(avg_white_spec, ps->patch[idx].spectrum);
            }
            free(avg_white_spec);
            avg_white_spec = s;
            spec_count++;
        }
    }
    if (spec_count > 0) {
        spec_scalar_multiply(avg_white_spec, 1.0 / spec_count);
    }
    avg_white = k_mul_v3(1.0 / tcs->white_count, avg_white);
    avg_white_xyz = k_mul_v3(1.0 / tcs->white_count, avg_white_xyz);
    v3 cp[3][tcs->white_count];
    for (int i = 0; i < tcs->white_count; i++) {
        int idx = tcs->white[i].idx;
        if (idx < 0 || idx >= ps->patch_count) {
            if (verbose || exit_on_failure) {
                elog("Error: %d is out of range (patch count %d).\n", idx, ps->patch_count);
            }
            if (exit_on_failure) {
                abort();
            } else {
                if (bad_correction != NULL) *bad_correction = true;
                free(avg_white_spec);
                return 0;
            }
        }
        cp[0][i] = cp[1][i] = cp[2][i] = (v3){{ tcs->white[i].r, 0, tcs->white[i].c }};
        cp[0][i].v[1] = avg_white.v[0] / ps->patch[idx].cam.v[0];
        cp[1][i].v[1] = avg_white.v[1] / ps->patch[idx].cam.v[1];
        cp[2][i].v[1] = avg_white.v[2] / ps->patch[idx].cam.v[2];
        //elog("%f %f %d\n", tcs->white[i].r, tcs->white[i].c, idx);
        //elog("%f %f %f\n", cp[0][i].v[1], cp[1][i].v[1], cp[2][i].v[1]);
        ps->patch[idx].xyz = avg_white_xyz;
        if (ps->patch[idx].spectrum != NULL) {
            free(ps->patch[idx].spectrum);
            ps->patch[idx].spectrum = spec_copy(avg_white_spec);
        }
    }
    free(avg_white_spec);
    tps_t *tps[3] = { NULL, NULL, NULL };
    if (tcs->white_count < 3) {
        if (verbose) elog("Warning: too few whites for flatfield correction (skipping).\n");
    } else {
        for (int i = 0; i < 3; i++) {
            tps[i] = tps_new(cp[i], tcs->white_count, 0);
            if (tps[i] == NULL) {
                tps[i] = tps_new(cp[i], tcs->white_count, 0.000000000001);
                if (tps[i] == NULL) {
                    if (exit_on_failure || verbose) {
                        elog("Error: could not make TPS.\n");
                    }
                    if (exit_on_failure) {
                        exit(EXIT_FAILURE);
                    } else {
                        tps_delete(tps[0]);
                        tps_delete(tps[1]);
                        tps_delete(tps[2]);
                        if (bad_correction != NULL) *bad_correction = true;
                        return 0;
                    }
                }
            }
        }
    }
    double max_dist = 0;
    double max_corr = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        double min_dist = 100000;
        for (int j = 0; j < tcs->white_count; j++) {
            double dist = test_chart_dist(&tcs->layout, i, tcs->white[j].idx);
            if (dist < min_dist) min_dist = dist;
        }
        if (min_dist > max_dist) max_dist = min_dist;

        double r, c;
        test_chart_rc(&tcs->layout, i, &r, &c);
        for (int j = 0; j < 3; j++) {
            double corr = tps[j] == NULL ? 1.0 : tps_interpolate(tps[j], r, c);
            if (corr < 0 || corr > 10000) {
                if (exit_on_failure || verbose) {
                    elog("Error: TPS returned bad correction value (%f).\n", corr);
                }
                if (exit_on_failure) {
                    exit(EXIT_FAILURE);
                } else {
                    tps_delete(tps[0]);
                    tps_delete(tps[1]);
                    tps_delete(tps[2]);
                    if (bad_correction != NULL) *bad_correction = true;
                    return 0;
                }
            }
            ps->patch[i].cam.v[j] *= corr;
            if (ps_ex) ps_ex->patch[i].cam.v[j] *= corr;
            if (corr < 1.0) corr = 1.0 / corr;
            if (corr > max_corr) max_corr = corr;
        }
        if (ps->patch[i].spectrum != NULL) {
            //ps->patch[i].xyz = spec2xyz_ill(ps->patch[i].spectrum, &dcam_ssf->ssf, illuminant, 1);
        }
    }
    tps_delete(tps[0]);
    tps_delete(tps[1]);
    tps_delete(tps[2]);
    if (verbose) {
        elog("  Max distance to nearest white (center to center): %.2f patch diagonals.\n"
             "  Largest flatfield correction %.2f%%.\n",
             max_dist / sqrt(tcs->layout.row_height*tcs->layout.row_height + tcs->layout.col_width*tcs->layout.col_width), 100.0 * (max_corr - 1.0));
    }
    return 100.0 * (max_corr - 1.0);
}

bool
glare_match(const struct test_chart_spec *tcs, // should contain neutrals
            struct patch_set *ps, // should be flatfield corrected
            const struct observer *obs,
            const spectrum_t *illuminant,
            const v3 *illuminant_wp,
            const char report_dir[],
            bool verbose,
            bool exit_on_failure)
{
    if (tcs->patch_count != ps->patch_count) {
        if (exit_on_failure || verbose) {
            elog("Error: target layout patch count not same as .ti3.\n");
        }
        if (exit_on_failure) {
            exit(EXIT_FAILURE);
        } else {
            return false;
        }
    }

    v3 wp;
    if (illuminant_wp == NULL || (illuminant != NULL && obs != NULL)) {
        if (illuminant == NULL || obs == NULL) {
            if (exit_on_failure || verbose) {
                elog("Error: illuminant required to model glare.\n");
            }
            if (exit_on_failure) {
                exit(EXIT_FAILURE);
            } else {
                return false;
            }
        }
        wp = v3_norm(spec2xyz(illuminant, obs));
    } else {
        wp = *illuminant_wp;
    }

    // find neutral patches
    const int neutrals_count = tcs->white_count + tcs->black_count + tcs->gray_count;
    struct {
        int idx;
        double lum_g;
        double lum_y;
        double r, c;
        bool is_white;
        bool is_black;
    } neutrals[neutrals_count];
    const int white_count = tcs->white_count;
    const int black_count = tcs->black_count;
    {
        double lightest_neutral_y = 0;
        double lightest_neutral_g = 0;
        for (int i = 0; i < neutrals_count; i++) {
            neutrals[i].idx = tcs->data[i].idx;
            neutrals[i].r = tcs->data[i].r;
            neutrals[i].c = tcs->data[i].c;
            assert(tcs->data[i].idx >= 0 && tcs->data[i].idx < ps->patch_count);
            neutrals[i].lum_g = ps->patch[tcs->data[i].idx].cam.v[1];
            neutrals[i].lum_y = ps->patch[tcs->data[i].idx].xyz.v[1];
            if (neutrals[i].lum_g > lightest_neutral_g) lightest_neutral_g = neutrals[i].lum_g;
            if (neutrals[i].lum_y > lightest_neutral_y) lightest_neutral_y = neutrals[i].lum_y;
            neutrals[i].is_white = tcs->data[i].is_white;
            neutrals[i].is_black = tcs->data[i].is_black;
        }
        for (int i = 0; i < neutrals_count; i++) {
            neutrals[i].lum_g /= lightest_neutral_g;
            neutrals[i].lum_y /= lightest_neutral_y;
        }
    }
    if (white_count == 0) {
        if (exit_on_failure || verbose) {
            elog("Error: no pure white patches found in chart, required for glare modeling.\n");
        }
        if (exit_on_failure) {
            exit(EXIT_FAILURE);
        } else {
            return false;
        }
    }

    struct {
        double handles[2][3];
        int handle_count;
        double *tc;
        int tc_len;
        double r, c;
    } curves[neutrals_count];
    int curves_count = 0;
    const double gamma = 2.2;
    const double patch_diag = sqrt(tcs->layout.row_height*tcs->layout.row_height + tcs->layout.col_width*tcs->layout.col_width);
    for (int i = 0; i < neutrals_count; i++) {
        int white_idx = -1, black_idx = -1;
        double endpoint_dist = 0;
        if (white_count > black_count) {
            if (!neutrals[i].is_white) {
                continue;
            }
            // find closest black
            white_idx = i;
            double black_dist = 1e9;
            for (int j = 0; j < neutrals_count; j++) {
                if (i == j || !neutrals[j].is_black) continue;
                double dist = test_chart_dist(&tcs->layout, neutrals[i].idx, neutrals[j].idx);
                if (dist < black_dist) {
                    black_dist = dist;
                    black_idx = j;
                }
            }
            endpoint_dist = black_dist;
        } else {
            if (!neutrals[i].is_black) {
                continue;
            }
            // find closest white
            black_idx = i;
            double white_dist = 1e9;
            for (int j = 0; j < neutrals_count; j++) {
                if (i == j || !neutrals[j].is_white) continue;
                double dist = test_chart_dist(&tcs->layout, neutrals[i].idx, neutrals[j].idx);
                if (dist < white_dist) {
                    white_dist = dist;
                    white_idx = j;
                }
            }
            endpoint_dist = white_dist;
        }
        if (white_idx == -1 || black_idx == -1) {
            continue;
        }

        // find gray patches close (enough) to closest white, if any

        // We could possibly get more handles if the target contains a detailed step wedge, but using only one middle gray gives us a more smooth and robust model
        struct {
            double lum_y;
            double lum_g;
            double dist;
            int idx;
        } grays[neutrals_count];
        int gray_count = 0;
        for (int j = 0; j < neutrals_count; j++) {
            if (i == j || neutrals[j].is_white || neutrals[j].is_black) continue;
            double dist = test_chart_dist(&tcs->layout, neutrals[i].idx, neutrals[j].idx);
            if (dist < endpoint_dist + 3 * patch_diag) {
                grays[gray_count].lum_y = neutrals[j].lum_y;
                grays[gray_count].lum_g = neutrals[j].lum_g;
                grays[gray_count].dist = dist;
                grays[gray_count].idx = neutrals[j].idx;
                gray_count++;
            }
        }
        // find the gray patch closest to middle gray
        int gray_idx = -1;
        double middle_gray = neutrals[black_idx].lum_y + 0.25 * (neutrals[white_idx].lum_y - neutrals[black_idx].lum_y);
        double min_diff = 1e9;
        for (int j = 0; j < gray_count; j++) {
            double lum_y = grays[j].lum_y;
            double diff = fabs(pow(lum_y, 1/gamma) - pow(middle_gray, 1/gamma));
            // weighting: prefer close patches over those further away
            diff *= grays[j].dist / patch_diag;
            if (tcs->gray_specified) {
                // specified grays, then we don't look if it's mid-gray, we only take the closest
                diff = grays[j].dist;
            }
            if (diff < min_diff && lum_y > neutrals[black_idx].lum_y && lum_y < neutrals[white_idx].lum_y) {
                min_diff = diff;
                gray_idx = j;
            }
        }

        {
            { // assign handles
                curves[curves_count].handles[0][0] = pow(neutrals[black_idx].lum_y, 1/gamma);
                curves[curves_count].handles[1][0] = neutrals[black_idx].lum_g - neutrals[black_idx].lum_y;
                int k;
                if (gray_idx != -1) {
                    curves[curves_count].handles[0][1] = pow(grays[gray_idx].lum_y, 1/gamma);
                    curves[curves_count].handles[1][1] = grays[gray_idx].lum_g - grays[gray_idx].lum_y;
                    curves[curves_count].handle_count = 3;
                    k = 2;
                } else {
                    curves[curves_count].handle_count = 2;
                    curves[curves_count].handles[0][2] = -1;
                    curves[curves_count].handles[1][2] = -1;
                    k = 1;
                }
                curves[curves_count].handles[0][k] = pow(neutrals[white_idx].lum_y, 1/gamma);
                curves[curves_count].handles[1][k] = neutrals[white_idx].lum_g - neutrals[white_idx].lum_y;
            }

            { // calculate curve position in chart as average of its handles positions
                double sum_r = 0, sum_c = 0, r, c;
                int cc = 0;
                if (gray_idx != -1) {
                    test_chart_rc(&tcs->layout, grays[gray_idx].idx, &r, &c);
                    sum_r += r;
                    sum_c += c;
                    cc++;
                }
                test_chart_rc(&tcs->layout, neutrals[white_idx].idx, &r, &c);
                sum_r += r;
                sum_c += c;
                cc++;
                test_chart_rc(&tcs->layout, neutrals[black_idx].idx, &r, &c);
                sum_r += r;
                sum_c += c;
                cc++;
                curves[curves_count].r = sum_r / cc;
                curves[curves_count].c = sum_c / cc;
            }

            if (0) {
                elog("curve %2d: %4.1f %4.1f | %f %f | %f %f | %f %f\n",
                     curves_count,
                     curves[curves_count].r,
                     curves[curves_count].c,
                     curves[curves_count].handles[0][0],
                     curves[curves_count].handles[1][0],
                     curves[curves_count].handles[0][1],
                     curves[curves_count].handles[1][1],
                     curves[curves_count].handles[0][2],
                     curves[curves_count].handles[1][2]);
            }
            curves_count++;
        }
    }

    if (curves_count == 0) {
        if (exit_on_failure || verbose) {
            elog("Error: did not find any handles for glare curves.\n");
        }
        if (exit_on_failure) {
            exit(EXIT_FAILURE);
        } else {
            return false;
        }
    }

    { // generate glare model curves

        /*
          Glare curve:
            - How many percent Y to add depending on input Y
            - Pre-defined shape inspired by actual CCSG shots
            - Always ends and 0% at max (white)
            - The pre-defined shape can only be adjusted if we have middle gray handles
            - White handles are ignored, as the model always ends at 0% for robustness

         */

        const int tc_len = 1024;
        double *outx = malloc(tc_len * sizeof(outx[0]));
        for (int i = 0; i < tc_len; i++) {
            outx[i] = (double)i / (tc_len - 1);
        }

        // with the overall curve as guide, make individual curves based on the handles
        FILE *stream = open_log_file(report_dir, "glare-curves.dat");
        for (int i = 0; i < curves_count; i++) {
            double *tc = malloc(sizeof(tc[0]) * tc_len);
            double inx[5];
            double iny[5];
            inx[0] = 0;
            iny[0] = curves[i].handles[1][0];
            for (int j = 0; j < curves[i].handle_count; j++) {
                inx[j+1] = curves[i].handles[0][j];
                iny[j+1] = curves[i].handles[1][j];
            }
            // lock last (white) to 1/0
            inx[curves[i].handle_count] = 1;
            iny[curves[i].handle_count] = 0;
            rounded_linear_interpolation(inx, iny, 1+curves[i].handle_count, outx, tc, tc_len);
            if (stream) {
                for (int j = 0; j < tc_len; j++) {
                    fprintf(stream, "%f %f\n", outx[j], tc[j]);
                }
                fprintf(stream, "\n");
            }
            curves[i].tc = tc;
            curves[i].tc_len = tc_len;
        }
        if (stream) {
            fclose(stream);
        }
        free(outx);
    }

    struct patch_error_report *pe = malloc(sizeof(*pe) + sizeof(pe->e[0]) * ps->patch_count);
    memset(pe, 0, sizeof(*pe));
    pe->patch_count = ps->patch_count;
    pe->ps = ps;
    bool spectrum_model_used = false;
    bool bad = false;
    { // apply glare to XYZ reference values
#pragma omp parallel for
        for (int i = 0; i < ps->patch_count; i++) {
            v3 cxyz;
            double r, c;
            test_chart_rc(&tcs->layout, i, &r, &c);
            if (ps->patch[i].spectrum == NULL) {
                // apply glare model in ProPhoto RGB space
                m3x3 xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, wp);
                v3 rgb = m3x3_mul_v3(xyz2rgb, ps->patch[i].xyz);
                for (int k = 0; k < 3; k++) {
                    // generate TPS with anchors in all glare curves
                    v3 cp[curves_count];
                    for (int j = 0; j < curves_count; j++) {
                        cp[j].v[0] = curves[j].r;
                        cp[j].v[1] = curve_y(pow(rgb.v[k], 1/gamma), curves[j].tc, curves[j].tc_len);
                        cp[j].v[2] = curves[j].c;
                    }
                    double y;
                    if (curves_count >= 3) {
                        tps_t *tps = tps_new(cp, curves_count,  0.000000000001);
                        if (tps == NULL) {
                            if ((verbose || exit_on_failure) && !bad) {
                                elog("Error: failed to create TPS.\n");
                            }
                            bad = true;
                            y = 0;
                        } else {
                            y = tps_interpolate(tps, r, c);
                            tps_delete(tps);
                        }
                    } else if (curves_count == 2) {
                        double d1 = sqrt((r-cp[0].v[0])*(r-cp[0].v[0]) + (c-cp[0].v[2])*(c-cp[0].v[2]));
                        double d2 = sqrt((r-cp[1].v[0])*(r-cp[1].v[0]) + (c-cp[1].v[2])*(c-cp[1].v[2]));
                        double y1 = cp[0].v[1];
                        double y2 = cp[1].v[1];
                        // longer distance, less weight so swap places on distances for right weight
                        y = d2 * y1 + d1 * y2 / (d1 + d2);
                    } else {
                        y = cp[0].v[1];
                    }
                    rgb.v[k] += y;
                }
                // use glare model's chroma and lightness, but keep original hue for robustness
                v3 jch = xyz2jch(m3x3_mul_v3(m3x3_invert(xyz2rgb), rgb), wp);
                if (!v3_isfinite(jch)) {
                    if ((verbose || exit_on_failure) && !bad) {
                        elog("Error: JCh value out of range. Too saturated reference XYZ values in target?\n");
                    }
                    bad = true;
                    exit(EXIT_FAILURE);
                }
                jch.v[2] = xyz2jch(ps->patch[i].xyz, wp).v[2];
                cxyz = jch2xyz(jch, wp);
            } else {
                spectrum_model_used = true;
                if (illuminant == NULL || obs == NULL) {
                    if ((verbose || exit_on_failure) && !bad) {
                        elog("Error: observer and illuminant required to model glare as the target has spectra.\n");
                    }
                    bad = true;
                    exit(EXIT_FAILURE);
                }
                ps->patch[i].xyz = spec2xyz_ill(ps->patch[i].spectrum, obs, illuminant, 1);
                double x = ps->patch[i].xyz.v[1];

                // generate TPS with anchors in all glare curves
                v3 cp[curves_count];
                for (int j = 0; j < curves_count; j++) {
                    cp[j].v[0] = curves[j].r;
                    cp[j].v[1] = curve_y(pow(x, 1/gamma), curves[j].tc, curves[j].tc_len);
                    cp[j].v[2] = curves[j].c;
                }
                double y;
                if (curves_count >= 3) {
                    tps_t *tps = tps_new(cp, curves_count,  0.000000000001);
                    if (tps == NULL) {
                        if ((verbose || exit_on_failure) && !bad) {
                            elog("Error: failed to create TPS.\n");
                        }
                        bad = true;
                        y = 0;
                    } else {
                        y = tps_interpolate(tps, r, c);
                        tps_delete(tps);
                    }
                } else if (curves_count == 2) {
                    double d1 = sqrt((r-cp[0].v[0])*(r-cp[0].v[0]) + (c-cp[0].v[2])*(c-cp[0].v[2]));
                    double d2 = sqrt((r-cp[1].v[0])*(r-cp[1].v[0]) + (c-cp[1].v[2])*(c-cp[1].v[2]));
                    double y1 = cp[0].v[1];
                    double y2 = cp[1].v[1];
                    // longer distance, less weight so swap places on distances for right weight
                    y = d2 * y1 + d1 * y2 / (d1 + d2);
                } else {
                    y = cp[0].v[1];
                }

                spectrum_t *s = ps->patch[i].spectrum;
                for (int j = 0; j < spec_band_count(s); j++) {
                    s->v[j].e += y;
                }

                cxyz = spec2xyz_ill(ps->patch[i].spectrum, obs, illuminant, 1);
            }
            const v3 xyz = ps->patch[i].xyz;
            pe->e[i].idx = i;
            pe->e[i].cxyz = cxyz;
            pe->e[i].errors[0] = ciede2000(xyz, cxyz, wp);
            pe->e[i].errors[1] = ciede2000_lch(xyz, cxyz, wp, 1, 100000, 100000);
            if (xyz2lab(xyz, wp).v[0] > xyz2lab(cxyz, wp).v[0]) {
                pe->e[i].errors[1] = -pe->e[i].errors[1];
            }
            pe->e[i].errors[2] = ciede2000_lch(xyz, cxyz, wp, 100000, 1, 100000);
            pe->e[i].errors[3] = ciede2000_lch(xyz, cxyz, wp, 100000, 100000, 1);
        }
    }
    if (report_dir) {
        render_pm_error_image(report_dir, "glare-match.tif", wp, pe, false);
    }

    if (!bad) {
        double min_y_old = 1e9;
        double min_y = 1e9;
        for (int i = 0; i < ps->patch_count; i++) {
            if (ps->patch[i].xyz.v[1] < min_y_old) min_y_old = ps->patch[i].xyz.v[1];
            ps->patch[i].xyz = pe->e[i].cxyz;
            if (ps->patch[i].xyz.v[1] < min_y) min_y = ps->patch[i].xyz.v[1];
        }
        if (verbose) {
            elog("  Minimum Y changed from %f to %f.%s\n", min_y_old, min_y,
                 spectrum_model_used ? " Glare was modeled into the spectra." : " Glare was modeled in RGB space.");
        }
    } else if (exit_on_failure) {
        exit(EXIT_FAILURE);
    }

    free(pe);
    for (int j = 0; j < curves_count; j++) {
        free(curves[j].tc);
    }

    return !bad;
}
