/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <strings.h>
#include <ctype.h>

#include <lcms2.h>

#include <nmsimplex.h>
#include <bisection.h>
#include <interp.h>
#include <target.h>
#include <colmath.h>
#include <matmath.h>
#include <matopt.h>
#include <profio.h>
#include <target.h>
#include <gamut.h>
#include <glare.h>
#include <spectraldb.h>
#include <xyz2spec.h>
#include <observers.h>
#include <argyllio.h>
#include <jsonio.h>
#include <dcamprof.h>
#include <wcompat.h>
#include <dngref.h>
#include <icclut.h>
#include <dnglut.h>
#include <tifio.h>
#include <look.h>
#include <elog.h>
#include <lut.h>
#include <tps.h>

#define UNUSED(x) (void)(x)
#define debug(...)
//#define debug(...) elog(__VA_ARGS__)

static struct {
    v3 d50;
    int pe_cmp_idx;
} glob;

FILE *
open_log_file(const char dir[],
              const char name[])
{
    if (dir == NULL) {
        return NULL;
    }
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    mkdir(dir);
#else
    mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);
#endif
    char path[strlen(dir) + strlen(name) + 32];
    sprintf(path, "%s/%s", dir, name);
    FILE *stream = utf8_fopen(path, "w");
    if (stream == NULL) {
        elog("Error: could not open \"%s\" for writing: %s.\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return stream;
}

struct find_distant_light_color_simplex_fun_arg {
    const struct patch_set *ps;
    m3x3 rgb2xyz;
    v3 wp;
};

static double
find_distant_light_color_simplex_fun(double m[],
                                     void *arg)
{
    //elog("%f %f\n", m[0], m[1]);
    if (m[0] < 0 || m[0] > 1 || m[1] < 0.0 || m[1] > 1.0) {
        return 1000000;
    }
    v3 hsl = {{ m[0], 1.0, 0.5 + 0.5 * m[1] }};
    v3 rgb = hsl2rgb(hsl);
    struct find_distant_light_color_simplex_fun_arg *a = (struct find_distant_light_color_simplex_fun_arg *)arg;

    v3 xyz = m3x3_mul_v3(a->rgb2xyz, rgb);
    double min_de = 1000000;
    for (int i = 0; i < a->ps->patch_count; i++) {
        double de = ciede2000_lch(xyz, a->ps->patch[i].xyz, a->wp, 1000000, 1, 1);
        //double de = euclidian_dist(rgb, a->ps->patch[i].cam);
        //double de = euclidian_dist(hsl, rgb2hsl(a->ps->patch[i].cam));
        if (de < min_de) min_de = de;
    }
    //elog("%f %f %f\n", min_de, m[0], m[1]);
    return 1.0 / min_de;
}

static v3
find_distant_light_color(const struct patch_set *ps,
                         m3x3 rgb2xyz)
{
    struct find_distant_light_color_simplex_fun_arg arg;
    arg.ps = ps;
    arg.rgb2xyz = rgb2xyz;
    arg.wp = m3x3_mul_v3(rgb2xyz, (v3){{1,1,1}});
    int hcount = 100;
    int scount = 100;
    double max_de = 0;
    v3 max_hsl = {{0,0,0}};
#pragma omp parallel for
    for (int h = 0; h < hcount; h++) {
        for (int s = scount - 1; s >= 0; s--) {
            v3 hsl = {{ (double)h / hcount, 1.0, 0.5 + 0.5 * (double)s / scount }};
            v3 rgb = hsl2rgb(hsl);
            v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
            double min_de = 1000000;
            for (int i = 0; i < ps->patch_count; i++) {
                double de = ciede2000_lch(xyz, ps->patch[i].xyz, arg.wp, 1000000, 1, 1);
                //double de = euclidian_dist(rgb, ps->patch[i].cam);
                //double de = euclidian_dist(hsl, rgb2hsl(ps->patch[i].cam));
                if (de < min_de) min_de = de;
            }
#pragma omp critical
            if (min_de > max_de) {
                max_de = min_de;
                max_hsl = hsl;
            }
        }
    }
    double m[2] = { max_hsl.v[0], 2.0 * (max_hsl.v[2] - 0.5) };
    double min = 1000, old_min;
    do {
        old_min = min;
        min = simplex(find_distant_light_color_simplex_fun, &arg, m, 2, 1.0e-7, 0.1, NULL);
    } while (min < old_min);
    //elog("%f %f %f %f\n", max_hsl.v[0], max_hsl.v[2], 1.0 / max_de, min);
    return hsl2rgb((v3){{ m[0], 1.0, 0.5 + 0.5 * m[1] }});
}

static struct patch_set *
exclude_darker_patches(struct patch_set *ps,
                       double uv_de_lim)
{
    double max_y = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].xyz.v[1] > max_y) {
            max_y = ps->patch[i].xyz.v[1];
        }
    }
    for (int i = 0; i < ps->patch_count; i++) {
        for (int j = 0; j < ps->patch_count; j++) {
            if (i == j) continue;
            v3 xyz1 = ps->patch[i].xyz;
            v3 xyz2 = ps->patch[j].xyz;
            for (int k = 0; k < 3; k++) {
                if (xyz1.v[k] != 0) xyz1.v[k] /= max_y;
                if (xyz2.v[k] != 0) xyz2.v[k] /= max_y;
            }
            v3 ltu1 = xyz2lutspace(xyz1);
            v3 ltu2 = xyz2lutspace(xyz2);
            ltu1.v[0] = ltu2.v[0];
            double uv_de = euclidean_dist(ltu1, ltu2);
            if (uv_de < uv_de_lim) {
                if (xyz1.v[0]+xyz1.v[1]+xyz1.v[2] < xyz2.v[0]+xyz2.v[1]+xyz2.v[2]) {
                    ps->patch[i].class = -1;
                } else {
                    ps->patch[j].class = -1;
                }
            }
        }
    }
    int orig_count = ps->patch_count;
    ps->patch_count = 0;
    for (int i = 0; i < orig_count; i++) {
        if (ps->patch[i].class >= 0) {
            if (i != ps->patch_count) {
                ps->patch[ps->patch_count] = ps->patch[i];
            }
            ps->patch_count++;
        }
    }
    if (ps->patch_count < orig_count) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Excluded %d patches from target due to lighter patch with same chroma.\n"
            "  %d patches left\n",
            orig_count - ps->patch_count, ps->patch_count);
    }
    return ps;
}

static struct patch_set *
exclude_weak_signal_patches(struct patch_set *ps,
                            bool verbose)
{
    double max_cam = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        for (int j = 0; j < 3; j++) {
            if (ps->patch[i].cam.v[j] > max_cam) {
                max_cam = ps->patch[i].cam.v[j];
            }
        }
    }
    int orig_count = ps->patch_count;
    ps->patch_count = 0;
    for (int i = 0; i < orig_count; i++) {
        double cam[3];
        cam[0] = ps->patch[i].cam.v[0] / max_cam;
        cam[1] = ps->patch[i].cam.v[1] / max_cam;
        cam[2] = ps->patch[i].cam.v[2] / max_cam;
        if (cam[0] >= 0.0001 || cam[1] >= 0.0001 || cam[2] >= 0.0001) {
            if (i != ps->patch_count) {
                ps->patch[ps->patch_count] = ps->patch[i];
            }
            ps->patch_count++;
        }
    }
    if (ps->patch_count < orig_count && verbose) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Excluded %d patches from target due to weak camera signal. %d patches left\n", orig_count - ps->patch_count, ps->patch_count);
    }
    return ps;
}

static void
assign_optimization_weights(struct patch_set *ps,
                            const double class_mtx_weight[],
                            const double (*class_lut_de_w)[7],
                            const double patch_mtx_weight[],
                            const double (*patch_lut_de_w)[7],
                            const double (*class_mtx_de_w)[7],
                            const double (*patch_mtx_de_w)[7],
                            bool skip_base_weighting,
                            const v3 illuminant_xyz)
{
    double wsum = 0;
#pragma omp parallel for
    for (int i = 0; i < ps->patch_count; i++) {
        if (!skip_base_weighting) {
            double de = 0.5; // start with some minimum
            for (int j = 0; j < ps->patch_count; j++) {
                if (i == j) continue;
                if (class_mtx_weight[ps->patch[j].class] == 0) continue;
                de += ciede2000_lch(ps->patch[j].xyz, ps->patch[i].xyz, illuminant_xyz, 1, 1, 1);
            }
            ps->patch[i].mtx_w = de;
        } else {
            ps->patch[i].mtx_w = 1.0;
        }
        ps->patch[i].mtx_w *= class_mtx_weight[ps->patch[i].class];
        ps->patch[i].mtx_w *= patch_mtx_weight[i];
#pragma omp critical
        {
            wsum += ps->patch[i].mtx_w;
        }
        {
            // auto: no lightness correction, and chroma-dependent hue/chroma correction
            const double chroma = xyz2jch(ps->patch[i].xyz, illuminant_xyz).v[1];
            const double lo_lim = 30;
            const double hi_lim = 70;
            double x;
            if (chroma < lo_lim) {
                x = 0;
            } else if (chroma < hi_lim) {
                x = (chroma - lo_lim) / (hi_lim - lo_lim);
            } else {
                x = 1;
            }
            ps->patch[i].lut_de_w[0] = -1000;
            ps->patch[i].lut_de_w[1] = +1000;
            ps->patch[i].lut_de_w[2] = -1.0*(1-x) + -4.0*x;
            ps->patch[i].lut_de_w[3] = 0.7*(1-x) + 1.5*x;
            ps->patch[i].lut_de_w[4] = -(0.5*(1-x) + 2.0*x);
            ps->patch[i].lut_de_w[5] = -ps->patch[i].lut_de_w[4];
            /*
            char buf[64];
            elog("%f %s %s %f %f %f\n", chroma, ps->patch[i].str_id, xyz2name(ps->patch[i].xyz, illuminant_xyz, buf),
                 ps->patch[i].de_w[2], ps->patch[i].de_w[3], ps->patch[i].de_w[4]);
            */
        }
        if (class_lut_de_w != NULL && class_lut_de_w[ps->patch[i].class][6] != 0) {
            for (int j = 0; j < 6; j++) ps->patch[i].lut_de_w[j] = class_lut_de_w[ps->patch[i].class][j];
        }
        if (patch_lut_de_w != NULL && patch_lut_de_w[i][6] != 0) {
            for (int j = 0; j < 6; j++) ps->patch[i].lut_de_w[j] = patch_lut_de_w[i][j];
        }
        ps->patch[i].use_mtx_de_w = false;
        if (class_mtx_de_w != NULL && class_mtx_de_w[ps->patch[i].class][6] != 0) {
            ps->patch[i].use_mtx_de_w = true;
            for (int j = 0; j < 6; j++) ps->patch[i].mtx_de_w[j] = class_mtx_de_w[ps->patch[i].class][j];
        }
        if (patch_mtx_de_w != NULL && patch_mtx_de_w[i][6] != 0) {
            ps->patch[i].use_mtx_de_w = true;
            for (int j = 0; j < 6; j++) ps->patch[i].mtx_de_w[j] = patch_mtx_de_w[i][j];
        }
    }
    if (wsum != 0) {
        for (int i = 0; i < ps->patch_count; i++) ps->patch[i].mtx_w /= wsum;
    }
}

static int
find_whitest_patch(const struct patch_set *ps,
                   v3 wp,
                   bool *found_lightest)
{
    int white_class = -1;
    int w_idx = -1;
    for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
        if (strcasecmp(ps->class_names[i], "white") == 0) {
            white_class = i;
            break;
        }
    }
    if (white_class != -1) {
        for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
            if (strcasecmp(ps->class_names[i], "illuminant") == 0) {
                white_class = i;
                break;
            }
        }
    }
    if (white_class != -1) {
        for (int i = 0; i < ps->patch_count; i++) {
            if (ps->patch[i].class == white_class) {
                w_idx = i;
                break;
            }
        }
    } else {
        bool *cand = calloc(1, sizeof(*cand) * ps->patch_count);
        bool has_cand = false;
        for (int i = 0; i < ps->patch_count; i++) {
            v3 xyz = ps->patch[i].xyz;
            if (xyz.v[1] < 0.1) {  // skip very dark patches
                continue;
            }
            v3 lab = xyz2lab(v3_norm(xyz), wp);
            double ab = sqrt(lab.v[1] * lab.v[1] + lab.v[2] * lab.v[2]);
            if (ab <= 2.0) {
                has_cand = true;
                cand[i] = true;
            }
        }
        double min_ab = 100000000;
        double min_ab_all = 100000000;
        double max_y = 0;
        int w_idx_all = -1;
        for (int i = 0; i < ps->patch_count; i++) {
            v3 xyz = ps->patch[i].xyz;
            if (has_cand) {
                if (!cand[i]) continue;
                if (xyz.v[1] > max_y) {
                    max_y = xyz.v[1];
                    w_idx = i;
                }
            } else {
                v3 lab = xyz2lab(v3_norm(xyz), wp);
                double ab = sqrt(lab.v[1] * lab.v[1] + lab.v[2] * lab.v[2]);
                if (xyz.v[1] >= 0.1) {
                    if (ab < min_ab) {
                        min_ab = ab;
                        w_idx = i;
                    }
                }
                if (ab < min_ab_all) {
                    min_ab_all = ab;
                    w_idx_all = i;
                }
            }
        }
        free(cand);
        if (w_idx == -1) w_idx = w_idx_all;
        if (found_lightest != NULL) {
            *found_lightest = true;
            for (int i = 0; i < ps->patch_count; i++) {
                if (ps->patch[i].xyz.v[1] > ps->patch[w_idx].xyz.v[1]) {
                    v3 lab = xyz2lab(v3_norm(ps->patch[i].xyz), wp);
                    double ab = sqrt(lab.v[1] * lab.v[1] + lab.v[2] * lab.v[2]);
                    if (ab < 6.0) { // reasonably white
                        *found_lightest = false;
                        break;
                    }
                }
            }
        }
    }
    return w_idx;
}

static int
find_neutral_and_rebalance(struct patch_set *ps,
                           const v3 ref_wp,
                           const char neutral_patch[],
                           bool skip_rebalancing)
{
    if (ps == NULL || ps->patch_count == 0) {
        return -1;
    }
    bool found_lightest = true;
    int w_idx = (neutral_patch == NULL) ? find_whitest_patch(ps, ref_wp, &found_lightest) + 1 : atoi(neutral_patch);
    if (w_idx == 0) {
        for (int i = 0; i < ps->patch_count; i++) {
            if (strcmp(ps->patch[i].str_id, neutral_patch) == 0) {
                w_idx = i + 1;
                break;
            }
        }
        if (w_idx == 0) {
            elog("Error: did not find neutral patch named \"%s\" in target.\n", neutral_patch);
            exit(EXIT_FAILURE);
        }
    } else if (w_idx < 0 || w_idx > ps->patch_count) {
        elog("Error: bad neutral patch index %s for target.\n", neutral_patch);
        exit(EXIT_FAILURE);
    }
    w_idx--;
    char w_name[sizeof(ps->patch[w_idx].str_id)+32];
    if (ps->patch[w_idx].str_id[0] == '\0') {
        sprintf(w_name, "index %d", w_idx+1);
    } else {
        strcpy(w_name, ps->patch[w_idx].str_id);
    }
    v3 old_wp = v3_norm(ps->patch[w_idx].xyz);
    double de = ciede2000(old_wp, ref_wp, ref_wp);
    if (!skip_rebalancing && !found_lightest) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: auto-selected neutral patch (%s) is not the lightest, as the\n"
            "  lightest patch is considerably off-white. That is if you later use the target\n"
            "  for white balancing you should use the indicated patch instead of the\n"
            "  lightest.\n",
            w_name);
    }
    if (skip_rebalancing || de == 0.0) {
        return w_idx;
    }
    if (de >= 0.01) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "The most neutral patch (%s) differs %.2f DE from actual neutral,\n"
            "  transforming target reference XYZ values to match, using CAT02.\n",
            w_name, de);
    }
    for (int i = 0; i < ps->patch_count; i++) {
        ps->patch[i].xyz = chromatic_adaptation_transform(CAT_CAT02, ps->patch[i].xyz, old_wp, ref_wp);
    }

    // CAT won't make the whitepoint match exactly, here we make sure that it does:
    double sum[2] = { 0, 0 };
    for (int i = 0; i < 3; i++) {
        sum[0] += ps->patch[w_idx].xyz.v[i];
        sum[1] += ref_wp.v[i];
    }
    ps->patch[w_idx].xyz = k_mul_v3(sum[0] / sum[1], ref_wp);

    return w_idx;
}

static struct patch_set *
normalize_patch_set(const struct patch_set *ps_,
                    int w_idx)
{
    if (ps_ == NULL || ps_->patch_count == 0) {
        return NULL;
    }
    struct patch_set *ps = malloc(sizeof(*ps) + sizeof(struct cam2xyz_patch) * ps_->patch_count);
    memcpy(ps, ps_, sizeof(*ps) + sizeof(struct cam2xyz_patch) * ps_->patch_count);

    double max_rgb = 0;
    double max_xyz = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        for (int j = 0; j < 3; j++) {
            if (ps->patch[i].cam.v[j] > max_rgb) max_rgb = ps->patch[i].cam.v[j];
            if (ps->patch[i].xyz.v[j] > max_xyz) max_xyz = ps->patch[i].xyz.v[j];
        }
    }
    for (int i = 0; i < ps->patch_count; i++) {
        ps->patch[i].spectrum = spec_copy(ps_->patch[i].spectrum);
        for (int j = 0; j < 3; j++) {
            ps->patch[i].xyz.v[j] /= max_xyz;
            ps->patch[i].cam.v[j] /= max_rgb;
        }
    }

    /* make sure that max(cam) == max(xyz) for white, which it probably already is if the white patch is the
       brightest in set, otherwise there's probably a small difference. */
    if (w_idx >= 0) {
        double cmax = v3_max(ps->patch[w_idx].cam);
        double rmax = v3_max(ps->patch[w_idx].xyz);
        if (cmax > rmax) {
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].cam = k_mul_v3(rmax / cmax, ps->patch[i].cam);
            }
        } else {
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].xyz = k_mul_v3(cmax / rmax, ps->patch[i].xyz);
            }
        }
    }
    return ps;
}

static bool
write_dcam_profile(const char filename[],
                   const struct dcam_profile *prof,
                   bool error_is_fatal)
{
    FILE *stream = utf8_fopen(filename, "w");
    if (stream == NULL) {
        elog("Error: could not open \"%s\" for writing: %s.\n", filename, strerror(errno));
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    fprintf(stream, "{\n");
    if (prof->camera_name != NULL) {
        fprintf(stream, "  \"CameraName\": \"%s\",\n", prof->camera_name);
    }
    fprintf(stream, "  \"CalibrationIlluminant1\": \"%s\",\n", exif_lightsource_ntoa(prof->illuminant_name));
    fprintf(stream, "  \"CalibrationIlluminantWhitepoint1\": [ %f, %f, %f ],\n", prof->illuminant_wp.v[0], prof->illuminant_wp.v[1], prof->illuminant_wp.v[2]);
    if (prof->illuminant_spec != NULL) {
        jsonio_spectrum_print(stream, "  ", "CalibrationIlluminantSpectrum1", prof->illuminant_spec);
        fprintf(stream, ",\n");
    }
    jsonio_3x3matrix_print(stream, "  ", "ColorMatrix1", prof->color_matrix);
    fprintf(stream, ",\n");
    fprintf(stream, "  \"ForwardMatrixWhitebalance1\": [ %f, %f, %f ],\n", prof->forward_matrix_wb.v[0], prof->forward_matrix_wb.v[1], prof->forward_matrix_wb.v[2]);
    jsonio_3x3matrix_print(stream, "  ", "ForwardMatrix1", prof->forward_matrix);
    fprintf(stream, ",\n");
    if (!m3x3_isequal(prof->forward_matrix, prof->lut_matrix)) {
        jsonio_3x3matrix_print(stream, "  ", "LUTMatrix1", prof->lut_matrix);
        fprintf(stream, ",\n");
    }
    chromalut_json_print(stream, "  ", "ChromaLUT1", prof->lut);
    fprintf(stream, "\n}\n");
    fclose(stream);

    return true;
}

static cmsHTRANSFORM
icc_lookup_fast_init(const struct profio_icc *icc)
{
    size_t sz;
    void *raw = profio_icc_toraw(icc, &sz);
    cmsHPROFILE h_in = cmsOpenProfileFromMem(raw, sz);
    free(raw);
    cmsHTRANSFORM tfm = cmsCreateTransform(h_in, TYPE_RGB_DBL, NULL, icc->pcs_is_lab ? TYPE_Lab_DBL : TYPE_XYZ_DBL, 0, 0);
    assert(tfm != NULL);
    cmsCloseProfile(h_in);
    return tfm;
}

static v3
icc_lookup_fast(cmsHTRANSFORM tfm,
                bool pcs_is_lab,
                v3 cam)
{
    double f[3] = { cam.v[0], cam.v[1], cam.v[2] };
    cmsDoTransform(tfm, f, f, 1);
    if (pcs_is_lab) {
        return lab2xyz((v3){{f[0],f[1],f[2]}}, icc_d50());
    }
    return (v3){{ f[0], f[1], f[2] }};
}

/*
static v3
icc_lookup(const struct profio_icc *icc,
           v3 cam)
{
    size_t sz;
    void *raw = profio_icc_toraw(icc, &sz);
    cmsHPROFILE h_in = cmsOpenProfileFromMem(raw, sz);
    free(raw);
    cmsHTRANSFORM tfm = cmsCreateTransform(h_in, TYPE_RGB_DBL, NULL, icc->pcs_is_lab ? TYPE_Lab_DBL : TYPE_XYZ_DBL, 0, 0);
    assert(tfm != NULL);
    double f[3] = { cam.v[0], cam.v[1], cam.v[2] };
    cmsDoTransform(tfm, f, f, 1);
    cmsDeleteTransform(tfm);
    cmsCloseProfile(h_in);
    if (icc->pcs_is_lab) {
        return lab2xyz((v3){{f[0],f[1],f[2]}}, icc_d50());
    }
    return (v3){{ f[0], f[1], f[2] }};
}
*/

struct icc_reverse_lookup_simplex_fun_arg {
    cmsHTRANSFORM tfm;
    bool pcs_is_lab;
    v3 xyz;
};

static double
icc_reverse_lookup_simplex_fun(double m[],
                               void *arg)
{
    for (int i = 0; i < 3; i++) {
        if (m[i] < 0 || m[i] > 1) {
            return 100000000;
        }
    }
    struct icc_reverse_lookup_simplex_fun_arg *a = (struct icc_reverse_lookup_simplex_fun_arg *)arg;
    double f[3];
    cmsDoTransform(a->tfm, m, f, 1);
    v3 xyz;
    if (a->pcs_is_lab) {
        xyz = lab2xyz((v3){{f[0],f[1],f[2]}}, icc_d50());
    } else {
        xyz = (v3){{f[0],f[1],f[2]}};
    }
    return euclidean_dist(v3_norm(xyz), v3_norm(a->xyz));
}

static v3
icc_reverse_lookup(const struct profio_icc *icc,
                   v3 xyz)
{
    size_t sz;
    void *raw = profio_icc_toraw(icc, &sz);
    cmsHPROFILE h_in = cmsOpenProfileFromMem(raw, sz);
    assert(h_in != NULL);
    free(raw);

    cmsHTRANSFORM tfm = cmsCreateTransform(h_in, TYPE_RGB_DBL, NULL, icc->pcs_is_lab ? TYPE_Lab_DBL : TYPE_XYZ_DBL, 0, 0);
    assert(tfm != NULL);

    struct icc_reverse_lookup_simplex_fun_arg arg;
    arg.tfm = tfm;
    arg.xyz = xyz;
    arg.pcs_is_lab = icc->pcs_is_lab;
    double m[3] = { 0.5, 0.5, 0.5 };
    double min = simplex(icc_reverse_lookup_simplex_fun, &arg, m, 3, 1.0e-7, 0.01, NULL);
    cmsDeleteTransform(tfm);
    cmsCloseProfile(h_in);
    if (min > 0.0001) {
        return (v3){{ -1, -1, -1 }};
    }
    return (v3){{ m[0], m[1], m[2] }};
}

void
dcam_profile_delete(struct dcam_profile *prof)
{
    if (prof == NULL) {
        return;
    }
    free(prof->illuminant_spec);
    chromalut_delete(prof->lut);
    free(prof->camera_name);
    free(prof);
}

struct dcp_lut {
    v3 *hsm;
    uint32_t *hsmdims;
    bool srgb_gamma;
};

static int
pe_cmp(const void *a_, const void *b_)
{
    const struct patch_error *a = a_;
    const struct patch_error *b = b_;
    if (fabs(a->errors[glob.pe_cmp_idx]) < fabs(b->errors[glob.pe_cmp_idx])) {
        return -1;
    }
    if (fabs(a->errors[glob.pe_cmp_idx]) > fabs(b->errors[glob.pe_cmp_idx])) {
        return 1;
    }
    return 0;
}

static void
render_test_tiff(const char report_dir[],
                 const char output_image_filename[],
                 const m3x3 *cam2xyz,
                 chromalut_t *lut,
                 struct dcp_lut *dcp_lut,
                 struct dcp_lut *dcp_lut_look,
                 const struct profio_icc *icc,
                 const double tc[],
                 int tc_len,
                 bool skip_lut,
                 bool skip_look_lut,
                 double *trc[3],
                 int trc_count,
                 const v3 *wb,
                 bool show_clip,
                 double *clip1,
                 double *clip2,
                 const struct rgbimg *custom_rgbimg)
{
    const double img_gamma = 461.0 / 256.0; // prophoto gamma, about 1.8
    float *p;
    int w, h;
    if (custom_rgbimg != NULL) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: interpretation of input image ICC not implemented, assuming input\n"
            "  image has ProPhoto ICC (linearizing 1.8 gamma).\n");
        w = custom_rgbimg->w;
        h = custom_rgbimg->h;
        p = calloc(1, w * h * 3 * sizeof(p[0]));
#pragma omp parallel for
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < 3; c++) {
                    float v = custom_rgbimg->p[3*(w*y+x)+c] / 65535.0;
                    p[3*(w*y+x)+c] = powf(v, img_gamma); // linearize
                }
            }
        }
    } else {
        const int s_count = 10;
        int w_base = 500;
        int w_extra = 100;
        w = 5 * (w_base + w_extra);
        h = w_base * s_count;
        p = calloc(1, w * h * 3 * sizeof(p[0]));

#pragma omp parallel for
        for (int s = 0; s < s_count; s++) {
            for (int l = 0; l < h/s_count; l++) {
                for (int h_ = 0; h_ < w_base + w_extra; h_++) {
                    const double sat[] = { 0.1, 0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.8, 1.0 };
                    double s_ = sat[s];
                    double l_ = (double)l / (h/s_count-1);
                    l_ = pow(l_, 2.2); // apply a gamma to not make too much focus on highlights
                    double hue = (double)h_ / (w_base-1);
                    bool grayscale = false;
                    if (h_ > w_base + 3 * w_extra / 4) {
                        grayscale = true;
                    } else if (hue > 1) {
                        hue -= 1;
                    }
                    v3 rgb;
                    if (grayscale) {
                        rgb = (v3){{ l_, l_, l_ }};
                    } else {
                        v3 hsl = {{ hue, s_, l_ }};
                        rgb = hsl2rgb(hsl);
                    }
                    int y = s * (h/s_count) + l;
                    p[3*(w*y+h_)+0] = rgb.v[0];
                    p[3*(w*y+h_)+1] = rgb.v[1];
                    p[3*(w*y+h_)+2] = rgb.v[2];
                }
            }
        }
#pragma omp parallel for
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w_base + w_extra; x++) {
                int x1 = x + w_base + w_extra;
                double m;
                m = 0.50;
                p[3*(w*y+x1)+0] = p[3*(w*y+x)+0] * m;
                p[3*(w*y+x1)+1] = p[3*(w*y+x)+1] * m;
                p[3*(w*y+x1)+2] = p[3*(w*y+x)+2] * m;
                x1 += w_base + w_extra;
                m = 0.25;
                p[3*(w*y+x1)+0] = p[3*(w*y+x)+0] * m;
                p[3*(w*y+x1)+1] = p[3*(w*y+x)+1] * m;
                p[3*(w*y+x1)+2] = p[3*(w*y+x)+2] * m;
                x1 += w_base + w_extra;
                m = 0.125;
                p[3*(w*y+x1)+0] = p[3*(w*y+x)+0] * m;
                p[3*(w*y+x1)+1] = p[3*(w*y+x)+1] * m;
                p[3*(w*y+x1)+2] = p[3*(w*y+x)+2] * m;
                x1 += w_base + w_extra;
                m = 0.0625;
                p[3*(w*y+x1)+0] = p[3*(w*y+x)+0] * m;
                p[3*(w*y+x1)+1] = p[3*(w*y+x)+1] * m;
                p[3*(w*y+x1)+2] = p[3*(w*y+x)+2] * m;
            }
        }

    }

    struct rgbimg *img = malloc(sizeof(*img) + w * h * 3 * sizeof(img->p[0]));
    img->w = w;
    img->h = h;

//#warning "ntro testing enabled!"
    const bool test_ntro = false; // compile time developer only option, used when testing neutral tone reproduction operator
#pragma omp parallel for
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            for (int c = 0; c < 3; c++) {
                float v = p[3*(w*y+x)+c];
                if (test_ntro) {
                    v = v * pow(2, 0.7); // +0.7 stops suitable for ACR curve A/B testing
                    if (v > 1) v = 0;
                    v = powf(v, 1/img_gamma);
                    img->p[3*(w*y+x)+c] = (uint16_t)roundf(v * 65535);
                } else {
                    img->p[3*(w*y+x)+c] = (uint16_t)roundf(v * 65535); // apply no gamma, actual RGB sample values in file = actual raw values
                }
            }
        }
    }
    if (report_dir != NULL) {
        char fname[strlen(report_dir) + 100];
        if (custom_rgbimg != NULL) {
            sprintf(fname, "%s/test-image-ref.tif", report_dir);
        } else {
            sprintf(fname, "%s/gradient-ref.tif", report_dir);
        }
        tifio_rgbimg_save(fname, img, test_ntro ? COLORSPACE_PROPHOTO : COLORSPACE_NONE, true);
    }

    if (trc_count > 0) {
        const int out_len = 4096;
        double *out_trc[3] = { NULL, NULL, NULL };
        if (trc_count > 1) {
            for (int j = 0; j < 3; j++) {
                out_trc[j] = malloc(out_len * sizeof(out_trc[j][0]));
#pragma omp parallel for
                for (int i = 0; i < out_len; i++) {
                    out_trc[j][i] = inverse_curve((double)i / (out_len-1), trc[j], trc_count);
                }
            }
        }
#pragma omp parallel for
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < 3; c++) {
                    float v = p[3*(w*y+x)+c];
                    if (trc_count == 1) {
                        v = powf(v, 1.0/trc[c][0]);
                    } else {
                        v = curve_y(v, out_trc[c], out_len);
                    }
                    p[3*(w*y+x)+c] = v;
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            free(out_trc[i]);
        }
    }

    cmsHTRANSFORM tfm = NULL;
    if (icc && icc->lut) {
        tfm = icc_lookup_fast_init(icc);
    }
    int xyz_clip_count = 0, alt_clip_count = 0;
#pragma omp parallel for
    for (int y = 0; y < h; y++) {
        int xyz_clip_count_ = 0, alt_clip_count_ = 0;
        for (int x = 0; x < w; x++) {
            v3 cam = {{ p[3*(w*y+x)+0], p[3*(w*y+x)+1], p[3*(w*y+x)+2] }};
            v3 cxyz, rgb;
            if (tfm == NULL) {
                cam.v[0] *= wb->v[0];
                cam.v[1] *= wb->v[1];
                cam.v[2] *= wb->v[2];
                cxyz = m3x3_mul_v3(*cam2xyz, cam);
                bool base_xyz_clip = false;
                int pre_clip = xyz_clip_count_ + alt_clip_count_;
                for (int i = 0; i < 3; i++) {
                    if (cxyz.v[i] < 0 || cxyz.v[i] > 1) {
                        base_xyz_clip = true;
                        break;
                    }
                }
                if (!skip_lut) {
                    bool xyz_clip = false, alt_clip = false;
                    if (lut) {
                        cxyz = chromalut_lookup_noclip(lut, cxyz);
                        for (int i = 0; i < 3; i++) {
                            if (cxyz.v[i] < 0 || cxyz.v[i] > 1) {
                                xyz_clip = true;
                            }
                            if (cxyz.v[i] < 0) cxyz.v[i] = 0;
                        }
                    } else if (dcp_lut && dcp_lut->hsm) {
                        cxyz = dcp_lookup(dcp_lut->hsm, dcp_lut->hsmdims, dcp_lut->srgb_gamma, cxyz, 1, NULL, 0, &xyz_clip, &alt_clip);
                    }
                    if (!skip_look_lut && dcp_lut_look && dcp_lut_look->hsm) {
                        bool xyz_clip1 = false, alt_clip1 = false;
                        cxyz = dcp_lookup(dcp_lut_look->hsm, dcp_lut_look->hsmdims, dcp_lut_look->srgb_gamma, cxyz, 1, NULL, 0, &xyz_clip1, &alt_clip1);
                        xyz_clip |= xyz_clip1;
                        alt_clip |= alt_clip1;
                    }
                    xyz_clip_count_ += (int)xyz_clip;
                    alt_clip_count_ += (int)alt_clip;
                } else {
                    xyz_clip_count_ += (int)base_xyz_clip;
                }
                rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, cxyz);
                if (lut || skip_lut) {
                    for (int i = 0; i < 3; i++) {
                        if (rgb.v[i] < 0 || rgb.v[i] > 1) {
                            alt_clip_count_++;
                            break;
                        }
                    }
                }
                if (tc != NULL) {
                    rgb = dcp_apply_tonecurve(rgb, tc, tc_len);
                }
                int post_clip = xyz_clip_count_ + alt_clip_count_;
                if (show_clip && pre_clip != post_clip) {
                    rgb.v[0] = 1;
                    rgb.v[1] = 0;
                    rgb.v[2] = 0;
                }
            } else {
                cxyz = icc_lookup_fast(tfm, icc->pcs_is_lab, cam);
                rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, cxyz);
            }
            rgb = gamut_acr_rgb_clip(rgb);
            p[3*(w*y+x)+0] = rgb.v[0];
            p[3*(w*y+x)+1] = rgb.v[1];
            p[3*(w*y+x)+2] = rgb.v[2];
        }

#pragma omp critical
        {
            xyz_clip_count += xyz_clip_count_;
            alt_clip_count += alt_clip_count_;
        }
    }
    *clip1 = (double)xyz_clip_count / (w*h);
    *clip2 = (double)alt_clip_count / (w*h);
    if (tfm != NULL) {
        cmsDeleteTransform(tfm);
    }

#pragma omp parallel for
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            for (int c = 0; c < 3; c++) {
                float v = p[3*(w*y+x)+c];
                img->p[3*(w*y+x)+c] = (uint16_t)roundf(pow(v, 1.0/img_gamma) * 65535); // made to match prophoto so we get correct display
            }
        }
    }
    free(p);

    if (output_image_filename != NULL) {
        tifio_rgbimg_save(output_image_filename, img, COLORSPACE_PROPHOTO, true);
    } else if (report_dir != NULL) {
        char fname[strlen(report_dir) + 100];
        if (custom_rgbimg != NULL) {
            sprintf(fname, "%s/test-image.tif", report_dir);
        } else {
            sprintf(fname, "%s/gradient.tif", report_dir);
        }
        tifio_rgbimg_save(fname, img, COLORSPACE_PROPHOTO, true);
    }
    free(img);
}

static struct patch_error_report *
test_patch_matching(const m3x3 *cam2xyz,
                    chromalut_t *lut,
                    struct dcp_lut *dcp_lut,
                    const struct profio_icc *icc,
                    const struct patch_set *ps,
                    const v3 ps_wp,
                    const v3 *wb)
{
    struct patch_error_report *pe = calloc(1, sizeof(*pe) + sizeof(pe->e[0]) * ps->patch_count);
    double sum[4] = {0,0,0,0};
    cmsHTRANSFORM tfm = NULL;
    if (lut) {
        pe->lookup_type = LOOKUP_LUT_NATIVE;
    } else if (dcp_lut && dcp_lut->hsm) {
        pe->lookup_type = LOOKUP_DCP_HSM;
    } else if (icc && icc->lut) {
        pe->lookup_type = LOOKUP_ICC_LUT;
        tfm = icc_lookup_fast_init(icc);
    } else {
        pe->lookup_type = LOOKUP_MATRIX;
    }
    for (int i = 0; i < ps->patch_count; i++) {
        v3 cxyz;
        if (pe->lookup_type != LOOKUP_ICC_LUT) {
            cxyz = m3x3_mul_v3(*cam2xyz, ps->patch[i].cam);
            if (lut) {
                cxyz = chromalut_lookup(lut, cxyz);
            } else if (dcp_lut && dcp_lut->hsm) {
                cxyz = dcp_lookup(dcp_lut->hsm, dcp_lut->hsmdims, dcp_lut->srgb_gamma, cxyz, 1, NULL, 0, NULL, NULL);
            }
        } else {
            v3 cam_wb = {{ ps->patch[i].cam.v[0] / wb->v[0], ps->patch[i].cam.v[1] / wb->v[1], ps->patch[i].cam.v[2] / wb->v[2] }};
            cxyz = icc_lookup_fast(tfm, icc->pcs_is_lab, cam_wb);
        }
        v3 xyz = ps->patch[i].xyz;
        v3 lch = lab2lch(xyz2lab(xyz, ps_wp));
        v3 clch = lab2lch(xyz2lab(cxyz, ps_wp));
        pe->e[i].idx = i;
        pe->e[i].cxyz = cxyz;
        pe->e[i].errors[0] = ciede2000(xyz, cxyz, ps_wp);
        pe->e[i].errors[1] = ciede2000_lch(xyz, cxyz, ps_wp, 1, 100000, 100000);
        if (lch.v[0] > clch.v[0]) pe->e[i].errors[1] = -pe->e[i].errors[1];
        pe->e[i].errors[2] = ciede2000_lch(xyz, cxyz, ps_wp, 100000, 1, 100000);
        if (lch.v[1] > clch.v[1]) pe->e[i].errors[2] = -pe->e[i].errors[2];
        pe->e[i].errors[3] = ciede2000_lch(xyz, cxyz, ps_wp, 100000, 100000, 1);
        if (lch.v[2] > clch.v[2]) pe->e[i].errors[3] = -pe->e[i].errors[3];
        for (int j = 0; j < 4; j++) {
            pe->e[i].errors[j] = round(1000 * pe->e[i].errors[j]) / 1000;
            if (pe->e[i].errors[j] == -0.0) {
                pe->e[i].errors[j] = +0.0;
            }
            if (fabs(pe->e[i].errors[j]) > pe->errors[j].max) {
                pe->errors[j].max = fabs(pe->e[i].errors[j]);
            }
            sum[j] += fabs(pe->e[i].errors[j]);
        }
    }
    if (tfm != NULL) {
        cmsDeleteTransform(tfm);
    }
    for (int j = 0; j < 4; j++) {
        pe->errors[j].avg = sum[j] / ps->patch_count;
    }
    pe->patch_count = ps->patch_count;
    pe->ps = ps;
    pe->ps_wp = ps_wp;
    struct patch_error *spe = malloc(sizeof(pe->e[0]) * ps->patch_count);
    memcpy(spe, pe->e, sizeof(pe->e[0]) * ps->patch_count);
    for (int j = 3; j >= 0; j--) {
        glob.pe_cmp_idx = j;
        qsort(spe, pe->patch_count, sizeof(spe[0]), pe_cmp);
        pe->errors[j].median = fabs(spe[pe->patch_count/2].errors[j]);
        pe->errors[j].p90 = fabs(spe[(9*pe->patch_count)/10].errors[j]);
        for (int k = 0; k < (int)(sizeof(pe->errors[j].worst)/sizeof(pe->errors[j].worst[0])) && k < ps->patch_count; k++) {
            pe->errors[j].worst[k] = spe[pe->patch_count - 1 - k].idx;
        }
        for (int k = 0; k < (int)(sizeof(pe->errors[j].best)/sizeof(pe->errors[j].best[0])) && k < ps->patch_count; k++) {
            pe->errors[j].best[k] = spe[k].idx;
        }
    }
    free(spe);
    return pe;
}

static void
print_profile_gamut_(FILE *stream,
                     const m3x3 *cam2xyz,
                     chromalut_t *lut,
                     struct dcp_lut *dcp_lut,
                     struct dcp_lut *dcp_lut_look,
                     const struct profio_icc *icc,
                     const v3 ps_wp,
                     const v3 *wb)
{
    UNUSED(ps_wp);
    cmsHTRANSFORM tfm = NULL;
    if (icc && icc->lut) {
        tfm = icc_lookup_fast_init(icc);
    }
    const int cs = 33;
    p2d *v = malloc(cs*cs*cs*sizeof(v[0]));
    int vcount = 0;
#pragma omp parallel for
    for (int r_ = 0; r_ < cs; r_++) {
        for (int g_ = 0; g_ < cs; g_++) {
            for (int b_ = 0; b_ < cs; b_++) {
                if (r_ == 0 && g_ == 0 && b_ == 0) continue;
                double r = (double)r_ / (cs-1);
                double g = (double)g_ / (cs-1);
                double b = (double)b_ / (cs-1);
                v3 cam = {{ r, g, b}};
                v3 cxyz;
                if (tfm == NULL) {
                    cxyz = m3x3_mul_v3(*cam2xyz, cam);
                    if (lut) {
                        cxyz = chromalut_lookup(lut, cxyz);
                    } else if (dcp_lut && dcp_lut->hsm) {
                        cxyz = dcp_lookup(dcp_lut->hsm, dcp_lut->hsmdims, dcp_lut->srgb_gamma, cxyz, 1, NULL, 0, NULL, NULL);
                    }
                    if (dcp_lut_look && dcp_lut_look->hsm) {
                        bool xyz_clip1 = false, alt_clip1 = false;
                        cxyz = dcp_lookup(dcp_lut_look->hsm, dcp_lut_look->hsmdims, dcp_lut_look->srgb_gamma, cxyz, 1, NULL, 0, &xyz_clip1, &alt_clip1);
                    }
                } else {
                    v3 cam_wb = {{ cam.v[0] / wb->v[0], cam.v[1] / wb->v[1], cam.v[2] / wb->v[2] }};
                    cxyz = icc_lookup_fast(tfm, icc->pcs_is_lab, cam_wb);
                }
                if (v3_isfinite(cxyz) && cxyz.v[0] >= 0 && cxyz.v[1] >= 0 && cxyz.v[2] >= 0) {
                    v3 ltu = xyz2lutspace(cxyz);
                    ltu.v[1] = ((int)(100000 * ltu.v[1])) / 100000.0;
                    ltu.v[2] = ((int)(100000 * ltu.v[2])) / 100000.0;
#pragma omp critical
                    {
                        v[vcount].x = ltu.v[1];
                        v[vcount].y = ltu.v[2];
                        vcount++;
                    }
                }
            }
        }
    }
    int hcount;
    p2d *h = convex_hull(v, vcount, &hcount);
    free(v);
    if (h != NULL) {
        for (int i = 0; i < hcount; i++) {
            fprintf(stream, "%f %f 0.0\n", h[i].x, h[i].y);
        }
        fprintf(stream, "%f %f 0.0\n", h[0].x, h[0].y);
        free(h);
    }
    if (tfm != NULL) {
        cmsDeleteTransform(tfm);
    }
}

static void
print_profile_gamut(const char report_dir[],
                    const m3x3 *cam2xyz,
                    chromalut_t *lut,
                    struct dcp_lut *dcp_lut,
                    struct dcp_lut *dcp_lut_look,
                    const struct profio_icc *icc,
                    const v3 ps_wp,
                    const v3 *wb)
{
    FILE *stream = open_log_file(report_dir, "gmt-prof.dat");
    if (stream == NULL) return;
    print_profile_gamut_(stream, cam2xyz, lut, dcp_lut, NULL, icc, ps_wp, wb);
    fclose(stream);

    if (dcp_lut_look && dcp_lut_look->hsm) {
        stream = open_log_file(report_dir, "gmt-prof-look.dat");
        if (stream == NULL) return;
        print_profile_gamut_(stream, cam2xyz, lut, dcp_lut, dcp_lut_look, icc, ps_wp, wb);
        fclose(stream);
    }
}

static void
print_matrices(const char report_dir[],
               const m3x3 *cm,
               const m3x3 *fm,
               const m3x3 *lm,
               const m3x3 *cm2,
               const m3x3 *fm2,
               const v3 ps_wp)
{
    const bool invert[5] = { true, false, false, true, false };
    const m3x3 *ms[5] = { cm, fm, lm, cm2, fm2 };
    const char *names[5] = { "gmt-cm.dat", "gmt-fm.dat", "gmt-lm.dat", "gmt-cm2.dat", "gmt-fm2.dat" };

    for (int i = 0; i < 5; i++) {
        if (ms[i] == NULL) continue;
        const m3x3 m = invert[i] ? m3x3_invert(*ms[i]) : *ms[i];
        FILE *stream = open_log_file(report_dir, names[i]);
        if (stream == NULL) return;
        v3 prim[3];
        prim[0] = m3x3_mul_v3(m, (v3){{1,0,0}});
        prim[1] = m3x3_mul_v3(m, (v3){{0,1,0}});
        prim[2] = m3x3_mul_v3(m, (v3){{0,0,1}});
        bool isbad = false;
        for (int i = 0; i < 3; i++) {
            if (prim[i].v[0] < 0 || prim[i].v[1] < 0 || prim[i].v[2] < 0) {
                isbad = true;
                break;
            }
        }
        if (isbad) {
            print_profile_gamut_(stream, &m, NULL, NULL, NULL, NULL, ps_wp, NULL);
        } else {
            prim[0] = xyz2lutspace(prim[0]);
            prim[1] = xyz2lutspace(prim[1]);
            prim[2] = xyz2lutspace(prim[2]);
            fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
            fprintf(stream, "%f %f 0.0\n", prim[1].v[1], prim[1].v[2]);
            fprintf(stream, "%f %f 0.0\n", prim[2].v[1], prim[2].v[2]);
            fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
        }
        fclose(stream);
    }
}

static void
print_common_spectra(const char report_dir[],
                     const struct observer *obs,
                     const spectrum_t *ill,
                     const struct dcam_ssf *ssf)
{
    if (obs != NULL) {
        FILE *stream1 = open_log_file(report_dir, "cmf-x.dat");
        FILE *stream2 = open_log_file(report_dir, "cmf-y.dat");
        FILE *stream3 = open_log_file(report_dir, "cmf-z.dat");
        if (stream1 == NULL) return;
        spec_print(obs->cmf[0], stream1);
        spec_print(obs->cmf[1], stream2);
        spec_print(obs->cmf[2], stream3);
        fclose(stream1);
        fclose(stream2);
        fclose(stream3);
    }
    if (ssf != NULL) {
        FILE *stream1 = open_log_file(report_dir, "ssf-r.dat");
        FILE *stream2 = open_log_file(report_dir, "ssf-g.dat");
        FILE *stream3 = open_log_file(report_dir, "ssf-b.dat");
        if (stream1 == NULL) return;
        spec_print(ssf->ssf.cmf[0], stream1);
        spec_print(ssf->ssf.cmf[1], stream2);
        spec_print(ssf->ssf.cmf[2], stream3);
        fclose(stream1);
        fclose(stream2);
        fclose(stream3);
    }
    if (ill != NULL) {
        FILE *stream1 = open_log_file(report_dir, "illuminant.dat");
        if (stream1 == NULL) return;
        spec_print(ill, stream1);
        fclose(stream1);
    }
    {
        FILE *stream1 = open_log_file(report_dir, "illuminant-d50.dat");
        if (stream1 == NULL) return;
        spec_print(spectraldb_illuminant(lsD50), stream1);
        fclose(stream1);
    }
}

static void
print_transfer_function(const char report_dir[],
                        double *trc[3],
                        int trc_count)
{
    if (trc_count > 0) {
        FILE *stream[3];
        stream[0] = open_log_file(report_dir, "tf-r.dat");
        stream[1] = open_log_file(report_dir, "tf-g.dat");
        stream[2] = open_log_file(report_dir, "tf-b.dat");
        if (stream[0] == NULL) return;
        for (int i = 0; i < 3; i++) {
            if (trc_count == 1) {
                for (int j = 0; j < 1024; j++) {
                    double x = (double)j / 1023;
                    fprintf(stream[i], "%f %f\n", x, pow(x, trc[i][0]));
                }
            } else {
                for (int j = 0; j < trc_count; j++) {
                    double x = (double)j / (trc_count-1);
                    fprintf(stream[i], "%f %f\n", x, trc[i][j]);
                }
            }
            fclose(stream[i]);
        }
    }
}

static void
print_tone_curve(const char report_dir[],
                 const double tc[],
                 int tc_len,
                 bool has_xy)
{
    if (tc != NULL) {
        FILE *stream[2];
        stream[0] = open_log_file(report_dir, "tc.dat");
        stream[1] = open_log_file(report_dir, "tc-srgb.dat");
        if (stream[0] == NULL) return;
        for (int j = 0; j < tc_len; j++) {
            if (has_xy) {
                fprintf(stream[0], "%f %f\n", tc[2*j+0], tc[2*j+1]);
                fprintf(stream[1], "%f %f\n", srgb_gamma_forward(tc[2*j+0]), srgb_gamma_forward(tc[2*j+1]));
            } else {
                double x = (double)j / (tc_len-1);
                fprintf(stream[0], "%f %f\n", x, tc[j]);
                fprintf(stream[1], "%f %f\n", srgb_gamma_forward(x), srgb_gamma_forward(tc[j]));
            }
        }
        fclose(stream[0]);
        fclose(stream[1]);
    }
}

static const char *
patch_name(char buf[PATCH_STRID_SIZE],
           const struct patch_set *ps,
           int idx)
{
    return target_patch_name(buf, ps, idx);
}

void
render_pm_error_image(const char report_dir[],
                      const char filename[],
                      const v3 ps_wp,
                      const struct patch_error_report *pe,
                      bool reverse_order)
{
    if (report_dir == NULL || filename == NULL) {
        return;
    }
    int margin = 5;
    int bs = 38;
    int w = 3 * margin + bs + 9 * 40;
    int h = margin + (bs+margin) * pe->patch_count + margin;
    struct rgbimg *img = malloc(sizeof(*img) + w * h * 3 * sizeof(img->p[0]));
    img->w = w;
    img->h = h;
    uint16_t bg = 32768;
    for (int y = 0; y < img->h; y++) {
        for (int x = 0; x < img->w; x++) {
            img->p[3*(img->w*y+x)+0] = bg;
            img->p[3*(img->w*y+x)+1] = bg;
            img->p[3*(img->w*y+x)+2] = bg;
        }
    }
    int offset = margin;
    for (int i = 0; i < pe->patch_count; i++) {
        const struct patch_error *e = &pe->e[reverse_order ? pe->patch_count - i - 1 : i];
        const struct cam2xyz_patch *p = &pe->ps->patch[e->idx];
        v3 ref = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, chromatic_adaptation_transform(CAT_CAT02, p->xyz, ps_wp, dcp_d50()));
        v3 cam = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, chromatic_adaptation_transform(CAT_CAT02, e->cxyz, ps_wp, dcp_d50()));
        for (int i = 0; i < 3; i++) {
            if (ref.v[i] < 0) ref.v[i] = 0;
            if (ref.v[i] > 1) ref.v[i] = 1;
            if (cam.v[i] < 0) cam.v[i] = 0;
            if (cam.v[i] > 1) cam.v[i] = 1;
        }
        for (int y = offset; y < offset + bs; y++) {
            for (int x = margin; x < margin+bs; x++) {
                v3 v;
                if (y - offset <= bs - (x - margin) - 1) {
                    v = ref;
                } else {
                    v = cam;
                }
                img->p[3*(img->w*y+x)+0] = pow(v.v[0], 1/1.8) * 65535;
                img->p[3*(img->w*y+x)+1] = pow(v.v[1], 1/1.8) * 65535;
                img->p[3*(img->w*y+x)+2] = pow(v.v[2], 1/1.8) * 65535;
            }
        }
        char msg[256], idbuf[PATCH_STRID_SIZE];
        sprintf(msg,
                "%5s DE %.2f DE LCh %+.2f %+.2f %+.2f",
                patch_name(idbuf, pe->ps, e->idx),
                e->errors[0],
                e->errors[1],
                e->errors[2],
                e->errors[3]);
        tifio_type(img, msg, bg, 0, bs + 2 * margin, offset + bs / 2 - 7);
        offset += bs + margin;
    }
    char fname[strlen(report_dir) + 32 + strlen(filename)];
    sprintf(fname, "%s/%s", report_dir, filename);
    tifio_rgbimg_save(fname, img, COLORSPACE_PROPHOTO, true);
    free(img);
}

static void
print_patch_matching_errors(const struct patch_error_report *pe,
                            const v3 ps_wp,
                            bool print_worst_best,
                            const char report_dir[],
                            const char report_filename[],
                            FILE *stream)
{
    const char *name = "Matrix";
    const char *pad  = "      ";
    if (pe->lookup_type == LOOKUP_LUT_NATIVE) {
        name = "Native LUT";
        pad  = "          ";
    }
    if (pe->lookup_type == LOOKUP_DCP_HSM) {
        name = "DCP HSM LUT";
        pad  = "           ";
    }
    if (pe->lookup_type == LOOKUP_ICC_LUT) {
        name = "ICC LUT";
        pad  = "       ";
    }
    fprintf(stream, "%s patch match average DE %.2f, DE LCh %.2f %.2f %.2f\n",
            name, pe->errors[0].avg, pe->errors[1].avg, pe->errors[2].avg, pe->errors[3].avg);
    fprintf(stream, "%s              median DE %.2f, DE LCh %.2f %.2f %.2f\n",
            pad, pe->errors[0].median, pe->errors[1].median, pe->errors[2].median, pe->errors[3].median);
    fprintf(stream, "%s                 p90 DE %.2f, DE LCh %.2f %.2f %.2f\n",
            pad, pe->errors[0].p90, pe->errors[1].p90, pe->errors[2].p90, pe->errors[3].p90);
    fprintf(stream, "%s                 max DE %.2f, DE LCh %.2f %.2f %.2f\n",
            pad, pe->errors[0].max, pe->errors[1].max, pe->errors[2].max, pe->errors[3].max);

    const char error_fmt[] = "%5s RGB %.3f %.3f %.3f XYZref %.3f %.3f %.3f XYZcam %.3f %.3f %.3f sRGB #%06X #%06X DE %.2f DE LCh %+.2f %+.2f %+.2f (%s)\n";
    char cnbuf[128], idbuf[PATCH_STRID_SIZE];
    if (report_filename != NULL) {
        char fname[strlen(report_filename)+100];
        sprintf(fname, "%s.txt", report_filename);
        FILE *stream1 = open_log_file(report_dir, fname);
        if (stream1 != NULL) {
            for (int i = 0; i < pe->patch_count; i++) {
                const struct patch_error *e = &pe->e[i];
                const struct cam2xyz_patch *p = &pe->ps->patch[e->idx];
                fprintf(stream1, error_fmt, patch_name(idbuf, pe->ps, e->idx),
                        p->cam.v[0], p->cam.v[1], p->cam.v[2], p->xyz.v[0], p->xyz.v[1], p->xyz.v[2],
                        e->cxyz.v[0], e->cxyz.v[1], e->cxyz.v[2],
                        xyz2srgb_u24(p->xyz, ps_wp),
                        xyz2srgb_u24(e->cxyz, ps_wp),
                        e->errors[0],
                        e->errors[1],
                        e->errors[2],
                        e->errors[3],
                        xyz2name(p->xyz, pe->ps_wp, cnbuf));
            }
            fclose(stream1);
        }
        sprintf(fname, "%s.tif", report_filename);
        render_pm_error_image(report_dir, fname, ps_wp, pe, false);

        struct patch_error_report *pe1 = malloc(sizeof(*pe) + sizeof(pe->e[0]) * pe->ps->patch_count);
        *pe1 = *pe;
        struct patch_error *spe = pe1->e;
        memcpy(spe, pe->e, sizeof(pe->e[0]) * pe->ps->patch_count);
        for (int j = 3; j >= 0; j--) {
            glob.pe_cmp_idx = j;
            qsort(spe, pe->patch_count, sizeof(spe[0]), pe_cmp);
            const char *dim[4] = { "-DE", "-DE-L", "-DE-C", "-DE-h" };
            sprintf(fname, "%s%s.tif", report_filename, dim[j]);
            render_pm_error_image(report_dir, fname, ps_wp, pe1, true);
        }
        free(pe1);
    }

    if (!print_worst_best) {
        return;
    }
    int worst_count = sizeof(pe->errors[0].worst)/sizeof(pe->errors[0].worst[0]);
    if (worst_count > pe->patch_count) worst_count = pe->patch_count;
    for (int j = 0; j < 4; j++) {
        const char *error_type[] = {
            "Overall DE",
            "Lightness DE",
            "Chroma DE",
            "Hue DE",
        };
        fprintf(stream, "%d worst patches for %s:\n", worst_count, error_type[j]);
        for (int k = 0; k < worst_count; k++) {
            const struct patch_error *e = &pe->e[pe->errors[j].worst[k]];
            const struct cam2xyz_patch *p = &pe->ps->patch[e->idx];
            fprintf(stream, error_fmt, patch_name(idbuf, pe->ps, e->idx),
                    p->cam.v[0], p->cam.v[1], p->cam.v[2], p->xyz.v[0], p->xyz.v[1], p->xyz.v[2],
                    e->cxyz.v[0], e->cxyz.v[1], e->cxyz.v[2],
                    xyz2srgb_u24(p->xyz, ps_wp),
                    xyz2srgb_u24(e->cxyz, ps_wp),
                    e->errors[0],
                    e->errors[1],
                    e->errors[2],
                    e->errors[3],
                    xyz2name(p->xyz, pe->ps_wp, cnbuf));
        }
    }
    int best_count = sizeof(pe->errors[0].best)/sizeof(pe->errors[0].best[0]);
    if (best_count > pe->patch_count) best_count = pe->patch_count;
    fprintf(stream, "%d best patches for Overall DE:\n", best_count);
    for (int k = 0; k < best_count; k++) {
        const struct patch_error *e = &pe->e[pe->errors[0].best[k]];
        const struct cam2xyz_patch *p = &pe->ps->patch[e->idx];
        fprintf(stream, error_fmt, patch_name(idbuf, pe->ps, e->idx),
                p->cam.v[0], p->cam.v[1], p->cam.v[2], p->xyz.v[0], p->xyz.v[1], p->xyz.v[2],
                e->cxyz.v[0], e->cxyz.v[1], e->cxyz.v[2],
                xyz2srgb_u24(p->xyz, ps_wp),
                xyz2srgb_u24(e->cxyz, ps_wp),
                e->errors[0],
                e->errors[1],
                e->errors[2],
                e->errors[3],
                xyz2name(p->xyz, pe->ps_wp, cnbuf));
    }
}

static void
print_pm_errors(const m3x3 cam2xyz,
                chromalut_t *lut,
                const struct patch_set *ps,
                const v3 ps_wp,
                bool print_worst_best,
                bool print_per_class,
                const char report_dir[],
                const char report_filename[])
{
    if (print_per_class && ps->class_names[1][0] != '\0') {
        for (int cl = 0; cl < TARGET_MAX_CLASS_COUNT && ps->class_names[cl][0] != '\0'; cl++) {
            struct patch_set *ps1 = malloc(sizeof(*ps1) + ps->patch_count * sizeof(ps1->patch[0]));
            *ps1 = *ps;
            int k = 0;
            for (int i = 0; i < ps->patch_count; i++) {
                if (ps->patch[i].class == cl) {
                    ps1->patch[k] = ps->patch[i];
                    ps1->patch[k].class = 0;
                    ps1->patch[k].spectrum = NULL;
                    k++;
                }
            }
            ps1->patch_count = k;
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "\n"
                "--------------------------------------------------------------------------------\n");
            elog("Patch matching for class \"%s\" (%d patches):\n\n", ps->class_names[cl], k);
            memset(ps1->class_names, 0, sizeof(ps1->class_names));
            strcpy(ps1->class_names[0], ps->class_names[cl]);
            print_pm_errors(cam2xyz, lut, ps1, ps_wp, print_worst_best, false, NULL, NULL);
            free(ps1);
        }
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "\n"
            "--------------------------------------------------------------------------------\n");
        elog("Patch matching for all classes (%d patches):\n\n", ps->patch_count);
    }
    struct patch_error_report *pe = test_patch_matching(&cam2xyz, lut, NULL, NULL, ps, ps_wp, NULL);
    print_patch_matching_errors(pe, ps_wp, print_worst_best, report_dir, report_filename, stderr);
    free(pe);
}

static double
test_matrix_lsq(const m3x3 *cam2xyz,
                const v3 *wp,
                const struct patch_set *ps,
                void *arg)
{
    if (arg != NULL) {
        const double *smallest = (const double *)arg;
        if (cam2xyz->v[0][2] < smallest[0] || cam2xyz->v[0][1] < smallest[0] || cam2xyz->v[0][0] < smallest[0]) return 1e9;
        if (cam2xyz->v[1][2] < smallest[1] || cam2xyz->v[1][1] < smallest[1] || cam2xyz->v[1][0] < smallest[1]) return 1e9;
        if (cam2xyz->v[2][2] < smallest[2] || cam2xyz->v[2][1] < smallest[2] || cam2xyz->v[2][0] < smallest[2]) return 1e9;
    }
    double sum = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].mtx_w == 0) {
            continue;
        }
        v3 cxyz = m3x3_mul_v3(*cam2xyz, ps->patch[i].cam);
        double de = ciede2000(ps->patch[i].xyz, cxyz, *wp);
        double wde = ps->patch[i].mtx_w * de;
        sum += wde * wde;
    }
    return sqrt(sum); // optimizer works better with slower approach to zero
}

struct test_matrix_lsq_lut_arg {
    double sum1, sum2;
    double sum1_mul, sum2_mul;
    double poly_x[3];
    double poly_y[3];
    double smallest;
};

static double
test_matrix_lsq_lut(const m3x3 *cam2xyz,
                    const v3 *wp,
                    const struct patch_set *ps,
                    void *arg)
{
    struct test_matrix_lsq_lut_arg *a = (struct test_matrix_lsq_lut_arg *)arg;

    if (cam2xyz->v[0][2] < a->smallest || cam2xyz->v[0][1] < a->smallest || cam2xyz->v[0][0] < a->smallest) return 1e9;
    if (cam2xyz->v[1][2] < a->smallest || cam2xyz->v[1][1] < a->smallest || cam2xyz->v[1][0] < a->smallest) return 1e9;
    if (cam2xyz->v[2][2] < a->smallest || cam2xyz->v[2][1] < a->smallest || cam2xyz->v[2][0] < a->smallest) return 1e9;

    double sum1;
    {
        double sum = 0;
        v3 uvY[3];
        uvY[0] = xyz2uvY(m3x3_mul_v3(*cam2xyz, (v3){{1,0,0}}));
        uvY[1] = xyz2uvY(m3x3_mul_v3(*cam2xyz, (v3){{0,1,0}}));
        uvY[2] = xyz2uvY(m3x3_mul_v3(*cam2xyz, (v3){{0,0,1}}));

        for (int i = 0; i < 3; i++) {
            int j = i+1;
            if (j == 3) j = 0;
            double de = 100*point_to_line_distance(a->poly_x[i], a->poly_y[i], a->poly_x[j], a->poly_y[j], uvY[i].v[0], uvY[i].v[1]);
            if (!point_inside_polygon(a->poly_x, a->poly_y, 3, uvY[i].v[0], uvY[i].v[1])) {
                de *= 3;
            }
            sum += de * de;
        }
        sum /= 3;
        sum1 = sqrt(sum);
        sum1 *= a->sum1_mul;
    }

    double sum2;
    {
        double sum = 0;
        int count = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            if (ps->patch[i].mtx_w == 0) {
                continue;
            }
            v3 cxyz = m3x3_mul_v3(*cam2xyz, ps->patch[i].cam);
            double de = ciede2000(ps->patch[i].xyz, cxyz, *wp);
            sum += de * de;
            count++;
        }
        sum /= count;
        sum2 = sqrt(sum);
        sum2 *= a->sum2_mul;
    }
    a->sum1 = sum1;
    a->sum2 = sum2;
    return sum1 + sum2;
}

static double
test_matrix_lsq_refinement(const m3x3 *cam2xyz,
                           const v3 *wp,
                           const struct patch_set *ps,
                           void *arg)
{
    int w_idx = (int)((size_t)arg);

    // white must keep reasonably close
    double sum = 0;
    if (w_idx >= 0) {
        v3 cxyz = m3x3_mul_v3(*cam2xyz, ps->patch[w_idx].cam);
        double de = ciede2000(ps->patch[w_idx].xyz, cxyz, *wp);
        if (de > 0.5) {
            return 1e9;
        }
        sum += de * de;
    }

    for (int i = 0; i < ps->patch_count; i++) {
        if (!ps->patch[i].use_mtx_de_w || i == w_idx) continue;
        v3 cxyz = m3x3_mul_v3(*cam2xyz, ps->patch[i].cam);
        // if zero DE
        if (ps->patch[i].mtx_de_w[0] == 0 && ps->patch[i].mtx_de_w[3] == 0 &&
            ps->patch[i].mtx_de_w[1] == 0 && ps->patch[i].mtx_de_w[4] == 0 &&
            ps->patch[i].mtx_de_w[2] == 0 && ps->patch[i].mtx_de_w[5] == 0)
        {
            double de = ciede2000(ps->patch[i].xyz, cxyz, *wp);
            sum += de * de;
            continue;
        }
        // if same range on all dimensions, check total DE instead (will lead to smaller error due to shorter euclidean distance)
        if (ps->patch[i].mtx_de_w[0] == -ps->patch[i].mtx_de_w[1] &&
            ps->patch[i].mtx_de_w[2] == -ps->patch[i].mtx_de_w[3] &&
            ps->patch[i].mtx_de_w[4] == -ps->patch[i].mtx_de_w[5] &&
            ps->patch[i].mtx_de_w[0] == ps->patch[i].mtx_de_w[2] &&
            ps->patch[i].mtx_de_w[0] == ps->patch[i].mtx_de_w[4])
        {
            double de = ciede2000(ps->patch[i].xyz, cxyz, *wp);
            de -= ps->patch[i].mtx_de_w[1];
            if (de < 0) de = 0;
            sum += de * de;
            continue;
        }
        v3 clch = lab2lch(xyz2lab(cxyz, *wp));
        v3 lch = lab2lch(xyz2lab(ps->patch[i].xyz, *wp));
        double de[3];
        de[0] = ciede2000_lch(ps->patch[i].xyz, cxyz, *wp, 1, 1000000, 1000000);
        de[1] = ciede2000_lch(ps->patch[i].xyz, cxyz, *wp, 1000000, 1, 1000000);
        de[2] = ciede2000_lch(ps->patch[i].xyz, cxyz, *wp, 1000000, 1000000, 1);
        double wde = 0;
        for (int j = 0; j < 3; j++) {
            if (clch.v[j] < lch.v[j]) de[j] = -de[j];
            double lo = ps->patch[i].mtx_de_w[2*j+0];
            double hi = ps->patch[i].mtx_de_w[2*j+1];
            double err;
            if (de[j] < lo) {
                err = de[j] - lo;
            } else if (de[j] > hi) {
                err = de[j] - hi;
            } else {
                err = 0;
            }
            if (fabs(err) < 0.01) {
                err = 0;
            }
            wde += err * err;
        }
        wde = sqrt(wde);
        sum += wde * wde;
    }
    return sqrt(sum);
}

static m3x3
find_cam2xyz_matrix(const char illname[],
                    enum exif_lightsource calib_ls,
                    const v3 calib_wp,
                    const int w_idx,
                    const struct patch_set *ps,
                    bool use_slow_matrix_optimizer,
                    const double *smallest_per_row,
                    bool lut_matrix,
                    bool apply_refinement,
                    bool exit_on_refinement_error)
{
    elog("Finding a camera raw RGB to CIE XYZ matrix for %s illuminant %s...\n", illname, exif_lightsource_ntoa(calib_ls));
    if (smallest_per_row != NULL && (smallest_per_row[0] > -10 || smallest_per_row[1] > -10 || smallest_per_row[2] > -10)) {
        if (smallest_per_row[0] <= -10 && smallest_per_row[2] <= -10) {
            elog("  Y row limit set to %g.\n", smallest_per_row[1]);
        } else {
            elog("  X,Y,Z row limits set to %g,%g,%g.\n", smallest_per_row[0], smallest_per_row[1], smallest_per_row[2]);
        }
    }

    // preserve white point
    double de = ciede2000(v3_norm(ps->patch[w_idx].xyz), v3_norm(calib_wp), calib_wp);
    debug("whitest patch difference from %s illuminant whitepoint %f\n", illname, de);
    char w_name[sizeof(ps->patch[w_idx].str_id)+32];
    if (ps->patch[w_idx].str_id[0] == '\0') {
        sprintf(w_name, "index %d", w_idx+1);
    } else {
        strcpy(w_name, ps->patch[w_idx].str_id);
    }
    if (de > 2) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: whitest (most neutral) patch in target (%s) differs DE %.2f\n"
            "  from %s illuminant, matrix precision may suffer.\n", w_name, de, illname);
    } else if (de > 1e-5) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Whitest patch in target (%s) differs DE %.2f from %s illuminant,\n"
            "  close enough to calculate whitepoint preservation.\n", w_name, de, illname);
    }
    struct patch_set *wb_ps = calloc(1, sizeof(*wb_ps) + sizeof(wb_ps->patch[0]));
    wb_ps->patch_count = 1;
    wb_ps->patch[0].mtx_w = 1.0;
    wb_ps->patch[0].xyz = ps->patch[w_idx].xyz;
    wb_ps->patch[0].cam = ps->patch[w_idx].cam;

    m3x3 cam2xyz;
    if (lut_matrix) {
        if (use_slow_matrix_optimizer) {
            elog("LUT matrix only needs to be approximate, slow optimizer not used.\n");
        }
        double smallest[3] = { 0.0002, 0.0002, 0.0002 };
        cam2xyz = matopt_find_matrix(calib_wp, ps, smallest, test_matrix_lsq, smallest, NULL, NULL);

        struct test_matrix_lsq_lut_arg arg;
        memset(&arg, 0, sizeof(arg));

        m3x3 rgb2xyz = m3x3_invert(dcp_xyz_D50_to_prophoto_rgb);
        v3 uvY[3];
        uvY[0] = xyz2uvY(m3x3_mul_v3(rgb2xyz, (v3){{1,0,0}}));
        uvY[1] = xyz2uvY(m3x3_mul_v3(rgb2xyz, (v3){{0,1,0}}));
        uvY[2] = xyz2uvY(m3x3_mul_v3(rgb2xyz, (v3){{0,0,1}}));
        for (int i = 0; i < 3; i++) {
            arg.poly_x[i] = uvY[i].v[0];
            arg.poly_y[i] = uvY[i].v[1];
        }
        arg.smallest = smallest[1];
        arg.sum1_mul = 1;
        arg.sum2_mul = 1;
        test_matrix_lsq_lut(&cam2xyz, &calib_wp, ps, &arg);
        arg.sum1_mul = arg.sum2 / (arg.sum1 + arg.sum2);
        arg.sum2_mul = arg.sum1 / (arg.sum1 + arg.sum2);

        matopt_refine_matrix(&cam2xyz, calib_wp, ps, test_matrix_lsq_lut, &arg);
        matopt_refine_matrix(&cam2xyz, calib_wp, wb_ps, test_matrix_lsq_lut, &arg);
    } else {
        if (use_slow_matrix_optimizer) {
            elog("Warning: using slow optimizer, this may take a while...\n");
            cam2xyz = matopt_find_matrix_brute_force(calib_wp, ps, test_matrix_lsq, (void *)smallest_per_row);
        } else {
            cam2xyz = matopt_find_matrix(calib_wp, ps, smallest_per_row, test_matrix_lsq, (void *)smallest_per_row, NULL, NULL);
        }
        matopt_refine_matrix(&cam2xyz, calib_wp, wb_ps, test_matrix_lsq, (void *)smallest_per_row);

        int has_mtx_de_w_count = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            if (ps->patch[i].use_mtx_de_w) has_mtx_de_w_count++;
        }
        if (apply_refinement && has_mtx_de_w_count > 0) {
            elog("Refining %d patch%s according to specified error ranges...\n", has_mtx_de_w_count, has_mtx_de_w_count > 1 ? "es" : "");
            if (test_matrix_lsq_refinement(&cam2xyz, &calib_wp, ps, (void *)((size_t)w_idx)) > 0.00001) {
                matopt_refine_matrix(&cam2xyz, calib_wp, ps, test_matrix_lsq_refinement, (void *)((size_t)w_idx));
                matopt_refine_matrix(&cam2xyz, calib_wp, wb_ps, test_matrix_lsq, (void *)smallest_per_row);
            }

            // test if desired refinements are met
            char idbuf[PATCH_STRID_SIZE];
            const char *dimname[3] = { "lightness", "chroma", "hue" };
            bool acceptable_de = true;
            for (int i = 0; i < ps->patch_count; i++) {
                if (!ps->patch[i].use_mtx_de_w) continue;
                v3 cxyz = m3x3_mul_v3(cam2xyz, ps->patch[i].cam);
                v3 clch = lab2lch(xyz2lab(cxyz, calib_wp));
                v3 lch = lab2lch(xyz2lab(ps->patch[i].xyz, calib_wp));
                double de[3];
                de[0] = ciede2000_lch(ps->patch[i].xyz, cxyz, calib_wp, 1, 1000000, 1000000);
                de[1] = ciede2000_lch(ps->patch[i].xyz, cxyz, calib_wp, 1000000, 1, 1000000);
                de[2] = ciede2000_lch(ps->patch[i].xyz, cxyz, calib_wp, 1000000, 1000000, 1);
                for (int j = 0; j < 3; j++) {
                    if (clch.v[j] < lch.v[j]) de[j] = -de[j];
                    if (de[j] < ps->patch[i].mtx_de_w[2*j+0] - 0.1) {
                        elog("%s: %s refinement DE lower limit exceeded by %s (%.2f < %.2f)\n", exit_on_refinement_error ? "Error": "Warning", dimname[j], patch_name(idbuf, ps, i), de[j], ps->patch[i].mtx_de_w[2*j+0]);
                        acceptable_de = false;
                        break;
                    }
                    if (de[j] > ps->patch[i].mtx_de_w[2*j+1] + 0.1) {
                        elog("%s: %s refinement DE lower limit exceeded by %s (%.2f > %.2f)\n", exit_on_refinement_error ? "Error": "Warning", dimname[j], patch_name(idbuf, ps, i), de[j], ps->patch[i].mtx_de_w[2*j+1]);
                        acceptable_de = false;
                        break;
                    }
                }
            }
            if (exit_on_refinement_error && !acceptable_de) {
                elog("Error: there is no linear matrix solution to meet the DE requirements set.\n"
                     "  Relax the requirements and retry.\n");
                exit(EXIT_FAILURE);
            }
            if (acceptable_de) {
                elog("Refinements successful.\n");
            }
        }
    }

    v3 calcxyz = m3x3_mul_v3(cam2xyz, wb_ps->patch[0].cam);
    de = ciede2000(wb_ps->patch[0].xyz, calcxyz, calib_wp);
    if (de > 1e-5) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: did not succeeed to match whitest patch, DE %.2f.\n"
            "  Precision may suffer.\n", de);
    }
    debug("DE %f\n", de);
    debug("%+.4f %+.4f %+.4f\n", cam2xyz.v[0][0], cam2xyz.v[0][1], cam2xyz.v[0][2]);
    debug("%+.4f %+.4f %+.4f\n", cam2xyz.v[1][0], cam2xyz.v[1][1], cam2xyz.v[1][2]);
    debug("%+.4f %+.4f %+.4f\n\n", cam2xyz.v[2][0], cam2xyz.v[2][1], cam2xyz.v[2][2]);
    free(wb_ps);
    return cam2xyz;
}

static double
norm_l_de(v3 xyz1, v3 xyz2, v3 wp)
{
    v3 lab1 = xyz2lab(xyz1, wp);
    v3 lab2 = xyz2lab(xyz2, wp);
    lab2.v[0] = lab1.v[0];
    xyz1 = lab2xyz(lab1, wp);
    xyz2 = lab2xyz(lab2, wp);
    return ciede2000(xyz1, xyz2, wp);
}

static int
double_cmp(const void *a_, const void *b_)
{
    double a = *(const double *)a_;
    double b = *(const double *)b_;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

static v3
de_vector(v3 ltu, v3 cltu, v3 wp, bool normalize)
{
    v3 xyz = lutspace2xyz(ltu);
    v3 cxyz = lutspace2xyz(cltu);
    v3 v;
    double len = 0;
    for (int i = 0; i < 3; i++) {
        v.v[i] = cltu.v[i]-ltu.v[i];
        len += v.v[i] * v.v[i];
    }
    len = sqrt(len);
    if (normalize) {
        double mx = v3_max(xyz);
        if (v3_max(cxyz) > mx) mx = v3_max(cxyz);
        cxyz = k_mul_v3(1/mx, cxyz);
        xyz = k_mul_v3(1/mx, xyz);
    }
    double de = ciede2000(xyz, cxyz, wp) * 0.01;
    for (int i = 0; i < 3; i++) {
        v.v[i] *= de / len;
    }
    return v;
}

static void
print_hsm(const char report_dir[],
          const char prefix[],
          const v3 *hsm,
          int hcount,
          int scount,
          int vcount,
          bool srgb_gamma)
{
    if (report_dir == NULL) {
        return;
    }
    char name[64];
    char name1[64];
    char name2[64];
    sprintf(name, "%s-ref.dat", prefix);
    sprintf(name1, "%s-lut.dat", prefix);
    sprintf(name2, "%s-lutv.dat", prefix);
    FILE *stream = open_log_file(report_dir, name);
    if (stream == NULL) {
        return;
    }
    FILE *stream1 = open_log_file(report_dir, name1);
    FILE *stream2 = open_log_file(report_dir, name2);
    m3x3 rgb2xyz = m3x3_invert(dcp_xyz_D50_to_prophoto_rgb);
    if (vcount == 1) {
        for (int h = 0; h < hcount; h++) {
            for (int s = 0; s < scount; s++) {
                const v3 *entry = &hsm[scount*h+s];
                v3 hsv = {{ h / (double)hcount * 6.0, s / (double)(scount-1), 0.5 }};
                v3 rgb = dcp_hsv2rgb(hsv);
                v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
                v3 ltu = xyz2lutspace(xyz);
                fprintf(stream, "%f %f 0.0\n", ltu.v[1], ltu.v[2]);
                v3 hsv1 = hsv;
                hsv1.v[0] += (entry->v[0] / 360.0 * 6.0);
                if (hsv1.v[0] < 0) hsv1.v[0] += 6.0;
                hsv1.v[1] *= entry->v[1];
                hsv1.v[2] *= entry->v[2];
                rgb = dcp_hsv2rgb(hsv1);
                xyz = m3x3_mul_v3(rgb2xyz, rgb);
                v3 ltu1 = xyz2lutspace(xyz);
                double ldiff = (ltu1.v[0]/ltu.v[0]-1)/10.0;
                fprintf(stream1, "%f %f %f\n", ltu1.v[1], ltu1.v[2], ldiff);
                fprintf(stream2, "%f %f %f %f %f %f\n", ltu.v[1], ltu.v[2], 0.0, ltu1.v[1] - ltu.v[1], ltu1.v[2] - ltu.v[2], ldiff);
            }
        }
    } else {
        for (int v = 0; v < vcount; v++) {
            char fname[64];
            sprintf(fname, "%s-lut%02d.dat", prefix, v);
            FILE *stream3 = open_log_file(report_dir, fname);
            for (int h = 0; h < hcount; h++) {
                for (int s = 0; s < scount; s++) {
                    const v3 *entry = &hsm[(hcount*scount)*v + scount*h + s];
                    v3 hsv = {{ h / (double)hcount * 6.0, s / (double)(scount-1), v / (double)(vcount-1) }};
                    if (srgb_gamma) {
                        hsv.v[2] = srgb_gamma_inverse(hsv.v[2]);
                    }
                    v3 rgb = dcp_hsv2rgb(hsv);
                    v3 xyz = m3x3_mul_v3(rgb2xyz, rgb);
                    v3 ltu = xyz2lutspace(xyz);
                    fprintf(stream, "%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
                    fprintf(stream3, "%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
                    v3 hsv1 = hsv;
                    hsv1.v[0] += (entry->v[0] / 360.0 * 6.0);
                    if (hsv1.v[0] < 0) hsv1.v[0] += 6.0;
                    hsv1.v[1] *= entry->v[1];
                    if (srgb_gamma) {
                        if (hsv.v[2] != 0) {
                            hsv1.v[2] *= srgb_gamma_inverse(entry->v[2] * srgb_gamma_forward(hsv.v[2])) / hsv.v[2];
                        }
                    } else {
                        hsv1.v[2] *= entry->v[2];
                    }
                    rgb = dcp_hsv2rgb(hsv1);
                    xyz = m3x3_mul_v3(rgb2xyz, rgb);
                    v3 ltu1 = xyz2lutspace(xyz);
                    fprintf(stream1, "%f %f %f\n", ltu1.v[1], ltu1.v[2], ltu1.v[0]);
                    fprintf(stream2, "%f %f %f %f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0], ltu1.v[1] - ltu.v[1], ltu1.v[2] - ltu.v[2], ltu1.v[0] - ltu.v[0]);
                }
            }
            fclose(stream3);
        }
    }
    if (stream) {
        fclose(stream);
        fclose(stream1);
        fclose(stream2);
    }
}

static void
print_lookup_reports(const char report_dir[],
                     chromalut_t *lut,
                     struct dcp_lut *dcp_lut,
                     const struct profio_icc *icc,
                     const struct dcam_ssf *dcam_ssf,
                     const m3x3 *cam2xyz,
                     const struct patch_set *ps_,
                     const v3 *ps_wp,
                     const int ps_white_idx,
                     const v3 *wb)
{
    if (report_dir == NULL) {
        return;
    }
    if (lut != NULL) {
        FILE *stream1 = open_log_file(report_dir, "nve-lut.dat");
        FILE *stream2 = open_log_file(report_dir, "nve-ref.dat");
        FILE *stream3 = open_log_file(report_dir, "nve-lutv.dat");
        if (stream1 == NULL) {
            return;
        }
        for (double u = 0.0; u <= 0.6; u += 0.01) {
            for (double v = 0.0; v <= 0.6; v += 0.01) {
                v3 diff = chromalut_lookup_diff(lut, u, v);
                double Ldiff = diff.v[0];
                double udiff = diff.v[1];
                double vdiff = diff.v[2];
                fprintf(stream1, "%f %f %f\n", u+udiff, v+vdiff, Ldiff);
                fprintf(stream2, "%f %f %f\n", u, v, 0.0);
                fprintf(stream3, "%f %f %f %f %f %f\n", u, v, 0.0, udiff, vdiff, Ldiff);
            }
            fprintf(stream1, "\n");
            fprintf(stream2, "\n");
            fprintf(stream3, "\n");
        }
        fclose(stream1);
        fclose(stream2);
        fclose(stream3);

        stream1 = open_log_file(report_dir, "nve-lutd.dat");
        for (double u = 0.0; u <= 0.6; u += 0.01 / 3.0) {
            for (double v = 0.0; v <= 0.6; v += 0.01 / 3.0) {
                v3 diff = chromalut_lookup_diff(lut, u, v);
                double Ldiff = diff.v[0];
                double udiff = diff.v[1];
                double vdiff = diff.v[2];
                fprintf(stream1, "%f %f %f\n", u+udiff, v+vdiff, Ldiff);
            }
            fprintf(stream1, "\n");
        }
        fclose(stream1);
    }
    if (dcp_lut && dcp_lut->hsm) {
        print_hsm(report_dir, "hsm", dcp_lut->hsm, dcp_lut->hsmdims[0], dcp_lut->hsmdims[1], dcp_lut->hsmdims[2], dcp_lut->srgb_gamma);
    }
    if (icc != NULL && icc->lut != NULL) {
        const int slice_count = 20;
        FILE *slices[slice_count];
        FILE *stream1 = open_log_file(report_dir, "icc-lut.dat");
        if (stream1 == NULL) {
            return;
        }
        double slice_lim[slice_count+1];
        for (int i = 0; i < slice_count; i++) {
            char name[64];
            sprintf(name, "icc-lut%02d.dat", i);
            slices[i] = open_log_file(report_dir, name);
            slice_lim[i] = (double)i / slice_count;
        }
        slice_lim[slice_count] = 100000000;
        const int cs = icc->lut->clut_side;
        for (int r_ = 0; r_ < cs; r_++) {
            for (int g_ = 0; g_ < cs; g_++) {
                for (int b_ = 0; b_ < cs; b_++) {
                    int oft = 3*(cs*cs*r_+cs*g_+b_);
                    v3 v = {{ icc->lut->clut[oft+0], icc->lut->clut[oft+1], icc->lut->clut[oft+2] }};
                    if (icc->pcs_is_lab) {
                        v = lab2xyz(v, icc_d50());
                    }
                    v3 ltu1 = xyz2lutspace(v);
                    fprintf(stream1, "%f %f %f\n", ltu1.v[1], ltu1.v[2], ltu1.v[0]);
                    for (int i = 0; i < slice_count; i++) {
                        if (ltu1.v[0] >= slice_lim[i] && ltu1.v[0] < slice_lim[i+1]) {
                            fprintf(slices[i], "%f %f %f\n", ltu1.v[1], ltu1.v[2], ltu1.v[0]);
                            break;
                        }
                    }
                }
            }
        }
        fclose(stream1);
        for (int i = 0; i < slice_count; i++) {
            fclose(slices[i]);
        }
    }

    if (ps_ != NULL) {
        FILE *stream1 = NULL, *stream2 = NULL, *stream3 = NULL, *stream3b = NULL, *stream3c = NULL, *stream3d = NULL, *stream4 = NULL, *stream5 = NULL, *stream5b = NULL, *stream5c = NULL, *stream5d = NULL,*stream6 = NULL, *stream7b = NULL, *stream7c = NULL, *stream7d = NULL, *stream1b = NULL, *stream1c = NULL, *stream1d = NULL;
        if (cam2xyz != NULL) {
            stream1 = open_log_file(report_dir, "target-mtx.dat");
            stream1b = open_log_file(report_dir, "target-mtxve.dat");
            stream1c = open_log_file(report_dir, "target-mtxve2.dat");
            stream1d = open_log_file(report_dir, "target-mtxve3.dat");
        }
        if (lut != NULL) {
            stream2 = open_log_file(report_dir, "target-nve-lut.dat");
            stream3 = open_log_file(report_dir, "target-nve-lutvm.dat");
            stream3b = open_log_file(report_dir, "target-nve-lutve.dat");
            stream3c = open_log_file(report_dir, "target-nve-lutve2.dat");
            stream3d = open_log_file(report_dir, "target-nve-lutve3.dat");
        }
        if (dcp_lut && dcp_lut->hsm) {
            stream4 = open_log_file(report_dir, "target-hsm-lut.dat");
            stream5 = open_log_file(report_dir, "target-hsm-lutvm.dat");
            stream5b = open_log_file(report_dir, "target-hsm-lutve.dat");
            stream5c = open_log_file(report_dir, "target-hsm-lutve2.dat");
            stream5d = open_log_file(report_dir, "target-hsm-lutve3.dat");
        }
        cmsHTRANSFORM tfm = NULL;
        if (icc != NULL && icc->lut != NULL) {
            stream6 = open_log_file(report_dir, "target-icc-lut.dat");
            stream7b = open_log_file(report_dir, "target-icc-lutve.dat");
            stream7c = open_log_file(report_dir, "target-icc-lutve2.dat");
            stream7d = open_log_file(report_dir, "target-icc-lutve3.dat");
            tfm = icc_lookup_fast_init(icc);
        }
        struct patch_set *ps = normalize_patch_set(ps_, ps_white_idx);
        for (int i = 0; i < ps->patch_count; i++) {
            v3 xyz = ps->patch[i].xyz;
            v3 cxyz = cam2xyz == NULL ? (v3){{0,0,0}} : m3x3_mul_v3(*cam2xyz, ps->patch[i].cam);
            v3 mtx = xyz2lutspace(cxyz);
            v3 ref = xyz2lutspace(xyz);
            v3 de;
            if (stream1) {
                fprintf(stream1, "%f %f %f\n", mtx.v[1], mtx.v[2], mtx.v[0]);
                fprintf(stream1b, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], mtx.v[1]-ref.v[1], mtx.v[2]-ref.v[2], mtx.v[0]-ref.v[0]);
                de = de_vector(ref, mtx, *ps_wp, false);
                fprintf(stream1c, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
                de = de_vector(ref, mtx, *ps_wp, true);
                fprintf(stream1d, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
            }
            if (stream2) {
                v3 ltu = xyz2lutspace(chromalut_lookup(lut, cxyz));
                fprintf(stream2, "%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
                fprintf(stream3, "%f %f %f %f %f %f\n", mtx.v[1], mtx.v[2], mtx.v[0], ltu.v[1]-mtx.v[1], ltu.v[2]-mtx.v[2], ltu.v[0]-mtx.v[0]);
                fprintf(stream3b, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], ltu.v[1]-ref.v[1], ltu.v[2]-ref.v[2], ltu.v[0]-ref.v[0]);
                de = de_vector(ref, ltu, *ps_wp, false);
                fprintf(stream3c, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
                de = de_vector(ref, ltu, *ps_wp, true);
                fprintf(stream3d, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
            }
            if (stream4) {
                v3 ltu = xyz2lutspace(dcp_lookup(dcp_lut->hsm, dcp_lut->hsmdims, dcp_lut->srgb_gamma, cxyz, 1, NULL, 0, NULL, NULL));
                fprintf(stream4, "%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
                fprintf(stream5, "%f %f %f %f %f %f\n", mtx.v[1], mtx.v[2], mtx.v[0], ltu.v[1]-mtx.v[1], ltu.v[2]-mtx.v[2], ltu.v[0]-mtx.v[0]);
                fprintf(stream5b, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], ltu.v[1]-ref.v[1], ltu.v[2]-ref.v[2], ltu.v[0]-ref.v[0]);
                de = de_vector(ref, ltu, *ps_wp, false);
                fprintf(stream5c, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
                de = de_vector(ref, ltu, *ps_wp, true);
                fprintf(stream5d, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
            }
            if (stream6) {
                v3 cam_wb = {{ ps->patch[i].cam.v[0] / wb->v[0], ps->patch[i].cam.v[1] / wb->v[1], ps->patch[i].cam.v[2] / wb->v[2] }};
                v3 ltu = xyz2lutspace(icc_lookup_fast(tfm, icc->pcs_is_lab, cam_wb));
                fprintf(stream6, "%f %f %f\n", ltu.v[1], ltu.v[2], ltu.v[0]);
                fprintf(stream7b, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], ltu.v[1]-ref.v[1], ltu.v[2]-ref.v[2], ltu.v[0]-ref.v[0]);
                de = de_vector(ref, ltu, *ps_wp, false);
                fprintf(stream7c, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
                de = de_vector(ref, ltu, *ps_wp, true);
                fprintf(stream7d, "%f %f %f %f %f %f\n", ref.v[1], ref.v[2], ref.v[0], de.v[1], de.v[2], de.v[0]);
            }
        }
        if (tfm != NULL) {
            cmsDeleteTransform(tfm);
        }
        if (stream1) {
            fclose(stream1);
            fclose(stream1b);
            fclose(stream1c);
            fclose(stream1d);
        }
        if (stream2) {
            fclose(stream2);
            fclose(stream3);
            fclose(stream3b);
            fclose(stream3c);
            fclose(stream3d);
        }
        if (stream4) {
            fclose(stream4);
            fclose(stream5);
            fclose(stream5b);
            fclose(stream5c);
            fclose(stream5d);
        }
        if (stream6) {
            fclose(stream6);
            fclose(stream7b);
            fclose(stream7c);
            fclose(stream7d);
        }

        if (dcam_ssf != NULL) {
            /*
unset key
set palette rgbformula 30,31,32
set cbrange [0:500]
plot 'ssf-csep.dat' pt 7 ps 2 lt palette, 'gmt-locus.dat' using 1:2:4 w l lw 4 lc rgb var, 'gmt-srgb.dat' w l, 'gmt-prophoto.dat' w l, 'gmt-pointer.dat' using 1:2:4 w l lw 2 lc rgb var
            */
            double max_uvdist;
            {
                double *uvdists = malloc(ps->patch_count * sizeof(*uvdists));
                for (int i = 0; i < ps->patch_count; i++) {
                    v3 xyz = ps->patch[i].xyz;
                    v3 ltu = xyz2lutspace(xyz);
                    double mindist = 10000;
                    for (int j = 0; j < ps->patch_count; j++) {
                        if (i == j) continue;
                        v3 ltu2 = xyz2lutspace(ps->patch[j].xyz);
                        double uvdist = sqrt((ltu2.v[2]-ltu.v[2])*(ltu2.v[2]-ltu.v[2]) + (ltu2.v[1]-ltu.v[1])*(ltu2.v[1]-ltu.v[1]));
                        if (uvdist < mindist) mindist = uvdist;
                    }
                    uvdists[i] = mindist;
                }
                qsort(uvdists, ps->patch_count, sizeof(double), double_cmp);
                double median_min_dist = uvdists[ps->patch_count/2];
                free(uvdists);
                max_uvdist = median_min_dist * 1.9;
            }

            FILE *stream1 = open_log_file(report_dir, "ssf-csep.dat");
            for (int i = 0; i < ps->patch_count; i++) {
                v3 xyz = ps->patch[i].xyz;
                v3 ltu = xyz2lutspace(xyz);
                v3 cam = ps->patch[i].cam;
                double minde = 1000;
                double mindiff = 10000000000;
                int test_count = 0;
                for (int j = 0; j < ps->patch_count; j++) {
                    if (i == j) continue;
                    {
                        v3 ltu2 = xyz2lutspace(ps->patch[j].xyz);
                        double uvdist = sqrt((ltu2.v[2]-ltu.v[2])*(ltu2.v[2]-ltu.v[2]) + (ltu2.v[1]-ltu.v[1])*(ltu2.v[1]-ltu.v[1]));
                        if (uvdist > max_uvdist) {
                            continue;
                        }
                    }
                    double diff = 0;
                    v3 cam1 = v3_norm(cam);
                    v3 cam2 = v3_norm(ps->patch[j].cam);
                    for (int j = 0; j < 3; j++) {
                        diff += fabs(cam1.v[j] - cam2.v[j]);
                    }
                    diff *= v3_max(cam);
                    diff /= 3;
                    double de = norm_l_de(xyz, ps->patch[j].xyz, *ps_wp);
                    if (de < minde) {
                        minde = de;
                    }
                    diff /= de;
                    if (diff < mindiff) {
                        mindiff = diff;
                    }
                    test_count++;
                }
                if (minde == 1000) {
                    continue;
                }
                // for one delta E chromaticity difference, the signal changes this much
                double sig_diff_per_de = mindiff;

                double csep = sig_diff_per_de * 65536;

                fprintf(stream1, "%f %f %f\n", ltu.v[1], ltu.v[2], csep);
            }
            fclose(stream1);
        }
        target_delete_patch_set(ps);
    }
}

static struct profio_dcp *
make_dcp(struct dcam_profile *prof,
         const char *camera_model,
         bool make_lut,
         int hcount,
         int scount,
         int vcount,
         bool enable_3d_hsm,
         bool skip_fwd,
         bool skip_lut_mtx,
         bool skip_observer_mapping,
         const double tc[],
         int tc_len,
         enum gc_type gc_type,
         enum tc_type tc_type,
         const struct tone_rep_op_config *ntro_conf,
         bool skip_lut_gamma,
         bool allow_discontinuity_hue_shifts,
         const char report_dir[])
{
    // test whitepoint position
    const v3 ref_d50 = dcp_d50();
    v3 map_d50 = ref_d50;
    m3x3 forward_matrix = (make_lut && !skip_lut_mtx) ? prof->lut_matrix : prof->forward_matrix;
    m3x3 map_white = v3_asdiagonal((v3){{1,1,1}});
    {
        v3 prof_d50 = m3x3_mul_v3(forward_matrix, (v3){{1,1,1}});
        double de = ciede2000(prof_d50, ref_d50, ref_d50);
        if (de > 0.1 && skip_observer_mapping) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Warning: forward matrix D50 whitepoint differs %.2f DE from the 1931_2 D50 used\n"
                "  by most DCP compatible products. This probably means that an observer other\n"
                "  than 1931_2 was used for making the profile. The profile may therefore not\n"
                "  work in some DCP products, and may not produce the expected colors.\n",
                de);
        }
        // if de is small enough it's probably just rounding errors and we make a final adjustment to "perfect" anyway
        if (!skip_observer_mapping || de <= 0.1) {
            if (de > 0.1) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "The D50 whitepoint in the forward matrix differs %.2f DE from the 1931_2 D50\n"
                    "  used by most DCP compatible products. This probably means that an observer\n"
                    "  other than 1931_2 was used for making the profile. Applying a Bradford\n"
                    "  transform to compensate.\n", de);
            }
            map_d50 = prof_d50;
            map_white = bradford_map_white(map_d50, ref_d50);
        }
    }

    struct profio_dcp *dcp = calloc(1, sizeof(*dcp));
    // we don't remap color matrix, not used for LUT and too small difference to care, seems to make more damage than improvement
    dcp->cm[0] = prof->color_matrix;
    dcp->fm[0] = m3x3_mul_m3x3(map_white, forward_matrix);
    dcp->has.fm = !skip_fwd;
    if (make_lut && chromalut_isnoop(prof->lut)) {
        if (enable_3d_hsm) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Native LUT has no corrections, but a HueSatMap is still made as it's 3D and can\n"
                "  thus encode extreme value gamut compression.\n");
        } else if (!skip_lut_mtx && !m3x3_isequal(prof->lut_matrix, prof->forward_matrix)) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Native LUT has no corrections, but a HueSatMap is still made as it's used for\n"
                "  expanding the LUT matrix.\n");
        } else {
            // makes no sense to make a HueSatMap
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Native LUT has no corrections, skipping HueSatMap.\n");
            make_lut = false;
        }
    }
    if (make_lut) {
        const int hsm_vcount = enable_3d_hsm ? vcount : 1;
        if (hsm_vcount == 1) {
            elog("Generating 2.5D HueSatMap with %dx%d = %d entries...", hcount, scount, hcount*scount);
        } else {
            elog("Generating 3D HueSatMap with %dx%dx%d = %d entries...", hcount, scount, hsm_vcount, hcount*scount*hsm_vcount);
        }
        dcp->huesatmap[0] = dnglut_huesatmap_new(dcp->hsmdims, hcount, scount, hsm_vcount, !skip_lut_gamma, prof->lut, map_d50, prof->forward_matrix, forward_matrix);
        elog("done!\n");
        if (!skip_lut_gamma && dcp->hsmdims[2] > 1) {
            dcp->huesatmap_srgb_gamma = true;
            dcp->has.huesatmap_srgb_gamma = true;
        }
        dnglut_test_discontinuity(dcp->huesatmap[0], dcp->hsmdims[0], dcp->hsmdims[1], dcp->hsmdims[2], "HueSatMap", allow_discontinuity_hue_shifts);
    }
    if (camera_model != NULL) {
        dcp->camera_model = strdup(camera_model);
    } else if (prof->camera_name != NULL) {
        dcp->camera_model = strdup(prof->camera_name);
    } else {
        dcp->camera_model = strdup("Unknown Camera");
        elog("Warning: unique camera name should be set to exif make/model\n");
    }
    dcp->profile_name = strdup("Generated by DCamProf v" DCAMPROF_VERSION);
    dcp->copyright = strdup("Copyright, the creator of this profile (generated by DCamProf v" DCAMPROF_VERSION ")");
    dcp->embed_policy = pepNoRestrictions;
    dcp->has.embed_policy = true;
    dcp->ill[0] = prof->illuminant_name;
    dcp->ill_count = 1;
    dcp->black_render_none = true;
    dcp->has.black_render_none = true;

    if (tc != NULL) {
        switch (tc_type) {
        case TC_NEUTRAL:
            if ((tc_len > 2 || gc_type != GC_NONE) && make_lut) {
                if (!skip_lut_gamma) {
                    dcp->looktable_srgb_gamma = true;
                    dcp->has.looktable_srgb_gamma = true;
                }
                double cv, cs;
                look_neutral_tone_rep_op_tc_analysis(tc, tc_len, &cv, &cs);
                elog("The tone curve's contrast value is %.2f (=> auto chroma scaling value %.3f)\n", cv, cs);
                dcp->looktable = dnglut_looktable_new(dcp->lookdims, hcount, scount, vcount, tc, tc_len, true, dcp->looktable_srgb_gamma, gc_type, ntro_conf);
                dnglut_test_discontinuity(dcp->looktable, hcount, scount, vcount, "LookTable", allow_discontinuity_hue_shifts);
                print_hsm(report_dir, "lkt", dcp->looktable, dcp->lookdims[0], dcp->lookdims[1], dcp->lookdims[2], dcp->looktable_srgb_gamma);
            }
            dcp->tonecurve = malloc(2 * tc_len * sizeof(dcp->tonecurve[0]));
            for (int i = 0; i < tc_len; i++) {
                dcp->tonecurve[2*i+0] = i / (double)(tc_len - 1);
                dcp->tonecurve[2*i+1] = tc[i];
            }
            dcp->tonecurve_len = tc_len;
            break;
        case TC_STANDARD:
            dcp->tonecurve = malloc(2 * tc_len * sizeof(dcp->tonecurve[0]));
            for (int i = 0; i < tc_len; i++) {
                dcp->tonecurve[2*i+0] = i / (double)(tc_len - 1);
                dcp->tonecurve[2*i+1] = tc[i];
            }
            dcp->tonecurve_len = tc_len;
            break;
        }
    }

    struct dcp_lut dcp_lut;
    dcp_lut.hsm = dcp->huesatmap[0];
    dcp_lut.hsmdims = dcp->hsmdims;
    dcp_lut.srgb_gamma = dcp->huesatmap_srgb_gamma;
    print_lookup_reports(report_dir, prof->lut, &dcp_lut, NULL, NULL, NULL, NULL, NULL, -1, NULL);
    print_tone_curve(report_dir, tc, tc_len, false);

    return dcp;
}

enum icc_profile_type {
    ICC_PROFILE_TYPE_UNDEFINED = 0,
    ICC_PROFILE_TYPE_LABLUT,
    ICC_PROFILE_TYPE_XYZLUT,
    ICC_PROFILE_TYPE_MATRIX
};

static struct profio_icc *
make_icc(struct dcam_profile *prof,
         const char *description,
         enum icc_profile_type icc_type,
         bool apply_wb,
         int clut_side,
         double *trc_[3],
         int trc_count_,
         const double tc[],
         int tc_len,
         enum gc_type gc_type,
         enum tc_type tc_type,
         bool skip_tc_apply,
         const struct tone_rep_op_config *ntro_conf,
         const char report_dir[])
{
    const v3 ref_d50 = icc_d50();
    m3x3 forward_matrix = prof->forward_matrix;
    v3 prof_d50 = m3x3_mul_v3(forward_matrix, (v3){{1,1,1}});
    double de = ciede2000(prof_d50, ref_d50, ref_d50);
    if (de > 0.1) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "The D50 whitepoint in the forward matrix differs %.2f DE from the ICC D50\n"
            "  This probably means that an observer other than 1931_2 was used for making\n"
            "  the profile. Applying a Bradford transform to compensate.\n", de);
    }

    double *trc[3] = { NULL, NULL, NULL };
    int trc_count = 0;
    if (trc_count_ > 0) {
        trc_count = 4096;
        if (trc_count_ > trc_count) {
            elog("Provided transfer function has %d entries, downsampling to %d (as ICCv2 can't handle more).\n", trc_count_, trc_count);
            for (int i = 0; i < 3; i++) {
                trc[i] = malloc(trc_count * sizeof(trc[i][0]));
                for (int j = 0; j < trc_count; j++) {
                    double x = (double)j / (trc_count-1);
                    trc[i][j] = trc_[i][(int)round(x * (trc_count_-1))];
                }
            }
        } else {
            trc_count = trc_count_;
            for (int i = 0; i < 3; i++) {
                trc[i] = malloc(trc_count * sizeof(trc[i][0]));
                for (int j = 0; j < trc_count; j++) {
                    trc[i][j] = trc_[i][j];
                }
            }
        }
        double mx = trc[0][trc_count-1];
        if (trc[1][trc_count-1] > mx) mx = trc[1][trc_count-1];
        if (trc[2][trc_count-1] > mx) mx = trc[2][trc_count-1];
        if (mx != 1.0) {
            elog("Warning: provided transfer function ends at %f rather than 1.0, compensating.\n", mx);
            mx = 1.0 / mx;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < trc_count; j++) {
                    trc[i][j] *= mx;
                }
            }
        }
        double mn = trc[0][0];
        if (trc[1][0] < mn) mn = trc[1][0];
        if (trc[2][0] < mn) mn = trc[2][0];
        if (mn != 0.0) {
            elog("Warning: provided transfer function starts at %f rather than 0.0.\n", mn);
        }
    }

    struct profio_icc *icc = calloc(1, sizeof(*icc));
    if (description != NULL) {
        icc->description = strdup(description);
    } else if (prof->camera_name != NULL) {
        icc->description = strdup(prof->camera_name);
    } else {
        icc->description = strdup("Unknown Camera");
    }
    icc->copyright = strdup("Copyright, the creator of this profile (generated by DCamProf v" DCAMPROF_VERSION ")");
    icc->wp = v3_norm(prof->illuminant_wp);
    m3x3 fm = m3x3_mul_m3x3(bradford_map_white(prof_d50, ref_d50), forward_matrix);
    if (icc_type != ICC_PROFILE_TYPE_MATRIX) {
        icc->pcs_is_lab = (icc_type == ICC_PROFILE_TYPE_LABLUT);
        v3 wb = apply_wb ? prof->forward_matrix_wb : (v3){{1,1,1}};
        icc->lut = icclut_new(fm, prof_d50, ref_d50, wb, trc, trc_count, tc, tc_len, gc_type, tc_type, skip_tc_apply, prof->lut, icc->pcs_is_lab, clut_side, ntro_conf);
    } else {
        if (tc != NULL) {
            elog("Error: tone-curve not supported for matrix-only profiles.\n");
            exit(EXIT_FAILURE);
        }
        icc->pcs_is_lab = false;
        icc->fm = fm;
        if (apply_wb) {
            v3 wb_mul = v3_invert(prof->forward_matrix_wb);
            icc->fm = m3x3_mul_m3x3(icc->fm, v3_asdiagonal(wb_mul));
            // make sure sum of Y is still 1.0 after rebalancing
            icc->fm = k_mul_m3x3(1.0 / (icc->fm.v[1][0] + icc->fm.v[1][1] + icc->fm.v[1][2]), icc->fm);
        }
        if (trc_count > 0) {
            for (int i = 0; i < 3; i++) {
                icc->trc[i] = malloc(sizeof(icc->trc[i][0]) * trc_count);
                for (int j = 0; j < trc_count; j++) {
                    icc->trc[i][j] = trc[i][j];
                }
            }
            icc->trc_len = trc_count;
        } else {
            for (int i = 0; i < 3; i++) {
                icc->trc[i] = malloc(sizeof(icc->trc[i][0]) * 1);
                icc->trc[i][0] = 1;
            }
            icc->trc_len = 1;
        }
    }
    free(trc[0]);
    free(trc[1]);
    free(trc[2]);
    print_lookup_reports(report_dir, prof->lut, NULL, icc, NULL, NULL, NULL, NULL, -1, NULL);
    print_tone_curve(report_dir, tc, tc_len, false);
    return icc;
}

static void
elog_matrix(const char indent[],
            const char name[],
            m3x3 m)
{
    elog("%s\"%s\": [\n"
         "%s  [ %9.6f, %9.6f, %9.6f ],\n"
         "%s  [ %9.6f, %9.6f, %9.6f ],\n"
         "%s  [ %9.6f, %9.6f, %9.6f ]\n"
         "%s]",
         indent, name,
         indent, m.v[0][0], m.v[0][1], m.v[0][2],
         indent, m.v[1][0], m.v[1][1], m.v[1][2],
         indent, m.v[2][0], m.v[2][1], m.v[2][2],
         indent);
}

static struct dcam_profile *
make_profile(enum exif_lightsource calib_ls,
             const v3 calib_wp,
             const spectrum_t *calib_spec,
             const m3x3 *precalc_cm,
             const m3x3 *precalc_fm,
             const m3x3 *precalc_lm,
             bool use_slow_matrix_optimizer,
             double min_chroma_dist,
             double compress_factor,
             const int white_idx,
             const struct patch_set *ps_,
             const struct patch_set *psD50,
             const struct observer *obs,
             const struct dcam_ssf *dcam_ssf,
             const double smallest_fm[3],
             bool skip_lut_in_report,
             bool refine_only_forward_matrix,
             bool exit_on_mtx_refine_failure,
             const char report_dir[])
{
    const v3 d50 = glob.d50;
    v3 wb, inv_wb;
    m3x3 cm;
    {
        if (precalc_cm != NULL) {
            cm = *precalc_cm;
            elog("Imported ColorMatrix:\n  {\n");
            jsonio_3x3matrix_print(stderr, "    ", "ColorMatrix1", cm);
            elog("\n  }\n");
        } else {
            struct patch_set *ps = normalize_patch_set(ps_, white_idx);
            m3x3 cam2xyz = find_cam2xyz_matrix("calibration", calib_ls, calib_wp, white_idx, ps, use_slow_matrix_optimizer, NULL, false, !refine_only_forward_matrix, exit_on_mtx_refine_failure);

            elog("Inverting to get ColorMatrix:\n  {\n");
            m3x3 xyz2cam = m3x3_invert(cam2xyz);
            elog_matrix("    ", "ColorMatrix1", xyz2cam);
            elog("\n  }\n");
            cm = xyz2cam;
            print_pm_errors(m3x3_invert(cm), NULL, ps, calib_wp, false, false, report_dir, "cm-patch-errors");
            target_delete_patch_set(ps);
        }
    }

    wb = v3_norm(m3x3_mul_v3(cm, calib_wp));
    inv_wb = (v3){{ 1.0/wb.v[0], 1.0/wb.v[1], 1.0/wb.v[2] }};
    elog("ColorMatrix optimal white balance for target: %g,%g,%g (m%g,%g,%g)\n", wb.v[0], wb.v[1], wb.v[2], inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);

    m3x3 lm;
    {
        if (precalc_lm != NULL) {
            elog("Imported LUTMatrix:\n");
            lm = *precalc_lm;
            elog("  {\n");
            elog_matrix("    ", "LUTMatrix1", lm);
            elog("\n  }\n");
        } else {
            struct patch_set *ps = normalize_patch_set(psD50, white_idx);
            m3x3 cam2xyz = find_cam2xyz_matrix("connection space", lsD50, d50, white_idx, ps, use_slow_matrix_optimizer, NULL, true, false, false);

            v3 camd50 = m3x3_mul_v3(m3x3_invert(cam2xyz), d50);
            wb = v3_norm(camd50);
            v3 mul;
            mul.v[0] = 1.0 / camd50.v[0];
            mul.v[1] = 1.0 / camd50.v[1];
            mul.v[2] = 1.0 / camd50.v[2];
            cam2xyz = m3x3_mul_m3x3(cam2xyz, m3x3_invert(v3_asdiagonal(mul)));
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].cam.v[0] /= camd50.v[0];
                ps->patch[i].cam.v[1] /= camd50.v[1];
                ps->patch[i].cam.v[2] /= camd50.v[2];
            }

            elog("  {\n");
            elog_matrix("    ", "LUTMatrix1", cam2xyz);
            elog("\n  }\n");

            lm = cam2xyz;

            inv_wb = (v3){{ 1.0/wb.v[0], 1.0/wb.v[1], 1.0/wb.v[2] }};
            elog("LUTMatrix optimal white balance for target: %g,%g,%g (m%g,%g,%g)\n", wb.v[0], wb.v[1], wb.v[2], inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);
            print_pm_errors(lm, NULL, ps, d50, false, false, report_dir, "lm-patch-errors");
            target_delete_patch_set(ps);
        }
    }

    m3x3 fm;
    {
        if (precalc_fm != NULL) {
            elog("Imported ForwardMatrix:\n");
            fm = *precalc_fm;
            elog("  {\n");
            elog_matrix("    ", "ForwardMatrix1", fm);
            elog("\n  }\n");
        } else {
            struct patch_set *ps = normalize_patch_set(psD50, white_idx);
            m3x3 cam2xyz = find_cam2xyz_matrix("connection space", lsD50, d50, white_idx, ps, use_slow_matrix_optimizer, smallest_fm, false, true, exit_on_mtx_refine_failure);

            elog("Applying white-balance to get ForwardMatrix:\n");

            v3 camd50 = m3x3_mul_v3(m3x3_invert(cam2xyz), d50);
            wb = v3_norm(camd50);
            v3 mul;
            mul.v[0] = 1.0 / camd50.v[0];
            mul.v[1] = 1.0 / camd50.v[1];
            mul.v[2] = 1.0 / camd50.v[2];
            cam2xyz = m3x3_mul_m3x3(cam2xyz, m3x3_invert(v3_asdiagonal(mul)));
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].cam.v[0] /= camd50.v[0];
                ps->patch[i].cam.v[1] /= camd50.v[1];
                ps->patch[i].cam.v[2] /= camd50.v[2];
            }

            elog("  {\n");
            elog_matrix("    ", "ForwardMatrix1", cam2xyz);
            elog("\n  }\n");

            fm = cam2xyz;

            inv_wb = (v3){{ 1.0/wb.v[0], 1.0/wb.v[1], 1.0/wb.v[2] }};
            elog("ForwardMatrix optimal white balance for target: %g,%g,%g (m%g,%g,%g)\n", wb.v[0], wb.v[1], wb.v[2], inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);

            print_pm_errors(fm, NULL, ps, d50, false, false, report_dir, "fm-patch-errors");
            target_delete_patch_set(ps);
        }
    }

    const m3x3 cam2xyzD50 = m3x3_mul_m3x3(fm, m3x3_invert(v3_asdiagonal(wb)));

    struct patch_set *ps1 = normalize_patch_set(psD50, white_idx);
    chromalut_t *lut = chromalut_new(fm, wb, d50, ps1, min_chroma_dist, compress_factor);

    {
        print_pm_errors(cam2xyzD50, skip_lut_in_report ? NULL : lut, ps1, d50, true, true, report_dir, "patch-errors");

        print_lookup_reports(report_dir, skip_lut_in_report ? NULL : lut, NULL, NULL, dcam_ssf, &cam2xyzD50, psD50, &d50, white_idx, NULL);
        target_print("targetd50", ps1, d50, obs, report_dir);
    }
    target_delete_patch_set(ps1);

    struct dcam_profile *prof = calloc(1, sizeof(*prof));
    prof->illuminant_name = calib_ls;
    prof->illuminant_spec = spec_copy(calib_spec);
    prof->illuminant_wp = calib_wp;
    prof->color_matrix = cm;
    prof->forward_matrix = fm;
    prof->lut_matrix = lm;
    prof->forward_matrix_wb = wb;
    prof->lut = lut;

    print_matrices(report_dir, &cm, &fm, &lm, NULL, NULL, d50);
    print_profile_gamut(report_dir, &fm, lut, NULL, NULL, NULL, d50, NULL);
    return prof;
}

struct find_cct_fun_arg {
    const struct observer *obs;
    double u, v;
};

double
find_cct_fun(double temp,
             void *arg)
{
    struct find_cct_fun_arg *a = (struct find_cct_fun_arg *)arg;
    const struct observer *obs = a->obs;
    spectrum_t *bb = spectraldb_render_blackbody(temp, 300, 830, 5);
    v3 uvY = xyz2uvY(spec2xyz(bb, obs));
    free(bb);
    double de = sqrt((a->u-uvY.v[0])*(a->u-uvY.v[0]) + (a->v-uvY.v[1])*(a->v-uvY.v[1]));
    return de;
}

static double
find_cct(const struct observer *obs,
         double u,
         double v,
         double *diff_de)
{
    struct find_cct_fun_arg arg = { .obs = obs, .u = u, .v = v };
    double best_temp = interval_halving_min(find_cct_fun, &arg, 500, 100000, 0.5, 1000);
    if (diff_de != NULL) {
        spectrum_t *bb = spectraldb_render_blackbody(best_temp, 300, 830, 5);
        v3 xyz1 = v3_norm(spec2xyz(bb, obs));
        v3 xyz2 = v3_norm(uv2xyz_lightest(u, v, (v3){{1,1,1}}));
        *diff_de = ciede2000(xyz1, xyz2, xyz1);
        free(bb);
    }
    return best_temp;
}

enum prof_type {
    PROF_NATIVE = 0,
    PROF_ICC,
    PROF_DCP,
};

struct prof_wrap {
    enum prof_type type;
    union {
        struct profio_dcp *dcp;
        struct dcam_profile *prof;
        struct profio_icc *icc;
    } p;
};

static void
test_profile(const struct prof_wrap *pw,
             const v3 *custom_wb,
             bool skip_lut,
             bool skip_look_lut,
             const v3 test_wp,
             const spectrum_t *test_spec,
             const int w_idx,
             const struct patch_set *ps_,
             const struct patch_set *psD50_,
             const struct observer *obs,
             const struct dcam_ssf *dcam_ssf,
             double *trc[3],
             int trc_count,
             const struct rgbimg *rgbimg,
             const char output_image_filename[],
             const char report_dir[])
{
    const v3 d50 = glob.d50;
    uint32_t hsmdims[3];
    chromalut_t *lut = NULL;
    struct dcp_lut dcp_lut;
    memset(&dcp_lut, 0, sizeof(dcp_lut));
    bool no_wb = false;
    m3x3 cam2xyzD50;
    v3 wb;
    if (custom_wb != NULL) {
        wb = v3_norm(*custom_wb);
        v3 inv_wb = v3_invert(wb);
        elog("Using user-provided white balance for target: %g,%g,%g (m%g,%g,%g)\n", wb.v[0], wb.v[1], wb.v[2], inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);
    } else if (psD50_ != NULL) {
        wb.v[0] = 1;
        wb.v[1] = 1;
        wb.v[2] = 1;
        elog("Whitest patch in target has index %d.\n", w_idx);
        wb = v3_norm(psD50_->patch[w_idx].cam);
        v3 cam_wb;
        v3 xyz_approx_D50 = v3_norm(psD50_->patch[w_idx].xyz);
        switch (pw->type) {
        case PROF_NATIVE: {
            //test_wb = v3_norm(m3x3_mul_v3(pw->p.prof->color_matrix, test_wp));
            m3x3 xyz2cam_wb = m3x3_invert(pw->p.prof->forward_matrix);
            cam_wb = v3_norm(m3x3_mul_v3(xyz2cam_wb, xyz_approx_D50));
            wb.v[0] = psD50_->patch[w_idx].cam.v[0] / cam_wb.v[0];
            wb.v[1] = psD50_->patch[w_idx].cam.v[1] / cam_wb.v[1];
            wb.v[2] = psD50_->patch[w_idx].cam.v[2] / cam_wb.v[2];
            wb = v3_norm(wb);
            break;
        }
        case PROF_DCP: {
            v3 old_wb = wb;
            for (int i = 0; i < 1000; i++) {
                m3x3 cam2xyzD50 = dcp_make_cam2xyz_d50(pw->p.dcp, NULL, NULL, NULL, NULL, wb);
                v3 invwb = {{ 1.0/wb.v[0], 1.0/wb.v[1], 1.0/wb.v[2] }};
                m3x3 fm = m3x3_mul_m3x3(cam2xyzD50, m3x3_invert(v3_asdiagonal(invwb)));
                //elog_matrix("", "", fm);
                m3x3 xyz2cam_wb = m3x3_invert(fm);
                cam_wb = v3_norm(m3x3_mul_v3(xyz2cam_wb, xyz_approx_D50));
                wb.v[0] = psD50_->patch[w_idx].cam.v[0] / cam_wb.v[0];
                wb.v[1] = psD50_->patch[w_idx].cam.v[1] / cam_wb.v[1];
                wb.v[2] = psD50_->patch[w_idx].cam.v[2] / cam_wb.v[2];
                wb = v3_norm(wb);
                //elog("%f\n", euclidean_dist(old_wb, wb));
                if (euclidean_dist(old_wb, wb) < 1e-7) {
                    break;
                }
                old_wb = wb;
            }
            break;
        }
        case PROF_ICC: {
            const struct profio_icc *icc = pw->p.icc;
            if (icc->lut == NULL) {
                m3x3 xyz2cam_wb = m3x3_invert(icc->fm);
                cam_wb = v3_norm(m3x3_mul_v3(xyz2cam_wb, xyz_approx_D50));
                wb.v[0] = psD50_->patch[w_idx].cam.v[0] / cam_wb.v[0];
                wb.v[1] = psD50_->patch[w_idx].cam.v[1] / cam_wb.v[1];
                wb.v[2] = psD50_->patch[w_idx].cam.v[2] / cam_wb.v[2];
                wb = v3_norm(wb);
            } else {
                v3 cam_wb = icc_reverse_lookup(icc, xyz_approx_D50);
                if (cam_wb.v[0] < 0) {
                    elog("Error: failed to find target white balance through a reverse ICC lookup. Provide white balance manually instead.\n");
                    exit(EXIT_FAILURE);
                }
                wb.v[0] = psD50_->patch[w_idx].cam.v[0] / cam_wb.v[0];
                wb.v[1] = psD50_->patch[w_idx].cam.v[1] / cam_wb.v[1];
                wb.v[2] = psD50_->patch[w_idx].cam.v[2] / cam_wb.v[2];
                wb = v3_norm(wb);
            }
            break;
        }
        }
        v3 inv_wb = v3_invert(wb);
        elog("Optimal white balance to match target: %g,%g,%g (m%g,%g,%g)\n", wb.v[0], wb.v[1], wb.v[2], inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);
    } else if (test_spec != NULL && dcam_ssf != NULL) {
        wb = v3_norm(spec2xyz(test_spec, &dcam_ssf->ssf));
        elog("Calculated white balance from illuminant spectra and camera SSFs: %g,%g,%g\n", wb.v[0], wb.v[1], wb.v[2]);
    } else {
        wb.v[0] = 1;
        wb.v[1] = 1;
        wb.v[2] = 1;
        no_wb = true;
    }
    switch (pw->type) {
    case PROF_NATIVE:
        if (!skip_lut) {
            lut = pw->p.prof->lut;
        }
        cam2xyzD50 = m3x3_mul_m3x3(pw->p.prof->forward_matrix, m3x3_invert(v3_asdiagonal(wb)));
        print_matrices(report_dir, &pw->p.prof->color_matrix, &pw->p.prof->forward_matrix, &pw->p.prof->lut_matrix, NULL, NULL, d50);
        break;
    case PROF_DCP:
        cam2xyzD50 = dcp_make_cam2xyz_d50(pw->p.dcp, NULL, NULL, NULL, NULL, wb);
        if (pw->p.dcp->ill_count == 2 && no_wb) {
            elog("Error: could not derive white balance which is required for the dual-illuminant DCP to work.\n");
            exit(EXIT_FAILURE);
        }
        if (!skip_lut) {
            dcp_lut.hsm = dcp_make_lut(pw->p.dcp, NULL, NULL, NULL, NULL, wb);
            if (dcp_lut.hsm != NULL) {
                memcpy(hsmdims, pw->p.dcp->hsmdims, sizeof(hsmdims));
                dcp_lut.hsmdims = hsmdims;
                dcp_lut.srgb_gamma = pw->p.dcp->huesatmap_srgb_gamma;
                if (!skip_look_lut && pw->p.dcp->looktable != NULL && psD50_ != NULL) {
                    elog("The DCP has both HueSatMap and LookTable, skipping LookTable in target matching.\n");
                }
            } else if (!skip_look_lut && pw->p.dcp->looktable != NULL && psD50_ != NULL) {
                elog("The DCP has no HueSatMap but has a LookTable, applying that instead for target matching.\n");
                size_t count = pw->p.dcp->lookdims[0]*pw->p.dcp->lookdims[1]*pw->p.dcp->lookdims[2];
                dcp_lut.hsm = malloc(count * sizeof(*dcp_lut.hsm));
                memcpy(dcp_lut.hsm, pw->p.dcp->looktable, count * sizeof(*dcp_lut.hsm));
                memcpy(hsmdims, pw->p.dcp->lookdims, sizeof(hsmdims));
                dcp_lut.hsmdims = hsmdims;
                dcp_lut.srgb_gamma = pw->p.dcp->looktable_srgb_gamma;
            }
        }
        if (pw->p.dcp->tonecurve != NULL) {
            if (psD50_ != NULL) {
                elog("The DCP has a tone curve (not used in target matching).\n");
            }
            print_tone_curve(report_dir, pw->p.dcp->tonecurve, pw->p.dcp->tonecurve_len, true);
        }
        if (pw->p.dcp->has.baseline_exposure_offset) {
            if (psD50_ != NULL) {
                elog("The DCP has a baseline exposure offset (not used in target matching).\n");
            }
        }
        print_matrices(report_dir, &pw->p.dcp->cm[0], pw->p.dcp->has.fm ? &pw->p.dcp->fm[0] : NULL, NULL, pw->p.dcp->ill_count == 2 ? &pw->p.dcp->cm[1] : NULL, (pw->p.dcp->has.fm && pw->p.dcp->ill_count == 2) ? &pw->p.dcp->fm[1] : NULL, dcp_d50());
        break;
    case PROF_ICC:
        if (skip_lut && pw->p.icc->lut != NULL) {
            elog("Error: cannot skip LUT with ICC Lab LUT profile.\n");
            exit(EXIT_FAILURE);
        }
        if (pw->p.icc->lut == NULL) {
            cam2xyzD50 = m3x3_mul_m3x3(pw->p.icc->fm, m3x3_invert(v3_asdiagonal(wb)));
        }
        break;
    }

    double cct = 0;
    if (pw->type != PROF_ICC && !no_wb) {
        // white balance white point
        v3 wb_wp;
        switch (pw->type) {
        case PROF_NATIVE:
            wb_wp = v3_norm(m3x3_mul_v3(m3x3_invert(pw->p.prof->color_matrix), wb));
            break;
        case PROF_DCP: {
            double xy[2];
            dcp_neutral2xy(pw->p.dcp, NULL, NULL, NULL, NULL, wb, xy);
            wb_wp = xy2xyz_lightest(xy[0], xy[1], (v3){{1,1,1}});
            break;
        }
        case PROF_ICC:
            abort();
            break;
        }
        v3 xyY = xyz2xyY(wb_wp);
        double temp, tint;
        double whiteXY[2] = { xyY.v[0], xyY.v[1] };
        dcp_xy2temp(whiteXY, &temp, &tint);
        v3 uvY = xyz2uvY(wb_wp);
        v3 invwb = v3_norm(v3_invert(wb));
        double cct_de;
        cct = find_cct(obs, uvY.v[0], uvY.v[1], &cct_de);
        elog("White balance whitepoint according to profile:\n"
             "  raw RGB: %.4f, %.4f, %.4f\n"
             "      XYZ: %.4f, %.4f, %.4f\n"
             "       xy: %.4f %.4f\n"
             "     u'v': %.4f %.4f\n"
             " DNG temp: %uK (tint %d)\n"
             "      CCT: %uK (%.2f DE off-axis)\n",
             invwb.v[0], invwb.v[1], invwb.v[2],
             wb_wp.v[0], wb_wp.v[1], wb_wp.v[2], xyY.v[0], xyY.v[1], uvY.v[0], uvY.v[1],
             (int)temp, (int)tint,
             (int)cct, cct_de);
    }

    if (pw->type != PROF_ICC) {
        // test illuminant white balance
        v3 xyY = xyz2xyY(test_wp);
        double whiteXY[2] = { xyY.v[0], xyY.v[1] };
        v3 test_wb;
        switch (pw->type) {
        case PROF_NATIVE:
            test_wb = v3_norm(m3x3_mul_v3(pw->p.prof->color_matrix, test_wp));
            break;
        case PROF_DCP: {
            m3x3 xyz2cam = dcp_find_xyz2cam(pw->p.dcp, whiteXY, NULL, NULL, NULL, NULL, NULL, NULL);
            test_wb = v3_norm(m3x3_mul_v3(xyz2cam, test_wp));
            break;
        }
        case PROF_ICC:
            abort();
            break;
        }
        v3 inv_wb = v3_invert(test_wb);

        double temp, tint;
        dcp_xy2temp(whiteXY, &temp, &tint);
        v3 uvY = xyz2uvY(test_wp);
        double test_cct_de;
        v3 test_uvY = xyz2uvY(v3_norm(test_wp));
        double test_cct = find_cct(obs, test_uvY.v[0], test_uvY.v[1], &test_cct_de);
        elog("Test illuminant:\n"
             "      XYZ: %.4f, %.4f, %.4f\n"
             "       xy: %.4f %.4f\n"
             "     u'v': %.4f %.4f\n"
             " DNG temp: %uK (tint %d)\n"
             "      CCT: %uK (%.2f DE off-axis)\n"
             "Test illuminant white balance according to profile: %g,%g,%g (m%g,%g,%g)\n",
             test_wp.v[0], test_wp.v[1], test_wp.v[2], xyY.v[0], xyY.v[1], uvY.v[0], uvY.v[1],
             (int)temp, (int)tint,
             (int)test_cct, test_cct_de,
             test_wb.v[0], test_wb.v[1], test_wb.v[2],
             inv_wb.v[0], inv_wb.v[1], inv_wb.v[2]);
        if (cct != 0 && fabs(test_cct - cct) > 200 && dcam_ssf == NULL) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Warning: a large difference between test illuminant CCT and white balance CCT\n"
                "  may mean that the target's RGB values were generated with a different\n"
                "  illuminant than the the provided test illuminant. Either make sure the test\n"
                "  illuminant matches the RGB values in the file, or provide camera SSFs so the\n"
                "  RGB values can be regenerated from scratch.\n");
        }
    }

    struct patch_set *ps = NULL, *psD50 = NULL;
    if (ps_) ps = normalize_patch_set(ps_, w_idx);
    if (psD50_) {
        psD50 = normalize_patch_set(psD50_, w_idx);
        for (int i = 0; i < ps->patch_count; i++) {
            psD50->patch[i].cam = ps->patch[i].cam;
        }
    }

    target_print("target", ps, test_wp, obs, report_dir);
    target_print("targetd50", psD50, d50, obs, report_dir);

    uint32_t hsmdims_look[3];
    struct dcp_lut dcp_lut_look;
    memset(&dcp_lut_look, 0, sizeof(dcp_lut_look));
    if (pw->type == PROF_DCP && pw->p.dcp->looktable != NULL) {
        print_hsm(report_dir, "lkt", pw->p.dcp->looktable, pw->p.dcp->lookdims[0], pw->p.dcp->lookdims[1], pw->p.dcp->lookdims[2], pw->p.dcp->looktable_srgb_gamma);
        size_t count = pw->p.dcp->lookdims[0]*pw->p.dcp->lookdims[1]*pw->p.dcp->lookdims[2];
        dcp_lut_look.hsm = malloc(count * sizeof(*dcp_lut_look.hsm));
        memcpy(dcp_lut_look.hsm, pw->p.dcp->looktable, count * sizeof(*dcp_lut_look.hsm));
        memcpy(hsmdims_look, pw->p.dcp->lookdims, sizeof(hsmdims_look));
        dcp_lut_look.hsmdims = hsmdims_look;
        dcp_lut_look.srgb_gamma = pw->p.dcp->looktable_srgb_gamma;
    }
    if (report_dir != NULL) {
        elog("Writing lookup plots to \"%s\"...", report_dir);
        print_lookup_reports(report_dir,
                             lut,
                             pw->type == PROF_DCP ? &dcp_lut : NULL,
                             pw->type == PROF_ICC ? pw->p.icc : NULL,
                             dcam_ssf,
                             pw->type != PROF_ICC || pw->p.icc->lut == NULL ? &cam2xyzD50 : NULL,
                             psD50, &d50, w_idx, &wb);
        elog("\n");

        print_profile_gamut(report_dir, pw->type != PROF_ICC || pw->p.icc->lut == NULL ? &cam2xyzD50 : NULL,
                            lut,
                            pw->type == PROF_DCP ? &dcp_lut : NULL,
                            pw->type == PROF_DCP ? &dcp_lut : NULL,
                            pw->type == PROF_ICC ? pw->p.icc : NULL,
                            d50, &wb);
    }

    if (output_image_filename != NULL || report_dir != NULL) {

        double *tc = NULL;
        int tc_len = 0;
        if (pw->type == PROF_DCP && pw->p.dcp->tonecurve != NULL) {
            tc_len = pw->p.dcp->tonecurve_len;
            tc = malloc(tc_len * sizeof(tc[0]));
            for (int i = 0; i < tc_len; i++) {
                tc[i] = pw->p.dcp->tonecurve[2*i+1];
            }
        }
        double xyz_clip, alt_clip;
        elog("Writing test image...\n");
        render_test_tiff(report_dir,
                         output_image_filename,
                         pw->type != PROF_ICC || pw->p.icc->lut == NULL ? &cam2xyzD50 : NULL,
                         lut,
                         pw->type == PROF_DCP && pw->p.dcp->huesatmap[0] != NULL ? &dcp_lut : NULL,
                         pw->type == PROF_DCP && pw->p.dcp->looktable != NULL ? &dcp_lut_look : NULL,
                         pw->type == PROF_ICC ? pw->p.icc : NULL,
                         tc, tc_len,
                         skip_lut,
                         skip_look_lut,
                         trc, trc_count,
                         &wb,
                         false,
                         &xyz_clip, &alt_clip,
                         rgbimg);
        if (pw->type == PROF_DCP) {
            elog("  When writing test image %g%% clipped in ProPhoto matrix conversion, and %g%% clipped in LUT.\n", 100*xyz_clip, 100*alt_clip);
        } else if (pw->type == PROF_NATIVE) {
            elog("  When writing test image %g%% was mapped to clipping XYZ coordinates, %g%% clipped ProPhoto.\n", 100*xyz_clip, 100*alt_clip);
        }
        free(tc);
    }
    free(dcp_lut_look.hsm);

    if (psD50 != NULL) {
        elog("Testing target matching...");
        struct patch_error_report *pe =
            test_patch_matching(pw->type != PROF_ICC || pw->p.icc->lut == NULL ? &cam2xyzD50 : NULL,
                                lut, pw->type == PROF_DCP ? &dcp_lut : NULL,
                                pw->type == PROF_ICC ? pw->p.icc : NULL,
                                psD50, d50, &wb);
        elog("\n");
        print_patch_matching_errors(pe, d50, true, report_dir, "patch-errors", stderr);
        free(pe);
    }

    target_delete_patch_set(ps);
    target_delete_patch_set(psD50);
    free(dcp_lut.hsm);
}

static struct rgbimg *
spb_render(const struct spb_image *spb,
           const spectrum_t *illuminant,
           const struct observer *observer,
           const m3x3 *xyz2rgb,
           const enum ca_transform *cat,
           const double gamma,
           bool white_balance)
{
    // FIXME: we rotate and mirror here based on trial-error from spectral sample images.
    // Don't really know where the flip/rotate error is
    struct rgbimg *img = malloc(sizeof(*img) + 3*spb->w*spb->h*sizeof(img->p[0]));
    img->w = spb->h;
    img->h = spb->w;
    double bands[spb->band_count];
    for (uint32_t i = 0; i < spb->band_count; i++) {
        bands[i] = spb->band_start + spb->band_spacing * i;
    }
    v3 wb = {{1,1,1}};
    if (white_balance) {
        wb = spec2xyz(illuminant, observer);
        wb.v[0] = 1/wb.v[0];
        wb.v[1] = 1/wb.v[1];
        wb.v[2] = 1/wb.v[2];
    }
    v3 src_wp = spec2xyz(illuminant, observer);
    v3 dst_wp = icc_d50();
    if (xyz2rgb != NULL) {
        dst_wp = m3x3_mul_v3(m3x3_invert(*xyz2rgb), (v3){{1,1,1}});
    }
#pragma omp parallel for
    for (uint32_t y = 0; y < spb->h; y++) {
        for (uint32_t x = 0; x < spb->w; x++) {
            double amp[spb->band_count];
            for (uint32_t i = 0; i < spb->band_count; i++) {
                amp[i] = spb->p[spb->w * spb->h * i + (spb->w * y + x)];
            }
            spectrum_t *s = spec_alloc1(bands, amp, spb->band_count);
            v3 xyz = spec2xyz_ill(s, observer, illuminant, 1);
            if (cat != NULL) {
                xyz = chromatic_adaptation_transform(*cat, xyz, src_wp, dst_wp);
            }
            v3 rgb = (xyz2rgb == NULL) ? xyz : m3x3_mul_v3(*xyz2rgb, xyz);
            free(s);
            rgb.v[0] *= wb.v[0];
            rgb.v[1] *= wb.v[1];
            rgb.v[2] *= wb.v[2];
            rgb = v3_clip01(rgb);
            for (int i = 0; i < 3; i++) {
                img->p[3*(img->w*x+y)+i] = (uint16_t)roundf(65535 * pow(rgb.v[i], gamma));
            }
        }
    }
    return img;
}

static void
linearize_rgb(struct patch_set *ps,
              double *trc[3],
              int trc_count)
{
    if (trc_count == 0) {
        return;
    }
    for (int i = 0; i < ps->patch_count; i++) {
        v3 cam = ps->patch[i].cam;
        for (int j = 0; j < 3; j++) {
            double v = cam.v[j];
            if (v < 0) v = 0;
            if (v > 1) v = 1;
            if (trc_count == 1) {
                cam.v[j] = pow(cam.v[j], trc[j][0]);
            } else {
                int idx = floor(v * (trc_count-1));
                double d = v * (trc_count-1) - idx;
                if (d > 1e-7) {
                    cam.v[j] = (1.0 - d) * trc[j][idx] + d * trc[j][idx+1];
                } else {
                    cam.v[j] = trc[j][idx];
                }
            }
        }
        ps->patch[i].cam = cam;
    }
}

static void
delinearize_rgb(struct patch_set *ps,
                double *trc[3],
                int trc_count)
{
    if (trc_count == 0) {
        return;
    }
    if (trc_count == 1) {
        double g[3] = { 1.0/trc[0][0], 1.0/trc[1][0], 1.0/trc[2][0] }, *trc_[3] = { &g[0], &g[1], &g[2] };
        linearize_rgb(ps, trc_, 1);
        return;
    }
    const int out_len = 4096;
    double *out_trc[3];
    for (int j = 0; j < 3; j++) {
        out_trc[j] = malloc(out_len * sizeof(out_trc[j][0]));
#pragma omp parallel for
        for (int i = 0; i < out_len; i++) {
            out_trc[j][i] = inverse_curve((double)i / (out_len-1), trc[j], trc_count);
        }
    }
    /*
    FILE *stream2 = utf8_fopen("out", "w");
    FILE *stream1 = utf8_fopen("in", "w");
    for (int i = 0; i < trc_count; i++) {
        fprintf(stream1, "%f %f\n", (double)i / (trc_count-1), trc[0][i]);
    }
    for (int i = 0; i < out_len; i++) {
        fprintf(stream2, "%f %f\n", (double)i / (out_len-1), out_trc[0][i]);
    }
    fclose(stream1);
    fclose(stream2);
    */

    linearize_rgb(ps, out_trc, out_len);
    for (int i = 0; i < 3; i++) {
        free(out_trc[i]);
    }
}

static void
regenerate_patches(const char filename[],
                   struct patch_set *ps,
                   bool regenerate_xyz,
                   bool prefer_cat,
                   const struct dcam_ssf *dcam_ssf,
                   const struct observer *observer,
                   const spectrum_t *xyz_illuminant,
                   const v3 xyz_illuminant_wp,
                   enum exif_lightsource xyz_illuminant_ls,
                   const spectrum_t *cam_illuminant,
                   const v3 cam_illuminant_wp,
                   enum exif_lightsource cam_illuminant_ls)
{
    if (regenerate_xyz && dcam_ssf != NULL) {
        if (xyz_illuminant_ls != cam_illuminant_ls) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Re-generating target reference XYZ for illuminant %s and camera raw \n"
                "  RGB values for illuminant %s...\n",
                exif_lightsource_ntoa(xyz_illuminant_ls), exif_lightsource_ntoa(cam_illuminant_ls));
        } else {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Re-generating target reference XYZ and camera raw RGB values for\n"
                "  illuminant %s...\n", exif_lightsource_ntoa(xyz_illuminant_ls));
        }
    } else if (regenerate_xyz) {
        elog("Re-generating target reference XYZ values for illuminant %s...\n", exif_lightsource_ntoa(xyz_illuminant_ls));
    } else if (dcam_ssf) {
        elog("Re-generating target camera raw RGB values for illuminant %s...\n", exif_lightsource_ntoa(cam_illuminant_ls));
    } else {
        // do nothing
        return;
    }

    if (((cam_illuminant == NULL || xyz_illuminant == NULL) && v3_isequal(xyz_illuminant_wp, cam_illuminant_wp)) ||
        ((cam_illuminant != NULL && xyz_illuminant != NULL) && spec_equal(xyz_illuminant, cam_illuminant)))
    {
        prefer_cat = false;
    }
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].spectrum == NULL) {
            elog("Error: target \"%s\" lacks spetral data, cannot regenerate XYZ or RGB values.\n", filename);
            exit(EXIT_FAILURE);
        }
        if (dcam_ssf != NULL) {
            if (ps->patch[i].emissive) {
                ps->patch[i].cam = spec2xyz(ps->patch[i].spectrum, &dcam_ssf->ssf);
            } else {
                ps->patch[i].cam = spec2xyz_ill(ps->patch[i].spectrum, &dcam_ssf->ssf, cam_illuminant, 1);
            }
        }
        if (regenerate_xyz) {
            if (ps->patch[i].emissive) {
                ps->patch[i].xyz = spec2xyz(ps->patch[i].spectrum, observer);
            } else {
                if (prefer_cat) {
                    ps->patch[i].xyz = spec2xyz_ill(ps->patch[i].spectrum, observer, cam_illuminant, 1);
                    ps->patch[i].xyz = chromatic_adaptation_transform(CAT_CAT02, ps->patch[i].xyz, cam_illuminant_wp, xyz_illuminant_wp);
                } else {
                    ps->patch[i].xyz = spec2xyz_ill(ps->patch[i].spectrum, observer, xyz_illuminant, 1);
                }
            }
        }
    }
    if (prefer_cat) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "CAT02 was applied on XYZ values to transform illuminant \"%s\" to \"%s\".\n",
            exif_lightsource_ntoa(cam_illuminant_ls), exif_lightsource_ntoa(xyz_illuminant_ls));
    }
}

static struct patch_set *
regenerate_target_values(struct patch_set **ps_,
                         const char filename[],
                         enum exif_lightsource calib_ls,
                         const spectrum_t *calib_spec,
                         const v3 calib_wp,
                         const spectrum_t *xyz_ref_spec,
                         const v3 xyz_ref_wp,
                         const struct observer *observer,
                         const struct dcam_ssf *dcam_ssf,
                         bool prefer_cat,
                         bool regenerate_spectra)
{
    struct patch_set *ps = *ps_;

    // copy xyz to temp buffer if we should need the original xyz values later
    for (int i = 0; i < ps->patch_count; i++) {
        ps->patch[i].tmp = ps->patch[i].xyz;
    }

    bool relight = false;
    if (ps->patch[0].spectrum != NULL) {
        if (calib_spec == NULL) {
            if (dcam_ssf != NULL) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "Error: spectrum for calibration illuminant is not known, cannot regenerate\n"
                    "  camera raw RGB values. Either change calibration illuminant or don't provide\n"
                    "  camera SSFs (that is use the RGB values in the target file as is).\n");
                exit(EXIT_FAILURE);
            } else if (memcmp(&calib_wp, &xyz_ref_wp, sizeof(calib_wp)) != 0) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "Warning: the target XYZ reference values' illuminant is different from the\n"
                    "  calibration illuminant, and while the target has spectra to allow\n"
                    "  regeneration from scratch the spectrum for the calibration illuminant is\n"
                    "  not known. A relighting transform is required, which is less precise.\n");
                relight = true;
            }
        } else {
            regenerate_patches(filename, ps, true, false, dcam_ssf, observer, calib_spec, calib_wp, calib_ls, calib_spec, calib_wp, calib_ls);
        }
    } else if (memcmp(&calib_wp, &xyz_ref_wp, sizeof(calib_wp)) != 0) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Warning: the target XYZ reference values illuminant is different from the\n"
            "  calibration illuminant and the target lacks spectra so XYZ values cannot be\n"
            "  regenerated from scratch. A relighting transform is required, which will not\n"
            "  produce as precise results.\n");
        relight = true;
    }

    if (relight) {
        if (regenerate_spectra && xyz_ref_spec != NULL) {
            elog("Reconstructing virtual spectra to allow relighting...");
#pragma omp parallel for
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].spectrum = xyz2spec(ps->patch[i].xyz, observer, xyz_ref_spec);
                if (ps->patch[i].spectrum == NULL) {
                    elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                        "\n"
                        "Error: could not generate spectra for patch \"%s\". This is normal for\n"
                        "  ultra-saturated patches. Disable spectral reconstruction for this target\n"
                        "  and try again.\n", ps->patch[i].str_id);
                    exit(EXIT_FAILURE);
                }
                elog(".");
            }
            elog("\n");
        } else {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Using Bradford CAT as relighting transform, which generally is less accurate\n"
                "  than spectral reconstruction, consider using that instead.\n");
            for (int i = 0; i < ps->patch_count; i++) {
                ps->patch[i].xyz = chromatic_adaptation_transform(CAT_BRADFORD, ps->patch[i].xyz, xyz_ref_wp, calib_wp);
            }
        }
    }

    ps = exclude_weak_signal_patches(ps, true);

    struct patch_set *psD50 = target_copy(ps);

    if (calib_ls == lsD50) {
        // do nothing
    } else if (ps->patch[0].spectrum != NULL) {
        const spectrum_t *d50 = spectraldb_illuminant(lsD50);
        // we don't provide camera ssf here even if we have it, as we want the psD50 target to have the same cam values as the original
        const v3 d50_wp = v3_norm(spec2xyz(d50, observer));
        regenerate_patches(filename, psD50, true, prefer_cat, NULL, observer, d50, d50_wp, lsD50, calib_spec, calib_wp, calib_ls);
    } else if (prefer_cat) {
        for (int i = 0; i < psD50->patch_count; i++) {
            psD50->patch[i].xyz = chromatic_adaptation_transform(CAT_CAT02, ps->patch[i].xyz, calib_wp, glob.d50);
        }
        elog("CAT02 was applied on XYZ values to go from illuminant \"%s\" to \"%s\".\n",
                exif_lightsource_ntoa(calib_ls), exif_lightsource_ntoa(lsD50));
    } else {
        for (int i = 0; i < psD50->patch_count; i++) {
            psD50->patch[i].xyz = chromatic_adaptation_transform(CAT_BRADFORD, psD50->patch[i].tmp, xyz_ref_wp, glob.d50);
        }
    }
    *ps_ = ps;
    if (ps->patch_count == 0 || psD50->patch_count == 0) {
        elog("Error: no patches left.\n");
        exit(EXIT_FAILURE);
    }
    return psD50;
}

static v3
parse_arg_illuminant(const char arg[],
                     const struct observer *observer,
                     bool *used_observer,
                     spectrum_t **ill_spec)
{
    if (ill_spec != NULL) {
        *ill_spec = NULL;
    }
    v3 p_ill = {{ 0,0,0 }};
    if (strchr(arg, ',') != NULL) {
        // x,y,z syntax
        p_ill.v[0] = atof(arg);
        char *p = strchr(arg, ',')+1;
        p_ill.v[1] = atof(p);
        p = strchr(p, ',');
        if (p != NULL) {
            p_ill.v[2] = atof(p+1);
        } else {
            v3 xyY = {{p_ill.v[0], p_ill.v[1], 1.0}};
            p_ill = xyY2xyz(xyY);
        }
    } else {
        // spectrum
        spectrum_t *s = jsonio_illuminant_parse(arg, false);
        if (s != NULL) {
            p_ill = spec2xyz(s, observer);
            *used_observer = true;
            if (ill_spec != NULL) {
                *ill_spec = s;
            } else {
                free(s);
            }
        } else if ((int)exif_lightsource_aton(arg) >= 0) {
            elog("Error: illuminant \"%s\" is a known EXIF illuminant but there is no spectrum\n"
                 "  for it in the built-in database, so no XYZ coordinates can be found.\n"
                 "  Please choose a different illuminant, or provide illuminant spectrum\n"
                 "  as a file.\n", arg);
            exit(EXIT_FAILURE);
        }
    }
    p_ill = v3_norm(p_ill);
    if (p_ill.v[0] <= 0.0 || p_ill.v[1] <= 0.0 || p_ill.v[2] <= 0.0) {
        elog("Error: invalid illuminant \"%s\".\n", arg);
        exit(EXIT_FAILURE);
    }
    return p_ill;
}

static void
print_usage_and_exit(void)
{
    char *exif_ls = exif_lightsource_list("    ", 78);
    const char *usage_fmt =
        "DCamProf v" DCAMPROF_VERSION " - a digital camera profiling tool.\n"
        "\n"
        "Usage: dcamprof <command> [command-specific flags] <command args>\n"
        "\n"
        "Commands:\n"
        "\n"
        "  make-target [flags] <output.ti3>\n"
        "    flags:\n"
        "      -c <ssf.json> camera's spectral sensitivity functions\n"
        "      -o <observer> (default: 1931_2)\n"
        "      -i <target/calibration illuminant> RGB values illuminant (default: D50)\n"
        "      -I <reference illuminant> XYZ values illuminant (default: same as target illuminant)\n"
        "      -C don't model color inconstancy (use relighting instead of CAT)\n"
        "      -p <patches.ti3> include patch set\n"
        "      -a <name> assign (new) class name to previously included patch set (-p)\n"
        "      -f <file.tif | tf.json> linearize imported RGB values to match transfer function\n"
        "         in provided tiff/json\n"
        "      -S render spectra for inputs that lacks it\n"
        //        "      -s <low-limit> guess spectra below given nm limit down to 380nm\n"
        "      -g <generated grid spacing> (default: 0.03)\n"
        "      -d <distance> mimimum distance between patches of different classes (default: 0.02)\n"
        "      -b <distance> exclude patch if there is a lighter patch with same chromaticity\n"
        "         suggested chromaticity distance 0.004 (default: not active)\n"
        "      -x <exclude.txt> text file with sample id to exclude from target\n"
        "      -k <keep.txt> text file with sample id to (force) keep in the target\n"
        "      -X don't regenerate XYZ values of imported patch sets\n"
        "      -R don't regenerate RGB values of imported patch sets\n"
        "      -n exclude spectra in output (default: include if all inputs has it)\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "    built-in patch sets:\n"
        "      from measured data: %s\n"
        "      generated borders: %s\n"
        "      generated grids: <border>-grid\n"
        "\n"
        "  make-profile [flags] <target.ti3> <output.json | output.dcp | output.icc>\n"
        "    flags:\n"
        "      -n <camera name>\n"
        "      -w <patch class|name> <matrix opt weight>\n"
        "         Matrix optimizer patch class relative weights\n"
        "      -W skip patch normalization weighting\n"
        "      -v <patch class|name> <max DE | -DEL,+DEL,-DEC,+DEC,DEh|-DEh,+DEh>\n"
        "         Matrix optimizer refinement error ranges (default: none)\n"
        "      -V Apply refinement error ranges only to the forward matrix\n"
        "      -l <patch class|name | all> <max DE | -DEL,+DEL,-DEC,+DEC,DEh|-DEh,+DEh>\n"
        "         LUT optimizer maximum accepted error ranges, per class or name (default: auto)\n"
        "      -a <target-adjustment.json> apply target adjustment configuration\n"
        "      -y <Y | X,Y,Z> smallest allowed Y (or X,Y,Z) row value in forward matrix\n"
        "         optimization (default: -0.2 on Y only)\n"
        "      -k <LUT compress factor> (default: 0.7)\n"
        "      -d <distance> minimum patch chromaticity distance in LUT (default: 0.02)\n"
        "      -g <target-layout.json> provide target layout for glare matching and/or\n"
        "         flatfield correction\n"
        "      -o <observer> (default: 1931_2)\n"
        "      -c <ssf.json> camera SSFs, used to re-generate target RGB values\n"
        "      -i <calibration illuminant> (default: D50)\n"
        "      -I <target XYZ reference values illuminant> (default: same as calibration illuminant)\n"
        "      -C don't model color inconstancy (use relighting instead of CAT)\n"
        "      -S render spectra if target lacks it\n"
        "      -B don't re-balance target so most neutral patch becomes 100%% neutral\n"
        "      -b <patch name or index> manually point out most neutral patch in target\n"
        "      -x <exclude.txt> text file with sample id to exclude from target\n"
        "      -p <matrix.json | dcp.json | profile.json> provide pre-calculated color matrix\n"
        "      -f <matrix.json | dcp.json | profile.json> provide pre-calculated forward matrix\n"
        "      -e <matrix.json | dcp.json | profile.json> provide pre-calculated LUT matrix\n"
        "      -m <profile.json> import all matrices\n"
        "      -s use (very) slow but more robust and possibly more exact matrix optimizer\n"
        "      -t <linear | none | acr | custom.json> embed a tone-curve in output DCP or ICC.\n"
        "      -L skip LUT in informational report (LUT is always generated anyway)\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  test-profile [flags] [target.ti3 | test.tif] <profile.json|.dcp|.icc> [out.tif]\n"
        "    flags:\n"
        "      -o <observer> (default: 1931_2)\n"
        "      -c <ssf.json> camera SSFs, used to re-generate target RGB values\n"
        "      -i <test illuminant> (default: same as profile's calibration illuminant)\n"
        "      -I <target XYZ reference values illuminant> (default: same as test illuminant)\n"
        "      -C don't model color inconstancy (use relighting instead of CAT)\n"
        "      -S render spectra if target lacks it\n"
        "      -B don't re-balance target so most neutral patch becomes 100%% neutral\n"
        "      -b <patch name or index> manually point out most neutral patch in target\n"
        "      -w <r,g,b> | m<r,g,b> provide camera WB (default: derived from target)\n"
        "      -L skip LUT(s)\n"
        "      -P skip LookTable (DCP only)\n"
        "      -T don't add default ACR curve if DCP lacks curve\n"
        "      -f <file.tif | tf.json> de-linearize RGB values in target, that is run\n"
        "         provided transfer function backwards\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  make-dcp [flags] <profile.json> [profile2.json] <output.dcp>\n"
        "    flags:\n"
        "      -n <unique camera name>\n"
        "      -d <profile name>\n"
        "      -c <copyright>\n"
        "      -s <calibration signature>\n"
        "      -b <baseline exposure offset>\n"
        "      -B don't include default black render = none tag\n"
        "      -i <calibration illuminant1> (default: from source profile)\n"
        "      -I <calibration illuminant2> (default: from source profile)\n"
        "      -m <other.dcp> copy illuminant(s) and color matrices from the provided DCP\n"
        "      -h <hdiv,sdiv,vdiv> hue, saturation and value divisions of LUTs (default: 90,30,30)\n"
        "      -v <max curve matching error> (default: 0.0019)\n"
        "      -F skip forward matrix (use only color matrix)\n"
        "      -E don't use LUT matrix as forward matrix for a LUT profile\n"
        "      -L skip LUT (= matrix-only profile)\n"
        "      -O disable forward matrix whitepoint remapping\n"
        "      -G skip gamma-encoding of 3D LUTs\n"
        "      -D make the HueSatMap LUT 3D instead of 2.5D\n"
        "      -H allow discontinuity in LUT hue shifts\n"
        "      -t <linear | none | acr | custom.json> embed a tone-curve (default: linear)\n"
        "      -o <neutral | standard | custom.json> tone reproduction operator (default: neutral)\n"
        "      -g <none | srgb | adobergb | srgb-strong | adobergb-strong> gamut compression\n"
        "         (default: none / from TRO file)\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  make-icc [flags] <profile.json> <output.icc>\n"
        "    flags:\n"
        "      -n <camera name> (profile description)\n"
        "      -c <copyright>\n"
        "      -s <clut divisions> divisions of side of LUT cube (default: 33)\n"
        "      -p <lablut | xyzlut | matrix> profile type\n"
        "         (default: lablut if input has LUT otherwise matrix)\n"
        "      -L skip the input profile's LUT\n"
        "      -W apply white balance correction\n"
        "      -f <file.tif | tf.json> adapt ICC to match transfer function in provided tiff/json\n"
        "      -t <none | acr | custom.json> embed a tone-curve (default: none)\n"
        "      -o <neutral | standard | custom.json> tone reproduction operator (default: neutral)\n"
        "      -g <none | srgb | adobergb | srgb-strong | adobergb-strong> gamut compression\n"
        "         (default: none / from TRO file)\n"
        "      -T don't apply tone curve to LUT\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  dcp2json <input.dcp> [<output.json>]\n"
        "  json2dcp <input.json> <output.dcp>\n"
        "  icc2json <input.icc> [<output.json>]\n"
        "  json2icc <input.json> <output.icc>\n"
        "  average-targets <input1.ti3> [<input2.ti3> ...] <output.ti3>\n"
        "\n"
        "  tiff-tf <input.tif> [<output.json>]\n"
        "    flags:\n"
        "      -R skip reconstruction of truncated values\n"
        "      -f <linear.tif> reference tiff with transfer function\n"
        "\n"
        "  txt2ti3 [flags] <input.txt> <output.ti3>\n"
        "    flags:\n"
        "      -f <low,high,spacing> (default: 400,700,5)\n"
        "      -s <scale> (default: 1.0)\n"
        "      -m <clip level after scaling> (default: unbounded)\n"
        "      -l <rows|cols> row or column layout (default: cols)\n"
        "      -a <name> assign class name (default: \"import\")\n"
        "\n"
        "  make-testchart [flags] <output.ti1>\n"
        "    flags:\n"
        "      -p <patch count> (default: 100)\n"
        "      -w <percentage white patches> (default: 20)\n"
        "      -b <black patches> (default: 5)\n"
        "      -g <inbetween gray steps> (default: 0)\n"
        "      -l <layout row count> (default: unspecified)\n"
        "      -d <layout row relative height>,<column relative width> (default: 1,1)\n"
        "      -O layout has even columns offset a half patch\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  testchart-ff <input.tif> <flatfield.tif> <corrected.tif>\n"
        "  testchart-ff [flags] <input.ti1 | layout.json> <input.ti3> [input2.ti3] <output.ti3>\n"
        "    flags:\n"
        "      -l, -d, -O flags from make-testchart to specify chart layout.\n"
        "      -L enable glare matching\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  match-spectra [flags] <reference.ti3> <match.ti3> <output-match.ti3> [<output-ref.ti3>]\n"
        "    flags:\n"
        "      -o <observer> observer for DE comparison (default: 1931_2)\n"
        "      -i <illuminant> illuminant for DE comparison (default: D50)\n"
        "      -c <ssf.json> camera SSFs, used instead of observer for DE comparison\n"
        "      -S scale spectra in output for best match\n"
        "      -N normalize before comparison\n"
        "      -U don't allow repeats of same spectrum in output\n"
        "      -E consider all spectra as emissive\n"
        "      -e <max DE> maximum acceptable delta E for match (default: infinite)\n"
        "      -r <dir> directory to save informational reports and plots\n"
        "\n"
        "  si-render [flags] <spectral_image.spb> <output.tif>"
        "    flags:\n"
        "      -i <illuminant> (default: D50)\n"
        "      -o <observer> (default: 1931_2)\n"
        "      -c <ssf.json> camera SSFs, use instead of observer\n"
        "      -g <gamma> gamma in output (default: 1.0)\n"
        "      -W apply white balance\n"
        "      -b <base band>[,<band width>] specify base band and width for indexed files\n"
        "      -P enable ProPhotoRGB output\n"
        "      -a <bradford | cat02> chromatic adaptation transform (default: cat02 when needed)\n"
        "\n"
        "Illuminants:\n"
        "  <x,y> whitepoint coordinate in xyY color space\n"
        "  <X,Y,Z> whitepoint coordinate in XYZ color space\n"
        "  <exif lightsource> integer or string.\n"
        "  exif lightsource names and their temperatures:\n"
        "%s\n"
        "  <spectrum.json | spectrum.sp> (JSON or Argyll spectrum file)\n"
        "  built-in illuminant spectra: %s, and <temp>K for blackbodies\n"
        "\n"
        "Observers:\n"
        "  built-in observers: %s\n"
        "\n"
        "Program version: " DCAMPROF_VERSION "\n"
        ;
    elog(usage_fmt, target_measured_list(), target_generated_list(),
            exif_ls,
            spectraldb_illuminant_list(), observer_list());
    free(exif_ls);
    exit(EXIT_FAILURE);
}

static void
init_observer_dependent_globals(const struct observer *obs)
{
    gamut_init(obs, spectraldb_illuminant(lsD50));
    glob.d50 = spec2xyz(spectraldb_illuminant(lsD50), obs);
    glob.d50 = v3_norm(glob.d50);
}

static bool
has_suffix(const char s[],
           const char suffix[])
{
    const char *p = strstr(s, suffix);
    if (!p) return false;
    return p[strlen(suffix)] == '\0';
}

struct compare_tristimulus_values_simplex_fun_arg {
    v3 xyz_ref;
    v3 xyz;
    v3 wp;
    bool use_ciede2000;
};

static double
compare_tristimulus_values_simplex_fun(double m[],
                                     void *arg)
{
    if (m[0] < 0) {
        return 100000000;
    }
    struct compare_tristimulus_values_simplex_fun_arg *a = (struct compare_tristimulus_values_simplex_fun_arg *)arg;
    v3 test_xyz = k_mul_v3(m[0], a->xyz);
    double de = a->use_ciede2000 ? ciede2000(test_xyz, a->xyz_ref, a->wp) : euclidean_dist(test_xyz, a->xyz_ref);
    return de;
}

static double
compare_tristimulus_values(v3 xyz,
                           v3 xyz_ref,
                           v3 wp,
                           bool normalize,
                           bool use_ciede2000,
                           double *scale_factor)
{
    if (normalize) {
        xyz_ref = v3_norm(xyz_ref);
        xyz = v3_norm(xyz);
    }
    struct compare_tristimulus_values_simplex_fun_arg arg;
    arg.xyz_ref = xyz_ref;
    arg.xyz = xyz;
    arg.wp = wp;
    arg.use_ciede2000 = use_ciede2000;
    double m = 1.0;
    double de = simplex(compare_tristimulus_values_simplex_fun, &arg, &m, 1, 1.0e-7, 0.01, NULL);
    *scale_factor = m;
    if (normalize) {
        return de;
    }
    return use_ciede2000 ? ciede2000(xyz, xyz_ref, wp) : euclidean_dist(xyz, xyz_ref);
}

static v3
apply_gamma(v3 rgb,
            bool is_srgb_gamma,
            double gamma,
            bool is_forward)
{
    for (int i = 0; i < 3; i++) {
        if (is_forward) {
            rgb.v[i] = (is_srgb_gamma) ? srgb_gamma_forward(rgb.v[i]) : pow(rgb.v[i], 1.0/gamma);
        } else {
            rgb.v[i] = (is_srgb_gamma) ? srgb_gamma_inverse(rgb.v[i]) : pow(rgb.v[i], gamma);
        }
    }
    return rgb;
}

static void
parse_tonecurve(const char arg[],
                double **tc_,
                int *tc_len_)
{
    double *old_tc = *tc_;
    int old_tc_len = *tc_len_;
    double *tc;
    int tc_len = 0;
    if (strcasecmp(arg, "linear") == 0) {
        tc_len = 2;
        tc = malloc(2 * sizeof(tc[0]));
        tc[0] = 0.0;
        tc[1] = 1.0;
    } else if (strcasecmp(arg, "acr") == 0) {
        tc = dcp_acr_tonecurve(&tc_len);
    } else if (strcasecmp(arg, "none") == 0) {
        // do nothing
        return;
    } else {
        jsonio_tonecurve_parse(arg, &tc, &tc_len, true);
        if (tc_len == 1) {
            double gamma = tc[0];
            free(tc);
            tc = malloc(4096 * sizeof(tc[0]));
            for (int i = 0; i < 4096; i++) {
                tc[i] = pow(i / 4095.0, gamma);
            }
            tc_len = 4096;
        }
    }
    reconstruct_tonecurve(tc, tc_len);
    if (old_tc == NULL) {
        *tc_ = tc;
        *tc_len_ = tc_len;
        return;
    }
    elog("Additional tone-curve specified, merging.\n");
    const int mtc_len = 65536;
    double *mtc = malloc(mtc_len * sizeof(mtc[0]));
    for (int i = 0; i < mtc_len; i++) {
        double x = (double)i / (mtc_len - 1);
        x = curve_y(x, old_tc, old_tc_len);
        mtc[i] = curve_y(x, tc, tc_len);
    }
    *tc_ = mtc;
    *tc_len_ = mtc_len;
    free(old_tc);
}

static enum gc_type
parse_gc_type(const char arg[])
{
    enum gc_type gc_type;
    if (strcasecmp(arg, "none") == 0) {
        gc_type = GC_NONE;
    } else if (strcasecmp(arg, "srgb") == 0) {
        gc_type = GC_SRGB;
    } else if (strcasecmp(arg, "adobergb") == 0) {
        gc_type = GC_ADOBERGB;
    } else if (strcasecmp(arg, "srgb-strong") == 0) {
        gc_type = GC_SRGB_STRONG;
    } else if (strcasecmp(arg, "adobergb-strong") == 0) {
        gc_type = GC_ADOBERGB_STRONG;
    } else {
        elog("Error: unknown gamut compression type \"%s\".\n", arg);
        exit(EXIT_FAILURE);
    }
    return gc_type;
}

static void
test_gc_type_setting(enum gc_type gc_type,
                     enum tc_type tc_type,
                     double *tc,
                     bool skip_lut)
{
    if (gc_type == GC_FROM_FILE || gc_type == GC_NONE) {
        return;
    }
    if (tc_type != TC_NEUTRAL) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Error: gamut compression will not be applied as it requires the neutral tone\n"
            "  reproduction operator ('-o neutral').\n");
        exit(EXIT_FAILURE);
    }
    if (tc == NULL) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Error: gamut compression will not be applied as it requires a tone curve\n"
            "  ('-t linear' is ok).\n");
        exit(EXIT_FAILURE);
    }
    if (skip_lut) {
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Error: gamut compression will not be applied as it requires a LUT (skip '-L').\n");
        exit(EXIT_FAILURE);
    }
}

static struct test_chart_spec *
ti1_to_tcs(const struct test_chart_layout *layout,
           const char filename[])
{
    struct patch_set *ps = argyllio_ti3_parse(filename, true, true);
    struct test_chart_spec *tcs = calloc(1, sizeof(*tcs) + ps->patch_count * sizeof(tcs->data[0]));
    tcs->patch_count = ps->patch_count;
    tcs->layout = *layout;
    tcs->col_count = tcs->patch_count / tcs->layout.row_count;
    if (tcs->col_count * tcs->layout.row_count != tcs->patch_count) {
        elog("Error: row count %d not divisable with .ti1 patch count %d.\n", tcs->layout.row_count, tcs->patch_count);
        exit(EXIT_FAILURE);
    }
    double max_ti1_g = 0;
    double min_ti1_g = 1e9;
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].cam.v[0] != ps->patch[i].cam.v[1] || ps->patch[i].cam.v[0] != ps->patch[i].cam.v[2]) {
            // not neutral
            continue;
        }
        if (ps->patch[i].cam.v[1] > max_ti1_g) {
            max_ti1_g = ps->patch[i].cam.v[1];
        }
        if (ps->patch[i].cam.v[1] < min_ti1_g) {
            min_ti1_g = ps->patch[i].cam.v[1];
        }
    }
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].cam.v[0] != ps->patch[i].cam.v[1] || ps->patch[i].cam.v[0] != ps->patch[i].cam.v[2]) {
            continue;
        }
        if (ps->patch[i].cam.v[1] == max_ti1_g) {
            tcs->white_count++;
        } else if (max_ti1_g != min_ti1_g && ps->patch[i].cam.v[1] == min_ti1_g) {
            tcs->black_count++;
        } else {
            tcs->gray_count++;
        }
    }
    tcs->black = (tcs->black_count > 0) ? &tcs->data[0] : NULL;
    tcs->white = (tcs->white_count > 0) ? &tcs->data[tcs->black_count] : NULL;
    tcs->gray = (tcs->gray_count > 0) ? &tcs->data[tcs->black_count+tcs->white_count] : NULL;
    tcs->white_count = tcs->black_count = tcs->gray_count = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        if (ps->patch[i].cam.v[0] != ps->patch[i].cam.v[1] || ps->patch[i].cam.v[0] != ps->patch[i].cam.v[2]) {
            continue;
        }
        struct test_chart_patch *p;
        if (ps->patch[i].cam.v[1] == max_ti1_g) {
            p = &tcs->white[tcs->white_count++];
            p->is_white = true;
        } else if (max_ti1_g != min_ti1_g && ps->patch[i].cam.v[1] == min_ti1_g) {
            p = &tcs->black[tcs->black_count++];
            p->is_black = true;
        } else {
            p = &tcs->gray[tcs->gray_count++];
            p->is_gray = true;
        }
        strcpy(p->str_id, ps->patch[i].str_id);
        p->idx = i;
        test_chart_rc(&tcs->layout, i, &p->r, &p->c);
    }
    free(ps);
    return tcs;
}

static void
match_tcs_with_ps(struct test_chart_spec *tcs,
                  const struct patch_set *ps)
{
    for (int i = 0; i < ps->patch_count; i++) {
        bool found = false;
        for (int j = 0; j < tcs->white_count + tcs->black_count + tcs->gray_count; j++) {
            if (tcs->data[j].str_id[0] == '\0') {
                elog("Error: empty patch name in target layout.\n");
                exit(EXIT_FAILURE);
            }
            if (strcasecmp(tcs->data[j].str_id, ps->patch[i].str_id) == 0) {
                if (found) {
                    elog("Error: duplicate patch name \"%s\" in target layout JSON file.\n", tcs->data[j].str_id);
                    exit(EXIT_FAILURE);
                }
                if (tcs->data[j].idx != -1) {
                    elog("Error: duplicate patch name \"%s\" in TI3 file.\n", tcs->data[j].str_id);
                    exit(EXIT_FAILURE);
                }
                tcs->data[j].idx = i;
                test_chart_rc(&tcs->layout, i, &tcs->data[j].r, &tcs->data[j].c);
                found = true;
            }
        }
    }
    for (int j = 0; j < tcs->white_count + tcs->black_count + tcs->gray_count; j++) {
        if (tcs->data[j].idx == -1) {
            elog("Error: patch \"%s\" in target layout JSON not found in TI3.\n", tcs->data[j].str_id);
            exit(EXIT_FAILURE);
        }
        for (int k = 0; k < tcs->white_count + tcs->black_count + tcs->gray_count; k++) {
            if (k != j && tcs->data[k].idx == tcs->data[j].idx) {
                elog("Error: duplicate patch index in target layout JSON for \"%s\".\n", tcs->data[j].str_id);
                exit(EXIT_FAILURE);
            }
        }
    }
}

int
main(int argc,
     char *argv[])
{
    if (argc < 3) { // minimum three: dcamprof <command> <arg>
        if (argc == 2 && strcmp(argv[1], "-v") == 0) {
            elog("Program version: " DCAMPROF_VERSION "\n");
            exit(EXIT_SUCCESS);
        }
        print_usage_and_exit();
    }

    bool parse_flags = true;
    int ai = 2;
    if (strcmp(argv[1], "dcp2json") == 0) {

        /******************
          DCP to JSON
         */

        if (ai != argc - 2 && ai != argc - 1) {
            print_usage_and_exit();
        }
        struct profio_dcp *dcp = profio_dcp_parse(argv[ai+0], true);
        FILE *out = stdout;
        if (ai == argc - 2) {
            out = utf8_fopen(argv[ai+1], "w");
            if (out == NULL) {
                elog("Error: could not open \"%s\" for writing: %s.\n", argv[ai+1], strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
        profio_dcp_json(dcp, out);
        if (out != stdout) {
            fclose(out);
        }
        profio_dcp_delete(dcp);

    } else if (strcmp(argv[1], "icc2json") == 0) {

        /******************
          ICC to JSON
         */
        if (ai != argc - 2 && ai != argc - 1) {
            print_usage_and_exit();
        }
        struct profio_icc *icc = profio_icc_parse(argv[ai+0], true);
        FILE *out = stdout;
        if (ai == argc - 2) {
            out = utf8_fopen(argv[ai+1], "w");
            if (out == NULL) {
                elog("Error: could not open \"%s\" for writing: %s.\n", argv[ai+1], strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
        profio_icc_json(icc, out);
        if (out != stdout) {
            fclose(out);
        }
        profio_icc_delete(icc);

    } else if (strcmp(argv[1], "json2dcp") == 0) {

        /******************
          JSON to DCP
         */
        if (ai != argc - 2) print_usage_and_exit();
        struct profio_dcp *dcp = jsonio_dcp_parse(argv[ai+0], NULL, true);
        profio_dcp_write(dcp, argv[ai+1], true);
        profio_dcp_delete(dcp);

    } else if (strcmp(argv[1], "json2icc") == 0) {

        /******************
          JSON to ICC
         */
        if (ai != argc - 2) print_usage_and_exit();
        struct profio_icc *icc = jsonio_icc_parse(argv[ai+0], NULL, true);
        profio_icc_write(icc, argv[ai+1], true);
        profio_icc_delete(icc);

    } else if (strcmp(argv[1], "make-target") == 0) {

        /******************
          Make target
         */
        const char *ps_name[argc];
        const char *ps_class_name[argc];
        memset(ps_class_name, 0, sizeof(ps_class_name));
        int ps_count = 0;
        const char *report_dir = NULL;
        struct dcam_ssf *dcam_ssf = NULL;
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        spectrum_t *illuminant = spec_copy(spectraldb_illuminant(lsD50));
        enum exif_lightsource illuminant_ls = lsD50;
        enum exif_lightsource ref_illuminant_ls = lsD50;
        spectrum_t *ref_illuminant = NULL;
        double min_chroma_dist = 0.02;
        double min_chroma_dist_lighter = -1;
        const char *exclude_list = NULL;
        const char *keep_list = NULL;
        bool exclude_spectra = false;
        bool regenerate_xyz = true;
        bool regenerate_cam = true;
        bool prefer_cat = true;
        bool used_observer = false;
        bool render_spectra = false;
        double grid_spacing = 0.03;
        //        double guess_below_nm = 0;
        double *trc[3] = { NULL, NULL, NULL };
        int trc_count = 0;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'c': dcam_ssf = jsonio_dcam_ssf_parse(argv[++ai], true); break;
            case 'i':
                free(illuminant);
                illuminant = jsonio_illuminant_parse(argv[++ai], true);
                if ((int)exif_lightsource_aton(argv[ai]) >= 0) {
                    illuminant_ls = exif_lightsource_aton(argv[ai]);
                } else {
                    illuminant_ls = lsOther;
                }
                break;
            case 'I':
                free(ref_illuminant);
                ref_illuminant = jsonio_illuminant_parse(argv[++ai], true);
                if ((int)exif_lightsource_aton(argv[ai]) >= 0) {
                    ref_illuminant_ls = exif_lightsource_aton(argv[ai]);
                } else {
                    ref_illuminant_ls = lsOther;
                }
                break;
            case 'C': prefer_cat = false; break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminants on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'p': ps_name[ps_count++] = argv[++ai]; break;
            case 'a': {
                if (ps_count == 0) {
                    elog("Error: cannot assign class name, no patch set included.\n");
                    exit(EXIT_FAILURE);
                }
                ai++;
                if (strcasecmp(argv[ai], "all") == 0 ||
                    strcasecmp(argv[ai], "white") == 0 ||
                    strcasecmp(argv[ai], "illuminant") == 0) {
                    elog("Error: \"%s\" is a reserved class name.\n", argv[ai]);
                    exit(EXIT_FAILURE);
                }
                ps_class_name[ps_count-1] = argv[ai];
                break;
            }
            case 'b': min_chroma_dist_lighter = atof(argv[++ai]); break;
            case 'd': min_chroma_dist = atof(argv[++ai]); break;
            case 'g': grid_spacing = atof(argv[++ai]); break;
            case 'x': exclude_list = argv[++ai]; break;
            case 'k': keep_list = argv[++ai]; break;
            case 'X': regenerate_xyz = false; break;
            case 'R': regenerate_cam = false; break;
            case 'S': render_spectra = true; break;
            case 'n': exclude_spectra = true; break;
                //case 's': guess_below_nm = atof(argv[++ai]); break;
            case 'f': jsonio_transfer_function_parse(argv[++ai], trc, &trc_count, true); break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer);
        if (ai != argc - 1) print_usage_and_exit();
        if (ps_count == 0) {
            elog("No patch sets included (-p), nothing to do.\n");
            exit(EXIT_FAILURE);
        }
        if (ref_illuminant == NULL) {
            ref_illuminant = spec_copy(illuminant);
            ref_illuminant_ls = illuminant_ls;
        }
        const char *load_ps[ps_count];
        const char *load_ps_names[ps_count];
        int load_ps_idx[ps_count];
        int load_ps_count = 0;
        const char *meas_ps[ps_count];
        const char *meas_ps_names[ps_count];
        int meas_ps_idx[ps_count];
        int meas_ps_count = 0;
        const char *gene_ps[ps_count];
        const char *gene_ps_names[ps_count];
        int gene_ps_idx[ps_count];
        int gene_ps_count = 0;
        const char *grid_ps[ps_count];
        const char *grid_ps_names[ps_count];
        int grid_ps_idx[ps_count];
        int grid_ps_count = 0;
        for (int i = 0; i < ps_count; i++) {
            bool is_measured;
            if (target_isvalid_target_name(ps_name[i], &is_measured)) {
                if (is_measured) {
                    meas_ps_names[meas_ps_count] = ps_class_name[i];
                    meas_ps_idx[meas_ps_count] = i;
                    meas_ps[meas_ps_count++] = ps_name[i];
                } else {
                    gene_ps_names[gene_ps_count] = ps_class_name[i];
                    gene_ps_idx[gene_ps_count] = i;
                    gene_ps[gene_ps_count++] = ps_name[i];
                }
            } else {
                bool found = false;
                const char *p = strstr(ps_name[i], "-grid");
                if (p == NULL) {
                    p = strstr(ps_name[i], "-GRID");
                }
                if (p != NULL && p[5] == '\0') {
                    char name[strlen(ps_name[i])];
                    strcpy(name, ps_name[i]);
                    name[strlen(name)-5] = '\0';
                    if (target_isvalid_target_name(name, &is_measured)) {
                        if (is_measured) {
                            elog("Error: cannot generate grid on built-in measurement based patch set \"%s\".\n", ps_name[i]);
                            exit(EXIT_FAILURE);
                        }
                        found = true;
                        grid_ps_names[grid_ps_count] = ps_class_name[i];
                        grid_ps_idx[grid_ps_count] = i;
                        grid_ps[grid_ps_count++] = ps_name[i];
                    }
                }
                if (!found) {
                    load_ps_names[load_ps_count] = ps_class_name[i];
                    load_ps_idx[load_ps_count] = i;
                    load_ps[load_ps_count++] = ps_name[i];
                }
            }
        }
        print_common_spectra(report_dir, observer, illuminant, dcam_ssf);
        struct patch_set *inps[ps_count];
        for (int i = 0; i < load_ps_count; i++) {
            elog("Loading patch set from \"%s\"...\n", load_ps[i]);
            struct patch_set *ps = argyllio_ti3_parse(load_ps[i], true, true);
            if (render_spectra) {
                bool need_spectra = false;
                for (int j = 0; j < ps->patch_count; j++) {
                    if (ps->patch[j].spectrum == NULL) {
                        need_spectra = true;
                        break;
                    }
                }
                if (need_spectra) {
                    elog("Patch set lacks spectra, rendering new...");
#pragma omp parallel for
                    for (int j = 0; j < ps->patch_count; j++) {
                        if (ps->patch[j].spectrum == NULL) {
                            ps->patch[j].spectrum = xyz2spec(ps->patch[j].xyz, observer, illuminant);
                            if (ps->patch[j].spectrum == NULL) {
                                elog("Error: failed to render spectra for patch %s.\n", ps->patch[j].str_id);
                                exit(EXIT_FAILURE);
                            }
                        }
                        elog(".");
                    }
                    elog("\n");
                }
            }
            if (!exclude_spectra) {
                for (int j = 0; j < ps->patch_count; j++) {
                    if (ps->patch[j].spectrum == NULL) {
                        if (regenerate_xyz) {
                            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                                "Error: patch set \"%s\" lacks spectral information so cannot regenerate XYZ\n"
                                "  values. You must use the -X flag to process this target.\n", load_ps[i]);
                            exit(EXIT_FAILURE);
                        }
                        if (regenerate_cam && dcam_ssf != NULL) {
                            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                                "Error: patch set \"%s\" lacks spectral information so cannot regenerate RGB\n"
                                "  values. You must use the -R flag to process this target.\n", load_ps[i]);
                            exit(EXIT_FAILURE);
                        }
                        elog("Patch set \"%s\" lacks spectral information; disabling spectra in output.\n", load_ps[i]);
                        exclude_spectra = true;
                        break;
                    }/* else if (guess_below_nm > 380) {
                        spectrum_t *s = ps->patch[j].spectrum;
                        ps->patch[j].spectrum = spec_guess_lows(s, guess_below_nm);
                        free(s);
                        }*/
                }
            }
            if (trc_count > 0) {
                elog("Linearizing RGB values using provided transfer function.\n");
                linearize_rgb(ps, trc, trc_count);
            }
            if (load_ps_names[i] != NULL) {
                target_set_name(ps, load_ps_names[i]);
            } else {
                target_set_name(ps, "import");
            }
            regenerate_patches(load_ps[i], ps, regenerate_xyz, prefer_cat, regenerate_cam ? dcam_ssf : NULL, observer, ref_illuminant, v3_norm(spec2xyz(ref_illuminant, observer)), ref_illuminant_ls, illuminant, v3_norm(spec2xyz(illuminant, observer)), illuminant_ls);
            inps[load_ps_idx[i]] = ps;
        }
        for (int i = 0; i < meas_ps_count; i++) {
            elog("Loading built-in measurement-based patch set \"%s\"...\n", meas_ps[i]);
            inps[meas_ps_idx[i]] = target_make_virtual(meas_ps[i], dcam_ssf, observer, illuminant, ref_illuminant, prefer_cat, 0.0, false, false, NULL, report_dir);
            if (meas_ps_names[i] != NULL) {
                target_set_name(inps[meas_ps_idx[i]], meas_ps_names[i]);
            }
        }
        for (int i = 0; i < gene_ps_count; i++) {
            elog("Generating chromaticity border patch set \"%s\"...\n", gene_ps[i]);
            inps[gene_ps_idx[i]] = target_make_virtual(gene_ps[i], dcam_ssf, observer, illuminant, ref_illuminant, prefer_cat, grid_spacing, false, false, NULL, report_dir);
            if (gene_ps_names[i] != NULL) {
                target_set_name(inps[gene_ps_idx[i]], gene_ps_names[i]);
            }
        }
        for (int i = 0; i < grid_ps_count; i++) {
            elog("Generating chromaticity grid patch set \"%s\"...\n", grid_ps[i]);
            char name[strlen(grid_ps[i])+1];
            strcpy(name, grid_ps[i]);
            name[strlen(name)-5] = '\0';
            inps[grid_ps_idx[i]] = target_make_virtual(name, dcam_ssf, observer, illuminant, ref_illuminant, prefer_cat, grid_spacing, true, false, NULL, report_dir);
            if (grid_ps_names[i] != NULL) {
                target_set_name(inps[grid_ps_idx[i]], grid_ps_names[i]);
            }
        }
        if (keep_list != NULL) {
            int keep_count = 0;
            for (int i = 0; i < ps_count; i++) {
                keep_count += textio_process_keep_list(inps[i], keep_list);
            }
            if (keep_count > 0) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "%d patches matched the provided keep list and will therefore not be excluded.\n", keep_count);
            } else {
                elog("No patches matched the provided keep list.\n");
            }
        }
        struct patch_set *ps;
        if (ps_count > 1) {
            int tot_count = 0;
            for (int i = 0; i < ps_count; i++) {
                tot_count += inps[i]->patch_count;
            }
            elog("Merging %d patch sets with %d patches in total...\n", ps_count, tot_count);
            ps = target_merge(inps, ps_count, min_chroma_dist, true);
            if (ps->patch_count != tot_count) {
                int class_count = 0;
                bool has_class[TARGET_MAX_CLASS_COUNT];
                memset(has_class, 0, sizeof(has_class));
                for (int i = 0; i < ps->patch_count; i++) {
                    if (!has_class[ps->patch[i].class]) {
                        has_class[ps->patch[i].class] = true;
                        class_count++;
                    }
                }
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "%d patches were excluded due to nearby neighbors of different class.\n"
                    "  %d patches in %d class%s left in merged set.\n",
                    tot_count - ps->patch_count, ps->patch_count, class_count, class_count == 1 ? "" : "es");
            }
        } else {
            ps = inps[0];
        }
        if (min_chroma_dist_lighter >= 0) {
            ps = exclude_darker_patches(ps, min_chroma_dist_lighter);
        }
        if (exclude_list != NULL) {
            int exclude_count = textio_process_exclude_list(ps, exclude_list, NULL);
            if (exclude_count > 0) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "%d patches matched provided exclude list and was therefore excluded.\n"
                    "  %d patches left\n",
                    exclude_count, ps->patch_count);
            } else {
                elog("No patches matched the provided exclude list, all patches kept.\n");
            }
        }
        if (ps->patch_count == 0) {
            elog("Error: no patches left in target.\n");
            exit(EXIT_FAILURE);
        }
        const v3 wp = v3_norm(spec2xyz(illuminant, observer));
        elog("Writing output to \"%s\"...\n", argv[ai]);
        argyllio_ti3_write(argv[ai], !exclude_spectra, ps, true);
        target_print("target", ps, wp, observer, report_dir);
        print_transfer_function(report_dir, trc, trc_count);
        target_delete_patch_set(ps);
        jsonio_dcam_ssf_delete(dcam_ssf);
        free(illuminant);
        free(ref_illuminant);
        free(trc[0]);
        free(trc[1]);
        free(trc[2]);
        elog("Complete!\n");

    } else if (strcmp(argv[1], "make-profile") == 0) {

        /******************
          Make profile
         */
        const char *camera_name = NULL;
        const char *report_dir = NULL;
        enum exif_lightsource calib_ls = lsD50;
        spectrum_t *calib_spec = spec_copy(spectraldb_illuminant(lsD50));
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        bool skip_lut_in_report = false;
        bool used_observer = false;
        bool use_slow_matrix_optimizer = false;
        bool auto_relax = true;
        v3 calib_wp = {{0,0,0}};
        v3 xyz_ref_wp = {{0,0,0}};
        spectrum_t *xyz_ref_spec = NULL;
        double min_chroma_dist = 0.02;
        double compress_factor = 0.7;
        double smallest_fm[3] = { -1000, -0.2, -1000 };
        bool skip_base_weighting = false;
        struct dcam_ssf *dcam_ssf = NULL;
        const char *exclude_list = NULL;
        const char *neutral_patch = NULL;
        const char *target_layout = NULL;
        bool refine_only_forward_matrix = false;
        bool skip_rebalancing = false;
        bool prefer_cat = true;
        bool render_spectra = false;
        m3x3 *precalc_fm = NULL;
        m3x3 *precalc_cm = NULL;
        m3x3 *precalc_lm = NULL;
        char all_name[4] = "all";
        int w_class_count = 0;
        int v_class_count = 0;
        int l_class_count = 0;
        struct target_adjustment_config *tadj = NULL;
        char *w_class_name[TARGET_MAX_CLASS_COUNT];
        char *v_class_name[TARGET_MAX_CLASS_COUNT];
        char *l_class_name[TARGET_MAX_CLASS_COUNT];
        const int de_w_count = 6;
        double w_class_weights[TARGET_MAX_CLASS_COUNT];
        double v_class_weights[TARGET_MAX_CLASS_COUNT][de_w_count];
        double l_class_weights[TARGET_MAX_CLASS_COUNT][de_w_count];
        memset(w_class_name, 0, sizeof(w_class_name));
        memset(v_class_name, 0, sizeof(v_class_name));
        memset(l_class_name, 0, sizeof(l_class_name));
        double *tc = NULL;
        int tc_len = 0;
        for (int j = 0; j < TARGET_MAX_CLASS_COUNT; j++) {
            w_class_weights[j] = 1.0;
            for (int k = 0; k < de_w_count; k++) v_class_weights[j][k] = 0.0;
            for (int k = 0; k < de_w_count; k++) l_class_weights[j][k] = 0.0;
        }
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'W': skip_base_weighting = true; break;
            case 'w': {
                if (w_class_count == TARGET_MAX_CLASS_COUNT) {
                    elog("Error: too many patch classes.\n");
                    exit(EXIT_FAILURE);
                }
                w_class_name[w_class_count] = argv[++ai];
                w_class_weights[w_class_count] = atof(argv[++ai]);
                w_class_count++;
                if (strchr(argv[ai], ',') != NULL) {
                    elog("Error: only one matrix optimizer weight per class can be specified (-w).\n");
                    exit(EXIT_FAILURE);
                }
                break;
            }
            case 'v': {
                if (v_class_count == TARGET_MAX_CLASS_COUNT) {
                    elog("Error: too many patch classes.\n");
                    exit(EXIT_FAILURE);
                }
                char *p = argv[++ai];
                v_class_name[v_class_count] = p;
                p = argv[++ai];
                double w = atof(p);
                if ((p = strchr(p, ',')) != NULL) {
                    v_class_weights[v_class_count][0] = w;
                    v_class_weights[v_class_count][1] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) v_class_weights[v_class_count][2] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) v_class_weights[v_class_count][3] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) v_class_weights[v_class_count][4] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) v_class_weights[v_class_count][5] = atof(++p);
                    else {
                        v_class_weights[v_class_count][5] = v_class_weights[v_class_count][4];
                        v_class_weights[v_class_count][4] = -v_class_weights[v_class_count][5];
                    }
                } else {
                    if (w < 0) {
                        elog("Error: single DE cannot be negative (-v).\n");
                        exit(EXIT_FAILURE);
                    }
                    v_class_weights[v_class_count][0] = -w;
                    v_class_weights[v_class_count][1] = w;
                    v_class_weights[v_class_count][2] = -w;
                    v_class_weights[v_class_count][3] = w;
                    v_class_weights[v_class_count][4] = -w;
                    v_class_weights[v_class_count][5] = w;
                }
                if (v_class_weights[v_class_count][0] > 0 ||
                    v_class_weights[v_class_count][1] < 0 ||
                    v_class_weights[v_class_count][2] > 0 ||
                    v_class_weights[v_class_count][3] < 0 ||
                    v_class_weights[v_class_count][4] > 0 ||
                    v_class_weights[v_class_count][5] < 0)
                {
                    elog("Error: bad DE range specified (-v).\n");
                    exit(EXIT_FAILURE);
                }
                v_class_count++;
                break;
            }
            case 'l': {
                auto_relax = false;
                if (l_class_count == TARGET_MAX_CLASS_COUNT) {
                    elog("Error: too many patch classes.\n");
                    exit(EXIT_FAILURE);
                }
                char *p = argv[++ai];
                if (isdigit(p[0]) || ((p[0] == '-' || p[0] == '+') && isdigit(p[1]))) {
                    l_class_name[l_class_count] = all_name;
                } else {
                    l_class_name[l_class_count] = p;
                    p = argv[++ai];
                }
                double w = atof(p);
                if ((p = strchr(p, ',')) != NULL) {
                    l_class_weights[l_class_count][0] = w;
                    l_class_weights[l_class_count][1] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) l_class_weights[l_class_count][2] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) l_class_weights[l_class_count][3] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) l_class_weights[l_class_count][4] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) l_class_weights[l_class_count][5] = atof(++p);
                    else {
                        l_class_weights[l_class_count][5] = l_class_weights[l_class_count][4];
                        l_class_weights[l_class_count][4] = -l_class_weights[l_class_count][5];
                    }
                } else {
                    if (w < 0) {
                        elog("Error: single DE cannot be negative (-l).\n");
                        exit(EXIT_FAILURE);
                    }
                    l_class_weights[l_class_count][0] = -w;
                    l_class_weights[l_class_count][1] = w;
                    l_class_weights[l_class_count][2] = -w;
                    l_class_weights[l_class_count][3] = w;
                    l_class_weights[l_class_count][4] = -w;
                    l_class_weights[l_class_count][5] = w;
                }
                if (l_class_weights[l_class_count][0] > 0 ||
                    l_class_weights[l_class_count][1] < 0 ||
                    l_class_weights[l_class_count][2] > 0 ||
                    l_class_weights[l_class_count][3] < 0 ||
                    l_class_weights[l_class_count][4] > 0 ||
                    l_class_weights[l_class_count][5] < 0)
                {
                    elog("Error: bad DE range specified (-l).\n");
                    exit(EXIT_FAILURE);
                }
                l_class_count++;
                break;
            }
            case 'V': refine_only_forward_matrix = true; break;
            case 'a': tadj = jsonio_target_adjustment_parse(argv[++ai], true); break;
            case 't': parse_tonecurve(argv[++ai], &tc, &tc_len); break;
            case 'd': min_chroma_dist = atof(argv[++ai]); break;
            case 'c': dcam_ssf = jsonio_dcam_ssf_parse(argv[++ai], true); break;
            case 'i':
                free(calib_spec);
                calib_spec = NULL;
                calib_wp = parse_arg_illuminant(argv[++ai], observer, &used_observer, &calib_spec);
                if ((int)exif_lightsource_aton(argv[ai]) >= 0) {
                    calib_ls = exif_lightsource_aton(argv[ai]);
                } else {
                    calib_ls = lsOther;
                }
                break;
            case 'I': xyz_ref_wp = parse_arg_illuminant(argv[++ai], observer, &used_observer, &xyz_ref_spec); break;
            case 'C': prefer_cat = false; break;
            case 'S': render_spectra = true; break;
            case 'B': skip_rebalancing = true; break;
            case 'b': neutral_patch = argv[++ai]; break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminants on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'n': camera_name = argv[++ai]; break;
            case 'p': free(precalc_cm); precalc_cm = jsonio_3x3matrix_parse(argv[++ai], "ColorMatrix1", NULL, true); break;
            case 'f': free(precalc_fm); precalc_fm = jsonio_3x3matrix_parse(argv[++ai], "ForwardMatrix1", NULL, true); break;
            case 'e': free(precalc_lm); precalc_lm = jsonio_3x3matrix_parse(argv[++ai], "ForwardMatrix1", "LUTMatrix1", true); break;
            case 'm': {
                ai++;
                free(precalc_cm);
                free(precalc_lm);
                free(precalc_fm);
                struct dcam_profile *prof = jsonio_profile_parse(argv[ai], true);
                precalc_cm = malloc(sizeof(*precalc_cm));
                precalc_lm = malloc(sizeof(*precalc_lm));
                precalc_fm = malloc(sizeof(*precalc_fm));
                *precalc_cm = prof->color_matrix;
                *precalc_lm = prof->lut_matrix;
                *precalc_fm = prof->forward_matrix;
                dcam_profile_delete(prof);
                break;
            }
            case 's': use_slow_matrix_optimizer = true; break;
            case 'L': skip_lut_in_report = true; break;
            case 'x': exclude_list = argv[++ai]; break;
            case 'k': compress_factor = atof(argv[++ai]); break;
            case 'g': target_layout = argv[++ai]; break;
            case 'y': {
                char *p = argv[++ai];
                smallest_fm[0] = -1000;
                smallest_fm[1] = -1000;
                smallest_fm[2] = -1000;
                if (strchr(p, ',') != NULL) {
                    smallest_fm[0] = atof(p);
                    p = strchr(p, ',');
                    if (p) smallest_fm[1] = atof(++p);
                    if (p) p = strchr(p, ',');
                    if (p) smallest_fm[2] = atof(++p);
                } else {
                    smallest_fm[0] = -1000;
                    smallest_fm[1] = atof(p);
                    smallest_fm[2] = -1000;
                }
                break;
            }
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer);
        if (calib_wp.v[0] == 0 && calib_wp.v[1] == 0 && calib_wp.v[2] == 0) {
            calib_wp = glob.d50;
        }
        if (xyz_ref_wp.v[0] == 0 && xyz_ref_wp.v[1] == 0 && xyz_ref_wp.v[2] == 0) {
            xyz_ref_wp = calib_wp;
            xyz_ref_spec = spec_copy(calib_spec);
        }
        if (ai != argc - 2) print_usage_and_exit();
        elog("Reading target...\n");
        struct patch_set *ps = argyllio_ti3_parse(argv[ai], true, false);

        double class_mtx_wt[TARGET_MAX_CLASS_COUNT];
        double class_mtx_de_w[TARGET_MAX_CLASS_COUNT][de_w_count+1];
        double class_lut_de_w[TARGET_MAX_CLASS_COUNT][de_w_count+1];
        for (int j = 0; j < TARGET_MAX_CLASS_COUNT; j++) {
            class_mtx_wt[j] = 1.0;
            for (int k = 0; k < de_w_count+1; k++) class_mtx_de_w[j][k] = 0.0;
            for (int k = 0; k < de_w_count+1; k++) class_lut_de_w[j][k] = 0.0;
        }
        double patch_mtx_wt[ps->patch_count];
        double patch_mtx_de_w[ps->patch_count][de_w_count+1];
        double patch_lut_de_w[ps->patch_count][de_w_count+1];
        for (int j = 0; j < ps->patch_count; j++) {
            patch_mtx_wt[j] = 1.0;
            for (int k = 0; k < de_w_count+1; k++) patch_mtx_de_w[j][k] = 0.0;
            for (int k = 0; k < de_w_count+1; k++) patch_lut_de_w[j][k] = 0.0;
        }
        if (v_class_count > 0) {
            for (int i = 0; i < v_class_count; i++) {
                bool is_class;
                int j = target_class_or_patch_to_idx(ps, v_class_name[i], &is_class);
                if (j == -1) {
                    elog("Error: target \"%s\" doesn't contain class or patch \"%s\", cannot assign weight.\n", argv[ai], v_class_name[i]);
                    exit(EXIT_FAILURE);
                }
                if (is_class) {
                    for (int k = 0; k < de_w_count; k++) class_mtx_de_w[j][k] = v_class_weights[i][k];
                    class_mtx_de_w[j][de_w_count] = 1;
                } else {
                    for (int k = 0; k < de_w_count; k++) patch_mtx_de_w[j][k] = v_class_weights[i][k];
                    patch_mtx_de_w[j][de_w_count] = 1;
                }
            }
        }
        if (l_class_count > 0) {
            if (l_class_count == 1 && strcasecmp(l_class_name[0], "all") == 0) {
                for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
                    for (int k = 0; k < de_w_count; k++) class_lut_de_w[i][k] = l_class_weights[0][k];
                    class_lut_de_w[i][de_w_count] = 1;
                }
            } else {
                for (int i = 0; i < l_class_count; i++) {
                    bool is_class;
                    int j = target_class_or_patch_to_idx(ps, l_class_name[i], &is_class);
                    if (j == -1) {
                        elog("Error: target \"%s\" doesn't contain class or patch \"%s\", cannot assign weight.\n", argv[ai], l_class_name[i]);
                        exit(EXIT_FAILURE);
                    }
                    if (is_class) {
                        for (int k = 0; k < de_w_count; k++) class_lut_de_w[j][k] = l_class_weights[i][k];
                        class_lut_de_w[j][de_w_count] = 1;
                    } else {
                        for (int k = 0; k < de_w_count; k++) patch_lut_de_w[j][k] = l_class_weights[i][k];
                        patch_lut_de_w[j][de_w_count] = 1;
                    }
                }
            }
        }
        if (w_class_count > 0) {
            if (w_class_count == 1 && strcasecmp(w_class_name[0], "all") == 0) {
                for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
                    class_mtx_wt[i] = w_class_weights[0];
                }
            } else {
                for (int i = 0; i < w_class_count; i++) {
                    bool is_class;
                    int j = target_class_or_patch_to_idx(ps, l_class_name[i], &is_class);
                    if (j == -1) {
                        elog("Error: target \"%s\" doesn't contain class or patch \"%s\", cannot assign weight.\n", argv[ai], w_class_name[i]);
                        exit(EXIT_FAILURE);
                    }
                    if (is_class) {
                        class_mtx_wt[j] = w_class_weights[i];
                    } else {
                        patch_mtx_wt[j] = w_class_weights[i];
                    }
                }
            }
        }

        if (dcam_ssf == NULL) {
            if (target_layout != NULL) {
                struct glare_data glare_data;
                struct test_chart_spec *tcs = jsonio_test_chart_spec_parse(target_layout, true);
                match_tcs_with_ps(tcs, ps);

                if (tcs->white_count > 3) {
                    elog("Flat-field correcting target...\n");
                    flatfield_correct_targets(tcs, ps, NULL, true, true, NULL);
                }
                if (tcs->white_count > 0 && tcs->black_count > 0) {
                    elog("Glare test before glare matching...\n");
                    glare_test(ps, &glare_data);
                    elog("Glare-matching target...\n");
                    glare_match(tcs, ps, observer, calib_spec, &calib_wp, NULL, true, true);
                    elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                        "Testing glare after adjusting reference values (camera G and observer Y should\n"
                        "  be close).\n");
                    glare_test(ps, &glare_data);
                } else {
                    glare_test(ps, &glare_data);
                }
                free(tcs);
            } else {
                struct glare_data glare_data;
                glare_test(ps, &glare_data);
            }
        }

        if (exclude_list != NULL) {
            int patch_remapped_indexes[ps->patch_count];
            int orig_count = ps->patch_count;
            int exclude_count = textio_process_exclude_list(ps, exclude_list, patch_remapped_indexes);
            if (exclude_count > 0) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "%d patches was in the provided exclude list and was therefore excluded.\n"
                    "  %d patches left\n",
                    exclude_count, ps->patch_count);
                for (int i = 0; i < orig_count; i++) {
                    int j = patch_remapped_indexes[i]; // j is guaranteed to be equal or less than i
                    if (j < 0) {
                        continue;
                    }
                    if (j < i) {
                        patch_mtx_wt[j] = patch_mtx_wt[i];
                        for (int k = 0; k < de_w_count+1; k++) {
                            patch_mtx_de_w[j][k] = patch_mtx_de_w[i][k];
                            patch_lut_de_w[j][k] = patch_lut_de_w[i][k];
                        }
                    }
                }
            } else {
                elog("No patches matched the provided exclude list, all patches kept.\n");
            }
        }
        if (ps->patch_count < 3) {
            elog("Error: too few patches in target (only %d).\n", ps->patch_count);
            exit(EXIT_FAILURE);
        }

        if (calib_ls != lsD50) {
            elog("Generating values for the calibration illuminant %s and for D50...\n", exif_lightsource_ntoa(calib_ls));
        } else {
            elog("Generating values for the calibration illuminant %s...\n", exif_lightsource_ntoa(calib_ls));
        }
        struct patch_set *psD50 = regenerate_target_values(&ps, argv[ai], calib_ls, calib_spec, calib_wp, xyz_ref_spec, xyz_ref_wp, observer, dcam_ssf, prefer_cat, render_spectra);

        int white_idx = find_neutral_and_rebalance(psD50, glob.d50, neutral_patch, skip_rebalancing);

        print_common_spectra(report_dir, observer, calib_spec, dcam_ssf);
        assign_optimization_weights(ps, class_mtx_wt, (const double(*)[7])class_lut_de_w, patch_mtx_wt, (const double(*)[7])patch_lut_de_w, (const double(*)[7])class_mtx_de_w, (const double(*)[7])patch_mtx_de_w, skip_base_weighting, xyz_ref_wp);
        assign_optimization_weights(psD50, class_mtx_wt, auto_relax ? NULL: (const double(*)[7])class_lut_de_w, patch_mtx_wt, auto_relax ? NULL: (const double(*)[7])patch_lut_de_w, (const double(*)[7])class_mtx_de_w, (const double(*)[7])patch_mtx_de_w, skip_base_weighting, glob.d50);
        if (auto_relax) {
            elog("Automatic LUT relaxation weights assigned.\n");
        }

        if (tadj != NULL) {
            look_apply_target_adjustment(psD50, tadj, glob.d50);
            elog("Applied target adjustment configuration to D50 target.\n");
        }

        elog("Making camera profile...\n");
        struct dcam_profile *prof = make_profile(calib_ls, calib_wp, calib_spec, precalc_cm, precalc_fm, precalc_lm, use_slow_matrix_optimizer, min_chroma_dist, compress_factor, white_idx, ps, psD50, observer, dcam_ssf, smallest_fm, skip_lut_in_report, refine_only_forward_matrix, true, report_dir);
        if (camera_name != NULL) {
            prof->camera_name = strdup(camera_name);
        }

        free(precalc_cm);
        free(precalc_fm);
        free(precalc_lm);
        target_delete_patch_set(ps);
        target_delete_patch_set(psD50);
        elog("Writing output to \"%s\"...\n", argv[ai+1]);

        if (has_suffix(argv[ai+1], ".dcp") || has_suffix(argv[ai+1], ".DCP")) {
            if (tc == NULL) {
                tc = malloc(2 * sizeof(tc[0]));
                tc[0] = 0.0;
                tc[1] = 1.0;
                tc_len = 2;
            }
            const bool enable_3d_hsm = false;
            struct profio_dcp *dcp = make_dcp(prof, prof->camera_name, !skip_lut_in_report, 90, 30, 30, enable_3d_hsm, false, false, false, tc, tc_len, GC_NONE, TC_NEUTRAL, NULL, false, false, NULL);
            profio_dcp_write(dcp, argv[ai+1], true);
            profio_dcp_delete(dcp);

        } else if (has_suffix(argv[ai+1], ".icc") || has_suffix(argv[ai+1], ".icm") || has_suffix(argv[ai+1], ".ICC") || has_suffix(argv[ai+1], ".ICM")) {
            struct profio_icc *icc = make_icc(prof, prof->camera_name, skip_lut_in_report ? ICC_PROFILE_TYPE_MATRIX : ICC_PROFILE_TYPE_LABLUT, false, 33, NULL, 0, tc, tc_len, GC_NONE, TC_NEUTRAL, false, NULL, NULL);
            profio_icc_write(icc, argv[ai+1], true);
            profio_icc_delete(icc);
        } else {
            if (tc != NULL) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "Warning: native format doesn't support embedding a tone-curve. Only provide a\n"
                    "  tone-curve if outputting a DNG (.dcp) or ICC (.icc) profile.\n");
            }
            write_dcam_profile(argv[ai+1], prof, true);
        }
        dcam_profile_delete(prof);
        jsonio_target_adjustment_delete(tadj);
        jsonio_dcam_ssf_delete(dcam_ssf);
        free(calib_spec);
        free(xyz_ref_spec);
        free(tc);
        elog("Complete!\n");

    } else if (strcmp(argv[1], "test-profile") == 0) {

        /******************
          Test profile
         */
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        bool used_observer = false;
        bool skip_lut = false;
        bool skip_look_lut = false;
        bool prefer_cat = true;
        bool render_spectra = false;
        bool add_acr_curve = true;
        v3 custom_wb_;
        v3 *custom_wb = NULL;
        struct dcam_ssf *dcam_ssf = NULL;
        v3 xyz_ref_wp = {{0,0,0}};
        v3 test_ill_wp = {{0,0,0}};
        const char *neutral_patch = NULL;
        bool skip_rebalancing = false;
        enum exif_lightsource test_ill_ls = lsOther;
        spectrum_t *test_ill_spec = NULL;
        spectrum_t *xyz_ref_spec = NULL;
        const char *report_dir = NULL;
        double *trc[3] = { NULL, NULL, NULL };
        int trc_count = 0;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminant on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'w': {
                char *p = argv[++ai];
                bool invert = false;
                if (*p == 'm') {
                    invert = true;
                    p++;
                }
                custom_wb_.v[0] = custom_wb_.v[1] = custom_wb_.v[2] = 0;
                custom_wb_.v[0] = atof(p);
                if ((p = strchr(p, ',')) != NULL) {
                    custom_wb_.v[1] = atof(++p);
                    if ((p = strchr(p, ',')) != NULL) {
                        custom_wb_.v[2] = atof(++p);
                    }
                }
                if (custom_wb_.v[0] <= 0 || custom_wb_.v[1] <= 0 || custom_wb_.v[2] <= 0) {
                    elog("Error: bad white balance value.\n");
                    exit(EXIT_FAILURE);
                }
                if (invert) {
                    custom_wb_.v[0] = 1.0 / custom_wb_.v[0];
                    custom_wb_.v[1] = 1.0 / custom_wb_.v[1];
                    custom_wb_.v[2] = 1.0 / custom_wb_.v[2];
                }
                custom_wb_ = v3_norm(custom_wb_);
                custom_wb = &custom_wb_;
                break;
            }
            case 'B': skip_rebalancing = true; break;
            case 'b': neutral_patch = argv[++ai]; break;
            case 'L': skip_lut = true; break;
            case 'P': skip_look_lut = true; break;
            case 'T': add_acr_curve = false; break;
            case 'c': dcam_ssf = jsonio_dcam_ssf_parse(argv[++ai], true); break;
            case 'i':
                test_ill_wp = parse_arg_illuminant(argv[++ai], observer, &used_observer, &test_ill_spec);
                if ((int)exif_lightsource_aton(argv[ai]) >= 0) {
                    test_ill_ls = exif_lightsource_aton(argv[ai]);
                } else {
                    test_ill_ls = lsOther;
                }
                break;
            case 'I': xyz_ref_wp = parse_arg_illuminant(argv[++ai], observer, &used_observer, &xyz_ref_spec); break;
            case 'C': prefer_cat = false; break;
            case 'S': render_spectra = true; break;
            case 'f': jsonio_transfer_function_parse(argv[++ai], trc, &trc_count, true); break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer);
        if (ai != argc - 2 && ai != argc - 3 && ai != argc - 1) print_usage_and_exit();
        bool has_target = (ai != argc - 1);
        bool has_image = false;
        if (has_target && (has_suffix(argv[ai], ".tif") || has_suffix(argv[ai], ".TIF") || has_suffix(argv[ai], ".tiff") || has_suffix(argv[ai], ".TIFF"))) {
            has_target = false;
            has_image = true;
        }
        bool has_output_image = (ai == argc - 3 || (ai == argc - 2 && (has_suffix(argv[ai+1], ".tif") || has_suffix(argv[ai+1], ".TIF") || has_suffix(argv[ai+1], ".tiff") || has_suffix(argv[ai+1], ".TIFF"))));
        const char *output_image_filename = NULL;
        if (has_output_image) {
            output_image_filename = argv[argc-1];
            if (ai == argc - 2) has_target = false;
        }
        struct prof_wrap pw;
        memset(&pw, 0, sizeof(pw));
        const char *prof_filename = argv[(has_target||has_image) ? ai+1 : ai];
        enum exif_lightsource calib_ls;
        elog("Reading profile...\n");
        if (has_suffix(prof_filename, ".dcp") || has_suffix(prof_filename, ".DCP")) {
            pw.type = PROF_DCP;
            pw.p.dcp = profio_dcp_parse(prof_filename, true);
            calib_ls = pw.p.dcp->ill[0];
            if (add_acr_curve && pw.p.dcp->tonecurve == NULL) {
                int tc_len;
                double *tc = dcp_acr_tonecurve(&tc_len);
                pw.p.dcp->tonecurve = malloc(2 * tc_len * sizeof(pw.p.dcp->tonecurve[0]));
                for (int i = 0; i < tc_len; i++) {
                    pw.p.dcp->tonecurve[2*i+0] = i / (double)(tc_len - 1);
                    pw.p.dcp->tonecurve[2*i+1] = tc[i];
                }
                pw.p.dcp->tonecurve_len = tc_len;
                free(tc);
            }
        } else if (has_suffix(prof_filename, ".icc") || has_suffix(prof_filename, ".ICC") ||
                   has_suffix(prof_filename, ".icm") || has_suffix(prof_filename, ".ICM")) {
            pw.type = PROF_ICC;
            pw.p.icc = profio_icc_parse(prof_filename, true);
            calib_ls = lsUnknown;
        } else {
            pw.type = PROF_NATIVE;
            pw.p.prof = jsonio_profile_parse(prof_filename, true);
            calib_ls = pw.p.prof->illuminant_name;
        }

        if (test_ill_wp.v[0] == 0 && test_ill_wp.v[1] == 0 && test_ill_wp.v[2] == 0) {
            test_ill_ls = calib_ls;
            if (pw.type == PROF_NATIVE && pw.p.prof->illuminant_spec != NULL) {
                test_ill_spec = spec_copy(pw.p.prof->illuminant_spec);
                test_ill_wp = v3_norm(spec2xyz(test_ill_spec, observer));
            } else if (calib_ls != lsUnknown) {
                test_ill_wp = parse_arg_illuminant(exif_lightsource_ntoa(test_ill_ls), observer, &used_observer, &test_ill_spec);
            } else if (pw.type == PROF_ICC) {
                test_ill_wp = pw.p.icc->wp;
            } else if (pw.type == PROF_NATIVE) {
                test_ill_wp = pw.p.prof->illuminant_wp;
            } else {
                elog("Error: could not figure out test illuminant white point.\n");
                exit(EXIT_FAILURE);
            }
        }
        if (xyz_ref_wp.v[0] == 0 && xyz_ref_wp.v[1] == 0 && xyz_ref_wp.v[2] == 0) {
            xyz_ref_wp = test_ill_wp;
        }
        print_common_spectra(report_dir, observer, test_ill_spec, dcam_ssf);
        struct patch_set *ps = NULL, *psD50 = NULL;
        struct rgbimg *rgbimg = NULL;
        int white_idx = -1;
        if (has_target) {
            elog("Reading target...\n");
            ps = argyllio_ti3_parse(argv[ai], true, true);
            psD50 = regenerate_target_values(&ps, argv[ai], test_ill_ls, test_ill_spec, test_ill_wp, xyz_ref_spec, xyz_ref_wp, observer, dcam_ssf, prefer_cat, render_spectra);
            if (trc_count > 0) {
                elog("De-linearizing RGB values by inverting provided transfer function.\n");
                delinearize_rgb(ps, trc, trc_count);
                delinearize_rgb(psD50, trc, trc_count);
            }
            if (ps->patch_count < 3) {
                elog("Error: too few patches in target (only %xd).\n", ps->patch_count);
                exit(EXIT_FAILURE);
            }

            white_idx = find_neutral_and_rebalance(psD50, glob.d50, neutral_patch, skip_rebalancing);
        } else if (has_image) {
            rgbimg = tifio_rgbimg_load(argv[ai], true);
            if (output_image_filename == NULL && report_dir == NULL) {
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "Error: when a test image is provided a an output image filename must also be\n"
                    "  provided, or alternatively a report directory (-r).\n");
                exit(EXIT_FAILURE);
            }
        }
        test_profile(&pw, custom_wb, skip_lut, skip_look_lut, test_ill_wp, test_ill_spec, white_idx, ps, psD50, observer, dcam_ssf, trc, trc_count, rgbimg, output_image_filename, report_dir);
        switch (pw.type) {
        case PROF_DCP: profio_dcp_delete(pw.p.dcp); break;
        case PROF_ICC: profio_icc_delete(pw.p.icc); break;
        case PROF_NATIVE: dcam_profile_delete(pw.p.prof); break;
        }
        free(trc[0]);
        free(trc[1]);
        free(trc[2]);
        free(xyz_ref_spec);
        free(test_ill_spec);
        free(rgbimg);
        jsonio_dcam_ssf_delete(dcam_ssf);
        target_delete_patch_set(ps);
        target_delete_patch_set(psD50);

    } else if (strcmp(argv[1], "make-dcp") == 0) {

        /******************
          Make DCP
         */
        const char *report_dir = NULL;
        const char *camera_name = NULL;
        const char *description = NULL;
        const char *copyright = NULL;
        const char *calibration_signature = NULL;
        double baseline_exposure_offset = 0;
        bool skip_black_render_none = false;
        bool enable_3d_hsm = false;
        bool has_beo = false;
        bool skip_lut = false;
        bool skip_lut_mtx = false;
        bool skip_fwd_mtx = false;
        bool skip_observer_mapping = false;
        bool skip_lut_gamma = false;
        bool allow_discontinuity_hue_shifts = false;
        int hcount = 90;
        int scount = 30;
        int vcount = 30;
        int illuminant_1 = -1;
        int illuminant_2 = -1;
        struct tone_rep_op_config *ntro_conf = NULL;
        struct profio_dcp *wb_dcp = NULL;
        enum gc_type gc_type = GC_FROM_FILE;
        enum tc_type tc_type = TC_NEUTRAL;
        double *tc = NULL;
        int tc_len = 0;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'n': camera_name = argv[++ai]; break;
            case 'd': description = argv[++ai]; break;
            case 'c': copyright = argv[++ai]; break;
            case 's': calibration_signature = argv[++ai]; break;
            case 'b': baseline_exposure_offset = atof(argv[++ai]); has_beo = true; break;
            case 'B': skip_black_render_none = true; break;
            case 'i':
                illuminant_1 = (int)exif_lightsource_aton(argv[++ai]);
                if (illuminant_1 < 0) {
                    elog("Error: unknown EXIF lightsource \"%s\".\n", argv[ai]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'I':
                illuminant_2 = (int)exif_lightsource_aton(argv[++ai]);
                if (illuminant_2 < 0) {
                    elog("Error: unknown EXIF lightsource \"%s\".\n", argv[ai]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                profio_dcp_delete(wb_dcp);
                wb_dcp = profio_dcp_parse(argv[++ai], true);
                break;
            case 'h': {
                char *p = argv[++ai];
                hcount = atoi(p);
                if ((p = strchr(p, ',')) == NULL) {
                    elog("Error: missing saturation division value (expected huediv,satdiv,valdiv example: \"-h 90,30,30\").\n");
                    exit(EXIT_FAILURE);
                }
                scount = atoi(++p);
                if ((p = strchr(p, ',')) == NULL) {
                    elog("Error: missing value division value (expected huediv,satdiv,valdiv example: \"-h 90,30,30\").\n");
                    exit(EXIT_FAILURE);
                }
                vcount = atoi(++p);
                if (hcount <= 0 || scount <= 0 || vcount < 0 || hcount > 10000 || scount > 10000 || vcount > 10000) {
                    elog("Error: invalid hue, saturation or value division.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            }
            case 'O': skip_observer_mapping = true; break;
            case 'F': skip_fwd_mtx = true; break;
            case 'E': skip_lut_mtx = true; break;
            case 'L': skip_lut = true; break;
            case 'H': allow_discontinuity_hue_shifts = true; break;
            case 'D': enable_3d_hsm = true; break;
            case 'r': report_dir = argv[++ai]; break;
            case 't': parse_tonecurve(argv[++ai], &tc, &tc_len); break;
            case 'o':
                ai++;
                if (strcasecmp(argv[ai], "neutral") == 0) {
                    tc_type = TC_NEUTRAL;
                } else if (strcasecmp(argv[ai], "standard") == 0) {
                    tc_type = TC_STANDARD;
                } else {
                    tc_type = TC_NEUTRAL;
                    ntro_conf = jsonio_ntro_conf_parse(argv[ai], true);
                }
                break;
            case 'g': gc_type = parse_gc_type(argv[++ai]); break;
            case 'G': skip_lut_gamma = true; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer_get(OBSERVER_1931_2));
        test_gc_type_setting(gc_type, tc_type, tc, skip_lut);
        if (tc == NULL) {
            tc = malloc(2 * sizeof(tc[0]));
            tc[0] = 0.0;
            tc[1] = 1.0;
            tc_len = 2;
        }
        if (ai != argc - 2 && ai != argc - 3) print_usage_and_exit();
        const char *outname = ai == argc - 3 ? argv[ai+2] : argv[ai+1];
        struct dcam_profile *prof = jsonio_profile_parse(argv[ai], true);
        struct profio_dcp *dcp = make_dcp(prof, camera_name, !skip_lut, hcount, scount, vcount, enable_3d_hsm, skip_fwd_mtx, skip_lut_mtx, skip_observer_mapping, tc, tc_len, gc_type, tc_type, ntro_conf, skip_lut_gamma, allow_discontinuity_hue_shifts, report_dir);
        free(tc);
        dcam_profile_delete(prof);
        if (ai == argc - 3) {
            prof = jsonio_profile_parse(argv[ai+1], true);
            struct profio_dcp *dcp2 = make_dcp(prof, camera_name, !skip_lut, hcount, scount, vcount, enable_3d_hsm, skip_fwd_mtx, skip_lut_mtx, skip_observer_mapping, NULL, 0, GC_NONE, TC_NEUTRAL, NULL, false, allow_discontinuity_hue_shifts, report_dir);
            dcam_profile_delete(prof);
            dcp->ill_count = 2;
            dcp->ill[1] = dcp2->ill[0];
            dcp->cm[1] = dcp2->cm[0];
            dcp->fm[1] = dcp2->fm[0];
            dcp->huesatmap[1] = dcp2->huesatmap[0];
            dcp2->huesatmap[0] = NULL;
            profio_dcp_delete(dcp2);
        }
        if (description != NULL) {
            free(dcp->profile_name);
            dcp->profile_name = strdup(description);
        } else {
            if (tc_type == TC_STANDARD) {
                dcp->profile_name = strdup("Standard");
            }
            if (tc_type == TC_STANDARD) {
               if (ntro_conf == NULL) {
                   dcp->profile_name = strdup("Neutral");
               } else {
                   dcp->profile_name = strdup("Custom");
               }
            }
        }
        if (copyright != NULL) {
            free(dcp->copyright);
            dcp->copyright = strdup(copyright);
        }
        if (calibration_signature != NULL) {
            free(dcp->calibration_signature);
            dcp->calibration_signature = strdup(calibration_signature);
        }
        if (has_beo) {
            dcp->baseline_exposure_offset = baseline_exposure_offset;
            dcp->has.baseline_exposure_offset = true;
        }
        if (skip_black_render_none) {
            dcp->black_render_none = false;
            dcp->has.black_render_none = false;
        }
        if (illuminant_1 != -1) {
            dcp->ill[0] = (enum exif_lightsource)illuminant_1;
        }
        if (illuminant_2 != -1 && dcp->ill_count > 1) {
            dcp->ill[1] = (enum exif_lightsource)illuminant_2;
        }
        if (wb_dcp != NULL) {
            if (dcp->ill_count == 2) {
                if (wb_dcp->ill_count == 2) {
                    // we need to copy illuminants as well to get the same mixing result
                    dcp->ill[0] = wb_dcp->ill[0];
                    dcp->ill[1] = wb_dcp->ill[1];
                    dcp->cm[0] = wb_dcp->cm[0];
                    dcp->cm[1] = wb_dcp->cm[1];
                } else {
                    // as both matrices are the same we can keep illuminants
                    dcp->cm[0] = wb_dcp->cm[0];
                    dcp->cm[1] = wb_dcp->cm[0];
                }
            } else {
                if (wb_dcp->ill_count == 2) {
                    // we need to make a fake dual illuminant profile to get the same WB result
                    elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                        "Warning: due to dual-illuminant WB profile and single-illuminant profile, the\n"
                        "  output profile is converted to (fake) dual-illuminant to match WB.\n");
                    dcp->ill_count = 2;
                    dcp->ill[0] = wb_dcp->ill[0];
                    dcp->ill[1] = wb_dcp->ill[1];
                    dcp->cm[0] = wb_dcp->cm[0];
                    dcp->cm[1] = wb_dcp->cm[1];
                    dcp->fm[1] = dcp->fm[0];
                    if (dcp->huesatmap[0] != NULL) {
                        size_t sz = dcp->hsmdims[0] * dcp->hsmdims[1] * dcp->hsmdims[2] * sizeof(dcp->huesatmap[0][0]);
                        dcp->huesatmap[1] = malloc(sz);
                        memcpy(dcp->huesatmap[1], dcp->huesatmap[0], sz);
                    }
                } else {
                    // as there's only one matrix, we can keep illuminant
                    dcp->cm[0] = wb_dcp->cm[0];
                }
            }
            profio_dcp_delete(wb_dcp);
        }
        jsonio_ntro_delete(ntro_conf);
        print_matrices(report_dir, &dcp->cm[0], dcp->has.fm ? &dcp->fm[0] : NULL, NULL, dcp->ill_count == 2 ? &dcp->cm[1] : NULL, (dcp->has.fm && dcp->ill_count == 2) ? &dcp->fm[1] : NULL, dcp_d50());
        elog("Writing output to \"%s\"...\n", outname);
        profio_dcp_write(dcp, outname, true);
        profio_dcp_delete(dcp);
        elog("Complete!\n");

    } else if (strcmp(argv[1], "make-icc") == 0) {

        /******************
          Make ICC
         */
        const char *report_dir = NULL;
        const char *camera_name = NULL;
        const char *copyright = NULL;
        bool apply_wb = false;
        bool skip_input_lut = false;
        bool skip_tc_apply = false;
        enum icc_profile_type icc_type = ICC_PROFILE_TYPE_UNDEFINED;
        int clut_side = 33;
        double *trc[3] = { NULL, NULL, NULL };
        int trc_count = 0;
        struct tone_rep_op_config *ntro_conf = NULL;
        enum gc_type gc_type = GC_FROM_FILE;
        enum tc_type tc_type = TC_NEUTRAL;
        double *tc = NULL;
        int tc_len = 0;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'n': camera_name = argv[++ai]; break;
            case 'c': copyright = argv[++ai]; break;
            case 's': clut_side = atoi(argv[++ai]); break;
            case 'L': skip_input_lut = true; break;
            case 'W': apply_wb = true; break;
            case 'f': jsonio_transfer_function_parse(argv[++ai], trc, &trc_count, true); break;
            case 'p':
                ai++;
                if (strcasecmp(argv[ai], "lablut") == 0) {
                    icc_type = ICC_PROFILE_TYPE_LABLUT;
                } else if (strcasecmp(argv[ai], "xyzlut") == 0) {
                    icc_type = ICC_PROFILE_TYPE_XYZLUT;
                } else if (strcasecmp(argv[ai], "matrix") == 0) {
                    icc_type = ICC_PROFILE_TYPE_MATRIX;
                } else {
                    elog("Error: unknown profile type \"%s\", expected \"lablut\", \"xyzlut\" or \"matrix\".\n", argv[ai]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 't': parse_tonecurve(argv[++ai], &tc, &tc_len); break;
            case 'T': skip_tc_apply = true; break;
            case 'o':
                ai++;
                if (strcasecmp(argv[ai], "neutral") == 0) {
                    tc_type = TC_NEUTRAL;
                } else if (strcasecmp(argv[ai], "standard") == 0) {
                    tc_type = TC_STANDARD;
                } else {
                    tc_type = TC_NEUTRAL;
                    ntro_conf = jsonio_ntro_conf_parse(argv[ai], true);
                }
                break;
            case 'g': gc_type = parse_gc_type(argv[++ai]); break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer_get(OBSERVER_1931_2));
        test_gc_type_setting(gc_type, tc_type, tc, skip_input_lut);
        if (ai != argc - 2) print_usage_and_exit();
        struct dcam_profile *prof = jsonio_profile_parse(argv[ai], true);
        if (skip_input_lut) {
            chromalut_delete(prof->lut);
            prof->lut = NULL;
        }
        if (icc_type == ICC_PROFILE_TYPE_UNDEFINED) {
            icc_type = prof->lut != NULL ? ICC_PROFILE_TYPE_LABLUT : ICC_PROFILE_TYPE_MATRIX;
        }
        struct profio_icc *icc = make_icc(prof, camera_name, icc_type, apply_wb, clut_side, trc, trc_count, tc, tc_len, gc_type, tc_type, skip_tc_apply, ntro_conf, report_dir);
        free(tc);
        print_transfer_function(report_dir, trc, trc_count);
        dcam_profile_delete(prof);
        if (copyright != NULL) {
            free(icc->copyright);
            icc->copyright = strdup(copyright);
        }
        elog("Writing output to \"%s\"...\n", argv[ai+1]);
        profio_icc_write(icc, argv[ai+1], true);
        elog("Complete!\n");
        profio_icc_delete(icc);
        free(trc[0]);
        free(trc[1]);
        free(trc[2]);
        jsonio_ntro_delete(ntro_conf);

    } else if (strcmp(argv[1], "make-testchart") == 0) {

        /******************
          Make test chart
         */
        const char *report_dir = NULL;
        int patch_count = 100;
        int black_count = 5;
        int gray_steps = 0;
        double white_percent = 20.0;
        struct test_chart_layout layout = {
            .row_count = -1,
            .has_col_offset = false,
            .row_height = 1.0,
            .col_width = 1.0
        };
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'p': patch_count = atoi(argv[++ai]); break;
            case 'w': white_percent = atof(argv[++ai]); break;
            case 'b': black_count = atoi(argv[++ai]); break;
            case 'g': gray_steps = atoi(argv[++ai]); break;
            case 'l': layout.row_count = atoi(argv[++ai]); break;
            case 'd': {
                char *p = argv[++ai];
                layout.row_height = atof(p);
                if ((p = strchr(p, ',')) == NULL) {
                    elog("Error: missing col width (expected row,col).\n");
                    exit(EXIT_FAILURE);
                }
                layout.col_width = atof(++p);
                break;
            }
            case 'O': layout.has_col_offset = true; break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        if (gray_steps < 0) {
            gray_steps = 0;
        }
        int white_count = (patch_count * white_percent) / 100;
        if (white_count < 1) white_count = 1;
        if (white_count >= patch_count) white_count = patch_count - 1;
        if (ai != argc - 1) print_usage_and_exit();

        elog("Generating %d patches...\n", patch_count);
        while (patch_count - black_count + gray_steps * black_count - white_count < black_count && black_count > 0) {
            black_count--;
        }
        if (black_count > white_count) black_count = white_count;
        int neutral_count = white_count + black_count + gray_steps * black_count;
        int nstep_count = 2 + gray_steps;
        int nstep_counts[nstep_count];
        nstep_counts[0] = white_count;
        for (int i = 1; i < nstep_count; i++) nstep_counts[i] = black_count;

        struct patch_set *ps = calloc(1, sizeof(*ps) + patch_count * sizeof(ps->patch[0]));
        const m3x3 rgb2xyz = {{
                { 0.4124564,  0.3575761,  0.1804375 },
                { 0.2126729,  0.7151522,  0.0721750 },
                { 0.0193339,  0.1191920,  0.9503041 }
            }};
        ps->patch[ps->patch_count].cam = (v3){{1,1,1}};
        ps->patch[ps->patch_count].xyz = m3x3_mul_v3(rgb2xyz, (v3){{1,1,1}});
        ps->patch_count++;
        for (int i = 1; i < patch_count - neutral_count + 1; i++) {
            v3 rgb = find_distant_light_color(ps, rgb2xyz);
            ps->patch[ps->patch_count].cam = rgb;
            ps->patch[ps->patch_count].xyz = m3x3_mul_v3(rgb2xyz, rgb);
            ps->patch_count++;
            elog(".");
        }
        elog("done!\n");
        if (layout.row_count <= 0) {
            layout.row_count = round(sqrt(patch_count));
            if (layout.row_count == 0) layout.row_count = 1;
        }
        {
            int color_idx[nstep_count][neutral_count];
            int c_count[nstep_count];
            memset(c_count, 0, sizeof(nstep_count));
            int has_color[patch_count]; // assigned color
            memset(has_color, 0, sizeof(has_color));
            has_color[0] = 1; // white
            color_idx[0][0] = 0; // assigned patch index
            c_count[0] = 1;
            if (black_count > 0) {
                for (int c = 1; c < nstep_count; c++) {
                    has_color[c] = c+1; // assign step wedge to first patches
                    color_idx[c][0] = c;
                    c_count[c] = 1;
                }
            }
            for (int c = 0; c < nstep_count; c++) {
                for (int i = 1; i < nstep_counts[c]; i++) {
                    double max_dist = 0;
                    int max_dist_idx = 0;
                    for (int j = 0; j < patch_count; j++) {
                        if (has_color[j] != 0) {
                            continue;
                        }
                        // find minimum distance from current position to already assigned patch of same color
                        double min_dist = -1;
                        for (int k = 0; k < c_count[c]; k++) {
                            double dist = test_chart_dist(&layout, j, color_idx[c][k]);
                            if (min_dist == -1 || dist < min_dist) {
                                min_dist = dist;
                            }
                        }
                        if (min_dist > max_dist) {
                            max_dist = min_dist;
                            max_dist_idx = j;
                        }
                    }
                    color_idx[c][c_count[c]] = max_dist_idx;
                    has_color[max_dist_idx] = c + 1;
                    c_count[c]++;
                }
            }
            struct patch_set *ps1 = calloc(1, sizeof(*ps1) + patch_count * sizeof(ps1->patch[0]));
            int k = 1;
            for (int i = 0; i < patch_count; i++) {
                if (has_color[i] == 0) {
                    ps1->patch[ps1->patch_count] = ps->patch[k++];
                } else {
                    double val = (double)(nstep_count - has_color[i]) / (nstep_count - 1);
                    ps1->patch[ps1->patch_count].cam = (v3){{val,val,val}};
                    ps1->patch[ps1->patch_count].xyz = m3x3_mul_v3(rgb2xyz, (v3){{val,val,val}});
                }
                ps1->patch_count++;
            }
            target_delete_patch_set(ps);
            ps = ps1;
        }
        argyllio_ti1_write(argv[ai], ps, true, true);
        target_print("target", ps, ps->patch[0].xyz, observer_get(OBSERVER_1931_2), report_dir);
        target_delete_patch_set(ps);

    } else if (strcmp(argv[1], "testchart-ff") == 0) {

        /******************
          Test chart flatfield
         */
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        bool illuminant_given = false;
        bool used_observer = false;
        spectrum_t *ill_spec = NULL;
        v3 ill_wp = {{0,0,0}};
        const char *report_dir = NULL;
        bool enable_linearization = false;
        struct test_chart_layout layout = {
            .row_count = -1,
            .has_col_offset = false,
            .row_height = 1.0,
            .col_width = 1.0
        };
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'l': layout.row_count = atoi(argv[++ai]); break;
            case 'd': {
                char *p = argv[++ai];
                layout.row_height = atof(p);
                if ((p = strchr(p, ',')) == NULL) {
                    elog("Error: missing col width (expected row,col).\n");
                    exit(EXIT_FAILURE);
                }
                layout.col_width = atof(++p);
                break;
            }
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminant on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'i':
                ill_wp = parse_arg_illuminant(argv[++ai], observer, &used_observer, &ill_spec);
                illuminant_given = true;
                break;
            case 'O': layout.has_col_offset = true; break;
            case 'L': enable_linearization = true; break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }

        if (has_suffix(argv[ai], ".tif") || has_suffix(argv[ai], ".TIF") || has_suffix(argv[ai], ".tiff") || has_suffix(argv[ai], ".TIFF")) {
            if (enable_linearization) {
                elog("Error: glare modeling of TIFF is not supported (remove the -L flag).\n");
                exit(EXIT_FAILURE);
            }
            if (ai != argc - 3) print_usage_and_exit();
            tifio_flatfield_correct(argv[ai], argv[ai+1], argv[ai+2], true);
            exit(EXIT_SUCCESS);
        }

        if (ai != argc - 3 && ai != argc - 4) print_usage_and_exit();

        struct test_chart_spec *tcs;
        bool tcs_from_json = false;
        if (has_suffix(argv[ai], ".ti1") || has_suffix(argv[ai], ".TI1")) {
            if (layout.row_count <= 0) {
                elog("Error: row count (-l) must be specified when reading .ti1.\n");
                exit(EXIT_FAILURE);
            }
            tcs = ti1_to_tcs(&layout, argv[ai]);
        } else {
            tcs = jsonio_test_chart_spec_parse(argv[ai], true);
            tcs_from_json = true;
        }
        elog("Found %d white patch%s, %d black patch%s and %d gray patch%s.\n", tcs->white_count, tcs->white_count != 1 ? "es" : "", tcs->black_count, tcs->black_count != 1 ? "es" : "", tcs->gray_count, tcs->gray_count != 1 ? "es" : "");
        struct patch_set *ps = argyllio_ti3_parse(argv[ai+1], true, true);
        if (ps->patch_count != tcs->patch_count) {
            elog("Error: mismatching patch count between \"%s\" and \"%s\", not the same target layout?\n", argv[ai], argv[ai+1]);
            exit(EXIT_FAILURE);
        }
        if (tcs_from_json) {
            match_tcs_with_ps(tcs, ps);
        }

        struct patch_set *ps_ex = NULL;
        char *out_filename;
        if (ai == argc - 4) {
            ps_ex = argyllio_ti3_parse(argv[ai+2], true, true);
            if (ps->patch_count != tcs->patch_count) {
                elog("Error: mismatching patch count between \"%s\" and \"%s\", not the same target layout?\n", argv[ai], argv[ai+2]);
                exit(EXIT_FAILURE);
            }
            out_filename = argv[ai+3];
        } else {
            out_filename = argv[ai+2];
        }

        if (tcs->white_count >= 3) {
            elog("Performing flatfield correction.\n");
            flatfield_correct_targets(tcs, ps, ps_ex, true, true, NULL);
        } else {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Warning: only %d white patches, flatfield correction skipped (requires at\n"
                "  least 3 whites).\n", tcs->white_count);
        }

        struct glare_data glare_data;
        elog("Testing for glare.\n");
        glare_test(ps, &glare_data);
        if (ps_ex != NULL) {
            target_delete_patch_set(ps);
            ps = normalize_patch_set(ps_ex, -1);
            target_delete_patch_set(ps_ex);
        } else {
            struct patch_set *ps1 = normalize_patch_set(ps, -1);
            target_delete_patch_set(ps);
            ps = ps1;
        }

        if (enable_linearization) {
            if (!illuminant_given) {
                elog("Error: illuminant must be provided to model glare (eg -i D50).\n");
                exit(EXIT_FAILURE);
            }
            elog("Adding glare to reference values to match glare in camera RGB values.\n");
            glare_match(tcs, ps, observer, ill_spec, &ill_wp, report_dir, true, true);
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Testing glare after adjusting reference values (camera G and observer Y should\n"
                "  be close).\n");
            glare_test(ps, &glare_data);
        }

        elog("Writing corrected output to \"%s\".\n", out_filename);
        argyllio_ti3_write(out_filename, ps->patch[0].spectrum != NULL, ps, true);
        target_delete_patch_set(ps);

    } else if (strcmp(argv[1], "tiff-tf") == 0) {

        /******************
          Tiff transfer functions
         */
        double *trc[3] = { NULL, NULL, NULL };
        int trc_count = 0;
        double *tf_trc[3] = { NULL, NULL, NULL };
        int tf_trc_count = 0;
        bool no_reconstruction = false;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'R': no_reconstruction = true; break;
            case 'f': jsonio_transfer_function_parse(argv[++ai], tf_trc, &tf_trc_count, true); break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        if (ai != argc - 2 && ai != argc - 1) print_usage_and_exit();
        tifio_get_transfer_function(argv[ai], trc, &trc_count, true);
        if (trc_count == 0) {
            elog("Error: no transfer function (TIFFTAG_TRANSFERFUNCTION) found in \"%s\".\n", argv[ai]);
            exit(EXIT_FAILURE);
        }
        if (!no_reconstruction) {
#pragma omp parallel
            for (int i = 0; i < 3; i++) {
                reconstruct_tonecurve(trc[i], trc_count);
            }
            if (tf_trc_count > 0) {
#pragma omp parallel
                for (int i = 0; i < 3; i++) {
                    reconstruct_tonecurve(tf_trc[i], tf_trc_count);
                }
            }
        }
        if (tf_trc_count > 0) {
            const char rgb[3] = { 'R', 'G', 'B' };
            for (int i = 0; i < 3; i++) {
                double *new_curve = malloc(sizeof(double) * trc_count);
                for (int j = 0; j < trc_count; j++) {
                    double x = (double)j / (trc_count - 1);
                    x = inverse_curve(x, trc[i], trc_count); // de-linerization curve, ie the tone curve including the reference transfer function
                    new_curve[j] = curve_y(x, tf_trc[i], tf_trc_count); // linearize the reference transfer function part, result tone curve alone
                }
                double mx = new_curve[trc_count - 1];
                double mn = new_curve[0];
                if (mx != 1.0) {
                    elog("Warning: curve %c ends at   %.10f, compensating so it will end at 1.0\n", rgb[i], mx);
                }
                if (mn != 0) {
                    if (mn > 0.01) {
                        elog("Warning: curve %c starts at %.10f, so large it's not a rounding error, will not compensate.\n", rgb[i], mn);
                        mn = 0;
                    } else {
                        elog("Warning: curve %c starts at %.10f, compensating so it will start at 0.0\n", rgb[i], mn);
                    }
                }
                mx = 1.0 / mx;
                mx = mx / (1.0 - mn);
                for (int j = 0; j < trc_count; j++) {
                    trc[i][j] = (new_curve[j] - mn) * mx;
                }
                reconstruct_tonecurve(trc[i], trc_count);
                free(new_curve);
            }
        }
        FILE *out = stdout;
        if (ai == argc - 2) {
            out = utf8_fopen(argv[ai+1], "w");
            if (out == NULL) {
                elog("Error: could not open \"%s\" for writing: %s.\n", argv[ai+1], strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
        fprintf(out, "{\n");
        const char *trcname[3] = { "RedTRC", "GreenTRC", "BlueTRC" };
        for (int i = 0; i < 3; i++) {
            fprintf(out,
                    "  \"%s\":", trcname[i]);
            if (trc_count > 1) {
                fprintf(out, " [\n");
                for (int j = 0; j < trc_count; j++) {
                    fprintf(out, "    %.8f", trc[i][j]);
                    if (j < trc_count-1) fprintf(out, ",\n");
                }
                fprintf(out, "\n  ]");
            } else {
                fprintf(out, " %f", trc[i][0]);
            }
            if (i < 2) fprintf(out, ",\n");
        }
        fprintf(out, "\n}\n");
        if (out != stdout) fclose(out);
        free(trc[0]);
        free(trc[1]);
        free(trc[2]);
        free(tf_trc[0]);
        free(tf_trc[1]);
        free(tf_trc[2]);

    } else if (strcmp(argv[1], "average-targets") == 0) {

        /******************
          Average targets
         */
        if (ai > argc - 2) {
            print_usage_and_exit();
        }
        elog("Reading \"%s\"...\n", argv[ai]);
        struct patch_set *ps = argyllio_ti3_parse(argv[ai], true, true);
        for (int i = 1; i < argc - 3; i++) {
            elog("Reading \"%s\"...\n", argv[ai+i]);
            struct patch_set *psin = argyllio_ti3_parse(argv[ai+i], true, true);
            if (psin->patch_count != ps->patch_count) {
                elog("Error: different patch count between \"%s\" and \"%s\". Different targets?\n", argv[ai], argv[ai+i]);
                exit(EXIT_FAILURE);
            }
            for (int j = 0; j < ps->patch_count; j++) {
                for (int k = 0; k < 3; k++) {
                    ps->patch[j].cam.v[k] += psin->patch[j].cam.v[k];
                    ps->patch[j].xyz.v[k] += psin->patch[j].xyz.v[k];
                }
                if (ps->patch[j].spectrum != NULL) {
                    if (psin->patch[j].spectrum == NULL) {
                        elog("Error: target \"%s\" has spectra, but not \"%s\".\n", argv[ai], argv[ai+i]);
                        exit(EXIT_FAILURE);
                    }
                    spectrum_t *s = spec_add(ps->patch[j].spectrum, psin->patch[j].spectrum);
                    free(ps->patch[j].spectrum);
                    ps->patch[j].spectrum = s;
                }
            }
            target_delete_patch_set(psin);
        }
        elog("Averaging all read targets...\n");
        double f = 1.0 / (argc - 3);
        for (int j = 0; j < ps->patch_count; j++) {
            for (int k = 0; k < 3; k++) {
                ps->patch[j].cam.v[k] *= f;
                ps->patch[j].xyz.v[k] *= f;
            }
            if (ps->patch[j].spectrum != NULL) {
                spec_scalar_multiply(ps->patch[j].spectrum, f);
            }
        }
        elog("Writing output to \"%s\"...\n", argv[argc-1]);
        argyllio_ti3_write(argv[argc-1], true, ps, true);
        target_delete_patch_set(ps);
        elog("Complete!\n");

    } else if (strcmp(argv[1], "match-spectra") == 0) {

        /******************
          Match spectra

          match-spectra [flags] <reference.ti3> <match.ti3> <output-match.ti3> [<output-ref.ti3>]
         */
        struct dcam_ssf *dcam_ssf = NULL;
        const char *report_dir = NULL;
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        double acceptable_error = -1;
        bool normalized_comparison = false;
        bool used_observer = false;
        bool scale_spectra = false;
        bool unique_output = false;
        bool force_emissive = false;
        spectrum_t *illuminant = spec_copy(spectraldb_illuminant(lsD50));
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'i':
                free(illuminant);
                illuminant = NULL;
                parse_arg_illuminant(argv[++ai], observer, &used_observer, &illuminant);
                if (illuminant == NULL) {
                    elog("Error: must provide illuminant with spectrum.\n");
                }
                break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminants on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'c': dcam_ssf = jsonio_dcam_ssf_parse(argv[++ai], true); break;
            case 'S': scale_spectra = true; break;
            case 'N': normalized_comparison = true; break;
            case 'U': unique_output = true; break;
            case 'E': force_emissive = true; break;
            case 'e': acceptable_error = atof(argv[++ai]); break;
            case 'r': report_dir = argv[++ai]; break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer);
        if (ai != argc - 3 && ai != argc - 4) print_usage_and_exit();
        struct patch_set *psref = argyllio_ti3_parse(argv[ai], true, true);
        struct patch_set *psm = argyllio_ti3_parse(argv[ai+1], true, true);
        struct patch_set *ps = calloc(1, sizeof(*ps) + sizeof(ps->patch[0]) * (psm->patch_count + psref->patch_count));
        struct patch_set *psrm = calloc(1, sizeof(*ps) + sizeof(ps->patch[0]) * (psm->patch_count + psref->patch_count));
        *ps = *psm;
        *psrm = *psref;
        ps->patch_count = 0;
        psrm->patch_count = 0;
        v3 illuminant_wp = v3_norm(spec2xyz(illuminant, observer));
        const struct observer *obs = dcam_ssf != NULL ? &dcam_ssf->ssf : observer;
        if (dcam_ssf != NULL) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Using camera SSFs instead of standard observer for comparisons, will also use\n"
                "  euclidean XYZ distance instead of CIEDE2000.\n");
        }
        struct {
            int idx;
            double de;
            v3 xyz_ref;
            v3 xyz;
            double scale_factor;
        } best_match[psref->patch_count];
        for (int i = 0; i < psref->patch_count; i++) {
            if (psref->patch[i].spectrum == NULL) {
                elog("Error: missing spectrum in reference target.\n");
                exit(EXIT_FAILURE);
            }
            v3 xyz_ref;
            if (force_emissive || psref->patch[i].emissive) {
                xyz_ref = spec2xyz(psref->patch[i].spectrum, obs);
            } else {
                xyz_ref = spec2xyz_ill(psref->patch[i].spectrum, obs, illuminant, 1);
            }
            double min_de = 100000;
            double min_sf = 1;
            int min_j = 0;
            v3 min_xyz = {{0,0,0}};
            for (int j = 0; j < psm->patch_count; j++) {
                if (psm->patch[j].spectrum == NULL) {
                    elog("Error: missing spectrum in matching target.\n");
                    exit(EXIT_FAILURE);
                }
                v3 xyz;
                if (force_emissive || psm->patch[j].emissive) {
                    xyz = spec2xyz(psm->patch[j].spectrum, obs);
                } else {
                    xyz = spec2xyz_ill(psm->patch[j].spectrum, obs, illuminant, 1);
                }
                double de, scale_factor;
                de = compare_tristimulus_values(xyz, xyz_ref, illuminant_wp, normalized_comparison, dcam_ssf == NULL, &scale_factor);
                if (de < min_de) {
                    min_j = j;
                    min_de = de;
                    min_xyz = xyz;
                    min_sf = scale_factor;
                }
            }
            best_match[i].idx = min_j;
            best_match[i].de = min_de;
            best_match[i].xyz = min_xyz;
            best_match[i].xyz_ref = xyz_ref;
            best_match[i].scale_factor = min_sf;
        }

        for (int i = 0; i < psref->patch_count; i++) {
            if (psref->patch[i].str_id[0] == '\0') sprintf(psref->patch[i].str_id, "%d", i+1);
            if (psm->patch[best_match[i].idx].str_id[0] == '\0') sprintf(psm->patch[best_match[i].idx].str_id, "%d", i+1);
            bool skip_patch = false;
            bool bad_match = false;
            if (acceptable_error >= 0 && best_match[i].de > acceptable_error) {
                skip_patch = true;
                bad_match = true;
            }
            if (!skip_patch && unique_output) {
                int best_i = i;
                for (int j = 0; j < psref->patch_count; j++) {
                    if (best_match[j].idx == best_match[i].idx) {
                        if (best_match[j].de < best_match[best_i].de) {
                            best_i = j;
                        }
                    }
                }
                if (i != best_i) {
                    skip_patch = true;
                }
            }
            elog("Best match for patch \"%s\" in reference target is patch \"%s\" in matching target. Matching result: %.2f DE%s\n", psref->patch[i].str_id, psm->patch[best_match[i].idx].str_id, best_match[i].de, bad_match ? " (not good enough, excluded)" : skip_patch ? " (better exists, excluded)" : "");
            if (!skip_patch) {
                ps->patch[ps->patch_count] = psm->patch[best_match[i].idx];
                ps->patch[ps->patch_count].spectrum = spec_copy(psm->patch[best_match[i].idx].spectrum);
                if (scale_spectra) {
                    spec_scalar_multiply(ps->patch[ps->patch_count].spectrum, best_match[i].scale_factor);
                    ps->patch[ps->patch_count].xyz = k_mul_v3(best_match[i].scale_factor, best_match[i].xyz);
                }
                ps->patch_count++;
                psrm->patch[psrm->patch_count] = psref->patch[i];
                psrm->patch[psrm->patch_count].spectrum = spec_copy(psref->patch[i].spectrum);
                psrm->patch_count++;
            }
        }
        free(illuminant);
        target_print("target-ref", psref, illuminant_wp, observer, report_dir);
        target_print("target-match", psm, illuminant_wp, observer, report_dir);
        target_print("target", ps, illuminant_wp, observer, report_dir);
        target_print("target-refm", psrm, illuminant_wp, observer, report_dir);
        target_delete_patch_set(psref);
        target_delete_patch_set(psm);
        argyllio_ti3_write(argv[ai+2], true, ps, true);
        if (ai == argc - 4) {
            argyllio_ti3_write(argv[ai+3], true, psrm, true);
        }
        target_delete_patch_set(psrm);
        target_delete_patch_set(ps);
        jsonio_dcam_ssf_delete(dcam_ssf);

    } else if (strcmp(argv[1], "txt2ti3") == 0) {

        /******************
          Text to ti3
         */
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        const spectrum_t *illuminant = spectraldb_illuminant(lsD50);
        const char *class_name = "import";
        double scale = 1;
        double clip_level = 0;
        int low = 300;
        int high = 700;
        int spacing = 2;
        bool column_layout = true;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'f': {
                char *p = argv[++ai];
                low = atoi(p);
                if ((p = strchr(p, ',')) != NULL) {
                    high = atoi(++p);
                    if ((p = strchr(p, ',')) != NULL) {
                        spacing = atoi(++p);
                    }
                }
                break;
            }
            case 's': scale = atof(argv[++ai]); break;
            case 'm': clip_level = atof(argv[++ai]); break;
            case 'a': class_name = argv[++ai]; break;
            case 'l':
                if (strcasecmp(argv[ai+1], "cols") == 0) {
                    column_layout = true;
                } else if (strcasecmp(argv[ai+1], "rows") == 0) {
                    column_layout = false;
                } else {
                    elog("Error: unknown layout \"%s\", expected \"cols\" or \"rows\".\n", argv[ai+1]);
                    exit(EXIT_FAILURE);
                }
                ai++;
                break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        if (ai != argc - 2) print_usage_and_exit();
        struct patch_set *ps = target_parse_text(argv[ai], class_name, observer, illuminant, scale, clip_level, low, high, spacing, column_layout, true);
        argyllio_ti3_write(argv[ai+1], true, ps, true);
        target_delete_patch_set(ps);

    } else if (strcmp(argv[1], "si-render") == 0) {

        /******************
          Spectral image render
         */
        bool userprovided_observer = false;
        bool used_observer = false;
        bool apply_wb = false;
        bool prophoto_enabled = false;
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        struct dcam_ssf *dcam_ssf = NULL;
        spectrum_t *ill_spec = NULL;
        int band_base = -1;
        int band_width = -1;
        m3x3 xyz2rgb_;
        m3x3 *xyz2rgb = NULL;
        enum ca_transform cat_;
        enum ca_transform *cat = NULL;
        double gamma = -1;
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminants on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                userprovided_observer = true;
                break;
            case 'c': dcam_ssf = jsonio_dcam_ssf_parse(argv[++ai], true); break;
            case 'g': gamma = atof(argv[++ai]); break;
            case 'i':
                free(ill_spec);
                ill_spec = NULL;
                parse_arg_illuminant(argv[++ai], observer, &used_observer, &ill_spec);
                if (ill_spec == NULL) {
                    elog("Error: illuminant with known spectrum is required.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'W': apply_wb = true; break;
            case 'b':
                ai++;
                band_base = atoi(argv[ai]);
                band_width = 10;
                if (strchr(argv[ai], ',') != NULL) {
                    char *p = strchr(argv[ai], ',') + 1;
                    band_width = atoi(p);
                }
                break;
            case 'a':
                ai++;
                if (strcasecmp(argv[ai], "bradford") == 0) {
                    cat_ = CAT_BRADFORD;
                } else if (strcasecmp(argv[ai], "cat02") == 0) {
                    cat_ = CAT_CAT02;
                } else {
                    elog("Error: CAT must be \"bradford\" or \"cat02\".\n");
                    exit(EXIT_FAILURE);
                }
                cat = &cat_;
                break;
            case 'P':
                xyz2rgb_ = dcp_xyz_D50_to_prophoto_rgb;
                xyz2rgb = &xyz2rgb_;
                if (gamma < 0) {
                    gamma = 461.0 / 256.0;
                }
                if (cat == NULL) {
                    cat_ = CAT_CAT02;
                    cat = &cat_;
                }
                prophoto_enabled = true;
                break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        init_observer_dependent_globals(observer);
        if (ill_spec == NULL) {
            ill_spec = spec_copy(spectraldb_illuminant(lsD50));
        }
        if (gamma <= 0) gamma = 1;
        if (userprovided_observer && dcam_ssf != NULL) {
            elog("Error: provide either camera SSF or observer, not both.\n");
            exit(EXIT_FAILURE);
        }
        if (ai != argc - 2) print_usage_and_exit();
        struct spb_image *spb = tifio_spb_load(argv[ai], band_base, band_width, true);
        struct rgbimg *img = spb_render(spb, ill_spec, dcam_ssf ? &dcam_ssf->ssf : observer, xyz2rgb, cat, 1/gamma, apply_wb);
        free(spb);
        tifio_rgbimg_save(argv[ai+1], img, prophoto_enabled ? COLORSPACE_PROPHOTO : COLORSPACE_NONE, true);
        free(img);
        free(ill_spec);
        jsonio_dcam_ssf_delete(dcam_ssf);

    } else if (strcmp(argv[1], "cspace") == 0) {

        /******************
          cspace - Transform coordinate
         */

        bool used_observer = false;
        bool custom_wp = false;
        const struct observer *observer = observer_get(OBSERVER_1931_2);
        m3x3 xyz2rgb = dcp_xyz_D50_to_prophoto_rgb;
        bool is_srgb_gamma = true;
        double gamma = 1.0;
        double unit_m = 100;
        v3 wp, in, out = {{0,0,0}};
        for (; ai < argc-1 && parse_flags && argv[ai][0] == '-'; ai++) {
            switch (argv[ai][1]) {
            case '\0': ai -= 2; parse_flags = false; break;
            case 'i':
                ai++;
                if (strcasecmp(argv[ai], "sRGB") == 0) {
                    wp = (v3){{0.95045471, 1.00000000, 1.08905029}};
                } else if (strcasecmp(argv[ai], "ProPhoto") == 0) {
                    wp = dcp_d50();
                } else if (strcasecmp(argv[ai], "AdobeRGB") == 0) {
                    wp = (v3){{0.95045471, 1.00000000, 1.08905029}};
                } else {
                    wp = parse_arg_illuminant(argv[ai], observer, &used_observer, NULL);
                }
                custom_wp = true;
                break;
            case 's':
                ai++;
                if (strcasecmp(argv[ai], "sRGB") == 0) {
                    xyz2rgb = xyz2rgb_srgb;
                } else if (strcasecmp(argv[ai], "ProPhoto") == 0) {
                    xyz2rgb = dcp_xyz_D50_to_prophoto_rgb;
                } else if (strcasecmp(argv[ai], "AdobeRGB") == 0) {
                    xyz2rgb = xyz2rgb_adobergb;
                } else {
                    elog("Unknown RGB color space \"%s\".\n", argv[ai]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'g':
                ai++;
                if (strcasecmp(argv[ai], "sRGB") == 0) {
                    is_srgb_gamma = true;
                } else if (strcasecmp(argv[ai], "ProPhoto") == 0) {
                    gamma = 1.8;
                } else if (strcasecmp(argv[ai], "AdobeRGB") == 0) {
                    gamma = 563.0/256.0;
                } else {
                    gamma = atof(argv[ai]);
                    if (gamma <= 0.0) {
                        elog("Bad RGB gamma \"%s\".\n", argv[ai]);
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'm': unit_m = atof(argv[++ai]); break;
            case 'o':
                observer = observer_get_byname(argv[++ai], true);
                if (used_observer) {
                    elog("Error: observer (-o) must be given before any illuminants on the command line.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                elog("unknown flag \"%s\".\n", argv[ai]);
                print_usage_and_exit();
                break;
            }
        }
        m3x3 rgb2xyz = m3x3_invert(xyz2rgb);
        init_observer_dependent_globals(observer);
        if (!custom_wp) {
            wp = glob.d50;
        }
        if (ai != argc - 4) print_usage_and_exit();
        const char *conv = argv[ai+0];
        if (strlen(conv) < 7) print_usage_and_exit();
        in.v[0] = atof(argv[ai+1]);
        in.v[1] = atof(argv[ai+2]);
        in.v[2] = atof(argv[ai+3]);
        char inspace[4];
        char outspace[4];
        memcpy(inspace, conv, 3);
        inspace[3] = '\0';
        memcpy(outspace, &conv[4], 3);
        outspace[3] = '\0';
        if (strcasecmp(inspace, "hsl") == 0 || strcasecmp(inspace, "hsv") == 0) {
            in.v[0] /= 360;
            in.v[1] /= unit_m;
            in.v[2] /= unit_m;
        }
        if (strcasecmp(inspace, "rgb") == 0 || strcasecmp(inspace, "xyz") == 0) {
            in.v[0] /= unit_m;
            in.v[1] /= unit_m;
            in.v[2] /= unit_m;
        }
        if (strcasecmp(conv, "lab2jch") == 0) {
            out = xyz2jch(lab2xyz(in, wp), wp);
        } else if (strcasecmp(conv, "jch2lab") == 0) {
            out = xyz2lab(jch2xyz(in, wp), wp);
        } else if (strcasecmp(conv, "xyz2lab") == 0) {
            out = xyz2lab(in, wp);
        } else if (strcasecmp(conv, "xyz2jch") == 0) {
            out = xyz2jch(in, wp);
        } else if (strcasecmp(conv, "xyz2xyY") == 0) {
            out = xyz2xyY(in);
        } else if (strcasecmp(conv, "xyz2uvY") == 0) {
            out = xyz2uvY(in);
        } else if (strcasecmp(conv, "lab2xyz") == 0) {
            out = lab2xyz(in, wp);
        } else if (strcasecmp(conv, "jch2xyz") == 0) {
            out = jch2xyz(in, wp);
        } else if (strcasecmp(conv, "xyY2xyz") == 0) {
            out = xyY2xyz(in);
        } else if (strcasecmp(conv, "uvY2xyz") == 0) {
            out = uvY2xyz(in);
        } else if (strcasecmp(conv, "rgb2hsl") == 0) {
            out = rgb2hsl(in);
        } else if (strcasecmp(conv, "rgb2hsv") == 0) {
            out = rgb2hsv(in);
        } else if (strcasecmp(conv, "hsl2rgb") == 0) {
            out = hsl2rgb(in);
        } else if (strcasecmp(conv, "hsl2jch") == 0) {
            out = xyz2jch(m3x3_mul_v3(rgb2xyz, apply_gamma(hsl2rgb(in), is_srgb_gamma, gamma, false)), wp);
        } else if (strcasecmp(conv, "hsl2jch") == 0) {
            out = xyz2jch(m3x3_mul_v3(rgb2xyz, apply_gamma(hsl2rgb(in), is_srgb_gamma, gamma, false)), wp);
        } else if (strcasecmp(conv, "hsv2rgb") == 0) {
            out = hsv2rgb(in);
        } else if (strcasecmp(conv, "hsv2hsl") == 0) {
            out = rgb2hsl(hsv2rgb(in));
        } else if (strcasecmp(conv, "hsl2hsv") == 0) {
            out = rgb2hsv(hsl2rgb(in));
        } else if (strcasecmp(conv, "rgb2lab") == 0) {
            out = xyz2lab(m3x3_mul_v3(rgb2xyz, apply_gamma(in, is_srgb_gamma, gamma, false)), wp);
        } else if (strcasecmp(conv, "rgb2jch") == 0) {
            out = xyz2jch(m3x3_mul_v3(rgb2xyz, apply_gamma(in, is_srgb_gamma, gamma, false)), wp);
        } else if (strcasecmp(conv, "rgb2xyz") == 0) {
            out = m3x3_mul_v3(rgb2xyz, apply_gamma(in, is_srgb_gamma, gamma, false));
        } else if (strcasecmp(conv, "lab2rgb") == 0) {
            out = apply_gamma(m3x3_mul_v3(xyz2rgb, lab2xyz(in, wp)), is_srgb_gamma, gamma, true);
        } else if (strcasecmp(conv, "jch2rgb") == 0) {
            out = apply_gamma(m3x3_mul_v3(xyz2rgb, jch2xyz(in, wp)), is_srgb_gamma, gamma, true);
        } else if (strcasecmp(conv, "xyz2rgb") == 0) {
            out = apply_gamma(m3x3_mul_v3(xyz2rgb, in), is_srgb_gamma, gamma, true);
        } else {
            elog("Unknown transformation \"%s\"\n", conv);
            exit(EXIT_FAILURE);
        }
        if (strcasecmp(outspace, "hsl") == 0 || strcasecmp(outspace, "hsv") == 0) {
            out.v[0] *= 360;
            out.v[1] *= unit_m;
            out.v[2] *= unit_m;
        }
        if (strcasecmp(outspace, "rgb") == 0 || strcasecmp(outspace, "xyz") == 0) {
            out.v[0] *= unit_m;
            out.v[1] *= unit_m;
            out.v[2] *= unit_m;
        }
        printf("%f %f %f\n", out.v[0], out.v[1], out.v[2]);
    } else {
        elog("Unknown command \"%s\". Run command without parameters to get usage documentation.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    return 0;
}
