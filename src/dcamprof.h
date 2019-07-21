/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef DCAMPROF_H_
#define DCAMPROF_H_

#include <stdio.h>
#include <stdbool.h>

#include <colmath.h>

struct cam2xyz_patch {
    spectrum_t *spectrum;
#define PATCH_STRID_SIZE 32
    char str_id[PATCH_STRID_SIZE];
    bool emissive;
    bool force_keep;
    bool use_mtx_de_w;
    int class;
    double lut_de_w[6];
    double mtx_de_w[6];
    double mtx_w;
    v3 cam;
    v3 xyz;
    v3 tmp;
};

struct patch_set {
    struct {
        int row_count;
        double patch_aspect_ratio;
        double patch_inner_ratio;
    } layout;
    enum exif_lightsource illuminant;
    v3 ill_wp; // set to D50 if illuminant is unknown
    int patch_count;
#define TARGET_MAX_CLASS_COUNT 128
#define TARGET_CLASS_NAME_SIZE 24
    char class_names[TARGET_MAX_CLASS_COUNT][TARGET_CLASS_NAME_SIZE];
    struct cam2xyz_patch patch[];
};

struct patch_error {
    int idx;
    v3 cxyz;
    double errors[4];
};

struct patch_error_report {
    enum {
        LOOKUP_MATRIX,
        LOOKUP_LUT_NATIVE,
        LOOKUP_DCP_HSM,
        LOOKUP_ICC_LUT
    } lookup_type;
    struct {
        double avg;
        double median;
        double p90;
        double max;
        int worst[5];
        int best[5];
    } errors[4];
    int patch_count;
    const struct patch_set *ps;
    v3 ps_wp;
    struct patch_error e[];
};

enum tc_type {
    TC_NEUTRAL,
    TC_STANDARD
};

enum gc_type {
    GC_FROM_FILE = 0,
    GC_NONE,
    GC_SRGB,
    GC_ADOBERGB,
    GC_SRGB_STRONG,
    GC_ADOBERGB_STRONG
};

enum lookop_type {
    LOOKOP_STRETCH,
    LOOKOP_ADDHUE,
    LOOKOP_ADDCHROMA,
    LOOKOP_ADDLIGHTNESS,
    LOOKOP_SCALECHROMA,
    LOOKOP_SCALELIGHTNESS,
    LOOKOP_SETTEMPERATURE,
    LOOKOP_CURVES
};

enum lookop_vs_x {
    LOOKOP_VS_X_LIGHTNESS,
    LOOKOP_VS_X_CHROMA,
    LOOKOP_VS_X_HUE,
    LOOKOP_VS_X_HSL_LIGHTNESS,
    LOOKOP_VS_X_HSL_SATURATION,
    LOOKOP_VS_X_HSL_HUE,
    LOOKOP_VS_X_HSV_VALUE,
    LOOKOP_VS_X_HSV_SATURATION,
    LOOKOP_VS_X_HSV_HUE
};

enum lookop_curve_type {
    LOOKOP_CURVE_SPLINE,
    LOOKOP_CURVE_ROUNDEDSTEP,
    LOOKOP_CURVE_LINEAR
};

struct lookop_curve {
    enum lookop_curve_type type;
    int handle_count;
    double handles[];
};

struct lookop_vs {
    enum lookop_vs_x x_type;
    double xrange[2];
    double yscale;
    struct lookop_curve *curve;
};

struct lookop {
    enum lookop_type type;
    union {
        struct {
            double value;
        } gen;
        struct {
            int dim_count;
            enum lookop_vs_x dims[3];
            double ranges[3][2];
            struct lookop_curve *curves[3];
        } stretch;
        struct {
            bool is_srgb_gamma;
            bool keep_lightness;
            double gamma;
            struct lookop_curve *curves[3];
        } curves;
        struct {
            double xy[2];
        } settemp;
    } p;
    bool blend_invert;
    bool blend_rgb;
    int vs_count;
    struct lookop_vs *vs;
};

struct gamut_spec {
    bool locus_enabled;
    bool xyz2rgb_enabled;
    m3x3 xyz2rgb;
    double chroma_scale;
};

struct gamut_compress_state_t_;
typedef struct gamut_compress_state_t_ gamut_compress_state_t;

struct gamut_compression {
    double hi_rgb_lim;
    double lo_rgb_lim;
    struct gamut_spec inner;
    struct gamut_spec dest;
    struct gamut_spec clip;
    struct {
        bool enabled;
        double lo_sat_blend_lim;
        double hi_sat_blend_lim;
        struct {
            double min_lift;
            double linear_limit;
            double cutoff_limit;
            double output_limit;
        } top;
        struct {
            double max_lift;
            double linear_limit;
            double cutoff_limit;
            double output_limit;
        } bottom;
    } hsv;
    gamut_compress_state_t *state;
};

struct tone_rep_op_config {
    double chroma_scaling;
    struct {
        double keep_factor;
        double low_chroma;
        double high_chroma;
    } curve;
    struct {
        double adjust_factor;
        double low_chroma;
        double high_chroma;
    } saturated;
    struct {
        double adjust_factor;
        double adjust_factor_high_chroma;
        double low_lightness;
        double high_lightness;
        double low_chroma;
        double high_chroma;
    } shadows;
    struct {
        struct lookop_curve *keep_factor;
        struct lookop_curve *low_satscale;
        struct lookop_curve *high_satscale;
    } rolloff;
    struct gamut_compression gamut;
    bool gamut_enabled;
    int lookop_count;
    struct lookop *lookops;
};

struct patch_adjustment {
    char class[TARGET_CLASS_NAME_SIZE];
    char name[8];
    double scale_rgb[3];
    double adjust_jch[3];
};

struct target_adjustment_config {
    struct {
        double scale_chroma;
        double scale_lightness;
        struct lookop_curve *scale_chroma_per_hue;
        struct lookop_curve *scale_lightness_per_hue;
        struct lookop_curve *add_hue_per_hue;
    } global;
    int patch_count;
    struct patch_adjustment patch[];
};

struct dcam_ssf {
    char *camera_name;
    struct observer ssf;
};

struct chromalut_t_;
typedef struct chromalut_t_ chromalut_t;

struct dcam_profile {
    enum exif_lightsource illuminant_name;
    spectrum_t *illuminant_spec;
    v3 illuminant_wp;
    m3x3 color_matrix;
    m3x3 forward_matrix;
    m3x3 lut_matrix;
    v3 forward_matrix_wb;
    char *camera_name;
    chromalut_t *lut;
};

struct test_chart_layout {
    int row_count;
    bool has_col_offset;
    double row_height;
    double col_width;
};

struct test_chart_patch {
    char str_id[PATCH_STRID_SIZE];
    bool is_black, is_white, is_gray;
    double r, c;
    int idx;
};

struct test_chart_spec {
    struct test_chart_layout layout;
    int col_count;
    int white_count;
    int black_count;
    int gray_count;
    int patch_count;
    bool gray_specified;
    struct test_chart_patch *white;
    struct test_chart_patch *black;
    struct test_chart_patch *gray;
    struct test_chart_patch data[];
};

static inline void
test_chart_rc(const struct test_chart_layout *layout,
              int idx,
              double *r, double *c)
{
    int row = idx % layout->row_count;
    int col = idx / layout->row_count;
    *r = (double)row + (layout->has_col_offset ? ((col&1)==1 ? 0.5 : 0) : 0);
    *c = (double)col;
    *r *= layout->row_height;
    *c *= layout->col_width;
}

static inline double
test_chart_dist(const struct test_chart_layout *layout,
                int idx1, int idx2)
{
    double r1, r2, c1, c2;
    test_chart_rc(layout, idx1, &r1, &c1);
    test_chart_rc(layout, idx2, &r2, &c2);
    double r = (r1 - r2);
    double c = (c1 - c2);
    double dist = sqrt(r*r + c*c);
    return dist;
}

void
dcam_profile_delete(struct dcam_profile *prof);

FILE *
open_log_file(const char dir[],
              const char name[]);

void
render_pm_error_image(const char report_dir[],
                      const char filename[],
                      const v3 ps_wp,
                      const struct patch_error_report *pe,
                      bool reverse_order);

#endif
