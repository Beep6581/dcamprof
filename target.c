/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdbool.h>
#include <inttypes.h>
#include <strings.h>
#include <assert.h>

#include <target.h>
#include <elog.h>
#include <nmsimplex.h>
#include <colmath.h>
#include <dcamprof.h>
#include <jsonio.h>
#include <lut.h>
#include <interp.h>
#include <spectraldb.h>
#include <xyz2spec.h>

enum vtarget_type {
    VTARGET_MUNSELL = 0,
    VTARGET_MUNSELL_BRIGHT,
    VTARGET_KUOPIO_NATURAL,
    VTARGET_CC24,
    VTARGET_LOCUS,
    VTARGET_POINTER,
    VTARGET_SRGB,
    VTARGET_ADOBERGB,
    VTARGET_PROPHOTO,
    VTARGET_ILLUMINANT,
    VTARGET_WHITE,
    /*
    VTARGET_JOENSUU_BIRCH,
    VTARGET_JOENSUU_SPRUCE,
    VTARGET_JOENSUU_PINE,
    */
};

static const struct {
    const char *name;
    enum vtarget_type en;
    bool measured;
} targ_list[] = {
    { "cc24", VTARGET_CC24, true },
    { "munsell", VTARGET_MUNSELL, true },
    { "munsell-bright", VTARGET_MUNSELL_BRIGHT, true },
    { "kuopio-natural", VTARGET_KUOPIO_NATURAL, true },
    { "locus", VTARGET_LOCUS, false },
    { "pointer", VTARGET_POINTER, false },
    { "srgb", VTARGET_SRGB, false },
    { "adobergb", VTARGET_ADOBERGB, false },
    { "prophoto", VTARGET_PROPHOTO, false },
    { "illuminant", VTARGET_ILLUMINANT, true },
    { "white", VTARGET_WHITE, true },
/*
    { "joensuu-birch", VTARGET_JOENSUU_BIRCH, true },
    { "joensuu-spruce", VTARGET_JOENSUU_SPRUCE, true },
    { "joensuu-pine", VTARGET_JOENSUU_PINE, true },
*/
    { NULL, 0, false }
};

static const char *munsell_brightest_names[] = {
    "ABB5010", "ABB7006", "ABB7008", "ABB8004", "ABB9002", "ABG6010", "ABG7008", "ABG8004", "ABG8006", "ABG9002", "AGG5012", "AGG7010", "AGG8006", "AGG8008", "AGG9002", "AGG9004", "AGY8008", "AGY8010", "AGY8012", "AGY9002", "AGY9004", "AGY9006", "APB5012", "APB6010", "APB7008", "APB8004", "APB8006", "APB9002", "APP4012", "APP5010", "APP7006", "APP7008", "APP8004", "APP9002", "ARP6012", "ARP7008", "ARP7010", "ARP8004", "ARP8006", "ARP9002", "ARR5014", "ARR6010", "ARR6012", "ARR7008", "ARR8004", "ARR8006", "ARR9002", "AYR6014", "AYR6016", "AYR7008", "AYR7010", "AYR7012", "AYR8004", "AYR8006", "AYR9002", "AYY8014", "AYY8016", "AYY8506", "AYY8508", "AYY8510", "AYY8512", "AYY9002", "AYY9004", "BBB6010", "BBB7006", "BBB7008", "BBB8004", "BBB9002", "BBG6010", "BBG7006", "BBG7008", "BBG8004", "BBG9002", "BGG7008", "BGG7010", "BGG8004", "BGG8006", "BGG9002", "BGY7012", "BGY8006", "BGY8008", "BGY8010", "BGY9002", "BGY9004", "BPB5012", "BPB6010", "BPB7008", "BPB8004", "BPB8006", "BPB9002", "BPP4012", "BPP5010", "BPP7006", "BPP7008", "BPP8004", "BPP9002", "BRP6012", "BRP7008", "BRP7010", "BRP8004", "BRP8006", "BRP9002", "BRR5014", "BRR6012", "BRR7008", "BRR7010", "BRR8004", "BRR8006", "BRR9002", "BYR7010", "BYR7012", "BYR7014", "BYR8004", "BYR8006", "BYR8008", "BYR9002", "BYY8508", "BYY8510", "BYY8512", "BYY8514", "BYY9002", "BYY9004", "BYY9006", "CBB6010", "CBB7006", "CBB7008", "CBB8004", "CBB9002", "CBG5010", "CBG7006", "CBG7008", "CBG8004", "CBG9002", "CGG6010", "CGG7008", "CGG8004", "CGG8006", "CGG9002", "CGY6012", "CGY8006", "CGY8008", "CGY8010", "CGY9002", "CGY9004", "CPB4012", "CPB6010", "CPB7006", "CPB7008", "CPB8004", "CPB9002", "CPP4012", "CPP6010", "CPP7008", "CPP8004", "CPP8006", "CPP9002", "CRP5014", "CRP6010", "CRP6012", "CRP7008", "CRP8004", "CRP8006", "CRP9002", "CRR5014", "CRR5016", "CRR6012", "CRR7008", "CRR7010", "CRR8004", "CRR8006", "CRR9002", "CYR7012", "CYR7014", "CYR7016", "CYR8004", "CYR8006", "CYR8008", "CYR8010", "CYR9002", "CYY8508", "CYY8510", "CYY8512", "CYY9002", "CYY9004", "CYY9006", "DBB5012", "DBB6010", "DBB7008", "DBB8004", "DBB8006", "DBB9002", "DBG5010", "DBG7006", "DBG7008", "DBG8004", "DBG9002", "DGG6010", "DGG7008", "DGG8004", "DGG8006", "DGG9002", "DGY6012", "DGY7010", "DGY8006", "DGY8008", "DGY9002", "DGY9004", "DPB4012", "DPB6010", "DPB7006", "DPB7008", "DPB8004", "DPB9002", "DPP5012", "DPP6010", "DPP7008", "DPP8004", "DPP8006", "DPP9002", "DRP5014", "DRP6010", "DRP6012", "DRP7008", "DRP8004", "DRP8006", "DRP9002", "DRR5016", "DRR6012", "DRR6014", "DRR7008", "DRR7010", "DRR8004", "DRR8006", "DRR9002", "DYR8006", "DYR8008", "DYR8010", "DYR8012", "DYR8014", "DYR9002", "DYR9004", "DYY8508", "DYY8510", "DYY8512", "DYY9002", "DYY9004", "DYY9006", "BBB9001", "BBG9001", "BPB9001", "DBB9001", "DBG9001", "DPB9001", "BGY9001", "DGY9001", "BGG9001", "DGG9001", "BPP9001", "DPP9001", "BRP9001", "DRP9001", "BRR9001", "DRR9001", "BYR9001", "DYR9001", "BYY9001", "DYY9001", "ERR5014", "ERR6012", "EYR6012", "EYR6014", "FRR5014", "FRR6012", "FYR6014", "FYR7012", "GRR5014", "GRR6012", "GYR7012", "GYR7014", "GYR7016", "HRR5016", "HRR6012", "HRR6014", "HYR7012", "HYR7014", "EYY8012", "EYY8014", "EYY8016", "EGY8012", "EGG6012", "FYY8014", "FYY8512", "GYY8512", "HYY8512", "HGY6012", "FPB5012", "GPB4014", "GPB5012", "HPB4012", "ERP4012", "FRP6012", "GRP5014", "GRP6012", "HRP5014", "HRP6012", "NEUT950"
};

static const double pointer_x[36] = {
    0.508,0.538,0.588,0.637,0.659,0.634,0.594,0.557,0.523,0.482,0.444,0.409,0.371,0.332,0.288,0.242,0.202,0.177,
    0.151,0.151,0.162,0.157,0.159,0.142,0.141,0.129,0.138,0.145,0.145,0.161,0.188,0.252,0.324,0.393,0.451,0.508
};
static const double pointer_y[36] = {
    0.226,0.258,0.280,0.298,0.316,0.351,0.391,0.427,0.462,0.491,0.515,0.546,0.558,0.573,0.584,0.576,0.530,0.454,
    0.389,0.330,0.295,0.266,0.245,0.214,0.195,0.168,0.141,0.129,0.106,0.094,0.084,0.104,0.127,0.165,0.199,0.226
};

static void
dump_common_gamuts(const char report_dir[],
                   const struct observer *obs)
{
    const v3 wp = {{0.950456,1.0,1.088754}};
    v3 prim[3];
    FILE *stream = open_log_file(report_dir, "gmt-prophoto.dat");
    if (stream == NULL) {
        return;
    }
    prim[0] = xyz2lutspace(xyY2xyz((v3){{0.7347,0.2653,100}}));
    prim[1] = xyz2lutspace(xyY2xyz((v3){{0.0366,0.0001,100}}));
    prim[2] = xyz2lutspace(xyY2xyz((v3){{0.1596,0.8404,100}}));
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[1].v[1], prim[1].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[2].v[1], prim[2].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fclose(stream);

    stream = open_log_file(report_dir, "gmt-adobergb.dat");
    prim[0] = xyz2lutspace(xyY2xyz((v3){{0.6400,0.3300,100}}));
    prim[1] = xyz2lutspace(xyY2xyz((v3){{0.2100,0.7100,100}}));
    prim[2] = xyz2lutspace(xyY2xyz((v3){{0.1500,0.0600,100}}));
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[1].v[1], prim[1].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[2].v[1], prim[2].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fclose(stream);

    stream = open_log_file(report_dir, "gmt-srgb.dat");
    prim[0] = xyz2lutspace(xyY2xyz((v3){{0.6400,0.3300,100}}));
    prim[1] = xyz2lutspace(xyY2xyz((v3){{0.3000,0.6000,100}}));
    prim[2] = xyz2lutspace(xyY2xyz((v3){{0.1500,0.0600,100}}));
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[1].v[1], prim[1].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[2].v[1], prim[2].v[2]);
    fprintf(stream, "%f %f 0.0\n", prim[0].v[1], prim[0].v[2]);
    fclose(stream);

    stream = open_log_file(report_dir, "gmt-pointer.dat");
    for (int i = 0; i < 36; i++) {
        v3 xyz = xy2xyz_lightest(pointer_x[i], pointer_y[i], wp);
        v3 ltu = xyz2lutspace(xyz);
        uint32_t rgb = xyz2srgb_u24(xyz, wp);
        fprintf(stream, "%f %f 0.0 0x%06x\n", ltu.v[1], ltu.v[2], rgb);
    }
    {
        v3 xyz = xy2xyz_lightest(pointer_x[0], pointer_y[0], wp);
        v3 ltu = xyz2lutspace(xyz);
        uint32_t rgb = xyz2srgb_u24(xyz, wp);
        fprintf(stream, "%f %f 0.0 0x%06x\n", ltu.v[1], ltu.v[2], rgb);
    }
    fclose(stream);

    stream = open_log_file(report_dir, "gmt-locus.dat");
    for (int i = 0; i < spec_band_count(obs->cmf[0]); i++) {
        v3 xyz = {{ spec_at_idx(obs->cmf[0], i), spec_at_idx(obs->cmf[1], i), spec_at_idx(obs->cmf[2], i) }};
        v3 ltu = xyz2lutspace(xyz);
        uint32_t rgb = xyz2srgb_u24(xyz, wp);
        fprintf(stream, "%f %f 0.0 0x%06x\n", ltu.v[1], ltu.v[2], rgb);
    }
    fclose(stream);
}

static enum vtarget_type
vtarget_aton(const char name[],
             bool *is_measured_data)
{
    int i;
    for (i = 0; targ_list[i].name != NULL; i++) {
        if (strcasecmp(name, targ_list[i].name) == 0) {
            break;
        }
    }
    if (targ_list[i].name == NULL) {
        return -1;
    }
    if (is_measured_data != NULL) {
        *is_measured_data = targ_list[i].measured;
    }
    return targ_list[i].en;
}

const char *
target_measured_list(void)
{
    static char *list = NULL;
    if (list != NULL) {
        return list;
    }
    int len = 0;
    for (int i = 0; targ_list[i].name != NULL; i++) {
        if (targ_list[i].measured) len += strlen(targ_list[i].name) + 3;
    }
    char *str = malloc(len);
    str[0] = '\0';
    int offset = 0;
    for (int i = 0; targ_list[i].name != NULL; i++) {
        if (targ_list[i].measured && strcasecmp(targ_list[i].name, "illuminant") != 0 && strcasecmp(targ_list[i].name, "white") != 0) {
            offset += sprintf(&str[offset], "%s, ", targ_list[i].name);
        }
    }
    str[offset-2] = '\0';
    list = str;
    return str;
}

const char *
target_generated_list(void)
{
    static char *list = NULL;
    if (list != NULL) {
        return list;
    }
    int len = 0;
    for (int i = 0; targ_list[i].name != NULL; i++) {
        if (!targ_list[i].measured) len += strlen(targ_list[i].name) + 3;
    }
    char *str = malloc(len);
    str[0] = '\0';
    int offset = 0;
    for (int i = 0; targ_list[i].name != NULL; i++) {
        if (!targ_list[i].measured) offset += sprintf(&str[offset], "%s, ", targ_list[i].name);
    }
    str[offset-2] = '\0';
    list = str;
    return str;
}

struct patch_set *
target_make_virtual(const char name[],
                    const struct dcam_ssf *ssf,
                    const struct observer *obs,
                    const spectrum_t *cam_illuminant,
                    const spectrum_t *ref_illuminant,
                    bool prefer_cat,
                    double grid_spacing,
                    bool fill_grid,
                    bool skip_border,
                    struct patch_grid **patch_grid,
                    const char report_dir[])
{
    bool is_measured;
    enum vtarget_type vtarget = vtarget_aton(name, &is_measured);
    if ((int)vtarget == -1) {
        elog("Unknown target name \"%s\".\n", name);
        return NULL;
    }
    if (grid_spacing < 0.00001 && !is_measured) {
        elog("Grid spacing %g is too small.\n", grid_spacing);
        return NULL;
    }

    {
        FILE *stream = open_log_file(report_dir, "live-illuminant.dat");
        if (stream) {
            spec_print(cam_illuminant, stream);
            fclose(stream);
        }
    }

    dump_common_gamuts(report_dir, obs);
    // gnuplot: splot 'live-patches.dat' pt 7 ps 2 lc rgb variable
    FILE *stream1 = open_log_file(report_dir, "live-patches.dat");
    // gnuplot: plot 'live-spectra.dat' w l lc rgb variable
    FILE *stream2 = open_log_file(report_dir, "live-spectra.dat");

    if (spec_equal(ref_illuminant, cam_illuminant)) {
        prefer_cat = false;
    }
    struct patch_set *ps = NULL;
    int patch_count = 0;
    const char *db_name = NULL;
    const double *poly_x = NULL;
    const double *poly_y = NULL;
    spectrum_t **locus_border = NULL;
    double *locus_poly_x = NULL;
    double *locus_poly_y = NULL;
    int poly_count = 0;
    const spectrum_t *gen_illuminant = prefer_cat ? cam_illuminant : ref_illuminant;
    const int lop_count = 20;
    switch (vtarget) {
    default:
        // measured target default case
        patch_count = spectraldb_group_count(name);
        db_name = name;
        break;
    case VTARGET_MUNSELL_BRIGHT: {
        patch_count = sizeof(munsell_brightest_names)/sizeof(munsell_brightest_names[0]);
        db_name = "munsell";
        break;
    }
    case VTARGET_ILLUMINANT:
    case VTARGET_WHITE:
        patch_count = 1;
        break;
    case VTARGET_LOCUS: {
        const double peak[3] = { 0, 1, 0 };
        int start = spec_band_start(obs->cmf[0]);
        int end = spec_band_end(obs->cmf[0]);
        int i;
        for (i = start; i <= end-10; i += 5) {
            spectrum_t *s = spec_alloc(peak, i, i+10, 5);
            v3 xyz = spec2xyz_ill(s, obs, gen_illuminant, 1);
            free(s);
            if (xyz.v[1] >= 0.005) {
                break;
            }
        }
        start = i;
        for (i = end-10; i >= start; i -= 5) {
            spectrum_t *s = spec_alloc(peak, i, i+10, 5);
            v3 xyz = spec2xyz_ill(s, obs, gen_illuminant, 1);
            free(s);
            if (xyz.v[1] >= 0.005) {
                break;
            }
        }
        end = i;
        poly_count = (end - start) / 5 + 1 + lop_count;
        locus_border = malloc(poly_count * sizeof(*locus_border));
        locus_poly_x = malloc(poly_count * sizeof(*locus_poly_x));
        locus_poly_y = malloc(poly_count * sizeof(*locus_poly_y));
        poly_x = locus_poly_x;
        poly_y = locus_poly_y;
        for (int i = 0; i < poly_count; i++) {
            if (i < poly_count - lop_count) {
                locus_border[i] = spec_alloc(peak, start+i*5, start+i*5+10, 5);
            } else {
                // line of purples
                double freq[6] = { start, start+5, start+10,  end, end+5, end+10 };
                int k = i - (poly_count - lop_count);
                double a = (double)(k+1) / (lop_count+1);
                a = pow(a, 2);
                double dual_peak[6] = { 0,a,0, 0,1-a,0 };
                locus_border[i] = spec_alloc1(freq, dual_peak, 6);
            }
            v3 xyz = spec2xyz_ill(locus_border[i], obs, gen_illuminant, 1);
            v3 xyY = xyz2xyY(xyz);
            locus_poly_x[i] = xyY.v[0];
            locus_poly_y[i] = xyY.v[1];
        }
        break;
    }
    case VTARGET_POINTER: {
        poly_x = pointer_x;
        poly_y = pointer_y;
        poly_count = sizeof(pointer_x) / sizeof(pointer_x[0]);
        break;
    }
    case VTARGET_SRGB: {
        static const double cx[3] = { 0.6400, 0.3000, 0.1500 };
        static const double cy[3] = { 0.3300, 0.6000, 0.0600 };
        poly_x = cx;
        poly_y = cy;
        poly_count = sizeof(cx) / sizeof(cx[0]);
        break;
    }
    case VTARGET_ADOBERGB: {
        static const double cx[3] = { 0.6400, 0.2100, 0.1500 };
        static const double cy[3] = { 0.3300, 0.7100, 0.0600 };
        poly_x = cx;
        poly_y = cy;
        poly_count = sizeof(cx) / sizeof(cx[0]);
        break;
    }
    case VTARGET_PROPHOTO: {
        static const double cx[3] = { 0.7347, 0.0366, 0.1596 };
        static const double cy[3] = { 0.2653, 0.0001, 0.8404 };
        poly_x = cx;
        poly_y = cy;
        poly_count = sizeof(cx) / sizeof(cx[0]);
        break;
    }
    }

    const v3 ref_wp = v3_norm(spec2xyz(ref_illuminant, obs));
    const v3 cam_wp = v3_norm(spec2xyz(cam_illuminant, obs));
    const v3 gen_wp = prefer_cat ? cam_wp : ref_wp;
    if (is_measured) {
        ps = calloc(1, sizeof(*ps) + patch_count * sizeof(ps->patch[0]));
        for (int i = 0; i < patch_count; i++) {
            const spectrum_t *s;
            spectrum_t *free_s = NULL;
            bool emissive = false;
            char color_name[COLOR_NAME_MAX+20];
            sprintf(color_name, "%d", i+1);
            switch (vtarget) {
            case VTARGET_ILLUMINANT:
                s = cam_illuminant;
                emissive = true;
                sprintf(color_name, "ill%d", i+1);
                break;
            case VTARGET_WHITE:
                free_s = spec_copy(spectraldb_illuminant_byname("stdE", true));
                spec_scalar_multiply(free_s, 1.0 / spec_max(free_s));
                s = free_s;
                sprintf(color_name, "wht%d", i+1);
                break;
            case VTARGET_MUNSELL_BRIGHT:
                s = spectraldb_color(db_name, munsell_brightest_names[i]);
                strcpy(color_name, munsell_brightest_names[i]);
                break;
            default:
                s = spectraldb_color_byidx(db_name, i, color_name);
                break;
            }
            assert(s != NULL);
            v3 xyz;
            if (prefer_cat) {
                xyz = spec2xyz_ill(s, obs, cam_illuminant, 1);
                xyz = chromatic_adaptation_transform(CAT_CAT02, xyz, cam_wp, ref_wp);
            } else {
                xyz = spec2xyz_ill(s, obs, ref_illuminant, 1);
            }
            ps->patch[ps->patch_count].cam = ssf == NULL ? (v3){{0,0,0}} : spec2xyz_ill(s, &ssf->ssf, cam_illuminant, 1);
            ps->patch[ps->patch_count].xyz = xyz;
            ps->patch[ps->patch_count].spectrum = spec_copy(s);
            ps->patch[ps->patch_count].emissive = emissive;
            strncpy(ps->patch[ps->patch_count].str_id, color_name, sizeof(ps->patch[ps->patch_count].str_id)-1);
            if (stream1 != NULL) {
                uint32_t rgb = xyz2srgb_u24(xyz, ref_wp);
                v3 ltu = xyz2lutspace(xyz);
                fprintf(stream1, "%f %f %f 0x%06x %s\n", ltu.v[1], ltu.v[2], ltu.v[0], rgb, ps->patch[ps->patch_count].str_id);
                spec_print_gpcolor(s, rgb, stream2);
                fprintf(stream2, "\n");
            }
            ps->patch_count++;
            free(free_s);
        }
    }

    if (poly_x != NULL && poly_y != NULL) {

        double poly_u[poly_count];
        double poly_v[poly_count];
        for (int i = 0; i < poly_count; i++) {
            v3 xyY = {{ poly_x[i], poly_y[i], 1 }};
            v3 uvY = xyz2uvY(xyY2xyz(xyY));
            poly_u[i] = uvY.v[0];
            poly_v[i] = uvY.v[1];
        }

        double *edge_coords = NULL;
        int edge_coord_count = 0;
        if (!skip_border) { // fill polygon border
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < poly_count; i++) {
                    double u1 = poly_u[i];
                    double v1 = poly_v[i];
                    double u2 = poly_u[i==poly_count-1?0:i+1];
                    double v2 = poly_v[i==poly_count-1?0:i+1];
                    double d = sqrtf((u2-u1)*(u2-u1) + (v2-v1)*(v2-v1));
                    // if distance between points larger than spacing, split it in equal sections and set patch coordinates there
                    int count = (int)(d / grid_spacing);
                    if (vtarget == VTARGET_LOCUS) count = 0; // special case
                    if (k == 0) {
                        edge_coord_count += count + 1;
                    } else {
                        edge_coords[2*edge_coord_count+0] = u1;
                        edge_coords[2*edge_coord_count+1] = v1;
                        edge_coord_count++;
                        if (count > 0) {
                            double s = d / (count + 1);
                            for (int j = 1; j <= count; j++) {
                                double offset = j * s;
                                double up = u1 + (u2 - u1) * offset / d;
                                double vp = v1 + (v2 - v1) * offset / d;
                                edge_coords[2*edge_coord_count+0] = up;
                                edge_coords[2*edge_coord_count+1] = vp;
                                edge_coord_count++;
                            }
                        }
                    }
                }
                if (k == 0) {
                    edge_coords = malloc(edge_coord_count * 2 * sizeof(edge_coords[0]));
                    edge_coord_count = 0;
                }
            }
        }

        struct patch_grid *pg = NULL;
        struct fill_coords {
            double u, v;
            int u_, v_;
        } *fill_coords = NULL;
        int fill_coord_count = 0;
        if (fill_grid) { // fill inside polygon

            const int mul = 0x7FFFFFF;
            int u_step = (int)(grid_spacing * mul);
            const int u_max = (int)(0.65 * mul);
            int side = u_max / u_step + 1;
            fill_coords = malloc(side * side * sizeof(fill_coords[0]));
            v3 lut_wp = xyz2lutspace(ref_wp);
            int wp_offset_u = (int)(lut_wp.v[1] * mul) % u_step;
            int wp_offset_v = (int)(lut_wp.v[2] * mul) % u_step;
            int lo_u = side;
            int lo_v = side;
            int hi_u = 0;
            int hi_v = 0;
#pragma omp parallel for
            for (int u_ = 0; u_ < side; u_++) {
                double u = ((u_*u_step)+wp_offset_u) * (1.0/mul);
                for (int v_ = 0; v_ < side; v_++) {
                    double v = ((v_*u_step)+wp_offset_v) * (1.0/mul);
                    if (point_inside_polygon(poly_u, poly_v, poly_count, u, v)) {
#pragma omp critical
                        {
                            if (u_ > hi_u) hi_u = u_;
                            if (u_ < lo_u) lo_u = u_;
                            if (v_ > hi_v) hi_v = v_;
                            if (v_ < lo_v) lo_v = v_;
                            fill_coords[fill_coord_count].u_ = u_;
                            fill_coords[fill_coord_count].v_ = v_;
                            fill_coords[fill_coord_count].u = u;
                            fill_coords[fill_coord_count].v = v;
                            fill_coord_count++;
                        }
                    }
                }
            }
            fill_coords = realloc(fill_coords, fill_coord_count * sizeof(fill_coords[0]));
            if (patch_grid) {
                int rows = hi_v - lo_v + 1;
                int cols = hi_u - lo_u + 1;
                pg = calloc(1, sizeof(*pg) + sizeof(pg->p[0]) * rows * cols);
                pg->rows = rows;
                pg->cols = cols;
                for (int v_ = 0; v_ < pg->rows; v_++) {
                    for (int u_ = 0; u_ < pg->cols; u_++) {
                        double u = (((lo_u+u_)*u_step)+wp_offset_u) * (1.0/mul);
                        double v = (((lo_v+v_)*u_step)+wp_offset_v) * (1.0/mul);
                        pg->p[pg->cols * v_ + u_].u = u;
                        pg->p[pg->cols * v_ + u_].v = v;
                        pg->p[pg->cols * v_ + u_].p = NULL;
                    }
                }
                for (int i = 0; i < fill_coord_count; i++) {
                    fill_coords[i].u_ -= lo_u;
                    fill_coords[i].v_ -= lo_v;
                }
                *patch_grid = pg;
            }
        }

        int coord_count = fill_coord_count + edge_coord_count;
        if (fill_coord_count > 0) {
            elog("Generating spectra for %d patches in a grid...\n", coord_count);
        } else {
            elog("Generating spectra for %d patches along border...\n", coord_count);
        }
        ps = calloc(1, sizeof(*ps) + coord_count * sizeof(ps->patch[0]));
        int nosolution_count = 0;
        int next_print = coord_count / 10;
        int next_print_step = next_print;
        int complete_count = 0;
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < coord_count; i++) {
            double u, v;
            int u_ = -1, v_ = -1;
            if (i < edge_coord_count) {
                u = edge_coords[2*i+0];
                v = edge_coords[2*i+1];
            } else {
                int j = i - edge_coord_count;
                u = fill_coords[j].u;
                v = fill_coords[j].v;
                u_ = fill_coords[j].u_;
                v_ = fill_coords[j].v_;
            }
            v3 xyz = uv2xyz_lightest(u, v, gen_wp);
            spectrum_t *s;
            if (i < edge_coord_count && vtarget == VTARGET_LOCUS) {
                // locus special case
                s = locus_border[i];
            } else {
                s = xyz2spec(xyz, obs, gen_illuminant);
            }
#pragma omp atomic
            complete_count++;
            if (s == NULL) {
#pragma omp atomic
                nosolution_count++;
                elog("_");
            } else {
                double mx = spec_max(s);
                if (mx > 1.0) {
                    spec_scalar_multiply(s, 1.0 / mx);
                }
                elog(".");
                int pc;
                v3 cam = ssf == NULL ? (v3){{0,0,0}} : spec2xyz_ill(s, &ssf->ssf, cam_illuminant, 1);
                xyz = spec2xyz_ill(s, obs, gen_illuminant, 1);
                if (prefer_cat) {
                    xyz = chromatic_adaptation_transform(CAT_CAT02, xyz, cam_wp, ref_wp);
                }
#pragma omp critical
                {
                    if (stream1 != NULL) {
                        uint32_t rgb = xyz2srgb_u24(xyz, ref_wp);
                        v3 ltu = xyz2lutspace(xyz);
                        fprintf(stream1, "%f %f %f 0x%06x\n", ltu.v[1], ltu.v[2], ltu.v[0], rgb);
                        fflush(stream1);
                        spec_print_gpcolor(s, rgb, stream2);
                        fprintf(stream2, "\n");
                        fflush(stream2);
                    }
                    pc = ps->patch_count;
                    ps->patch_count++;
                    if (complete_count > next_print) {
                        next_print += next_print_step;
                        int percent = (100 * complete_count) / coord_count;
                        if (percent < 100) {
                            elog("%d%%", percent);
                        }
                    }
                }
                ps->patch[pc].cam = cam;
                ps->patch[pc].xyz = xyz;
                ps->patch[pc].spectrum = s;
                if (pg && u_ != -1) {
                    pg->p[pg->cols * v_ + u_].u = u;
                    pg->p[pg->cols * v_ + u_].v = v;
                    pg->p[pg->cols * v_ + u_].p = &ps->patch[pc];
                }
            }
        }
        elog("100%%\n");
        if (nosolution_count > 0) {
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "Could not find any spectrum for %d patches (%.2f%% of total amount). Those\n"
                "  patches were probably close or outside the observer's limits and have been\n"
                "  excluded from the set. Remaining patch count: %d.\n",
                nosolution_count, (double)100.0 * nosolution_count / coord_count, ps->patch_count);
        }
        free(edge_coords);
        free(fill_coords);
    }
    free(locus_poly_x);
    free(locus_poly_y);
    free(locus_border);

    //normalize_range(ps, 0, ps->patch_count);
    strcpy(ps->class_names[0], name);

    if (stream1 != NULL) {
        fclose(stream1);
        fclose(stream2);
    }
    return ps;
}

struct patch_set *
target_copy(const struct patch_set *ps)
{
    if (ps == NULL) {
        return NULL;
    }
    struct patch_set *ps1 = malloc(sizeof(*ps1) + ps->patch_count * sizeof(ps1->patch[0]));
    memcpy(ps1, ps, sizeof(*ps1) + ps->patch_count * sizeof(ps1->patch[0]));
    for (int i = 0; i < ps1->patch_count; i++) {
        if (ps->patch[i].spectrum != NULL) {
            ps1->patch[i].spectrum = spec_copy(ps->patch[i].spectrum);
        }
    }
    return ps1;
}

bool
target_isvalid_target_name(const char name[],
                           bool *is_measured_data)
{
    return (int)vtarget_aton(name, is_measured_data) != -1;
}

void
target_delete_patch_set(struct patch_set *ps)
{
    if (ps == NULL) {
        return;
    }
    for (int i = 0; i < ps->patch_count; i++) {
        free(ps->patch[i].spectrum);
    }
    free(ps);
}

struct patch_set *
target_merge(struct patch_set *inps[], // will be deleted
             int ps_count,
             double min_chroma_distance,
             bool error_is_fatal)
{
    int tot_count = 0;
    struct patch_set base;
    memset(&base, 0, sizeof(base));
    for (int i = 0; i < ps_count; i++) {
        tot_count += inps[i]->patch_count;
        int transtab[TARGET_MAX_CLASS_COUNT];
        bool need_translation = false;
        for (int j = 0; inps[i]->class_names[j][0] != '\0' && j < TARGET_MAX_CLASS_COUNT; j++) {
            transtab[j] = j;
            int k;
            bool found = false;
            for (k = 0; base.class_names[k][0] != '\0' && k < TARGET_MAX_CLASS_COUNT; k++) {
                if (strcasecmp(base.class_names[k], inps[i]->class_names[j]) == 0) {
                    transtab[j] = k;
                    found = true;
                    break;
                }
            }
            if (transtab[j] != j) {
                need_translation = true;
            } else if (!found) {
                if (base.class_names[TARGET_MAX_CLASS_COUNT-1][0] != '\0') {
                    elog("Too many classes, could not merge patch sets!\n");
                    if (error_is_fatal) exit(EXIT_FAILURE);
                    return NULL;
                }
                strcpy(base.class_names[k], inps[i]->class_names[j]);
                if (k != j) {
                    transtab[j] = k;
                    need_translation = true;
                }
            }
        }
        if (need_translation) {
            for (int j = 0; j < inps[i]->patch_count; j++) {
                inps[i]->patch[j].class = transtab[inps[i]->patch[j].class];
            }
        }
    }
    struct patch_set *ps = calloc(1, sizeof(*ps) + tot_count * sizeof(ps->patch[0]));
    *ps = base;
    ps->patch_count = tot_count;
    for (int i = 0, pc = 0; i < ps_count; i++) {
        for (int j = 0; j < inps[i]->patch_count; j++) {
            ps->patch[pc++] = inps[i]->patch[j];
            inps[i]->patch[j].spectrum = NULL;
        }
        target_delete_patch_set(inps[i]);
    }
    int illum_class = -1;
    int white_class = -1;
    for (int i = 0; ps->class_names[i][0] != '\0' && i < TARGET_MAX_CLASS_COUNT; i++) {
        if (strcmp(ps->class_names[i], "illuminant") == 0) {
            illum_class = i;
        }
        if (strcmp(ps->class_names[i], "white") == 0) {
            white_class = i;
        }
    }
    if (min_chroma_distance > 0.0) {
        int skip_count = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            if (ps->patch[i].class == illum_class || ps->patch[i].class == white_class || ps->patch[i].force_keep) {
                continue;
            }
            for (int j = 0; j < ps->patch_count; j++) {
                if (ps->patch[j].class < ps->patch[i].class) {
                    v3 ltu1, ltu2;
                    ltu1 = xyz2lutspace(ps->patch[i].xyz);
                    ltu2 = xyz2lutspace(ps->patch[j].xyz);
                    ltu1.v[0] = ltu2.v[0];
                    double uv_de = euclidean_dist(ltu1, ltu2);
                    if (uv_de < min_chroma_distance) {
                        ps->patch[i].class = -1;
                        skip_count++;
                    }
                }
            }
        }
        if (skip_count > 0) {
            ps->patch_count = 0;
            for (int i = 0; i < tot_count; i++) {
                if (ps->patch[i].class >= 0) {
                    ps->patch[ps->patch_count++] = ps->patch[i];
                }
            }
        }
    }
    return ps;
}

void
target_set_name(struct patch_set *ps,
                const char class_name[])
{
    strncpy(ps->class_names[0], class_name, TARGET_CLASS_NAME_SIZE - 1);
    if (ps->class_names[1][0] == '\0') {
        return;
    }
    for (int i = 1; i < TARGET_MAX_CLASS_COUNT; i++) {
        memset(ps->class_names[i], 0, sizeof(ps->class_names[i]));
    }
    for (int i = 0; i < ps->patch_count; i++) {
        ps->patch[i].class = 0;
    }
}

int
target_class_or_patch_to_idx(const struct patch_set *ps,
                             const char name[],
                             bool *is_class)
{
    if (ps == NULL || name == NULL) {
        return -1;
    }
    for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
        if (strcasecmp(ps->class_names[i], name) == 0) {
            if (is_class != NULL) *is_class = true;
            return i;
        }
    }
    for (int i = 0; i < ps->patch_count; i++) {
        if (strcasecmp(ps->patch[i].str_id, name) == 0) {
            if (is_class != NULL) *is_class = false;
            return i;
        }
    }
    return -1;
}

int
target_name2idx(const struct patch_set *ps,
                const char name[],
                const char class_name[])
{
    if (ps == NULL || name == NULL) {
        return -1;
    }
    for (int i = 0; i < ps->patch_count; i++) {
        if (strcasecmp(ps->patch[i].str_id, name) == 0) {
            if (class_name == NULL || class_name[0] == '\0' || strcasecmp(ps->class_names[ps->patch[i].class], class_name) == 0) {
                return i;
            }
        }
    }
    return -1;
}

const char *
target_patch_name(char buf[PATCH_STRID_SIZE],
                  const struct patch_set *ps,
                  int idx)
{
    const char *id;
    if (ps->patch[idx].str_id[0] != '\0') {
        id = ps->patch[idx].str_id;
    } else {
        sprintf(buf, "%d", idx + 1);
        id = buf;
    }
    return id;
}


struct patch_set *
target_parse_text(const char filename[],
                  const char class_name[],
                  const struct observer *obs,
                  const spectrum_t *illuminant,
                  double scale,
                  double clip_level,
                  int low,
                  int high,
                  int spacing,
                  bool column_layout,
                  bool error_is_fatal)
{
    struct patch_set *ps = NULL;
    double *buf = NULL;
    char *str = fileio_file2string(filename);
    if (str == NULL) {
        goto fail;
    }
    if (strchr(str, '.') == NULL && strchr(str, ',') != NULL) {
        for (char *p = str; *p != '\0'; p++) {
            if (*p == ',') *p = '.';
        }
    }
    size_t bufsize = 4096;
    size_t count = 0;
    buf = (double *)malloc(bufsize * sizeof(*buf));
    char *p = str, *endp;
    for (;;) {
        double val = 0;
        do {
            endp = NULL;
            val = strtod(p, &endp);
            if (p == endp) {
                p++;
            } else {
                break;
            }
        } while (*p != '\0');
        if (*p == '\0') {
            break;
        }
        p = endp;
        buf[count++] = val;
        if (count == bufsize) {
            bufsize *= 2;
            buf = (double *)realloc(buf, bufsize * sizeof(*buf));
        }
    }
    buf = (double *)realloc(buf, count * sizeof(*buf));
    free(str);
    str = NULL;
    const int band_count = (high - low) / spacing + 1;
    if (count % band_count != 0) {
        elog("Found %d values in the file, which is not divisable by the band count %d\n", (int)count, band_count);
        goto fail;
    }
    const int patch_count = count / band_count;
    ps = calloc(1, sizeof(*ps) + patch_count * sizeof(ps->patch[0]));
    ps->patch_count = patch_count;
    for (int i = 0; i < patch_count; i++) {
        double amp[band_count];
        if (column_layout) {
            for (int j = 0; j < band_count; j++) {
                amp[j] = buf[i + j * patch_count];
            }
        } else {
            for (int j = 0; j < band_count; j++) {
                amp[j] = buf[i * band_count + j];
            }
        }
        for (int j = 0; j < band_count; j++) {
            amp[j] *= scale;
            if (clip_level > 0 && amp[j] > clip_level) {
                amp[j] = clip_level;
            }
        }
        spectrum_t *s = spec_alloc(amp, low, high, spacing);
        v3 xyz = spec2xyz_ill(s, obs, illuminant, 1);
        ps->patch[i].xyz = xyz;
        ps->patch[i].spectrum = s;
    }
    if (class_name != NULL) {
        strncpy(ps->class_names[0], class_name, TARGET_CLASS_NAME_SIZE - 1);
    }
    free(buf);
    return ps;
fail:
    free(str);
    free(buf);
    target_delete_patch_set(ps);
    if (error_is_fatal) {
        exit(EXIT_FAILURE);
    }
    return NULL;
}

void
target_print(const char prefix[],
             const struct patch_set *ps,
             const v3 wp,
             const struct observer *obs,
             const char report_dir[])
{
    if (ps == NULL || ps->patch_count == 0) {
        return;
    }
    char filename[256];
    sprintf(filename, "%s-xyz.dat", prefix);
    FILE *stream = open_log_file(report_dir, filename);
    if (stream == NULL) {
        return;
    }
    bool print_spectra = (ps->patch[0].spectrum != NULL);
    sprintf(filename, "%s-spectra.dat", prefix);
    FILE *stream_s = print_spectra ? open_log_file(report_dir, filename) : NULL;
    FILE *cstream[TARGET_MAX_CLASS_COUNT];
    FILE *cstream_s[TARGET_MAX_CLASS_COUNT];
    memset(cstream, 0, sizeof(cstream));
    memset(cstream_s, 0, sizeof(cstream_s));
    for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
        if (ps->class_names[i][0] != '\0') {
            char name[TARGET_CLASS_NAME_SIZE+256];
            sprintf(name, "%s-xyz-%s.dat", prefix, ps->class_names[i]);
            cstream[i] = open_log_file(report_dir, name);
            if (print_spectra) {
                sprintf(name, "%s-spectra-%s.dat", prefix, ps->class_names[i]);
                cstream_s[i] = open_log_file(report_dir, name);
            }
        }
    }
    for (int i = 0; i < ps->patch_count; i++) {
        v3 xyz = ps->patch[i].xyz;
        char line[128+sizeof(ps->patch[i].str_id)];
        uint32_t rgb = xyz2srgb_u24(xyz, wp);
        v3 ltu = xyz2lutspace(xyz);
        const spectrum_t *s = ps->patch[i].spectrum;
        char id_[PATCH_STRID_SIZE];
        const char *id;
        if (ps->patch[i].str_id[0] != '\0') {
            id = ps->patch[i].str_id;
        } else {
            sprintf(id_, "%d", i + 1);
            id = id_;
        }
        sprintf(line, "%f %f %f 0x%06x %s\n", ltu.v[1], ltu.v[2], ltu.v[0], rgb, id);
        fprintf(stream, "%s", line);
        if (cstream[ps->patch[i].class] != NULL) {
            fprintf(cstream[ps->patch[i].class], "%s", line);
        }
        if (s != NULL && print_spectra) {
            spec_print_gpcolor(s, rgb, stream_s);
            fprintf(stream_s, "\n");
            if (cstream_s[ps->patch[i].class] != NULL) {
                spec_print_gpcolor(s, rgb, cstream_s[ps->patch[i].class]);
                fprintf(cstream_s[ps->patch[i].class], "\n");
            }
        }
    }
    fclose(stream);
    if (stream_s != NULL) {
        fclose(stream_s);
    }
    for (int i = 0; i < TARGET_MAX_CLASS_COUNT; i++) {
        if (cstream[i] != NULL) {
            fclose(cstream[i]);
        }
        if (cstream_s[i] != NULL) {
            fclose(cstream_s[i]);
        }
    }
    dump_common_gamuts(report_dir, obs);
}

#if 0
#include <observers.h>

static void __attribute__((constructor))
precalc_grid(void)
{
    struct patch_grid *pg;
    struct patch_set *ps = target_make_virtual("locus", NULL, observer_get(OBSERVER_2006_10), spectraldb_illuminant(lsD50), 0.01, true, true, &pg, NULL);
    printf("#include <spectrum.h>\n");
    printf("#include <dcamprof.h>\n");
    printf("#include <target.h>\n");
    printf("struct patch_set *allocate_precalc_patch_grid(struct patch_grid **pg_) {\n");
    printf("struct patch_grid *pg; struct patch_set *ps;\n");
    printf("pg = calloc(1, sizeof(*pg) + %d * %d * sizeof(pg->p[0]));\n", pg->rows, pg->cols);
    printf("pg->rows = %d; pg->cols = %d;\n", pg->rows, pg->cols);
    printf("ps = calloc(1, sizeof(*ps) + %zd * sizeof(ps->patch[0]));\n", ps->patch_count);
    printf("ps->patch_count = %zd;\n", ps->patch_count);
    for (int i = 0; i < ps->patch_count; i++) {
        printf("ps->patch[%d].xyz = (v3){{%g,%g,%g}};\n", i, ps->patch[i].xyz.v[0], ps->patch[i].xyz.v[1], ps->patch[i].xyz.v[2]);
        printf("ps->patch[%d].spectrum = malloc(sizeof(spectrum_t)* %d * sizeof(struct spectrum_band_));\n", i, ps->patch[i].spectrum->band_count);
        printf("ps->patch[%d].spectrum->band_count = %d;\n", i, ps->patch[i].spectrum->band_count);
        for (int j = 0; j < ps->patch[i].spectrum->band_count; j++) {
            printf("ps->patch[%d].spectrum->v[%d].f = %g; ps->patch[%d].spectrum->v[%d].e = %g;\n",
                   i, j, ps->patch[i].spectrum->v[j].f, i, j, ps->patch[i].spectrum->v[j].e);
        }
    }
    for (int i = 0; i < pg->rows; i++) {
        for (int j = 0; j < pg->cols; j++) {
            int idx = pg->cols * i + j;
            printf("pg->p[%d].u = %g; pg->p[%d].v = %g; ", idx, pg->p[idx].u, idx, pg->p[idx].v);
            if (pg->p[idx].p == NULL) {
                printf("pg->p[%d].p = NULL;\n", idx);
            } else {
                size_t pidx = ((uintptr_t)pg->p[idx].p - (uintptr_t)ps->patch) / sizeof(ps->patch[0]);
                printf("pg->p[%d].p = &ps->patch[%zd];\n", idx, pidx);
            }
        }
    }
    printf("*pg_ = pg;\nreturn ps;\n}\n");
    exit(1);
}

#endif
