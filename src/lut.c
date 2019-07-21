/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include <lcms2.h>

#include <target.h>
#include <interp.h>
#include <lut.h>
#include <tps.h>
#include <elog.h>
#include <nmsimplex.h>
#include <bisection.h>
#include <dngref.h>
#include <gamut.h>

struct chromalut_t_ {
    tps_t *tps[3];
    double regularization[3];
    double compress_factor;
    v3 *cp[3];
    int cp_count;
    m3x3 fm;
    m3x3 ifm;
};

v3
xyz2lutspace(v3 xyz)
{
    double x = xyz.v[0], y = xyz.v[1], z = xyz.v[2];
    // u' v' and L*
    double up = 4*x / (x + 15*y + 3*z);
    double vp = 9*y / (x + 15*y + 3*z);
    double L = 116*lab_ft_forward(y) - 16;
    if (!isfinite(up)) up = 0;
    if (!isfinite(vp)) vp = 0;

    return (v3){{ L*0.01, up, vp }};
}

v3
lutspace2xyz(v3 lutspace)
{
    double L = lutspace.v[0]*100.0, up = lutspace.v[1], vp = lutspace.v[2];
    double y = (L + 16)/116;
    y = lab_ft_inverse(y);
    double x = y*9*up / (4*vp);
    double z = y * (12 - 3*up - 20*vp) / (4*vp);
    if (!isfinite(x)) x = 0;
    if (!isfinite(z)) z = 0;

    return (v3){{ x, y, z }};
}

static double
tps_interpolate_(tps_t *tps,
                 double x,
                 double z)
{
    if (tps == NULL) {
        return 0;
    }
    return tps_interpolate(tps, x, z);
}

struct extreme_compress_hsv_ge_arg {
    v3 hsv;
    const m3x3 *fm;
    v3 last_valid_hsv;
};

static double
extreme_compress_hsv_ge(double x, void *arg)
{
    struct extreme_compress_hsv_ge_arg *a = (struct extreme_compress_hsv_ge_arg *)arg;
    v3 hsv = a->hsv;
    hsv.v[1] = x;
    v3 xyz = m3x3_mul_v3(*a->fm, hsv2rgb(hsv));
    if (!gamut_isvalid(xyz)) {
        return 1 + x;
    }
    a->last_valid_hsv = hsv;
    return 1 - x;
}

static v3
extreme_compress(const m3x3 fm,
                 const m3x3 ifm,
                 const double compress_factor,
                 const v3 orig_xyz)
{
    if (compress_factor <= 0) {
        return orig_xyz;
    }
    /* Why compressing in HSV-Saturation? More tonality is kept than in HSL-Saturation
       (as we go from saturated to white instead of gray)
       In the extreme range keeping nice gradients with visible tonality is more
       important than correct color, as we don't have correct color anyway in this range.

       Scaling raw RGB channels with constant HSV-hue should be reasonably hue-stable.
    */
    v3 cam = m3x3_mul_v3(ifm, orig_xyz);
    cam = gamut_acr_rgb_clip(cam);
    v3 hsv = rgb2hsv(cam);

    double hsv_chroma_transition;
    {
        // principle: longer out-of-gamut range => longer transition
        struct extreme_compress_hsv_ge_arg arg = {
            .hsv = hsv,
            .fm = &fm
        };
        arg.hsv.v[1] = 1.0;
        interval_halving_min(extreme_compress_hsv_ge, &arg, 0, 1, 1e-6, 100);
        const v3 hsv_edge = arg.last_valid_hsv;
        hsv_chroma_transition = compress_factor * hsv_edge.v[1];
    }

    if (gamut_isvalid(orig_xyz)) {
        v3 test_hsv = hsv;
        test_hsv.v[1] /= hsv_chroma_transition;
        if (test_hsv.v[1] > 1) test_hsv.v[1] = 1;
        if (gamut_isvalid(m3x3_mul_v3(fm, hsv2rgb(test_hsv)))) {
            v3 xyz = m3x3_mul_v3(fm, cam);
            return xyz;
        }
    }

    struct extreme_compress_hsv_ge_arg arg = {
        .hsv = hsv,
        .fm = &fm
    };
    interval_halving_min(extreme_compress_hsv_ge, &arg, 0, 1, 1e-6, 100);
    const v3 hsv_edge = arg.last_valid_hsv;
    v3 hsv_desat = hsv_edge;
    hsv_desat.v[1] *= hsv_chroma_transition;
    double dist = (hsv.v[1] - hsv_desat.v[1]) / (1 - hsv_desat.v[1]);

    // gradually increasing
    double diff = (1 - hsv_desat.v[1]) / (hsv_edge.v[1] - hsv_desat.v[1]);
    double ndist = dist * (1 + (diff - 1) * (1 - pow(dist, 0.2)));
    dist = ndist;

    hsv.v[1] = hsv_desat.v[1] + dist * (hsv_edge.v[1] - hsv_desat.v[1]);
    v3 xyz = m3x3_mul_v3(fm, hsv2rgb(hsv));
    return xyz;
}

static v3
chromalut_lookup_(tps_t *tps[3],
                  const m3x3 fm,
                  const m3x3 ifm,
                  const double compress_factor,
                  v3 xyz,
                  bool do_clip)
{
    xyz = extreme_compress(fm, ifm, compress_factor, xyz);
    int zero_count = 0;
    for (int i = 0; i < 3; i++) {
        if (xyz.v[i] < 0) {
            xyz.v[i] = 0;
            zero_count++;
        }
    }
    if (zero_count == 3) {
        return (v3){{ 0,0,0 }};
    }
    v3 lutspace = xyz2lutspace(xyz);
    v3 lso = lutspace;
    lso.v[0] *= 1.0 + 10.0 * tps_interpolate_(tps[0], lutspace.v[1], lutspace.v[2]);
    if (do_clip) {
        if (lso.v[0] < 0) lso.v[0] = 0;
    }
    lso.v[1] += tps_interpolate_(tps[1], lutspace.v[1], lutspace.v[2]);
    lso.v[2] += tps_interpolate_(tps[2], lutspace.v[1], lutspace.v[2]);
    v3 out = lutspace2xyz(lso);
    if (do_clip) {
        for (int i = 0; i < 3; i++) {
            if (out.v[i] < 0) out.v[i] = 0;
        }
    }
    if (!v3_isfinite(out)) {
        v3_print("xyz  input", xyz);
        v3_print("xyz output", out);
        v3_print("lutspace input", lutspace);
        v3_print("lutspace output", lso);
        elog("Bug: non-finite number in lookup.\n");
        abort();
    }
    return out;
}

static p2d *
edge_coords(const p2d v[],
            int nvert,
            double spacing,
            int *count)
{
    // this edge getter will prioritize spacing, it will thus cut corners
    *count = 0;
    if (nvert == 0) {
        return NULL;
    }
    p2d *ec = NULL;
    for (int k = 0; k < 2; k++) {
        int j = 0;
        int i = 1;
        double xp = v[0].x;
        double yp = v[0].y;
        while (i <= nvert) {
            if (k == 1) {
                ec[j].x = xp;
                ec[j].y = yp;
            }
            double x = xp;
            double y = yp;
            double x1 = i == nvert ? v[0].x : v[i].x;
            double y1 = i == nvert ? v[0].y : v[i].y;
            double d = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
            if (d < spacing) {
                i++;
                continue;
            }
            double x0 = v[i-1].x;
            double y0 = v[i-1].y;

            /*
              sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp)) == spacing;
              xp = x0 + (x1 - x0) * p;
              yp = y0 + (y1 - y0) * p;
            */

            const double s = spacing;
            double p = (sqrt((-x0*x0+2*x*x0-x*x+s*s)*y1*y1+(((2*x0-2*x)*x1-2*x*x0+2*x*x-2*s*s)*y0+((2*x-2*x0)*x1+2*x0*x0-2*x*x0)*y)*y1+(-x1*x1+2*x*x1-x*x+s*s)*y0*y0+(2*x1*x1+(-2*x0-2*x)*x1+2*x*x0)*y*y0+(-x1*x1+2*x0*x1-x0*x0)*y*y+s*s*x1*x1-2*s*s*x0*x1+s*s*x0*x0)+(y-y0)*y1+y0*y0-y*y0+(x-x0)*x1+x0*x0-x*x0)/(y1*y1-2*y0*y1+y0*y0+x1*x1-2*x0*x1+x0*x0);

            xp = x0 + (x1 - x0) * p;
            yp = y0 + (y1 - y0) * p;

            j++;
            //((x-(x0+(x1-x0)*p))*(x-(x0+(x1-x0)*p))+(y-(y0+(y1-y0)*p))*(y-(y0+(y1-y0)*p)))^0.5=s,X=,Y=
        }
        if (k == 0) {
            j++;
            *count = j;
            ec = malloc(j * sizeof(ec[0]));
        }
    }
    return ec;
}

struct find_nolref_simplex_fun_arg {
    double y;
    v3 ref_xyz;
    v3 d50;
};

static double
find_nolref_simplex_fun(double m[],
                               void *arg)
{
    for (int i = 0; i < 2; i++) {
        if (m[i] < 0 || m[i] > 1) {
            return 100000000;
        }
    }
    struct find_nolref_simplex_fun_arg *a = (struct find_nolref_simplex_fun_arg *)arg;
    v3 apx = {{ m[0], a->y, m[1] }};
    double deC = ciede2000_lch(a->ref_xyz, apx, a->d50, 100000, 1, 100000);
    double deH = ciede2000_lch(a->ref_xyz, apx, a->d50, 100000, 100000, 1);
    return deC+deH;
}

static v3
find_nolightness_reference(v3 cam_xyz,
                           v3 ref_xyz,
                           v3 d50)
{
    struct find_nolref_simplex_fun_arg arg;
    arg.y = cam_xyz.v[1];
    arg.ref_xyz = ref_xyz;
    arg.d50 = d50;
    double m[2] = { ref_xyz.v[0], ref_xyz.v[2] };
    simplex(find_nolref_simplex_fun, &arg, m, 2, 1.0e-7, 0.01, NULL);
    return (v3){{ m[0], arg.y, m[1] }};
}

static double
worst_stretch(const m3x3 cam2xyz,
              const struct patch_set *ps,
              v3 **cp,
              int cp_count,
              int j) // index of patch to test (reference patch)
{
    double max_stretch = 0;
    v3 ltu = xyz2lutspace(m3x3_mul_v3(cam2xyz, ps->patch[j].cam));
    double c1m[3] = { 0, ltu.v[1], ltu.v[2] }; // ref patch matrix position
    ltu.v[0]  = cp[0][j].v[1];
    ltu.v[1] += cp[1][j].v[1];
    ltu.v[2] += cp[2][j].v[1];
    double c1L[3] = { ltu.v[0], ltu.v[1], ltu.v[2] }; // ref patch LUT position
    //double sum = 0;
    for (int i = 0; i < cp_count; i++) {
        if (i == j) continue;

        ltu = xyz2lutspace(m3x3_mul_v3(cam2xyz, ps->patch[i].cam));
        double c2m[3] = { 0, ltu.v[1], ltu.v[2] }; // matrix only
        ltu.v[0]  = cp[0][i].v[1];
        ltu.v[1] += cp[1][i].v[1];
        ltu.v[2] += cp[2][i].v[1];
        double c2L[3] = { ltu.v[0], ltu.v[1], ltu.v[2] }; // with LUT

        double cdm[3] = { c2m[0]-c1m[0], c2m[1]-c1m[1], c2m[2]-c1m[2] };
        double cdL[3] = { c2L[0]-c1L[0], c2L[1]-c1L[1], c2L[2]-c1L[2] };
        double d1m = sqrt(cdm[0]*cdm[0] + cdm[1]*cdm[1] + cdm[2]*cdm[2]);
        double d1L = sqrt(cdL[0]*cdL[0] + cdL[1]*cdL[1] + cdL[2]*cdL[2]);

        double adiff1 = atan2(cdm[1], cdm[2]) - atan2(cdL[1], cdL[2]);
        double adiff2 = atan2(d1m, cdm[0]) - atan2(d1L, cdL[0]);
        double adiff = fabs(adiff1) + fabs(adiff2);

        double stretch;
        if (d1L < d1m) {
            stretch = d1m / d1L;
        } else {
            stretch = d1L / d1m;
        }
        stretch *= adiff;
        //sum += stretch;
        if (stretch > max_stretch) {
            max_stretch = stretch;
        }
    }
    //elog("%f\n", sum);
    //return sum;
    return max_stretch;
}

chromalut_t *
chromalut_new(const m3x3 fm,
              const v3 wb,
              const v3 d50,
              const struct patch_set *normalized_ps,
              const double min_chroma_distance,
              const double compress_factor)
{
    const m3x3 cam2xyz = m3x3_mul_m3x3(fm, m3x3_invert(v3_asdiagonal(wb)));
    const m3x3 ifm = m3x3_invert(fm);
    struct patch_set *ps = target_copy(normalized_ps);

    elog("Making 2.5D chromaticity-addressed lookup table for XYZ correction...\n");

    if (min_chroma_distance > 0.0) {
        const int orig_count = ps->patch_count;
        const v3 ltu_wp = xyz2lutspace(d50);
        struct patch_set *ps1 = malloc(sizeof(*ps) + ps->patch_count * sizeof(ps->patch[0]));
        *ps1 = *ps;
        ps1->patch_count = 0;
        int *cgrp = calloc(1, sizeof(*cgrp) * ps->patch_count);
        bool *used = calloc(1, sizeof(*used) * ps->patch_count);
        int grouped = 0;
        int white_skip = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            v3 ltu1 = xyz2lutspace(m3x3_mul_v3(cam2xyz, ps->patch[i].cam));
            ltu1.v[0] = ltu_wp.v[0];
            double uv_de = euclidean_dist(ltu1, ltu_wp);
            if (uv_de < min_chroma_distance) {
                white_skip++;
                continue;
            }
            ps1->patch[ps1->patch_count++] = ps->patch[i];
        }
        if (white_skip > 0) {
            struct patch_set *tmp = ps;
            ps = ps1;
            ps1 = tmp;
        }
        ps1->patch_count = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            if (used[i]) continue;
            ps1->patch[ps1->patch_count++] = ps->patch[i];
            used[i] = true;
            const int base_i = ps1->patch_count - 1;
            cgrp[base_i] = 1;
            const v3 ltu1 = xyz2lutspace(m3x3_mul_v3(cam2xyz, ps->patch[i].cam));
            for (int j = i+1; j < ps->patch_count; j++) {
                if (used[j]) continue;
                v3 ltu2 = xyz2lutspace(m3x3_mul_v3(cam2xyz, ps->patch[j].cam));
                ltu2.v[0] = ltu1.v[0];
                double uv_de = euclidean_dist(ltu1, ltu2);
                if (uv_de < min_chroma_distance) {
                    cgrp[base_i]++;
                    ps1->patch[ps1->patch_count++] = ps->patch[j];
                    used[j] = true;
                    grouped++;
                }
            }
        }
        free(used);
        if (grouped > 0) {
            ps->patch_count = 0;
            int largest_grp = 0;
            const int de_w_count = sizeof(ps->patch[0].lut_de_w)/sizeof(ps->patch[0].lut_de_w[0]);
            for (int i = 0; i < ps1->patch_count;  i += cgrp[i]) {
                if (cgrp[i] > 1) {
                    if (cgrp[i] > largest_grp) largest_grp = cgrp[i];
                    int cc[TARGET_MAX_CLASS_COUNT];
                    memset(cc, 0, sizeof(cc));
                    struct cam2xyz_patch patch;
                    memset(&patch, 0, sizeof(patch));
                    for (int j = i; j < i + cgrp[i]; j++) {
                        for (int k = 0; k < 3; k++) {
                            patch.cam.v[k] += ps1->patch[j].cam.v[k];
                            patch.xyz.v[k] += ps1->patch[j].xyz.v[k];
                        }
                        for (int k = 0; k < de_w_count; k++) {
                            patch.lut_de_w[k] += ps1->patch[j].lut_de_w[k];
                        }
                        cc[ps1->patch[j].class]++;
                    }
                    for (int k = 0; k < 3; k++) {
                        patch.cam.v[k] /= cgrp[i];
                        patch.xyz.v[k] /= cgrp[i];
                    }
                    for (int k = 0; k < de_w_count; k++) {
                        patch.lut_de_w[k] /= cgrp[i];
                    }
                    int max_class = 0;
                    for (int k = 0; k < TARGET_MAX_CLASS_COUNT; k++) {
                        if (cc[k] > cc[max_class]) {
                            max_class = k;
                        }
                    }
                    for (int j = i; j < i + cgrp[i]; j++) {
                        if (ps1->patch[j].class == cc[max_class]) {
                            ps->patch[ps->patch_count] = ps1->patch[j];
                            break;
                        }
                    }
                    ps->patch[ps->patch_count].cam = patch.cam;
                    ps->patch[ps->patch_count].xyz = patch.xyz;
                    for (int k = 0; k < de_w_count; k++) {
                        ps->patch[ps->patch_count].lut_de_w[k] = patch.lut_de_w[k];
                    }
                    ps->patch_count++;
                } else {
                    ps->patch[ps->patch_count] = ps1->patch[i];
                    ps->patch_count++;
                }
            }
            elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                "%.2f%% of the patches was put in a chromaticity group due to nearby neighbor.\n"
                "  %.2f%% of the patches was removed due to being nearby the whitepoint.\n"
                "  Largest chromaticity group contains %d patches. Patch count reduced from\n"
                "  %d to %d. Note that patch matching cannot reach 100%% when chromaticity\n"
                "  groups are formed, as the LUT matches the average within a group.\n",
                100.0 * (double)grouped / orig_count, 100.0 * (double)white_skip / orig_count, largest_grp, orig_count, ps->patch_count);
        }
        free(ps1);
        free(cgrp);
    }

    v3 *anchor_cp = NULL;
    int anchor_count = 0;
    {
        const int cs = 63;
        p2d *lp = malloc(cs*cs*cs*sizeof(lp[0]));
        p2d *v = malloc(cs*cs*cs*sizeof(v[0]));
        int vcount = 0;
        int lp_count = 0;
#pragma omp parallel for
        for (int r_ = 0; r_ < cs; r_++) {
            for (int g_ = 0; g_ < cs; g_++) {
                for (int b_ = 0; b_ < cs; b_++) {
                    if (r_ == 0 && g_ == 0 && b_ == 0) continue;
                    double r = (double)r_ / (cs-1);
                    double g = (double)g_ / (cs-1);
                    double b = (double)b_ / (cs-1);
                    v3 cam = {{ r, g, b}};
                    v3 xyz = m3x3_mul_v3(cam2xyz, cam);
                    if (gamut_isvalid(xyz)) {
                        v3 ltu = xyz2lutspace(xyz);
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
        lp = convex_hull(v, vcount, &lp_count);
        free(v);
        int edge_coord_count;
        p2d *edge = edge_coords(lp, lp_count, 0.04, &edge_coord_count);
        free(lp);

        for (int kk = 0; kk < 3; kk++) {
            double anchor_chroma_distance = 0.08 + kk * 0.01;
            p2d *p = malloc(edge_coord_count * sizeof(*p));
            int p_count = 0;
            for (int i = 0; i < edge_coord_count; i++) {
                double x = edge[i].x;
                double y = edge[i].y;
                double nx = i == edge_coord_count - 1 ? edge[0].x : edge[i+1].x;
                double ny = i == edge_coord_count - 1 ? edge[0].y : edge[i+1].y;
                double nx1 = i == 0 ? edge[edge_coord_count-1].x : edge[i-1].x;
                double ny1 = i == 0 ? edge[edge_coord_count-1].y : edge[i-1].y;
                double rad = atan2(ny1 - ny, nx1 - nx) + M_PI / 2.0;
                p[p_count].x = x + anchor_chroma_distance * cos(rad);
                p[p_count].y = y + anchor_chroma_distance * sin(rad);
                p_count++;
            }
            anchor_cp = realloc(anchor_cp, (anchor_count + edge_coord_count) * sizeof(anchor_cp[0]));
            for (int i = 0; i < edge_coord_count; i++) {
                anchor_cp[anchor_count++] = (v3){{ p[i].x, 0, p[i].y }};
                //printf("%f %f %f\n", p[i].x, p[i].y, 0.0);
            }
            free(p);
        }
    }

    v3 *cp[3];
    for (int i = 0; i < 3; i++) {
        cp[i] = malloc((ps->patch_count + anchor_count + 1) * sizeof(v3));
    }

    /*
    double csum = 0, rsum = 0;
    for (int i = 0; i < ps->patch_count; i++) {
        v3 cam = m3x3_mul_v3(cam2xyz, ps->patch[i].cam);
        v3 ref = ps->patch[i].xyz;
        cam = xyz2lutspace(cam);
        ref = xyz2lutspace(ref);
        csum += cam.v[0];
        rsum += ref.v[0];
    }
    double f = rsum / csum;
    */
    for (int i = 0; i < ps->patch_count; i++) {
        v3 cam = extreme_compress(fm, ifm, compress_factor, m3x3_mul_v3(cam2xyz, ps->patch[i].cam));
        v3 ref = ps->patch[i].xyz;
        cam = xyz2lutspace(cam);
        ref = xyz2lutspace(ref);
        cp[0][i] = cp[1][i] = cp[2][i] = (v3){{ cam.v[1], 0, cam.v[2] }};
        //cp[0][i].v[1] = ref.v[0] - f * cam.v[0];
        cp[0][i].v[1] = (ref.v[0] / cam.v[0] - 1.0) * 0.1;
        cp[1][i].v[1] = ref.v[1] - cam.v[1];
        cp[2][i].v[1] = ref.v[2] - cam.v[2];
        //elog("cam %f %f %f\n", cam.v[0],cam.v[1],cam.v[2]);
        //elog("ref %f %f %f\n", ref.v[0],ref.v[1],ref.v[2]);
        //printf("%f %f %f\n", cp[0][i].v[0]+cp[1][i].v[1], cp[0][i].v[2]+cp[2][i].v[1], cp[0][i].v[1]);
    }
    for (int i = ps->patch_count; i < ps->patch_count + anchor_count; i++) {
        int j = i - ps->patch_count;
        cp[0][i] = cp[1][i] = cp[2][i] = anchor_cp[j];
    }
    free(anchor_cp);

    { // add extra patch for whitepoint preservation
        const int i = ps->patch_count + anchor_count;
        v3 cam = extreme_compress(fm, ifm, compress_factor, m3x3_mul_v3(cam2xyz, wb));
        v3 ref = d50;
        cam = xyz2lutspace(cam);
        ref = xyz2lutspace(ref);
        cp[0][i] = cp[1][i] = cp[2][i] = (v3){{ cam.v[1], 0, cam.v[2] }};
        cp[0][i].v[1] = (ref.v[0] / cam.v[0] - 1.0) * 0.1;
        cp[1][i].v[1] = ref.v[1] - cam.v[1];
        cp[2][i].v[1] = ref.v[2] - cam.v[2];
    }

    const int cp_count = ps->patch_count + anchor_count + 1;

    double max_de_lim = 0.0;
    double max_de_lim_print = 0.0;
    bool axis_disabled[3] = { true, true, true };
    for (int j = 0; j < ps->patch_count; j++) {
        for (int k = 0; k < 6; k++) {
            const int w2a[] = { 0, 0, 1, 1, 2, 2 };
            double de = fabs(ps->patch[j].lut_de_w[k]);
            if (de < 50) {
                axis_disabled[w2a[k]] = false;
            }
            if (de > max_de_lim) {
                max_de_lim = de;
            }
            if (de < 50 && de > max_de_lim_print) {
                max_de_lim_print = de;
            }
        }
    }

    double relax[3] = { 0.0, 0.0, 0.0 };
    { // workaround: test if we need to add regularization
        if (axis_disabled[0]) {
            relax[0] = -1;
        }
        if (axis_disabled[1] && axis_disabled[2]) {
            // chroma/hue can't be disabled separately due to chromaticity indexing, so either both or none
            relax[1] = -1;
            relax[2] = -1;
        }
        tps_t *tps_exact[3];
        for (int i = 0; i < 3; i++) {
            if (relax[i] < 0.0) {
                tps_exact[i] = NULL;
                continue;
            }
            tps_exact[i] = tps_new(cp[i], cp_count, relax[i]);
            if (tps_exact[i] == NULL) {
                // tps seems to sometimes have some issue with 0 regularization
                relax[i] = 0.000000000001;
                tps_exact[i] = tps_new(cp[i], cp_count, relax[i]);
            }
            if (tps_exact[i] == NULL) {
                relax[i] = 0.0;
                do {
                    relax[i] += 0.0005;
                    tps_exact[i] = tps_new(cp[i], cp_count, relax[i]);
                } while (tps_exact[i] == NULL);
                elog("Exact TPS solution for axis %d was not possible, relaxed with %f.\n", i, relax[i]);
            }
        }
        for (int i = 0; i < 3; i++) {
            tps_delete(tps_exact[i]);
        }
    }

    if (axis_disabled[0] && !axis_disabled[1]) {
        // Common to disable lightness, make adjustments for better corrections
        elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
            "Lightness axis is disabled. Since lightness affects chroma, the LUT chroma\n"
            "  control points are recalculated to better match the uncorrected lightness.\n"
            "  A residual error of up to about 0.2 DE is expected.\n");
#pragma omp parallel for
        for (int i = 0; i < ps->patch_count; i++) {
            v3 cam = extreme_compress(fm, ifm, compress_factor, m3x3_mul_v3(cam2xyz, ps->patch[i].cam));
            v3 ref = ps->patch[i].xyz;

            ref = find_nolightness_reference(cam, ref, d50);

            cam = xyz2lutspace(cam);
            ref = xyz2lutspace(ref);
            cp[0][i] = cp[1][i] = cp[2][i] = (v3){{ cam.v[1], 0, cam.v[2] }};
            cp[0][i].v[1] = 0; // disabled
            cp[1][i].v[1] = ref.v[1] - cam.v[1];
            cp[2][i].v[1] = ref.v[2] - cam.v[2];
        }
    }

    if (max_de_lim > 0.0) { // relax LUT stretching according to provide max acceptable DE values

        v3 *cpa[3];
        for (int i = 0; i < 3; i++) {
            cpa[i] = malloc(cp_count * sizeof(v3));
            memcpy(cpa[i], cp[i], cp_count * sizeof(v3));
        }

        elog("Relaxing LUT stretch with up to %.2f DE. Iterating over %d patches%s...\n", max_de_lim_print, ps->patch_count,
                (ps->patch_count > 100) ? ", this can take very long time" : "");
        if (axis_disabled[0]) {
            elog("  Lightness correction is disabled.\n");
        }
        if (axis_disabled[1]) {
            elog("  Chroma correction is disabled.\n");
        }
        if (axis_disabled[2]) {
            elog("  Hue correction is disabled.\n");
        }
        struct timeval tv, tvnext;
        gettimeofday(&tv, NULL);
        tvnext = tv;
        tvnext.tv_sec += 3;
        double de_sum_max = 0;
        double de_sum_min = 1e9;
        int de_sum_not_increase_count = 0;
        int de_sum_decrease_count = 0;
        for (;;) {
            tps_t *tps_approx[3];
            for (int i = 0; i < 3; i++) {
                if (relax[i] < 0.0 || (i == 0 && axis_disabled[0])) { // fast past for most common axis disabling (lightness)
                    tps_approx[i] = NULL;
                } else {
                    tps_approx[i] = tps_new(cpa[i], cp_count, relax[i]);
                }
            }

            double de_sum = 0;
            for (int j = 0; j < ps->patch_count; j++) {
                v3 cam = m3x3_mul_v3(cam2xyz, ps->patch[j].cam);
                v3 ref = ps->patch[j].xyz;
                v3 approx = chromalut_lookup_(tps_approx, fm, ifm, compress_factor, cam, true);
                double de = ciede2000(approx, ref, d50);
                de_sum += de;
            }
            if (de_sum > de_sum_max) {
                de_sum_max = de_sum;
                de_sum_not_increase_count = 0;
            } else if (de_sum < de_sum_min) {
                de_sum_min = de_sum;
                de_sum_not_increase_count = 0;
                de_sum_decrease_count++;
            } else {
                de_sum_not_increase_count++;
                if (de_sum_decrease_count > 2) {
                    de_sum_max = de_sum;
                }
                de_sum_decrease_count = 0;
            }
            //elog("de_sum: %f %d\n", de_sum, de_sum_not_increase_count);
            if (de_sum_not_increase_count > 10) {
                for (int i = 0; i < 3; i++) {
                    memcpy(cp[i], cpa[i], cp_count * sizeof(v3));
                    free(cpa[i]);
                    tps_delete(tps_approx[i]);
                }
                int zero_count = 0;
                bool non_zero[3] = { false, false, false };
                for (int j = 0; j < ps->patch_count; j++) {
                    if (cp[0][j].v[1] == 0 && cp[1][j].v[1] == 0 && cp[2][j].v[1] == 0) {
                        zero_count++;
                    }
                    for (int i = 0; i < 3; i++) {
                        if (cp[i][j].v[1] != 0) {
                            non_zero[i] = true;
                        }
                    }
                }
                elog(//2345678901234567890123456789012345678901234567890123456789012345678901234567890
                    "Average DE for the %d tested patches increased to %.2f after LUT relax.\n"
                    "  %.2f%% could do without LUT correction.\n",
                    ps->patch_count, de_sum / ps->patch_count, (100.0 * zero_count) / ps->patch_count);
                if (zero_count == ps->patch_count) {
                    elog("  LUT was disabled as there is no correction need.\n");
                    axis_disabled[0] = true;
                    axis_disabled[1] = true;
                    axis_disabled[2] = true;
                    relax[0] = -1;
                    relax[1] = -1;
                    relax[2] = -1;
                } else {
                    for (int i = 0; i < 3; i++) {
                        if (!axis_disabled[i] && !non_zero[i]) {
                            axis_disabled[i] = true;
                            elog("  %s LUT axis was disabled as there is no correction need.\n",
                                    i == 0 ? "Lightness" : i == 1 ? "Chromaticity u'" : "Chromaticity v'");
                            relax[i] = -1;
                        }
                    }
                }
                break;
            }
            gettimeofday(&tv, NULL);
            if (tv.tv_sec >= tvnext.tv_sec) {
                elog(".");
                tvnext = tv;
                tvnext.tv_sec += 1;
            }
            for (int j = 0; j < ps->patch_count; j++) {
                v3 cam = m3x3_mul_v3(cam2xyz, ps->patch[j].cam);
                v3 ref = ps->patch[j].xyz;
                v3 apx = chromalut_lookup_(tps_approx, fm, ifm, compress_factor, cam, true);

                v3 reflch = xyz2lch(ref, d50);
                v3 apxlch = xyz2lch(apx, d50);
                double deL = ciede2000_lch(ref, apx, d50, 1, 100000, 100000);
                if (reflch.v[0] > apxlch.v[0]) {
                    deL = -deL;
                }
                double deC = ciede2000_lch(ref, apx, d50, 100000, 1, 100000);
                if (reflch.v[1] > apxlch.v[1]) {
                    deC = -deC;
                }
                double deH = ciede2000_lch(ref, apx, d50, 100000, 100000, 1);
                if (reflch.v[2] > apxlch.v[2]) {
                    deH = -deH;
                }

                apx = xyz2lutspace(apx);
                cam = xyz2lutspace(cam);
                ref = xyz2lutspace(ref);
                for (int i = 0; i < 3; i++) {
                    double pre_relax[3] = { cpa[0][j].v[1], cpa[1][j].v[1], cpa[2][j].v[1] };
                    double stretch_before = worst_stretch(cam2xyz, ps, cpa, ps->patch_count, j);
                    bool rel = false;
                    const double *de_w = ps->patch[j].lut_de_w;
                    if ((i == 0 && (deL < de_w[0]-0.1 || deL > de_w[1]+0.1)) ||
                        (i != 0 && (deC < de_w[2]-0.1 || deC > de_w[3]+0.1 || (deH < de_w[4]-0.1 || deH > de_w[5]+0.1))))
                    {
                        // increase correction
                        if (fabs(cpa[i][j].v[1] - cp[i][j].v[1]) <= 0.00005) {
                            cpa[i][j].v[1] = cp[i][j].v[1]; // full correction
                        } else if (cpa[i][j].v[1] < cp[i][j].v[1]) {
                            cpa[i][j].v[1] += 0.00005;
                        } else {
                            cpa[i][j].v[1] -= 0.00005;
                        }
                    } else {
                        // decrease correction
                        if (fabs(cpa[i][j].v[1]) <= 0.00005) {
                            cpa[i][j].v[1] = 0.0; // no correction
                        } else if (cpa[i][j].v[1] < 0) {
                            cpa[i][j].v[1] += 0.00005;
                        } else {
                            cpa[i][j].v[1] -= 0.00005;
                        }
                        rel = true;
                    }
                    double stretch_after = worst_stretch(cam2xyz, ps, cpa, ps->patch_count, j);
                    if (stretch_after > stretch_before * 1.02) {
                        //elog("w");
                        cpa[0][j].v[1] = pre_relax[0];
                        cpa[1][j].v[1] = pre_relax[1];
                        cpa[2][j].v[1] = pre_relax[2];
                        if (rel) {
                            //elog("worse stretch during relax %d (%s) %f => %f\n", i, ps->patch[j].str_id, stretch_before, stretch_after);
                        } else {
                            //elog("worse stretch during relax %d (%s) %f => %f\n", i, ps->patch[j].str_id, stretch_before, stretch_after);
                        }
                    }
                }
            }
            if (0) { // test approx
                for (int j = 0; j < ps->patch_count; j++) {
                    v3 cxyz = m3x3_mul_v3(cam2xyz, ps->patch[j].cam);
                    v3 xyz = ps->patch[j].xyz;
                    cxyz = chromalut_lookup_(tps_approx, fm, ifm, compress_factor, cxyz, true);

                    double errors[4];
                    errors[0] = ciede2000(xyz, cxyz, d50);
                    errors[1] = ciede2000_lch(xyz, cxyz, d50, 1, 100000, 100000);
                    if (xyz2lab(xyz, d50).v[0] > xyz2lab(cxyz, d50).v[0]) {
                        errors[1] = -errors[1];
                    }
                    errors[2] = ciede2000_lch(xyz, cxyz, d50, 100000, 1, 100000);
                    {
                        v3 lab = xyz2lab(xyz, d50);
                        double c1 = sqrt(lab.v[1] * lab.v[1] + lab.v[2] * lab.v[2]);
                        lab = xyz2lab(cxyz, d50);
                        double c2 = sqrt(lab.v[1] * lab.v[1] + lab.v[2] * lab.v[2]);
                        if (c1 > c2) {
                            errors[2] = -errors[2];
                        }
                    }
                    errors[3] = ciede2000_lch(xyz, cxyz, d50, 100000, 100000, 1);
                    elog("%d %s DE %.2f DE LCh %+.2f %+.2f %.2f\n", j, ps->patch[j].str_id, errors[0], errors[1], errors[2], errors[3]);
                }
                getchar();
            }
            for (int i = 0; i < 3; i++) {
                tps_delete(tps_approx[i]);
            }
        }

        /*
        for (;;) {
            for (int i = 0; i < 3; i++) {
                tps_approx[i] = tps_new(cpa[i], ps->patch_count, relax[i]);
            }
            bool complete = true;

            { // test approx
                for (int j = 0; j < ps->patch_count; j++) {
                    v3 cam = m3x3_mul_v3(cam2xyz, ps->patch[j].cam);
                    v3 ref = ps->patch[j].xyz;
                    v3 approx = chromalut_lookup_(tps_approx, cam);
                    double de = ciede2000(approx, ref, d50);
                    //elog("%d DE %f %d\n", j, de, rdir[j]);
                }
                //getchar();
            }

            {
                complete = true;
                for (int j = 0; j < ps->patch_count; j++) {
                    v3 cam = m3x3_mul_v3(cam2xyz, ps->patch[j].cam);
                    v3 ref = ps->patch[j].xyz;
                    v3 apx = chromalut_lookup_(tps_approx, cam);
                    double de = ciede2000(apx, ref, d50);
                    apx = xyz2lutspace(apx);
                    cam = xyz2lutspace(cam);
                    ref = xyz2lutspace(ref);
                    if (de > ps->patch[j].de && rdir[j] != +2) {
                        for (int i = 0; i < 3; i++) {
                            if (cpa[i][j].v[1] < 0) {
                                cpa[i][j].v[1] -= 0.00005;
                            } else {
                                cpa[i][j].v[1] += 0.00005;
                            }
                        }
                        rdir[j] = -1;
                        complete = false;
                    } else if (de < ps->patch[j].de - 0.1 && rdir[j] != +2) {
                        for (int i = 0; i < 3; i++) {
                            if (fabs(cpa[i][j].v[1]) <= 0.00005) {
                                cpa[i][j].v[1] = 0.0;
                            } else if (cpa[i][j].v[1] < 0) {
                                cpa[i][j].v[1] += 0.00005;
                            } else {
                                cpa[i][j].v[1] -= 0.00005;
                            }
                        }
                        bool finished = true;
                        for (int i = 0; i < 3; i++) {
                            if (cpa[i][j].v[1] != 0.0) {
                                finished = false;
                            }
                        }
                        if (!finished) {
                            rdir[j] = +1;
                            complete = false;
                        } else {
                            rdir[j] = +2;
                        }
                    } else {
                        rdir[j] = +2;
                    }
                }
            }

            static int cc = 0;
            if (++cc == 1000) break;
            if (complete) {
                break;
            }
            for (int i = 0; i < 3; i++) {
                tps_delete(tps_approx[i]);
            }
        }
        */
    }

    chromalut_t *lut = calloc(1, sizeof(*lut));
    for (int i = 0; i < 3; i++) {
        lut->regularization[i] = relax[i];
        if (relax[i] < 0.0) {
            lut->tps[i] = NULL;
            lut->cp[i] = NULL;
            free(cp[i]);
        } else {
            lut->tps[i] = tps_new(cp[i], cp_count, relax[i]);
            if (lut->tps[i] == NULL) {
                elog("Failed to generate TPS! Aborting\n");
                abort();
            }
            lut->cp[i] = cp[i];
        }
    }
    lut->cp_count = cp_count;
    lut->fm = fm;
    lut->ifm = ifm;
    lut->compress_factor = compress_factor;

    return lut;
}

void
chromalut_delete(chromalut_t *lut)
{
    if (lut == NULL) {
        return;
    }
    for (int i = 0; i < 3; i++) {
        tps_delete(lut->tps[i]);
        free(lut->cp[i]);
    }
    free(lut);
}

v3
chromalut_lookup_diff(chromalut_t *lut,
                      double u,
                      double v)
{
    double Ldiff = tps_interpolate_(lut->tps[0], u, v);
    double udiff = tps_interpolate_(lut->tps[1], u, v);
    double vdiff = tps_interpolate_(lut->tps[2], u, v);
    return (v3){{ Ldiff, udiff, vdiff }};
}

v3
chromalut_lookup(chromalut_t *lut,
                 v3 xyz)
{
    if (lut == NULL) {
        return xyz;
    }
    xyz = chromalut_lookup_(lut->tps, lut->fm, lut->ifm, lut->compress_factor, xyz, true);
    return xyz;
}

v3
chromalut_lookup_noclip(chromalut_t *lut,
                        v3 xyz)
{
    if (lut == NULL) {
        return xyz;
    }
    xyz = chromalut_lookup_(lut->tps, lut->fm, lut->ifm, lut->compress_factor, xyz, false);
    return xyz;
}

/*
struct chromalut_reverse_lookup_simplex_fun_arg {
    chromalut_t *lut;
    v3 xyz;
};

static double
chromalut_reverse_lookup_simplex_fun(double m[],
                                     void *arg)
{
    struct chromalut_reverse_lookup_simplex_fun_arg *a = (struct chromalut_reverse_lookup_simplex_fun_arg *)arg;
    v3 xyz = {{ m[0], m[1], m[2] }};
    xyz = chromalut_lookup(a->lut, xyz);
    return euclidean_dist(xyz, a->xyz);
}

v3
chromalut_reverse_lookup(chromalut_t *lut,
                         v3 xyz)
{
    double m[3] = { xyz.v[0], xyz.v[1], xyz.v[2] };
    if (m[0] < 0 || m[1] < 0 || m[2] < 0) {
        return (v3){{ -1, -1, -1 }};
    }

    struct chromalut_reverse_lookup_simplex_fun_arg arg;
    arg.lut = lut;
    arg.xyz = xyz;
    double min = simplex(chromalut_reverse_lookup_simplex_fun, &arg, m, 3, 1.0e-7, 0.01, NULL);
    if (min > 0.0001) {
        //elog("%f %f %f %f\n", min, m[0], m[1], m[2]);
        return (v3){{ -1, -1, -1 }};
    }
    return (v3){{ m[0], m[1], m[2] }};

}
*/

static void
print_cp_array(FILE *stream,
               const v3 *cp,
               int cp_count)
{
    fprintf(stream, "[ ");
    if (cp != NULL) {
        for (int i = 0; i < cp_count; i++) {
            fprintf(stream, "%g, %g, %g%s", cp[i].v[0], cp[i].v[2], cp[i].v[1], i < cp_count - 1 ? ", " : "");
        }
    }
    fprintf(stream, " ]");
}

void
chromalut_json_print(FILE *stream,
                     const char indent[],
                     const char name[],
                     chromalut_t *lut)
{
    fprintf(stream, "%s\"%s\": {\n", indent, name);
    fprintf(stream, "%s  \"compressFactor\": %g,\n", indent, lut->compress_factor);
    fprintf(stream, "%s  \"regularization\": [ %g, %g, %g ],\n", indent, lut->regularization[0], lut->regularization[1], lut->regularization[2]);
    fprintf(stream, "%s  \"uvL\": ", indent);
    print_cp_array(stream, lut->cp[0], lut->cp_count);
    fprintf(stream, ",\n%s  \"uvU\": ", indent);
    print_cp_array(stream, lut->cp[1], lut->cp_count);
    fprintf(stream, ",\n%s  \"uvV\": ", indent);
    print_cp_array(stream, lut->cp[2], lut->cp_count);
    fprintf(stream, "\n%s}", indent);
}

chromalut_t *
chromalut_new_from_data(v3 *cp[3],
                        int cp_count,
                        double regularization[3],
                        const m3x3 fm,
                        double compress_factor)
{
    chromalut_t *lut = calloc(1, sizeof(*lut));
    for (int i = 0; i < 3; i++) {
        lut->regularization[i] = regularization[i];
        if (lut->regularization[i] < 0) {
            lut->tps[i] = NULL;
            lut->cp[i] = NULL;
        } else {
            lut->tps[i] = tps_new(cp[i], cp_count, lut->regularization[i]);
            if (lut->tps[i] == NULL) {
                elog("Failed to generate TPS! Aborting\n");
                abort();
            }
            lut->cp[i] = malloc(cp_count * sizeof(v3));
            memcpy(lut->cp[i], cp[i], cp_count * sizeof(v3));
        }
    }
    lut->cp_count = cp_count;
    lut->fm = fm;
    lut->ifm = m3x3_invert(lut->fm);
    lut->compress_factor = compress_factor;
    return lut;
}

bool
chromalut_isnoop(chromalut_t *lut)
{
    if (lut == NULL) {
        return true;
    }
    for (int i = 0; i < 3; i++) {
        if (lut->tps[i] != NULL) return false;
    }
    return true;
}
