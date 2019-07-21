/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdbool.h>
#include <assert.h>

#include <elog.h>
#include <nmsimplex.h>
#include <xyz2spec.h>

enum xyz2spec_mode {
    XYZ2SPEC_MODE_GAUSSIAN_BANDS = 0,
    XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_LEFT,
    XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_RIGHT,
    XYZ2SPEC_MODE_GAUSSIAN,
    XYZ2SPEC_MODE_2X_GAUSSIAN,
    XYZ2SPEC_MODE_4X_GAUSSIAN,
#define XYZ2SPEC_MODE_FIRST XYZ2SPEC_MODE_GAUSSIAN_BANDS
#define XYZ2SPEC_MODE_FIRST_EXTREME XYZ2SPEC_MODE_4X_GAUSSIAN
#define XYZ2SPEC_MODE_LAST XYZ2SPEC_MODE_4X_GAUSSIAN
};

struct xyz2spec_simplex_fun_arg {
    double *obs_max;
    double *bands;
    int band_count;
    int m_count;
    enum xyz2spec_mode mode;
    const struct observer *obs;
    const spectrum_t *ill;
    double xyz_norm;
    v3 xyz;
    spectrum_t *last;
    bool smoothness_exactness_tradeoff;
    double freq_lo;
    double freq_hi;
    double last_sf;
};

static double
cmp_xyz(const spectrum_t *s,
        struct xyz2spec_simplex_fun_arg *arg)
{
    const v3 xyz_ref = arg->xyz;
    v3 xyz = spec2xyz_ill(s, arg->obs, arg->ill, 1);
    if (!isfinite(xyz.v[0]) || !isfinite(xyz.v[1]) || !isfinite(xyz.v[2])) {
        xyz = (v3){{0,0,0}};
    }
    /*
    if (arg->ill != NULL) {
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < 3; i++) {
            sum1 += xyz.v[i];
            sum2 += xyz_ref.v[i];
        }
        double f = sum2 / sum1;
        if (!isfinite(f)) {
            f = 1.0;
        }
        if (f > 1) {
            f = 1;
        }
        for (int i = 0; i < 3; i++) {
            xyz.v[i] *= f;
        }
    }
    */
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        double diff = xyz.v[i] - xyz_ref.v[i];
        diff *= 100;
        sum += diff * diff;
    }
    return sqrt(sum);
}

static double
xyz2spec_rough_simplex_fun(double m[],
                           void *arg)
{
    struct xyz2spec_simplex_fun_arg *a = (struct xyz2spec_simplex_fun_arg *)arg;
    double amp[a->band_count];
    double range = a->freq_hi - a->freq_lo;
    const double worst = 100000000;
    if (m[0] < 0 || m[1] < 0 || m[2] < 0) {
        return worst;
    }
    for (int i = 0; i < a->band_count; i++) {
        if (a->bands[i] < a->freq_lo) {
            amp[i] = 0;
        } else if (a->bands[i] < a->freq_hi) {
            double x = (a->bands[i] - a->freq_lo) / range;
            if (x < 0.3333) {
                amp[i] = m[0];
            } else if (x < 0.6666) {
                amp[i] = m[1];
            } else {
                amp[i] = m[2];
            }
        } else {
            amp[i] = 0;
        }
    }
    spectrum_t *s = spec_alloc1(a->bands, amp, a->band_count);
    double de = cmp_xyz(s, a);
    if (0) {
        FILE *stream = fopen("dump1/r.dat", "w");
        spec_print(s, stream);
        fclose(stream);
    }
    //elog("DE %f (%f %f %f)\n", de, m[0], m[1], m[2]);

    free(s);
    return de;
}

static double
xyz2spec_simplex_fun(double m[],
                     void *arg)
{
    struct xyz2spec_simplex_fun_arg *a = (struct xyz2spec_simplex_fun_arg *)arg;
    const double worst = 100000000;
    double amp[a->band_count];

    // render spectrum
    switch (a->mode) {
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS:
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_LEFT:
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_RIGHT:
    {
        const int gauss_count = a->m_count;
        for (int i = 0; i < gauss_count; i++) {
            if (m[i] < 0) {
                return worst;
            }
        }
        double offset;
        switch (a->mode) {
        default: offset = 0.0; break;
        case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_LEFT: offset = -0.05; break;
        case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_RIGHT: offset = +0.05; break;
        }
        double band_lo = a->freq_lo;
        double fspan = a->freq_hi - a->freq_lo;
        for (int i = 0; i < a->band_count; i++) {
            double x = (double)i / (a->band_count-1);
            double f = a->bands[i] - band_lo;
            f /= fspan;
            x = f;
            amp[i] = 0;
            for (int j = 0; j < gauss_count; j++) {
                double p0 = (j-1) * 0.1 + offset;
                double xp = (x - p0) * 10;
                amp[i] += m[j] * exp(-(xp*xp));
            }
        }
        /*
        FILE *stream = fopen("dump1/s.dat", "w");
        for (int i = 0; i < a->band_count; i++) {
            double x = (double)i / (a->band_count-1);
            //fprintf(stream, "%f %f\n", x, amp[i]);
            fprintf(stream, "%f %f\n", a->bands[i], amp[i]);
        }
        fclose(stream);

        exit(1);
        */
        break;
    }
    case XYZ2SPEC_MODE_GAUSSIAN:
    case XYZ2SPEC_MODE_2X_GAUSSIAN:
    case XYZ2SPEC_MODE_4X_GAUSSIAN: {
        const int gauss_count = a->m_count / 3;
        for (int j = 0; j < gauss_count; j++) {
            double *p = &m[j*3];
            if (p[0] < 0.0 || p[1] <= 0 || p[2] <= 0) {
                return worst;
            }
        }
        for (int i = 0; i < a->band_count; i++) {
            double x = (double)i / (a->band_count-1);
            amp[i] = 0;
            for (int j = 0; j < gauss_count; j++) {
                double *p = &m[j*3];
                double xp = (x - p[0]) / p[2];
                amp[i] += p[1]*exp(-(xp*xp));
            }
        }
        break;
    }
    }

    double sf = 0;
    { // calculate smoothness factor
        sf = 0;
        double der[a->band_count];
        for (int i = 1; i < a->band_count; i++) {
            der[i] = amp[i] - amp[i-1];
        }
        for (int i = 2; i < a->band_count; i++) {
            double der2 = der[i] - der[i-1];
            double f = fabs(der2);
            { // we don't like rise on the sides
                if (a->bands[i] < a->freq_lo + 25) {
                    if (der[i] < 0) {
                        f *= 2;
                    }
                } else if (a->bands[i] > a->freq_hi - 25) {
                    if (der[i] > 0) {
                        f *= 2;
                    }
                }
            }
            sf += f;
        }
        sf /= a->band_count - 2;
        sf *= 100;
    }

    // test spectrum
    free(a->last);
    a->last = spec_alloc1(a->bands, amp, a->band_count);

    /*
    if (0) {
        if (a->mode == XYZ2SPEC_MODE_4X_GAUSSIAN) {
            FILE *stream = fopen("dump1/s.dat", "w");
            spec_print(a->last, stream);
            fclose(stream);
            exit(1);
        }
    }
    */

    assert(a->last);
    double de = cmp_xyz(a->last, a);

    a->last_sf = sf;
    double f;
    if (a->smoothness_exactness_tradeoff) {
        f = de + sf;
    } else {
        if (de < 0.000001) {
            f = de;
        } else if (de < 1.0) {
            f = de + sf * pow(de, 0.25);
            if (f < 0.000015) {
                f = de + sf * de;
            }
        } else {
            f = de + sf * de;
        }
    }
    if (0) {
        elog("F %f DE %f (SF %f):", f, de, sf);
        for (int j = 0; j < a->m_count; j++) {
            elog(" %.4f", m[j]);
        }
        elog("\n");
    }
    return f;
}

static void
init_xyz2spec_simplex(double m[],
                      struct xyz2spec_simplex_fun_arg *arg,
                      enum xyz2spec_mode mode)
{
    double pm[3] = {0,0,0};
    simplex(xyz2spec_rough_simplex_fun, arg, pm, 3, 1.0e-5, 0.08, NULL);
    arg->mode = mode;
    switch (mode) {
    case XYZ2SPEC_MODE_GAUSSIAN: {
        int maxi = 0;
        for (int i = 1; i < 3; i++) if (pm[i] > pm[maxi]) maxi = i;
        m[0] = arg->freq_lo + (1.0/6.0 + maxi * 1.0/3.0) * (arg->freq_hi - arg->freq_lo);
        m[0] -= arg->bands[0];
        m[0] /= arg->bands[arg->band_count-1] - arg->bands[0];
        m[1] = pm[maxi];
        m[2] = 0.2;
        arg->m_count = 3;
        break;
    }
    case XYZ2SPEC_MODE_2X_GAUSSIAN: {
        int maxi = 0;
        for (int i = 1; i < 3; i++) if (pm[i] > pm[maxi]) maxi = i;
        m[0] = arg->freq_lo + (1.0/6.0 + maxi * 1.0/3.0) * (arg->freq_hi - arg->freq_lo);
        m[0] -= arg->bands[0];
        m[0] /= arg->bands[arg->band_count-1] - arg->bands[0];
        m[1] = pm[maxi];
        m[2] = 0.15;

        pm[maxi] = 0.0;
        for (int i = 0; i < 3; i++) if (pm[i] > pm[maxi]) maxi = i;
        m[3] = arg->freq_lo + (1.0/6.0 + maxi * 1.0/3.0) * (arg->freq_hi - arg->freq_lo);
        m[3] -= arg->bands[0];
        m[3] /= arg->bands[arg->band_count-1] - arg->bands[0];
        m[4] = pm[maxi];
        m[5] = 0.15;
        arg->m_count = 6;
        break;
    }
    case XYZ2SPEC_MODE_4X_GAUSSIAN: {
        int maxi = 0;
        for (int i = 1; i < 3; i++) if (pm[i] > pm[maxi]) maxi = i;
        m[0] = arg->freq_lo + (1.0/6.0 + maxi * 1.0/3.0) * (arg->freq_hi - arg->freq_lo);
        m[0] -= arg->bands[0];
        m[0] /= arg->bands[arg->band_count-1] - arg->bands[0];
        m[1] = pm[maxi] / 2;
        m[2] = 0.1;

        pm[maxi] = 0.0;
        for (int i = 0; i < 3; i++) if (pm[i] > pm[maxi]) maxi = i;
        m[3] = arg->freq_lo + (1.0/6.0 + maxi * 1.0/3.0) * (arg->freq_hi - arg->freq_lo);
        m[3] -= arg->bands[0];
        m[3] /= arg->bands[arg->band_count-1] - arg->bands[0];
        m[4] = pm[maxi] / 2;
        m[5] = 0.1;

        m[6] = m[0];
        m[7] = m[1];
        m[8] = m[2];
        m[9] = m[3];
        m[10] = m[4];
        m[11] = m[5];
        arg->m_count = 12;
        break;
    }
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS:
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_LEFT:
    case XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_RIGHT: {
        double f = 0.57;
        m[0] = pm[0] * f / 4;
        m[1] = pm[0] * f;
        m[2] = pm[0] * f;
        m[3] = pm[0] * f;
        m[4] = (pm[0]*0.66+pm[1]*0.33) * f;
        m[5] = pm[1] * f;
        m[6] = pm[1] * f;
        m[7] = pm[1] * f;
        m[8] = (pm[1]*0.33+pm[2]*0.66) * f;
        m[9] = pm[2] * f;
        m[10] = pm[2] * f;
        m[11] = pm[2] * f;
        m[12] = pm[2] * f / 4;
        arg->m_count = 13;
        break;
    }
    }
}

static double
run_xyz2spec_simplex(double m[],
                     struct xyz2spec_simplex_fun_arg *arg)
{
    double min, oldmin = -1;
    for (int i = 0; i < 10000; i++) {
        min = simplex(xyz2spec_simplex_fun, arg, m, arg->m_count, 1.0e-7, 0.01, NULL);
        if (oldmin + 0.000001 > min && oldmin - 0.000001 < min) break;
        oldmin = min;
    }
    return min;
}

static void
refine_xyz2spec_simplex(double m[],
                        struct xyz2spec_simplex_fun_arg *arg)
{
    const double max_good_de = 0.00001;
    double last_success[arg->m_count];
    double best_sf = arg->last_sf;
    spectrum_t *best_spectrum = spec_copy(arg->last);
    memcpy(last_success, m, arg->m_count * sizeof(m[0]));

    arg->smoothness_exactness_tradeoff = true;
    double min = run_xyz2spec_simplex(m, arg);
    if (0) {
        FILE *stream = fopen("dump1/s.dat", "w");
        spectrum_t *s = spec_copy(arg->last);
        spec_scalar_multiply(s, 1.0 / arg->xyz_norm);
        spec_print(s, stream);
        free(s);
        fclose(stream);
        elog("smoothness tradeoff %f\n", min);
        getchar();
    }
    arg->smoothness_exactness_tradeoff = false;
    min = run_xyz2spec_simplex(m, arg);
    if (0) {
        FILE *stream = fopen("dump1/s.dat", "w");
        spec_print(arg->last, stream);
        fclose(stream);
        elog("refine exact %f (%f %f)\n", min, arg->last_sf, best_sf);
        getchar();
    }
    if (min <= max_good_de && arg->last_sf < best_sf) {
        free(best_spectrum);
        best_spectrum = arg->last;
    } else {
        free(arg->last);
    }
    arg->smoothness_exactness_tradeoff = false;
    arg->last = best_spectrum;
}

static double
midpoint(double f1,
         double mid_f,
         double f2,
         double e1,
         double e2)
{
    if (mid_f == f1) return e1;
    if (mid_f == f2) return e2;
    double b = e2 - e1;
    double a = f2 - f1;
    double a2 = mid_f - f1;
    return e1 + a2*(b/a);
}

static void
optimize_observer(struct observer *out_obs,
                  double *freq_lo_,
                  double *freq_hi_,
                  const struct observer *obs)
{
    // find side frequencies where observer response has risen with a meaningful amaount
    double freq_lo = 500;
    double freq_hi = 600;
    double freq_min = freq_lo;
    double freq_max = freq_hi;
    int bsp = 5;
    const double limit = 0.05;
    for (int i = 0; i < 3; i++) {
        int bc;
        double *bands, *amp;
        spec_getraw(obs->cmf[i], &bc, &bands, &amp);
        bsp = (bands[1] - bands[0]);
        if (bands[0] < freq_min) freq_min = bands[0];
        if (bands[bc-1] > freq_max) freq_max = bands[bc-1];
        double lo = -1, hi = -1;
        for (int j = 0; j < bc; j++) {
            if (amp[j] >= 0.1) {
                if (j > 0) {
                    lo = midpoint(amp[j-1], limit, amp[j], bands[j-1], bands[j]);
                } else {
                    lo = bands[j];
                }
                break;
            }
        }
        for (int j = bc-1; j >= 0; j--) {
            if (amp[j] >= 0.1) {
                if (j < bc-1) {
                    hi = midpoint(amp[j], limit, amp[j+1], bands[j], bands[j+1]);
                } else {
                    hi = bands[j];
                }
                break;
            }
        }
        if (lo != -1 && lo < freq_lo) {
            freq_lo = lo;
        }
        if (hi != -1 && hi > freq_hi) {
            freq_hi = hi;
        }
        free(bands);
        free(amp);
    }
    if (0) {
        FILE *stream = fopen("dump1/s1.dat", "w");
        fprintf(stream, "%f %f\n", freq_lo, 1.0);
        fprintf(stream, "%f %f\n", freq_hi, 1.0);
        fclose(stream);
    }
    *freq_lo_ = freq_lo;
    *freq_hi_ = freq_hi;
    if (bsp < 5) bsp = 5;

    { // generate new observer
        int lo = (int)(freq_lo - 100);
        while (lo % bsp != 0) lo--;
        while (lo > freq_min) lo -= bsp;
        if (lo < 0) lo = 0;
        int hi = (int)(freq_hi + 100);
        while (hi % bsp != 0) hi++;
        while (hi < freq_max) hi += bsp;
        int bc = (hi - lo) / bsp + 1;
        for (int j = 0; j < 3; j++) {
            double *amp = malloc(bc * sizeof(double));;
            for (int i = 0; i < bc; i++) {
                amp[i] = spec_at(obs->cmf[j], (double)(lo + i*bsp));
            }
            out_obs->cmf[j] = spec_alloc(amp, lo, hi, bsp);
        }
    }
}

spectrum_t *
xyz2spec(const v3 xyz,
         const struct observer *obs,
         const spectrum_t *illuminant)
{
    if ((xyz.v[0] == 0 && xyz.v[1] == 0 && xyz.v[2] == 0) ||
        (xyz.v[0] < 0 || xyz.v[1] < 0 || xyz.v[2] < 0) ||
        (!isfinite(xyz.v[0]) || !isfinite(xyz.v[1]) || !isfinite(xyz.v[2])))
    {
        return NULL;
    }
    const double max_good_de = 0.00001;
    double min;
    double m[64];
    struct observer opt_obs;
    struct xyz2spec_simplex_fun_arg arg;
    memset(&arg, 0, sizeof(arg));

    optimize_observer(&opt_obs, &arg.freq_lo, &arg.freq_hi, obs);

    spec_getraw(opt_obs.cmf[0], &arg.band_count, &arg.bands, NULL);
    arg.obs = &opt_obs;
    arg.ill = illuminant;

    arg.obs_max = malloc(arg.band_count * sizeof(double));
    for (int i = 0; i < arg.band_count; i++) {
        arg.obs_max[i] = 0;
        for (int j = 0; j < 3; j++) {
            double e = spec_at(obs->cmf[j], arg.bands[i]);
            if (e > arg.obs_max[i]) {
                arg.obs_max[i] = e;
            }
        }
    }

    // normalize XYZ response so the spectrum shape finder gets to work on a predictable range
    arg.xyz_norm = 1.0 / ((xyz.v[0]+xyz.v[1]+xyz.v[2]) / 3.0);
    arg.xyz = (v3){{xyz.v[0]*arg.xyz_norm,xyz.v[1]*arg.xyz_norm,xyz.v[2]*arg.xyz_norm}};

    enum xyz2spec_mode best_mode = -1;
    struct {
        spectrum_t *spectrum;
        double peak;
        double sf;
    } result[XYZ2SPEC_MODE_LAST+1];
    memset(&result, 0, sizeof(result));
    for (enum xyz2spec_mode mode = XYZ2SPEC_MODE_FIRST; mode <= XYZ2SPEC_MODE_LAST; mode++) {
        init_xyz2spec_simplex(m, &arg, mode);
        min = run_xyz2spec_simplex(m, &arg);
        if (0) {
            if (arg.last) {
                FILE *stream = fopen("dump1/s.dat", "w");
                spectrum_t *s = spec_copy(arg.last);
                spec_scalar_multiply(s, 1.0 / arg.xyz_norm);
                spec_print(s, stream);
                free(s);
                fclose(stream);
                elog("first exact (%d) %f\n", mode, min);
                getchar();
            }
        }
        if (min <= max_good_de) {
            refine_xyz2spec_simplex(m, &arg);
            if ((int)best_mode == -1) {
                best_mode = mode;
            }
            result[mode].spectrum = spec_copy(arg.last);
            result[mode].peak = spec_max(arg.last);
            result[mode].sf = arg.last_sf;
        }
        if (0) {
            if (result[mode].spectrum) {
                FILE *stream = fopen("dump1/s.dat", "w");
                spectrum_t *s = spec_copy(arg.last);
                spec_scalar_multiply(s, 1.0 / arg.xyz_norm);
                spec_print(s, stream);
                free(s);
                fclose(stream);
                elog("mode %d refined exact\n", mode);
                getchar();
            } else {
                FILE *stream = fopen("dump1/s.dat", "w");
                spectrum_t *s = spec_copy(arg.last);
                spec_scalar_multiply(s, 1.0 / arg.xyz_norm);
                spec_print(s, stream);
                free(s);
                fclose(stream);
                elog("mode %d no solution\n", mode);
                getchar();
            }
        }
        if (mode+1 == XYZ2SPEC_MODE_FIRST_EXTREME && (int)best_mode != -1) {
            // only do the last types (designed for extremes) if we haven't found a match yet
            break;
        }
    }
    // gaussian bands: pick the one with lowest sf
    double sf = -1;
    for (enum xyz2spec_mode mode = XYZ2SPEC_MODE_GAUSSIAN_BANDS; mode <= XYZ2SPEC_MODE_GAUSSIAN_BANDS_OFFSET_RIGHT; mode++) {
        if (result[mode].spectrum != NULL && (result[mode].sf < sf || sf == -1)) {
            sf = result[mode].sf;
            best_mode = mode;
        }
    }
    if (result[XYZ2SPEC_MODE_GAUSSIAN].spectrum != NULL &&
        result[XYZ2SPEC_MODE_GAUSSIAN].peak * 0.8 < result[best_mode].peak)
    {
        best_mode = XYZ2SPEC_MODE_GAUSSIAN;
    }
    if (result[XYZ2SPEC_MODE_2X_GAUSSIAN].spectrum != NULL &&
        result[XYZ2SPEC_MODE_2X_GAUSSIAN].peak * 1.1 < result[best_mode].peak)
    {
        best_mode = XYZ2SPEC_MODE_2X_GAUSSIAN;
    }

    spectrum_t *best_spectrum = NULL;
    for (enum xyz2spec_mode mode = XYZ2SPEC_MODE_FIRST; mode <= XYZ2SPEC_MODE_LAST; mode++) {
        if (mode != best_mode) {
            free(result[mode].spectrum);
        } else {
            best_spectrum = result[mode].spectrum;
        }
    }
    if (best_spectrum != NULL) {
        // revert normalization
        spec_scalar_multiply(best_spectrum, 1.0 / arg.xyz_norm);
    }

    if (0) {
        if (best_spectrum) {
            FILE *stream = fopen("dump1/s.dat", "w");
            spec_print(best_spectrum, stream);
            fclose(stream);
        }
        elog("%sbest (%d)\n", best_spectrum == NULL ? "no " : "", best_mode);
    }
    free(opt_obs.cmf[0]);
    free(opt_obs.cmf[1]);
    free(opt_obs.cmf[2]);
    free(arg.bands);
    free(arg.obs_max);
    free(arg.last);
    return best_spectrum;
}
