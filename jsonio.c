/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <jsonio.h>
#include <lut.h>
#include <tifio.h>
#include <dngref.h>
#include <observers.h>
#include <nmsimplex.h>
#include <spectraldb.h>
#include <argyllio.h>
#include <wcompat.h>
#include <interp.h>
#include <gamut.h>
#include <look.h>
#include <elog.h>
#include <cJSON.h>

static uint8_t *
file2mem_nullterm(const char filename[],
                  const char mode[],
                  size_t *size)
{
    FILE *stream = utf8_fopen(filename, mode);
    if (stream == NULL) {
        elog("Could not open file \"%s\": %s\n", filename, strerror(errno));
        return NULL;
    }
    size_t bufsize = 4096;
    size_t datasize = 0, ret;
    uint8_t *buf = (uint8_t *)malloc(bufsize);
    while ((ret = fread(&buf[datasize], 1, bufsize - datasize, stream)) != 0) {
        datasize += ret;
        if (datasize == bufsize) {
            bufsize *= 2;
            buf = (uint8_t *)realloc(buf, bufsize);
        }
    }
    if (!feof(stream)) {
        fclose(stream);
        free(buf);
        elog("Failed to read file \"%s\"\n", filename);
        return NULL;
    }
    fclose(stream);
    buf = (uint8_t *)realloc(buf, datasize + 1);
    buf[datasize] = '\0';
    *size = datasize;
    return buf;
}

uint8_t *
fileio_file2mem(const char filename[],
                size_t *size)
{
    return file2mem_nullterm(filename, "rb", size);
}

char *
fileio_file2string(const char filename[])
{
    size_t size;
    return (char *)file2mem_nullterm(filename, "r", &size);
}

static bool
test_if_json(const char filename[])
{
    char *str = fileio_file2string(filename);
    if (str == NULL) {
        return false;
    }
    cJSON_Minify(str);
    cJSON *js = cJSON_Parse(str);
    free(str);
    if (js != NULL) {
        cJSON_Delete(js);
        return true;
    }
    return false;
}

static cJSON *
json_parse(const char filename[])
{
    char *str = fileio_file2string(filename);
    if (str == NULL) {
        return NULL;
    }
    cJSON_Minify(str); // cJSON requires comments to be removed before parsing
    cJSON *js = cJSON_Parse(str);
    free(str);
    if (js == NULL) {
        elog("Invalid JSON data in file \"%s\"\n", filename);
        return NULL;
    }
    return js;
}

static bool
get_and_test(cJSON **js_,
             cJSON *jsroot,
             const char name[],
             int type,
             bool is_required)
{
    cJSON *js = cJSON_GetObjectItem(jsroot, name);
    *js_ = NULL;
    if (js == NULL) {
        if (is_required) {
            elog("Required field \"%s\" is missing.\n", name);
            return false;
        }
        return true;
    }
    if (type >= 0 && js->type != type) {
        elog("Type mismatch for field \"%s\".\n", name);
        return false;
    }
    *js_ = js;
    return true;
}

static double *
json_numberarray(cJSON *js,
                 int *count_)
{
    if (js == NULL) {
        elog("Missing array.\n");
        return NULL;
    }
    if (js->type != cJSON_Array) {
        elog("Expected array, got other type.\n");
        return NULL;
    }
    cJSON *ji;
    int count = 0;
    for (ji = js->child; ji != NULL; ji = ji->next) {
        if (ji->type != cJSON_Number) {
            elog("Expected Number in array, got other type.\n");
            return NULL;
        }
        count++;
    }
    if (count == 0) {
        elog("Zero length array, expected more.\n");
        return NULL;
    }
    int i;
    double *a = malloc(sizeof(*a) * count);
    for (ji = js->child, i = 0; ji != NULL; ji = ji->next, i++) {
        a[i] = ji->valuedouble;
    }
    *count_ = count;
    return a;
}

static bool
get_fixed_numberarray(double *out,
                      cJSON *js,
                      int count)
{
    int c;
    double *a = json_numberarray(js, &c);
    if (a == NULL) {
        return false;
    }
    if (c != count) {
        free(a);
        elog("Error: expected array length %d, got %d.\n", count, c);
        return false;
    }
    memcpy(out, a, count * sizeof(a[0]));
    free(a);
    return true;
}

static double *
json_spectrum_bands_array(cJSON *js,
                          int *count_)
{
    int count;
    double *bands = json_numberarray(js, &count);
    if (bands == NULL) {
        elog("Could not parse JSON number array.\n");
        return NULL;
    }
    if (count != 3 || bands[2] > bands[1]) {
        for (int i = 0; i < count; i++) {
            if (bands[i] < 0 || (i > 0 && bands[i] <= bands[i-1])) {
                elog("Invalid bands for spectrum.\n");
                free(bands);
                return NULL;
            }
        }
        *count_ = count;
        return bands;
    }
    if (bands[0] >= bands[1]) {
        elog("Invalid bands for spectrum.\n");
        free(bands);
        return NULL;
    }
    int bspec[3] = { (int)bands[0], (int)bands[1], (int)bands[2] };
    free(bands);
    if (bspec[2] == 0 || (bspec[1] - bspec[0]) % bspec[2] != 0) {
        elog("Invalid bands for spectrum.\n");
        return NULL;
    }
    count = (bspec[1] - bspec[0]) / bspec[2] + 1;
    if (count <= 2 || count > 1000000) {
        elog("Invalid bands for spectrum.\n");
        return NULL;
    }
    double *ob = malloc(count * sizeof(ob[0]));
    for (int i = 0; i < count; i++) {
        ob[i] = bspec[0] + i * bspec[2];
    }
    ob[count-1] = bspec[1];
    if (ob[count-2] > ob[count-1]) {
        elog("Invalid bands for spectrum.\n");
        free(ob);
        return NULL;
    }
    *count_ = count;
    return ob;
}

struct dcam_ssf *
jsonio_dcam_ssf_parse(const char filename[],
                      bool parse_error_is_fatal)
{
    cJSON *js = json_parse(filename);
    if (js == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    // json_print(js);
    int band_count;
    double *band_f = NULL;
    spectrum_t *ssfs[3] = { NULL, NULL, NULL };
    cJSON *camera_name = cJSON_GetObjectItem(js, "camera_name");
    cJSON *bands = cJSON_GetObjectItem(js, "ssf_bands");
    cJSON *spec[3];
    spec[0] = cJSON_GetObjectItem(js, "red_ssf");
    spec[1] = cJSON_GetObjectItem(js, "green_ssf");
    spec[2] = cJSON_GetObjectItem(js, "blue_ssf");

    if (camera_name == NULL || camera_name->type != cJSON_String) {
        goto fail;
    }
    band_f = json_spectrum_bands_array(bands, &band_count);
    if (band_f == NULL) {
        goto fail;
    }
    for (int i = 0; i < 3; i++) {
        int amp_count;
        double *amp = json_numberarray(spec[i], &amp_count);
        if (amp == NULL || amp_count != band_count) {
            free(amp);
            goto fail;
        }
        ssfs[i] = spec_alloc1(band_f, amp, band_count);
        free(amp);
        if (ssfs[i] == NULL) {
            goto fail;
        }
    }
    free(band_f);
    double maxa = 0;
    for (int i = 0; i < 3; i++) {
        if (spec_max(ssfs[i]) > maxa) maxa = spec_max(ssfs[i]);
    }
    if (maxa != 1.0) {
        for (int i = 0; i < 3; i++) {
            spec_scalar_multiply(ssfs[i], 1.0 / maxa);
        }
    }
    struct dcam_ssf *ssf = calloc(1, sizeof(*ssf));
    ssf->camera_name = strdup(camera_name->valuestring);
    ssf->ssf.cmf[0] = ssfs[0];
    ssf->ssf.cmf[1] = ssfs[1];
    ssf->ssf.cmf[2] = ssfs[2];
    cJSON_Delete(js);
    return ssf;
fail:
    free(ssfs[0]);
    free(ssfs[1]);
    free(ssfs[2]);
    free(band_f);
    cJSON_Delete(js);
    elog("Could not parse Digital Camera SSF file \"%s\".\n", filename);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

void
jsonio_dcam_ssf_delete(struct dcam_ssf *ssf)
{
    if (ssf == NULL) {
        return;
    }
    free(ssf->camera_name);
    free(ssf->ssf.cmf[0]);
    free(ssf->ssf.cmf[1]);
    free(ssf->ssf.cmf[2]);
    free(ssf);
}

spectrum_t *
jsonio_spectrum_parse(const char filename[],
                      bool parse_error_is_fatal)
{
    cJSON *js = json_parse(filename);
    if (js == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    cJSON *js_bands = cJSON_GetObjectItem(js, "bands");
    cJSON *js_scale = cJSON_GetObjectItem(js, "scale");
    cJSON *js_spectrum = cJSON_GetObjectItem(js, "spectrum");
    int band_count;
    double *bands = NULL;
    double scale = 1.0;
    if (js_scale != NULL) {
        if (js_scale->type != cJSON_Number) {
            elog("\"scale\" is not a number.\n");
            goto fail;
        }
        scale = js_scale->valuedouble;
        if (scale <= 0) {
            elog("Bad \"scale\" value.\n");
            goto fail;
        }
    }
    bands = json_spectrum_bands_array(js_bands, &band_count);
    if (bands == NULL) {
        goto fail;
    }
    int amp_count;
    double *amp = json_numberarray(js_spectrum, &amp_count);
    if (amp == NULL || amp_count != band_count) {
        free(amp);
        goto fail;
    }
    if (scale != 1.0) {
        scale = 1.0 / scale;
        for (int i = 0; i < band_count; i++) {
            amp[i] *= scale;
        }
    }
    spectrum_t *s = spec_alloc1(bands, amp, band_count);
    if (s == NULL) {
        free(amp);
        goto fail;
    }
    free(amp);
    cJSON_Delete(js);
    return s;
fail:
    free(bands);
    elog("Could not parse Spectrum file \"%s\".\n", filename);
    cJSON_Delete(js);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

static double *
parse_dcp_tonecurve(cJSON *tonecurve,
                    int *tc_len)
{
    double *out_tc = NULL;
    int tonecurve_len = 0;
    for (cJSON *ji = tonecurve->child; ji != NULL; ji = ji->next) {
        if (ji->type != cJSON_Array) {
            elog("Expected array in array, got other type.\n");
            goto fail;
        }
        tonecurve_len++;
    }
    if (tonecurve_len < 2) {
        elog("Too short ProfileToneCurve.\n");
        goto fail;
    }
    out_tc = malloc(2 * tonecurve_len * sizeof(out_tc[0]));
    int count = 0;
    for (cJSON *ji = tonecurve->child; ji != NULL; ji = ji->next) {
        int i = 0;
        for (cJSON *ji1 = ji->child; ji1 != NULL; ji1 = ji1->next) {
            if (ji1->type != cJSON_Number) {
                elog("Expected Number in inner array, got other type.\n");
                goto fail;
            }
            if (i > 1) {
                elog("Expected length 2 of array element in ProfileToneCurve.\n");
                goto fail;
            }
            out_tc[2*count+i] = ji1->valuedouble;
            i++;
        }
        count++;
    }
    if (out_tc[0+0] != 0 || out_tc[0+1] != 0) {
        elog("Expected [0,0] as first element in ProfileToneCurve.\n");
        goto fail;
    }
    if (out_tc[2*(tonecurve_len-1)+0] != 1.0 || out_tc[2*(tonecurve_len-1)+1] != 1.0) {
        elog("Expected [1,1] as last element in ProfileToneCurve.\n");
        goto fail;
    }
    *tc_len = tonecurve_len;
    return out_tc;
fail:
    free(out_tc);
    *tc_len = 0;
    return NULL;
}

static double *
parse_single_trc(cJSON *trc,
                 int *count_)
{
    double *out_trc = NULL;
    int trc_len = 0;
    if (trc->type == cJSON_Number) {
        trc_len = 1;
        out_trc = malloc(sizeof(out_trc[0]));
        out_trc[0] = trc->valuedouble;
    } else if (trc->type == cJSON_Array) {
        trc_len = cJSON_GetArraySize(trc);
        if (trc_len < 2) {
            elog("Too short TRC array.\n");
            goto fail;
        }
        out_trc = malloc(trc_len * sizeof(out_trc[0]));
        if (!get_fixed_numberarray(out_trc, trc, trc_len)) {
            elog("Bad TRC array.\n");
            goto fail;
        }
    } else {
        elog("Expected array or number in TRC.\n");
        goto fail;
    }
    *count_ = trc_len;
    return out_trc;
fail:
    free(out_trc);
    *count_ = 0;
    return NULL;
}

static bool
parse_trc(cJSON *trc[3],
          double *out_trc[3],
          int *count_)
{
    out_trc[0] = NULL;
    out_trc[1] = NULL;
    out_trc[2] = NULL;
    int trc_len = 0;
    if (trc[0] != NULL || trc[1] != NULL || trc[2] != NULL) {
        if (trc[0] == NULL || trc[1] == NULL || trc[2] == NULL) {
            elog("Either all three or none of the TRCs must be set.\n");
            goto fail;
        }
        if (trc[0]->type == cJSON_Number || trc[1]->type == cJSON_Number || trc[2]->type == cJSON_Number) {
            if (trc[0]->type != cJSON_Number || trc[1]->type != cJSON_Number || trc[2]->type != cJSON_Number) {
                elog("TRC content mismatch.\n");
                goto fail;
            }
            trc_len = 1;
            for (int i = 0; i < 3; i++) {
                out_trc[i] = malloc(sizeof(out_trc[i][0]));
                out_trc[i][0] = trc[i]->valuedouble;
            }
        } else if (trc[0]->type == cJSON_Array || trc[1]->type == cJSON_Array || trc[2]->type == cJSON_Array) {
            if (trc[0]->type != cJSON_Array || trc[1]->type != cJSON_Array || trc[2]->type != cJSON_Array) {
                elog("TRC content mismatch.\n");
                goto fail;
            }
            trc_len = cJSON_GetArraySize(trc[0]);
            if (trc_len < 2) {
                elog("Too short TRC array.\n");
                goto fail;
            }
            for (int i = 0; i < 3; i++) {
                out_trc[i] = malloc(trc_len * sizeof(out_trc[i][0]));
                if (!get_fixed_numberarray(out_trc[i], trc[i], trc_len)) {
                    elog("Bad TRC array.\n");
                    goto fail;
                }
            }
        } else {
            elog("Expected array or number in TRC.\n");
            goto fail;
        }
    }
    *count_ = trc_len;
    return true;
fail:
    free(out_trc[0]);
    free(out_trc[1]);
    free(out_trc[2]);
    *count_ = 0;
    return false;
}

static bool
test_if_tiff(const char filename[])
{
    FILE *stream = utf8_fopen(filename, "rb");
    if (stream == NULL) {
        return false;
    }
    bool istiff = false;
    uint8_t bo[2];
    if (fread(bo, 1, 2, stream) == 2 && ((bo[0] == 'I' && bo[1] == 'I') || (bo[0] == 'M' && bo[1] == 'M'))) {
        istiff = true;
    }
    fclose(stream);
    return istiff;
}

bool
jsonio_transfer_function_parse(const char filename[],
                               double *trc[3],
                               int *count_,
                               bool parse_error_is_fatal)
{
    *count_ = 0;
    trc[0] = NULL;
    trc[1] = NULL;
    trc[2] = NULL;
    if (test_if_tiff(filename)) {
        bool success = tifio_get_transfer_function(filename, trc, count_, parse_error_is_fatal);
        if (!success) return false;
        if (*count_ == 0) {
            elog("Error: no transfer function (TIFFTAG_TRANSFERFUNCTION) found in TIFF \"%s\".\n", filename);
            if (parse_error_is_fatal) exit(EXIT_FAILURE);
            return false;
        }
        return true;
    }
    cJSON *js = json_parse(filename);
    if (js == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    cJSON *jstrc[3];
    if (!get_and_test(&jstrc[0], js, "RedTRC", -1, true)) goto fail;
    if (!get_and_test(&jstrc[1], js, "GreenTRC", -1, true)) goto fail;
    if (!get_and_test(&jstrc[2], js, "BlueTRC", -1, true)) goto fail;
    if (!parse_trc(jstrc, trc, count_)) {
        goto fail;
    }
    cJSON_Delete(js);
    return true;

fail:
    elog("Could not parse transfer function file \"%s\".\n", filename);
    cJSON_Delete(js);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return false;
}

static struct lookop_curve *
parse_lookop_curve(cJSON *js);

bool
jsonio_tonecurve_parse(const char filename[],
                       double **trc,
                       int *count_,
                       bool error_is_fatal)
{
    if (file_exists(filename)) {
        char s[256];
        FILE *stream = utf8_fopen(filename, "r");
        if (stream != NULL) {
            if (fgets(s, sizeof(s), stream) != NULL && (strstr(s, "Spline") == s || strstr(s, "Linear") == s)) {
                bool is_linear = (strstr(s, "Linear") == s);
                int max_handles = is_linear ? 100000 : 1000;
                double *inx, *iny;
                inx = malloc(max_handles * sizeof(inx[0]));
                iny = malloc(max_handles * sizeof(iny[0]));
                int count = 0;
                while (fgets(s, sizeof(s), stream) != NULL) {
                    float x, y;
                    if (sscanf(s, "%f %f", &x, &y) != 2) {
                        elog("Error: could not parse line in %s curve file \"%s\".\n", is_linear ? "Linear" : "Spline", filename);
                        free(inx);
                        free(iny);
                        fclose(stream);
                        if (error_is_fatal) exit(EXIT_FAILURE);
                        return false;
                    }
                    inx[count] = x;
                    iny[count] = y;
                    count++;
                    if (count == max_handles) {
                        elog("Error: too many lines in %s curve file \"%s\".\n", is_linear ? "Linear" : "Spline", filename);
                        free(inx);
                        free(iny);
                        fclose(stream);
                        if (error_is_fatal) exit(EXIT_FAILURE);
                        return false;
                    }
                }
                fclose(stream);
                if (count == 0) {
                    elog("Error: too few lines in %s curve file \"%s\".\n", is_linear ? "Linear" : "Spline", filename);
                    if (error_is_fatal) exit(EXIT_FAILURE);
                    free(inx);
                    free(iny);
                    return false;
                }
                if (inx[count-1] != 1.0) {
                    for (int i = 0; i < count; i++) {
                        iny[i] /= inx[count-1];
                        inx[i] /= inx[count-1];
                    }
                }
                const int tc_len = 65536;
                double *outx = malloc(sizeof(outx[0]) * tc_len);
                *trc = malloc(sizeof((*trc)[0]) * tc_len);
                for (int i = 0; i < tc_len; i++) {
                    outx[i] = (double)i / (tc_len - 1);
                    outx[i] = srgb_gamma_forward(outx[i]);
                }
                *count_ = tc_len;
                if (is_linear) {
                    linear_interpolation(inx, iny, count, outx, *trc, *count_);
                    reconstruct_tonecurve(*trc, *count_);
                } else {
                    cubic_spline(inx, iny, count, outx, *trc, *count_);
                }
                free(outx);
                free(inx);
                free(iny);
                for (int i = 0; i < tc_len; i++) {
                    (*trc)[i] = srgb_gamma_inverse((*trc)[i]);
                }
                return true;
            }
            fclose(stream);
        }
    }
    cJSON *js = json_parse(filename);
    if (js == NULL) {
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    bool did_parse = true;
    cJSON *jstrc[5], *tonecurve;
    if (!get_and_test(&tonecurve, js, "ProfileToneCurve", cJSON_Array, false)) goto fail;
    if (tonecurve != NULL) {
        double *tc = parse_dcp_tonecurve(tonecurve, count_);
        if (tc == NULL) goto fail;
        *trc = malloc(*count_ * sizeof((*trc)[0]));
        double step = tc[2] - tc[0];
        for (int i = 0; i < *count_; i++) {
            (*trc)[i] = tc[2*i+1];
            if (i > 0) {
                double step1 = tc[2*i+0] - tc[2*i-2];
                if (step1 / step > 1.05 || step / step1 > 1.05) {
                    elog("ProfileToneCurve has uneven spacing, not supported for tonecurves.\n");
                    free(tc);
                    goto fail;
                }
            }
        }
        free(tc);
        return true;
    }
    cJSON *curvtype;
    if (!get_and_test(&curvtype, js, "CurveType", cJSON_String, false)) goto fail;
    if (curvtype != NULL) {
        struct lookop_curve *lc = parse_lookop_curve(js);
        if (lc == NULL) goto fail;
        if (lc->type != LOOKOP_CURVE_SPLINE && lc->type != LOOKOP_CURVE_LINEAR && lc->type != LOOKOP_CURVE_ROUNDEDSTEP) {
            elog("The given curve type is not allowed for tone curves.\n");
            goto fail;
        }
        bool is_srgb = false;
        double gamma_num = 1.0;
        cJSON *gamma;
        if (!get_and_test(&gamma, js, "CurveGamma", -1, false)) goto fail;
        if (gamma == NULL) {
            gamma_num = 1.0;
        } else if (gamma->type == cJSON_String) {
            if (strcasecmp(gamma->valuestring, "sRGB") == 0) {
                is_srgb = true;
            } else {
                elog("Unknown \"CurveGamma\" \"%s\".\n", gamma->valuestring);
                goto fail;
            }
        } else if (gamma->type == cJSON_Number) {
            gamma_num = gamma->valuedouble;
        } else {
            elog("Type mismatch for \"CurveGamma\".\n");
            goto fail;
        }
        double inx[lc->handle_count];
        double iny[lc->handle_count];
        for (int i = 0; i < lc->handle_count; i++) {
            inx[i] = lc->handles[2*i+0];
            iny[i] = lc->handles[2*i+1];
        }
        const int tc_len = 65536;
        double *outx = malloc(sizeof(outx[0]) * tc_len);
        *trc = malloc(sizeof((*trc)[0]) * tc_len);
        for (int i = 0; i < tc_len; i++) {
            outx[i] = (double)i / (tc_len - 1);
            if (is_srgb || gamma_num != 1.0) {
                if (is_srgb) {
                    outx[i] = srgb_gamma_forward(outx[i]);
                } else {
                    outx[i] = pow(outx[i], 1.0 / gamma_num);
                }
            }
        }
        *count_ = tc_len;
        switch (lc->type) {
        case LOOKOP_CURVE_SPLINE:
            cubic_spline(inx, iny, lc->handle_count, outx, *trc, *count_);
            break;
        case LOOKOP_CURVE_LINEAR:
            linear_interpolation(inx, iny, lc->handle_count, outx, *trc, *count_);
            reconstruct_tonecurve(*trc, *count_);
            break;
        case LOOKOP_CURVE_ROUNDEDSTEP:
            rounded_linear_interpolation(inx, iny, lc->handle_count, outx, *trc, *count_);
            break;
        }
        free(outx);
        if (is_srgb || gamma_num != 1.0) {
            for (int i = 0; i < tc_len; i++) {
                if (is_srgb) {
                    (*trc)[i] = srgb_gamma_inverse((*trc)[i]);
                } else {
                    (*trc)[i] = pow((*trc)[i], gamma_num);
                }
            }
        }
        cJSON_Delete(js);
        free(lc);
        return true;
    }

    if (!get_and_test(&jstrc[0], js, "TRC", -1, false)) goto fail;
    if (!get_and_test(&jstrc[1], js, "GrayTRC", -1, false)) goto fail;
    if (!get_and_test(&jstrc[2], js, "GreenTRC", -1, false)) goto fail;
    if (!get_and_test(&jstrc[3], js, "RedTRC", -1, false)) goto fail;
    if (!get_and_test(&jstrc[4], js, "BlueTRC", -1, false)) goto fail;
    *trc = NULL;
    *count_ = 0;
    did_parse = false;
    for (size_t i = 0; i < sizeof(jstrc)/sizeof(jstrc[0]); i++) {
        if (jstrc[i] == NULL) {
            continue;
        }
        did_parse = true;
        *trc = parse_single_trc(jstrc[i], count_);
        break;
    }
fail:
    cJSON_Delete(js);
    if (*trc == NULL) {
        if (!did_parse) {
            elog("Error: did not find any \"TRC\" entry in \"%s\".\n", filename);
        }
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    return true;
}

void
jsonio_spectrum_print(FILE *stream,
                      const char indent[],
                      const char name[],
                      const spectrum_t *s)
{
    fprintf(stream,
            "%s\"%s\": {\n"
            "%s  \"bands\": [ %g",
            indent, name, indent, s->v[0].f);
    for (int i = 1; i < s->band_count; i++) {
        fprintf(stream, ", %g", s->v[i].f);
    }
    fprintf(stream,
            "],\n"
            "%s  \"spectrum\": [ %g",
            indent, s->v[0].e);
    for (int i = 1; i < s->band_count; i++) {
        fprintf(stream, ", %g", s->v[i].e);
    }
    fprintf(stream,
            "]\n"
            "%s}",
            indent);
}

spectrum_t *
jsonio_illuminant_parse(const char filename[],
                        bool parse_error_is_fatal)
{
    const spectrum_t *ill = NULL;
    enum exif_lightsource ls = exif_lightsource_aton(filename);
    if ((int)ls > 0) {
        ill = spectraldb_illuminant(ls);
    } else {
        ill = spectraldb_illuminant_byname(filename, false);
    }
    if (ill != NULL) {
        return spec_copy(ill);
    }
    if (file_exists(filename)) {
        if (test_if_json(filename)) {
            return jsonio_spectrum_parse(filename, parse_error_is_fatal);
        } else {
            // assume it's a CGATS file
            struct patch_set *ps = argyllio_ti3_parse(filename, parse_error_is_fatal, false);
            if (ps == NULL || ps->patch_count < 1 || ps->patch[0].spectrum == NULL) {
                target_delete_patch_set(ps);
                if (parse_error_is_fatal) {
                    elog("Error: unexpected content in \"%s\".\n", filename);
                    exit(EXIT_FAILURE);
                }
                return NULL;
            }
            spectrum_t *s = spec_copy(ps->patch[0].spectrum);
            spec_scalar_multiply(s, 1.0 / spec_at(s, 560));
            target_delete_patch_set(ps);
            return s;
        }
    }
    if (parse_error_is_fatal) {
        // run this to get error message
        spectraldb_illuminant_byname(filename, true);
        exit(EXIT_FAILURE);
    }
    return NULL;
}

static bool
matrix3x3_parse(m3x3 *m_,
                cJSON *mjs)
{
    if (mjs->type != cJSON_Array) {
        elog("Expected array type, got other.\n");
        goto fail;
    }
    int row = 0;
    m3x3 m;
    for (cJSON *ji = mjs->child; ji != NULL; ji = ji->next) {
        if (ji->type != cJSON_Array) {
            elog("Expected array in array, got other type.\n");
            goto fail;
        }
        if (++row > 3) goto dimension_fail;
        int col = 0;
        for (cJSON *ji1 = ji->child; ji1 != NULL; ji1 = ji1->next) {
            if (ji1->type != cJSON_Number) {
                elog("Expected Number in inner array, got other type.\n");
                goto fail;
            }
            if (++col > 3) goto dimension_fail;
            m.v[row-1][col-1] = ji1->valuedouble;
        }
        if (col != 3) goto dimension_fail;
    }
    if (row != 3) goto dimension_fail;
    *m_ = m;
    return true;
dimension_fail:
    elog("Expected 3x3 matrix, got other dimensions.\n");
    return false;
fail:
    return false;
}

m3x3 *
jsonio_3x3matrix_parse(const char filename[],
                       const char optional_matrix_name[],
                       const char optional_matrix_name2[],
                       bool parse_error_is_fatal)
{
    cJSON *js = json_parse(filename);
    if (js == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    cJSON *mjs = NULL;
    if (optional_matrix_name2 != NULL) {
        mjs = cJSON_GetObjectItem(js, optional_matrix_name2);
        if (mjs == NULL && optional_matrix_name != NULL) {
            mjs = cJSON_GetObjectItem(js, optional_matrix_name);
        }
    } else if (optional_matrix_name != NULL) {
        mjs = cJSON_GetObjectItem(js, optional_matrix_name);
    }
    if (mjs == NULL) {
        mjs = cJSON_GetObjectItem(js, "matrix");
        if (mjs == NULL) {
            mjs = js;
        }
    }
    m3x3 m;
    if (!matrix3x3_parse(&m, mjs)) {
        goto fail;
    }
    cJSON_Delete(js);
    m3x3 *m_ = malloc(sizeof(*m_));
    *m_ = m;
    return m_;

fail:
    elog("Could not parse Matrix file \"%s\".\n", filename);
    cJSON_Delete(js);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

void
jsonio_3x3matrix_print(FILE *stream,
                       const char indent[],
                       const char name[],
                       m3x3 m)
{
    fprintf(stream,
            "%s\"%s\": [\n"
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

static v3 *
parse_dcp_huesatmap(const uint32_t dims[],
                    cJSON *js)
{
    // { "HueDiv:"  0, "SatDiv:"  0, "ValDiv:"  0, "HueShift:" +0.000000, "SatScale:" 1.000000, "ValScale": 1.000000 },
    double dimsize = (double)dims[0] * (double)dims[1] * (double)dims[2];
    if (dimsize > 100 * 1024 * 1024 || dimsize < 1) {
        elog("Bad table dimensions.\n");
        return NULL;
    }
    unsigned divs[3] = { 0,0,0 };
    v3 *hsv = malloc(dims[0]*dims[1]*dims[2] * sizeof(*hsv));
    for (cJSON *ji = js->child; ji != NULL; ji = ji->next) {
        if (ji->type != cJSON_Object) {
            elog("Expected object array, got other type.\n");
            goto fail;
        }
        cJSON *fields[6];
        if (!get_and_test(&fields[0], ji, "HueDiv", cJSON_Number, true)) goto fail;
        if (!get_and_test(&fields[1], ji, "SatDiv", cJSON_Number, true)) goto fail;
        if (!get_and_test(&fields[2], ji, "ValDiv", cJSON_Number, true)) goto fail;
        if (!get_and_test(&fields[3], ji, "HueShift", cJSON_Number, true)) goto fail;
        if (!get_and_test(&fields[4], ji, "SatScale", cJSON_Number, true)) goto fail;
        if (!get_and_test(&fields[5], ji, "ValScale", cJSON_Number, true)) goto fail;
        const int k = (dims[0]*dims[1])*divs[2] + dims[1]*divs[0] + divs[1];
        for (int i = 0; i < 3; i++) {
            if (fields[i]->valueint != (int)divs[i]) {
                elog("Bad HueSatDiv indexing.\n");
                goto fail;
            }
            hsv[k].v[i] = fields[3+i]->valuedouble;
        }
        divs[2]++;
        if (divs[2] == dims[2]) {
            divs[2] = 0;
            divs[1]++;
            if (divs[1] == dims[1]) {
                divs[1] = 0;
                divs[0]++;
                if (divs[0] == dims[0] && ji->next != NULL) {
                    elog("Expected end of array.\n");
                    goto fail;
                }
            }
        }
    }
    return hsv;
fail:
    free(hsv);
    return NULL;
}

struct profio_dcp *
jsonio_dcp_parse(const char filename[],
                 cJSON *root,
                 bool parse_error_is_fatal)
{
    struct profio_dcp *dcp = NULL;
    cJSON *jsroot = root != NULL ? root : json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }

    cJSON *camera_name, *profile_name, *copyright, *embed_policy, *calibration_sig, *ill[2], *cm[2], *fm[2], *black_render, *baseline_exp;
    cJSON *look_dims, *look_enc, *hsm_dims, *hsm_enc, *hsm[2], *look, *tonecurve;
    if (!get_and_test(&camera_name, jsroot, "UniqueCameraModel", cJSON_String, true)) goto fail;
    if (!get_and_test(&profile_name, jsroot, "ProfileName", cJSON_String, false)) goto fail;
    if (!get_and_test(&copyright, jsroot, "ProfileCopyright", cJSON_String, false)) goto fail;
    if (!get_and_test(&embed_policy, jsroot, "ProfileEmbedPolicy", cJSON_String, false)) goto fail;
    if (!get_and_test(&calibration_sig, jsroot, "ProfileCalibrationSignature", cJSON_String, false)) goto fail;
    if (calibration_sig == NULL) {
        // old name
        if (!get_and_test(&calibration_sig, jsroot, "CalibrationSignature", cJSON_String, false)) goto fail;
    }
    if (!get_and_test(&ill[0], jsroot, "CalibrationIlluminant1", cJSON_String, true)) goto fail;
    if (!get_and_test(&ill[1], jsroot, "CalibrationIlluminant2", cJSON_String, false)) goto fail;
    if (!get_and_test(&cm[0], jsroot, "ColorMatrix1", cJSON_Array, true)) goto fail;
    if (!get_and_test(&cm[1], jsroot, "ColorMatrix2", cJSON_Array, ill[1] != NULL)) goto fail;
    if (!get_and_test(&fm[0], jsroot, "ForwardMatrix1", cJSON_Array, false)) goto fail;
    if (!get_and_test(&fm[1], jsroot, "ForwardMatrix2", cJSON_Array, fm[0] != NULL && ill[1] != NULL)) goto fail;
    if (!get_and_test(&black_render, jsroot, "DefaultBlackRender", cJSON_String, false)) goto fail;
    if (!get_and_test(&baseline_exp, jsroot, "BaselineExposureOffset", cJSON_Number, false)) goto fail;
    if (!get_and_test(&tonecurve, jsroot, "ProfileToneCurve", cJSON_Array, false)) goto fail;
    if (!get_and_test(&hsm[0], jsroot, "ProfileHueSatMap1", cJSON_Array, false)) goto fail;
    if (!get_and_test(&hsm[1], jsroot, "ProfileHueSatMap2", cJSON_Array, hsm[0] != NULL && ill[1] != NULL)) goto fail;
    if (!get_and_test(&hsm_dims, jsroot, "ProfileHueSatMapDims", cJSON_Array, hsm[0] != NULL)) goto fail;
    if (!get_and_test(&look, jsroot, "ProfileLookTable", cJSON_Array, false)) goto fail;
    if (!get_and_test(&look_dims, jsroot, "ProfileLookTableDims", cJSON_Array, look != NULL)) goto fail;
    if (!get_and_test(&hsm_enc, jsroot, "ProfileHueSatMapEncoding", cJSON_String, false)) goto fail;
    if (!get_and_test(&look_enc, jsroot, "ProfileLookTableEncoding", cJSON_String, false)) goto fail;

    dcp = calloc(1, sizeof(*dcp));
    dcp->camera_model = strdup(camera_name->valuestring);
    if (profile_name) dcp->profile_name = strdup(profile_name->valuestring);
    if (copyright) dcp->copyright = strdup(copyright->valuestring);
    if (calibration_sig) dcp->calibration_signature = strdup(calibration_sig->valuestring);
    dcp->ill_count = (ill[1] != NULL) ? 2 : 1;
    if (hsm_dims != NULL) {
        int count;
        double *array = json_numberarray(hsm_dims, &count);
        if (count != 3) {
            elog("Expected three integer array for ProfileHueSatMapDims.\n");
            free(array);
            goto fail;
        }
        dcp->hsmdims[0] = (int)array[0];
        dcp->hsmdims[1] = (int)array[1];
        dcp->hsmdims[2] = (int)array[2];
        free(array);
    }
    if (look_dims != NULL) {
        int count;
        double *array = json_numberarray(look_dims, &count);
        if (count != 3) {
            elog("Expected three integer array for ProfileLookTableDims.\n");
            free(array);
            goto fail;
        }
        dcp->lookdims[0] = (int)array[0];
        dcp->lookdims[1] = (int)array[1];
        dcp->lookdims[2] = (int)array[2];
        free(array);
    }
    for (int i = 0; i < dcp->ill_count; i++) {
        dcp->ill[i] = exif_lightsource_aton(ill[i]->valuestring);
        if ((int)dcp->ill[i] < 0) {
            elog("Unknown CalibrationIlluminant%d \"%s\".\n", i+1, ill[i]->valuestring);
            goto fail;
        }
        if (!matrix3x3_parse(&dcp->cm[i], cm[i])) {
            elog("Failed to parse ColorMatrix%d.\n", i+1);
            goto fail;
        }
        if (fm[i] != NULL) {
            dcp->has.fm = true;
            if (!matrix3x3_parse(&dcp->fm[i], fm[i])) {
                elog("Failed to parse ForwardMatrix%d.\n", i+1);
                goto fail;
            }
        }
        if (hsm[i] != NULL) {
            if ((dcp->huesatmap[i] = parse_dcp_huesatmap(dcp->hsmdims, hsm[i])) == NULL) {
                elog("Failed to parse ProfileHueSatMap%d.\n", i+1);
                goto fail;
            }
        }
    }
    if (look != NULL) {
        if ((dcp->looktable = parse_dcp_huesatmap(dcp->lookdims, look)) == NULL) {
            elog("Failed to parse ProfileLookTable.\n");
            goto fail;
        }
    }
    if (tonecurve) {
        dcp->tonecurve = parse_dcp_tonecurve(tonecurve, &dcp->tonecurve_len);
        if (dcp->tonecurve == NULL) {
            goto fail;
        }
    }
    if (look_enc) {
        dcp->has.looktable_srgb_gamma = true;
        if (strcasecmp(look_enc->valuestring, "sRGB") == 0) {
            dcp->looktable_srgb_gamma = true;
        } else if (strcasecmp(look_enc->valuestring, "Linear") == 0) {
            dcp->looktable_srgb_gamma = false;
        } else {
            elog("Unknown ProfileLookTableEncoding \"%s\".\n", look_enc->valuestring);
            goto fail;
        }
    }
    if (hsm_enc) {
        dcp->has.huesatmap_srgb_gamma = true;
        if (strcasecmp(hsm_enc->valuestring, "sRGB") == 0) {
            dcp->huesatmap_srgb_gamma = true;
        } else if (strcasecmp(hsm_enc->valuestring, "Linear") == 0) {
            dcp->huesatmap_srgb_gamma = false;
        } else {
            elog("Unknown ProfileHueSatMapEncoding \"%s\".\n", hsm_enc->valuestring);
            goto fail;
        }
    }
    if (baseline_exp) {
        dcp->has.baseline_exposure_offset = true;
        dcp->baseline_exposure_offset = baseline_exp->valuedouble;
    }
    if (black_render) {
        dcp->has.black_render_none = true;
        if (strcasecmp(black_render->valuestring, "None") == 0) {
            dcp->black_render_none = true;
        } else if (strcasecmp(black_render->valuestring, "Auto") == 0) {
            dcp->black_render_none = false;
        } else {
            elog("Unknown DefaultBlackRender \"%s\".\n", black_render->valuestring);
            goto fail;
        }
    }
    if (embed_policy) {
        dcp->has.embed_policy = true;
        if (strcasecmp(embed_policy->valuestring, "Allow copying") == 0 ||
            strcasecmp(embed_policy->valuestring, "AllowCopying") == 0)
        {
            dcp->embed_policy = pepAllowCopying;
        } else if (strcasecmp(embed_policy->valuestring, "Embed if used") == 0 ||
                   strcasecmp(embed_policy->valuestring, "EmbedIfUsed") == 0)
        {
            dcp->embed_policy = pepEmbedIfUsed;
        } else if (strcasecmp(embed_policy->valuestring, "Embed never") == 0 ||
                   strcasecmp(embed_policy->valuestring, "EmbedNever") == 0)
        {
            dcp->embed_policy = pepEmbedNever;
        } else if (strcasecmp(embed_policy->valuestring, "No restrictions") == 0 ||
                   strcasecmp(embed_policy->valuestring, "NoRestrictions") == 0)
        {
            dcp->embed_policy = pepNoRestrictions;
        } else {
            elog("Unknown ProfileEmbedPolicy \"%s\".\n", embed_policy->valuestring);
            goto fail;
        }
    }

    cJSON_Delete(jsroot);
    return dcp;

fail:
    if (filename != NULL) {
        elog("Could not parse DCP JSON file \"%s\".\n", filename);
    } else {
        elog("Could not parse DCP JSON.\n");
    }
    cJSON_Delete(jsroot);
    profio_dcp_delete(dcp);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

struct dcam_profile *
jsonio_profile_parse(const char filename[],
                     bool parse_error_is_fatal)
{
    struct dcam_profile *prof = NULL;
    cJSON *jsroot = json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    cJSON *cm, *fm, *lm, *fmwb, *ill, *lut, *ill_spec, *ill_wp, *camname;

    if (!get_and_test(&camname, jsroot, "CameraName", cJSON_String, false)) goto fail;
    if (!get_and_test(&ill, jsroot, "CalibrationIlluminant1", cJSON_String, true)) goto fail;
    if (!get_and_test(&ill_wp, jsroot, "CalibrationIlluminantWhitepoint1", cJSON_Array, false)) goto fail;
    if (!get_and_test(&ill_spec, jsroot, "CalibrationIlluminantSpectrum1", cJSON_Object, false)) goto fail;
    if (!get_and_test(&cm, jsroot, "ColorMatrix1", cJSON_Array, true)) goto fail;
    if (!get_and_test(&fm, jsroot, "ForwardMatrix1", cJSON_Array, true)) goto fail;
    if (!get_and_test(&lm, jsroot, "LUTMatrix1", cJSON_Array, false)) goto fail;
    if (!get_and_test(&fmwb, jsroot, "ForwardMatrixWhitebalance1", cJSON_Array, false)) goto fail;
    if (!get_and_test(&lut, jsroot, "ChromaLUT1", cJSON_Object, false)) goto fail;

    prof = calloc(1, sizeof(*prof));
    if (camname != NULL) {
        prof->camera_name = strdup(camname->valuestring);
    }
    prof->illuminant_name = exif_lightsource_aton(ill->valuestring);
    if ((int)prof->illuminant_name < 0) {
        elog("Unknown CalibrationIlluminant1 \"%s\".\n", ill->valuestring);
        goto fail;
    }
    const spectrum_t *ill_spectrum = spectraldb_illuminant(prof->illuminant_name);
    if (ill_wp == NULL) {
        if (ill_spectrum == NULL) {
            elog("Spectrum for CalibrationIlluminant1 \"%s\" is not known.\n", ill->valuestring);
            goto fail;
        }
        prof->illuminant_wp = spec2xyz(ill_spectrum, observer_get(OBSERVER_1931_2));
    } else {
        double a[3];
        if (!get_fixed_numberarray(a, ill_wp, 3)) {
            elog("Failed to parse CalibrationIlluminantWhitepoint1.\n");
            goto fail;
        }
        prof->illuminant_wp.v[0] = a[0];
        prof->illuminant_wp.v[1] = a[1];
        prof->illuminant_wp.v[2] = a[2];
    }
    if (fmwb == NULL) {
        prof->forward_matrix_wb = (v3){{1,1,1}};
    } else {
        double a[3];
        if (!get_fixed_numberarray(a, fmwb, 3)) {
            elog("Failed to parse ForwardMatrixWhitebalance1.\n");
            goto fail;
        }
        prof->forward_matrix_wb.v[0] = a[0];
        prof->forward_matrix_wb.v[1] = a[1];
        prof->forward_matrix_wb.v[2] = a[2];
    }
    if (!matrix3x3_parse(&prof->color_matrix, cm)) {
        elog("Failed to parse ColorMatrix1.\n");
        goto fail;
    }
    if (!matrix3x3_parse(&prof->forward_matrix, fm)) {
        elog("Failed to parse ForwardMatrix1.\n");
        goto fail;
    }
    if (lm != NULL) {
        if (!matrix3x3_parse(&prof->lut_matrix, lm)) {
            elog("Failed to parse LUTMatrix1.\n");
            goto fail;
        }
    } else {
        prof->lut_matrix = prof->forward_matrix;
    }
    if (ill_spec != NULL) {
        cJSON *jbands, *jspec;
        if (!get_and_test(&jbands, ill_spec, "bands", cJSON_Array, true)) goto fail;
        if (!get_and_test(&jspec, ill_spec, "spectrum", cJSON_Array, true)) goto fail;
        int band_count;
        double *band_f = json_spectrum_bands_array(jbands, &band_count);
        if (band_f == NULL) {
            goto fail;
        }
        int amp_count;
        double *amp = json_numberarray(jspec, &amp_count);
        if (amp == NULL || amp_count != band_count) {
            free(amp);
            free(band_f);
            goto fail;
        }
        prof->illuminant_spec = spec_alloc1(band_f, amp, band_count);
        free(band_f);
        free(amp);
        if (prof->illuminant_spec == NULL) {
            goto fail;
        }
    } else if (ill_spectrum != NULL) {
        prof->illuminant_spec = spec_copy(ill_spectrum);
    }

    if (lut != NULL) {
        cJSON *reg, *cf;
        cJSON *uv[3];
        if (!get_and_test(&reg, lut, "regularization", cJSON_Array, true)) goto fail;
        if (!get_and_test(&uv[0], lut, "uvL", cJSON_Array, true)) goto fail;
        if (!get_and_test(&uv[1], lut, "uvU", cJSON_Array, true)) goto fail;
        if (!get_and_test(&uv[2], lut, "uvV", cJSON_Array, true)) goto fail;
        if (!get_and_test(&cf, lut, "compressFactor", cJSON_Number, false)) goto fail;

        double regularization[3];
        {
            int count;
            double *relax = json_numberarray(reg, &count);
            if (relax == NULL || count != 3) {
                elog("Unexpected regularization content in LUT.\n");
                free(relax);
                goto fail;
            }
            regularization[0] = relax[0];
            regularization[1] = relax[1];
            regularization[2] = relax[2];
            free(relax);
        }
        double compress_factor = 0.7;
        if (cf != NULL) {
            compress_factor = cf->valuedouble;
        }
        int count[3];
        double *da[3] = { NULL, NULL, NULL };
        for (int i = 0; i < 3; i++) {
            if (regularization[i] < 0) {
                da[i] = NULL;
                count[i] = 0;
            } else {
                da[i] = json_numberarray(uv[i], &count[i]);
                if (da[i] == NULL) {
                    free(da[0]);free(da[1]);free(da[2]);
                    goto fail;
                }
            }
        }
        int ref_i = 0;
        for (int i = 0; i < 3; i++) {
            if (da[i] != NULL) {
                ref_i = i;
                break;
            }
        }
        for (int i = 0; i < 3; i++) {
            if (da[i] == NULL) continue;
            if (count[i] != count[ref_i] || count[i] % 3 != 0) {
                elog("Unexpected array length in LUT.\n");
                free(da[0]);free(da[1]);free(da[2]);
                goto fail;
            }
        }
        int cp_count = count[ref_i]/3;
        v3 *cp[3];
        for (int i = 0; i < 3; i++) {
            if (da[i] == NULL) {
                cp[i] = NULL;
                continue;
            }
            cp[i] = malloc(cp_count * sizeof(*cp[i]));
            for (int j = 0; j < cp_count; j++) {
                cp[i][j].v[0] = da[i][j*3+0];
                cp[i][j].v[1] = da[i][j*3+2]; // yes it is this order, 012 vs 021
                cp[i][j].v[2] = da[i][j*3+1];
            }
            free(da[i]);
        }
        prof->lut = chromalut_new_from_data(cp, cp_count, regularization, prof->forward_matrix, compress_factor);
        for (int i = 0; i < 3; i++) {
            free(cp[i]);
        }
        if (prof->lut == NULL) {
            elog("Invalid LUT data.\n");
            goto fail;
        }
    }

    cJSON_Delete(jsroot);
    return prof;

fail:
    elog("Could not parse camera profile in file \"%s\".\n", filename);
    cJSON_Delete(jsroot);
    if (prof != NULL) {
        free(prof->illuminant_spec);
        chromalut_delete(prof->lut);
        free(prof);
    }
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

struct profio_icc *
jsonio_icc_parse(const char filename[],
                     cJSON *root,
                 bool parse_error_is_fatal)
{
    struct profio_icc *icc = NULL;
    cJSON *jsroot = root != NULL ? root : json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }

    cJSON *description, *copyright, *fm, *trc[3], *point[2], *clut, *pcs;
    if (!get_and_test(&description, jsroot, "Description", cJSON_String, false)) goto fail;
    if (!get_and_test(&copyright, jsroot, "Copyright", cJSON_String, false)) goto fail;
    if (!get_and_test(&point[0], jsroot, "Whitepoint", cJSON_Array, false)) goto fail;
    if (!get_and_test(&point[1], jsroot, "Blackpoint", cJSON_Array, false)) goto fail;
    if (!get_and_test(&pcs, jsroot, "ProfileConnectionSpace", cJSON_String, false)) goto fail;
    if (!get_and_test(&fm, jsroot, "ForwardMatrix", cJSON_Array, false)) goto fail;
    if (!get_and_test(&trc[0], jsroot, "RedTRC", -1, false)) goto fail;
    if (!get_and_test(&trc[1], jsroot, "GreenTRC", -1, false)) goto fail;
    if (!get_and_test(&trc[2], jsroot, "BlueTRC", -1, false)) goto fail;
    if (!get_and_test(&clut, jsroot, "CLUT", cJSON_Object, false)) goto fail;

    icc = calloc(1, sizeof(*icc));
    if (description) icc->description = strdup(description->valuestring);
    if (copyright) icc->copyright = strdup(copyright->valuestring);
    for (int i = 0; i < 2; i++) {
        if (point[i] != NULL) {
            int count;
            double *array = json_numberarray(point[i], &count);
            if (count != 3) {
                elog("Expected three number array for %s.\n", i == 0 ? "Whitepoint" : "Blackpoint");
                free(array);
                goto fail;
            }
            v3 p;
            p.v[0] = array[0];
            p.v[1] = array[1];
            p.v[2] = array[2];
            free(array);
            if (i == 0) icc->wp = p; else icc->bp = p;
        }
    }
    if (fm != NULL) {
        if (!matrix3x3_parse(&icc->fm, fm)) {
            elog("Failed to parse ForwardMatrix.\n");
            goto fail;
        }
    }
    if (pcs != NULL) {
        if (strcasecmp(pcs->valuestring, "Lab") == 0) {
            icc->pcs_is_lab = true;
        } else if (strcasecmp(pcs->valuestring, "XYZ") != 0) {
            elog("Expected 'Lab' or 'XYZ' as ProfileConnectionSpace.\n");
            goto fail;
        }
    }
    if (!parse_trc(trc, icc->trc, &icc->trc_len)) {
        goto fail;
    }

    if (clut != NULL) {
        cJSON *im, *curv[2], *lt;
        if (!get_and_test(&im, clut, "InputMatrix", cJSON_Array, true)) goto fail;
        if (!get_and_test(&curv[0], clut, "InputCurves", cJSON_Array, true)) goto fail;
        if (!get_and_test(&curv[1], clut, "OutputCurves", cJSON_Array, true)) goto fail;
        if (!get_and_test(&lt, clut, "CLUT", cJSON_Array, true)) goto fail;

        int curv_sz[2] = { cJSON_GetArraySize(curv[0]), cJSON_GetArraySize(curv[1]) };
        if (curv_sz[0] < 2 || curv_sz[1] < 2) {
            elog("Bad Curves array length.\n");
            goto fail;
        }
        int clut_sz = cJSON_GetArraySize(lt);
        const int cs = (int)round(pow(clut_sz, 1.0/3.0));
        if (cs*cs*cs != clut_sz) {
            elog("Bad CLUT array length (%d).\n", clut_sz);
            goto fail;
        }

        icc->lut = calloc(1, sizeof(*icc->lut)+sizeof(icc->lut->data[0])*(3*curv_sz[0] + 3*curv_sz[1] + 3*cs*cs*cs));
        icc->lut->in_len = curv_sz[0];
        icc->lut->out_len = curv_sz[1];
        icc->lut->clut_side = cs;

        if (!matrix3x3_parse(&icc->lut->cm, im)) {
            elog("Failed to parse CLUT InputMatrix.\n");
            goto fail;
        }
        int base = 0;
        for (int i = 0; i < 2; i++) {
            double *cv[3];
            cv[0] = &icc->lut->data[base+0];
            cv[1] = &icc->lut->data[base+curv_sz[i]];
            cv[2] = &icc->lut->data[base+2*curv_sz[i]];
            base += 3 * curv_sz[i];
            int j = 0;
            for (cJSON *ji = curv[i]->child; ji != NULL; ji = ji->next) {
                double a[3];
                if (!get_fixed_numberarray(a, ji, 3)) {
                    elog("Bad inner array in Curves.\n");
                    goto fail;
                }
                cv[0][j] = a[0];
                cv[1][j] = a[1];
                cv[2][j] = a[2];
                j++;
            }
            memcpy(i == 0 ? icc->lut->in : icc->lut->out, cv, sizeof(cv));
        }

        icc->lut->clut = &icc->lut->data[base];
        int divs[3] = { 0, 0, 0 };
        for (cJSON *ji = lt->child; ji != NULL; ji = ji->next) {
            if (ji->type != cJSON_Object) {
                elog("Expected object array for CLUT, got other type.\n");
                goto fail;
            }
            cJSON *fields[6];
            if (!get_and_test(&fields[0], ji, "Div1", cJSON_Number, true)) goto fail;
            if (!get_and_test(&fields[1], ji, "Div2", cJSON_Number, true)) goto fail;
            if (!get_and_test(&fields[2], ji, "Div3", cJSON_Number, true)) goto fail;
            if (!get_and_test(&fields[3], ji, "Ch1", cJSON_Number, true)) goto fail;
            if (!get_and_test(&fields[4], ji, "Ch2", cJSON_Number, true)) goto fail;
            if (!get_and_test(&fields[5], ji, "Ch3", cJSON_Number, true)) goto fail;
            for (int i = 0; i < 3; i++) {
                if (fields[i]->valueint != divs[i]) {
                    elog("Bad Div indexing in CLUT.\n");
                    goto fail;
                }
                icc->lut->clut[3*(cs*cs*divs[0]+cs*divs[1]+divs[2])+i] = fields[3+i]->valuedouble;
            }
            divs[2]++;
            if (divs[2] == cs) {
                divs[2] = 0;
                divs[1]++;
                if (divs[1] == cs) {
                    divs[1] = 0;
                    divs[0]++;
                }
            }
        }
    }
    cJSON_Delete(jsroot);
    return icc;

fail:
    if (filename != NULL) {
        elog("Error: could not parse ICC JSON file \"%s\".\n", filename);
    } else {
        elog("Could not parse ICC JSON.\n");
    }
    cJSON_Delete(jsroot);
    profio_icc_delete(icc);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

static struct lookop_curve *
parse_lookop_curve(cJSON *js)
{
    struct lookop_curve *c = NULL;
    if (js->type != cJSON_Object) {
        elog("Expected Object type of \"Curve\", got other type.\n");
        goto fail;
    }
    double curve_max = 1.0;
    cJSON *type, *handles, *cmax;
    if (!get_and_test(&type, js, "CurveType", cJSON_String, true)) goto fail;
    if (!get_and_test(&handles, js, "CurveHandles", cJSON_Array, true)) goto fail;
    if (!get_and_test(&cmax, js, "CurveMax", cJSON_Number, false)) goto fail;
    int handle_count = cJSON_GetArraySize(handles);
    if (handle_count < 2) {
        elog("Too few handles in curve.\n");
        goto fail;
    }
    if (cmax != NULL) {
        curve_max = cmax->valuedouble;
        if (curve_max <= 0) {
            elog("Bad \"CurveMax\" value.\n");
            goto fail;
        }
    }
    c = calloc(1, sizeof(*c) + sizeof(c->handles[0]) * 2 * handle_count);
    if (strcasecmp(type->valuestring, "Spline") == 0) {
        c->type = LOOKOP_CURVE_SPLINE;
    } else if (strcasecmp(type->valuestring, "RoundedStep") == 0) {
        c->type = LOOKOP_CURVE_ROUNDEDSTEP;
    } else if (strcasecmp(type->valuestring, "Linear") == 0) {
        c->type = LOOKOP_CURVE_LINEAR;
    } else {
        elog("Unknown Curve type \"%s\".\n", type->valuestring);
        goto fail;
    }

    for (cJSON *ji = handles->child; ji != NULL; ji = ji->next) {
        if (ji->type != cJSON_Array) {
            elog("Expected array in array in Curve Handles, got other type.\n");
            goto fail;
        }
    }
    int count = 0;
    for (cJSON *ji = handles->child; ji != NULL; ji = ji->next) {
        int i = 0;
        for (cJSON *ji1 = ji->child; ji1 != NULL; ji1 = ji1->next) {
            if (ji1->type != cJSON_Number) {
                elog("Expected Number in Curve Handles inner array, got other type.\n");
                goto fail;
            }
            if (i > 1) {
                elog("Expected length 2 of array element in Curve Handles.\n");
                goto fail;
            }
            c->handles[2*count+i] = ji1->valuedouble / curve_max;
            i++;
        }
        count++;
    }
    for (int i = 1; i < handle_count; i++) {
        if (c->handles[2*(i-1)] > c->handles[2*i]) {
            elog("Expected X values of Curve Handles to be ordered from smallest to largest X.\n");
            goto fail;
        }
    }
    c->handle_count = handle_count;
    return c;

fail:
    free(c);
    return NULL;
}

static bool
parse_lookop_blend(struct lookop_vs *vs,
                   cJSON *js)
{
    if (js->type != cJSON_Object) {
        elog("Expected Object in Blend array, got other type.\n");
        return false;
    }
    cJSON *x, *xrange, *yscale, *curve;
    if (!get_and_test(&x, js, "X", cJSON_String, true)) return false;
    if (!get_and_test(&xrange, js, "XRange", cJSON_Array, false)) return false;
    if (!get_and_test(&yscale, js, "YScale", cJSON_Number, false)) return false;
    if (!get_and_test(&curve, js, "Curve", cJSON_Object, true)) return false;

    if (strcasecmp(x->valuestring, "Lightness") == 0) {
        vs->x_type = LOOKOP_VS_X_LIGHTNESS;
    } else if (strcasecmp(x->valuestring, "Chroma") == 0) {
        vs->x_type = LOOKOP_VS_X_CHROMA;
    } else if (strcasecmp(x->valuestring, "Hue") == 0) {
        vs->x_type = LOOKOP_VS_X_HUE;
    } else if (strcasecmp(x->valuestring, "HSL-Lightness") == 0) {
        vs->x_type = LOOKOP_VS_X_HSL_LIGHTNESS;
    } else if (strcasecmp(x->valuestring, "HSL-Saturation") == 0) {
        vs->x_type = LOOKOP_VS_X_HSL_SATURATION;
    } else if (strcasecmp(x->valuestring, "HSL-Hue") == 0) {
        vs->x_type = LOOKOP_VS_X_HSL_HUE;
    } else if (strcasecmp(x->valuestring, "HSV-Value") == 0) {
        vs->x_type = LOOKOP_VS_X_HSV_VALUE;
    } else if (strcasecmp(x->valuestring, "HSV-Saturation") == 0) {
        vs->x_type = LOOKOP_VS_X_HSV_SATURATION;
    } else if (strcasecmp(x->valuestring, "HSV-Hue") == 0) {
        vs->x_type = LOOKOP_VS_X_HSV_HUE;
    } else {
        elog("Unknown Blend curve X axis \"%s\".\n", x->valuestring);
        return false;
    }
    if (xrange != NULL) {
        if (!get_fixed_numberarray(vs->xrange, xrange, 2)) {
            elog("Bad Blend XRange array.\n");
            return false;
        }
        if (vs->xrange[0] >= vs->xrange[1]) {
            elog("Bad range specified in Blend XRange.\n");
            return false;
        }
    }
    if (yscale != NULL) {
        vs->yscale = yscale->valuedouble;
    } else {
        vs->yscale = 1.0;
    }
    if ((vs->curve = parse_lookop_curve(curve)) == NULL) return false;

    return true;
}

static bool
parse_lookop(struct lookop *lo,
             cJSON *js)
{
    if (js->type != cJSON_Object) {
        elog("Expected Object in LookOperators array, got other type.\n");
        return false;
    }
    cJSON *operator;
    if (!get_and_test(&operator, js, "Operator", cJSON_String, true)) return false;
    if (strcasecmp(operator->valuestring, "Stretch") == 0) {
        lo->type = LOOKOP_STRETCH;
    } else if (strcasecmp(operator->valuestring, "AddHue") == 0) {
        lo->type = LOOKOP_ADDHUE;
    } else if (strcasecmp(operator->valuestring, "AddChroma") == 0) {
        lo->type = LOOKOP_ADDCHROMA;
    } else if (strcasecmp(operator->valuestring, "AddLightness") == 0) {
        lo->type = LOOKOP_ADDLIGHTNESS;
    } else if (strcasecmp(operator->valuestring, "ScaleChroma") == 0) {
        lo->type = LOOKOP_SCALECHROMA;
    } else if (strcasecmp(operator->valuestring, "ScaleLightness") == 0) {
        lo->type = LOOKOP_SCALELIGHTNESS;
    } else if (strcasecmp(operator->valuestring, "SetTemperature") == 0) {
        lo->type = LOOKOP_SETTEMPERATURE;
    } else if (strcasecmp(operator->valuestring, "Curves") == 0) {
        lo->type = LOOKOP_CURVES;
    } else {
        elog("Unknown Look Operator \"%s\".\n", operator->valuestring);
        return false;
    }

    cJSON *blend = NULL, *blendinvert = NULL, *blendrgb = NULL;
    if (!get_and_test(&blend, js, "Blend", cJSON_Array, false)) return false;
    if (!get_and_test(&blendinvert, js, "BlendInvert", -1, false)) return false;
    if (!get_and_test(&blendrgb, js, "BlendRGB", -1, false)) return false;

    if (blend != NULL && cJSON_GetArraySize(blend) > 0) {
        lo->vs_count = cJSON_GetArraySize(blend);
        lo->vs = calloc(1, sizeof(lo->vs[0]) * lo->vs_count);
        int i = 0;
        for (cJSON *ji = blend->child; ji != NULL; ji = ji->next, i++) {
            if (!parse_lookop_blend(&lo->vs[i], ji)) {
                elog("Could not parse \"Blend\" element %d.\n", i + 1);
                return false;
            }
        }
    }
    if (blendinvert != NULL && blendinvert->type == cJSON_True) {
        lo->blend_invert = true;
    }
    if (blendrgb != NULL && blendrgb->type == cJSON_True) {
        lo->blend_rgb = true;
    }

    switch (lo->type) {
    case LOOKOP_STRETCH: {
        cJSON *stretch;
        if (!get_and_test(&stretch, js, "Stretch", cJSON_Array, true)) return false;
        int dim_count = cJSON_GetArraySize(stretch);
        if (dim_count <= 0 || dim_count > 3) {
            elog("Too few or too many \"Stretch\" dimensions.\n");
            return false;
        }
        lo->p.stretch.dim_count = dim_count;
        int i = 0;
        for (cJSON *ji = stretch->child; ji != NULL; ji = ji->next, i++) {
            if (ji->type != cJSON_Object) {
                elog("Expected Object in Stretch array, got other type.\n");
                return false;
            }
            cJSON *x, *xrange, *curve;
            if (!get_and_test(&x, ji, "X", cJSON_String, true)) return false;
            if (!get_and_test(&xrange, ji, "XRange", cJSON_Array, false)) return false;
            if (!get_and_test(&curve, ji, "Curve", cJSON_Object, true)) return false;

            if (strcasecmp(x->valuestring, "Lightness") == 0) {
                lo->p.stretch.dims[i] = LOOKOP_VS_X_LIGHTNESS;
            } else if (strcasecmp(x->valuestring, "Chroma") == 0) {
                lo->p.stretch.dims[i] = LOOKOP_VS_X_CHROMA;
            } else if (strcasecmp(x->valuestring, "Hue") == 0) {
                lo->p.stretch.dims[i] = LOOKOP_VS_X_HUE;
            } else {
                elog("Unknown Stretch curve axis \"%s\".\n", x->valuestring);
                return false;
            }
            for (int j = 0; j < i; j++) {
                if (lo->p.stretch.dims[j] == lo->p.stretch.dims[i]) {
                    elog("Duplicate Stretch curve axis \"%s\".\n", x->valuestring);
                    return false;
                }
            }
            if (xrange != NULL) {
                if (!get_fixed_numberarray(lo->p.stretch.ranges[i], xrange, 2)) {
                    elog("Bad Stretch XRange array.\n");
                    return false;
                }
                if (lo->p.stretch.ranges[i][0] >= lo->p.stretch.ranges[i][1]) {
                    elog("Bad range specified in Stretch XRange.\n");
                    return false;
                }
            }
            if ((lo->p.stretch.curves[i] = parse_lookop_curve(curve)) == NULL) return false;
            int hc = lo->p.stretch.curves[i]->handle_count;
            if (lo->p.stretch.curves[i]->handles[0] != 0 || lo->p.stretch.curves[i]->handles[1] != 0 ||
                lo->p.stretch.curves[i]->handles[2*hc-2] != 1 || lo->p.stretch.curves[i]->handles[2*hc-1] != 1)
            {
                elog("Stretch curve must be diagonal (start at 0,0 and end at 1,1)!\n");
                return false;
            }
        }
        break;
    }
    case LOOKOP_ADDHUE:
    case LOOKOP_ADDCHROMA:
    case LOOKOP_ADDLIGHTNESS:
    case LOOKOP_SCALECHROMA:
    case LOOKOP_SCALELIGHTNESS: {
        cJSON *value;
        if (!get_and_test(&value, js, "Value", cJSON_Number, true)) return false;
        lo->p.gen.value = value->valuedouble;
        break;
    }
    case LOOKOP_SETTEMPERATURE: {
        cJSON *temp;
        bool is_xy = false;
        if (!get_and_test(&temp, js, "TempTint", cJSON_Array, false)) return false;
        if (temp == NULL) {
            if (!get_and_test(&temp, js, "XY", cJSON_Array, false)) return false;
            if (temp == NULL) {
                elog("Either \"TempTint\" or \"XY\" must be specified for the \"SetTemperature\" Look Operator.\n");
                return false;
            }
            is_xy = true;
        }
        double t2[2];
        if (!get_fixed_numberarray(t2, temp, 2)) {
            elog("Bad SetTemperature value array.\n");
            return false;
        }
        if (!is_xy) {
            v3 wp = dcp_d50();
            v3 xyY = xyz2xyY(wp);
            double xy[2] = { xyY.v[0], xyY.v[0] };
            dcp_temp2xy(t2[0], t2[1], xy);
            lo->p.settemp.xy[0] = xy[0];
            lo->p.settemp.xy[1] = xy[1];
        } else {
            lo->p.settemp.xy[0] = t2[0];
            lo->p.settemp.xy[1] = t2[1];
        }
        break;
    }
    case LOOKOP_CURVES: {
        cJSON *gamma, *curves, *kl;
        if (!get_and_test(&curves, js, "Curves", cJSON_Array, true)) return false;
        if (!get_and_test(&gamma, js, "CurvesGamma", -1, true)) return false;
        if (!get_and_test(&kl, js, "KeepLightness", -1, false)) return false;
        if (gamma->type == cJSON_String) {
            if (strcasecmp(gamma->valuestring, "sRGB") == 0) {
                lo->p.curves.is_srgb_gamma = true;
            } else {
                elog("Unknown \"CurvesGamma\" \"%s\".\n", gamma->valuestring);
                return false;
            }
        } else if (gamma->type == cJSON_Number) {
            lo->p.curves.gamma = gamma->valuedouble;
        } else {
            elog("Type mismatch for \"CurvesGamma\".\n");
            return false;
        }

        lo->p.curves.keep_lightness = (kl != NULL && kl->type == cJSON_True);

        if (cJSON_GetArraySize(curves) != 3) {
            elog("Bad \"Curves\" count (should be 3).\n");
            return false;
        }
        int i = 0;
        for (cJSON *ji = curves->child; ji != NULL; ji = ji->next, i++) {
            lo->p.curves.curves[i] = parse_lookop_curve(ji);
            if (lo->p.curves.curves[i] == NULL) {
                return false;
            }
        }
        break;
    }
    }
    return true;
}

static bool
parse_gamut_matrix(m3x3 *m_,
                   cJSON *gam)
{
    if (gam->type == cJSON_String) {
        if (strcasecmp(gam->valuestring, "sRGB") == 0) {
            *m_ = xyz2rgb_srgb;
        } else if (strcasecmp(gam->valuestring, "AdobeRGB") == 0) {
            *m_ = xyz2rgb_adobergb;
        } else if (strcasecmp(gam->valuestring, "ProphotoRGB") == 0) {
            *m_ = dcp_xyz_D50_to_prophoto_rgb, look_wp();
        } else {
            elog("Unknown gamut name \"%s\".\n", gam->valuestring);
            return false;
        }
    } else {
        if (!matrix3x3_parse(m_, gam)) {
            return false;
        }
    }
    *m_ = gamut_remap_xyz2rgb(*m_, look_wp());
    return true;
}

struct tone_rep_op_config *
jsonio_ntro_conf_parse(const char filename[],
                       bool parse_error_is_fatal)
{
    struct tone_rep_op_config *conf = NULL;
    cJSON *jsroot = json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    cJSON *chroma_scaling, *saturated, *shadows, *rolloff, *curve, *gamut, *lookop;
    if (!get_and_test(&chroma_scaling, jsroot, "ChromaScaling", cJSON_Number, true)) goto fail;
    if (!get_and_test(&saturated, jsroot, "Saturated", cJSON_Object, true)) goto fail;
    if (!get_and_test(&shadows, jsroot, "Shadows", cJSON_Object, true)) goto fail;
    if (!get_and_test(&curve, jsroot, "Curve", cJSON_Object, true)) goto fail;
    if (!get_and_test(&rolloff, jsroot, "Rolloff", cJSON_Object, true)) goto fail;
    if (!get_and_test(&gamut, jsroot, "GamutCompression", cJSON_Object, false)) goto fail;
    if (!get_and_test(&lookop, jsroot, "LookOperators", cJSON_Array, false)) goto fail;

    cJSON *crv_kf, *crv_lo, *crv_hi;
    if (!get_and_test(&crv_kf, curve, "KeepFactor", cJSON_Number, true)) goto fail;
    if (!get_and_test(&crv_lo, curve, "LowChroma", cJSON_Number, true)) goto fail;
    if (!get_and_test(&crv_hi, curve, "HighChroma", cJSON_Number, true)) goto fail;

    cJSON *sat_af, *sat_lo, *sat_hi;
    if (!get_and_test(&sat_af, saturated, "AdjustFactor", cJSON_Number, true)) goto fail;
    if (!get_and_test(&sat_lo, saturated, "LowChroma", cJSON_Number, true)) goto fail;
    if (!get_and_test(&sat_hi, saturated, "HighChroma", cJSON_Number, true)) goto fail;

    cJSON *shw_af, *shw_afhc, *shw_lo, *shw_hi, *shw_loc, *shw_hic;
    if (!get_and_test(&shw_af, shadows, "AdjustFactor", cJSON_Number, true)) goto fail;
    if (!get_and_test(&shw_afhc, shadows, "AdjustFactorHighChroma", cJSON_Number, true)) goto fail;
    if (!get_and_test(&shw_lo, shadows, "LowLightness", cJSON_Number, true)) goto fail;
    if (!get_and_test(&shw_hi, shadows, "HighLightness", cJSON_Number, true)) goto fail;
    if (!get_and_test(&shw_loc, shadows, "LowChroma", cJSON_Number, true)) goto fail;
    if (!get_and_test(&shw_hic, shadows, "HighChroma", cJSON_Number, true)) goto fail;

    cJSON *rlf_kf, *rlf_lo, *rlf_hi;
    if (!get_and_test(&rlf_kf, rolloff, "KeepFactorHueCurve", cJSON_Object, true)) goto fail;
    if (!get_and_test(&rlf_lo, rolloff, "LowSatScaleHueCurve", cJSON_Object, true)) goto fail;
    if (!get_and_test(&rlf_hi, rolloff, "HighSatScaleHueCurve", cJSON_Object, true)) goto fail;

    conf = calloc(1, sizeof(*conf));

    conf->chroma_scaling = chroma_scaling->valuedouble;

    conf->curve.keep_factor = crv_kf->valuedouble;
    conf->curve.low_chroma = crv_lo->valuedouble;
    conf->curve.high_chroma = crv_hi->valuedouble;

    conf->saturated.adjust_factor = sat_af->valuedouble;
    conf->saturated.low_chroma = sat_lo->valuedouble;
    conf->saturated.high_chroma = sat_hi->valuedouble;

    conf->shadows.adjust_factor = shw_af->valuedouble;
    conf->shadows.adjust_factor_high_chroma = shw_afhc->valuedouble;
    conf->shadows.low_lightness = shw_lo->valuedouble;
    conf->shadows.high_lightness = shw_hi->valuedouble;
    conf->shadows.low_chroma = shw_loc->valuedouble;
    conf->shadows.high_chroma = shw_hic->valuedouble;

    if ((conf->rolloff.keep_factor = parse_lookop_curve(rlf_kf)) == NULL) goto fail;
    if ((conf->rolloff.low_satscale = parse_lookop_curve(rlf_lo)) == NULL) goto fail;
    if ((conf->rolloff.high_satscale = parse_lookop_curve(rlf_hi)) == NULL) goto fail;

    if (gamut != NULL) {
        cJSON *ig, *igcs, *dg, *dgcs, *cg, *cgcs, *hilim, *lolim, *preproc;
        if (!get_and_test(&ig, gamut, "InnerGamut", -1, true)) goto fail;
        if (!get_and_test(&igcs, gamut, "InnerGamutChromaScale", cJSON_Number, true)) goto fail;
        if (!get_and_test(&dg, gamut, "DestinationGamut", -1, true)) goto fail;
        if (!get_and_test(&dgcs, gamut, "DestinationGamutChromaScale", cJSON_Number, true)) goto fail;
        if (!get_and_test(&cg, gamut, "ClipGamut", -1, true)) goto fail;
        if (!get_and_test(&cgcs, gamut, "ClipGamutChromaScale", cJSON_Number, true)) goto fail;
        if (!get_and_test(&hilim, gamut, "RGBChannelHighLimit", cJSON_Number, true)) goto fail;
        if (!get_and_test(&lolim, gamut, "RGBChannelLowLimit", cJSON_Number, true)) goto fail;
        if (!get_and_test(&preproc, gamut, "HSVPreprocessing", cJSON_Object, false)) goto fail;
        conf->gamut.hi_rgb_lim = hilim->valuedouble;
        conf->gamut.lo_rgb_lim = lolim->valuedouble;
        for (int i = 0; i < 3; i++) {
            struct gamut_spec *gs = i == 0 ? &conf->gamut.inner : i == 1 ? &conf->gamut.dest : &conf->gamut.clip;
            cJSON *gam = i == 0 ? ig : i == 1 ? dg : cg;
            if (gam->type == cJSON_String) {
                if (strcasecmp(gam->valuestring, "Locus") == 0) {
                    gs->locus_enabled = true;
                } else if (strcasecmp(gam->valuestring, "ProphotoRGB-Locus-Intersection") == 0) {
                    gs->locus_enabled = true;
                    gs->xyz2rgb_enabled = true;
                    gs->xyz2rgb = gamut_remap_xyz2rgb(dcp_xyz_D50_to_prophoto_rgb, look_wp());
                } else {
                    gs->xyz2rgb_enabled = true;
                    if (!parse_gamut_matrix(&gs->xyz2rgb, gam)) goto fail;
                }
            } else {
                gs->xyz2rgb_enabled = true;
                if (!parse_gamut_matrix(&gs->xyz2rgb, gam)) goto fail;
            }
            gs->chroma_scale = i == 0 ? igcs->valuedouble : i == 1 ? dgcs->valuedouble : cgcs->valuedouble;
        }
        conf->gamut_enabled = true;

        if (preproc != NULL) {
            cJSON *tll, *tol, *tcl, *bll, *bol, *bcl, *tl, *bl, *br;
            if (!get_and_test(&tll, preproc, "TopLinearLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&tol, preproc, "TopOutputLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&tcl, preproc, "TopCutoffLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&bll, preproc, "BottomLinearLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&bol, preproc, "BottomOutputLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&bcl, preproc, "BottomCutoffLimit", cJSON_Number, true)) goto fail;
            if (!get_and_test(&tl, preproc, "TopRGBMinLift", cJSON_Number, true)) goto fail;
            if (!get_and_test(&bl, preproc, "BottomRGBMaxLift", cJSON_Number, true)) goto fail;
            if (!get_and_test(&br, preproc, "InnerGamutSatBlendRange", cJSON_Array, true)) goto fail;

            conf->gamut.hsv.top.linear_limit = tll->valuedouble;
            conf->gamut.hsv.top.output_limit = tol->valuedouble;
            conf->gamut.hsv.top.cutoff_limit = tcl->valuedouble;
            conf->gamut.hsv.top.min_lift = tl->valuedouble;

            conf->gamut.hsv.bottom.linear_limit = bll->valuedouble;
            conf->gamut.hsv.bottom.output_limit = bol->valuedouble;
            conf->gamut.hsv.bottom.cutoff_limit = bcl->valuedouble;
            conf->gamut.hsv.bottom.max_lift = bl->valuedouble;
            double a[2];
            if (!get_fixed_numberarray(a, br, 2)) goto fail;
            conf->gamut.hsv.lo_sat_blend_lim = a[0];
            conf->gamut.hsv.hi_sat_blend_lim = a[1];
            conf->gamut.hsv.enabled = true;
        }
    }

    if (lookop != NULL && cJSON_GetArraySize(lookop) > 0) {
        conf->lookop_count = cJSON_GetArraySize(lookop);
        conf->lookops = calloc(1, sizeof(conf->lookops[0]) * conf->lookop_count);
        int i = 0;
        for (cJSON *ji = lookop->child; ji != NULL; ji = ji->next, i++) {
            if (!parse_lookop(&conf->lookops[i], ji)) {
                goto fail;
            }
        }
    }

    cJSON_Delete(jsroot);

    return conf;
fail:
    elog("Error: could not parse tone reproduction configuration file \"%s\".\n", filename);
    cJSON_Delete(jsroot);
    jsonio_ntro_delete(conf);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

void
jsonio_ntro_delete(struct tone_rep_op_config *conf)
{
    if (conf == NULL) {
        return;
    }
    if (conf->lookops != NULL) {
        for (int i = 0; i < conf->lookop_count; i++) {
            for (int j = 0; j < conf->lookops[i].vs_count; j++) {
                free(conf->lookops[i].vs[j].curve);
            }
            free(conf->lookops[i].vs);
            if (conf->lookops[i].type == LOOKOP_STRETCH) {
                free(conf->lookops[i].p.stretch.curves[0]);
                free(conf->lookops[i].p.stretch.curves[1]);
                free(conf->lookops[i].p.stretch.curves[2]);
            }
            if (conf->lookops[i].type == LOOKOP_CURVES) {
                free(conf->lookops[i].p.curves.curves[0]);
                free(conf->lookops[i].p.curves.curves[1]);
                free(conf->lookops[i].p.curves.curves[2]);
            }
        }
        free(conf->lookops);
    }
    free(conf->rolloff.keep_factor);
    free(conf->rolloff.low_satscale);
    free(conf->rolloff.high_satscale);
    free(conf);
}

struct target_adjustment_config *
jsonio_target_adjustment_parse(const char filename[],
                               bool parse_error_is_fatal)
{
    struct target_adjustment_config *conf = NULL;
    cJSON *jsroot = json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }

    cJSON *global, *patches;
    if (!get_and_test(&global, jsroot, "GlobalAdjustments", cJSON_Object, false)) goto fail;
    if (!get_and_test(&patches, jsroot, "PatchAdjustments", cJSON_Array, false)) goto fail;

    int patch_count = 0;
    if (patches != NULL) {
        patch_count = cJSON_GetArraySize(patches);
    }

    conf = calloc(1, sizeof(*conf) + patch_count * sizeof(conf->patch[0]));
    conf->global.scale_chroma = 1.0;
    conf->global.scale_lightness = 1.0;

    if (global != NULL) {
        cJSON *lph, *cph, *hph, *sl, *sc;
        if (!get_and_test(&sc, global, "ScaleChroma", cJSON_Number, false)) goto fail;
        if (!get_and_test(&sl, global, "ScaleLightness", cJSON_Number, false)) goto fail;
        if (!get_and_test(&lph, global, "ScaleLightnessPerHue", cJSON_Object, false)) goto fail;
        if (!get_and_test(&cph, global, "ScaleChromaPerHue", cJSON_Object, false)) goto fail;
        if (!get_and_test(&hph, global, "AddHuePerHue", cJSON_Object, false)) goto fail;

        if (lph != NULL && (conf->global.scale_lightness_per_hue = parse_lookop_curve(lph)) == NULL) goto fail;
        if (cph != NULL && (conf->global.scale_chroma_per_hue = parse_lookop_curve(cph)) == NULL) goto fail;
        if (hph != NULL && (conf->global.add_hue_per_hue = parse_lookop_curve(hph)) == NULL) goto fail;
        if (sl != NULL) conf->global.scale_lightness = sl->valuedouble;
        if (sc != NULL) conf->global.scale_chroma = sc->valuedouble;
    }

    if (patches != NULL) {
        conf->patch_count = cJSON_GetArraySize(patches);
        int i = 0;
        for (cJSON *js = patches->child; js != NULL; js = js->next, i++) {
            if (js->type != cJSON_Object) {
                elog("Expected Object in PatchAdjustments array, got other type.\n");
                goto fail;
            }
            conf->patch[i].scale_rgb[0] = 1;
            conf->patch[i].scale_rgb[1] = 1;
            conf->patch[i].scale_rgb[2] = 1;
            conf->patch[i].adjust_jch[0] = 1;
            conf->patch[i].adjust_jch[1] = 1;
            conf->patch[i].adjust_jch[2] = 0;
            cJSON *class, *name, *scale_rgb, *adjust_jch;
            if (!get_and_test(&name, js, "Name", cJSON_String, true)) goto fail;
            if (!get_and_test(&class, js, "Class", cJSON_String, false)) goto fail;
            if (!get_and_test(&scale_rgb, js, "ScaleRGB", cJSON_Array, false)) goto fail;
            if (!get_and_test(&adjust_jch, js, "AdjustJCh", cJSON_Array, false)) goto fail;
            if (scale_rgb != NULL && !get_fixed_numberarray(conf->patch[i].scale_rgb, scale_rgb, 3)) goto fail;
            if (adjust_jch != NULL && !get_fixed_numberarray(conf->patch[i].adjust_jch, adjust_jch, 3)) goto fail;
            strncpy(conf->patch[i].name, name->valuestring, sizeof(conf->patch[i].name)-1);
            if (class != NULL) strncpy(conf->patch[i].class, class->valuestring, sizeof(conf->patch[i].class)-1);
        }
    }

    cJSON_Delete(jsroot);
    return conf;

fail:
    elog("Error: could not parse target adjustment configuration file \"%s\".\n", filename);
    cJSON_Delete(jsroot);
    jsonio_target_adjustment_delete(conf);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}

void
jsonio_target_adjustment_delete(struct target_adjustment_config *conf)
{
    if (conf == NULL) {
        return;
    }
    free(conf->global.scale_chroma_per_hue);
    free(conf->global.scale_lightness_per_hue);
    free(conf->global.add_hue_per_hue);
    free(conf);
}

struct test_chart_spec *
jsonio_test_chart_spec_parse(const char filename[],
                             bool parse_error_is_fatal)
{
    struct test_chart_spec *spec = NULL;
    cJSON *jsroot = json_parse(filename);
    if (jsroot == NULL) {
        if (parse_error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    cJSON *rows, *cols, *rh, *cw, *co, *wp, *bp, *gp;
    if (!get_and_test(&rows, jsroot, "RowCount", cJSON_Number, true)) goto fail;
    if (!get_and_test(&cols, jsroot, "ColCount", cJSON_Number, true)) goto fail;
    if (!get_and_test(&rh, jsroot, "RelativeRowHeight", cJSON_Number, false)) goto fail;
    if (!get_and_test(&cw, jsroot, "RelativeColWidth", cJSON_Number, false)) goto fail;
    if (!get_and_test(&co, jsroot, "HasColOffset", -1, false)) goto fail;
    if (!get_and_test(&wp, jsroot, "WhitePatches", cJSON_Array, false)) goto fail;
    if (!get_and_test(&bp, jsroot, "BlackPatches", cJSON_Array, false)) goto fail;
    if (!get_and_test(&gp, jsroot, "MiddleGrayPatches", cJSON_Array, false)) goto fail;

    int col_count = cols->valueint;
    int row_count = rows->valueint;
    int patch_count = col_count * row_count;
    if (patch_count <= 0 || patch_count >= 1000000) {
        elog("Bad row/col count.\n");
        goto fail;
    }
    spec = calloc(1, sizeof(*spec) + patch_count * sizeof(spec->data[0]));
    spec->gray_specified = true;
    spec->patch_count = patch_count;
    spec->col_count = col_count;
    spec->layout.row_count = row_count;
    spec->layout.has_col_offset = (co != NULL && co->type == cJSON_True);
    spec->layout.row_height = (rh != NULL) ? rh->valuedouble : 1;
    spec->layout.col_width = (cw != NULL) ? cw->valuedouble : 1;
    spec->black_count = (bp != NULL) ? cJSON_GetArraySize(bp) : 0;
    spec->white_count = (wp != NULL) ? cJSON_GetArraySize(wp) : 0;
    spec->gray_count = (gp != NULL) ? cJSON_GetArraySize(gp) : 0;
    if (spec->black_count + spec->white_count + spec->gray_count > patch_count) {
        elog("Too many black/white/gray patches to fit on chart.\n");
        goto fail;
    }
    spec->black = (bp != NULL) ? &spec->data[0] : NULL;
    spec->white = (wp != NULL) ? &spec->data[spec->black_count] : NULL;
    spec->gray = (gp != NULL) ? &spec->data[spec->black_count+spec->white_count] : NULL;
    for (int i = 0; i < 3; i++) {
        cJSON *array = i == 0 ? wp : i == 1 ? bp : gp;
        struct test_chart_patch *patch = i == 0 ? spec->white : i == 1 ? spec->black : spec->gray;
        if (array == NULL || patch == NULL) {
            continue;
        }
        int k = 0;
        for (cJSON *js = array->child; js != NULL; js = js->next, k++) {
            if (js->type != cJSON_String) {
                elog("Expected string patch name.\n");
                goto fail;
            }
            if (strlen(js->valuestring) > PATCH_STRID_SIZE-1) {
                elog("Too long patch name.\n");
                goto fail;
            }
            if (strlen(js->valuestring) == 0) {
                elog("Empty patch name.\n");
                goto fail;
            }
            strcpy(patch[k].str_id, js->valuestring);
            patch[k].is_white = i == 0;
            patch[k].is_black = i == 1;
            patch[k].is_gray = i == 2;
            patch[k].idx = -1;
            patch[k].r = -1;
            patch[k].c = -1;
        }
    }
    for (int i = 0; i < spec->black_count + spec->white_count + spec->gray_count; i++) {
        for (int j = 0; j < spec->black_count + spec->white_count + spec->gray_count; j++) {
            if (i != j && strcasecmp(spec->data[i].str_id, spec->data[j].str_id) == 0) {
                elog("Duplicate patch name \"%s\".\n", spec->data[j].str_id);
                goto fail;
            }
        }
    }
    cJSON_Delete(jsroot);
    return spec;

fail:
    elog("Error: could not parse test chart layout file \"%s\".\n", filename);
    cJSON_Delete(jsroot);
    free(spec);
    if (parse_error_is_fatal) exit(EXIT_FAILURE);
    return NULL;
}
