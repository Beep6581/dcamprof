/*
 * (c) Copyright 2015 - 2018 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>

#include <argyllio.h>
#include <observers.h>
#include <spectraldb.h>
#include <wcompat.h>
#include <elog.h>

// like fgets() but skips empty lines and can handle any combination of \r and \n
static char *
fgets_data(char *s, size_t sz, char *data, size_t *data_offset, size_t data_size, bool to_upper)
{
    if (*data_offset == data_size) {
        return NULL;
    }
    size_t i, j = 0;
    for (i = *data_offset; i < data_size && j < sz - 1; i++) {
        s[j] = to_upper ? toupper(data[i]) : data[i];
        if (s[j] == '\n' || s[j] == '\r') {
            if (j > 0) {
                s[j] = '\n';
                j++;
                break;
            }
        } else {
            j++;
        }
    }
    s[j] = '\0';
    *data_offset = i;
    return s;
}

static double
atof_flex(const char *nptr)
{
    char *ep = NULL;
    double val = strtod(nptr, &ep);
    if (ep == nptr) {
        return val;
    }
    if (*ep == ',') {
        int len = (int)((uintptr_t)ep - (uintptr_t)nptr);
        char s[len+32];
        strncpy(s, nptr, sizeof(s));
        s[sizeof(s)-1] = '\0';
        char *p = strchr(s, '.');
        if (p != NULL && (uintptr_t)p < (uintptr_t)ep) {
            return val;
        }
        s[len] = '.';
        val = strtod(s, &ep);
    }
    return val;
}

static bool
is_empty_str(const char str[])
{
    if (str == NULL) return true;
    const char *p = str;
    while (*p != '\0') {
        if (!isblank(*p) && *p != '\n') {
            return false;
        }
        p++;
    }
    return true;
}

struct patch_set *
argyllio_ti3_parse(const char filename[],
                   bool error_is_fatal,
                   bool verbose)
{
    double *bands = NULL;
    double *fields = NULL;
    char **strings = NULL;
    struct patch_set *ti3 = NULL;
    size_t data_offset = 0, data_size;
    char *data = NULL;
    {
        data = (char *)fileio_file2mem(filename, &data_size);
        if (data == NULL) {
            goto fail;
        }

        // line-feed pre-processing: if a line starts with space or tab then remove previous line break
        ssize_t linefeed = -1;
        for (ssize_t i = 0; i < (ssize_t)data_size; i++) {
            if (data[i] == '\n' || data[i] == '\r') {
                if (linefeed == -1) linefeed = i;
            } else {
                if (linefeed != -1 && isblank(data[i])) {
                    while (linefeed != i) data[linefeed++] = ' ';
                }
                linefeed = -1;
            }
        }
    }
    char s[8192];
    bool found_begin = false;
    int id_offset = -1;
    int class_offset = -1;
    int stype_offset = -1;
    int xyz_offset[3] = { -1, -1, -1 };
    int xyY_offset[3] = { -1, -1, -1 };
    int rgb_offset[3] = { -1, -1, -1 };
    int lab_offset[3] = { -1, -1, -1 };
    int lch_offset[3] = { -1, -1, -1 };
    int spec_offset = -1;
    int spec_band_count = -1;
    int row_count = 0;
    double patch_aspect_ratio = 0;
    double patch_inner_ratio = 0;
    double spectral_scaling = 1.0;
    bool spec_nm_specified = false;
    bool spec_pct = false;
    // normally illuminant is either D50 or D65, that's the only two we try to identify
    enum exif_lightsource illuminant = lsUnknown;
    const struct observer *obs = observer_get(OBSERVER_1931_2);
    const spectrum_t *ill_spec = spectraldb_illuminant(lsD50);
    v3 ill_wp = (v3){{ 0.9642, 1.0000, 0.8249 }}; // D50, assume that if we don't find anything else
    const char *spec_prefixes[] = {
        "SPEC_",
        "SPECTRAL_NM_",
        "NM",
        "R_",
        NULL
    };
    struct field_info {
        const char *name;
        int offset;
    };
    const char *ignore_fields[] = {
        //"SAMPLE_ID",
        //"RGB_R",
        //"RGB_G",
        //"RGB_B",
        //"SPECTRAL_NM",
        //"SPECTRAL_PCT",
        //"SPECTRAL_DEC",
        //"XYZ_X",
        //"XYZ_Y",
        //"XYZ_Z",
        //"XYY_X",
        //"XYY_Y",
        //"XYY_CAPY",
        //"LAB_L",
        //"LAB_A",
        //"LAB_B",
        //"LAB_C",
        //"LAB_H",
        "STRING",
        "CMYK_C",
        "CMYK_M",
        "CMYK_Y",
        "CMYK_K",
        "D_RED",
        "D_GREEN",
        "D_BLUE",
        "D_VIS",
        "D_MAJOR_FILTER",
        "LAB_DE",
        "LAB_DE_94",
        "LAB_DE_CMC",
        "LAB_DE_2000",
        "MEAN_DE",
        "STDEV_X",
        "STDEV_Y",
        "STDEV_Z",
        "STDEV_L",
        "STDEV_A",
        "STDEV_B",
        "STDEV_DE",
        "CHI_SQD_PAR",
        NULL
    };
    while (fgets_data(s, sizeof(s), data, &data_offset, data_size, true) != NULL) {
        if (strstr(s, "CTI3") == s || strstr(s, "IT8.7/2") == s) {
            spectral_scaling = 100.0;
        } else if (strstr(s, "SPEC_SCALE") == s) {
            spectral_scaling = atof_flex(&s[strlen("SPEC_SCALE")]);
        } else if (strstr(s, "LGOROWLENGTH") == s) { // non-standard but defacto way to specify row count
            row_count = atoi(&s[strlen("LGOROWLENGTH")]);
        } else if (strstr(s, "DCAMPROF_PATCH_ASPECT_RATIO") == s) { // DCamProf tag
            patch_aspect_ratio = atof_flex(&s[strlen("DCAMPROF_PATCH_ASPECT_RATIO")]);
        } else if (strstr(s, "DCAMPROF_PATCH_INNER_RATIO") == s) { // DCamProf tag
            patch_inner_ratio = atof_flex(&s[strlen("DCAMPROF_PATCH_INNER_RATIO")]);
        } else if (strstr(s, "DESCRIPTOR") == s ||
                   strstr(s, "MEASUREMENT_CONDITION") == s ||
                   strstr(s, "MEASUREMENT_SOURCE") == s ||
                   strstr(s, "ILLUMINATION_NAME") == s ||
                   strstr(s, "ILLUMINANT_NAME") == s ||
                   strstr(s, "ILLUMINANT NAME") == s ||
                   strstr(s, "WEIGHTING_FUNCTION") == s)
        {
            char *p = s;
            while (!isblank(*p)) p++;
            if (strstr(p, "D50") != NULL) {
                illuminant = lsD50;
                ill_spec = spectraldb_illuminant(lsD50);
                ill_wp = (v3){{ 0.9642, 1.0000, 0.8249 }}; // D50
            } else if (strstr(p, "D65") != NULL) {
                illuminant = lsD65;
                ill_spec = spectraldb_illuminant(lsD65);
                ill_wp = (v3){{ 0.950456, 1.0, 1.088754 }}; // D65
            }
        } else if (strstr(s, "SPECTRAL_NORM") == s) {
            char *p = &s[strlen("SPECTRAL_NORM")];
            if (strchr(p, '\"') != NULL) p = strchr(p, '\"') + 1;
            spectral_scaling = atof_flex(p);
        } else if (strstr(s, "BEGIN_DATA_FORMAT") == s && is_empty_str(&s[17])) {
            if (fgets_data(s, sizeof(s), data, &data_offset, data_size, true) == NULL) {
                elog("Could not read data format\n");
                goto fail;
            }
            char *p = s;
            int field = 0;
            for (;;) {
                while (isblank(*p)) p++;
                if (*p == '\n' || *p == '\0') {
                    break;
                }
                bool ignore_field = false;
                for (int i = 0; ignore_fields[i] != NULL; i++) {
                    const int field_len = strlen(ignore_fields[i]);
                    if (strncmp(p, ignore_fields[i], field_len) == 0) {
                        ignore_field = true;
                        p += field_len;
                        field++;
                    }
                }
                if (ignore_field) {
                    continue;
                }
                if (strncmp(p, "RGB_R", 5) == 0) {
                    rgb_offset[0] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "RGB_G", 5) == 0) {
                    rgb_offset[1] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "RGB_B", 5) == 0) {
                    rgb_offset[2] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYZ_X", 5) == 0) {
                    xyz_offset[0] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYZ_Y", 5) == 0) {
                    xyz_offset[1] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYZ_Z", 5) == 0) {
                    xyz_offset[2] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYY_X", 5) == 0) {
                    xyY_offset[0] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYY_Y", 5) == 0) {
                    xyY_offset[1] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "XYY_CAPY", 8) == 0) {
                    xyY_offset[2] = field;
                    p += 8;
                    field++;
                } else if (strncmp(p, "LAB_L", 5) == 0) {
                    lab_offset[0] = field;
                    lch_offset[0] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "LAB_A", 5) == 0) {
                    lab_offset[1] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "LAB_B", 5) == 0) {
                    lab_offset[2] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "LAB_C", 5) == 0) {
                    lch_offset[1] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "LAB_H", 5) == 0) {
                    lch_offset[2] = field;
                    p += 5;
                    field++;
                } else if (strncmp(p, "SAMPLE_CLASS", 12) == 0) {
                    class_offset = field;
                    p += 12;
                    field++;
                } else if (strncmp(p, "SAMPLE_ID", 9) == 0) {
                    if (id_offset == -1) {
                        id_offset = field;
                    }
                    p += 9;
                    field++;
                } else if (strncmp(p, "SAMPLE_NAME", 11) == 0) { // overrides SAMPLE_ID
                    id_offset = field;
                    p += 11;
                    field++;
                } else if (strncmp(p, "SAMPLE_LOC", 10) == 0) { // overrides SAMPLE_ID
                    id_offset = field;
                    p += 10;
                    field++;
                } else if (strncmp(p, "SAMPLE_TYPE", 11) == 0) {
                    stype_offset = field;
                    p += 11;
                    field++;
                } else if (strncmp(p, "SPECTRAL_NM", 11) == 0 && isblank(p[11])) {
                    // assume SPECTRAL_NM SPECTRAL_PCT special case
                    spec_offset = field;
                    p += 11;
                    field++;
                    int bc = 0;
                    bool has_pct = false, has_dec = false;
                    for (;;) {
                        while (isblank(*p)) p++;
                        if (strncmp(p, "SPECTRAL_PCT", 12) == 0) {
                            has_pct = true;
                        } else if (strncmp(p, "SPECTRAL_DEC", 12) == 0) {
                            has_dec = true;
                        } else {
                            elog("Expected SPECTRAL_PCT or SPECTRAL_DEC column after SPECTRAL_NM.\n");
                            goto fail;
                        }
                        bc++;
                        field++;
                        p += 12;
                        while (isblank(*p)) p++;
                        if (strncmp(p, "SPECTRAL_NM", 11) != 0) {
                            break;
                        }
                        field++;
                        p += 11;
                        while (isblank(*p)) p++;
                    }
                    if (has_pct && has_dec) {
                        elog("Cannot mix SPECTRAL_PCT and SPECTRAL_DEC in the same file.\n");
                        goto fail;
                    }
                    spec_pct = has_pct;
                    spec_band_count = 2 * bc;
                    spec_nm_specified = true;
                    bands = malloc(bc * sizeof(bands[0]));
                } else {
                    int sp;
                    for (sp = 0; spec_prefixes[sp] != NULL; sp++) {
                        const int prefix_len = strlen(spec_prefixes[sp]);
                        if (strncmp(p, spec_prefixes[sp], prefix_len) == 0) {
                            spec_offset = field;
                            bands = malloc(strlen(s) * sizeof(bands[0]));
                            int bc = 0;
                            for (;;) {
                                p += prefix_len;
                                bands[bc++] = (double)atoi(p);
                                while (!isblank(*p) && *p != '\0') p++;
                                while (isblank(*p)) p++;
                                if (strncmp(p, spec_prefixes[sp], prefix_len) != 0) {
                                    break;
                                }
                            }
                            for (int i = 0; i < bc; i++) {
                                if (bands[i] <= 0 || (i > 0 && bands[i] <= bands[i-1])) {
                                    elog("Bad spectral bands (%f).\n", bands[i]);
                                    goto fail;
                                }
                            }
                            if (bc < 3) {
                                elog("Too few spectral bands.\n");
                                goto fail;
                            }
                            spec_band_count = bc;
                            field += bc;
                            break;
                        }
                    }
                    if (spec_prefixes[sp] == NULL) {
                        // ignore unknown field
                        while (!isblank(*p) && *p != '\0') p++;
                        field++;
                    }
                }
            }
        } else if (strstr(s, "BEGIN_DATA") == s && is_empty_str(&s[10])) {
            found_begin = true;
            break;
        }
    }
    const bool has_xyz = xyz_offset[0] >= 0 && xyz_offset[1] >= 0 && xyz_offset[2] >= 0;
    const bool has_xyY = xyY_offset[0] >= 0 && xyY_offset[1] >= 0 && xyY_offset[2] >= 0;
    const bool has_lab = lab_offset[0] >= 0 && lab_offset[1] >= 0 && lab_offset[2] >= 0;
    const bool has_lch = lch_offset[0] >= 0 && lch_offset[1] >= 0 && lch_offset[2] >= 0;
    const bool has_rgb = rgb_offset[0] >= 0 && rgb_offset[1] >= 0 && rgb_offset[2] >= 0;
    if (!has_xyz && !has_xyY) {
        if (!has_lab && !has_lch) {
            if (verbose) {
                elog("Missing XYZ fields, creating dummy values.\n");
            }
        } else {
            if (verbose) {
                elog("Missing XYZ fields, creating through conversion from provided LAB values.\n");
            }
        }
    }
    if (!has_rgb) {
        if (verbose) {
            elog("Missing RGB fields, creating dummy values.\n");
        }
    }
    if (!found_begin) {
        elog("BEGIN_DATA not found.\n");
        goto fail;
    }
    int alloc_size = 1;
    ti3 = calloc(1, sizeof(*ti3) + alloc_size * sizeof(struct cam2xyz_patch));
    fields = malloc(sizeof(s) * sizeof(fields[0]));
    strings = malloc(sizeof(s) * sizeof(strings[0]));
    for (;;) {
        if (fgets_data(s, sizeof(s), data, &data_offset, data_size, false) == NULL) {
            break;
        }
        if (strncasecmp(s, "END_DATA", 8) == 0) {
            break;
        }
        char *p = s;
        while (isblank(*p)) p++;
        int fi = 0;
        const int i = ti3->patch_count;
        while (*p != '\0' && *p != '\n') {
            if (*p == '\"') {
                fields[fi] = 0;
                strings[fi++] = ++p;
                if ((p = strchr(p, '\"')) == NULL) {
                    elog("Bad string on sample data line %d.\n", i + 1);
                    goto fail;
                }
                *p = '\0'; p++;
            } else {
                if (fi == id_offset) {
                    fields[fi] = 0;
                    strings[fi++] = p;
                    while (!isblank(*p) && *p != '\0') p++;
                    if (*p != '\0') {
                        *p = '\0'; p++;
                    }
                } else {
                    strings[fi] = NULL;
                    fields[fi++] = atof_flex(p);
                    while (!isblank(*p) && *p != '\0') p++;
                }
            }
            while (isblank(*p)) p++;
        }
        if (fi <= xyz_offset[0] || fi <= xyz_offset[1] || fi <= xyz_offset[2] ||
            fi <= xyY_offset[0] || fi <= xyY_offset[1] || fi <= xyY_offset[2] ||
            fi <= lab_offset[0] || fi <= lab_offset[1] || fi <= lab_offset[2] ||
            fi <= rgb_offset[0] || fi <= rgb_offset[1] || fi <= rgb_offset[2] ||
            fi < spec_offset + spec_band_count)
        {
            elog("Too few data fields on sample data line %d.\n", i + 1);
            goto fail;
        }
        if (id_offset >= 0) {
            if (strlen(strings[id_offset]) > sizeof(ti3->patch[i].str_id)-1) {
                elog("Patch name \"%s\" is longer than the maximum supported length which is %d characters. Edit the file and replace with shorter patch names.\n", strings[id_offset], (int)(sizeof(ti3->patch[i].str_id)-1));
                goto fail;
            }
            strncpy(ti3->patch[i].str_id, strings[id_offset], sizeof(ti3->patch[i].str_id)-1);
        } else {
            sprintf(ti3->patch[i].str_id, "%d", i + 1);
        }
        if (has_xyz || has_xyY) {
            v3 xyz;
            if (has_xyz) {
                xyz.v[0] = fields[xyz_offset[0]] / 100.0;
                xyz.v[1] = fields[xyz_offset[1]] / 100.0;
                xyz.v[2] = fields[xyz_offset[2]] / 100.0;
            } else { // has_xyY
                v3 xyY;
                xyY.v[0] = fields[xyY_offset[0]];
                xyY.v[1] = fields[xyY_offset[1]];
                xyY.v[2] = fields[xyY_offset[2]] / 100.0;
                xyz = xyY2xyz(xyY);
            }
            ti3->patch[i].xyz = xyz;
        }
        if (has_rgb) {
            ti3->patch[i].cam.v[0] = fields[rgb_offset[0]] / 100.0;
            ti3->patch[i].cam.v[1] = fields[rgb_offset[1]] / 100.0;
            ti3->patch[i].cam.v[2] = fields[rgb_offset[2]] / 100.0;
        }
        if ((has_lab || has_lch) && (!has_xyz && !has_xyY)) {
            v3 lab;
            if (has_lab) {
                lab.v[0] = fields[lab_offset[0]];
                lab.v[1] = fields[lab_offset[1]];
                lab.v[2] = fields[lab_offset[2]];
            } else { // has_lch
                v3 lch;
                lch.v[0] = fields[lch_offset[0]];
                lch.v[1] = fields[lch_offset[1]];
                lch.v[2] = fields[lch_offset[2]];
                lab = lch2lab(lch);
            }
            ti3->patch[i].xyz = lab2xyz(lab, ill_wp);
        }
        if (stype_offset >= 0) {
            ti3->patch[i].emissive = strings[stype_offset] != NULL && strcmp(strings[stype_offset], "E") == 0;
        }
        if (class_offset >= 0) {
            if (strings[class_offset] == NULL) {
                elog("Bad class name on sample data line %d.\n", i + 1);
                goto fail;
            }
            ti3->patch[i].class = -1;
            int j;
            for (j = 0; j < TARGET_MAX_CLASS_COUNT && ti3->class_names[j][0] != '\0'; j++) {
                if (strcasecmp(ti3->class_names[j], strings[class_offset]) == 0) {
                    ti3->patch[i].class = j;
                    break;
                }
            }
            if (ti3->patch[i].class == -1) {
                if (j == TARGET_MAX_CLASS_COUNT) {
                    elog("Too many class names, max %d supported.\n", TARGET_MAX_CLASS_COUNT);
                    goto fail;
                }
                if (strlen(strings[class_offset]) >= TARGET_CLASS_NAME_SIZE) {
                    elog("Too long class name on sample data line %d.\n", i+1);
                    goto fail;
                }
                strcpy(ti3->class_names[j], strings[class_offset]);
                ti3->patch[i].class = j;
            }
        }
        if (spec_offset >= 0) {
            if (spectral_scaling != 1.0 || spec_pct) {
                double mul = 1.0 / (spectral_scaling * (spec_pct ? 100.0 : 1.0));
                for (int j = 0; j < spec_band_count; j++) {
                    if (spec_nm_specified) j++;
                    fields[spec_offset+j] *= mul;
                }
            }
            if (spec_nm_specified) {
                double amps[spec_band_count/2];
                for (int j = 0; j < spec_band_count; j += 2) {
                    bands[j/2] = fields[spec_offset+j];
                    amps[j/2] = fields[spec_offset+j+1];
                }
                for (int j = 0; j < spec_band_count/2; j++) {
                    if (bands[j] <= 0 || (j > 0 && bands[j] <= bands[j-1])) {
                        elog("Bad spectral bands (%f) at sample data line %d.\n", bands[j], i + 1);
                        goto fail;
                    }
                }
                ti3->patch[i].spectrum = spec_alloc1(bands, amps, spec_band_count / 2);
            } else {
                ti3->patch[i].spectrum = spec_alloc1(bands, &fields[spec_offset], spec_band_count);
            }
            if (ti3->patch[i].spectrum == NULL) {
                elog("Bad spectral data at sample data line %d.\n", i + 1);
                goto fail;
            }
            if (!has_xyz) {
                ti3->patch[i].xyz = spec2xyz_ill(ti3->patch[i].spectrum, obs, ill_spec, 1);
            }
        }
        ti3->patch_count++;
        if (ti3->patch_count == alloc_size) {
            alloc_size += 256;
            ti3 = realloc(ti3, sizeof(*ti3) + alloc_size * sizeof(struct cam2xyz_patch));
            memset(&ti3->patch[ti3->patch_count], 0, 256 * sizeof(struct cam2xyz_patch));
        }
    }
    free(data);
    free(bands);
    free(fields);
    free(strings);
    ti3->layout.row_count = row_count;
    ti3->layout.patch_aspect_ratio = patch_aspect_ratio;
    ti3->layout.patch_inner_ratio = patch_inner_ratio;
    ti3->illuminant = illuminant;
    ti3->ill_wp = ill_wp;
    ti3 = realloc(ti3, sizeof(*ti3) + ti3->patch_count * sizeof(struct cam2xyz_patch));
    return ti3;
fail:
    free(data);
    free(ti3);
    free(bands);
    free(fields);
    free(strings);
    elog("Failed to parse \"%s\".\n", filename);
    if (error_is_fatal) {
        exit(EXIT_FAILURE);
    }
    return NULL;
}

bool
argyllio_ti3_write(const char filename[],
                   bool print_spectral_data,
                   const struct patch_set *ps,
                   bool error_is_fatal)
{
    int *bands = NULL;
    int band_count = 0;
    bool matching_bands[ps->patch_count];
    if (print_spectral_data) {
        int band_start = 10000, band_end = 0;
        for (int i = 0; i < ps->patch_count; i++) {
            int bs = round(spec_band_start(ps->patch[i].spectrum));
            int be = round(spec_band_end(ps->patch[i].spectrum));
            if (bs < band_start) band_start = bs;
            if (be > band_end) band_end = be;
        }
        if (band_end < band_start) { // should not happen
            band_end = band_start;
        }
        bool *has_f = calloc(1, ((band_end - band_start) + 1) * sizeof(has_f[0]));
        for (int i = 0; i < ps->patch_count; i++) {
            const spectrum_t *s = ps->patch[i].spectrum;
            if (s == NULL) {
                elog("At least one patch lacks spectrum information, cannot write spectral data.\n");
                if (error_is_fatal) exit(EXIT_FAILURE);
                return false;
            }
            for (int i = 0; i < spec_band_count(s); i++) {
                has_f[(int)round(s->v[i].f) - band_start] = true;
            }
        }
        bands = malloc((band_end - band_start + 1) * sizeof(bands[0]));
        band_count = 0;
        for (int i = 0; i < band_end - band_start + 1; i++) {
            if (has_f[i]) {
                bands[band_count++] = band_start + i;
            }
        }
        for (int i = 0; i < ps->patch_count; i++) {
            const spectrum_t *s = ps->patch[i].spectrum;
            if (spec_band_count(s) != band_count) {
                matching_bands[i] = false;
            } else {
                matching_bands[i] = true;
                for (int i = 0; i < band_count; i++) {
                    if (round(s->v[i].f) != bands[i]) {
                        matching_bands[i] = false;
                        break;
                    }
                }
            }
        }
    }
    FILE *stream = utf8_fopen(filename, "w");
    if (stream == NULL) {
        elog("Could not open \"%s\" for writing: %s.\n", filename, strerror(errno));
        if (error_is_fatal) exit(EXIT_FAILURE);
        free(bands);
        return false;
    }

    struct tm tm;
    time_t tt = time(NULL);
    char timestr[128];
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    tm = *localtime(&tt);
    strcpy(timestr, asctime(&tm));
#else
    localtime_r(&tt, &tm);
    asctime_r(&tm, timestr);
#endif
    timestr[strlen(timestr)-1] = '\0';
    fprintf(stream,
            "CTI3   \n"
            "\n"
            "DESCRIPTOR \"Argyll Calibration Target chart information 3\"\n"
            "ORIGINATOR \"DCamProf\"\n"
            "CREATED \"%s\"\n"
            "KEYWORD \"DEVICE_CLASS\"\n"
            "DEVICE_CLASS \"INPUT\"\n"
            "KEYWORD \"COLOR_REP\"\n"
            "COLOR_REP \"RGB_XYZ\"\n",
            timestr);
    if (print_spectral_data) {
        fprintf(stream,
                "KEYWORD \"SPECTRAL_BANDS\"\n"
                "SPECTRAL_BANDS \"%d\"\n"
                "KEYWORD \"SPECTRAL_START_NM\"\n"
                "SPECTRAL_START_NM \"%f\"\n"
                "KEYWORD \"SPECTRAL_END_NM\"\n"
                "SPECTRAL_END_NM \"%f\"\n",
                band_count,
                (double)bands[0],
                (double)bands[band_count-1]);
    }
    fprintf(stream,
            "\n"
            "KEYWORD \"SAMPLE_CLASS\"\n"
        );
    if (print_spectral_data) {
        for (int i = 0; i < band_count; i++) {
            fprintf(stream, "KEYWORD \"SPEC_%d\"\n", bands[i]);
        }
    }
    fprintf(stream,
            "NUMBER_OF_FIELDS %d\n"
            "BEGIN_DATA_FORMAT\n"
            "SAMPLE_ID SAMPLE_CLASS SAMPLE_TYPE RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z",
            8 + (print_spectral_data ? band_count : 0));
    if (print_spectral_data) {
        for (int i = 0; i < band_count; i++) {
            fprintf(stream, " SPEC_%d", bands[i]);
        }
    }
    fprintf(stream, "\nEND_DATA_FORMAT\n\n");

    fprintf(stream, "NUMBER_OF_SETS %d\n"
            "BEGIN_DATA\n", ps->patch_count);
    for (int i = 0; i < ps->patch_count; i++) {
        const struct cam2xyz_patch *p = &ps->patch[i];
        const char *id;
        char id_[16];
        if (p->str_id[0] != '\0') {
            id = p->str_id;
        } else {
            id = id_;
            sprintf(id_, "%d", i+1);
        }
        fprintf(stream, "%s \"%s\" \"%c\" %.6G %.6G %.6G %.6G %.6G %.6G",
                id,
                ps->class_names[p->class],
                p->emissive ? 'E' : 'R',
                100*p->cam.v[0], 100*p->cam.v[1], 100*p->cam.v[2],
                100*p->xyz.v[0], 100*p->xyz.v[1], 100*p->xyz.v[2]);
        if (print_spectral_data) {
            for (int j = 0; j < band_count; j++) {
                if (p->spectrum != NULL) {
                    double e;
                    if (matching_bands[i]) {
                        e = 100.0 * p->spectrum->v[j].e;
                    } else {
                        e = 100.0 * spec_at(p->spectrum, bands[j]);
                    }
                    e *= 1e6;
                    e = round(e);
                    e /= 1e6;
                    char buf[32];
                    sprintf(buf, " %.6G", e);
                    if (strchr(buf, 'E') != NULL) {
                        fprintf(stream, " %.6f", e);
                    } else {
                        fprintf(stream, " %.6G", e);
                    }
                } else {
                    fprintf(stream, " 0");
                }
            }
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "END_DATA\n");
    fclose(stream);
    free(bands);

    return true;
}

bool
argyllio_ti1_write(const char filename[],
                   struct patch_set *ps,
                   bool fill_in_xyz,
                   bool error_is_fatal)
{
    const m3x3 beta_rgb2xyz = {{
            { 0.6712537,  0.1745834,  0.1183829 },
            { 0.3032726,  0.6637861,  0.0329413 },
            { 0.0000000,  0.0407010,  0.7845090 }
        }};
    const v3 wp = m3x3_mul_v3(beta_rgb2xyz, (v3){{1,1,1}});

    FILE *stream = utf8_fopen(filename, "w");
    if (stream == NULL) {
        elog("Could not open \"%s\" for writing: %s.\n", filename, strerror(errno));
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    struct tm tm;
    time_t tt = time(NULL);
    char timestr[128];
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    tm = *localtime(&tt);
    strcpy(timestr, asctime(&tm));
#else
    localtime_r(&tt, &tm);
    asctime_r(&tm, timestr);
#endif
    timestr[strlen(timestr)-1] = '\0';
    fprintf(stream,
            "CTI1   \n"
            "\n"
            "DESCRIPTOR \"Argyll Calibration Target chart information 1\"\n"
            "ORIGINATOR \"DCamProf\"\n"
            "CREATED \"%s\"\n"
            "KEYWORD \"APPROX_WHITE_POINT\"\n"
            "APPROX_WHITE_POINT \"%g %g %g\"\n"
            "KEYWORD \"COLOR_REP\"\n"
            "COLOR_REP \"iRGB\"\n",
            timestr,
            wp.v[0] * 100, wp.v[1] * 100, wp.v[2] * 100);

    fprintf(stream, "\n");

    fprintf(stream,
            "NUMBER_OF_FIELDS 7\n"
            "BEGIN_DATA_FORMAT\n"
            "SAMPLE_ID RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z\n"
            "END_DATA_FORMAT\n\n");

    fprintf(stream, "NUMBER_OF_SETS %d\n"
            "BEGIN_DATA\n", ps->patch_count);
    for (int i = 0; i < ps->patch_count; i++) {
        struct cam2xyz_patch *p = &ps->patch[i];
        v3 xyz = m3x3_mul_v3(beta_rgb2xyz, p->cam);;
        fprintf(stream, "%d %.6G %.6G %.6G %.6G %.6G %.6G\n",
                i+1,
                100*p->cam.v[0], 100*p->cam.v[1], 100*p->cam.v[2],
                100*xyz.v[0], 100*xyz.v[1], 100*xyz.v[2]);
        if (fill_in_xyz) {
            p->xyz = xyz;
        }
    }
    fprintf(stream, "END_DATA\n");

    fprintf(stream,
            "CTI1   \n"
            "\n"
            "DESCRIPTOR \"Argyll Calibration Target chart information 1\"\n"
            "ORIGINATOR \"DCamProf\"\n"
            "KEYWORD \"DENSITY_EXTREME_VALUES\"\n"
            "DENSITY_EXTREME_VALUES \"8\"\n"
            "CREATED \"%s\"\n"
            "\n"
            "KEYWORD \"INDEX\"\n"
            "NUMBER_OF_FIELDS 7\n"
            "BEGIN_DATA_FORMAT\n"
            "INDEX RGB_R RGB_G RGB_B XYZ_X XYZ_Y XYZ_Z \n"
            "END_DATA_FORMAT\n"
            "\n"
            "NUMBER_OF_SETS 8\n"
            "BEGIN_DATA\n",
            timestr);

    v3 val[8] = {
        {{ 1.0, 1.0, 1.0}},
        {{ 1.0, 0.47, 1.0}},
        {{ 1.0, 0.0, 0.79}},
        {{ 0.0, 0.0, 0.58}},
        {{ 1.0, 0.66, 0.0}},
        {{ 0.0, 0.35, 0.0}},
        {{ 0.84, 0.0, 0.0}},
        {{ 0.0, 0.0, 0.0}},
    };
    for (int i = 0; i < 8; i++) {
        v3 xyz = m3x3_mul_v3(beta_rgb2xyz, val[i]);;
        fprintf(stream, "%d %.6G %.6G %.6G %.6G %.6G %.6G\n",
                i+1,
                100*val[i].v[0], 100*val[i].v[1], 100*val[i].v[2],
                100*xyz.v[0], 100*xyz.v[1], 100*xyz.v[2]);
    }

    fprintf(stream, "END_DATA\n");
    fclose(stream);
    return true;
}

static int
parse_patch_name_list(bool selected[], // patch_count length
                      struct patch_set *ps,
                      const char filename[])
{
    size_t data_offset = 0, data_size;
    char *data = (char *)fileio_file2mem(filename, &data_size);
    if (data == NULL) {
        exit(EXIT_FAILURE);
    }

    char s[8192];
    int selected_count = 0;
    memset(selected, 0, sizeof(selected[0]) * ps->patch_count);
    while (fgets_data(s, sizeof(s), data, &data_offset, data_size, false) != NULL) {
        char *p = s;
        while (isblank(*p)) p++;
        if (*p == '\n' || *p == '\r' || *p == '\0') {
            continue;
        }
        char *class = NULL;
        char *name = p;
        if (*name == '#') {
            continue;
        }
        while (!isblank(*p) && *p != '\0' && *p != '\n' && *p != '\r') p++;
        if (isblank(*p)) {
            // class + name
            *p = '\0';
            p++;
            while (isblank(*p)) p++;
            if (*p != '\n' && *p != '\0' && *p != '\r') {
                class = name;
                name = p;
                while (!isblank(*p) && *p != '\0' && *p != '\n' && *p != '\r') p++;
                *p = '\0';
            }
        } else {
            // only name
            *p = '\0';
        }
        int class_idx = -1;
        if (class != NULL) {
            class_idx = TARGET_MAX_CLASS_COUNT;
            for (int j = 0; j < TARGET_MAX_CLASS_COUNT && ps->class_names[j][0] != '\0'; j++) {
                if (strcasecmp(class, ps->class_names[j]) == 0) {
                    class_idx = j;
                    break;
                }
            }
        }
        if (isdigit(name[0])) {
            int idx = atoi(name);
            if (idx > 0 && idx <= ps->patch_count && (class_idx == -1 || ps->patch[idx-1].class == class_idx) && !selected[idx-1]) {
                selected[idx-1] = true;
                selected_count++;
            }
        } else {
            for (int i = 0; i < ps->patch_count; i++) {
                if (strcmp(ps->patch[i].str_id, name) == 0 && (class_idx == -1 || ps->patch[i].class == class_idx) && !selected[i]) {
                    selected[i] = true;
                    selected_count++;
                }
            }
        }
    }
    free(data);
    return selected_count;
}

int
textio_process_exclude_list(struct patch_set *ps,
                            const char filename[],
                            int patch_remapped_indexes[])
{
    bool exclude[ps->patch_count];
    int exclude_count = parse_patch_name_list(exclude, ps, filename);
    if (patch_remapped_indexes != NULL) {
        for (int i = 0; i < ps->patch_count; i++) {
            patch_remapped_indexes[i] = i;
        }
    }
    if (exclude_count > 0) {
        int orig_count = ps->patch_count;
        ps->patch_count = 0;
        for (int i = 0; i < orig_count; i++) {
            if (exclude[i] && !ps->patch[i].force_keep) {
                free(ps->patch[i].spectrum);
                if (patch_remapped_indexes != NULL) {
                    patch_remapped_indexes[i] = -1;
                }
                continue;
            }
            if (patch_remapped_indexes != NULL) {
                patch_remapped_indexes[i] = ps->patch_count;
            }
            if (i != ps->patch_count) {
                ps->patch[ps->patch_count] = ps->patch[i];
            }
            ps->patch_count++;
        }
    }
    return exclude_count;
}

int
textio_process_keep_list(struct patch_set *ps,
                         const char filename[])
{
    bool keep[ps->patch_count];
    int keep_count = parse_patch_name_list(keep, ps, filename);
    if (keep_count > 0) {
        for (int i = 0; i < ps->patch_count; i++) {
            if (keep[i]) {
                ps->patch[i].force_keep = true;
            }
        }
    }
    return keep_count;
}
