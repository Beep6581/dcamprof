/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef PROFIO_H_
#define PROFIO_H_

#include <stdbool.h>
#include <stdio.h>

#include <colmath.h>

enum dng_embed_policy {

        // Freely embedable and copyable into installations that encounter this
        // profile, so long as the profile is only used to process DNG files.

        pepAllowCopying                         = 0,

        // Can be embeded in a DNG for portable processing, but cannot be used
        // to process other files that the profile is not embedded in.

        pepEmbedIfUsed                          = 1,

        // Can only be used if installed on the machine processing the file.
        // Note that this only applies to stand-alone profiles.  Profiles that
        // are already embedded inside a DNG file allowed to remain embedded
        // in that DNG, even if the DNG is resaved.

        pepEmbedNever                           = 2,

        // No restricts on profile use or embedding.

        pepNoRestrictions                       = 3

};

struct profio_dcp {
    char *camera_model;
    char *profile_name;
    char *copyright;
    char *calibration_signature;
    enum dng_embed_policy embed_policy;
    enum exif_lightsource ill[2];
    m3x3 cm[2];
    m3x3 fm[2];
    int ill_count;
    uint32_t hsmdims[3];
    uint32_t lookdims[3];
    bool huesatmap_srgb_gamma;
    bool looktable_srgb_gamma;
    bool black_render_none;
    v3 *huesatmap[2];
    v3 *looktable;
    double baseline_exposure_offset;
    double *tonecurve;
    int tonecurve_len;
    // bools for optional tags
    struct {
        bool fm;
        bool huesatmap_srgb_gamma;
        bool looktable_srgb_gamma;
        bool black_render_none;
        bool baseline_exposure_offset;
        bool embed_policy;
    } has;
};

struct icc_lut {
    m3x3 cm;
    double *in[3];
    int in_len;
    double *out[3];
    int out_len;
    double *clut;
    int clut_side;
    double data[];
};

struct profio_icc {
    char *description;
    char *copyright;
    bool pcs_is_lab;
    v3 wp;
    v3 bp;
    m3x3 fm;
    double *trc[3];
    int trc_len;
    struct icc_lut *lut;
};

const char *
exif_lightsource_ntoa(enum exif_lightsource ls);

enum exif_lightsource
exif_lightsource_aton(const char str[]);

double
exif_lightsource_temp(enum exif_lightsource ls);

char *
exif_lightsource_list(const char prefix[],
                      int lw);

static inline v3
icc_d50(void)
{
    return (v3){{ 0.9642, 1.0000, 0.8249 }};
}

struct profio_dcp *
profio_dcp_parse(const char filename[],
                 bool error_is_fatal);

void
profio_dcp_json(struct profio_dcp *dcp,
                FILE *stream);

char *
profio_dcp_tojson(struct profio_dcp *dcp);

struct profio_dcp *
profio_dcp_fromraw(const uint8_t data[],
                   size_t size);

void *
profio_dcp_toraw(const struct profio_dcp *dcp,
                 size_t *size);

struct profio_dcp *
profio_dcp_copy(const struct profio_dcp *dcp);

bool
profio_dcp_write(const struct profio_dcp *dcp,
                 const char filename[],
                 bool error_is_fatal);

void
profio_dcp_delete(struct profio_dcp *dcp);

struct profio_icc *
profio_icc_parse(const char filename[],
                 bool error_is_fatal);

struct profio_icc *
profio_icc_fromraw(const uint8_t data[],
                   size_t size);

void *
profio_icc_toraw(const struct profio_icc *icc,
                 size_t *size);

bool
profio_icc_write(const struct profio_icc *icc,
                 const char filename[],
                 bool error_is_fatal);

void
profio_icc_json(struct profio_icc *icc,
                FILE *stream);

char *
profio_icc_tojson(struct profio_icc *icc);

void
profio_icc_delete(struct profio_icc *icc);

#endif
