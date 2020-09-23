/*
 * (c) Copyright 2015 - 2018 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <pthread.h>
#include <inttypes.h>
#include <time.h>

#include <wcompat.h>
#include <profio.h>
#include <jsonio.h>
#include <strbuf.h>
#include <elog.h>

#define debug(...)
//#define debug(...) elog(__VA_ARGS__)

#define GET4(dst) do {                             \
        if (data_offset + 4 > data_len) goto fail; \
        memcpy(dst, &data[data_offset], 4);        \
        if (do_swap) *dst = bit_swap32(*(dst));    \
        data_offset += 4;                          \
    } while (0)
#define GET2(dst) do {                             \
        if (data_offset + 2 > data_len) goto fail; \
        memcpy(dst, &data[data_offset], 2);        \
        if (do_swap) *dst = bit_swap16(*(dst));    \
        data_offset += 2;                          \
    } while (0)

typedef	enum {
	TIFF_NOTYPE	= 0,	/* placeholder */
	TIFF_BYTE	= 1,	/* 8-bit unsigned integer */
	TIFF_ASCII	= 2,	/* 8-bit bytes w/ last byte null */
	TIFF_SHORT	= 3,	/* 16-bit unsigned integer */
	TIFF_LONG	= 4,	/* 32-bit unsigned integer */
	TIFF_RATIONAL	= 5,	/* 64-bit unsigned fraction */
	TIFF_SBYTE	= 6,	/* !8-bit signed integer */
	TIFF_UNDEFINED	= 7,	/* !8-bit untyped data */
	TIFF_SSHORT	= 8,	/* !16-bit signed integer */
	TIFF_SLONG	= 9,	/* !32-bit signed integer */
	TIFF_SRATIONAL	= 10,	/* !64-bit signed fraction */
	TIFF_FLOAT	= 11,	/* !32-bit IEEE floating point */
	TIFF_DOUBLE	= 12,	/* !64-bit IEEE floating point */
	TIFF_IFD	= 13	/* %32-bit unsigned integer (offset) */
} TIFFDataType;

#define TIFFTAG_UNIQUECAMERAMODEL 50708
#define TIFFTAG_COLORMATRIX1      50721
#define TIFFTAG_COLORMATRIX2      50722
#define TIFFTAG_REDUCTIONMATRIX1  50725 //
#define TIFFTAG_REDUCTIONMATRIX2  50726 //
#define TIFFTAG_CALIBRATIONILLUMINANT1 50778
#define TIFFTAG_CALIBRATIONILLUMINANT2 50779
#define TIFFTAG_PROFILECALIBRATIONSIGNATURE 50932
#define TIFFTAG_PROFILENAME 50936
#define TIFFTAG_PROFILEHUESATMAPDIMS 50937
#define TIFFTAG_PROFILEHUESATMAPDATA1 50938
#define TIFFTAG_PROFILEHUESATMAPDATA2 50939
#define TIFFTAG_PROFILETONECURVE 50940
#define TIFFTAG_PROFILEEMBEDPOLICY 50941
#define TIFFTAG_PROFILECOPYRIGHT 50942
#define TIFFTAG_FORWARDMATRIX1 50964
#define TIFFTAG_FORWARDMATRIX2 50965
#define TIFFTAG_PROFILELOOKTABLEDIMS 50981
#define TIFFTAG_PROFILELOOKTABLEDATA 50982
#define TIFFTAG_PROFILEHUESATMAPENCODING 51107 // long 1
#define TIFFTAG_PROFILELOOKTABLEENCODING 51108
#define TIFFTAG_BASELINEEXPOSUREOFFSET 51109 // rational 1
#define TIFFTAG_DEFAULTBLACKRENDER 51110 // long 1


struct tiff_tag {
    uint16_t id;
    uint16_t type;
    uint32_t count;
    uint32_t data;
    bool parsed_subdir;
    uint32_t ext_data_size;
    uint8_t *ext_data;
};

static int
tag_value_size(uint16_t type)
{
    int vsize = 1;
    switch (type) {
    case TIFF_UNDEFINED:
    case TIFF_SBYTE:
    case TIFF_BYTE:
    case TIFF_ASCII:
        vsize = 1;
        break;
    case TIFF_SHORT:
    case TIFF_SSHORT:
        vsize = 2;
        break;
    case TIFF_LONG:
    case TIFF_SLONG:
    case TIFF_FLOAT:
        vsize = 4;
        break;
    case TIFF_RATIONAL:
    case TIFF_SRATIONAL:
    case TIFF_DOUBLE:
        vsize = 8;
        break;
    case TIFF_IFD:
        vsize = 4;
        break;
    }
    return vsize;
}

static bool
memiszero(const void *mem,
          size_t size)
{
    const uint8_t *p = (const uint8_t *)mem;
    for (size_t i = 0; i < size; i++) {
        if (p[i] != 0) return false;
    }
    return true;
}

struct profio_dcp *
profio_dcp_parse(const char filename[],
                 bool error_is_fatal)
{
    size_t data_len;
    uint8_t *data = fileio_file2mem(filename, &data_len);
    if (data == NULL) {
        elog("Could not read DCP in file \"%s\".\n", filename);
        if (error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    struct profio_dcp *dcp = profio_dcp_fromraw(data, data_len);
    free(data);
    if (dcp == NULL) {
        elog("Could not parse DCP in file \"%s\".\n", filename);
        if (error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    return dcp;
}

struct profio_dcp *
profio_dcp_copy(const struct profio_dcp *dcp)
{
    if (dcp == NULL) {
        return NULL;
    }
    struct profio_dcp *dcp1 = malloc(sizeof(*dcp1));
    memcpy(dcp1, dcp, sizeof(*dcp1));
    if (dcp->camera_model != NULL) dcp1->camera_model = strdup(dcp->camera_model);
    if (dcp->profile_name != NULL) dcp1->profile_name = strdup(dcp->profile_name);
    if (dcp->copyright != NULL) dcp1->copyright = strdup(dcp->copyright);
    if (dcp->calibration_signature != NULL) dcp1->calibration_signature = strdup(dcp->calibration_signature);
    if (dcp->tonecurve != NULL) {
        dcp1->tonecurve = malloc(2 * dcp->tonecurve_len * sizeof(dcp1->tonecurve[0]));
        memcpy(dcp1->tonecurve, dcp->tonecurve, 2 * dcp->tonecurve_len * sizeof(dcp1->tonecurve[0]));
    }
    if (dcp->looktable != NULL) {
        size_t table_size =  dcp->lookdims[0]*dcp->lookdims[1]*dcp->lookdims[2] * sizeof(dcp->looktable[0]);
        dcp1->looktable = malloc(table_size);
        memcpy(dcp1->looktable, dcp->looktable, table_size);
    }
    for (int i = 0; i < 2; i++) {
        if (dcp->huesatmap[i] != NULL) {
            size_t table_size =  dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2] * sizeof(dcp->huesatmap[i][0]);
            dcp1->huesatmap[i] = malloc(table_size);
            memcpy(dcp1->huesatmap[i], dcp->huesatmap[i], table_size);
        }
    }
    return dcp1;
}

struct profio_dcp *
profio_dcp_fromraw(const uint8_t data[],
                   size_t data_len)
{
    uint8_t *tag_ext_data = NULL;
    struct profio_dcp *dcp = NULL;
    size_t data_offset = 0;
    if (data == NULL || data_len < 8) {
        debug("too little data\n");
        goto fail;
    }
    // TIFF directory format
    uint8_t bo[2];
    bool little_endian;
    memcpy(bo, data, 2);
    if (bo[0] == 'I' && bo[1] == 'I') {
        little_endian = true;
    } else if (bo[0] == 'M' && bo[1] == 'M') {
        little_endian = false;
    } else {
        debug("unknown tiff header %d %d\n", bo[0], bo[1]);
        goto fail;
    }

    bool do_swap;
    { // test if we need to swap data to get host order
        uint16_t t = 1;
        memcpy(bo, &t, 2);
        if ((bo[0] == 1 && !little_endian) || (bo[1] == 1 && little_endian)) {
            do_swap = true;
        } else {
            do_swap = false;
        }
    }

    uint32_t ifdoffset;
    { // goto first IFDOffset
        data_offset = 4;
        GET4(&ifdoffset);
        debug("IFD offset %d\n", ifdoffset);
        if (ifdoffset > data_len) {
            debug("IFD Offset seek to %d failed\n", ifdoffset);
            goto fail;
        }
        data_offset = ifdoffset;
    }

    dcp = calloc(1, sizeof(*dcp));
    dcp->ill_count = 1;
    uint16_t tag_count;
    uint32_t hsd_count[3] = { 0,0,0 };
    GET2(&tag_count);
    for (int i = 0; i < (int)tag_count; i++) {
        uint16_t tag_id;
        uint16_t tag_type;
        uint32_t tag_count;
        uint32_t tag_data;
        uint32_t tag_ext_data_size;
        GET2(&tag_id);
        GET2(&tag_type);
        GET4(&tag_count);
        GET4(&tag_data);
        const int vsize = tag_value_size(tag_type);
        debug("tag %05d|0x%04X\t%d\t%d\t0x%08X\t(%d)\t%d [%d - %d]\n", tag_id, tag_id, tag_type, tag_count, tag_data, vsize * tag_count, (int)data_offset, tag_data, tag_data + vsize * tag_count);
        if (vsize * tag_count > 4) {
            tag_ext_data_size = vsize * tag_count + (tag_type == TIFF_ASCII ? 1 : 0);
            tag_ext_data = malloc(tag_ext_data_size);
            if (tag_data + vsize * tag_count > data_len) {
                debug("Failed to read %d bytes (%d / %d) at %d, returned\n", vsize * tag_count, vsize, tag_count, tag_data);

                goto fail;
            }
            memcpy(tag_ext_data, &data[tag_data], vsize * tag_count);
            if (tag_type == TIFF_ASCII) {
                // ascii strings are generally zero-terminated, but here we guarantee it
                tag_ext_data[tag_count] = '\0';
            }
        } else {
            uint32_t raw_data = do_swap ? bit_swap32(tag_data) : tag_data;
            tag_ext_data_size = vsize * tag_count + (tag_type == TIFF_ASCII ? 1 : 0);
            if (tag_ext_data_size < 4) {
                tag_ext_data_size = 4;
            }
            tag_ext_data = malloc(tag_ext_data_size);
            memcpy(tag_ext_data, &raw_data, 4);
        }

        switch (tag_id) {
        case TIFFTAG_UNIQUECAMERAMODEL:
            dcp->camera_model = strdup((char *)tag_ext_data);
            break;
        case TIFFTAG_REDUCTIONMATRIX1:
        case TIFFTAG_REDUCTIONMATRIX2:
            elog("ReductionMatrix support not implemented!\n");
            goto fail;
            break;
        case TIFFTAG_COLORMATRIX1:
        case TIFFTAG_COLORMATRIX2:
        case TIFFTAG_FORWARDMATRIX1:
        case TIFFTAG_FORWARDMATRIX2: {
            if (tag_type != TIFF_SRATIONAL || tag_count != 9) {
                debug("unexpected matrix content");
                goto fail;
            }
            int32_t raw[9*2];
            m3x3 cm;
            memcpy(raw, tag_ext_data, sizeof(raw));
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    if (do_swap) {
                        raw[2*(row*3+col)+0] = bit_swap32(raw[2*(row*3+col)+0]);
                        raw[2*(row*3+col)+1] = bit_swap32(raw[2*(row*3+col)+1]);
                    }
                    cm.v[row][col] = (double)raw[2*(row*3+col)+0] / raw[2*(row*3+col)+1];
                }
            }
            switch (tag_id) {
            case TIFFTAG_COLORMATRIX1: memcpy(&dcp->cm[0], &cm, sizeof(cm)); break;
            case TIFFTAG_COLORMATRIX2: memcpy(&dcp->cm[1], &cm, sizeof(cm)); break;
            case TIFFTAG_FORWARDMATRIX1: memcpy(&dcp->fm[0], &cm, sizeof(cm)); dcp->has.fm = true; break;
            case TIFFTAG_FORWARDMATRIX2: memcpy(&dcp->fm[1], &cm, sizeof(cm)); break;
            }
            break;
        }
        case TIFFTAG_CALIBRATIONILLUMINANT1:
        case TIFFTAG_CALIBRATIONILLUMINANT2: {
            if (tag_type != TIFF_SHORT) {
                debug("unexpected illuminant content");
                goto fail;
            }
            int ill_idx = 0;
            if (tag_id == TIFFTAG_CALIBRATIONILLUMINANT2) {
                ill_idx = 1;
                dcp->ill_count = 2;
            }
            uint16_t ill;
            memcpy(&ill, tag_ext_data, 2);
            dcp->ill[ill_idx] = do_swap ? bit_swap16(ill) : ill;
            break;
        }
        case TIFFTAG_PROFILECALIBRATIONSIGNATURE:
            dcp->calibration_signature = strdup((char *)tag_ext_data);
            break;
        case TIFFTAG_PROFILENAME:
            dcp->profile_name = strdup((char *)tag_ext_data);
            break;
        case TIFFTAG_PROFILEHUESATMAPDIMS:
        case TIFFTAG_PROFILELOOKTABLEDIMS:
            if (tag_type != TIFF_LONG || tag_count != 3) {
                debug("unexpected map dims content");
                goto fail;
            }
            uint32_t dims[3];
            for (int i = 0; i < 3; i++) {
                dims[i] = do_swap ? bit_swap32(((uint32_t *)tag_ext_data)[i]) : ((uint32_t *)tag_ext_data)[i];
            }
            if (tag_id == TIFFTAG_PROFILEHUESATMAPDIMS) {
                memcpy(dcp->hsmdims, dims, sizeof(dims));
            } else {
                memcpy(dcp->lookdims, dims, sizeof(dims));
            }
            break;
        case TIFFTAG_PROFILEHUESATMAPDATA1:
        case TIFFTAG_PROFILEHUESATMAPDATA2:
        case TIFFTAG_PROFILELOOKTABLEDATA:
            if (tag_type != TIFF_FLOAT || tag_count == 0 || tag_count % 3 != 0) {
                debug("unexpected map content");
                goto fail;
            }
            v3 *map = malloc((tag_count / 3) * sizeof(*map));
            for (int i = 0; i < (int)tag_count / 3; i++) {
                float row[3];
                if (do_swap) {
                    uint32_t raw[3];
                    memcpy(raw, &((uint32_t *)tag_ext_data)[3*i], 12);
                    for (int j = 0; j < 3; j++) {
                        raw[j] = bit_swap32(raw[j]);
                    }
                    memcpy(row, raw, 12);
                } else {
                    memcpy(row, &((uint32_t *)tag_ext_data)[3*i], 12);
                }
                map[i].v[0] = row[0];
                map[i].v[1] = row[1];
                map[i].v[2] = row[2];
            }
            switch (tag_id) {
            case TIFFTAG_PROFILEHUESATMAPDATA1: dcp->huesatmap[0] = map; hsd_count[0] = tag_count; break;
            case TIFFTAG_PROFILEHUESATMAPDATA2: dcp->huesatmap[1] = map; hsd_count[1] = tag_count;  break;
            case TIFFTAG_PROFILELOOKTABLEDATA: dcp->looktable = map; hsd_count[2] = tag_count; break;
            }
            break;
        case TIFFTAG_PROFILEHUESATMAPENCODING:
        case TIFFTAG_PROFILELOOKTABLEENCODING:
            if (tag_type != TIFF_LONG || tag_count != 1) {
                debug("unexpected table encoding content\n");
                goto fail;
            }
            if (tag_id == TIFFTAG_PROFILEHUESATMAPENCODING) {
                dcp->huesatmap_srgb_gamma = !!(*(uint32_t *)tag_ext_data);
                dcp->has.huesatmap_srgb_gamma = true;
            } else {
                dcp->looktable_srgb_gamma = !!(*(uint32_t *)tag_ext_data);
                dcp->has.looktable_srgb_gamma = true;
            }
            break;
        case TIFFTAG_PROFILETONECURVE:
            if (tag_type != TIFF_FLOAT || tag_count < 4 || tag_count % 2 != 0) {
                debug("unexpected table encoding content\n");
                goto fail;
            }
            dcp->tonecurve = malloc(tag_count * sizeof(double));
            for (int i = 0; i < (int)tag_count; i++) {
                float val;
                if (do_swap) {
                    uint32_t raw = bit_swap32(((uint32_t *)tag_ext_data)[i]);
                    memcpy(&val, &raw, 4);
                } else {
                    val = ((float *)tag_ext_data)[i];
                }
                dcp->tonecurve[i] = val;
            }
            dcp->tonecurve_len = tag_count / 2;
            break;
        case TIFFTAG_PROFILEEMBEDPOLICY:
            if (tag_type != TIFF_LONG || tag_count != 1) {
                debug("unexpected embed policy content\n");
                goto fail;
            }
            dcp->embed_policy = do_swap ? bit_swap32(*(uint32_t *)tag_ext_data) : *(uint32_t *)tag_ext_data;
            dcp->has.embed_policy = true;
            if ((unsigned)dcp->embed_policy > 3) {
                debug("unexpected embed policy content %d\n", (int)dcp->embed_policy);
                goto fail;
            }
            break;
        case TIFFTAG_PROFILECOPYRIGHT:
            dcp->copyright = strdup((char *)tag_ext_data);
            break;
        case TIFFTAG_BASELINEEXPOSUREOFFSET:
            if (tag_type != TIFF_SRATIONAL || tag_count != 1) {
                debug("unexpected baseline exposure offset content\n");
                goto fail;
            }
            int32_t raw[2];
            memcpy(raw, tag_ext_data, sizeof(raw));
            if (do_swap) {
                raw[0] = bit_swap32(raw[0]);
                raw[1] = bit_swap32(raw[1]);
            }
            dcp->baseline_exposure_offset = (double)raw[0] / raw[1];
            dcp->has.baseline_exposure_offset = true;
            break;
        case TIFFTAG_DEFAULTBLACKRENDER:
            if (tag_type != TIFF_LONG || tag_count != 1) {
                debug("unexpected default black render content\n");
                goto fail;
            }
            dcp->black_render_none = !!(*(uint32_t *)tag_ext_data);
            dcp->has.black_render_none = true;
            break;
        default:
            debug("skipping unknown tag id %05d\n", tag_id);
            break;
        }
        free(tag_ext_data);
        tag_ext_data = NULL;
    }
    GET4(&ifdoffset);
    if (ifdoffset != 0) {
        debug("expected ifd termination, got %u\n", ifdoffset);
        goto fail;
    }
    if (hsd_count[0]/3 != dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2] ||
        (dcp->huesatmap[1] != NULL && hsd_count[1]/3 != dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]) ||
        hsd_count[2]/3 != dcp->lookdims[0]*dcp->lookdims[1]*dcp->lookdims[2])
    {
        debug("mismatching table dimensions\n");
        goto fail;
    }
    debug("success!\n");
    return dcp;

fail:
    free(tag_ext_data);
    profio_dcp_delete(dcp);
    return NULL;
}

static uint8_t *
tiff_tag_toraw(uint8_t *p,
               uint16_t tag_id,
               uint16_t tag_type,
               uint32_t tag_count,
               const void *tag_data,
               uint32_t *data_offset)
{
    memcpy(p, &tag_id, 2); p+= 2;
    memcpy(p, &tag_type, 2); p+= 2;
    memcpy(p, &tag_count, 4); p+= 4;
    const int vsize = tag_value_size(tag_type);
    if (vsize * tag_count > 4) {
        memcpy(p, data_offset, 4); p+= 4;
        *data_offset += (vsize * tag_count + 3) & ~0x3;
    } else {
        uint32_t data = 0;;
        if (tag_type == TIFF_BYTE || tag_type == TIFF_ASCII) {
            for (uint32_t i = 0; i < tag_count; i++) {
                data |= ((const uint8_t *)tag_data)[i] << (i*8);
            }
        } else if (tag_type == TIFF_SHORT) {
            for (uint32_t i = 0; i < tag_count; i++) {
                data |= ((const uint16_t *)tag_data)[i] << (i*16);
            }
        } else {
            memcpy(&data, tag_data, 4);
        }
        memcpy(p, &data, 4); p+= 4;
    }
    return p;
}

static uint8_t *
tiff_data_toraw(uint8_t *p,
                const void *data,
                size_t size)
{
    if (size <= 4) {
        return p;
    }
    memcpy(p, data, size);
    p += size;
    size_t ps = ((size + 3) & ~0x3) - size;
    if (ps > 0) {
        uint32_t pad = 0;
        memcpy(p, &pad, ps);
        p += ps;
    }
    return p;
}

static void
make_srational_matrix(const m3x3 cm,
                      int32_t srational[18])
{
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            srational[2*(row*3+col)+0] = (int32_t)round(10000 * cm.v[row][col]);
            srational[2*(row*3+col)+1] = 10000;
        }
    }
}

static float *
make_float_hsmap(const uint32_t dims[],
                 const v3 *hsv)
{
    float *m = malloc(dims[0]*dims[1]*dims[2]*sizeof(*m)*3);
    for (uint32_t i = 0; i < dims[0]*dims[1]*dims[2]; i++) {
        m[3*i+0] = hsv[i].v[0];
        m[3*i+1] = hsv[i].v[1];
        m[3*i+2] = hsv[i].v[2];
    }
    return m;
}

void *
profio_dcp_toraw(const struct profio_dcp *dcp,
                 size_t *size)
{
    const int max_tag_count = 32;
    uint8_t hdr[8+2+max_tag_count * 12 + 4];
    uint8_t *p = hdr;
    {
        uint8_t bo[2];
        uint16_t t = 1;
        memcpy(bo, &t, 2);
        if (bo[1] == 1) {
            memcpy(p, "MMCR", 4);
        } else {
            memcpy(p, "IIRC", 4);
        }
        p+= 4;
    }
    uint32_t ifdoffset = 8;
    memcpy(p, &ifdoffset, 4);
    p += 4;
    uint16_t tag_count = 0;
    if (dcp->camera_model) tag_count++;
    tag_count += dcp->ill_count;
    tag_count += dcp->ill_count;
    if (dcp->calibration_signature) tag_count++;
    if (dcp->profile_name) tag_count++;
    if (dcp->huesatmap[0] != NULL) tag_count += 2;
    if (dcp->huesatmap[1] != NULL) tag_count++;
    if (dcp->has.huesatmap_srgb_gamma) tag_count++;
    if (dcp->tonecurve != NULL) tag_count++;
    if (dcp->has.embed_policy) tag_count++;
    if (dcp->copyright) tag_count++;
    if (dcp->has.fm) tag_count += dcp->ill_count;
    if (dcp->looktable != NULL) tag_count += 2;
    if (dcp->has.looktable_srgb_gamma) tag_count++;
    if (dcp->has.baseline_exposure_offset) tag_count++;
    if (dcp->has.black_render_none) tag_count++;
    debug("tag count: %d\n", tag_count);

    memcpy(p, &tag_count, 2);
    p += 2;

    uint32_t data_offset = 8 + 2 + tag_count * 12 + 4;
    if (dcp->camera_model) p = tiff_tag_toraw(p, TIFFTAG_UNIQUECAMERAMODEL, TIFF_ASCII, strlen(dcp->camera_model) + 1, dcp->camera_model, &data_offset);
    p = tiff_tag_toraw(p, TIFFTAG_COLORMATRIX1, TIFF_SRATIONAL, 9, NULL, &data_offset);
    if (dcp->ill_count == 2) {
        p = tiff_tag_toraw(p, TIFFTAG_COLORMATRIX2, TIFF_SRATIONAL, 9, NULL, &data_offset);
    }
    p = tiff_tag_toraw(p, TIFFTAG_CALIBRATIONILLUMINANT1, TIFF_SHORT, 1, &dcp->ill[0], NULL);
    if (dcp->ill_count == 2) {
        p = tiff_tag_toraw(p, TIFFTAG_CALIBRATIONILLUMINANT2, TIFF_SHORT, 1, &dcp->ill[1], NULL);
    }
    if (dcp->calibration_signature) p = tiff_tag_toraw(p, TIFFTAG_PROFILECALIBRATIONSIGNATURE, TIFF_ASCII, strlen(dcp->calibration_signature) + 1, dcp->calibration_signature, &data_offset);
    if (dcp->profile_name) p = tiff_tag_toraw(p, TIFFTAG_PROFILENAME, TIFF_ASCII, strlen(dcp->profile_name) + 1, dcp->profile_name, &data_offset);
    if (dcp->huesatmap[0] != NULL) {
        p = tiff_tag_toraw(p, TIFFTAG_PROFILEHUESATMAPDIMS, TIFF_LONG, 3, NULL, &data_offset);
        p = tiff_tag_toraw(p, TIFFTAG_PROFILEHUESATMAPDATA1, TIFF_FLOAT, dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]*3, NULL, &data_offset);
        if (dcp->huesatmap[1] != NULL) {
            p = tiff_tag_toraw(p, TIFFTAG_PROFILEHUESATMAPDATA2, TIFF_FLOAT, dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]*3, NULL, &data_offset);
        }
    }
    if (dcp->tonecurve) p = tiff_tag_toraw(p, TIFFTAG_PROFILETONECURVE, TIFF_FLOAT, dcp->tonecurve_len*2, NULL, &data_offset);
    uint32_t u32;
    if (dcp->has.embed_policy) {
        u32 = dcp->embed_policy;
        p = tiff_tag_toraw(p, TIFFTAG_PROFILEEMBEDPOLICY, TIFF_LONG, 1, &u32, &data_offset);
    }
    if (dcp->copyright) p = tiff_tag_toraw(p, TIFFTAG_PROFILECOPYRIGHT, TIFF_ASCII, strlen(dcp->copyright) + 1, dcp->copyright, &data_offset);
    if (dcp->has.fm) {
        p = tiff_tag_toraw(p, TIFFTAG_FORWARDMATRIX1, TIFF_SRATIONAL, 9, NULL, &data_offset);
        if (dcp->ill_count == 2) {
            p = tiff_tag_toraw(p, TIFFTAG_FORWARDMATRIX2, TIFF_SRATIONAL, 9, NULL, &data_offset);
        }
    }
    if (dcp->looktable != NULL) {
        p = tiff_tag_toraw(p, TIFFTAG_PROFILELOOKTABLEDIMS, TIFF_LONG, 3, NULL, &data_offset);
        p = tiff_tag_toraw(p, TIFFTAG_PROFILELOOKTABLEDATA, TIFF_FLOAT, dcp->lookdims[0]*dcp->lookdims[1]*dcp->lookdims[2]*3, NULL, &data_offset);
    }
    if (dcp->has.huesatmap_srgb_gamma) {
        u32 = !!dcp->huesatmap_srgb_gamma;
        p = tiff_tag_toraw(p, TIFFTAG_PROFILEHUESATMAPENCODING, TIFF_LONG, 1, &u32, NULL);
    }
    if (dcp->has.looktable_srgb_gamma) {
        u32 = !!dcp->looktable_srgb_gamma;
        p = tiff_tag_toraw(p, TIFFTAG_PROFILELOOKTABLEENCODING, TIFF_LONG, 1, &u32, NULL);
    }
    if (dcp->has.baseline_exposure_offset) {
        p = tiff_tag_toraw(p, TIFFTAG_BASELINEEXPOSUREOFFSET, TIFF_SRATIONAL, 1, NULL, &data_offset);
    }
    if (dcp->has.black_render_none) {
        u32 = !!dcp->black_render_none;
        p = tiff_tag_toraw(p, TIFFTAG_DEFAULTBLACKRENDER, TIFF_LONG, 1, &u32, &data_offset);
    }
    ifdoffset = 0;
    memcpy(p, &ifdoffset, 4); p+= 4;

    uint8_t *dcp_data = (uint8_t *)malloc(data_offset);
    *size = (size_t)data_offset;
    memcpy(dcp_data, hdr, (uintptr_t)p - (uintptr_t)hdr);
    p = &dcp_data[(uintptr_t)p - (uintptr_t)hdr];

    int32_t matrix[9*2];
    if (dcp->camera_model) p = tiff_data_toraw(p, dcp->camera_model, strlen(dcp->camera_model) + 1);
    make_srational_matrix(dcp->cm[0], matrix);
    p = tiff_data_toraw(p, matrix, 9*2*4);
    if (dcp->ill_count == 2) {
        make_srational_matrix(dcp->cm[1], matrix);
        p = tiff_data_toraw(p, matrix, 9*2*4);
    }
    if (dcp->calibration_signature) p = tiff_data_toraw(p, dcp->calibration_signature, strlen(dcp->calibration_signature) + 1);
    if (dcp->profile_name) p = tiff_data_toraw(p, dcp->profile_name, strlen(dcp->profile_name) + 1);
    if (dcp->huesatmap[0] != NULL) {
        float *hsmap = make_float_hsmap(dcp->hsmdims, dcp->huesatmap[0]);
        p = tiff_data_toraw(p, dcp->hsmdims, 3*4);
        p = tiff_data_toraw(p, hsmap, dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]*3*4);
        free(hsmap);
        if (dcp->huesatmap[1] != NULL) {
            hsmap = make_float_hsmap(dcp->hsmdims, dcp->huesatmap[1]);
            p = tiff_data_toraw(p, hsmap, dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]*3*4);
            free(hsmap);
        }
    }
    if (dcp->tonecurve) {
        float *tc = malloc(sizeof(*tc)*2*dcp->tonecurve_len);
        for (int i = 0; i < 2*dcp->tonecurve_len; i++) {
            tc[i] = dcp->tonecurve[i];
        }
        p = tiff_data_toraw(p, tc, dcp->tonecurve_len*2*4);
        free(tc);
    }
    if (dcp->copyright) p = tiff_data_toraw(p, dcp->copyright, strlen(dcp->copyright) + 1);
    if (dcp->has.fm) {
        make_srational_matrix(dcp->fm[0], matrix);
        p = tiff_data_toraw(p, matrix, 9*2*4);
        if (dcp->ill_count == 2) {
            make_srational_matrix(dcp->fm[1], matrix);
            p = tiff_data_toraw(p, matrix, 9*2*4);
        }
    }
    if (dcp->looktable != NULL) {
        float *hsmap = make_float_hsmap(dcp->lookdims, dcp->looktable);
        p = tiff_data_toraw(p, dcp->lookdims, 3*4);
        p = tiff_data_toraw(p, hsmap, dcp->lookdims[0]*dcp->lookdims[1]*dcp->lookdims[2]*3*4);
        free(hsmap);
    }
    if (dcp->has.baseline_exposure_offset) {
        int32_t srat[2] = { 10000 * dcp->baseline_exposure_offset, 10000 };
        p = tiff_data_toraw(p, srat, 2*4);
    }

    return dcp_data;
}

bool
profio_dcp_write(const struct profio_dcp *dcp,
                 const char filename[],
                 bool error_is_fatal)
{
    FILE *stream = utf8_fopen(filename, "wb");
    if (stream == NULL) {
        elog("Failed to open \"%s\" for writing: %s\n", filename, strerror(errno));
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    size_t size;
    void *data = profio_dcp_toraw(dcp, &size);
    fwrite(data, size, 1, stream);
    free(data);
    fclose(stream);
    return true;
}

const char *
exif_lightsource_ntoa(enum exif_lightsource ls)
{
    switch (ls) {
    case lsUnknown: return "Unknown";
    case lsDaylight: return "Daylight";
    case lsFluorescent: return "Fluorescent";
    case lsTungsten: return "Tungsten";
    case lsFlash: return "Flash";
    case lsFineWeather: return "FineWeather";
    case lsCloudyWeather: return "CloudyWeather";
    case lsShade: return "Shade";
    case lsDaylightFluorescent: return "DaylightFluorescent";	// D  5700 - 7100K
    case lsDayWhiteFluorescent: return "DayWhiteFluorescent";	// N  4600 - 5500K
    case lsCoolWhiteFluorescent: return "CoolWhiteFluorescent"; // W  3800 - 4500K
    case lsWhiteFluorescent: return "WhiteFluorescent"; // WW 3250 - 3800K
    case lsWarmWhiteFluorescent: return "WarmWhiteFluorescent"; // L  2600 - 3250K
    case lsStandardLightA: return "StdA";
    case lsStandardLightB: return "StdB";
    case lsStandardLightC: return "StdC";
    case lsD55: return "D55";
    case lsD65: return "D65";
    case lsD75: return "D75";
    case lsD50: return "D50";
    case lsISOStudioTungsten: return "ISOStudioTungsten";
    case lsOther: return "Other";
    default: return NULL;
    }
}

enum exif_lightsource
exif_lightsource_aton(const char str[])
{
    if (strcasecmp(str, "Unknown") == 0) {
        return lsUnknown;
    } else if (strcasecmp(str, "Daylight") == 0) {
        return lsDaylight;
    } else if (strcasecmp(str, "Fluorescent") == 0) {
        return lsFluorescent;
    } else if (strcasecmp(str, "Tungsten") == 0) {
        return lsTungsten;
    } else if (strcasecmp(str, "Flash") == 0) {
        return lsFlash;
    } else if (strcasecmp(str, "FineWeather") == 0) {
        return lsFineWeather;
    } else if (strcasecmp(str, "CloudyWeather") == 0) {
        return lsCloudyWeather;
    } else if (strcasecmp(str, "Shade") == 0) {
        return lsShade;
    } else if (strcasecmp(str, "DaylightFluorescent") == 0) {
        return lsDaylightFluorescent;
    } else if (strcasecmp(str, "DayWhiteFluorescent") == 0) {
        return lsDayWhiteFluorescent;
    } else if (strcasecmp(str, "CoolWhiteFluorescent") == 0) {
        return lsCoolWhiteFluorescent;
    } else if (strcasecmp(str, "WhiteFluorescent") == 0) {
        return lsWhiteFluorescent;
    } else if (strcasecmp(str, "WarmWhiteFluorescent") == 0) {
        return lsWarmWhiteFluorescent;
    } else if (strcasecmp(str, "Std A") == 0 || strcasecmp(str, "StdA") == 0) {
        return lsStandardLightA;
    } else if (strcasecmp(str, "Std B") == 0 || strcasecmp(str, "StdB") == 0) {
        return lsStandardLightB;
    } else if (strcasecmp(str, "Std C") == 0 || strcasecmp(str, "StdC") == 0) {
        return lsStandardLightC;
    } else if (strcasecmp(str, "D55") == 0) {
        return lsD55;
    } else if (strcasecmp(str, "D65") == 0) {
        return lsD65;
    } else if (strcasecmp(str, "D75") == 0) {
        return lsD75;
    } else if (strcasecmp(str, "D50") == 0) {
        return lsD50;
    } else if (strcasecmp(str, "ISOStudioTungsten") == 0) {
        return lsISOStudioTungsten;
    } else if (strcasecmp(str, "Other") == 0) {
        return lsOther;
    } else {
        char *ep = NULL;
        int ill = (int)strtol(str, &ep, 10);
        if (*ep != '\0' || exif_lightsource_ntoa(ill) == NULL) {
            return -1;
        }
        return (enum exif_lightsource)ill;
    }
    return -1;
}

double
exif_lightsource_temp(enum exif_lightsource ls)
{
    switch (ls) {
    case lsStandardLightA:
    case lsTungsten:
        return 2850.0;
    case lsISOStudioTungsten:
        return 3200.0;
    case lsD50:
        return 5000.0;
    case lsD55:
    case lsDaylight:
    case lsFineWeather:
    case lsFlash:
    case lsStandardLightB:
        return 5500.0;
    case lsD65:
    case lsStandardLightC:
    case lsCloudyWeather:
        return 6500.0;
    case lsD75:
    case lsShade:
        return 7500.0;
    case lsDaylightFluorescent:
        return (5700.0 + 7100.0) * 0.5;
    case lsDayWhiteFluorescent:
        return (4600.0 + 5500.0) * 0.5;
    case lsCoolWhiteFluorescent:
    case lsFluorescent:
        return (3800.0 + 4500.0) * 0.5;
    case lsWhiteFluorescent:
        return (3250.0 + 3800.0) * 0.5;
    case lsWarmWhiteFluorescent:
        return (2600.0 + 3250.0) * 0.5;
    default:
        return 0.0;
    }
}

char *
exif_lightsource_list(const char prefix[],
                      int lw)
{
    // in temperature order:
    const enum exif_lightsource ls[] = {
        lsStandardLightA,
        lsTungsten,
        lsWarmWhiteFluorescent,
        lsISOStudioTungsten,
        lsWhiteFluorescent,
        lsCoolWhiteFluorescent,
        lsFluorescent,
        lsD50,
        lsDayWhiteFluorescent,
        lsD55,
        lsDaylight,
        lsFineWeather,
        lsFlash,
        lsStandardLightB,
        lsDaylightFluorescent,
        lsD65,
        lsStandardLightC,
        lsCloudyWeather,
        lsD75,
        lsShade,
        lsOther
    };
    int len = 0;
    for (int i = 0; ls[i] != lsOther; i++) {
        int namlen = strlen(exif_lightsource_ntoa(ls[i]));
        len += namlen + 16 + strlen(prefix);
        if (lw < namlen + 10) {
            lw = namlen + 10;
        }
    }
    char *str = malloc(len);
    str[0] = '\0';
    int offset = 0;
    int last_lb = 0;
    offset += sprintf(str, "%s", prefix);
    for (int i = 0; ls[i] != lsOther; i++) {
        int cur = sprintf(&str[offset], "%s (%.0fK), ", exif_lightsource_ntoa(ls[i]), exif_lightsource_temp(ls[i]));
        if (offset + cur - last_lb > lw && offset > 0) {
            last_lb = offset;
            offset--;
            offset += sprintf(&str[offset], "\n%s%s (%.0fK), ", prefix, exif_lightsource_ntoa(ls[i]), exif_lightsource_temp(ls[i]));
        } else {
            offset += cur;
        }
    }
    str[offset-2] = '\0';
    return str;
}

void
profio_dcp_json(struct profio_dcp *dcp,
                FILE *stream)
{
    if (stream == NULL) {
        return;
    }
    char *json = profio_dcp_tojson(dcp);
    if (json == NULL) {
        return;
    }
    fprintf(stream, "%s", json);
    free(json);
}

char *
profio_dcp_tojson(struct profio_dcp *dcp)
{
    if (dcp == NULL) {
        return NULL;
    }
    strbuf_t sb = STRBUF_INITIALIZER;
    const char *pre_str = "";
    strbuf_asprintf(&sb, "{\n");
    if (dcp->camera_model) {
        strbuf_asprintf(&sb, "  \"UniqueCameraModel\": \"%s\"", dcp->camera_model);
        pre_str = ",\n";
    }
    if (dcp->profile_name) {
        strbuf_asprintf(&sb, "%s  \"ProfileName\": \"%s\"", pre_str, dcp->profile_name);
        pre_str = ",\n";
    }
    if (dcp->copyright) {
        strbuf_asprintf(&sb, "%s  \"ProfileCopyright\": \"%s\"", pre_str, dcp->copyright);
        pre_str = ",\n";
    }
    if (dcp->has.embed_policy) {
        const char *ep = "Allow Copying";
        switch (dcp->embed_policy) {
        case pepAllowCopying: ep = "Allow copying"; break;
        case pepEmbedIfUsed: ep = "Embed if used"; break;
        case pepEmbedNever: ep = "Embed never"; break;
        case pepNoRestrictions: ep = "No restrictions"; break;
        }
        strbuf_asprintf(&sb, "%s  \"ProfileEmbedPolicy\": \"%s\"", pre_str, ep);
        pre_str = ",\n";
    }
    if (dcp->calibration_signature) {
        strbuf_asprintf(&sb, "%s  \"ProfileCalibrationSignature\": \"%s\"", pre_str, dcp->calibration_signature);
        pre_str = ",\n";
    }
    for (int i = 0; i < dcp->ill_count; i++) {
        strbuf_asprintf(&sb, "%s  \"CalibrationIlluminant%d\": \"", pre_str, i+1);
        const char *ill = exif_lightsource_ntoa(dcp->ill[i]);
        if (ill == NULL) {
            strbuf_asprintf(&sb, "Undefined %d", dcp->ill[i]);
        } else {
            strbuf_asprintf(&sb, "%s", ill);
        }
        strbuf_asprintf(&sb, "\"");
        pre_str = ",\n";
    }
    for (int i = 0; i < dcp->ill_count; i++) {
        strbuf_asprintf(&sb,
                "%s"
                "  \"ColorMatrix%d\": [\n"
                "    [ %9.6f, %9.6f, %9.6f ],\n"
                "    [ %9.6f, %9.6f, %9.6f ],\n"
                "    [ %9.6f, %9.6f, %9.6f ]\n"
                "  ]",
                pre_str,
                i+1,
                dcp->cm[i].v[0][0], dcp->cm[i].v[0][1], dcp->cm[i].v[0][2],
                dcp->cm[i].v[1][0], dcp->cm[i].v[1][1], dcp->cm[i].v[1][2],
                dcp->cm[i].v[2][0], dcp->cm[i].v[2][1], dcp->cm[i].v[2][2]);
        pre_str = ",\n";
    }
    if (dcp->has.fm) {
        for (int i = 0; i < dcp->ill_count; i++) {
            strbuf_asprintf(&sb,
                    ",\n"
                    "  \"ForwardMatrix%d\": [\n"
                    "    [ %9.6f, %9.6f, %9.6f ],\n"
                    "    [ %9.6f, %9.6f, %9.6f ],\n"
                    "    [ %9.6f, %9.6f, %9.6f ]\n"
                    "  ]",
                    i+1,
                    dcp->fm[i].v[0][0], dcp->fm[i].v[0][1], dcp->fm[i].v[0][2],
                    dcp->fm[i].v[1][0], dcp->fm[i].v[1][1], dcp->fm[i].v[1][2],
                    dcp->fm[i].v[2][0], dcp->fm[i].v[2][1], dcp->fm[i].v[2][2]);
        }
    }
    if (dcp->has.black_render_none) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"DefaultBlackRender\": \"%s\"",
                dcp->black_render_none ? "None" : "Auto");
    }
    if (dcp->has.baseline_exposure_offset) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"BaselineExposureOffset\": %f",
                dcp->baseline_exposure_offset);
    }
    if (dcp->huesatmap[0] != NULL) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileHueSatMapDims\": [ %u, %u, %u ]",
                dcp->hsmdims[0], dcp->hsmdims[1], dcp->hsmdims[2]);
    }
    if (dcp->has.huesatmap_srgb_gamma) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileHueSatMapEncoding\": \"%s\"",
                dcp->huesatmap_srgb_gamma ? "sRGB" : "Linear");
    }
    if (dcp->looktable != NULL) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileLookTableDims\": [ %u, %u, %u ]",
                dcp->lookdims[0], dcp->lookdims[1], dcp->lookdims[2]);
    }
    if (dcp->has.looktable_srgb_gamma) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileLookTableEncoding\": \"%s\"",
                dcp->looktable_srgb_gamma ? "sRGB" : "Linear");
    }
    for (int i = 0; i < dcp->ill_count; i++) {
        if (dcp->huesatmap[i] == NULL) {
            continue;
        }
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileHueSatMap%d\": [\n",
                i+1);
        for (unsigned h = 0; h < dcp->hsmdims[0]; h++) {
            for (unsigned s = 0; s < dcp->hsmdims[1]; s++) {
                for (unsigned v = 0; v < dcp->hsmdims[2]; v++) {
                    const unsigned j = (dcp->hsmdims[0]*dcp->hsmdims[1])*v + dcp->hsmdims[1]*h + s;
                    strbuf_asprintf(&sb, "    { \"HueDiv\": %2u, \"SatDiv\": %2u, \"ValDiv\": %2u, \"HueShift\": %10.6f, \"SatScale\": %.6f, \"ValScale\": %.6f }",
                            h, s, v, dcp->huesatmap[i][j].v[0], dcp->huesatmap[i][j].v[1], dcp->huesatmap[i][j].v[2]);
                    if (j < dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2]-1) {
                        strbuf_asprintf(&sb, ",\n");
                    }
                }
            }
        }
        strbuf_asprintf(&sb, "\n  ]");
    }
    if (dcp->looktable != NULL) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileLookTable\": [\n");
        for (unsigned h = 0; h < dcp->lookdims[0]; h++) {
            for (unsigned s = 0; s < dcp->lookdims[1]; s++) {
                for (unsigned v = 0; v < dcp->lookdims[2]; v++) {
                    const unsigned j = (dcp->lookdims[0]*dcp->lookdims[1])*v + dcp->lookdims[1]*h + s;
                    strbuf_asprintf(&sb, "    { \"HueDiv\": %2u, \"SatDiv\": %2u, \"ValDiv\": %2u, \"HueShift\": %10.6f, \"SatScale\": %.6f, \"ValScale\": %.6f }",
                            h, s, v, dcp->looktable[j].v[0], dcp->looktable[j].v[1], dcp->looktable[j].v[2]);
                    if (j < dcp->lookdims[0]*dcp->lookdims[1]*dcp->lookdims[2]-1) {
                        strbuf_asprintf(&sb, ",\n");
                    }
                }
            }
        }
        strbuf_asprintf(&sb, "\n  ]");
    }
    if (dcp->tonecurve != NULL) {
        strbuf_asprintf(&sb,
                ",\n"
                "  \"ProfileToneCurve\": [\n");
        for (int j = 0; j < dcp->tonecurve_len; j++) {
            strbuf_asprintf(&sb, "    [ %f, %f ]", dcp->tonecurve[2*j+0], dcp->tonecurve[2*j+1]);
            if (j < dcp->tonecurve_len-1) {
                strbuf_asprintf(&sb, ",\n");
            }
        }
        strbuf_asprintf(&sb, "\n  ]");
    }
    strbuf_asprintf(&sb, "\n}\n");

    return strbuf_base(&sb);
}

void
profio_dcp_delete(struct profio_dcp *dcp)
{
    if (dcp == NULL) {
        return;
    }
    free(dcp->camera_model);
    free(dcp->profile_name);
    free(dcp->copyright);
    free(dcp->calibration_signature);
    free(dcp->huesatmap[0]);
    free(dcp->huesatmap[1]);
    free(dcp->looktable);
    free(dcp->tonecurve);
    free(dcp);
}

struct icctag {
    uint32_t sig;
    uint32_t offset;
    uint32_t size;
};

struct icchdr {
    uint32_t size;
    uint32_t cmm_sig;
    uint32_t version;
    uint32_t device_class;
    uint32_t input_space;
    uint32_t pcs;
    uint16_t create_date[6];
    uint32_t file_sig;
    uint32_t platform_sig;
    uint32_t cmm_flags;
    uint32_t dev_manufacturer;
    uint32_t dev_model;
    uint32_t dev_attr[2];
    uint32_t rendering_intent;
    uint32_t pcs_ill_xyz[3];
    uint32_t creator_sig;
    uint8_t reserved[44];
};

static void
labu16_tod(const uint16_t in[],
           double out[])
{
    out[0] = 100.0 * ((in[0] >> 8) / 255.0) + (in[0] & 0xFF) * 100.0 / 65280.0;
    out[1] = -128.0 + ((double)(in[1] >> 8)) + (in[1] & 0xFF) / 256.0;
    out[2] = -128.0 + ((double)(in[2] >> 8)) + (in[2] & 0xFF) / 256.0;
}

static void
dto_labu16(const double in[],
           uint16_t out[])
{
    int16_t num, frac;
    double w, val;
    val = in[0];
    if (val < 0) val = 0;
    if (val > 100.0 + 25500/65280.0) val = 100.0 + 25500/65280.0;
    val *= 2.55;
    w = floor(round(val * 100000) / 100000);
    num = (int16_t)w;
    frac = (int16_t)(((val - w) / 2.55 / (25500/65280.0)) * 256.0);
    if (frac > 0xFF) frac = 0xFF;
    out[0] = (((uint16_t)num) << 8) | frac;
    for (int i = 1; i < 3; i++) {
        val = in[i]+128.0;
        if (val < 0) val = 0;
        if (val > 255+255/256.0) val = 255+255/256.0;
        w = floor(val);
        num = (int16_t)w;
        frac = (int16_t)(((val - w) / (255.0/256.0)) * 256.0);
        if (frac > 0xFF) frac = 0xFF;
        out[i] = (((uint16_t)num) << 8) | frac;
    }
}

static double
u1fixed15_tod(uint16_t num)
{
    return (double)(num >> 15) + (num & 0x7FFF) / 32768.0;
}

static uint16_t
dto_u1fixed15(double val)
{
    double w = floor(val);
    int16_t num = (int16_t)w;
    int32_t frac = (int32_t)(((val - w) / (32767.0/32768.0)) * 32768.0);
    if (frac > 0x7FFF) frac = 0x7FFF;
    uint16_t out = (((uint16_t)num) << 15) | (uint16_t)frac;
    return out;
}

static double
s15fixed16_tod(uint32_t num)
{
    return (double)((int16_t)(num >> 16)) + (num & 0xFFFF) / 65536.0;
}

static uint32_t
dto_s15fixed16(double val)
{
    double w = floor(val);
    int16_t num = (int16_t)w;
    int32_t frac = (int32_t)(((val - w) / (65535.0/65536.0)) * 65536.0);
    if (frac > 0xFFFF) frac = 0xFFFF;
    uint32_t out = (((uint32_t)num) << 16) | frac;
    return out;
}

static uint16_t
dto_u8fixed8(double val)
{
    double w = floor(val);
    uint8_t num = (uint8_t)w;
    uint16_t frac = (uint16_t)(((val - w) / (255.0/256.0)) * 256.0);
    if (frac > 0xFF) frac = 0xFF;
    return (((uint16_t)num) << 8) | frac;
}

static bool
get_icc_xyztype(const uint8_t data[],
                const struct icctag tag,
                int count,
                v3 *dst)
{
    const uint8_t *p = &data[tag.offset];
    if (tag.size != 8 + 3 * 4 * (unsigned)count) {
        elog("Error: unexpected XYZType size.\n");
        return false;
    }
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t num;
            memcpy(&num, &p[8+12*i+j*4], 4);
            num = ntohl(num);
            dst[i].v[j] = s15fixed16_tod(num);
        }
    }
    return true;
}

static char *
get_icc_texttype(const uint8_t data[],
                 const struct icctag tag)
{
    const char *p = (char *)&data[tag.offset+8];
    if (tag.size < 8 + 1) {
        elog("Error: unexpected textType size.\n");
        return NULL;
    }
    if (p[tag.size-9] != '\0') {
        elog("Error: textType not null-terminated.\n");
        return NULL;
    }
    char *str = malloc(tag.size-8);
    memcpy(str, p, tag.size-8);
    return str;
}

static char *
get_icc_descriptiontype(const uint8_t data[],
                        const struct icctag tag)
{
    const char *p = (char *)&data[tag.offset+8];
    if (tag.size < 12 + 1) {
        elog("Error: unexpected textDescriptionType size.\n");
        return NULL;
    }
    uint32_t ascii_count;
    memcpy(&ascii_count, p, 4);
    ascii_count = ntohl(ascii_count);
    if (ascii_count == 0) {
        return strdup("");
    }
    if (tag.size < 12 + ascii_count) {
        elog("Error: unexpected textDescriptionType size.\n");
        return NULL;
    }
    if (p[4+ascii_count-1] != '\0') {
        elog("Error: ASCII not null-terminated.\n");
        return false;
    }
    char *str = malloc(ascii_count);
    memcpy(str, &p[4], ascii_count);
    return str;
}

static double *
get_icc_trc(const uint8_t data[],
            const struct icctag tag,
            int *len)
{
    const uint8_t *p = &data[tag.offset+8];
    if (tag.size < 12) {
        elog("Error: unexpected curveType size.\n");
        return NULL;
    }
    uint32_t count;
    memcpy(&count, p, 4);
    count = ntohl(count);

    double *tc;
    if (count == 0) {
        tc = malloc(2 * sizeof(tc[0]));
        tc[0] = 0;
        tc[1] = 1;
        *len = 2;
        return tc;
    }
    if (count == 1) {
        if (tag.size < 14) {
            elog("Error: unexpected curveType size.\n");
            return NULL;
        }
        // special case, gamma only;
        tc = malloc(1 * sizeof(tc[0]));
        uint16_t val;
        memcpy(&val, &p[4], 2);
        val = ntohs(val);
        tc[0] = (double)val / 256.0;
        *len = 1;
        return tc;
    }

    if (count * 2 > tag.size - 12) {
        elog("Error: unexpected curveType size.\n");
        return NULL;
    }
    tc = malloc(count * sizeof(tc[0]));
    const uint16_t *p16 = (const uint16_t *)&p[4];
    for (uint32_t j = 0; j < count; j++) {
        tc[j] = ntohs(p16[j]) / 65535.0;
    }
    *len = count;
    return tc;
}

static uint32_t
sig2num(const char sig[])
{
    return (((uint32_t)((uint8_t)sig[0]) << 24) | ((uint32_t)((uint8_t)sig[1]) << 16) | ((uint32_t)((uint8_t)sig[2]) << 8) | ((uint32_t)((uint8_t)sig[3]) << 0));
}

struct icc_lut *
get_icc_luttype(const uint8_t data[],
                const struct icctag tag,
                bool pcs_is_lab)
{
    const uint32_t *p32 = (const uint32_t *)&data[tag.offset];
    const uint16_t *p16 = (const uint16_t *)&data[tag.offset];
    const uint8_t *p8 = (const uint8_t *)&data[tag.offset];
    if (tag.size < 12) {
        elog("Error: unexpected lutType size.\n");
        return NULL;
    }
    if (ntohl(p32[0]) == sig2num("mft1")) {
        elog("Error: profile contains an 8 bit LUT, only supporting 16 bit LUT.\n");
        return NULL;
    }
    if (ntohl(p32[0]) != sig2num("mft2")) {
        elog("Error: only supporting 16 bit LUT, got some other type (0x%08X).\n", ntohl(p32[0]));
        return NULL;
    }
    unsigned in_ch = p8[8];
    unsigned out_ch = p8[9];
    if (in_ch != 3) {
        elog("Error: only exactly 3 LUT input channels supported, got %d.\n", in_ch);
        return NULL;
    }
    if (out_ch != 3) {
        elog("Error: only exactly 3 LUT output channels supported, got %d.\n", out_ch);
        return NULL;
    }
    int gridp = p8[10];
    unsigned in_tab = ntohs(p16[24]);
    unsigned out_tab = ntohs(p16[25]);
    if (in_tab < 2 || out_tab < 2 || gridp < 2) {
        elog("Error: unexpected tables sizes in lut16Type (%d %d %d).\n", in_tab, out_tab, gridp);
        return NULL;
    }
    if (tag.size != 52 + 2 * in_ch * in_tab + 2 * out_ch * out_tab + 2 * 3 * gridp * gridp * gridp) {
        elog("Error: unexpected lut16Type size.\n");
        return NULL;
    }
    struct icc_lut *lut = calloc(1, sizeof(*lut)+sizeof(lut->data[0])*(3*in_tab + 3*out_tab + 3*gridp*gridp*gridp));
    lut->cm.v[0][0] = s15fixed16_tod(ntohl(p32[3]));
    lut->cm.v[0][1] = s15fixed16_tod(ntohl(p32[4]));
    lut->cm.v[0][2] = s15fixed16_tod(ntohl(p32[5]));
    lut->cm.v[1][0] = s15fixed16_tod(ntohl(p32[6]));
    lut->cm.v[1][1] = s15fixed16_tod(ntohl(p32[7]));
    lut->cm.v[1][2] = s15fixed16_tod(ntohl(p32[8]));
    lut->cm.v[2][0] = s15fixed16_tod(ntohl(p32[9]));
    lut->cm.v[2][1] = s15fixed16_tod(ntohl(p32[10]));
    lut->cm.v[2][2] = s15fixed16_tod(ntohl(p32[11]));
    int idx = 0;
    lut->in_len = in_tab;
    for (unsigned i = 0; i < 3; i++) {
        lut->in[i] = &lut->data[idx];
        for (unsigned j = 0; j < in_tab; j++) {
            lut->data[idx] = ntohs(p16[26+idx]) / 65535.0;
            idx++;
        }
    }
    lut->clut_side = gridp;
    lut->clut = &lut->data[idx];
    for (int i = 0; i < 3*gridp*gridp*gridp; i++) {
        if (pcs_is_lab) {
            uint16_t in[3];
            in[0] = ntohs(p16[26+idx+0]);
            in[1] = ntohs(p16[26+idx+1]);
            in[2] = ntohs(p16[26+idx+2]);
            labu16_tod(in, &lut->data[idx]);
            i += 2;
            idx += 2;
        } else {
            // assume XYZ
            lut->data[idx] = u1fixed15_tod(ntohs(p16[26+idx]));
        }
        idx++;
    }
    lut->out_len = out_tab;
    for (unsigned i = 0; i < 3; i++) {
        lut->out[i] = &lut->data[idx];
        for (unsigned j = 0; j < out_tab; j++) {
            lut->data[idx] = ntohs(p16[26+idx]) / 65535.0;
            idx++;
        }
    }
    return lut;
}

struct profio_icc *
profio_icc_parse(const char filename[],
                 bool error_is_fatal)
{
    size_t data_len;
    uint8_t *data = fileio_file2mem(filename, &data_len);
    if (data == NULL) {
        elog("Error: could not read ICC profile \"%s\".\n", filename);
        if (error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    struct profio_icc *icc = profio_icc_fromraw(data, data_len);
    free(data);
    if (icc == NULL) {
        elog("Error: could not parse ICC profile \"%s\".\n", filename);
        if (error_is_fatal) exit(EXIT_FAILURE);
        return NULL;
    }
    return icc;
}

struct profio_icc *
profio_icc_fromraw(const uint8_t data[],
                   size_t len)
{
    struct profio_icc *icc = calloc(1, sizeof(*icc));
    if (data == NULL || len < 128 + 4) {
        elog("Error: file too small to contain an ICC profile.\n");
        goto fail;
    }
    uint32_t *swap = (uint32_t *)data;
    uint16_t *swap16 = (uint16_t *)data;
    for (int i = 0; i < 6; i++) {
        swap[i] = ntohl(swap[i]);
    }
    for (int i = 12; i < 12+6; i++) {
        swap16[i] = ntohs(swap16[i]);
    }
    for (int i = 9; i < 22; i++) {
        swap[i] = ntohl(swap[i]);
    }
    struct icchdr *hdr = (struct icchdr *)data;
    if ((hdr->version & 0xFF000000) != 0x02000000 || hdr->file_sig != 0x61637370) {
        if ((hdr->version & 0xFF000000) == 0x04000000) {
            elog("Error: ICC v4 profiles are not supported.\n");
            goto fail;
        }
        elog("Error: not an v2 ICC profile.\n");
        goto fail;
    }
    if (hdr->device_class != 0x73636E72) {
        elog("Error: the ICC is not for an input device, only camera profiles are supported.\n");
        goto fail;
    }
    if (hdr->input_space != 0x52474220) {
        elog("Error: the input color space is not RGB, only RGB camera profiles are supported.\n");
        goto fail;
    }
    if (hdr->pcs != 0x58595A20 && hdr->pcs != 0x4C616220) {
        elog("Error: PCS is not Lab or XYZ, which are the only PCSes supported.\n");
        goto fail;
    }
    if (hdr->size != len) {
        elog("Error: header size mismatch.\n");
        goto fail;
    }
    icc->pcs_is_lab = hdr->pcs == 0x4C616220;
    const uint8_t *p = &data[128]; // skip 128 byte header
    uint32_t tag_count;
    memcpy(&tag_count, p, 4);
    p += 4;
    tag_count = ntohl(tag_count);
    if (tag_count * 12 > len - 128 - 4 || tag_count > 1000) {
        elog("Error: bad tag count.\n");
        goto fail;
    }
    struct icctag *tags = (struct icctag *)p;
    for (unsigned i = 0; i < tag_count; i++, p += 12) {
        tags[i].sig = ntohl(tags[i].sig);
        tags[i].offset = ntohl(tags[i].offset);
        tags[i].size = ntohl(tags[i].size);
        if (tags[i].offset + tags[i].size > len) {
            elog("Error: bad tag offset.\n");
            goto fail;
        }
        if ((tags[i].offset & 3) != 0) {
            elog("Error: bad tag offset alignment.\n");
            goto fail;
        }
        char sig[5];
        sig[0] = (tags[i].sig >> 24) & 0xFF;
        sig[1] = (tags[i].sig >> 16) & 0xFF;
        sig[2] = (tags[i].sig >>  8) & 0xFF;
        sig[3] = (tags[i].sig >>  0) & 0xFF;
        sig[4] = '\0';
        if (strcmp(sig, "desc") == 0) {
            if ((icc->description = get_icc_descriptiontype(data, tags[i])) == NULL) {
                goto fail;
            }
        } else if (strcmp(sig, "cprt") == 0) {
            if ((icc->copyright = get_icc_texttype(data, tags[i])) == NULL) {
                goto fail;
            }
        } else if (strcmp(sig, "wtpt") == 0 || strcmp(sig, "bkpt") == 0) {
            v3 *v = sig[0] == 'w' ? &icc->wp : &icc->bp;
            if (!get_icc_xyztype(data, tags[i], 1, v)) {
                goto fail;
            }
        } else if (strcmp(sig, "A2B0") == 0) {
            if ((icc->lut = get_icc_luttype(data, tags[i], icc->pcs_is_lab)) == NULL) {
                goto fail;
            }
        } else if (strcmp(sig, "rXYZ") == 0 || strcmp(sig, "gXYZ") == 0 || strcmp(sig, "bXYZ") == 0) {
            int col = sig[0] == 'r' ? 0 : sig[0] == 'g' ? 1 : 2;
            v3 v;
            if (!get_icc_xyztype(data, tags[i], 1, &v)) {
                goto fail;
            }
            icc->fm.v[0][col] = v.v[0];
            icc->fm.v[1][col] = v.v[1];
            icc->fm.v[2][col] = v.v[2];
        } else if (strcmp(sig, "rTRC") == 0 || strcmp(sig, "gTRC") == 0 || strcmp(sig, "bTRC") == 0) {
            int idx = sig[0] == 'r' ? 0 : sig[0] == 'g' ? 1 : 2;
            int len = 0;
            double *trc = get_icc_trc(data, tags[i], &len);
            if (trc == NULL) {
                goto fail;
            }
            if (icc->trc_len != 0 && len != icc->trc_len) {
                elog("Error: varying TRC length, no support for that.\n");
                free(trc);
                goto fail;
            }
            if (icc->trc[idx] != NULL) {
                elog("Error: duplicate tag.\n");
                free(trc);
                goto fail;
            }
            icc->trc_len = len;
            icc->trc[idx] = trc;
        } else if (strcmp(sig, "kTRC") == 0) {
            elog("Error: found tag \"%s\" which means this is monochrome profile. Only color profiles are supported.\n", sig);
            goto fail;
        } else if (strcmp(sig, "vcgt") == 0 ||
                   strcmp(sig, "lumi") == 0 ||
                   strcmp(sig, "gamt") == 0 ||
                   strcmp(sig, "A2B2") == 0 ||
                   strcmp(sig, "A2B1") == 0 ||
                   strcmp(sig, "B2A2") == 0 ||
                   strcmp(sig, "B2A1") == 0 ||
                   strcmp(sig, "B2A0") == 0)
        {
            // tags that indicate this is not a camera profile
            elog("Error: found tag \"%s\" which means this is not a camera profile. Only ICCv2 camera profiles are supported.\n", sig);
            goto fail;
        } else if (strcmp(sig, "targ") == 0 ||
                   strcmp(sig, "DevD") == 0 ||
                   strcmp(sig, "CIED") == 0 ||
                   strcmp(sig, "SC20") == 0 ||
                   strcmp(sig, "tech") == 0 ||
                   strcmp(sig, "dmnd") == 0 ||
                   strcmp(sig, "dmdd") == 0 ||
                   strcmp(sig, "mmod") == 0 ||
                   strcmp(sig, "calt") == 0 ||
                   strcmp(sig, "clrt") == 0)
        {
            // silently ignore
        } else {
            elog("Ignoring unknown tag with signature \"%s\" (only a subset of standard tags is understood by this parser).\n", sig);
        }
    }
    if (icc->trc_len != 0 && (icc->trc[0] == NULL || icc->trc[1] == NULL || icc->trc[2] == NULL)) {
        elog("Error: missing TRC.\n");
        goto fail;
    }
    return icc;
fail:
    profio_icc_delete(icc);
    return NULL;
}

void
profio_icc_json(struct profio_icc *icc,
                FILE *stream)
{
    if (stream == NULL) {
        return;
    }
    char *json = profio_icc_tojson(icc);
    if (json == NULL) {
        return;
    }
    fprintf(stream, "%s", json);
    free(json);
}

char *
profio_icc_tojson(struct profio_icc *icc)
{
    if (icc == NULL) {
        return NULL;
    }
    const char *pre_str = "";
    strbuf_t sb = STRBUF_INITIALIZER;

    strbuf_asprintf(&sb, "{\n");
    if (icc->description) {
        strbuf_asprintf(&sb, "%s  \"Description\": \"%s\"", pre_str, icc->description);
        pre_str = ",\n";
    }
    if (icc->copyright) {
        strbuf_asprintf(&sb, "%s  \"Copyright\": \"%s\"", pre_str, icc->copyright);
        pre_str = ",\n";
    }
    if (!memiszero(&icc->wp, sizeof(icc->wp))) {
        strbuf_asprintf(&sb,
                "%s  \"Whitepoint\": [ %9.6f, %9.6f, %9.6f ]",
                pre_str,
                icc->wp.v[0], icc->wp.v[1], icc->wp.v[2]);
        pre_str = ",\n";
    }
    if (!memiszero(&icc->bp, sizeof(icc->bp))) {
        strbuf_asprintf(&sb,
                "%s  \"Blackpoint\": [ %9.6f, %9.6f, %9.6f ]",
                pre_str,
                icc->bp.v[0], icc->bp.v[1], icc->bp.v[2]);
        pre_str = ",\n";
    }
    if (!memiszero(&icc->fm, sizeof(icc->fm))) {
        strbuf_asprintf(&sb,
                "%s  \"ForwardMatrix\": [\n"
                "    [ %9.6f, %9.6f, %9.6f ],\n"
                "    [ %9.6f, %9.6f, %9.6f ],\n"
                "    [ %9.6f, %9.6f, %9.6f ]\n"
                "  ]",
                pre_str,
                icc->fm.v[0][0], icc->fm.v[0][1], icc->fm.v[0][2],
                icc->fm.v[1][0], icc->fm.v[1][1], icc->fm.v[1][2],
                icc->fm.v[2][0], icc->fm.v[2][1], icc->fm.v[2][2]);
        pre_str = ",\n";
    }
    {
        strbuf_asprintf(&sb, "%s  \"ProfileConnectionSpace\": \"%s\"", pre_str, icc->pcs_is_lab ? "Lab" : "XYZ");
        pre_str = ",\n";
    }
    if (icc->trc_len > 0) {
        const char *trcname[3] = { "RedTRC", "GreenTRC", "BlueTRC" };
        for (int i = 0; i < 3; i++) {
            strbuf_asprintf(&sb,
                    "%s  \"%s\":", pre_str, trcname[i]);
            if (icc->trc_len > 1) {
                strbuf_asprintf(&sb, " [\n");
                for (int j = 0; j < icc->trc_len; j++) {
                    strbuf_asprintf(&sb, "    %f", icc->trc[i][j]);
                    if (j < icc->trc_len-1) {
                        strbuf_asprintf(&sb, ",\n");
                    }
                }
                strbuf_asprintf(&sb, "\n  ]");
            } else {
                strbuf_asprintf(&sb, " %f", icc->trc[i][0]);
            }
            pre_str = ",\n";
        }
    }
    if (icc->lut != NULL) {
        const struct icc_lut *lut = icc->lut;
        strbuf_asprintf(&sb,
                "%s  \"CLUT\": {\n"
                "    \"InputMatrix\": [\n"
                "      [ %9.6f, %9.6f, %9.6f ],\n"
                "      [ %9.6f, %9.6f, %9.6f ],\n"
                "      [ %9.6f, %9.6f, %9.6f ]\n"
                "    ],\n",
                pre_str,
                lut->cm.v[0][0], lut->cm.v[0][1], lut->cm.v[0][2],
                lut->cm.v[1][0], lut->cm.v[1][1], lut->cm.v[1][2],
                lut->cm.v[2][0], lut->cm.v[2][1], lut->cm.v[2][2]);

        strbuf_asprintf(&sb, "    \"InputCurves\": [\n");
        for (int i = 0; i < lut->in_len; i++) {
            strbuf_asprintf(&sb, "      [ %f, %f, %f ]", lut->in[0][i], lut->in[1][i], lut->in[2][i]);
            if (i < lut->in_len-1) strbuf_asprintf(&sb, ",\n");
        }
        strbuf_asprintf(&sb, "\n    ],\n");

        strbuf_asprintf(&sb, "    \"CLUT\": [\n");
        const int cs = lut->clut_side;
        int j = 0;
        for (int h = 0; h < cs; h++) {
            for (int s = 0; s < cs; s++) {
                for (int v = 0; v < cs; v++, j++) {
                    strbuf_asprintf(&sb, "      { \"Div1\": %2u, \"Div2\": %2u, \"Div3\": %2u, \"Ch1\": %.6f, \"Ch2\": %.6f, \"Ch3\": %.6f }",
                            h, s, v, lut->clut[3*(cs*cs*h+cs*s+v)+0], lut->clut[3*(cs*cs*h+cs*s+v)+1], lut->clut[3*(cs*cs*h+cs*s+v)+2]);
                    if (j < cs*cs*cs-1) {
                        strbuf_asprintf(&sb, ",\n");
                    }
                }
            }
        }
        strbuf_asprintf(&sb, "\n    ],\n");

        strbuf_asprintf(&sb, "    \"OutputCurves\": [\n");
        for (int i = 0; i < lut->out_len; i++) {
            strbuf_asprintf(&sb, "      [ %f, %f, %f ]", lut->out[0][i], lut->out[1][i], lut->out[2][i]);
            if (i < lut->out_len-1) strbuf_asprintf(&sb, ",\n");
        }
        strbuf_asprintf(&sb, "\n    ]");
        strbuf_asprintf(&sb, "\n  }");
        pre_str = ",\n";
    }
    strbuf_asprintf(&sb, "\n}\n");

    return strbuf_base(&sb);
}

static uint16_t
dto_u16(double in)
{
    int32_t val = (int32_t)round(in * 65535.0);
    if (val < 0) val = 0;
    if (val > 65535) val = 65535;
    return (uint16_t)val;
}

void *
profio_icc_toraw(const struct profio_icc *icc,
                 size_t *size)
{
    struct icchdr hdr;
    memset(&hdr, 0, sizeof(hdr));

    hdr.cmm_sig = htonl(0);
    hdr.version = htonl(0x02100000);
    hdr.device_class = htonl(0x73636E72); // scnr
    hdr.input_space = htonl(0x52474220); // RGB
    hdr.pcs = icc->pcs_is_lab ? htonl(0x4C616220) : htonl(0x58595A20); // Lab or XYZ
    time_t tt = time(NULL);
    struct tm tm;
#ifdef __MINGW32__
    localtime_s(&tm, &tt);
#else
    localtime_r(&tt, &tm);
#endif
    hdr.create_date[0] = htons(tm.tm_year + 1900);
    hdr.create_date[1] = htons(tm.tm_mon+1);
    hdr.create_date[2] = htons(tm.tm_mday);
    hdr.create_date[3] = htons(tm.tm_hour);
    hdr.create_date[4] = htons(tm.tm_min);
    hdr.create_date[5] = htons(tm.tm_sec);
    hdr.file_sig = htonl(0x61637370);
    hdr.platform_sig = htonl(0);
    hdr.cmm_flags = htonl(0);
    hdr.dev_manufacturer = htonl(0);
    hdr.dev_model = htonl(0);
    hdr.dev_attr[0] = htonl(0);
    hdr.dev_attr[1] = htonl(0);
    hdr.rendering_intent = htonl(0);
    hdr.pcs_ill_xyz[0] = htonl(0x0000F6D6);
    hdr.pcs_ill_xyz[1] = htonl(0x00010000);
    hdr.pcs_ill_xyz[2] = htonl(0x0000D32D);
    hdr.creator_sig = htonl(0);

    int tag_count = 0;
    struct icctag tags[32];
    void *tag_data[32];
    uint32_t *td;
    memset(tags, 0, sizeof(tags));
    if (icc->description != NULL) {
        tags[tag_count].sig = sig2num("desc");
        tags[tag_count].size = 90 + strlen(icc->description) + 1;
        td = calloc(1, tags[tag_count].size);
        td[0] = htonl(sig2num("desc"));
        td[1] = htonl(0);
        td[2] = htonl(strlen(icc->description)+1);
        strcpy((char *)&td[3], icc->description);
        tag_data[tag_count] = td;
        tag_count++;
    }
    if (icc->copyright != NULL) {
        tags[tag_count].sig = sig2num("cprt");
        tags[tag_count].size = 8 + strlen(icc->copyright) + 1;
        td = calloc(1, tags[tag_count].size);
        td[0] = htonl(sig2num("text"));
        td[1] = htonl(0);
        strcpy((char *)&td[2], icc->copyright);
        tag_data[tag_count] = td;
        tag_count++;
    }
    if (!memiszero(&icc->wp, sizeof(icc->wp))) {
        tags[tag_count].sig = sig2num("wtpt");
        tags[tag_count].size = 20;
        td = calloc(1, tags[tag_count].size);
        td[0] = htonl(sig2num("XYZ "));
        td[1] = htonl(0);
        td[2] = htonl(dto_s15fixed16(icc->wp.v[0]));
        td[3] = htonl(dto_s15fixed16(icc->wp.v[1]));
        td[4] = htonl(dto_s15fixed16(icc->wp.v[2]));
        tag_data[tag_count] = td;
        tag_count++;
    }
    if (!memiszero(&icc->bp, sizeof(icc->bp))) {
        tags[tag_count].sig = sig2num("bkpt");
        tags[tag_count].size = 20;
        td = calloc(1, tags[tag_count].size);
        td[0] = htonl(sig2num("XYZ "));
        td[1] = htonl(0);
        td[2] = htonl(dto_s15fixed16(icc->bp.v[0]));
        td[3] = htonl(dto_s15fixed16(icc->bp.v[1]));
        td[4] = htonl(dto_s15fixed16(icc->bp.v[2]));
        tag_data[tag_count] = td;
        tag_count++;
    }
    if (!memiszero(&icc->fm, sizeof(icc->fm))) {
        for (int i = 0; i < 3; i++) {
            const char *sig = i == 0 ? "rXYZ" : i == 1 ? "gXYZ" : "bXYZ";
            tags[tag_count].sig = sig2num(sig);
            tags[tag_count].size = 20;
            td = calloc(1, tags[tag_count].size);
            td[0] = htonl(sig2num("XYZ "));
            td[1] = htonl(0);
            td[2] = htonl(dto_s15fixed16(icc->fm.v[0][i]));
            td[3] = htonl(dto_s15fixed16(icc->fm.v[1][i]));
            td[4] = htonl(dto_s15fixed16(icc->fm.v[2][i]));
            tag_data[tag_count] = td;
            tag_count++;
        }
    }
    if (icc->trc_len > 0) {
        for (int i = 0; i < 3; i++) {
            const char *sig = i == 0 ? "rTRC" : i == 1 ? "gTRC" : "bTRC";
            tags[tag_count].sig = sig2num(sig);
            uint16_t val;
            if (icc->trc_len == 1) {
                tags[tag_count].size = 14;
                td = calloc(1, tags[tag_count].size);
                td[0] = htonl(sig2num("curv"));
                td[1] = htonl(0);
                td[2] = htonl(1);
                val = htons(dto_u8fixed8(icc->trc[i][0]));
                memcpy(&td[3], &val, 2);
            } else {
                tags[tag_count].size = 12 + 2 * icc->trc_len;
                td = calloc(1, tags[tag_count].size);
                td[0] = htonl(sig2num("curv"));
                td[1] = htonl(0);
                td[2] = htonl(icc->trc_len);
                uint16_t *td16 = (uint16_t *)&td[3];
                for (int j = 0; j < icc->trc_len; j++) {
                    td16[j] = htons(dto_u16(icc->trc[i][j]));
                }
            }
            tag_data[tag_count] = td;
            tag_count++;
        }
    }
    if (icc->lut != NULL) {
        const struct icc_lut *lut = icc->lut;
        const int cs = lut->clut_side;
        tags[tag_count].sig = sig2num("A2B0");
        tags[tag_count].size = 52 + 2*3*lut->in_len + 2*3*lut->out_len + 2*3*cs*cs*cs;
        td = calloc(1, tags[tag_count].size);
        td[0] = htonl(sig2num("mft2"));
        td[1] = htonl(0);
        uint16_t *p16 = (uint16_t *)td;
        uint8_t *p8 = (uint8_t *)td;
        p8[8] = 3;
        p8[9] = 3;
        p8[10] = cs;
        td[3] =  htonl(dto_s15fixed16(lut->cm.v[0][0]));
        td[4] =  htonl(dto_s15fixed16(lut->cm.v[0][1]));
        td[5] =  htonl(dto_s15fixed16(lut->cm.v[0][2]));
        td[6] =  htonl(dto_s15fixed16(lut->cm.v[1][0]));
        td[7] =  htonl(dto_s15fixed16(lut->cm.v[1][1]));
        td[8] =  htonl(dto_s15fixed16(lut->cm.v[1][2]));
        td[9] =  htonl(dto_s15fixed16(lut->cm.v[2][0]));
        td[10] = htonl(dto_s15fixed16(lut->cm.v[2][1]));
        td[11] = htonl(dto_s15fixed16(lut->cm.v[2][2]));
        p16[24] = htons(lut->in_len);
        p16[25] = htons(lut->out_len);
        int idx = 26;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < lut->in_len; j++) {
                p16[idx++] = htons(dto_u16(lut->in[i][j]));
            }
        }
        for (int i = 0; i < 3*cs*cs*cs; i++) {
            if (icc->pcs_is_lab) {
                dto_labu16(&lut->clut[i], &p16[idx]);
                p16[idx+0] = htons(p16[idx+0]);
                p16[idx+1] = htons(p16[idx+1]);
                p16[idx+2] = htons(p16[idx+2]);
                idx += 3;
                i += 2;
            } else {
                p16[idx++] = htons(dto_u1fixed15(lut->clut[i]));
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < lut->out_len; j++) {
                p16[idx++] = htons(dto_u16(lut->out[i][j]));
            }
        }
        tag_data[tag_count] = td;
        tag_count++;
    }

    uint32_t offset = 128 + 4 + 12 * tag_count;
    for (int i = 0; i < tag_count; i++) {
        bool own_data = true;
        for (int j = 0; j < tag_count; j++) {
            if (i != j && tags[j].offset != 0 && tags[j].size == tags[i].size && memcmp(tag_data[j], tag_data[i], tags[i].size) == 0) {
                tags[i].offset = tags[j].offset;
                free(tag_data[i]);
                tag_data[i] = NULL;
                own_data = false;
                break;
            }
        }
        if (own_data) {
            tags[i].offset = offset;
            offset += tags[i].size;
            offset = (offset + 3) & ~3;
        }
    }
    for (int i = 0; i < tag_count; i++) {
        tags[i].sig = htonl(tags[i].sig);
        tags[i].offset = htonl(tags[i].offset);
        tags[i].size = htonl(tags[i].size);
    }

    hdr.size = htonl(offset);
    *size = offset;
    void *raw = malloc(offset);
    {
        uint32_t *p32 = (uint32_t *)raw;
        memcpy(raw, &hdr, sizeof(hdr));
        p32[32] = htonl(tag_count);
        memcpy(&p32[33], tags, 12 * tag_count);
        offset = 128 + 4 + 12 * tag_count;
        for (int i = 0; i < tag_count; i++) {
            uint32_t sz = ntohl(tags[i].size);
            if (tag_data[i] != NULL) {
                memcpy(&((uint8_t *)raw)[offset], tag_data[i], sz);
                offset += sz;
                int pad = ((sz + 3) & ~3) - sz;
                if (pad > 0) {
                    memset(&((uint8_t *)raw)[offset], 0, pad);
                    offset += pad;
                }
                free(tag_data[i]);
            }
        }
    }
    return raw;
}

bool
profio_icc_write(const struct profio_icc *icc,
                 const char filename[],
                 bool error_is_fatal)
{
    FILE *stream = utf8_fopen(filename, "wb");
    if (stream == NULL) {
        elog("Failed to open \"%s\" for writing: %s\n", filename, strerror(errno));
        if (error_is_fatal) exit(EXIT_FAILURE);
        return false;
    }
    size_t size;
    void *data = profio_icc_toraw(icc, &size);
    fwrite(data, size, 1, stream);
    free(data);
    fclose(stream);
    return true;
}

void
profio_icc_delete(struct profio_icc *icc)
{
    if (icc == NULL) {
        return;
    }
    free(icc->description);
    free(icc->copyright);
    free(icc->trc[0]);
    free(icc->trc[1]);
    free(icc->trc[2]);
    free(icc->lut);
    free(icc);
}
