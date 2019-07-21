/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef TIFIO_H_
#define TIFIO_H_

#include <stdbool.h>

#include <colmath.h>

struct rgbimg {
    int w;
    int h;
    uint16_t p[];
};

enum tifio_colorspace {
    COLORSPACE_NONE, // that is sRGB
    COLORSPACE_PROPHOTO
};

void
tifio_type(struct rgbimg *img,
           const char str[],
           uint16_t bg,
           uint16_t fg,
           int x,
           int y);

struct rgbimg *
tifio_rgbimg_load(const char filename[],
                  bool error_is_fatal);

bool
tifio_rgbimg_save(const char filename[],
                  const struct rgbimg *img,
                  enum tifio_colorspace colorspace,
                  bool error_is_fatal);

bool
tifio_get_transfer_function(const char filename[],
                            double *trc[3],
                            int *count_,
                            bool error_is_fatal);

bool
tifio_flatfield_correct(const char orig_filename[],
                        const char ff_filename[],
                        const char out_filename[],
                        bool error_is_fatal);


// SPB format
struct spb_image {
    uint32_t w;
    uint32_t h;
    uint32_t band_count;
    float band_start;
    float band_spacing;
    float band_end;
    float p[];
};

struct spb_image *
tifio_spb_load(char filename[],
               int band_base,
               int band_width,
               bool error_is_fatal);

#endif
