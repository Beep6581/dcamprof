/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef SPECTRALDB_H_
#define SPECTRALDB_H_

#include <stdbool.h>

#include <spectrum.h>
#include <colmath.h>

#define COLOR_NAME_MAX 8

const spectrum_t *
spectraldb_color(const char group[],
                 const char name[]);

const spectrum_t *
spectraldb_color_byidx(const char group[],
                       int idx,
                       char name[]); // optional, COLOR_NAME_MAX

int
spectraldb_group_count(const char group[]);

const spectrum_t *
spectraldb_illuminant(enum exif_lightsource ls);

const spectrum_t *
spectraldb_illuminant_byname(const char name[],
                             bool missing_is_fatal);

const char *
spectraldb_illuminant_list(void);

spectrum_t *
spectraldb_render_blackbody(double temp,
                            int band_start,
                            int band_end,
                            int band_spacing);

#endif
