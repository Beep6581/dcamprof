/*
 * (c) Copyright 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef LUT3D_H_
#define LUT3D_H_

#include <matmath.h>

double
lut3d_lookup_1d(const double *lut,
                int cs,
                v3 abc);

v3
lut3d_lookup(const double *lut,
             int cs,
             v3 abc);

#endif
