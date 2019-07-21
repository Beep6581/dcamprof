/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef XYZ2SPEC_H_
#define XYZ2SPEC_H_

#include <colmath.h>

spectrum_t *
xyz2spec(v3 xyz,
         const struct observer *obs,
         const spectrum_t *illuminant);

#endif
