/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef TPS_H_
#define TPS_H_

#include <matmath.h>

struct tps_t_;
typedef struct tps_t_ tps_t;

tps_t *
tps_new(const v3 control_points[],
        size_t control_points_len,
        double regularization);

void
tps_delete(tps_t *tps);

double
tps_bending_energy(tps_t *tps);

double
tps_interpolate(tps_t *tps,
                double x,
                double z);

#endif
