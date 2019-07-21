/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef BISECTION_H_
#define BISECTION_H_

double
bisection(double (*func)(double x, void *arg),
          void *arg,
          double low_lim,
          double high_lim,
          double tolerance,
          int max_iterations);

double
interval_halving_min(double (*func)(double x, void *arg),
                     void *arg,
                     double low_lim,
                     double high_lim,
                     double tolerance,
                     int max_iterations);

#endif

