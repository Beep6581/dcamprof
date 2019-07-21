/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <bisection.h>

double
bisection(double (*func)(double x, void *arg),
          void *arg,
          double a,
          double b,
          double tol,
          int nmax)
{
    double c = 0;
    for (int i = 0; i < nmax; i++) {
        c = (a + b) * 0.5;
        double f_c = func(c, arg);
        if (f_c == 0 || (b - a) * 0.5 < tol) {
            return c;
        }
        double f_a = func(a, arg);
        if ((f_c < 0 && f_a < 0) || (f_c > 0 && f_a > 0)) {
            a = c;
        } else {
            b = c;
        }
    }
    return c;
}

double
interval_halving_min(double (*func)(double x, void *arg),
                     void *arg,
                     double a,
                     double b,
                     double tol,
                     int nmax)
{
    double L = b - a;
    double x = (a + b) * 0.5;
    for (int i = 0; i < nmax; i++) {
        double f_x = func(x, arg);
        if ((b - a) * 0.5 < tol) {
            return x;
        }
        double x1 = a + L/4;
        double f_x1 = func(x1, arg);
        if (f_x1 < f_x) {
            b = x;
            x = x1;
        } else {
            double x2 = b - L/4;
            double f_x2 = func(x2, arg);
            if (f_x2 < f_x) {
                a = x;
                x = x2;
            } else {
                a = x1;
                b = x2;
            }
        }
        L = b - a;
    }
    return x;
}
