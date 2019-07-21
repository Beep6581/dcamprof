/*
 * (c) Copyright 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <lut3d.h>

// based on pseudo code from http://www.filmlight.ltd.uk/pdf/whitepapers/FL-TL-TN-0057-SoftwareLib.pdf
static double
tetrahedral_interpolation(v3 abc,
                          const double *lut,
                          const int output_channel_count,
                          const int output_channel_index,
                          const int cube_side)
{
    const int cs = cube_side;

    // truncate range to 0...1 to avoid going outside the table
    abc = v3_clip01(abc);

    // change range from 0...1 to 0...cs-1
    double a = abc.v[0] * (cs-1);
    double b = abc.v[1] * (cs-1);
    double c = abc.v[2] * (cs-1);

    // prev integer value, truncate down
    int pa = (int)a;
    int pb = (int)b;
    int pc = (int)c;

    // next integer value
    int na = pa+1;
    int nb = pb+1;
    int nc = pc+1;
    if (na > cs-1) na = cs-1;
    if (nb > cs-1) nb = cs-1;
    if (nc > cs-1) nc = cs-1;

    // fraction value, that is how close real value is to next integer value
    double fa = a - pa;
    double fb = b - pb;
    double fc = c - pc;

#define LUT(a,b,c) lut[output_channel_count*(cs*cs*a+cs*b+c)+output_channel_index]

    double res;
    double c000 = LUT(pa, pb, pc);
    double c111 = LUT(na, nb, nc);
    if (fa > fb) {
        if (fb > fc) {
            double c100 = LUT(na, pb, pc);
            double c110 = LUT(na, nb, pc);
            res = (1-fa) * c000 + (fa-fb) * c100 + (fb-fc) * c110 + (fc) * c111;
        } else if (fa > fc) {
            double c100 = LUT(na, pb, pc);
            double c101 = LUT(na, pb, nc);
            res = (1-fa) * c000 + (fa-fc) * c100 + (fc-fb) * c101 + (fb) * c111;
        } else {
            double c001 = LUT(pa, pb, nc);
            double c101 = LUT(na, pb, nc);
            res = (1-fc) * c000 + (fc-fa) * c001 + (fa-fb) * c101 + (fb) * c111;
        }
    } else {
        if (fc > fb) {
            double c001 = LUT(pa, pb, nc);
            double c011 = LUT(pa, nb, nc);
            res = (1-fc) * c000 + (fc-fb) * c001 + (fb-fa) * c011 + (fa) * c111;
        } else if (fc > fa) {
            double c010 = LUT(pa, nb, pc);
            double c011 = LUT(pa, nb, nc);
            res = (1-fb) * c000 + (fb-fc) * c010 + (fc-fa) * c011 + (fa) * c111;
        } else {
            double c010 = LUT(pa, nb, pc);
            double c110 = LUT(na, nb, pc);
            res = (1-fb) * c000 + (fb-fa) * c010 + (fa-fc) * c110 + (fc) * c111;
        }
    }
    return res;

#undef LUT
}

double
lut3d_lookup_1d(const double *lut,
                int cs,
                v3 abc)
{
    return tetrahedral_interpolation(abc, lut, 1, 0, cs);
}
