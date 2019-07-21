/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef COLMATH_H_
#define COLMATH_H_

#include <math.h>

#include <spectrum.h>
#include <matmath.h>

enum exif_lightsource {
    lsUnknown = 0,
    lsDaylight = 1,
    lsFluorescent = 2,
    lsTungsten = 3,
    lsFlash = 4,
    lsFineWeather = 9,
    lsCloudyWeather = 10,
    lsShade = 11,
    lsDaylightFluorescent = 12,	// D  5700 - 7100K
    lsDayWhiteFluorescent = 13,	// N  4600 - 5500K
    lsCoolWhiteFluorescent = 14, // W  3800 - 4500K
    lsWhiteFluorescent = 15, // WW 3250 - 3800K
    lsWarmWhiteFluorescent = 16, // L  2600 - 3250K
    lsStandardLightA = 17,
    lsStandardLightB = 18,
    lsStandardLightC = 19,
    lsD55 = 20,
    lsD65 = 21,
    lsD75 = 22,
    lsD50 = 23,
    lsISOStudioTungsten = 24,
    lsOther = 255
};

struct observer {
    spectrum_t *cmf[3];
};

/*
  rgb2xyz_adobergb
     +0.5767 +0.1856 +0.1882
     +0.2973 +0.6274 +0.0753
     +0.0270 +0.0707 +0.9913

  rgb2xyz_srgb
     +0.4124 +0.3576 +0.1805
     +0.2126 +0.7152 +0.0722
     +0.0193 +0.1192 +0.9505
*/
static const m3x3 xyz2rgb_adobergb = {{{  2.04159, -0.56501, -0.34473 },
                                       { -0.96924,  1.87597,  0.04156 },
                                       {  0.01344, -0.11836,  1.01517 }}};
static const m3x3 xyz2rgb_srgb = {{{  3.2406, -1.5372, -0.4986 },
                                   { -0.9689,  1.8758,  0.0415 },
                                   {  0.0557, -0.2040,  1.0570 }}};

static inline double
srgb_gamma_forward(const double x)
{
    if (x <= 0.0031308) {
        return x * 12.92;
    } else {
        return 1.055 * pow (x, 1.0 / 2.4) - 0.055;
    }
}

static inline double
srgb_gamma_inverse(const double y)
{
    if (y <= 0.0031308 * 12.92) {
        return y * (1.0 / 12.92);
    } else {
        return pow ((y + 0.055) * (1.0 / 1.055), 2.4);
    }
}

static inline v3
wp_from_rgb2xyz(const m3x3 rgb2xyz)
{
    return (v3){{ rgb2xyz.v[0][0] + rgb2xyz.v[0][1] + rgb2xyz.v[0][2], rgb2xyz.v[1][0] + rgb2xyz.v[1][1] + rgb2xyz.v[1][2], rgb2xyz.v[2][0] + rgb2xyz.v[2][1] + rgb2xyz.v[2][2] }};
}

v3
xyz2srgb(v3 xyz, v3 wp);

uint32_t
xyz2srgb_u24(v3 xyz, v3 wp);

static inline double
lab_ft_forward(double t)
{
    if (t >= 8.85645167903563082e-3) {
        return pow(t, 1.0/3.0);
    } else {
        return t * (841.0/108.0) + 4.0/29.0;
    }
}

static inline double
lab_ft_inverse(double t)
{
    if (t >= 0.206896551724137931) {
        return t*t*t;
    } else {
        return 108.0 / 841.0 * (t - 4.0/29.0);
    }
}

const char *
xyz2name(v3 xyz, v3 wp, char buf[]);

v3
xyz2lab(v3 xyz, v3 wp);

v3
lab2xyz(v3 lab, v3 wp);

v3
lab2lch(v3 lab);

v3
lch2lab(v3 lch);

v3
hsl2rgb(v3 hsl);

v3
rgb2hsl(v3 rgb);

v3
hsv2rgb(v3 hsv);

v3
rgb2hsv(v3 rgb);

static inline v3
xyz2lch(v3 xyz, v3 wp)
{
    return lab2lch(xyz2lab(xyz, wp));
}

static inline v3
lch2xyz(v3 lch, v3 wp)
{
    return lab2xyz(lch2lab(lch), wp);
}

static inline v3
xyz2xyY(const v3 xyz)
{
    double x = xyz.v[0], y = xyz.v[1], z = xyz.v[2];
    return (v3){{ x / (x+y+z), y / (x+y+z), y }};
}

static inline v3
xyY2xyz(const v3 xyY)
{
    double x = xyY.v[0], y = xyY.v[1], Y = xyY.v[2];
    return (v3){{ x*Y / y, Y, (1-x-y)*Y/y }};
}

static inline v3
xyz2uvY(const v3 xyz)
{
    double x = xyz.v[0], y = xyz.v[1], z = xyz.v[2];
    return (v3){{ 4*x / (x + 15*y + 3*z), 9*y / (x + 15*y + 3*z), y }};
}

static inline v3
uvY2xyz(const v3 uvY)
{
    double u = uvY.v[0], v = uvY.v[1], Y = uvY.v[2];
    return (v3){{ Y*9*u / (4*v), Y, Y * (12 - 3*u - 20*v) / (4*v) }};
}

v3
xyz2luv(v3 xyz, v3 wp);

v3
luv2xyz(v3 luv, v3 wp);

v3
xyz2ipt(v3 xyz, v3 wp);

v3
ipt2xyz(v3 ipt, v3 wp);

v3
xyz2ipt_polar(v3 xyz, v3 wp);

v3
ipt2xyz_polar(v3 ipt_polar, v3 wp);

static inline v3
jab2jch(v3 jab)
{
    return lab2lch(jab);
}

static inline v3
jch2jab(v3 jch)
{
    return lch2lab(jch);
}

v3
jch2ucs(v3 jch);

v3
ucs2jch(v3 jab);

v3
xyz2ucs(v3 xyz, v3 wp);

v3
ucs2xyz(v3 jab, v3 wp);

v3
xyz2ucs_s(v3 xyz, void *h02); // h02 with La == 20 is required

v3
ucs2xyz_s(v3 jab, void *h02); // h02 with La == 20 is required

v3
xyz2jch(v3 xyz, v3 wp);

v3
jch2xyz(v3 jch, v3 wp);

v3
xyz2jch_s(v3 xyz, void *h02);

v3
jch2xyz_s(v3 jch, void *h02);

static inline v3
rgb2jch_s(const v3 rgb, const m3x3 rgb2xyz, void *h02)
{
    return xyz2jch_s(m3x3_mul_v3(rgb2xyz, rgb), h02);
}

static inline v3
jch2rgb_s(const v3 jch, const m3x3 xyz2rgb, void *h02)
{
    return m3x3_mul_v3(xyz2rgb, jch2xyz_s(jch, h02));
}

v3
jch_blend(v3 jch_a, v3 jch_b, double vs); // vs == 1 => 100% of jch_a, vs == 0 => 100% of jch_b

static inline v3
lch_blend(v3 lch_a, v3 lch_b, double vs)
{
    return jch_blend(lch_a, lch_b, vs);
}

double
spec2xyz_weight(const struct observer *o, const spectrum_t *illuminant);

v3
spec2xyz_ill(const spectrum_t *s, const struct observer *o, const spectrum_t *illuminant, double weight);

v3
spec2xyz(const spectrum_t *s, const struct observer *o);

v3
xy2xyz_lightest(double x,
                double y,
                const v3 wp);

v3
uv2xyz_lightest(double u,
                double v,
                const v3 wp);

m3x3
bradford_map_white(v3 src_wp,
                   v3 dst_wp);

enum ca_transform {
    CAT_CAT02,
    CAT_BRADFORD
};

v3
chromatic_adaptation_transform(enum ca_transform cat,
                               const v3 xyz,
                               const v3 src_wp,
                               const v3 dst_wp);

double
ciede2000(v3 xyz1, v3 xyz2, v3 wp);

double
ciede2000_lch(v3 xyz1, v3 xyz2, v3 wp,
              double Kl, double Kc, double Kh);

double
ciede2000_lch_wlc(v3 xyz1, v3 xyz2, v3 wp,
                  double Kl, double Kc, double Kh,
                  double wl, double wc, double wh);

double
euclidean_dist(v3 v1, v3 v2);

#endif
