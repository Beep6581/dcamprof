/*
 * (c) Copyright 2015 - 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#include <inttypes.h>
#include <lcms2.h>

#include <colmath.h>
#include <dngref.h>
#include <interp.h>
#include <bisection.h>
#include <elog.h>

v3
xyz2srgb(v3 xyz,
         v3 wp)
{
    v3 lab = xyz2lab(xyz, wp);
    xyz = lab2xyz(lab, (v3){{0.950456,1.0,1.088754}}); // D65

    double X = xyz.v[0], Y = xyz.v[1], Z = xyz.v[2];

    // matrix xyz to srgb
    double r =  3.2406 * X - 1.5372 * Y - 0.4986 * Z;
    double g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
    double b =  0.0557 * X - 0.2040 * Y + 1.0570 * Z;

    double mn = r;
    if (g < mn) mn = g;
    if (b < mn) mn = b;

    // clip
    if (mn < 0) {
        r -= mn;
        g -= mn;
        b -= mn;
    }
    if (r > 1.0) r = 1.0;
    if (g > 1.0) g = 1.0;
    if (b > 1.0) b = 1.0;

    r = srgb_gamma_forward(r);
    g = srgb_gamma_forward(g);
    b = srgb_gamma_forward(b);
    return (v3){{ r, g, b }};
}

uint32_t
xyz2srgb_u24(v3 xyz,
             v3 wp)
{
    v3 rgb = xyz2srgb(xyz, wp);
    double m = rgb.v[0];
    if (rgb.v[1] > m) m = rgb.v[1];
    if (rgb.v[2] > m) m = rgb.v[2];
    if (m > 1.0) {
        rgb.v[0] /= m;
        rgb.v[1] /= m;
        rgb.v[2] /= m;
    }
    uint32_t r = (int)(255*rgb.v[0]);
    uint32_t g = (int)(255*rgb.v[1]);
    uint32_t b = (int)(255*rgb.v[2]);
    return (r << 16) | (g << 8) | b;
}

static const struct {
    const char *hue;
    double x, y;
} munsell_d50_hues[] = {
    { "BB", 0.217131, 0.315868 },{ "BB", 0.275236, 0.347227 },{ "BB", 0.255686, 0.341913 },{ "BB", 0.298299, 0.353967 },{ "BB", 0.323967, 0.357194 },{ "BG", 0.258165, 0.396637 },{ "BG", 0.279862, 0.393782 },{ "BG", 0.311970, 0.378013 },{ "BG", 0.297305, 0.386500 },{ "BG", 0.327192, 0.371531 },{ "GG", 0.278950, 0.490437 },{ "GG", 0.305379, 0.455827 },{ "GG", 0.324514, 0.419209 },{ "GG", 0.319540, 0.431871 },{ "GG", 0.339541, 0.380845 },{ "GG", 0.336353, 0.392202 },{ "GY", 0.407816, 0.470689 },{ "GY", 0.420957, 0.496964 },{ "GY", 0.427532, 0.512579 },{ "GY", 0.366782, 0.392454 },{ "GY", 0.382432, 0.418067 },{ "GY", 0.393575, 0.440349 },{ "PB", 0.207014, 0.243848 },{ "PB", 0.239980, 0.275875 },{ "PB", 0.265236, 0.300240 },{ "PB", 0.308592, 0.335909 },{ "PB", 0.296567, 0.325932 },{ "PB", 0.332922, 0.350082 },{ "PP", 0.304965, 0.241597 },{ "PP", 0.307228, 0.260952 },{ "PP", 0.325928, 0.309731 },{ "PP", 0.322069, 0.300891 },{ "PP", 0.332970, 0.329410 },{ "PP", 0.342517, 0.348298 },{ "RP", 0.402653, 0.289468 },{ "RP", 0.384608, 0.314777 },{ "RP", 0.390647, 0.308870 },{ "RP", 0.362413, 0.341006 },{ "RP", 0.372502, 0.331839 },{ "RP", 0.353167, 0.352562 },{ "RR", 0.508579, 0.321605 },{ "RR", 0.452868, 0.338031 },{ "RR", 0.480443, 0.334220 },{ "RR", 0.425255, 0.345933 },{ "RR", 0.384199, 0.354444 },{ "RR", 0.404360, 0.352322 },{ "RR", 0.361289, 0.357174 },{ "YR", 0.540115, 0.391552 },{ "YR", 0.551470, 0.390229 },{ "YR", 0.453526, 0.384401 },{ "YR", 0.477741, 0.388540 },{ "YR", 0.499720, 0.393156 },{ "YR", 0.398224, 0.375064 },{ "YR", 0.425159, 0.381713 },{ "YR", 0.370963, 0.367278 },{ "YY", 0.494669, 0.467295 },{ "YY", 0.505413, 0.471458 },{ "YY", 0.425934, 0.419658 },{ "YY", 0.445713, 0.432631 },{ "YY", 0.464976, 0.447532 },{ "YY", 0.472587, 0.452932 },{ "YY", 0.373952, 0.381755 },{ "YY", 0.400350, 0.402783 },{ "BB", 0.231807, 0.311689 },{ "BB", 0.276401, 0.336919 },{ "BB", 0.257974, 0.329090 },{ "BB", 0.299768, 0.347419 },{ "BB", 0.324169, 0.354651 },{ "BG", 0.254657, 0.377117 },{ "BG", 0.288186, 0.376536 },{ "BG", 0.271895, 0.379882 },{ "BG", 0.306849, 0.371851 },{ "BG", 0.325097, 0.367590 },{ "GG", 0.298283, 0.423702 },{ "GG", 0.292228, 0.432206 },{ "GG", 0.324970, 0.391943 },{ "GG", 0.313851, 0.409734 },{ "GG", 0.336205, 0.376167 },{ "GY", 0.403852, 0.533601 },{ "GY", 0.381755, 0.446133 },{ "GY", 0.393906, 0.472894 },{ "GY", 0.398337, 0.495523 },{ "GY", 0.361849, 0.393853 },{ "GY", 0.374441, 0.419863 },{ "PB", 0.227170, 0.246541 },{ "PB", 0.253714, 0.275055 },{ "PB", 0.272827, 0.295552 },{ "PB", 0.313073, 0.333215 },{ "PB", 0.302665, 0.323606 },{ "PB", 0.334894, 0.349297 },{ "PP", 0.325963, 0.248159 },{ "PP", 0.328672, 0.267426 },{ "PP", 0.334443, 0.313025 },{ "PP", 0.335695, 0.305557 },{ "PP", 0.339110, 0.330557 },{ "PP", 0.344491, 0.349427 },{ "RP", 0.433282, 0.302753 },{ "RP", 0.402140, 0.324577 },{ "RP", 0.407568, 0.321139 },{ "RP", 0.368779, 0.344680 },{ "RP", 0.384076, 0.339124 },{ "RP", 0.356081, 0.353957 },{ "RR", 0.533032, 0.332706 },{ "RR", 0.486339, 0.342869 },{ "RR", 0.433828, 0.352886 },{ "RR", 0.462123, 0.350787 },{ "RR", 0.388238, 0.357837 },{ "RR", 0.405594, 0.356815 },{ "RR", 0.362002, 0.359078 },{ "YR", 0.482709, 0.402199 },{ "YR", 0.508020, 0.405092 },{ "YR", 0.519706, 0.411761 },{ "YR", 0.399165, 0.380049 },{ "YR", 0.431022, 0.388896 },{ "YR", 0.453647, 0.397276 },{ "YR", 0.372653, 0.370282 },{ "YY", 0.438274, 0.448105 },{ "YY", 0.454033, 0.462071 },{ "YY", 0.465448, 0.474594 },{ "YY", 0.479370, 0.482212 },{ "YY", 0.374509, 0.387951 },{ "YY", 0.396849, 0.410628 },{ "YY", 0.416465, 0.428461 },{ "BB", 0.227870, 0.296933 },{ "BB", 0.277909, 0.331666 },{ "BB", 0.254419, 0.317190 },{ "BB", 0.302873, 0.344761 },{ "BB", 0.326861, 0.353112 },{ "BG", 0.238509, 0.356068 },{ "BG", 0.282545, 0.365287 },{ "BG", 0.268422, 0.366386 },{ "BG", 0.305710, 0.365400 },{ "BG", 0.325028, 0.362353 },{ "GG", 0.270407, 0.425938 },{ "GG", 0.289817, 0.412053 },{ "GG", 0.318631, 0.388274 },{ "GG", 0.306406, 0.400921 },{ "GG", 0.331668, 0.375429 },{ "GY", 0.359706, 0.564143 },{ "GY", 0.362019, 0.444232 },{ "GY", 0.365566, 0.467597 },{ "GY", 0.366143, 0.496220 },{ "GY", 0.353765, 0.389167 },{ "GY", 0.357179, 0.410872 },{ "PB", 0.244333, 0.228092 },{ "PB", 0.279303, 0.276817 },{ "PB", 0.304409, 0.310417 },{ "PB", 0.297416, 0.301458 },{ "PB", 0.319884, 0.329987 },{ "PB", 0.337229, 0.347599 },{ "PP", 0.349810, 0.257085 },{ "PP", 0.349427, 0.289906 },{ "PP", 0.352028, 0.304991 },{ "PP", 0.348610, 0.334147 },{ "PP", 0.353907, 0.323979 },{ "PP", 0.348424, 0.351140 },{ "RP", 0.475495, 0.297232 },{ "RP", 0.436511, 0.320651 },{ "RP", 0.447731, 0.314764 },{ "RP", 0.412527, 0.332089 },{ "RP", 0.375118, 0.349162 },{ "RP", 0.391302, 0.342148 },{ "RP", 0.358824, 0.354053 },{ "RR", 0.540646, 0.345592 },{ "RR", 0.562856, 0.340977 },{ "RR", 0.497544, 0.356299 },{ "RR", 0.442399, 0.363412 },{ "RR", 0.464682, 0.361874 },{ "RR", 0.390405, 0.363311 },{ "RR", 0.419314, 0.364436 },{ "RR", 0.362673, 0.360376 },{ "YR", 0.499934, 0.421473 },{ "YR", 0.513693, 0.427648 },{ "YR", 0.524820, 0.429330 },{ "YR", 0.400974, 0.386229 },{ "YR", 0.424770, 0.398347 },{ "YR", 0.451001, 0.408442 },{ "YR", 0.472411, 0.416351 },{ "YR", 0.373492, 0.373748 },{ "YY", 0.428938, 0.456054 },{ "YY", 0.443866, 0.471543 },{ "YY", 0.454617, 0.486150 },{ "YY", 0.373334, 0.390779 },{ "YY", 0.392448, 0.412029 },{ "YY", 0.411583, 0.435333 },{ "BB", 0.208824, 0.262791 },{ "BB", 0.233262, 0.287720 },{ "BB", 0.263465, 0.311240 },{ "BB", 0.308077, 0.341504 },{ "BB", 0.286892, 0.330279 },{ "BB", 0.332383, 0.353241 },{ "BG", 0.227651, 0.336552 },{ "BG", 0.278578, 0.355484 },{ "BG", 0.266784, 0.355264 },{ "BG", 0.300206, 0.360236 },{ "BG", 0.326522, 0.360626 },{ "GG", 0.267699, 0.410029 },{ "GG", 0.284616, 0.402156 },{ "GG", 0.316248, 0.382511 },{ "GG", 0.300779, 0.392824 },{ "GG", 0.330536, 0.373102 },{ "GY", 0.318129, 0.535048 },{ "GY", 0.333000, 0.490118 },{ "GY", 0.340818, 0.431695 },{ "GY", 0.339856, 0.449877 },{ "GY", 0.346798, 0.383689 },{ "GY", 0.346905, 0.406571 },{ "PB", 0.275211, 0.232007 },{ "PB", 0.295436, 0.275769 },{ "PB", 0.316057, 0.309503 },{ "PB", 0.308142, 0.299076 },{ "PB", 0.326840, 0.329737 },{ "PB", 0.341227, 0.349397 },{ "PP", 0.378492, 0.270042 },{ "PP", 0.372990, 0.291614 },{ "PP", 0.367307, 0.307914 },{ "PP", 0.356219, 0.337398 },{ "PP", 0.360752, 0.324247 },{ "PP", 0.351284, 0.351398 },{ "RP", 0.498562, 0.307463 },{ "RP", 0.440859, 0.329614 },{ "RP", 0.462048, 0.324584 },{ "RP", 0.423364, 0.339854 },{ "RP", 0.379474, 0.350973 },{ "RP", 0.395976, 0.345590 },{ "RP", 0.357988, 0.356193 },{ "RR", 0.574340, 0.362111 },{ "RR", 0.505622, 0.373079 },{ "RR", 0.528235, 0.373996 },{ "RR", 0.449983, 0.374560 },{ "RR", 0.472041, 0.376113 },{ "RR", 0.394610, 0.369676 },{ "RR", 0.422940, 0.374169 },{ "RR", 0.364045, 0.362187 },{ "YR", 0.426276, 0.411566 },{ "YR", 0.453148, 0.426090 },{ "YR", 0.475271, 0.433532 },{ "YR", 0.493931, 0.444480 },{ "YR", 0.510888, 0.450853 },{ "YR", 0.375179, 0.379737 },{ "YR", 0.401646, 0.396725 },{ "YY", 0.421603, 0.462009 },{ "YY", 0.433660, 0.481951 },{ "YY", 0.440810, 0.493134 },{ "YY", 0.372024, 0.392012 },{ "YY", 0.390211, 0.416857 },{ "YY", 0.404777, 0.441432 },{ "BB", 0.335194, 0.356566 },{ "BG", 0.335646, 0.363125 },{ "PB", 0.338992, 0.352641 },{ "BB", 0.337050, 0.354498 },{ "BG", 0.333896, 0.360057 },{ "PB", 0.342496, 0.352884 },{ "GY", 0.354238, 0.374797 },{ "GY", 0.345006, 0.373898 },{ "GG", 0.339079, 0.368827 },{ "GG", 0.337609, 0.365871 },{ "PP", 0.345490, 0.353237 },{ "PP", 0.348832, 0.355137 },{ "RP", 0.351078, 0.356184 },{ "RP", 0.352843, 0.356822 },{ "RR", 0.353762, 0.358448 },{ "RR", 0.356648, 0.361029 },{ "YR", 0.361153, 0.365548 },{ "YR", 0.362710, 0.370318 },{ "YY", 0.360116, 0.371916 },{ "YY", 0.354666, 0.370461 },{ "RR", 0.479434, 0.322451 },{ "RR", 0.467709, 0.330151 },{ "YR", 0.503992, 0.379380 },{ "YR", 0.537582, 0.380709 },{ "RR", 0.508005, 0.327769 },{ "RR", 0.480649, 0.339177 },{ "YR", 0.526492, 0.398552 },{ "YR", 0.501312, 0.399218 },{ "RR", 0.535365, 0.340137 },{ "RR", 0.505244, 0.351693 },{ "YR", 0.497134, 0.414701 },{ "YR", 0.515626, 0.418647 },{ "YR", 0.524312, 0.420654 },{ "RR", 0.550114, 0.351307 },{ "RR", 0.494330, 0.364766 },{ "RR", 0.534662, 0.363767 },{ "YR", 0.489765, 0.429750 },{ "YR", 0.503800, 0.434350 },{ "YY", 0.478396, 0.443929 },{ "YY", 0.489905, 0.451697 },{ "YY", 0.499474, 0.455609 },{ "GY", 0.431502, 0.493054 },{ "GG", 0.304899, 0.502606 },{ "YY", 0.477816, 0.467157 },{ "YY", 0.463155, 0.457451 },{ "YY", 0.452412, 0.470727 },{ "YY", 0.443618, 0.482227 },{ "GY", 0.341344, 0.519733 },{ "PB", 0.222957, 0.253551 },{ "PB", 0.215342, 0.213268 },{ "PB", 0.240200, 0.244630 },{ "PB", 0.263457, 0.234409 },{ "RP", 0.392725, 0.267129 },{ "RP", 0.418746, 0.297476 },{ "RP", 0.460474, 0.294397 },{ "RP", 0.440402, 0.309330 },{ "RP", 0.488310, 0.302537 },{ "RP", 0.459460, 0.318325 }
};

const char *
xyz2name(v3 xyz, v3 wp, char buf[])
{
    static char static_buf[64];
    if (buf == NULL) {
        buf = static_buf;
    }
    const char *hue_name = NULL;
    {
        v3 d50 = {{ cmsD50X, cmsD50Y, cmsD50Z }};
        xyz = m3x3_mul_v3(bradford_map_white(wp, d50), xyz);
        v3 xyY = xyz2xyY(xyz);
        int min_e = 0;
        double min_de = 1000000000;
        for (size_t i = 0; i < sizeof(munsell_d50_hues)/sizeof(munsell_d50_hues[0]); i++) {
            double x = xyY.v[0];
            double y = xyY.v[1];
            double x1 = munsell_d50_hues[i].x;
            double y1 = munsell_d50_hues[i].y;
            double de = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
            if (de < min_de) {
                min_e = i;
                min_de = de;
            }
        }
        const char *mhue = munsell_d50_hues[min_e].hue;
        if (strcmp(mhue, "RR") == 0) {
            hue_name = "red";
        } else if (strcmp(mhue, "YR") == 0) {
            hue_name = "orange";
        } else if (strcmp(mhue, "YY") == 0) {
            hue_name = "yellow";
        } else if (strcmp(mhue, "GY") == 0) {
            hue_name = "yellow-green";
        } else if (strcmp(mhue, "GG") == 0) {
            hue_name = "green";
        } else if (strcmp(mhue, "BG") == 0) {
            hue_name = "cyan";
        } else if (strcmp(mhue, "BB") == 0) {
            hue_name = "blue";
        } else if (strcmp(mhue, "PB") == 0) {
            hue_name = "purple-blue";
        } else if (strcmp(mhue, "PP") == 0) {
            hue_name = "purple";
        } else if (strcmp(mhue, "RP") == 0) {
            hue_name = "purple-red";
        } else {
            elog("unexpected name %s\n", mhue);
            abort();
        }
    }

    v3 lch;
    {
        v3 lab = xyz2lab(xyz, wp);
        cmsCIELab lab1 = { .L = lab.v[0], .a = lab.v[1], .b = lab.v[2] };
        cmsCIELCh lch1;
        cmsLab2LCh(&lch1, &lab1);
        lch = (v3){{ lch1.L, lch1.C, lch1.h }};
    }

    char *p = buf;
    *p = '\0';
    int lightness = (int)lch.v[0];
    int chroma = (int)lch.v[1];
    //int hue = (int)lch.v[2];
    //p += sprintf(p, "%3d %2d %3d ", lightness, chroma, hue);
    char *base = p;

    if (lightness > 70 && chroma > 5 && chroma <= 20) {
        // light grayish
        if (lightness >= 90) {
            p += sprintf(p, "whitish ");
        } else {
            p += sprintf(p, "pale ");
        }
    } else {
        if (chroma > 5) {
            if (lightness <= 15) {
            p += sprintf(p, "very dark ");
            } else if (lightness <= 40) {
                p += sprintf(p, "dark ");
            } else if (lightness <= 70) {
                // medium
            } else if (lightness <= 90) {
                p += sprintf(p, "light ");
            } else {
                p += sprintf(p, "very light ");
            }
        }
        if (chroma <= 5) {
            if (lightness > 90) {
                strcpy(p, "white");
                return buf;
            } else if (lightness <= 5) {
                strcpy(p, "black");
                return buf;
            }
            p += sprintf(p, "gray %d%%", 10 * (int)(round(lightness / 10.0)));
            return buf;
        } else if (chroma <= 20) {
            p += sprintf(p, "grayish ");
        } else if (chroma <= 60) {
            // moderate
            //p += sprintf(p, "moderate ");
        } else if (chroma <= 80) {
            p += sprintf(p, "strong ");
        } else if (chroma <= 100) {
            p += sprintf(p, "vivid ");
        } else {
            p += sprintf(p, "very vivid ");
        }
    }
    strcat(p, hue_name);
    // special cases
    if (strcmp(base, "pale red") == 0) {
        strcpy(base, "pink");
    } else if (strcmp(base, "whitish red") == 0) {
        strcpy(base, "whitish pink");
    } else if ((p = strstr(base, "grayish orange")) != NULL) {
        strcpy(p, "brown");
    }
    return buf;
}

v3
xyz2lab(const v3 xyz1_,
        const v3 wp_)
{
    cmsCIEXYZ xyz1 = { .X = xyz1_.v[0], .Y = xyz1_.v[1], .Z = xyz1_.v[2] };
    cmsCIELab lab1;
    cmsCIEXYZ wp = { .X = wp_.v[0], .Y = wp_.v[1], .Z = wp_.v[2] };
    cmsXYZ2Lab(&wp, &lab1, &xyz1);
    return (v3){{ lab1.L, lab1.a, lab1.b }};
}

v3
lab2xyz(const v3 lab1_,
        const v3 wp_)
{
    cmsCIELab lab1 = { .L = lab1_.v[0], .a = lab1_.v[1], .b = lab1_.v[2] };
    cmsCIEXYZ xyz1;
    cmsCIEXYZ wp = { .X = wp_.v[0], .Y = wp_.v[1], .Z = wp_.v[2] };
    cmsLab2XYZ(&wp, &xyz1, &lab1);
    return (v3){{ xyz1.X, xyz1.Y, xyz1.Z }};
}

v3
lab2lch(v3 lab)
{
    cmsCIELab lab1 = { .L = lab.v[0], .a = lab.v[1], .b = lab.v[2] };
    cmsCIELCh lch1;
    cmsLab2LCh(&lch1, &lab1);
    return (v3){{ lch1.L, lch1.C, lch1.h }};
}

v3
lch2lab(v3 lch)
{
    cmsCIELCh lch1 = { .L = lch.v[0], .C = lch.v[1], .h = lch.v[2] };
    cmsCIELab lab1;
    cmsLCh2Lab(&lab1, &lch1);
    return (v3){{ lab1.L, lab1.a, lab1.b }};
}

v3
xyz2jch_s(v3 xyz, void *h02)
{
    cmsCIEXYZ XYZ = (cmsCIEXYZ){ .X = xyz.v[0] * 100.0, .Y = xyz.v[1] * 100.0, .Z = xyz.v[2] * 100.0 };
    cmsJCh JCh;
    cmsCIECAM02Forward((cmsHANDLE *)h02, &XYZ, &JCh);
    return (v3){{ JCh.J, JCh.C, JCh.h }};
}

v3
jch2xyz_s(v3 jch, void *h02)
{
    //ciecam02 reverse can't handle 0 lightness (division by J=0 in the formula), we make a special case handling for that here
    if (jch.v[0] == 0 && jch.v[1] == 0) {
        return (v3){{ 0, 0, 0 }};
    }
    cmsJCh JCh = (cmsJCh){ jch.v[0], jch.v[1], jch.v[2] };
    cmsCIEXYZ XYZ;
    cmsCIECAM02Reverse(h02, &JCh, &XYZ);
    return (v3){{ XYZ.X / 100.0, XYZ.Y / 100.0, XYZ.Z / 100.0 }};
}

v3
jch_blend(v3 jch_a, v3 jch_b, double vs)
{
    if (vs >= 1) {
        return jch_a;
    }
    if (vs <= 0) {
        return jch_b;
    }
    v3 jch;
    jch.v[0] = jch_a.v[0] * vs + jch_b.v[0] * (1 - vs);
    jch.v[1] = jch_a.v[1] * vs + jch_b.v[1] * (1 - vs);
    double h1 = jch_a.v[2];
    double h2 = jch_b.v[2];
    double d = h2 - h1;
    vs = 1 - vs;
    if (h1 > h2) {
        double tmp = h1;
        h1 = h2;
        h2 = tmp;
        d = -d;
        vs = 1 - vs;
    }
    double h;
    if (d > 180) {
        h1 = h1 + 360;
        h = fmod((h1 + vs * (h2 - h1)), 360);
    } else {
        h = h1 + vs * d;
    }
    jch.v[2] = h;
    return jch;
}

v3
luv2xyz(const v3 luv,
        const v3 wp)
{
    double L = luv.v[0], u = luv.v[1], v = luv.v[2], x, y, z;
    y = (L + 16) / 116;
    y = wp.v[1] * lab_ft_inverse(y);
    if (L != 0) {
        u /= L;
        v /= L;
    }
    u = u/13 + 4*wp.v[0] / (wp.v[0] + 15*wp.v[1] + 3*wp.v[2]);
    v = v/13 + 9*wp.v[1] / (wp.v[0] + 15*wp.v[1] + 3*wp.v[2]);
    x = (y) * ((9*u)/(4*v));
    z = (y) * ((3 - 0.75*u)/v - 5);
    return (v3){{ x, y, z }};
}

v3
xyz2luv(const v3 xyz,
        const v3 wp)
{
    double x = xyz.v[0], y = xyz.v[1], z = xyz.v[2];
    double u1, v1, denom;

    denom = x + 15*y + 3*z;
    if (denom > 0) {
        u1 = (4*x) / denom;
        v1 = (9*y) / denom;
    } else {
        u1 = v1 = 0;
    }
    y /= wp.v[1];
    y = lab_ft_forward(y);
    v3 luv;
    luv.v[0] = 116*y - 16;
    luv.v[1] = 13*(luv.v[0])*(u1 - 4*wp.v[0] / (wp.v[0] + 15*wp.v[1] + 3*wp.v[2]));
    luv.v[2] = 13*(luv.v[0])*(v1 - 9*wp.v[1] / (wp.v[0] + 15*wp.v[1] + 3*wp.v[2]));
    return luv;
}

static double
hue2rgb(double v1,
        double v2,
        double vH)
{
    if (vH < 0) vH += 1;
    if (vH > 1) vH -= 1;
    if (6 * vH < 1) return v1 + (v2 - v1) * 6 * vH;
    if (2 * vH < 1) return v2;
    if (3 * vH < 2) return v1 + (v2 - v1) * (2.0/3.0 - vH) * 6;
    return v1;
}


v3
hsl2rgb(v3 hsl)
{
    double H = hsl.v[0], S = hsl.v[1], L = hsl.v[2];
    double R, G, B;
    if (S == 0) {
        R = G = B = 1;
    } else {
        double v1, v2;
        if (L < 0.5) {
            v2 = L * (1 + S);
        } else {
            v2 = (L + S) - (S * L);
        }
        v1 = 2 * L - v2;
        R = hue2rgb(v1, v2, H + 1.0/3.0);
        G = hue2rgb(v1, v2, H);
        B = hue2rgb(v1, v2, H - 1.0/3.0);
    }
    return (v3){{R,G,B}};
}

v3
rgb2hsl(v3 rgb)
{
    double R = rgb.v[0], G = rgb.v[1], B = rgb.v[2];
    double H, S, L;

    double mn = R;
    if (G < mn) mn = G;
    if (B < mn) mn = B;
    double mx = R;
    if (G > mx) mx = G;
    if (B > mx) mx = B;
    double delta = mx - mn;

    L = (mx + mn) * 0.5;
    if (delta == 0.0) {
        H = S = 0;
    } else {
        if (L < 0.5) S = delta / (mx + mn);
        else         S = delta / (2.0 - mx - mn);

        double delta_R = ((mx - R) / 6.0 + delta * 0.5) / delta;
        double delta_G = ((mx - G) / 6.0 + delta * 0.5) / delta;
        double delta_B = ((mx - B) / 6.0 + delta * 0.5) / delta;

        H = 0;
        if      (R == mx) H = delta_B - delta_G;
        else if (G == mx) H = 1.0/3.0 + delta_R - delta_B;
        else if (B == mx) H = 2.0/3.0 + delta_G - delta_R;

        if (H < 0.0) H += 1.0;
        if (H > 1.0) H -= 1.0;
    }
    return (v3){{ H, S, L }};
}

v3
hsv2rgb(v3 hsv)
{
    return dcp_hsv2rgb((v3){{ hsv.v[0] * 6, hsv.v[1], hsv.v[2] }});
}

v3
rgb2hsv(v3 rgb)
{
    v3 hsv = dcp_rgb2hsv(rgb);
    hsv.v[0] /= 6.0;
    return hsv;
}

double
spec2xyz_weight(const struct observer *o,
                const spectrum_t *illuminant)
{
    return 1.0 / (spec_mulint(o->cmf[1], illuminant, -1, -1) * 0.01);
}

v3
spec2xyz(const spectrum_t *s,
         const struct observer *o)
{
    return spec2xyz_ill(s, o, NULL, 1);
}

v3
spec2xyz_ill(const spectrum_t *s,
             const struct observer *o,
             const spectrum_t *illuminant,
             double weight)
{
    v3 xyz;
    for (int i = 0; i < 3; i++) {
        if (illuminant == NULL) {
            xyz.v[i] = spec_mulint(o->cmf[i], s, -1, -1) * 0.01;
        } else {
            spectrum_t *s_ = spec_multiply(s, illuminant);
            xyz.v[i] = spec_mulint(o->cmf[i], s_, -1, -1) * weight * 0.01;
            free(s_);
        }
    }
    return xyz;
}

v3
xy2xyz_lightest(double x,
                double y,
                const v3 wp)
{
    double X, Y, Z;
    Y = wp.v[1];
    Z = (1-x-y)*Y/y;
    if (Z > wp.v[2]) {
        Z = wp.v[2];
        Y = y*Z/(1-y-x);
    }
    X = x*Y / y;
    if (X > wp.v[0]) {
        X = wp.v[0];
        Y = y*X/x;
        Z = (1-x-y)*Y/y;
    }
    return (v3){{ X, Y, Z }};
}

v3
uv2xyz_lightest(double u,
                double v,
                const v3 wp)
{
    double X, Y, Z;
    Y = wp.v[1];
    Z = Y * (12-3*u-20*v) / (4*v);
    if (Z > wp.v[2]) {
        Z = wp.v[2];
        Y = 4*v*Z/(12-3*u-20*v);
    }
    X = Y*9*u / (4*v);
    if (X > wp.v[0]) {
        X = wp.v[0];
        Y = 4*v*X/(9*u);
        Z = Y * (12-3*u-20*v)/(4*v);
    }
    return (v3){{ X, Y, Z }};
}

m3x3
bradford_map_white(v3 src_wp,
                   v3 dst_wp)
{
    static const m3x3 bfd = {{ {  0.8951,  0.2664, -0.1614 },
                               { -0.7502,  1.7135,  0.0367 },
                               {  0.0389, -0.0685,  1.0296 } }};

    v3 w1 = m3x3_mul_v3(bfd, src_wp);
    v3 w2 = m3x3_mul_v3(bfd, dst_wp);

    for (int i = 0; i < 3; i++) {
        if (w1.v[i] < 0) w1.v[i] = 0.0;
        if (w2.v[i] < 0) w2.v[i] = 0.0;
    }

    m3x3 a = {{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} }};
    for (int i = 0; i < 3; i++) {
        double v = w1.v[i] == 0 ? 10.0 : w2.v[i] / w1.v[i];
        if (v < 0.1) v = 0.1;
        if (v > 10.0) v = 10.0;
        a.v[i][i] = v;
    }

    m3x3 b = m3x3_mul_m3x3(m3x3_mul_m3x3(m3x3_invert(bfd), a), bfd);
    return b;
}

/*
static double
ciecam02_la2fl(const double la)
{
    const double la5 = la * 5;
    const double k = 1.0 / (la5 + 1.0);
    const double k4 = k * k * k * k;
    return (0.2 * k4 * la5) + (0.1 * (1.0 - k4) * (1.0 - k4) * pow(la5, 1.0 / 3.0));
}
*/

// CAM02-UCS
v3
jch2ucs(const v3 jch)
{
    //const double La = 20;
    //const double fl_p025 = pow(ciecam02_la2fl(La), 0.25);
    const double fl_p025 = 0.825404;
    double M = jch.v[1] * fl_p025;

    const double c1 = 0.007;
    const double c2 = 0.0228;
    double Jp = (1.0 + 100.0 * c1) * jch.v[0] / (1.0 + c1 * jch.v[0]);
    double Mp = (1.0 / c2) * log(1.0 + c2 * M);
    double a = Mp * cos(jch.v[2] / 180.0 * M_PI);
    double b = Mp * sin(jch.v[2] / 180.0 * M_PI);
    v3 jab = (v3){{ Jp, a, b }};
    return jab;
}

v3
xyz2ucs(const v3 xyz,
        const v3 wp)
{
    return jch2ucs(xyz2jch(xyz, wp));
}

v3
ucs2xyz(const v3 jab,
        const v3 wp)
{
    return jch2xyz(ucs2jch(jab), wp);
}

v3
ucs2jch(const v3 jab)
{
    if (jab.v[0] == 0 && jab.v[1] == 0 && jab.v[2] == 0) {
        return (v3){{ 0, 0, 0 }};
    }

    const double fl_p025 = 0.825404;
    const double c1 = 0.007;
    const double c2 = 0.0228;

    double Jp = jab.v[0];
    double J = Jp / (1 - (Jp - 100) * c1);
    double Mp = sqrt(jab.v[1] * jab.v[1] + jab.v[2] * jab.v[2]);
    double M = (exp(Mp * c2) - 1.0) / c2;
    double C = M / fl_p025;
    double h = atan2(jab.v[2], jab.v[1]) * 180.0 / M_PI;
    if (h < 0) h += 360.0;
    return (v3){{ J, C, h }};
}

v3
xyz2ucs_s(v3 xyz, void *h02) // must have La == 20
{
    return jch2ucs(xyz2jch_s(xyz, h02));
}

v3
ucs2xyz_s(v3 jab, void *h02)
{
    return jch2xyz_s(ucs2jch(jab), h02);
}

v3
xyz2jch(const v3 xyz,
        const v3 wp)
{
    cmsViewingConditions vc;
    vc.whitePoint.X = wp.v[0] * 100;
    vc.whitePoint.Y = wp.v[1] * 100;
    vc.whitePoint.Z = wp.v[2] * 100;
    vc.Yb = 20;
    vc.La = 20;
    vc.surround = AVG_SURROUND;
    vc.D_value = 1.0;
    cmsHANDLE h = cmsCIECAM02Init(0, &vc);
    cmsCIEXYZ XYZ = { .X = xyz.v[0] * 100, .Y = xyz.v[1] * 100, .Z = xyz.v[2] * 100 };
    cmsJCh JCh;
    cmsCIECAM02Forward(h, &XYZ, &JCh);
    cmsCIECAM02Done(h);
    return (v3){{ JCh.J, JCh.C, JCh.h }};
}

v3
jch2xyz(const v3 jch,
        const v3 wp)
{
    if (jch.v[0] == 0 && jch.v[1] == 0) {
        return (v3){{ 0, 0, 0 }};
    }
    cmsViewingConditions vc;
    vc.whitePoint.X = wp.v[0] * 100;
    vc.whitePoint.Y = wp.v[1] * 100;
    vc.whitePoint.Z = wp.v[2] * 100;
    vc.Yb = 20;
    vc.La = 20;
    vc.surround = AVG_SURROUND;
    vc.D_value = 1.0;
    cmsHANDLE h = cmsCIECAM02Init(0, &vc);
    cmsJCh JCh = { .J = jch.v[0], .C = jch.v[1], .h = jch.v[2] };
    cmsCIEXYZ XYZ;
    cmsCIECAM02Reverse(h, &JCh, &XYZ);
    cmsCIECAM02Done(h);
    v3 xyz = (v3){{ XYZ.X / 100.0, XYZ.Y / 100.0, XYZ.Z / 100.0 }};
    return xyz;
}

v3
xyz2ipt(v3 xyz, v3 wp)
{
    const v3 d65 = {{0.950456, 1.0, 1.088754}};
    if (!v3_isequal(wp, d65)) {
        xyz = chromatic_adaptation_transform(CAT_BRADFORD, xyz, wp, d65);
    }
    const m3x3 xyz2lms = {{
            {  0.4002, 0.7075, -0.0807 },
            { -0.228,  1.1500,  0.0612 },
            {  0.0,    0.0,     0.9184 },
        }};
    v3 lms = m3x3_mul_v3(xyz2lms, xyz);
    double Lp = lms.v[0] >= 0 ? pow(lms.v[0], 0.43) : -pow(-lms.v[0], 0.43);
    double Mp = lms.v[1] >= 0 ? pow(lms.v[1], 0.43) : -pow(-lms.v[1], 0.43);
    double Sp = lms.v[2] >= 0 ? pow(lms.v[2], 0.43) : -pow(-lms.v[2], 0.43);

    const m3x3 lmsp2ipt = {{
            { 0.4000,  0.4000,  0.2000 },
            { 4.4550, -4.8510,  0.3960 },
            { 0.8056,  0.3572, -1.1628 },
        }};
    v3 ipt = m3x3_mul_v3(lmsp2ipt, (v3){{Lp,Mp,Sp}});
    return ipt;
}

v3
ipt2xyz(v3 ipt, v3 wp)
{
    const m3x3 ipt2lmsp = m3x3_invert((m3x3){{
                { 0.4000,  0.4000,  0.2000 },
                    { 4.4550, -4.8510,  0.3960 },
                        { 0.8056,  0.3572, -1.1628 },
                }});
    v3 lmsp = m3x3_mul_v3(ipt2lmsp, ipt);
    double L = lmsp.v[0] >= 0 ? pow(lmsp.v[0], 1/0.43) : -pow(-lmsp.v[0], 1/0.43);
    double M = lmsp.v[1] >= 0 ? pow(lmsp.v[1], 1/0.43) : -pow(-lmsp.v[1], 1/0.43);
    double S = lmsp.v[2] >= 0 ? pow(lmsp.v[2], 1/0.43) : -pow(-lmsp.v[2], 1/0.43);
    const m3x3 lms2xyz = m3x3_invert((m3x3){{
            {  0.4002, 0.7075, -0.0807 },
            { -0.228,  1.1500,  0.0612 },
            {  0.0,    0.0,     0.9184 },
                }});
    v3 xyz = m3x3_mul_v3(lms2xyz, (v3){{L,M,S}});
    const v3 d65 = {{0.950456, 1.0, 1.088754}};
    if (!v3_isequal(wp, d65)) {
        xyz = chromatic_adaptation_transform(CAT_BRADFORD, xyz, d65, wp);
    }
    for (int i = 0; i < 3; i++) {
        // minor precision errors, avoid negative xyz
        if (xyz.v[i] < 0 && xyz.v[i] > -1e-11) xyz.v[i] = 0;
    }
    return xyz;
}

v3
xyz2ipt_polar(v3 xyz, v3 wp)
{
    v3 ipt = xyz2ipt(xyz, wp);
    double lightness = ipt.v[0];
    double chroma = sqrt(ipt.v[1]*ipt.v[1] + ipt.v[2]*ipt.v[2]);
    double hue = atan2(ipt.v[2], ipt.v[1]) * 180.0 / M_PI;
    if (hue < 0) hue += 360.0;
    return (v3){{ lightness, chroma, hue }};
}

v3
ipt2xyz_polar(v3 ipt_polar, v3 wp)
{
    double rad = ipt_polar.v[2] / 180.0 * M_PI;
    double p = ipt_polar.v[1] * cos(rad);
    double t = ipt_polar.v[1] * sin(rad);
    v3 ipt = {{ipt_polar.v[0], p, t }};
    return ipt2xyz(ipt, wp);
}

/*
static double
get_hue_fade(v3 jch,
             double low,
             double high,
             double fade)
{
    double hue = jch.v[2] + 360;
    low += 360;
    high += 360;
    if (hue < low - fade || hue > high + fade) {
        return 0;
    }
    if (hue < low) {
        return 1.0 - scurve((low - hue) / fade);
    }
    if (hue > high) {
        return 1.0 - scurve((hue - high) / fade);
    }
    return 1;
}
*/

v3
chromatic_adaptation_transform(enum ca_transform cat,
                               const v3 xyz,
                               const v3 src_wp,
                               const v3 dst_wp)
{
    if (memcmp(&src_wp, &dst_wp, sizeof(src_wp)) == 0) {
        return xyz;
    }
    switch (cat) {
    case CAT_CAT02: {
        cmsViewingConditions vc1, vc2;
        cmsJCh Out;
        cmsCIEXYZ In;
        cmsHANDLE h1, h2;
        vc1.whitePoint.X = src_wp.v[0] * 100;
        vc1.whitePoint.Y = src_wp.v[1] * 100;
        vc1.whitePoint.Z = src_wp.v[2] * 100;
        vc1.Yb = 20;
        vc1.La = 20;
        vc1.surround = AVG_SURROUND;
        vc1.D_value = 1.0;
        h1 = cmsCIECAM02Init(0, &vc1);
        vc2.whitePoint.X = dst_wp.v[0] * 100;
        vc2.whitePoint.Y = dst_wp.v[1] * 100;
        vc2.whitePoint.Z = dst_wp.v[2] * 100;
        vc2.Yb = 20;
        vc2.La = 20;
        vc2.surround = AVG_SURROUND;
        vc2.D_value = 1.0;
        h2 = cmsCIECAM02Init(0, &vc2);
        In.X = xyz.v[0] * 100;
        In.Y = xyz.v[1] * 100;
        In.Z = xyz.v[2] * 100;
        cmsCIECAM02Forward(h1, &In, &Out);
        cmsCIECAM02Reverse(h2, &Out, &In);
        cmsCIECAM02Done(h1);
        cmsCIECAM02Done(h2);
        In.X /= 100;
        In.Y /= 100;
        In.Z /= 100;
        /*{

            elog("%f %f %f | %f %f %f\n", xyz.v[0],xyz.v[1],xyz.v[2],In.X, In.Y, In.Z);
            v3 out = chromatic_adaptation_transform(CAT_BRADFORD, xyz, src_wp, dst_wp);
            elog("1 %f %f %f | %f %f %f\n", xyz.v[0],xyz.v[1],xyz.v[2],out.v[0],out.v[1],out.v[2]);
            }*/
        return (v3){{ In.X, In.Y, In.Z }};
    }
        /*
          Experimental code, turned out that CAT02 is perfectly fine for StdA.
    case CAT_DCAMPROF: {
        // Here's a modified CAT02 to better simulate the color appearance in tungsten light, that is
        // the conversion from StdA to D50. This CAT is not reversible (in the current implementation).
        v3 cat02_xyz = chromatic_adaptation_transform(CAT_CAT02, xyz, src_wp, dst_wp);
        v3 src_xyY = xyz2xyY(src_wp);
        double src_whiteXY[2] = { src_xyY.v[0], src_xyY.v[1] }, src_temp, src_tint;
        dcp_xy2temp(src_whiteXY, &src_temp, &src_tint);
        double adjval;
        if (src_temp < 2900) {
            adjval = 1;
        } else if (src_temp > 4900) {
            return cat02_xyz;
        } else {
            adjval = 1.0 - (src_temp - 2900) / (4900 - 2900);
        }

        v3 dst_xyY = xyz2xyY(dst_wp);
        double dst_whiteXY[2] = { dst_xyY.v[0], dst_xyY.v[1] }, dst_temp, dst_tint;
        dcp_xy2temp(dst_whiteXY, &dst_temp, &dst_tint);
        if (dst_temp <= src_temp) {
            return cat02_xyz;
        } else if (dst_temp < 4900) {
            adjval *= (dst_temp - src_temp) / (4900 - src_temp);
        }
        v3 cat02_jch = xyz2jch(cat02_xyz, dst_wp);
        bool did_change = false;
        v3 jch = cat02_jch;
        double adj = get_hue_fade(jch, 6, 66, 36);
        if (adj > 0) {
            v3 alt_jch = jch;
            alt_jch.v[1] *= 1.05;
            jch = jch_blend(alt_jch, jch, adj);
            did_change = true;
        }
        adj = get_hue_fade(jch, 18, 32, 36);
        if (adj > 0) {
            v3 alt_jch = jch;
            alt_jch.v[1] *= 1.05;
            alt_jch.v[2] += 4;
            jch = jch_blend(alt_jch, jch, adj);
            did_change = true;
        }
        if (!did_change) {
            return cat02_xyz;
        }
        jch = jch_blend(jch, cat02_jch, adjval);
        v3 alt_xyz = jch2xyz(jch, dst_wp);
        if (!v3_isfinite(alt_xyz)) {
            return cat02_xyz;
        }
        return alt_xyz;
    }
        */
    case CAT_BRADFORD: {
        /*
        cmsCIEXYZ xyz_ = { .X = xyz.v[0], .Y = xyz.v[1], .Z = xyz.v[2] };
        cmsCIEXYZ src_wp_ = { .X = src_wp.v[0], .Y = src_wp.v[1], .Z = src_wp.v[2] };
        cmsCIEXYZ dst_wp_ = { .X = dst_wp.v[0], .Y = dst_wp.v[1], .Z = dst_wp.v[2] };
        cmsCIEXYZ res;
        cmsAdaptToIlluminant(&res, &src_wp_, &dst_wp_, &xyz_);
        return (v3){{ res.X, res.Y, res.Z }};
        */
        return m3x3_mul_v3(bradford_map_white(src_wp, dst_wp), xyz);
    }
    default: abort();
    }
    return xyz;
}

double
ciede2000(v3 xyz1_,
          v3 xyz2_,
          v3 wp_)
{
    cmsCIEXYZ xyz1 = { .X = xyz1_.v[0], .Y = xyz1_.v[1], .Z = xyz1_.v[2] };
    cmsCIEXYZ xyz2 = { .X = xyz2_.v[0], .Y = xyz2_.v[1], .Z = xyz2_.v[2] };
    cmsCIELab lab1, lab2;
    cmsCIEXYZ wp = { .X = wp_.v[0], .Y = wp_.v[1], .Z = wp_.v[2] };
    cmsXYZ2Lab(&wp, &lab1, &xyz1);
    cmsXYZ2Lab(&wp, &lab2, &xyz2);
    return cmsCIE2000DeltaE(&lab1, &lab2, 1, 1, 1);
}

double
ciede2000_lch(v3 xyz1_, v3 xyz2_, v3 wp_,
              double Kl, double Kc, double Kh)
{
    cmsCIEXYZ xyz1 = { .X = xyz1_.v[0], .Y = xyz1_.v[1], .Z = xyz1_.v[2] };
    cmsCIEXYZ xyz2 = { .X = xyz2_.v[0], .Y = xyz2_.v[1], .Z = xyz2_.v[2] };
    cmsCIELab lab1, lab2;
    cmsCIEXYZ wp = { .X = wp_.v[0], .Y = wp_.v[1], .Z = wp_.v[2] };
    cmsXYZ2Lab(&wp, &lab1, &xyz1);
    cmsXYZ2Lab(&wp, &lab2, &xyz2);
    return cmsCIE2000DeltaE(&lab1, &lab2, Kl, Kc, Kh);
}

double
ciede2000_lch_wlc(v3 xyz1, v3 xyz2, v3 wp,
                  double Kl, double Kc, double Kh,
                  double wl, double wc, double wh)
{
    double de = ciede2000_lch(xyz1, xyz2, wp, Kl, Kc, Kh);
    if (wl == 1.0 && wc == 1.0 && wh == 1.0) {
        return de;
    }
    double deL = ciede2000_lch(xyz1, xyz2, wp, 1, 1000000, 1000000);
    double deC = ciede2000_lch(xyz1, xyz2, wp, 1000000, 1, 1000000);
    double deH = ciede2000_lch(xyz1, xyz2, wp, 1000000, 1000000, 1);
    double desum = deL + deC + deH;
    v3 lab1 = xyz2lab(xyz1, wp);
    v3 lab2 = xyz2lab(xyz2, wp);
    double new_de = de;
    if (wl != 1.0) {
        if (wl > 0.0) {
            if (lab2.v[0] > lab1.v[0]) {
                deL = de * deL / desum;
                new_de -= deL;
                new_de += deL * wl;
            }
        } else {
            if (lab2.v[0] < lab1.v[0]) {
                deL = de * deL / desum;
                new_de -= deL;
                new_de += deL * -wl;
            }
        }
    }
    if (wc != 1.0) {
        double c1 = sqrt(lab1.v[1]*lab1.v[1] + lab1.v[2]*lab1.v[2]);
        double c2 = sqrt(lab2.v[1]*lab2.v[1] + lab2.v[2]*lab2.v[2]);
        if (wc > 0.0) {
            if (c2 > c1) {
                deC = de * deC / desum;
                new_de -= deC;
                new_de += deC * wc;
            }
        } else {
            if (c2 < c1) {
                deC = de * deC / desum;
                new_de -= deC;
                new_de += deC * -wc;
            }
        }
    }
    if (wh != 1.0) {
        deH = de * deH / desum;
        new_de -= deH;
        new_de += deH * wh;
    }
//    elog("%f %f (%f %f | %f %f)\n", de, new_de, lab1.v[0], c1, lab2.v[0], c2);
    if (new_de < 0) {
        new_de = 0;
    }
    return new_de;
}

double
euclidean_dist(v3 v1,
               v3 v2)
{
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        double diff = v2.v[i] - v1.v[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
