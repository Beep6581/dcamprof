/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */

// utility functions derived from DNG reference code
#include <stdbool.h>

#include <colmath.h>
#include <elog.h>
#include <dngref.h>

#define UNUSED(x) (void)(x)

static const double d50_xy[2] = { 0.3457, 0.3585 };

v3
dcp_lookup_rgb(const v3 *hsm,
               const uint32_t hsmdims[3],
               const bool srgb_gamma,
               const v3 rgb,
               bool *hsv_clip)
{
    v3 hsv = dcp_rgb2hsv(rgb);

    v3 hsv1;

    { // DCP ref code
        uint32_t hueDivisions = hsmdims[0];
        uint32_t satDivisions = hsmdims[1];
        uint32_t valDivisions = hsmdims[2];

        double hScale = (hueDivisions < 2) ? 0.0 : (hueDivisions * (1.0 / 6.0));
        double sScale = (double) (satDivisions - 1);
        double vScale = (double) (valDivisions - 1);
        int32_t maxHueIndex0 = hueDivisions - 1;
        int32_t maxSatIndex0 = satDivisions - 2;
        int32_t maxValIndex0 = valDivisions - 2;
        const v3 *tableBase = hsm;
        int32_t hueStep = satDivisions;
        int32_t valStep = hueDivisions * hueStep;

        double h, s, v;

        h = hsv.v[0];
        s = hsv.v[1];
        v = hsv.v[2];

        double vEncoded = v;

        double hueShift;
        double satScale;
        double valScale;
        if (valDivisions < 2) {// Optimize most common case of "2.5D" table.
            double hScaled = h * hScale;
            double sScaled = s * sScale;
            int32_t hIndex0 = (int32_t) hScaled;
            int32_t sIndex0 = (int32_t) sScaled;
            if (sIndex0 > maxSatIndex0) sIndex0 = maxSatIndex0;
            int32_t hIndex1 = hIndex0 + 1;
            if (hIndex0 >= maxHueIndex0) {
                hIndex0 = maxHueIndex0;
                hIndex1 = 0;
            }
            double hFract1 = hScaled - (double) hIndex0;
            double sFract1 = sScaled - (double) sIndex0;
            double hFract0 = 1.0f - hFract1;
            double sFract0 = 1.0f - sFract1;
            const v3 *entry00 = tableBase + hIndex0 * hueStep + sIndex0;
            const v3 *entry01 = entry00 + (hIndex1 - hIndex0) * hueStep;
            double hueShift0 = hFract0 * entry00->v[0] + hFract1 * entry01->v[0];
            double satScale0 = hFract0 * entry00->v[1] + hFract1 * entry01->v[1];
            double valScale0 = hFract0 * entry00->v[2] + hFract1 * entry01->v[2];

            entry00++;
            entry01++;

            double hueShift1 = hFract0 * entry00->v[0] + hFract1 * entry01->v[0];
            double satScale1 = hFract0 * entry00->v[1] + hFract1 * entry01->v[1];
            double valScale1 = hFract0 * entry00->v[2] + hFract1 * entry01->v[2];
            hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
            satScale = sFract0 * satScale0 + sFract1 * satScale1;
            valScale = sFract0 * valScale0 + sFract1 * valScale1;
        } else {
            double hScaled = h		  * hScale;
            double sScaled = s		  * sScale;
            if (srgb_gamma) {
                double val = v;
                if (val < 0.0) val = 0.0;
                if (val > 1.0) val = 1.0;
                vEncoded = srgb_gamma_forward(val);
            }
            double vScaled = vEncoded * vScale;
            int32_t hIndex0 = (int32_t) hScaled;
            int32_t sIndex0 = (int32_t) sScaled;
            int32_t vIndex0 = (int32_t) vScaled;
            if (sIndex0 > maxSatIndex0) sIndex0 = maxSatIndex0;
            if (vIndex0 > maxValIndex0) vIndex0 = maxValIndex0;
            int32_t hIndex1 = hIndex0 + 1;
            if (hIndex0 >= maxHueIndex0) {
                hIndex0 = maxHueIndex0;
                hIndex1 = 0;
            }
            double hFract1 = hScaled - (double) hIndex0;
            double sFract1 = sScaled - (double) sIndex0;
            double vFract1 = vScaled - (double) vIndex0;
            double hFract0 = 1.0 - hFract1;
            double sFract0 = 1.0 - sFract1;
            double vFract0 = 1.0 - vFract1;
            const v3 *entry00 = tableBase + vIndex0 * valStep +  hIndex0 * hueStep + sIndex0;
            const v3 *entry01 = entry00 + (hIndex1 - hIndex0) * hueStep;
            const v3 *entry10 = entry00 + valStep;
            const v3 *entry11 = entry01 + valStep;
            double hueShift0 = vFract0 * (hFract0 * entry00->v[0] + hFract1 * entry01->v[0]) + vFract1 * (hFract0 * entry10->v[0] + hFract1 * entry11->v[0]);
            double satScale0 = vFract0 * (hFract0 * entry00->v[1] + hFract1 * entry01->v[1]) + vFract1 * (hFract0 * entry10->v[1] + hFract1 * entry11->v[1]);
            double valScale0 = vFract0 * (hFract0 * entry00->v[2] + hFract1 * entry01->v[2]) + vFract1 * (hFract0 * entry10->v[2] + hFract1 * entry11->v[2]);
            entry00++;
            entry01++;
            entry10++;
            entry11++;

            double hueShift1 = vFract0 * (hFract0 * entry00->v[0] + hFract1 * entry01->v[0]) + vFract1 * (hFract0 * entry10->v[0] + hFract1 * entry11->v[0]);
            double satScale1 = vFract0 * (hFract0 * entry00->v[1] + hFract1 * entry01->v[1]) + vFract1 * (hFract0 * entry10->v[1] + hFract1 * entry11->v[1]);
            double valScale1 = vFract0 * (hFract0 * entry00->v[2] + hFract1 * entry01->v[2]) + vFract1 * (hFract0 * entry10->v[2] + hFract1 * entry11->v[2]);
            hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
            satScale = sFract0 * satScale0 + sFract1 * satScale1;
            valScale = sFract0 * valScale0 + sFract1 * valScale1;
        }

        bool did_clip = false;
        hueShift *= (6.0 / 360.0);	// Convert to internal hue range.
        h += hueShift;
        if (h > 6.0) h -= 6.0;
        if (h < 0) h += 6.0;
        s = s * satScale;
        if (s > 1.0) {
            s = 1.0;
            did_clip = true;
        }

        vEncoded = vEncoded * valScale;
        if (vEncoded < 0) {
            vEncoded = 0;
            did_clip = true;
        }
        if (vEncoded > 1) {
            vEncoded = 1;
            did_clip = true;
        }

        if (srgb_gamma) {
            v = srgb_gamma_inverse(vEncoded);
        } else {
            v = vEncoded;
        }

        hsv1 = (v3){{ h, s, v }};
        if (did_clip) {
            //elog("HSV clipping!\n");
            //v3_print("hsv", hsv);
        }
        if (hsv_clip) *hsv_clip = did_clip;
    }
    v3 rgb1 = dcp_hsv2rgb(hsv1);

    return rgb1;
}

v3
dcp_lookup(const v3 *hsm,
           const uint32_t hsmdims[3],
           const bool srgb_gamma,
           const v3 cxyz,
           double bse_scale,
           const double tc[],
           int tc_len,
           bool *xyz_clip,
           bool *hsv_clip)
{
    /*
      Note about a full DNG pipeline:
        1) convert to rgb via matrix
        2) apply huesatmap
        3) apply exposure curve including exposure adjustment and black subtraction
            - here we assume that there is no black subtraction and exposure is zero or pre-applied, we only apply
            - baseline exposure offset
        4) apply look table
        5) apply tone curve
     */
    v3 rgb = m3x3_mul_v3(dcp_xyz_D50_to_prophoto_rgb, cxyz);
    if (bse_scale != 0 && bse_scale != 1) {
        rgb = k_mul_v3(bse_scale, rgb);
    }
    bool did_clip = false;
    for (int i = 0; i < 3; i++) {
        if (rgb.v[i] < 0) {
            rgb.v[i] = 0;
            did_clip = true;
        } else if (rgb.v[i] > 1) {
            rgb.v[i] = 1;
            did_clip = true;
        }
    }
    if (xyz_clip) *xyz_clip = did_clip;
    if (did_clip) {
        //elog("Prophoto clipping!\n");
        //v3_print("rgb", rgb);
    }

    // output
    v3 rgb1 = dcp_lookup_rgb(hsm, hsmdims, srgb_gamma, rgb, hsv_clip);
    if (tc != NULL) {
        rgb1 = dcp_apply_tonecurve(rgb1, tc, tc_len);
    }
    const m3x3 rgb2xyz = m3x3_invert(dcp_xyz_D50_to_prophoto_rgb);
    v3 xyz1 = m3x3_mul_v3(rgb2xyz, rgb1);
    return xyz1;
}

// rgb range 0 - 1, hsv range 0 - 1 for s and v, 0 - 6 for h
v3
dcp_rgb2hsv(v3 rgb)
{
    double h, s, v;
    double largest = rgb.v[0];
    double smallest = rgb.v[0];
    for (int i = 1; i < 3; i++) {
        if (rgb.v[i] > largest) largest = rgb.v[i];
        if (rgb.v[i] < smallest) smallest = rgb.v[i];
    }
    v = largest;
    double span = largest - smallest;
    if (span > 0.0) {
        if (rgb.v[0] == v) {
            h = (rgb.v[1] - rgb.v[2]) / span;
            if (h < 0.0) {
                h += 6.0;
            }
        } else if (rgb.v[1] == v) {
            h = 2.0 + (rgb.v[2] - rgb.v[0]) / span;
        } else {
            h = 4.0 + (rgb.v[0] - rgb.v[1]) / span;
        }
        s = span / v;
    } else {
        h = s = 0.0;
    }
    return (v3){{ h, s, v }};
}

v3
dcp_hsv2rgb(const v3 hsv)
{
    double h = hsv.v[0], s = hsv.v[1], v = hsv.v[2], r, g, b;
    if (s > 0.0f) {
        if (h < 0.0f) {
            h += 6.0f;
        }
        if (h >= 6.0f) {
            h -= 6.0f;
        }
        int  i = (int)h;
        double f = h - (double)i;
        double p = v * (1.0f - s);
#define q       (v * (1.0f - s * f))
#define t       (v * (1.0f - s * (1.0f - f)))
        switch (i) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        case 5: r = v; g = p; b = q; break;
        default: r = g = b = 0; // never reached
        }
#undef q
#undef t
    } else {
        r = v;
        g = v;
        b = v;
    }
    return (v3){{ r, g, b }};
}

struct ruvt {
    double r;
    double u;
    double v;
    double t;
};

static const double kTintScale = -3000.0;
static const struct ruvt kTempTable [] =
	{
	{   0, 0.18006, 0.26352, -0.24341 },
	{  10, 0.18066, 0.26589, -0.25479 },
	{  20, 0.18133, 0.26846, -0.26876 },
	{  30, 0.18208, 0.27119, -0.28539 },
	{  40, 0.18293, 0.27407, -0.30470 },
	{  50, 0.18388, 0.27709, -0.32675 },
	{  60, 0.18494, 0.28021, -0.35156 },
	{  70, 0.18611, 0.28342, -0.37915 },
	{  80, 0.18740, 0.28668, -0.40955 },
	{  90, 0.18880, 0.28997, -0.44278 },
	{ 100, 0.19032, 0.29326, -0.47888 },
	{ 125, 0.19462, 0.30141, -0.58204 },
	{ 150, 0.19962, 0.30921, -0.70471 },
	{ 175, 0.20525, 0.31647, -0.84901 },
	{ 200, 0.21142, 0.32312, -1.0182 },
	{ 225, 0.21807, 0.32909, -1.2168 },
	{ 250, 0.22511, 0.33439, -1.4512 },
	{ 275, 0.23247, 0.33904, -1.7298 },
	{ 300, 0.24010, 0.34308, -2.0637 },
	{ 325, 0.24702, 0.34655, -2.4681 },
	{ 350, 0.25591, 0.34951, -2.9641 },
	{ 375, 0.26400, 0.35200, -3.5814 },
	{ 400, 0.27218, 0.35407, -4.3633 },
	{ 425, 0.28039, 0.35577, -5.3762 },
	{ 450, 0.28863, 0.35714, -6.7262 },
	{ 475, 0.29685, 0.35823, -8.5955 },
	{ 500, 0.30505, 0.35907, -11.324 },
	{ 525, 0.31320, 0.35968, -15.628 },
	{ 550, 0.32129, 0.36011, -23.325 },
	{ 575, 0.32931, 0.36038, -40.770 },
	{ 600, 0.33724, 0.36051, -116.45 }
	};

void
dcp_xy2temp(const double whiteXY[2],
            double *temp,
            double *tint)
{
    double fTemperature = 0;
    double fTint = 0;

    // Convert to uv space.
    double u = 2.0 * whiteXY[0] / (1.5 - whiteXY[0] + 6.0 * whiteXY[1]);
    double v = 3.0 * whiteXY[1] / (1.5 - whiteXY[0] + 6.0 * whiteXY[1]);

    // Search for line pair coordinate is between.
    double last_dt = 0.0;
    double last_dv = 0.0;
    double last_du = 0.0;

    for (uint32_t index = 1; index <= 30; index++) {
        // Convert slope to delta-u and delta-v, with length 1.
        double du = 1.0;
        double dv = kTempTable [index] . t;
        double len = sqrt (1.0 + dv * dv);
        du /= len;
        dv /= len;

        // Find delta from black body point to test coordinate.
        double uu = u - kTempTable [index] . u;
        double vv = v - kTempTable [index] . v;

        // Find distance above or below line.
        double dt = - uu * dv + vv * du;

        // If below line, we have found line pair.
        if (dt <= 0.0 || index == 30) {
            // Find fractional weight of two lines.
            if (dt > 0.0)
                dt = 0.0;
            dt = -dt;
            double f;
            if (index == 1)
            {
                f = 0.0;
            }
            else
            {
                f = dt / (last_dt + dt);
            }

            // Interpolate the temperature.
            fTemperature = 1.0E6 / (kTempTable [index - 1] . r * f +
                                    kTempTable [index    ] . r * (1.0 - f));

            // Find delta from black body point to test coordinate.
            uu = u - (kTempTable [index - 1] . u * f +
                      kTempTable [index    ] . u * (1.0 - f));
            vv = v - (kTempTable [index - 1] . v * f +
                      kTempTable [index    ] . v * (1.0 - f));
            // Interpolate vectors along slope.
            du = du * (1.0 - f) + last_du * f;
            dv = dv * (1.0 - f) + last_dv * f;
            len = sqrt (du * du + dv * dv);
            du /= len;
            dv /= len;

            // Find distance along slope.
            fTint = (uu * du + vv * dv) * kTintScale;
            break;
        }
        // Try next line pair.
        last_dt = dt;
        last_du = du;
        last_dv = dv;
    }
    if (temp != NULL)
        *temp = fTemperature;
    if (tint != NULL)
        *tint = fTint;
}

void
dcp_temp2xy(double temp,
            double tint,
            double xy[])
{
    const double fTemperature = temp;
    const double fTint = tint;

    // Find inverse temperature to use as index.
    double r = 1.0E6 / fTemperature;

    // Convert tint to offset is uv space.
    double offset = fTint * (1.0 / kTintScale);

    // Search for line pair containing coordinate.
    for (uint32_t index = 0; index <= 29; index++) {
        if (r < kTempTable [index + 1] . r || index == 29) {
            // Find relative weight of first line.
            double f = (kTempTable [index + 1] . r - r) / (kTempTable [index + 1] . r - kTempTable [index] . r);
            // Interpolate the black body coordinates.
            double u = kTempTable [index    ] . u * f + kTempTable [index + 1] . u * (1.0 - f);
            double v = kTempTable [index    ] . v * f + kTempTable [index + 1] . v * (1.0 - f);

            // Find vectors along slope for each line.
            double uu1 = 1.0;
            double vv1 = kTempTable [index] . t;
            double uu2 = 1.0;
            double vv2 = kTempTable [index + 1] . t;
            double len1 = sqrt (1.0 + vv1 * vv1);
            double len2 = sqrt (1.0 + vv2 * vv2);
            uu1 /= len1;
            vv1 /= len1;
            uu2 /= len2;
            vv2 /= len2;

            // Find vector from black body point.
            double uu3 = uu1 * f + uu2 * (1.0 - f);
            double vv3 = vv1 * f + vv2 * (1.0 - f);
            double len3 = sqrt (uu3 * uu3 + vv3 * vv3);
            uu3 /= len3;
            vv3 /= len3;

            // Adjust coordinate along this vector.
            u += uu3 * offset;
            v += vv3 * offset;

            // Convert to xy coordinates.
            xy[0] = 1.5 * u / (u - 4.0 * v + 2.0);
            xy[1] =       v / (u - 4.0 * v + 2.0);
            break;
        }
    }
}


static double
dngref_get_mix(const struct profio_dcp *dcp,
               const double whiteXY[2])
{
    const double fTemperature1 = exif_lightsource_temp(dcp->ill[0]);
    const double fTemperature2 = dcp->ill_count == 2 ? exif_lightsource_temp(dcp->ill[1]) : 0;
    // Convert to temperature/offset space.
    double wbtemp;
    dcp_xy2temp(whiteXY, &wbtemp, NULL);

    // Find fraction to weight the first calibration.
    double g;
    if (wbtemp <= fTemperature1 || dcp->ill_count == 1) {
        g = 1.0;
    } else if (wbtemp >= fTemperature2) {
        g = 0.0;
    } else {
        double invT = 1.0 / wbtemp;
        g = (invT                  - (1.0 / fTemperature2)) /
            ((1.0 / fTemperature1) - (1.0 / fTemperature2));
    }
    return g;
}

m3x3
dcp_find_xyz2cam(const struct profio_dcp *dcp,
                 const double white_xy[2],
                 m3x3 *fm,
                 m3x3 *rm,
                 const m3x3 *cc1,
                 const m3x3 *cc2,
                 const char cc_sig[],
                 m3x3 *cc)
{
    UNUSED(rm); // reduction matrix not implemented yet

    const double g = dngref_get_mix(dcp, white_xy);

    // Color Matrix
    m3x3 cm;
    if (g >= 1.0)
        cm = dcp->cm[0];
    else if (g <= 0.0)
        cm = dcp->cm[1];
    else {
        cm = m3x3_add_m3x3(k_mul_m3x3(g, dcp->cm[0]), k_mul_m3x3(1.0 - g, dcp->cm[1]));
    }

    // Forward Matrix
    if (fm) {
        bool has1 = dcp->has.fm;
        bool has2 = has1 && dcp->ill_count == 2;
        if (has1 && has2) {
            if (g >= 1.0)
                *fm = dcp->fm[0];
            else if (g <= 0.0)
                *fm = dcp->fm[1];
            else
                *fm = m3x3_add_m3x3(k_mul_m3x3(g, dcp->fm[0]), k_mul_m3x3(1.0 - g, dcp->fm[1]));
        } else if (has1) {
            *fm = dcp->fm[0];
        } else if (has2) {
            *fm = dcp->fm[1];
        } else {
            *fm = (m3x3){ .v = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } };
        }
    }

    // CameraCalibration Matrix
    if (cc) {
        *cc = (m3x3){ .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
        if (dcp->calibration_signature != NULL && cc_sig != NULL && strcmp(dcp->calibration_signature, cc_sig) == 0) {
            const m3x3 cc1_ = cc1 != NULL ? *cc1 : (m3x3){ .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
            const m3x3 cc2_ = cc2 != NULL ? *cc2 : (m3x3){ .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
            if (g >= 1.0) {
                *cc = cc1_;
            } else if (g <= 0.0) {
                *cc = cc2_;
            } else {
                *cc = m3x3_add_m3x3(k_mul_m3x3(g, cc1_), k_mul_m3x3(1.0 - g, cc2_));
            }
        }
    }

    /*
    // Interpolate reduction matrix, if any.
    if (reductionMatrix) {
        bool has1 = fReductionMatrix1.NotEmpty ();
        bool has2 = fReductionMatrix2.NotEmpty ();
        if (has1 && has2) {
            if (g >= 1.0)
                *reductionMatrix = fReductionMatrix1;
            else if (g <= 0.0)
                *reductionMatrix = fReductionMatrix2;
            else
                *reductionMatrix = (g      ) * fReductionMatrix1 +
                    (1.0 - g) * fReductionMatrix2;
        } else if (has1) {
            *reductionMatrix = fReductionMatrix1;
        } else if (has2) {
            *reductionMatrix = fReductionMatrix2;
        } else {
            reductionMatrix->Clear ();
        }
    }
    */

    return cm;
}

static void
dngref_xyz2xy(v3 xyz, double xy[2])
{
    const double total = xyz.v[0] + xyz.v[1] + xyz.v[2];
    if (total > 0.0) {
        xy[0] = xyz.v[0] / total;
        xy[1] = xyz.v[1] / total;
    } else {
        xy[0] = 0.3457;
        xy[1] = 0.3585;
    }
}

static v3
dngref_xy2xyz(const double xy[2])
{
    double temp[2] = { xy[0], xy[1] };
    // Restrict xy coord to someplace inside the range of real xy coordinates.
    // This prevents math from doing strange things when users specify
    // extreme temperature/tint coordinates.
    for (int i = 0; i < 2; i++) {
        if (temp[i] < 0.000001) temp[i] = 0.000001;
        if (temp[i] > 0.999999) temp[i] = 0.999999;
    }
    if (temp[0] + temp[1] > 0.999999) {
        double scale = 0.999999 / (temp[0] + temp[1]);
        temp[0] *= scale;
        temp[1] *= scale;
    }
    v3 xyz;
    xyz.v[0] = temp[0] / temp[1];
    xyz.v[1] = 1.0;
    xyz.v[2] = (1.0 - temp[0] - temp[1]) / temp[1];
    return xyz;
}

void
dcp_neutral2xy(const struct profio_dcp *dcp,
               const m3x3 *cc1,
               const m3x3 *cc2,
               const char cc_sig[],
               const v3 *ab_vector,
               const v3 neutral,
               double xy[2])
{
    const int max_passes = 30;
    double last_xy[2] = { 0.3457, 0.3585 }; // D50
    m3x3 ab = (ab_vector != NULL) ? v3_asdiagonal(*ab_vector) : (m3x3){ .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
    for (int pass = 0; pass < max_passes; pass++) {
        // from DNG ref code, but adding in use of CC and AB, as in the DNG spec 1.4
        m3x3 cc;
        m3x3 cm = dcp_find_xyz2cam(dcp, last_xy, NULL, NULL, cc1, cc2, cc_sig, &cc);
        m3x3 xyz2cam = m3x3_mul_m3x3(m3x3_mul_m3x3(ab, cc), cm);

        double next_xy[2];
        v3 next_xyz = m3x3_mul_v3(m3x3_invert(xyz2cam), neutral);
        dngref_xyz2xy(next_xyz, next_xy);

        if (fabs(next_xy[0] - last_xy[0]) + fabs(next_xy[1] - last_xy[1]) < 0.0000001) {
            xy[0] = next_xy[0];
            xy[1] = next_xy[1];
            return;
        }
        if (pass == max_passes - 1) {
            next_xy[0] = (last_xy[0] + next_xy[0]) * 0.5;
            next_xy[1] = (last_xy[1] + next_xy[1]) * 0.5;
        }
        last_xy[0] = next_xy[0];
        last_xy[1] = next_xy[1];
    }
    xy[0] = last_xy[0];
    xy[1] = last_xy[1];
}

m3x3
dcp_make_cam2xyz_d50(const struct profio_dcp *dcp,
                     const m3x3 *cc1,
                     const m3x3 *cc2,
                     const char cc_sig[],
                     const v3 *ab_vector,
                     const v3 wb)
{
    m3x3 cm, fm, cc;
    double xy[2];
    dcp_neutral2xy(dcp, cc1, cc2, cc_sig, ab_vector, wb, xy);
    cm = dcp_find_xyz2cam(dcp, xy, &fm, NULL, cc1, cc2, cc_sig, &cc);
    m3x3 ab = (ab_vector != NULL) ? v3_asdiagonal(*ab_vector) : (m3x3){ .v = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

    m3x3 cam2xyz_D50;
    if (!dcp->has.fm) {
        // doing according to DNG spec 1.4, include CC and AB, which is not used in ref code!

        // XYZtoCamera = AB * CC * CM
        m3x3 xyz2cam = m3x3_mul_m3x3(m3x3_mul_m3x3(ab, cc), cm);
        // First, invert the XYZtoCamera matrix. CameraToXYZ = Inverse (XYZtoCamera)
        // FIXME: add case with reduction matrix!
        m3x3 cam2xyz = m3x3_invert(xyz2cam);
        // CA [...] is a chromatic adaptation matrix that maps from the white balance xy value to
        // the D50 white point. The recommended method for computing this chromatic adaptation
        // matrix is to use the linear Bradford algorithm.
        v3 wp = dngref_xy2xyz(xy);
        m3x3 ca = bradford_map_white(wp, dcp_d50());
        // CameraToXYZ_D50 = CA * CameraToXYZ
        cam2xyz_D50 = m3x3_mul_m3x3(ca, cam2xyz);
    } else {
        // ReferenceNeutral = Inverse (AB * CC) * CameraNeutral
        m3x3 inv_abcc = m3x3_invert(m3x3_mul_m3x3(ab, cc));
        v3 ref_wb = m3x3_mul_v3(inv_abcc, wb);
        // D = Invert (AsDiagonalMatrix (ReferenceNeutral))
        m3x3 d = m3x3_invert(v3_asdiagonal(ref_wb));
        // CameraToXYZ_D50 = FM * D * Inverse (AB * CC)
        cam2xyz_D50 = m3x3_mul_m3x3(m3x3_mul_m3x3(fm, d), inv_abcc);

        // Note if cc is a diagonal matrix, which it almost always is, the above code is equivalent to
        //   cam2xyz_D50 = m3x3_mul_m3x3(fm, m3x3_invert(v3_asdiagonal(wb))));
    }
    return cam2xyz_D50;
}

v3 *
dcp_make_lut(const struct profio_dcp *dcp,
             const m3x3 *cc1,
             const m3x3 *cc2,
             const char cc_sig[],
             const v3 *ab,
             const v3 wb)
{
    if (dcp == NULL || dcp->huesatmap[0] == NULL) {
        return NULL;
    }
    size_t count = dcp->hsmdims[0]*dcp->hsmdims[1]*dcp->hsmdims[2];
    v3 *hsm = malloc(count * sizeof(*hsm));
    if (dcp->ill_count == 1 || dcp->huesatmap[1] == NULL) {
        memcpy(hsm, dcp->huesatmap[0], count * sizeof(*hsm));
        return hsm;
    }
    double xy[2];
    dcp_neutral2xy(dcp, cc1, cc2, cc_sig, ab, wb, xy);
    double mix = dngref_get_mix(dcp, xy);

    if (mix >= 1.0) {
        memcpy(hsm, dcp->huesatmap[0], count * sizeof(*hsm));
        return hsm;
    } else if (mix <= 0.0) {
        memcpy(hsm, dcp->huesatmap[1], count * sizeof(*hsm));
        return hsm;
    }

    // Interpolate between the tables.
    double w1 = mix;
    double w2 = 1.0 - mix;
    for (size_t i = 0; i < count; i++) {
        for (int k = 0; k < 3; k++) {
            hsm[i].v[k] = w1 * dcp->huesatmap[0][i].v[k] + w2 * dcp->huesatmap[1][i].v[k];
        }
    }
    return hsm;
}

v3
dcp_d50(void)
{
    return xy2xyz_lightest(d50_xy[0], d50_xy[1], (v3){{1,1,1}});
}

double *
dcp_acr_tonecurve(int *len)
{
    static const float curve[] = {
        0.00000f, 0.00078f, 0.00160f, 0.00242f,
        0.00314f, 0.00385f, 0.00460f, 0.00539f,
        0.00623f, 0.00712f, 0.00806f, 0.00906f,
        0.01012f, 0.01122f, 0.01238f, 0.01359f,
        0.01485f, 0.01616f, 0.01751f, 0.01890f,
        0.02033f, 0.02180f, 0.02331f, 0.02485f,
        0.02643f, 0.02804f, 0.02967f, 0.03134f,
        0.03303f, 0.03475f, 0.03648f, 0.03824f,
        0.04002f, 0.04181f, 0.04362f, 0.04545f,
        0.04730f, 0.04916f, 0.05103f, 0.05292f,
        0.05483f, 0.05675f, 0.05868f, 0.06063f,
        0.06259f, 0.06457f, 0.06655f, 0.06856f,
        0.07057f, 0.07259f, 0.07463f, 0.07668f,
        0.07874f, 0.08081f, 0.08290f, 0.08499f,
        0.08710f, 0.08921f, 0.09134f, 0.09348f,
        0.09563f, 0.09779f, 0.09996f, 0.10214f,
        0.10433f, 0.10652f, 0.10873f, 0.11095f,
        0.11318f, 0.11541f, 0.11766f, 0.11991f,
        0.12218f, 0.12445f, 0.12673f, 0.12902f,
        0.13132f, 0.13363f, 0.13595f, 0.13827f,
        0.14061f, 0.14295f, 0.14530f, 0.14765f,
        0.15002f, 0.15239f, 0.15477f, 0.15716f,
        0.15956f, 0.16197f, 0.16438f, 0.16680f,
        0.16923f, 0.17166f, 0.17410f, 0.17655f,
        0.17901f, 0.18148f, 0.18395f, 0.18643f,
        0.18891f, 0.19141f, 0.19391f, 0.19641f,
        0.19893f, 0.20145f, 0.20398f, 0.20651f,
        0.20905f, 0.21160f, 0.21416f, 0.21672f,
        0.21929f, 0.22185f, 0.22440f, 0.22696f,
        0.22950f, 0.23204f, 0.23458f, 0.23711f,
        0.23963f, 0.24215f, 0.24466f, 0.24717f,
        0.24967f, 0.25216f, 0.25465f, 0.25713f,
        0.25961f, 0.26208f, 0.26454f, 0.26700f,
        0.26945f, 0.27189f, 0.27433f, 0.27676f,
        0.27918f, 0.28160f, 0.28401f, 0.28641f,
        0.28881f, 0.29120f, 0.29358f, 0.29596f,
        0.29833f, 0.30069f, 0.30305f, 0.30540f,
        0.30774f, 0.31008f, 0.31241f, 0.31473f,
        0.31704f, 0.31935f, 0.32165f, 0.32395f,
        0.32623f, 0.32851f, 0.33079f, 0.33305f,
        0.33531f, 0.33756f, 0.33981f, 0.34205f,
        0.34428f, 0.34650f, 0.34872f, 0.35093f,
        0.35313f, 0.35532f, 0.35751f, 0.35969f,
        0.36187f, 0.36404f, 0.36620f, 0.36835f,
        0.37050f, 0.37264f, 0.37477f, 0.37689f,
        0.37901f, 0.38112f, 0.38323f, 0.38533f,
        0.38742f, 0.38950f, 0.39158f, 0.39365f,
        0.39571f, 0.39777f, 0.39982f, 0.40186f,
        0.40389f, 0.40592f, 0.40794f, 0.40996f,
        0.41197f, 0.41397f, 0.41596f, 0.41795f,
        0.41993f, 0.42191f, 0.42388f, 0.42584f,
        0.42779f, 0.42974f, 0.43168f, 0.43362f,
        0.43554f, 0.43747f, 0.43938f, 0.44129f,
        0.44319f, 0.44509f, 0.44698f, 0.44886f,
        0.45073f, 0.45260f, 0.45447f, 0.45632f,
        0.45817f, 0.46002f, 0.46186f, 0.46369f,
        0.46551f, 0.46733f, 0.46914f, 0.47095f,
        0.47275f, 0.47454f, 0.47633f, 0.47811f,
        0.47989f, 0.48166f, 0.48342f, 0.48518f,
        0.48693f, 0.48867f, 0.49041f, 0.49214f,
        0.49387f, 0.49559f, 0.49730f, 0.49901f,
        0.50072f, 0.50241f, 0.50410f, 0.50579f,
        0.50747f, 0.50914f, 0.51081f, 0.51247f,
        0.51413f, 0.51578f, 0.51742f, 0.51906f,
        0.52069f, 0.52232f, 0.52394f, 0.52556f,
        0.52717f, 0.52878f, 0.53038f, 0.53197f,
        0.53356f, 0.53514f, 0.53672f, 0.53829f,
        0.53986f, 0.54142f, 0.54297f, 0.54452f,
        0.54607f, 0.54761f, 0.54914f, 0.55067f,
        0.55220f, 0.55371f, 0.55523f, 0.55673f,
        0.55824f, 0.55973f, 0.56123f, 0.56271f,
        0.56420f, 0.56567f, 0.56715f, 0.56861f,
        0.57007f, 0.57153f, 0.57298f, 0.57443f,
        0.57587f, 0.57731f, 0.57874f, 0.58017f,
        0.58159f, 0.58301f, 0.58443f, 0.58583f,
        0.58724f, 0.58864f, 0.59003f, 0.59142f,
        0.59281f, 0.59419f, 0.59556f, 0.59694f,
        0.59830f, 0.59966f, 0.60102f, 0.60238f,
        0.60373f, 0.60507f, 0.60641f, 0.60775f,
        0.60908f, 0.61040f, 0.61173f, 0.61305f,
        0.61436f, 0.61567f, 0.61698f, 0.61828f,
        0.61957f, 0.62087f, 0.62216f, 0.62344f,
        0.62472f, 0.62600f, 0.62727f, 0.62854f,
        0.62980f, 0.63106f, 0.63232f, 0.63357f,
        0.63482f, 0.63606f, 0.63730f, 0.63854f,
        0.63977f, 0.64100f, 0.64222f, 0.64344f,
        0.64466f, 0.64587f, 0.64708f, 0.64829f,
        0.64949f, 0.65069f, 0.65188f, 0.65307f,
        0.65426f, 0.65544f, 0.65662f, 0.65779f,
        0.65897f, 0.66013f, 0.66130f, 0.66246f,
        0.66362f, 0.66477f, 0.66592f, 0.66707f,
        0.66821f, 0.66935f, 0.67048f, 0.67162f,
        0.67275f, 0.67387f, 0.67499f, 0.67611f,
        0.67723f, 0.67834f, 0.67945f, 0.68055f,
        0.68165f, 0.68275f, 0.68385f, 0.68494f,
        0.68603f, 0.68711f, 0.68819f, 0.68927f,
        0.69035f, 0.69142f, 0.69249f, 0.69355f,
        0.69461f, 0.69567f, 0.69673f, 0.69778f,
        0.69883f, 0.69988f, 0.70092f, 0.70196f,
        0.70300f, 0.70403f, 0.70506f, 0.70609f,
        0.70711f, 0.70813f, 0.70915f, 0.71017f,
        0.71118f, 0.71219f, 0.71319f, 0.71420f,
        0.71520f, 0.71620f, 0.71719f, 0.71818f,
        0.71917f, 0.72016f, 0.72114f, 0.72212f,
        0.72309f, 0.72407f, 0.72504f, 0.72601f,
        0.72697f, 0.72794f, 0.72890f, 0.72985f,
        0.73081f, 0.73176f, 0.73271f, 0.73365f,
        0.73460f, 0.73554f, 0.73647f, 0.73741f,
        0.73834f, 0.73927f, 0.74020f, 0.74112f,
        0.74204f, 0.74296f, 0.74388f, 0.74479f,
        0.74570f, 0.74661f, 0.74751f, 0.74842f,
        0.74932f, 0.75021f, 0.75111f, 0.75200f,
        0.75289f, 0.75378f, 0.75466f, 0.75555f,
        0.75643f, 0.75730f, 0.75818f, 0.75905f,
        0.75992f, 0.76079f, 0.76165f, 0.76251f,
        0.76337f, 0.76423f, 0.76508f, 0.76594f,
        0.76679f, 0.76763f, 0.76848f, 0.76932f,
        0.77016f, 0.77100f, 0.77183f, 0.77267f,
        0.77350f, 0.77432f, 0.77515f, 0.77597f,
        0.77680f, 0.77761f, 0.77843f, 0.77924f,
        0.78006f, 0.78087f, 0.78167f, 0.78248f,
        0.78328f, 0.78408f, 0.78488f, 0.78568f,
        0.78647f, 0.78726f, 0.78805f, 0.78884f,
        0.78962f, 0.79040f, 0.79118f, 0.79196f,
        0.79274f, 0.79351f, 0.79428f, 0.79505f,
        0.79582f, 0.79658f, 0.79735f, 0.79811f,
        0.79887f, 0.79962f, 0.80038f, 0.80113f,
        0.80188f, 0.80263f, 0.80337f, 0.80412f,
        0.80486f, 0.80560f, 0.80634f, 0.80707f,
        0.80780f, 0.80854f, 0.80926f, 0.80999f,
        0.81072f, 0.81144f, 0.81216f, 0.81288f,
        0.81360f, 0.81431f, 0.81503f, 0.81574f,
        0.81645f, 0.81715f, 0.81786f, 0.81856f,
        0.81926f, 0.81996f, 0.82066f, 0.82135f,
        0.82205f, 0.82274f, 0.82343f, 0.82412f,
        0.82480f, 0.82549f, 0.82617f, 0.82685f,
        0.82753f, 0.82820f, 0.82888f, 0.82955f,
        0.83022f, 0.83089f, 0.83155f, 0.83222f,
        0.83288f, 0.83354f, 0.83420f, 0.83486f,
        0.83552f, 0.83617f, 0.83682f, 0.83747f,
        0.83812f, 0.83877f, 0.83941f, 0.84005f,
        0.84069f, 0.84133f, 0.84197f, 0.84261f,
        0.84324f, 0.84387f, 0.84450f, 0.84513f,
        0.84576f, 0.84639f, 0.84701f, 0.84763f,
        0.84825f, 0.84887f, 0.84949f, 0.85010f,
        0.85071f, 0.85132f, 0.85193f, 0.85254f,
        0.85315f, 0.85375f, 0.85436f, 0.85496f,
        0.85556f, 0.85615f, 0.85675f, 0.85735f,
        0.85794f, 0.85853f, 0.85912f, 0.85971f,
        0.86029f, 0.86088f, 0.86146f, 0.86204f,
        0.86262f, 0.86320f, 0.86378f, 0.86435f,
        0.86493f, 0.86550f, 0.86607f, 0.86664f,
        0.86720f, 0.86777f, 0.86833f, 0.86889f,
        0.86945f, 0.87001f, 0.87057f, 0.87113f,
        0.87168f, 0.87223f, 0.87278f, 0.87333f,
        0.87388f, 0.87443f, 0.87497f, 0.87552f,
        0.87606f, 0.87660f, 0.87714f, 0.87768f,
        0.87821f, 0.87875f, 0.87928f, 0.87981f,
        0.88034f, 0.88087f, 0.88140f, 0.88192f,
        0.88244f, 0.88297f, 0.88349f, 0.88401f,
        0.88453f, 0.88504f, 0.88556f, 0.88607f,
        0.88658f, 0.88709f, 0.88760f, 0.88811f,
        0.88862f, 0.88912f, 0.88963f, 0.89013f,
        0.89063f, 0.89113f, 0.89163f, 0.89212f,
        0.89262f, 0.89311f, 0.89360f, 0.89409f,
        0.89458f, 0.89507f, 0.89556f, 0.89604f,
        0.89653f, 0.89701f, 0.89749f, 0.89797f,
        0.89845f, 0.89892f, 0.89940f, 0.89987f,
        0.90035f, 0.90082f, 0.90129f, 0.90176f,
        0.90222f, 0.90269f, 0.90316f, 0.90362f,
        0.90408f, 0.90454f, 0.90500f, 0.90546f,
        0.90592f, 0.90637f, 0.90683f, 0.90728f,
        0.90773f, 0.90818f, 0.90863f, 0.90908f,
        0.90952f, 0.90997f, 0.91041f, 0.91085f,
        0.91130f, 0.91173f, 0.91217f, 0.91261f,
        0.91305f, 0.91348f, 0.91392f, 0.91435f,
        0.91478f, 0.91521f, 0.91564f, 0.91606f,
        0.91649f, 0.91691f, 0.91734f, 0.91776f,
        0.91818f, 0.91860f, 0.91902f, 0.91944f,
        0.91985f, 0.92027f, 0.92068f, 0.92109f,
        0.92150f, 0.92191f, 0.92232f, 0.92273f,
        0.92314f, 0.92354f, 0.92395f, 0.92435f,
        0.92475f, 0.92515f, 0.92555f, 0.92595f,
        0.92634f, 0.92674f, 0.92713f, 0.92753f,
        0.92792f, 0.92831f, 0.92870f, 0.92909f,
        0.92947f, 0.92986f, 0.93025f, 0.93063f,
        0.93101f, 0.93139f, 0.93177f, 0.93215f,
        0.93253f, 0.93291f, 0.93328f, 0.93366f,
        0.93403f, 0.93440f, 0.93478f, 0.93515f,
        0.93551f, 0.93588f, 0.93625f, 0.93661f,
        0.93698f, 0.93734f, 0.93770f, 0.93807f,
        0.93843f, 0.93878f, 0.93914f, 0.93950f,
        0.93986f, 0.94021f, 0.94056f, 0.94092f,
        0.94127f, 0.94162f, 0.94197f, 0.94231f,
        0.94266f, 0.94301f, 0.94335f, 0.94369f,
        0.94404f, 0.94438f, 0.94472f, 0.94506f,
        0.94540f, 0.94573f, 0.94607f, 0.94641f,
        0.94674f, 0.94707f, 0.94740f, 0.94774f,
        0.94807f, 0.94839f, 0.94872f, 0.94905f,
        0.94937f, 0.94970f, 0.95002f, 0.95035f,
        0.95067f, 0.95099f, 0.95131f, 0.95163f,
        0.95194f, 0.95226f, 0.95257f, 0.95289f,
        0.95320f, 0.95351f, 0.95383f, 0.95414f,
        0.95445f, 0.95475f, 0.95506f, 0.95537f,
        0.95567f, 0.95598f, 0.95628f, 0.95658f,
        0.95688f, 0.95718f, 0.95748f, 0.95778f,
        0.95808f, 0.95838f, 0.95867f, 0.95897f,
        0.95926f, 0.95955f, 0.95984f, 0.96013f,
        0.96042f, 0.96071f, 0.96100f, 0.96129f,
        0.96157f, 0.96186f, 0.96214f, 0.96242f,
        0.96271f, 0.96299f, 0.96327f, 0.96355f,
        0.96382f, 0.96410f, 0.96438f, 0.96465f,
        0.96493f, 0.96520f, 0.96547f, 0.96574f,
        0.96602f, 0.96629f, 0.96655f, 0.96682f,
        0.96709f, 0.96735f, 0.96762f, 0.96788f,
        0.96815f, 0.96841f, 0.96867f, 0.96893f,
        0.96919f, 0.96945f, 0.96971f, 0.96996f,
        0.97022f, 0.97047f, 0.97073f, 0.97098f,
        0.97123f, 0.97149f, 0.97174f, 0.97199f,
        0.97223f, 0.97248f, 0.97273f, 0.97297f,
        0.97322f, 0.97346f, 0.97371f, 0.97395f,
        0.97419f, 0.97443f, 0.97467f, 0.97491f,
        0.97515f, 0.97539f, 0.97562f, 0.97586f,
        0.97609f, 0.97633f, 0.97656f, 0.97679f,
        0.97702f, 0.97725f, 0.97748f, 0.97771f,
        0.97794f, 0.97817f, 0.97839f, 0.97862f,
        0.97884f, 0.97907f, 0.97929f, 0.97951f,
        0.97973f, 0.97995f, 0.98017f, 0.98039f,
        0.98061f, 0.98082f, 0.98104f, 0.98125f,
        0.98147f, 0.98168f, 0.98189f, 0.98211f,
        0.98232f, 0.98253f, 0.98274f, 0.98295f,
        0.98315f, 0.98336f, 0.98357f, 0.98377f,
        0.98398f, 0.98418f, 0.98438f, 0.98458f,
        0.98478f, 0.98498f, 0.98518f, 0.98538f,
        0.98558f, 0.98578f, 0.98597f, 0.98617f,
        0.98636f, 0.98656f, 0.98675f, 0.98694f,
        0.98714f, 0.98733f, 0.98752f, 0.98771f,
        0.98789f, 0.98808f, 0.98827f, 0.98845f,
        0.98864f, 0.98882f, 0.98901f, 0.98919f,
        0.98937f, 0.98955f, 0.98973f, 0.98991f,
        0.99009f, 0.99027f, 0.99045f, 0.99063f,
        0.99080f, 0.99098f, 0.99115f, 0.99133f,
        0.99150f, 0.99167f, 0.99184f, 0.99201f,
        0.99218f, 0.99235f, 0.99252f, 0.99269f,
        0.99285f, 0.99302f, 0.99319f, 0.99335f,
        0.99351f, 0.99368f, 0.99384f, 0.99400f,
        0.99416f, 0.99432f, 0.99448f, 0.99464f,
        0.99480f, 0.99495f, 0.99511f, 0.99527f,
        0.99542f, 0.99558f, 0.99573f, 0.99588f,
        0.99603f, 0.99619f, 0.99634f, 0.99649f,
        0.99664f, 0.99678f, 0.99693f, 0.99708f,
        0.99722f, 0.99737f, 0.99751f, 0.99766f,
        0.99780f, 0.99794f, 0.99809f, 0.99823f,
        0.99837f, 0.99851f, 0.99865f, 0.99879f,
        0.99892f, 0.99906f, 0.99920f, 0.99933f,
        0.99947f, 0.99960f, 0.99974f, 0.99987f,
        1.00000f
    };

    const size_t sz = sizeof(curve)/sizeof (curve[0]);
    double *tc = malloc(sizeof(tc[0]) * sz);
    for (size_t i = 0; i < sz; i++) {
        tc[i] = (double)curve[i];
    }
    *len = (int)sz;
    return tc;
}

static inline double
curve_interpolate_point(double v,
                        const double tc[],
                        int tc_len)
{
    if (v <= 0) {
        v = 0;
    } else if (v >= 1.0) {
        v = 1.0;
    } else {
        // we assume the tc has evenly spaced points (it should have)
        int k = (int)((tc_len-1) * v);
        double diff = (tc_len-1) * v - (double)k;
        v = (1.0 - diff) * tc[k] + diff * tc[k+1];
    }
    return v;
}

#define adobe_rgb_tone(r, g, b, tc, tc_len)                           \
    do {                                                              \
        double r_ = curve_interpolate_point(*r, tc, tc_len);          \
        double b_ = curve_interpolate_point(*b, tc, tc_len);          \
        double g_ = b_ + ((r_ - b_) * (*g - *b) / (*r - *b));         \
        *r = r_;                                                      \
        *g = g_;                                                      \
        *b = b_;                                                      \
    } while (0);

v3
dcp_apply_tonecurve(v3 rgb,
                    const double tc[],
                    int tc_len)
{
    if (!v3_isfinite(rgb)) {
        elog("Bug: bad RGB input to tone curve.\n");
        v3_print("rgb", rgb);
        abort();
    }
    double *r = &rgb.v[0];
    double *g = &rgb.v[1];
    double *b = &rgb.v[2];

    if (*r >= *g) {
        if (*g > *b) {         // Case 1: r >= g >  b
            adobe_rgb_tone(r, g, b, tc, tc_len);
        } else if (*b > *r) {  // Case 2: b >  r >= g
            adobe_rgb_tone(b, r, g, tc, tc_len);
        } else if (*b > *g) {  // Case 3: r >= b >  g
            adobe_rgb_tone(r, b, g, tc, tc_len);
        } else {               // Case 4: r >= g == b
            *r = curve_interpolate_point(*r, tc, tc_len);
            *g = curve_interpolate_point(*g, tc, tc_len);
            *b = *g;
        }
    } else {
        if (*r >= *b) {        // Case 5: g >  r >= b
            adobe_rgb_tone(g, r, b, tc, tc_len);
        } else if (*b > *g) {  // Case 6: b >  g >  r
            adobe_rgb_tone(b, g, r, tc, tc_len);
        } else {               // Case 7: g >= b >  r
            adobe_rgb_tone(g, b, r, tc, tc_len);
        }
    }
    return (v3){{*r, *g, *b}};
}
