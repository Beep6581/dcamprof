/*
 * (c) Copyright 2015 - 2016 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */

/*
  The spectral data is copyright the respective owners. See the links
  embedded in the associated data files for full terms of use.
*/
#include <inttypes.h>
#include <assert.h>
#include <string.h>
#include <strings.h>

#include <elog.h>
#include <spectraldb.h>

struct color_group {
    const char *name;
    const void *data;
    int band_start;
    int band_end;
    int band_spacing;
    int count;
    spectrum_t **s;
};

static struct {
    spectrum_t *macbeth_cc24[24];
    spectrum_t *munsell_glossy[1600];
    spectrum_t *kuopio_natural[218];
    spectrum_t *joensuu_birch[337];
    spectrum_t *joensuu_spruce[349];
    spectrum_t *joensuu_pine[370];
} spectra;

extern const double spectraldb_macbeth_colorchecker[];
extern const double spectraldb_munsell380_780_1_glossy[];
extern const short spectraldb_kuopio_natural[];
extern const double spectraldb_joensuu_birch[];
extern const double spectraldb_joensuu_spruce[];
extern const double spectraldb_joensuu_pine[];
enum cg_type {
    // enum must match cg indexes
    CG_CC24 = 0,
    CG_MUNSELL = 1,
    CG_KUOPIO_NATURAL = 2,
    /*
    CG_JOENSUU_BIRCH = 3,
    CG_JOENSUU_SPRUCE = 4,
    CG_JOENSUU_PINE = 5,
    */
};
static const struct color_group glob_cg[] = {
    {
        .name = "cc24",
        .data = spectraldb_macbeth_colorchecker,
        .band_start = 380,
        .band_end = 730,
        .band_spacing = 10,
        .count = 24,
        .s = spectra.macbeth_cc24
    },
    {
        .name = "munsell",
        .data = spectraldb_munsell380_780_1_glossy,
        .band_start = 380,
        .band_end = 780,
        .band_spacing = 1,
        .count = 1600,
        .s = spectra.munsell_glossy
    },
    {
        .name = "kuopio-natural",
        .data = spectraldb_kuopio_natural,
        .band_start = 400,
        .band_end = 700,
        .band_spacing = 5,
        .count = 218,
        .s = spectra.kuopio_natural
    },
    /*
    {
        .name = "joensuu-birch",
        .data = spectraldb_joensuu_birch,
        .band_start = 390,
        .band_end = 850,
        .band_spacing = 5,
        .count = 337,
        .s = spectra.joensuu_birch
    },
    {
        .name = "joensuu-spruce",
        .data = spectraldb_joensuu_spruce,
        .band_start = 390,
        .band_end = 850,
        .band_spacing = 5,
        .count = 349,
        .s = spectra.joensuu_spruce
    },
    {
        .name = "joensuu-pine",
        .data = spectraldb_joensuu_pine,
        .band_start = 390,
        .band_end = 850,
        .band_spacing = 5,
        .count = 370,
        .s = spectra.joensuu_pine
    },
    */
    { NULL, NULL, 0,0,0,0, NULL }
};

static const char *munsell_glossy_names[] = {
    "ABB2002", "ABB2004", "ABB2006", "ABB3002", "ABB3004", "ABB3006", "ABB3008", "ABB4002", "ABB4004", "ABB4006", "ABB4008", "ABB5002", "ABB5004", "ABB5006", "ABB5008", "ABB5010", "ABB6002", "ABB6004", "ABB6006", "ABB6008", "ABB7002", "ABB7004", "ABB7006", "ABB7008", "ABB8002", "ABB8004", "ABB9002", "ABG2002", "ABG2004", "ABG2006", "ABG3002", "ABG3004", "ABG3006", "ABG3008", "ABG4002", "ABG4004", "ABG4006", "ABG4008", "ABG5002", "ABG5004", "ABG5006", "ABG5008", "ABG5010", "ABG6002", "ABG6004", "ABG6006", "ABG6008", "ABG6010", "ABG7002", "ABG7004", "ABG7006", "ABG7008", "ABG8002", "ABG8004", "ABG8006", "ABG9002", "AGG2002", "AGG2004", "AGG3002", "AGG3004", "AGG3006", "AGG3008", "AGG3010", "AGG4002", "AGG4004", "AGG4006", "AGG4008", "AGG4010", "AGG5002", "AGG5004", "AGG5006", "AGG5008", "AGG5010", "AGG5012", "AGG6002", "AGG6004", "AGG6006", "AGG6008", "AGG6010", "AGG7002", "AGG7004", "AGG7006", "AGG7008", "AGG7010", "AGG8002", "AGG8004", "AGG8006", "AGG8008", "AGG9002", "AGG9004", "AGY2002", "AGY3002", "AGY3004", "AGY4002", "AGY4004", "AGY4006", "AGY5002", "AGY5004", "AGY5006", "AGY5008", "AGY6002", "AGY6004", "AGY6006", "AGY6008", "AGY6010", "AGY7002", "AGY7004", "AGY7006", "AGY7008", "AGY7010", "AGY7012", "AGY8002", "AGY8004", "AGY8006", "AGY8008", "AGY8010", "AGY8012", "AGY9002", "AGY9004", "AGY9006", "APB2002", "APB2004", "APB2006", "APB3002", "APB3004", "APB3006", "APB3008", "APB3010", "APB4002", "APB4004", "APB4006", "APB4008", "APB4010", "APB5002", "APB5004", "APB5006", "APB5008", "APB5010", "APB5012", "APB6002", "APB6004", "APB6006", "APB6008", "APB6010", "APB7002", "APB7004", "APB7006", "APB7008", "APB8002", "APB8004", "APB8006", "APB9002", "APP2002", "APP2004", "APP2006", "APP2008", "APP2010", "APP3002", "APP3004", "APP3006", "APP3008", "APP3010", "APP4002", "APP4004", "APP4006", "APP4008", "APP4010", "APP4012", "APP5002", "APP5004", "APP5006", "APP5008", "APP5010", "APP6002", "APP6004", "APP6006", "APP6008", "APP7002", "APP7004", "APP7006", "APP7008", "APP8002", "APP8004", "APP9002", "ARP2002", "ARP2004", "ARP2006", "ARP2008", "ARP3002", "ARP3004", "ARP3006", "ARP3008", "ARP3010", "ARP4002", "ARP4004", "ARP4006", "ARP4008", "ARP4010", "ARP4012", "ARP5002", "ARP5004", "ARP5006", "ARP5008", "ARP5010", "ARP5012", "ARP6002", "ARP6004", "ARP6006", "ARP6008", "ARP6010", "ARP6012", "ARP7002", "ARP7004", "ARP7006", "ARP7008", "ARP7010", "ARP8002", "ARP8004", "ARP8006", "ARP9002", "ARR2002", "ARR2004", "ARR2006", "ARR2008", "ARR3002", "ARR3004", "ARR3006", "ARR3008", "ARR3010", "ARR4002", "ARR4004", "ARR4006", "ARR4008", "ARR4010", "ARR4012", "ARR4014", "ARR5002", "ARR5004", "ARR5006", "ARR5008", "ARR5010", "ARR5012", "ARR5014", "ARR6002", "ARR6004", "ARR6006", "ARR6008", "ARR6010", "ARR6012", "ARR7002", "ARR7004", "ARR7006", "ARR7008", "ARR8002", "ARR8004", "ARR8006", "ARR9002", "AYR2002", "AYR2004", "AYR3002", "AYR3004", "AYR3006", "AYR3008", "AYR4002", "AYR4004", "AYR4006", "AYR4008", "AYR4010", "AYR5002", "AYR5004", "AYR5006", "AYR5008", "AYR5010", "AYR5012", "AYR5014", "AYR6002", "AYR6004", "AYR6006", "AYR6008", "AYR6010", "AYR6012", "AYR6014", "AYR6016", "AYR7002", "AYR7004", "AYR7006", "AYR7008", "AYR7010", "AYR7012", "AYR8002", "AYR8004", "AYR8006", "AYR9002", "AYY2002", "AYY3002", "AYY3004", "AYY4002", "AYY4004", "AYY4006", "AYY5002", "AYY5004", "AYY5006", "AYY5008", "AYY6002", "AYY6004", "AYY6006", "AYY6008", "AYY6010", "AYY7002", "AYY7004", "AYY7006", "AYY7008", "AYY7010", "AYY7012", "AYY8002", "AYY8004", "AYY8006", "AYY8008", "AYY8010", "AYY8012", "AYY8014", "AYY8016", "AYY8502", "AYY8504", "AYY8506", "AYY8508", "AYY8510", "AYY8512", "AYY9002", "AYY9004", "BBB2002", "BBB2004", "BBB2006", "BBB3002", "BBB3004", "BBB3006", "BBB3008", "BBB4002", "BBB4004", "BBB4006", "BBB4008", "BBB4010", "BBB5002", "BBB5004", "BBB5006", "BBB5008", "BBB5010", "BBB6002", "BBB6004", "BBB6006", "BBB6008", "BBB6010", "BBB7002", "BBB7004", "BBB7006", "BBB7008", "BBB8002", "BBB8004", "BBB9002", "BBG2002", "BBG2004", "BBG2006", "BBG3002", "BBG3004", "BBG3006", "BBG3008", "BBG4002", "BBG4004", "BBG4006", "BBG4008", "BBG5002", "BBG5004", "BBG5006", "BBG5008", "BBG5010", "BBG6002", "BBG6004", "BBG6006", "BBG6008", "BBG6010", "BBG7002", "BBG7004", "BBG7006", "BBG7008", "BBG8002", "BBG8004", "BBG9002", "BGG2002", "BGG2004", "BGG2006", "BGG3002", "BGG3004", "BGG3006", "BGG3008", "BGG4002", "BGG4004", "BGG4006", "BGG4008", "BGG4010", "BGG5002", "BGG5004", "BGG5006", "BGG5008", "BGG5010", "BGG6002", "BGG6004", "BGG6006", "BGG6008", "BGG6010", "BGG7002", "BGG7004", "BGG7006", "BGG7008", "BGG7010", "BGG8002", "BGG8004", "BGG8006", "BGG9002", "BGY2002", "BGY3002", "BGY3004", "BGY3006", "BGY4002", "BGY4004", "BGY4006", "BGY4008", "BGY5002", "BGY5004", "BGY5006", "BGY5008", "BGY5010", "BGY6002", "BGY6004", "BGY6006", "BGY6008", "BGY6010", "BGY7002", "BGY7004", "BGY7006", "BGY7008", "BGY7010", "BGY7012", "BGY8002", "BGY8004", "BGY8006", "BGY8008", "BGY8010", "BGY9002", "BGY9004", "BPB2002", "BPB2004", "BPB2006", "BPB2008", "BPB3002", "BPB3004", "BPB3006", "BPB3008", "BPB3010", "BPB4002", "BPB4004", "BPB4006", "BPB4008", "BPB4010", "BPB4012", "BPB5002", "BPB5004", "BPB5006", "BPB5008", "BPB5010", "BPB5012", "BPB6002", "BPB6004", "BPB6006", "BPB6008", "BPB6010", "BPB7002", "BPB7004", "BPB7006", "BPB7008", "BPB8002", "BPB8004", "BPB8006", "BPB9002", "BPP2002", "BPP2004", "BPP2006", "BPP2008", "BPP3002", "BPP3004", "BPP3006", "BPP3008", "BPP3010", "BPP4002", "BPP4004", "BPP4006", "BPP4008", "BPP4010", "BPP4012", "BPP5002", "BPP5004", "BPP5006", "BPP5008", "BPP5010", "BPP6002", "BPP6004", "BPP6006", "BPP6008", "BPP7002", "BPP7004", "BPP7006", "BPP7008", "BPP8002", "BPP8004", "BPP9002", "BRP2002", "BRP2004", "BRP2006", "BRP2008", "BRP3002", "BRP3004", "BRP3006", "BRP3008", "BRP3010", "BRP4002", "BRP4004", "BRP4006", "BRP4008", "BRP4010", "BRP4012", "BRP5002", "BRP5004", "BRP5006", "BRP5008", "BRP5010", "BRP5012", "BRP6002", "BRP6004", "BRP6006", "BRP6008", "BRP6010", "BRP6012", "BRP7002", "BRP7004", "BRP7006", "BRP7008", "BRP7010", "BRP8002", "BRP8004", "BRP8006", "BRP9002", "BRR2002", "BRR2004", "BRR2006", "BRR2008", "BRR3002", "BRR3004", "BRR3006", "BRR3008", "BRR3010", "BRR4002", "BRR4004", "BRR4006", "BRR4008", "BRR4010", "BRR4012", "BRR4014", "BRR5002", "BRR5004", "BRR5006", "BRR5008", "BRR5010", "BRR5012", "BRR5014", "BRR6002", "BRR6004", "BRR6006", "BRR6008", "BRR6010", "BRR6012", "BRR7002", "BRR7004", "BRR7006", "BRR7008", "BRR7010", "BRR8002", "BRR8004", "BRR8006", "BRR9002", "BYR2002", "BYR2004", "BYR3002", "BYR3004", "BYR3006", "BYR4002", "BYR4004", "BYR4006", "BYR4008", "BYR5002", "BYR5004", "BYR5006", "BYR5008", "BYR5010", "BYR5012", "BYR6002", "BYR6004", "BYR6006", "BYR6008", "BYR6010", "BYR6012", "BYR6014", "BYR7002", "BYR7004", "BYR7006", "BYR7008", "BYR7010", "BYR7012", "BYR7014", "BYR8002", "BYR8004", "BYR8006", "BYR8008", "BYR9002", "BYY2002", "BYY3002", "BYY3004", "BYY4002", "BYY4004", "BYY4006", "BYY5002", "BYY5004", "BYY5006", "BYY5008", "BYY6002", "BYY6004", "BYY6006", "BYY6008", "BYY6010", "BYY7002", "BYY7004", "BYY7006", "BYY7008", "BYY7010", "BYY7012", "BYY8002", "BYY8004", "BYY8006", "BYY8008", "BYY8010", "BYY8012", "BYY8014", "BYY8502", "BYY8504", "BYY8506", "BYY8508", "BYY8510", "BYY8512", "BYY8514", "BYY9002", "BYY9004", "BYY9006", "CBB2002", "CBB2004", "CBB2006", "CBB3002", "CBB3004", "CBB3006", "CBB3008", "CBB4002", "CBB4004", "CBB4006", "CBB4008", "CBB4010", "CBB5002", "CBB5004", "CBB5006", "CBB5008", "CBB5010", "CBB6002", "CBB6004", "CBB6006", "CBB6008", "CBB6010", "CBB7002", "CBB7004", "CBB7006", "CBB7008", "CBB8002", "CBB8004", "CBB9002", "CBG2002", "CBG2004", "CBG2006", "CBG3002", "CBG3004", "CBG3006", "CBG3008", "CBG4002", "CBG4004", "CBG4006", "CBG4008", "CBG5002", "CBG5004", "CBG5006", "CBG5008", "CBG5010", "CBG6002", "CBG6004", "CBG6006", "CBG6008", "CBG7002", "CBG7004", "CBG7006", "CBG7008", "CBG8002", "CBG8004", "CBG9002", "CGG2002", "CGG2004", "CGG2006", "CGG3002", "CGG3004", "CGG3006", "CGG3008", "CGG4002", "CGG4004", "CGG4006", "CGG4008", "CGG4010", "CGG5002", "CGG5004", "CGG5006", "CGG5008", "CGG5010", "CGG6002", "CGG6004", "CGG6006", "CGG6008", "CGG6010", "CGG7002", "CGG7004", "CGG7006", "CGG7008", "CGG8002", "CGG8004", "CGG8006", "CGG9002", "CGY2002", "CGY2004", "CGY3002", "CGY3004", "CGY3006", "CGY4002", "CGY4004", "CGY4006", "CGY4008", "CGY5002", "CGY5004", "CGY5006", "CGY5008", "CGY5010", "CGY6002", "CGY6004", "CGY6006", "CGY6008", "CGY6010", "CGY6012", "CGY7002", "CGY7004", "CGY7006", "CGY7008", "CGY7010", "CGY8002", "CGY8004", "CGY8006", "CGY8008", "CGY8010", "CGY9002", "CGY9004", "CPB2002", "CPB2004", "CPB2006", "CPB2008", "CPB2010", "CPB3002", "CPB3004", "CPB3006", "CPB3008", "CPB3010", "CPB3012", "CPB4002", "CPB4004", "CPB4006", "CPB4008", "CPB4010", "CPB4012", "CPB5002", "CPB5004", "CPB5006", "CPB5008", "CPB5010", "CPB6002", "CPB6004", "CPB6006", "CPB6008", "CPB6010", "CPB7002", "CPB7004", "CPB7006", "CPB7008", "CPB8002", "CPB8004", "CPB9002", "CPP2002", "CPP2004", "CPP2006", "CPP3002", "CPP3004", "CPP3006", "CPP3008", "CPP3010", "CPP4002", "CPP4004", "CPP4006", "CPP4008", "CPP4010", "CPP4012", "CPP5002", "CPP5004", "CPP5006", "CPP5008", "CPP5010", "CPP6002", "CPP6004", "CPP6006", "CPP6008", "CPP6010", "CPP7002", "CPP7004", "CPP7006", "CPP7008", "CPP8002", "CPP8004", "CPP8006", "CPP9002", "CRP2002", "CRP2004", "CRP2006", "CRP2008", "CRP3002", "CRP3004", "CRP3006", "CRP3008", "CRP3010", "CRP4002", "CRP4004", "CRP4006", "CRP4008", "CRP4010", "CRP4012", "CRP5002", "CRP5004", "CRP5006", "CRP5008", "CRP5010", "CRP5012", "CRP5014", "CRP6002", "CRP6004", "CRP6006", "CRP6008", "CRP6010", "CRP6012", "CRP7002", "CRP7004", "CRP7006", "CRP7008", "CRP8002", "CRP8004", "CRP8006", "CRP9002", "CRR2002", "CRR2004", "CRR2006", "CRR2008", "CRR3002", "CRR3004", "CRR3006", "CRR3008", "CRR3010", "CRR3012", "CRR4002", "CRR4004", "CRR4006", "CRR4008", "CRR4010", "CRR4012", "CRR4014", "CRR4016", "CRR5002", "CRR5004", "CRR5006", "CRR5008", "CRR5010", "CRR5012", "CRR5014", "CRR5016", "CRR6002", "CRR6004", "CRR6006", "CRR6008", "CRR6010", "CRR6012", "CRR7002", "CRR7004", "CRR7006", "CRR7008", "CRR7010", "CRR8002", "CRR8004", "CRR8006", "CRR9002", "CYR2002", "CYR2004", "CYR3002", "CYR3004", "CYR3006", "CYR4002", "CYR4004", "CYR4006", "CYR4008", "CYR5002", "CYR5004", "CYR5006", "CYR5008", "CYR5010", "CYR6002", "CYR6004", "CYR6006", "CYR6008", "CYR6010", "CYR6012", "CYR6014", "CYR7002", "CYR7004", "CYR7006", "CYR7008", "CYR7010", "CYR7012", "CYR7014", "CYR7016", "CYR8002", "CYR8004", "CYR8006", "CYR8008", "CYR8010", "CYR9002", "CYY2002", "CYY3002", "CYY3004", "CYY4002", "CYY4004", "CYY4006", "CYY5002", "CYY5004", "CYY5006", "CYY5008", "CYY6002", "CYY6004", "CYY6006", "CYY6008", "CYY6010", "CYY7002", "CYY7004", "CYY7006", "CYY7008", "CYY7010", "CYY7012", "CYY8002", "CYY8004", "CYY8006", "CYY8008", "CYY8010", "CYY8012", "CYY8502", "CYY8504", "CYY8506", "CYY8508", "CYY8510", "CYY8512", "CYY9002", "CYY9004", "CYY9006", "DBB2002", "DBB2004", "DBB2006", "DBB3002", "DBB3004", "DBB3006", "DBB3008", "DBB3010", "DBB4002", "DBB4004", "DBB4006", "DBB4008", "DBB4010", "DBB5002", "DBB5004", "DBB5006", "DBB5008", "DBB5010", "DBB5012", "DBB6002", "DBB6004", "DBB6006", "DBB6008", "DBB6010", "DBB7002", "DBB7004", "DBB7006", "DBB7008", "DBB8002", "DBB8004", "DBB8006", "DBB9002", "DBG2002", "DBG2004", "DBG2006", "DBG3002", "DBG3004", "DBG3006", "DBG3008", "DBG4002", "DBG4004", "DBG4006", "DBG4008", "DBG5002", "DBG5004", "DBG5006", "DBG5008", "DBG5010", "DBG6002", "DBG6004", "DBG6006", "DBG6008", "DBG7002", "DBG7004", "DBG7006", "DBG7008", "DBG8002", "DBG8004", "DBG9002", "DGG2002", "DGG2004", "DGG2006", "DGG3002", "DGG3004", "DGG3006", "DGG3008", "DGG4002", "DGG4004", "DGG4006", "DGG4008", "DGG4010", "DGG5002", "DGG5004", "DGG5006", "DGG5008", "DGG5010", "DGG6002", "DGG6004", "DGG6006", "DGG6008", "DGG6010", "DGG7002", "DGG7004", "DGG7006", "DGG7008", "DGG8002", "DGG8004", "DGG8006", "DGG9002", "DGY2002", "DGY2004", "DGY3002", "DGY3004", "DGY3006", "DGY4002", "DGY4004", "DGY4006", "DGY4008", "DGY5002", "DGY5004", "DGY5006", "DGY5008", "DGY5010", "DGY5012", "DGY6002", "DGY6004", "DGY6006", "DGY6008", "DGY6010", "DGY6012", "DGY7002", "DGY7004", "DGY7006", "DGY7008", "DGY7010", "DGY8002", "DGY8004", "DGY8006", "DGY8008", "DGY9002", "DGY9004", "DPB2002", "DPB2004", "DPB2006", "DPB2008", "DPB2010", "DPB3002", "DPB3004", "DPB3006", "DPB3008", "DPB3010", "DPB4002", "DPB4004", "DPB4006", "DPB4008", "DPB4010", "DPB4012", "DPB5002", "DPB5004", "DPB5006", "DPB5008", "DPB5010", "DPB6002", "DPB6004", "DPB6006", "DPB6008", "DPB6010", "DPB7002", "DPB7004", "DPB7006", "DPB7008", "DPB8002", "DPB8004", "DPB9002", "DPP2002", "DPP2004", "DPP2006", "DPP3002", "DPP3004", "DPP3006", "DPP3008", "DPP3010", "DPP4002", "DPP4004", "DPP4006", "DPP4008", "DPP4010", "DPP4012", "DPP5002", "DPP5004", "DPP5006", "DPP5008", "DPP5010", "DPP5012", "DPP6002", "DPP6004", "DPP6006", "DPP6008", "DPP6010", "DPP7002", "DPP7004", "DPP7006", "DPP7008", "DPP8002", "DPP8004", "DPP8006", "DPP9002", "DRP2002", "DRP2004", "DRP2006", "DRP2008", "DRP3002", "DRP3004", "DRP3006", "DRP3008", "DRP3010", "DRP4002", "DRP4004", "DRP4006", "DRP4008", "DRP4010", "DRP4012", "DRP4014", "DRP5002", "DRP5004", "DRP5006", "DRP5008", "DRP5010", "DRP5012", "DRP5014", "DRP6002", "DRP6004", "DRP6006", "DRP6008", "DRP6010", "DRP6012", "DRP7002", "DRP7004", "DRP7006", "DRP7008", "DRP8002", "DRP8004", "DRP8006", "DRP9002", "DRR2002", "DRR2004", "DRR2006", "DRR3002", "DRR3004", "DRR3006", "DRR3008", "DRR3010", "DRR4002", "DRR4004", "DRR4006", "DRR4008", "DRR4010", "DRR4012", "DRR5002", "DRR5004", "DRR5006", "DRR5008", "DRR5010", "DRR5012", "DRR5014", "DRR5016", "DRR6002", "DRR6004", "DRR6006", "DRR6008", "DRR6010", "DRR6012", "DRR6014", "DRR7002", "DRR7004", "DRR7006", "DRR7008", "DRR7010", "DRR8002", "DRR8004", "DRR8006", "DRR9002", "DYR2002", "DYR3002", "DYR3004", "DYR3006", "DYR4002", "DYR4004", "DYR4006", "DYR4008", "DYR5002", "DYR5004", "DYR5006", "DYR5008", "DYR5010", "DYR6002", "DYR6004", "DYR6006", "DYR6008", "DYR6010", "DYR6012", "DYR7002", "DYR7004", "DYR7006", "DYR7008", "DYR7010", "DYR7012", "DYR7014", "DYR8002", "DYR8004", "DYR8006", "DYR8008", "DYR8010", "DYR8012", "DYR8014", "DYR9002", "DYR9004", "DYY2002", "DYY3002", "DYY3004", "DYY4002", "DYY4004", "DYY4006", "DYY5002", "DYY5004", "DYY5006", "DYY5008", "DYY6002", "DYY6004", "DYY6006", "DYY6008", "DYY6010", "DYY7002", "DYY7004", "DYY7006", "DYY7008", "DYY7010", "DYY7012", "DYY8002", "DYY8004", "DYY8006", "DYY8008", "DYY8010", "DYY8012", "DYY8502", "DYY8504", "DYY8506", "DYY8508", "DYY8510", "DYY8512", "DYY9002", "DYY9004", "DYY9006", "DPB7001", "BBB3001", "BBB4001", "BBB5001", "BBB6001", "BBB7001", "BBB8001", "BBB9001", "BBG2001", "BBG3001", "BBG4001", "BBG5001", "BBG6001", "BBG7001", "BBG8001", "BBG9001", "BPB2001", "BPB3001", "BPB4001", "BPB5001", "BPB6001", "BPB7001", "BPB8001", "BPB9001", "DBB2001", "DBB3001", "DBB4001", "DBB5001", "DBB6001", "DBB7001", "DBB8001", "DBB9001", "DBG2001", "DBG3001", "DBG4001", "DBG5001", "DBG6001", "DBG7001", "DBG8001", "DBG9001", "DPB2001", "DPB3001", "DPB4001", "DPB5001", "DPB6001", "BBB2001", "DPB8001", "DPB9001", "BGY2001", "BGY3001", "BGY4001", "BGY5001", "BGY6001", "BGY7001", "BGY8001", "BGY9001", "DGY2001", "DGY3001", "DGY4001", "DGY5001", "DGY6001", "DGY7001", "DGY8001", "DGY9001", "BGG2001", "BGG3001", "BGG4001", "BGG5001", "BGG6001", "BGG7001", "BGG8001", "BGG9001", "DGG2001", "DGG3001", "DGG4001", "DGG5001", "DGG6001", "DGG7001", "DGG8001", "DGG9001", "BPP2001", "BPP3001", "BPP4001", "BPP5001", "BPP6001", "BPP7001", "BPP8001", "BPP9001", "DPP2001", "DPP3001", "DPP4001", "DPP5001", "DPP6001", "DPP7001", "DPP8001", "DPP9001", "BRP2001", "BRP3001", "BRP4001", "BRP5001", "BRP6001", "BRP7001", "BRP8001", "BRP9001", "DRP2001", "DRP3001", "DRP4001", "DRP5001", "DRP6001", "DRP7001", "DRP8001", "DRP9001", "BRR2001", "BRR3001", "BRR4001", "BRR5001", "BRR6001", "BRR7001", "BRR8001", "BRR9001", "DRR2001", "DRR3001", "DRR4001", "DRR5001", "DRR6001", "DRR7001", "DRR8001", "DRR9001", "BYR2001", "BYR3001", "BYR4001", "BYR5001", "BYR6001", "BYR7001", "BYR8001", "BYR9001", "DYR2001", "DYR3001", "DYR4001", "DYR5001", "DYR6001", "DYR7001", "DYR8001", "DYR9001", "BYY2001", "BYY3001", "BYY4001", "BYY5001", "BYY6001", "BYY7001", "BYY8001", "BYY9001", "DYY2001", "DYY3001", "DYY4001", "DYY5001", "DYY6001", "DYY7001", "DYY8001", "DYY9001", "NEUT050", "NEUT075", "NEUT100", "NEUT125", "NEUT150", "NEUT175", "NEUT200", "NEUT225", "NEUT250", "NEUT275", "NEUT300", "NEUT325", "NEUT350", "NEUT375", "NEUT400", "NEUT425", "NEUT450", "NEUT475", "NEUT500", "NEUT525", "NEUT550", "NEUT575", "NEUT600", "NEUT625", "NEUT650", "NEUT675", "NEUT700", "NEUT725", "NEUT750", "NEUT775", "NEUT800", "NEUT825", "NEUT850", "NEUT875", "NEUT900", "NEUT925", "NEUT950", "ERR4012", "ERR4014", "ERR5012", "ERR5014", "ERR6012", "EYR4012", "EYR5012", "EYR6012", "EYR6014", "FRR4012", "FRR4014", "FRR5012", "FRR5014", "FRR6012", "FYR5012", "FYR6012", "FYR6014", "FYR7012", "GRR3012", "GRR4012", "GRR4014", "GRR5012", "GRR5014", "GRR6012", "GYR6012", "GYR7012", "GYR7014", "GYR7016", "HRR4012", "HRR4014", "HRR5012", "HRR5014", "HRR5016", "HRR6012", "HRR6014", "HYR6012", "HYR7012", "HYR7014", "EYY6012", "EYY7012", "EYY7014", "EYY8012", "EYY8014", "EYY8016", "EGY8012", "EGG5012", "EGG6012", "FYY7012", "FYY8012", "FYY8014", "FYY8512", "GYY7012", "GYY8012", "GYY8512", "HYY7012", "HYY8012", "HYY8512", "HGY6012", "FPB4012", "FPB5012", "GPB3012", "GPB4012", "GPB4014", "GPB5012", "HPB4012", "ERP4012", "FRP4012", "FRP5012", "FRP6012", "GRP4012", "GRP5012", "GRP5014", "GRP6012", "HRP4012", "HRP4014", "HRP5012", "HRP5014", "HRP6012"
};

static const struct color_group *
get_cg(const char group[])
{
    for (int i = 0; glob_cg[i].name != NULL; i++) {
        if (strcasecmp(group, glob_cg[i].name) == 0) {
            return &glob_cg[i];
        }
    }
    return NULL;
}

static enum cg_type
cg_type(const struct color_group *cg)
{
    return (enum cg_type)((uintptr_t)cg - (uintptr_t)glob_cg) / sizeof(*cg);
}

const spectrum_t *
spectraldb_color(const char group[],
                 const char name[])
{
    const struct color_group *cg = get_cg(group);
    if (cg == NULL) {
        return NULL;
    }
    // FIXME: make proper indexing rather than just scanning
    switch (cg_type(cg)) {
    case CG_MUNSELL:
        for (int i = 0; i < cg->count; i++) {
            if (strcmp(name, munsell_glossy_names[i]) == 0) {
                return spectraldb_color_byidx(group, i, NULL);
            }
        }
        break;
    case CG_CC24: {
        int col, row;
        switch (name[0]) {
        case 'A': case 'a': row = 0; break;
        case 'B': case 'b': row = 1; break;
        case 'C': case 'c': row = 2; break;
        case 'D': case 'd': row = 3; break;
        default: return NULL;
        }
        switch (name[2]) {
        case '1': col = 0; break;
        case '2': col = 1; break;
        case '3': col = 2; break;
        case '4': col = 3; break;
        case '5': col = 3; break;
        case '6': col = 3; break;
        default: return NULL;
        }
        if (name[1] != '0' || name[3] != '\0') return NULL;
        return spectraldb_color_byidx(group, col*4+row, NULL);
        break;
    }
    default:
        elog("No named indexing for color group \"%s\".\n", group);
        break;
    }
    return NULL;
}

const spectrum_t *
spectraldb_color_byidx(const char group[],
                       int idx,
                       char name[])
{
    const struct color_group *cg = get_cg(group);
    if (cg == NULL) {
        return NULL;
    }
    if (idx < 0 || idx >= cg->count) {
        return NULL;
    }

    if (name != NULL) {
        // names
        switch (cg_type(cg)) {
        case CG_MUNSELL:
            memcpy(name, munsell_glossy_names[idx], 8);
            break;
        case CG_CC24:
            name[0] = 'A' + idx / 6;
            name[1] = '0';
            name[2] = '1' + idx % 6;
            name[3] = '\0';
            break;
        default:
            name[0] = '\0';
            break;
        }
    }

    if (cg->s[idx] != NULL) {
        return cg->s[idx];
    }

    // load spectrum
    const int band_count = (cg->band_end - cg->band_start) / cg->band_spacing + 1;
    double amp_[band_count];
    const double *amp = amp_;
    const double *data_f64 = (const double *)cg->data;
    const short *data_i16 = (const short *)cg->data;
    switch (cg_type(cg)) {
    case CG_MUNSELL:
        for (int i = 0; i < band_count; i++) {
            amp_[i] = data_f64[idx + i * cg->count];
        }
        break;
    case CG_KUOPIO_NATURAL:
        for (int i = 0; i < band_count; i++) {
            short val = data_i16[idx * band_count + i];
            if (val < 0) val = 0;
            if (val > 4095) val = 4095;
            amp_[i] = (double)val / 4095.0;
        }
        break;
    case CG_CC24: {
        // transpose, data is in 4x6 rather than 6x4
        int tidx[24] = { 0, 4, 8,  12, 16, 20,
                         1, 5, 9,  13, 17, 21,
                         2, 6, 10, 14, 18, 22,
                         3, 7, 11, 15, 19, 23 };
        amp = &data_f64[tidx[idx] * band_count];
        break;
    }
    default:
        amp = &data_f64[idx * band_count];
        break;
    }

    cg->s[idx] = spec_alloc(amp, cg->band_start, cg->band_end, cg->band_spacing);

    return cg->s[idx];
}

int
spectraldb_group_count(const char group[])
{
    const struct color_group *cg = get_cg(group);
    return cg ? cg->count : 0;
}

static double
blackbody_f(double temp, // in K
            double wl) // in nm
{
    double L = (double)wl / 1000000000.0;
    double f = 3.7418e-16 / ((pow(L,5)) * (exp(0.014388/(L*temp))-1));
    f /= 1e8;
    return f;
}

spectrum_t *
spectraldb_render_blackbody(double temp,
                            int band_start,
                            int band_end,
                            int band_spacing)
{
    int band_count = (band_end - band_start) / band_spacing + 1;
    double amp[band_count];
    double ref = blackbody_f(temp, 560);
    for (int i = 0; i < band_count; i++) {
        double wl = (band_start + i * band_spacing);
        amp[i] = blackbody_f(temp, wl) / ref;
    }
    return spec_alloc(amp, band_start, band_end, band_spacing);
}

struct ill_spectrum {
    const char *name;
    enum exif_lightsource exif;
    int band_start;
    int band_end;
    int band_spacing;
    double scale;
#define ILL_SPECTRUM_MAX_BAND_COUNT 400
    const double data[ILL_SPECTRUM_MAX_BAND_COUNT];
};

static const struct ill_spectrum ill_list[] = {
    {
        .name = "StdA",
        .exif = lsStandardLightA,
        .band_start = 300,
        .band_end = 830,
        .band_spacing = 5,
        .scale = 100,
        .data = {
            0.930483, 1.128210, 1.357690, 1.622190, 1.925080,
            2.269800, 2.659810, 3.098610, 3.589680, 4.136480,
            4.742380, 5.410700, 6.144620, 6.947200, 7.821350,
            8.769800, 9.795100, 10.899600, 12.085300, 13.354300,
            14.708, 16.148, 17.6753, 19.2907, 20.995,
            22.7883, 24.6709, 26.6425, 28.7027, 30.8508,
            33.0859, 35.4068, 37.8121, 40.3002, 42.8693,
            45.5174, 48.2423, 51.0418, 53.9132, 56.8539,
            59.8611, 62.9320, 66.0635, 69.2525, 72.4959,
            75.7903, 79.1326, 82.5193,  85.947, 89.4124,
             92.912, 96.4423, 100.000, 103.582, 107.184,
            110.803, 114.436, 118.080, 121.731, 125.386,
            129.043, 132.697, 136.346, 139.988, 143.618,
            147.235, 150.836, 154.418, 157.979, 161.516,
            165.028, 168.510, 171.963, 175.383, 178.769,
            182.118, 185.429, 188.701, 191.931, 195.118,
            198.261, 201.359, 204.409, 207.411, 210.365,
            213.268, 216.120, 218.920, 221.667, 224.361,
            227.000, 229.585, 232.115, 234.589, 237.008,
            239.370, 241.675, 243.924, 246.116, 248.251,
            250.329, 252.350, 254.314, 256.221, 258.071,
            259.865, 261.602
        }
    },
    {
        .name = "D50",
        .exif = lsD50,
        .band_start = 300,
        .band_end = 830,
        .band_spacing = 5,
        .scale = 100,
        .data = {
            0.02, 1.03, 2.05, 4.91, 7.78, 11.26, 14.75, 16.35, 17.95, 19.48,
            21.01, 22.48, 23.94, 25.45, 26.96, 25.72, 24.49, 27.18, 29.87, 39.59,
            49.31, 52.91, 56.51, 58.27, 60.03, 58.93, 57.82, 66.32, 74.82, 81.04,
            87.25, 88.93, 90.61, 90.99, 91.37, 93.24, 95.11, 93.54, 91.96, 93.84,
            95.72, 96.17, 96.61, 96.87, 97.13, 99.61, 102.10, 101.43, 100.75, 101.54,
            102.32, 101.16, 100.00, 98.87, 97.74, 98.33, 98.92, 96.21, 93.50, 95.59,
            97.69, 98.48, 99.27, 99.16, 99.04, 97.38, 95.72, 97.29, 98.86, 97.26,
            95.67, 96.93, 98.19, 100.60, 103.00, 101.07, 99.13, 93.26, 87.38, 89.49,
            91.60, 92.25, 92.89, 84.87, 76.85, 81.68, 86.51, 89.55, 92.58, 85.40,
            78.23, 67.96, 57.69, 70.31, 82.92, 80.60, 78.27, 78.91, 79.55, 76.48,
            73.40, 68.66, 63.92, 67.35, 70.78, 72.61, 74.44
        }
    },
    {
        .name = "D55",
        .exif = lsD55,
        .band_start = 300,
        .band_end = 780,
        .band_spacing = 5,
        .scale = 100,
        .data = {
            0.024,1.048,2.072,6.648,11.224,15.936,20.647,22.266,23.885,25.851,27.817,29.219,30.621,32.464,34.308,33.446,32.584,35.335,38.087,49.518,60.949,64.751,68.554,70.065,71.577,69.746,67.914,76.760,85.605,91.799,97.993,99.228,100.463,100.188,99.913,101.326,102.739,100.409,98.078,99.379,100.680,100.688,100.695,100.341,99.987,102.098,104.210,103.156,102.102,102.535,102.968,101.484,100.000,98.608,97.216,97.482,97.749,94.590,91.432,92.926,94.419,94.780,95.140,94.680,94.220,92.334,90.448,91.389,92.330,90.592,88.854,89.586,90.317,92.133,93.950,91.953,89.956,84.817,79.677,81.258,82.840,83.842,84.844,77.539,70.235,74.768,79.301,82.147,84.993,78.437,71.880,62.337,52.793,64.360,75.927,73.872,71.818
        }
    },
    {
        .name = "D65",
        .exif = lsD65,
        .band_start = 300,
        .band_end = 830,
        .band_spacing = 5,
        .scale = 100,
        .data = {
                0.03410, 1.66430, 3.29450, 11.76520, 20.23600,
                28.64470, 37.05350, 38.50110, 39.94880, 42.43020,
                44.91170, 45.77500, 46.63830, 49.36370, 52.08910,
                51.03230, 49.97550, 52.31180, 54.64820, 68.70150,
                82.75490, 87.12040, 91.48600, 92.45890, 93.43180,
                90.05700, 86.68230, 95.77360, 104.86500, 110.93600,
                117.00800, 117.41000, 117.81200, 116.33600, 114.86100,
                115.39200, 115.92300, 112.36700, 108.81100, 109.08200,
                109.35400, 108.57800, 107.80200, 106.29600, 104.79000,
                106.23900, 107.68900, 106.04700, 104.40500, 104.22500,
                104.04600, 102.02300, 100.00000, 98.16710, 96.33420,
                96.06110, 95.78800, 92.23680, 88.68560, 89.34590,
                90.00620, 89.80260, 89.59910, 88.64890, 87.69870,
                85.49360, 83.28860, 83.49390, 83.69920, 81.86300,
                80.02680, 80.12070, 80.21460, 81.24620, 82.27780,
                80.28100, 78.28420, 74.00270, 69.72130, 70.66520,
                71.60910, 72.97900, 74.34900, 67.97650, 61.60400,
                65.74480, 69.88560, 72.48630, 75.08700, 69.33980,
                63.59270, 55.00540, 46.41820, 56.61180, 66.80540,
                65.09410, 63.38280, 63.84340, 64.30400, 61.87790,
                59.45190, 55.70540, 51.95900, 54.69980, 57.44060,
                58.87650, 60.31250
        }
    },
    {
        .name = "D75",
        .exif = lsD75,
        .band_start = 300,
        .band_end = 780,
        .band_spacing = 5,
        .scale = 100,
        .data = {
            0.043,2.588,5.133,17.470,29.808,42.369,54.930,56.095,57.259,60.000,62.740,62.861,62.982,66.647,70.312,68.507,66.703,68.333,69.963,85.946,101.929,106.911,111.894,112.346,112.798,107.945,103.092,112.145,121.198,127.104,133.010,132.682,132.355,129.838,127.322,127.061,126.800,122.291,117.783,117.186,116.589,115.146,113.702,111.181,108.659,109.552,110.445,108.367,106.289,105.596,104.904,102.452,100.000,97.808,95.616,94.914,94.213,90.605,86.997,87.112,87.227,86.684,86.140,84.861,83.581,81.164,78.747,78.587,78.428,76.614,74.801,74.562,74.324,74.873,75.422,73.499,71.576,67.714,63.852,64.464,65.076,66.573,68.070,62.256,56.443,60.343,64.242,66.697,69.151,63.890,58.629,50.623,42.617,51.985,61.352,59.838,58.324
        }
    },
    {
        .name = "StdE",
        .exif = -1,
        .band_start = 300,
        .band_end = 830,
        .band_spacing = 5,
        .scale = 100,
        .data = {
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100
        }
    },
    { .name = NULL,
      .exif = -1 }
};

static const spectrum_t *
spectraldb_illuminant_byidx(int i)

{
    static spectrum_t *ill[sizeof(ill_list)/sizeof(ill_list[0])];

    if (i < 0 || i >= (int)(sizeof(ill)/sizeof(ill[0]))) {
        return NULL;
    }
    if (ill[i] != NULL) {
        return ill[i];
    }
    ill[i] = spec_alloc(ill_list[i].data, ill_list[i].band_start, ill_list[i].band_end, ill_list[i].band_spacing);
    if (ill_list[i].scale != 1.0) {
        spec_scalar_multiply(ill[i], 1.0 / ill_list[i].scale);
    }
    return ill[i];
}

const spectrum_t *
spectraldb_illuminant(enum exif_lightsource ls)

{
    for (int i = 0; i < (int)(sizeof(ill_list)/sizeof(ill_list[0])); i++) {
        if (ill_list[i].exif == ls) {
            return spectraldb_illuminant_byidx(i);
        }
    }
    return NULL;
}

const spectrum_t *
spectraldb_illuminant_byname(const char name[],
                             bool missing_is_fatal)
{
    assert(name != NULL);
    double temp = atof(name);
    if ((strchr(name, 'K') != NULL || strchr(name, 'k') != NULL) && temp > 0) {
        // blackbody
        static struct {
            spectrum_t *s;
            double temp;
        } bb[100];
        int i = 0;
        for (i = 0; i < 100; i++) {
            if (bb[i].temp == 0) break;
            if (bb[i].temp == temp) return bb[i].s;
        }
        if (i == 100) {
            elog("Too many blackbodies\n");
            abort();
        }
        bb[i].s = spectraldb_render_blackbody(temp, 300, 830, 5);
        bb[i].temp = temp;
        return bb[i].s;
    }
    int i;
    for (i = 0; ill_list[i].name != NULL; i++) {
        if (strcasecmp(name, ill_list[i].name) == 0) {
            break;
        }
    }
    if (ill_list[i].name == NULL) {
        if (missing_is_fatal) {
            elog("Unknown spectral illuminant name \"%s\".\n", name);
            elog("  Known spectral illuminant names are: \"%s\"", ill_list[0].name);
            for (int i = 1; ill_list[i].name != NULL; i++) {
                elog(", \"%s\"", ill_list[i].name);
            }
            elog(", and blackbodies \"<temp>K\", eg \"2854K\".\n");
            exit(EXIT_FAILURE);
        }
        return NULL;
    }
    return spectraldb_illuminant_byidx(i);
}

const char *
spectraldb_illuminant_list(void)
{
    static char *list = NULL;
    if (list != NULL) {
        return list;
    }
    int len = 0;
    for (int i = 0; ill_list[i].name != NULL; i++) len += strlen(ill_list[i].name) + 3;
    char *str = malloc(len);
    str[0] = '\0';
    int offset = 0;
    for (int i = 0; ill_list[i].name != NULL; i++) {
        offset += sprintf(&str[offset], "%s, ", ill_list[i].name);
    }
    str[offset-2] = '\0';
    list = str;
    return str;
}

/*
static void __attribute__((constructor))
get_brightest_munsell(void)
{
    for (int i = 0; i < glob_cg[CG_MUNSELL].count; i++) {
        char n1[8];
        strcpy(n1, munsell_glossy_names[i]);
        n1[3] = n1[4] = ' ';
        bool brighter_exist = false;
        int v1 = (munsell_glossy_names[i][3]-'0')*10+(munsell_glossy_names[i][4]-'0');
        for (int j = 0; j < glob_cg[CG_MUNSELL].count; j++) {
            if (i == j) continue;
            char n2[8];
            strcpy(n2, munsell_glossy_names[j]);
            n2[3] = n2[4] = ' ';
            if (strcmp(n1, n2) == 0) {
                int v2 = (munsell_glossy_names[j][3]-'0')*10+(munsell_glossy_names[j][4]-'0');
                if (v2 > v1) {
                    brighter_exist = true;
                }
            }
        }
        if (!brighter_exist) {
            fprintf(stdout, "\"%s\", ", munsell_glossy_names[i]);
        }
    }
    fprintf(stdout, "\"NEUT950\"\n");
    exit(1);
}
*/

 /*
#include <colmath.h>
#include <observers.h>
static void __attribute__((constructor))
get_brightest_munsell(void)
{
    const struct observer *obs = observer_get(OBSERVER_1931_2);
    const spectrum_t *ill = spectraldb_illuminant(lsD50);
    for (int i = 0; i < glob_cg[CG_MUNSELL].count; i++) {
        char n1[8];
        strcpy(n1, munsell_glossy_names[i]);
        n1[3] = n1[4] = ' ';
        bool brighter_exist = false;
        int v1 = (munsell_glossy_names[i][3]-'0')*10+(munsell_glossy_names[i][4]-'0');
        for (int j = 0; j < glob_cg[CG_MUNSELL].count; j++) {
            if (i == j) continue;
            char n2[8];
            strcpy(n2, munsell_glossy_names[j]);
            n2[3] = n2[4] = ' ';
            if (strcmp(n1, n2) == 0) {
                int v2 = (munsell_glossy_names[j][3]-'0')*10+(munsell_glossy_names[j][4]-'0');
                if (v2 > v1) {
                    brighter_exist = true;
                }
            }
        }
        if (!brighter_exist && strstr(munsell_glossy_names[i], "NEUT") == NULL) {
            const spectrum_t *s = spectraldb_color("munsell", munsell_glossy_names[i]);
            v3 xyz = v3_norm(spec2xyz_ill(s, obs, ill, 1));
            v3 xyY = xyz2xyY(xyz);
            fprintf(stdout, "{ \"%c%c\", %f, %f },", munsell_glossy_names[i][1], munsell_glossy_names[i][2], xyY.v[0], xyY.v[1]);
        }
    }
    exit(1);
}
 */
