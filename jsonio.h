/*
 * (c) Copyright 2015, 2017 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef JSONIO_H_
#define JSONIO_H_

#include <stdbool.h>
#include <inttypes.h>

#include <spectrum.h>
#include <colmath.h>
#include <profio.h>
#include <dcamprof.h>
#include <cJSON.h>

uint8_t *
fileio_file2mem(const char filename[],
                size_t *size);

char *
fileio_file2string(const char filename[]);

struct target_adjustment_config *
jsonio_target_adjustment_parse(const char filename[],
                               bool parse_error_is_fatal);

void
jsonio_target_adjustment_delete(struct target_adjustment_config *conf);

struct tone_rep_op_config *
jsonio_ntro_conf_parse(const char filename[],
                       bool parse_error_is_fatal);

void
jsonio_ntro_delete(struct tone_rep_op_config *ntro);

struct dcam_ssf *
jsonio_dcam_ssf_parse(const char filename[],
                      bool parse_error_is_fatal);

void
jsonio_dcam_ssf_delete(struct dcam_ssf *ssf);

spectrum_t *
jsonio_spectrum_parse(const char filename[],
                      bool parse_error_is_fatal);

void
jsonio_spectrum_print(FILE *stream,
                      const char indent[],
                      const char name[],
                      const spectrum_t *s);

spectrum_t *
jsonio_illuminant_parse(const char filename[],
                        bool parse_error_is_fatal);

bool
jsonio_transfer_function_parse(const char filename[],
                               double *trc[3],
                               int *count_,
                               bool error_is_fatal);

bool
jsonio_tonecurve_parse(const char filename[],
                       double **trc,
                       int *count_,
                       bool error_is_fatal);

m3x3 *
jsonio_3x3matrix_parse(const char filename[],
                       const char optional_matrix_name[],
                       const char optional_matrix_name2[],
                       bool parse_error_is_fatal);

void
jsonio_3x3matrix_print(FILE *stream,
                       const char indent[],
                       const char name[],
                       m3x3 m);

struct profio_dcp *
jsonio_dcp_parse(const char filename[],
                 cJSON *root,
                 bool parse_error_is_fatal);

struct dcam_profile *
jsonio_profile_parse(const char filename[],
                     bool parse_error_is_fatal);

struct profio_icc *
jsonio_icc_parse(const char filename[],
                 cJSON *root,
                 bool parse_error_is_fatal);

struct test_chart_spec *
jsonio_test_chart_spec_parse(const char filename[],
                             bool parse_error_is_fatal);

#endif
