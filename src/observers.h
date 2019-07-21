/*
 * (c) Copyright 2015 -- Anders Torger
 *
 * This program is open source. For license terms, see the LICENSE file.
 *
 */
#ifndef OBSERVERS_H_
#define OBSERVERS_H_

#include <stdbool.h>
#include <colmath.h>

enum observer_name {
    OBSERVER_2006_10 = 0,
    OBSERVER_2006_2,
    OBSERVER_1964_10,
    OBSERVER_1978_JUDDVOS_2,
    OBSERVER_1931_2,
#define OBSERVER_MAX OBSERVER_1931_2
};

const struct observer *
observer_get(enum observer_name obs);

const struct observer *
observer_get_byname(const char name[],
                    bool missing_is_fatal);

const char*
observer_list(void);

#endif
