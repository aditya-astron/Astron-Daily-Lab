/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * Monte Carlo of a physical memory array between two scrub passes: upsets
 * land on physical cells, the real decoder runs on every word, and the
 * outcome is compared with the stored truth.
 */
#ifndef SCRUB_SIM_H
#define SCRUB_SIM_H

#include <stdint.h>

typedef struct {
    long words;      /* word-intervals simulated                      */
    long clean;      /* zero syndrome, data intact                    */
    long corrected;  /* single error corrected back to the truth      */
    long due;        /* detected uncorrectable                        */
    long sdc;        /* decoder returned wrong data without flagging  */
} scrub_tally;

uint64_t rng_next(uint64_t *s);
double   rng_uniform(uint64_t *s);

/* n_words must be a multiple of interleave. lamT = event count mean per
 * cell over one scrub interval; f = multi-cell (2-cell) fraction. */
void scrub_mc(int n_words, int interleave, double lamT, double f,
              int intervals, uint64_t seed, scrub_tally *t);

#endif
