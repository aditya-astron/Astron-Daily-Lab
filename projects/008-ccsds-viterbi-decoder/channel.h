/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
#ifndef CHANNEL_H
#define CHANNEL_H

#include <stddef.h>
#include <stdint.h>

enum rx_mode { RX_HARD, RX_Q3, RX_SOFT };

/* xoshiro256** so runs are reproducible across platforms. */
typedef struct { uint64_t s[4]; } rng_t;

void   rng_seed(rng_t *r, uint64_t seed);
uint64_t rng_next(rng_t *r);
double rng_gauss(rng_t *r);

/* BPSK over AWGN at the given Eb/N0 (dB) for a code of rate R, followed by
 * the receiver's quantizer. sym: 0/1 code symbols, out: metric inputs. */
void bpsk_awgn(const uint8_t *sym, size_t m, double ebn0_db, double rate,
               enum rx_mode mode, rng_t *r, double *out);

double q_func(double x);

#endif
