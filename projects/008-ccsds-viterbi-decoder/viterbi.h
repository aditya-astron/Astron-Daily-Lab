/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
#ifndef VITERBI_H
#define VITERBI_H

#include <stddef.h>
#include <stdint.h>

/* CCSDS 131.0-B basic convolutional code: K = 7, R = 1/2,
 * G1 = 171 (octal), G2 = 133 (octal), second symbol inverted on the wire. */
#define CC_K        7
#define CC_MEM      (CC_K - 1)
#define CC_STATES   (1 << CC_MEM)
#define CC_G1       0171
#define CC_G2       0133

/* Encode n info bits (0/1) plus CC_MEM zero tail bits.
 * Writes 2 * (n + CC_MEM) code symbols (0/1) to out. */
void cc_encode(const uint8_t *bits, size_t n, uint8_t *out);

/* Maximum-likelihood decode of a tail-terminated block.
 * rx holds 2 * (n + CC_MEM) real-valued soft symbols where a positive value
 * means "code bit 0 more likely" (BPSK map 0 -> +1, 1 -> -1). Hard or
 * quantized decisions are just a special case of the input values.
 * Returns 0 on success, -1 on allocation failure. */
int cc_viterbi_decode(const double *rx, size_t n, uint8_t *bits_out);

/* Distance spectrum of the code: for output weights d = 0..dmax, a[d] is the
 * number of error events of weight d and b[d] the total number of
 * information-bit errors in them. Returns d_free. */
int cc_distance_spectrum(int dmax, double *a, double *b);

#endif
