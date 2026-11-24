/*
 * Astron Built this AP
 * Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-24 IST
 */
#ifndef CIC_H
#define CIC_H

#include <stddef.h>
#include <stdint.h>

#define CIC_MAX_N 6
#define CIC_MAX_M 4

/*
 * Hogenauer CIC decimator: N integrators at the input rate, decimate by R,
 * N combs with differential delay M. All registers are unsigned and wrap
 * modulo 2^width; two's complement wraparound is harmless as long as the
 * last stage is wide enough for the true output range (Hogenauer 1981).
 *
 * Stage j = 1..2N (integrators 1..N, combs N+1..2N), j = 2N+1 is the output
 * register. prune[j] is the number of LSBs (of full precision) discarded at
 * the input of stage j; it must be non-decreasing in j. Register j holds
 * bits [prune[j], width_full).
 */
typedef struct {
    int N, R, M;
    int width_full;                 /* Bmax: bits of the unpruned output   */
    int prune[2 * CIC_MAX_N + 2];   /* index 0 = input (always 0)          */
    uint64_t mask[2 * CIC_MAX_N + 2];
    uint64_t integ[CIC_MAX_N];
    uint64_t comb_dly[CIC_MAX_N][CIC_MAX_M];
    int dly_pos;
    int phase;
} cic_t;

/* Full-precision width for a Bin-bit two's complement input. */
int cic_full_width(int N, int R, int M, int bin);

/*
 * width: register width of the unpruned design (normally cic_full_width).
 * prune: array of 2N+2 entries or NULL for no pruning.
 * Returns 0 on success, -1 on bad parameters.
 */
int cic_init(cic_t *c, int N, int R, int M, int width, const int *prune);

/*
 * Push one input sample. Returns 1 and writes *out when a decimated output
 * is produced. *out is a signed integer in units of 2^prune[2N+1].
 */
int cic_push(cic_t *c, int32_t x, int64_t *out);

/* Run a whole block; returns the number of outputs written. */
size_t cic_run(cic_t *c, const int32_t *x, size_t n, int64_t *out);

/* Same, for a +/-1 bitstream stored as int8. */
size_t cic_run_i8(cic_t *c, const int8_t *x, size_t n, int64_t *out);

/*
 * Hogenauer's variance gain F_j^2 from the input of stage j to the output
 * (sum of squared impulse-response coefficients), j = 1..2N+1.
 */
double cic_hogenauer_F2(int N, int R, int M, int j);

/*
 * Hogenauer pruning for an output word of bout bits:
 *   B_j = floor(-log2 F_j + log2 sigma_T + 0.5 log2(6/N)),
 *   sigma_T^2 = 2^(2 B_{2N+1}) / 12,  B_{2N+1} = width - bout.
 * Writes 2N+2 entries (prune[0] = 0), clipped to [0, B_{2N+1}] and made
 * non-decreasing.
 */
void cic_hogenauer_prune(int N, int R, int M, int width, int bout, int *prune);

/*
 * Second-order 1-bit sigma-delta modulator with NTF = (1 - z^-1)^2 and
 * STF = z^-2:
 *   y[n]    = sign(s2[n])
 *   s1[n+1] = s1[n] + x[n] - y[n]
 *   s2[n+1] = s2[n] + s1[n] - 2 y[n]
 * e[n] = y[n] - s2[n] is returned (if e != NULL) so the exact linear
 * decomposition Y = z^-2 X + (1 - z^-1)^2 E can be checked.
 */
void sdm2_run(const double *x, size_t n, int8_t *y, double *e);

#endif
