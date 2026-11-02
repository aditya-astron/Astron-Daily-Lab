/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
#include "viterbi.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

/* Register layout: reg = (u << 6) | state, where bit 6 is the current input
 * and bit 5 the previous one. The octal generator's MSB taps the current
 * input, so a generator applies directly as a mask on reg. */
static int parity7(unsigned x)
{
    x ^= x >> 4;
    x ^= x >> 2;
    x ^= x >> 1;
    return (int)(x & 1u);
}

static void branch_out(unsigned state, unsigned u, int *c1, int *c2)
{
    unsigned reg = (u << CC_MEM) | state;
    *c1 = parity7(reg & CC_G1);
    *c2 = parity7(reg & CC_G2) ^ 1;   /* CCSDS symbol inversion on G2 */
}

void cc_encode(const uint8_t *bits, size_t n, uint8_t *out)
{
    unsigned state = 0;
    for (size_t i = 0; i < n + CC_MEM; i++) {
        unsigned u = i < n ? (bits[i] & 1u) : 0u;
        int c1, c2;
        branch_out(state, u, &c1, &c2);
        out[2 * i] = (uint8_t)c1;
        out[2 * i + 1] = (uint8_t)c2;
        state = ((u << CC_MEM) | state) >> 1;
    }
}

int cc_viterbi_decode(const double *rx, size_t n, uint8_t *bits_out)
{
    size_t steps = n + CC_MEM;
    /* One survivor decision bit per state per step, packed into 64 bits. */
    uint64_t *dec = malloc(steps * sizeof *dec);
    if (!dec)
        return -1;

    /* Expected BPSK symbols for every (state, input) branch. */
    double e1[CC_STATES][2], e2[CC_STATES][2];
    for (unsigned s = 0; s < CC_STATES; s++)
        for (unsigned u = 0; u < 2; u++) {
            int c1, c2;
            branch_out(s, u, &c1, &c2);
            e1[s][u] = 1.0 - 2.0 * c1;
            e2[s][u] = 1.0 - 2.0 * c2;
        }

    /* Path metric = correlation sum(r * e), maximised. For AWGN this is the
     * log-likelihood up to constants; for +/-1 inputs it is Hamming distance. */
    double pm[CC_STATES], nm[CC_STATES];
    for (int s = 0; s < CC_STATES; s++)
        pm[s] = -DBL_MAX / 4;
    pm[0] = 0.0;   /* encoder starts in the zero state */

    for (size_t t = 0; t < steps; t++) {
        double r1 = rx[2 * t], r2 = rx[2 * t + 1];
        uint64_t d = 0;
        for (unsigned ns = 0; ns < CC_STATES; ns++) {
            /* ns = (u << 5) | (s >> 1): two predecessors differ in bit 0. */
            unsigned u = ns >> (CC_MEM - 1);
            unsigned s0 = (ns << 1) & (CC_STATES - 1);
            unsigned s1 = s0 | 1u;
            double m0 = pm[s0] + r1 * e1[s0][u] + r2 * e2[s0][u];
            double m1 = pm[s1] + r1 * e1[s1][u] + r2 * e2[s1][u];
            if (m1 > m0) {
                nm[ns] = m1;
                d |= 1ull << ns;
            } else {
                nm[ns] = m0;
            }
        }
        dec[t] = d;
        memcpy(pm, nm, sizeof pm);
    }

    /* Tail bits force the encoder back to state 0, so trace back from it. */
    unsigned state = 0;
    for (size_t t = steps; t-- > 0;) {
        unsigned u = state >> (CC_MEM - 1);
        if (t < n)
            bits_out[t] = (uint8_t)u;
        unsigned b = (unsigned)((dec[t] >> state) & 1u);
        state = ((state << 1) & (CC_STATES - 1)) | b;
    }
    free(dec);
    return 0;
}

int cc_distance_spectrum(int dmax, double *a, double *b)
{
    /* Enumerate every path that leaves state 0 and first returns to it.
     * live_n[s][w]: number of partial paths now in state s with output
     * weight w; live_i: their summed input weight. The code is not
     * catastrophic (no zero-weight cycle), so every surviving path gains
     * weight eventually and the loop terminates once all exceed dmax. */
    size_t W = (size_t)dmax + 1;
    double *ln = calloc(CC_STATES * W, sizeof *ln);
    double *li = calloc(CC_STATES * W, sizeof *li);
    double *tn = calloc(CC_STATES * W, sizeof *tn);
    double *ti = calloc(CC_STATES * W, sizeof *ti);
    if (!ln || !li || !tn || !ti) {
        free(ln); free(li); free(tn); free(ti);
        return -1;
    }
    for (int d = 0; d <= dmax; d++)
        a[d] = b[d] = 0.0;

    /* Distances are those of the linear code, so use uninverted outputs. */
#define RAW_W(s, u) (parity7((((u) << CC_MEM) | (s)) & CC_G1) + \
                     parity7((((u) << CC_MEM) | (s)) & CC_G2))

    /* Diverge: input 1 from state 0. */
    int w0 = RAW_W(0u, 1u);
    if (w0 <= dmax) {
        ln[(1u << (CC_MEM - 1)) * W + (size_t)w0] = 1.0;
        li[(1u << (CC_MEM - 1)) * W + (size_t)w0] = 1.0;
    }

    for (;;) {
        int alive = 0;
        memset(tn, 0, CC_STATES * W * sizeof *tn);
        memset(ti, 0, CC_STATES * W * sizeof *ti);
        for (unsigned s = 1; s < CC_STATES; s++)
            for (size_t w = 0; w < W; w++) {
                double cn = ln[s * W + w], ci = li[s * W + w];
                if (cn == 0.0)
                    continue;
                for (unsigned u = 0; u < 2; u++) {
                    size_t nw = w + (size_t)RAW_W(s, u);
                    unsigned ns = ((u << CC_MEM) | s) >> 1;
                    if (nw >= W)
                        continue;
                    double add_i = ci + cn * u;
                    if (ns == 0) {
                        a[nw] += cn;
                        b[nw] += add_i;
                    } else {
                        tn[ns * W + nw] += cn;
                        ti[ns * W + nw] += add_i;
                        alive = 1;
                    }
                }
            }
        double *x = ln; ln = tn; tn = x;
        x = li; li = ti; ti = x;
        if (!alive)
            break;
    }
#undef RAW_W
    free(ln); free(li); free(tn); free(ti);

    for (int d = 0; d <= dmax; d++)
        if (a[d] > 0.0)
            return d;
    return -1;
}
