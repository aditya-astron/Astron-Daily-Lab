/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 */
#include "secded.h"

#include <math.h>

uint8_t  secded_col[SECDED_N];
uint64_t secded_row[SECDED_CHECK_BITS];
static int8_t syn_to_pos[256];

static int popcount8(unsigned v) { return __builtin_popcount(v & 0xFFu); }

/*
 * Hsiao construction: every column has odd weight, so any single error gives
 * an odd syndrome and any double error an even, nonzero one (columns are
 * distinct). That is what separates "correct" from "detect" without the
 * overall parity bit of an extended Hamming code.
 *
 * 64 data columns = all C(8,3) = 56 weight-3 columns + 8 weight-5 columns.
 * The weight-3 set already loads each row with 21 ones; the 8 cyclic
 * rotations of 00011111 add exactly 5 more to every row, so all rows carry
 * 26 data ones. Equal row weight means equal XOR-tree depth per check bit,
 * which is the point of Hsiao's choice in hardware.
 */
void secded_init(void)
{
    int j = 0;
    for (unsigned v = 0; v < 256; v++)
        if (popcount8(v) == 3) secded_col[j++] = (uint8_t)v;
    for (int r = 0; r < 8; r++) {
        unsigned v = 0x1Fu;
        secded_col[j++] = (uint8_t)(((v << r) | (v >> (8 - r))) & 0xFFu);
    }
    for (int i = 0; i < SECDED_CHECK_BITS; i++) secded_col[64 + i] = (uint8_t)(1u << i);

    for (int i = 0; i < SECDED_CHECK_BITS; i++) {
        secded_row[i] = 0;
        for (int k = 0; k < 64; k++)
            if (secded_col[k] >> i & 1u) secded_row[i] |= 1ULL << k;
    }
    for (int s = 0; s < 256; s++) syn_to_pos[s] = -1;
    for (int k = 0; k < SECDED_N; k++) syn_to_pos[secded_col[k]] = (int8_t)k;
}

uint8_t secded_encode(uint64_t data)
{
    uint8_t c = 0;
    for (int i = 0; i < SECDED_CHECK_BITS; i++)
        c |= (uint8_t)((__builtin_popcountll(data & secded_row[i]) & 1) << i);
    return c;
}

int secded_locate(uint8_t syndrome) { return syn_to_pos[syndrome]; }

secded_status secded_decode(uint64_t *data, uint8_t *check, int *pos)
{
    uint8_t s = secded_encode(*data) ^ *check;
    *pos = -1;
    if (s == 0) return SECDED_CLEAN;
    if ((popcount8(s) & 1) == 0) return SECDED_DUE;  /* even: >= 2 errors */
    int k = syn_to_pos[s];
    if (k < 0) return SECDED_DUE;                    /* odd but not a column */
    if (k < 64) *data ^= 1ULL << k;
    else        *check ^= (uint8_t)(1u << (k - 64));
    *pos = k;
    return SECDED_CORRECTED;
}

/* Probability that a Poisson count with mean m is odd: (1 - e^{-2m}) / 2. */
static double odd_prob(double m) { return -0.5 * expm1(-2.0 * m); }

/*
 * A bit's final state after T is the XOR of the parities of every event
 * source that touches it, so only the odd-count probability of each source
 * matters (two hits on the same cell restore it).
 *
 * interleave >= 2: the neighbour of a cell always belongs to another word, so
 * each bit is driven by independent sources only: its own SBUs and MCUs and
 * the MCU spilling from its physical predecessor. Bits are i.i.d. with
 * per-bit event rate lam*(1+f), giving a binomial on 72 bits.
 *
 * interleave == 1: an MCU starting on bit c flips bits c and c+1 of the same
 * word, so bits are correlated. Walk the 72 bits carrying the parity of the
 * MCU source that straddles the boundary; state = (carry, flips so far
 * capped at 2). Bit 0 receives the spill from the previous word, and the
 * spill out of bit 71 lands in the next word, so for this word it acts as a
 * single-bit source.
 */
double word_fail_prob(double lam, double f, int interleave, double T,
                      double *p0, double *p1)
{
    double P[2][3] = {{0}};
    if (interleave >= 2) {
        double q = odd_prob(lam * (1.0 + f) * T);
        double a = 1.0 - q;
        double b0 = pow(a, 72), b1 = 72.0 * q * pow(a, 71), tail = 0.0;
        double term = b1;
        for (int k = 2; k <= 72; k++) {           /* C(72,k) q^k a^(72-k) */
            term *= (double)(72 - k + 1) / k * q / a;
            tail += term;
        }
        if (p0) *p0 = b0;
        if (p1) *p1 = b1;
        return tail;
    }
    double ps = odd_prob(lam * (1.0 - f) * T); /* single-cell events  */
    double pm = odd_prob(lam * f * T);         /* multi-cell events   */
    P[0][0] = 1.0 - pm;                        /* spill into bit 0    */
    P[1][0] = pm;
    for (int c = 0; c < 72; c++) {
        double N[2][3] = {{0}};
        for (int cin = 0; cin < 2; cin++)
            for (int n = 0; n < 3; n++) {
                double w = P[cin][n];
                if (w == 0.0) continue;
                for (int s = 0; s < 2; s++)
                    for (int m = 0; m < 2; m++) {
                        double pr = w * (s ? ps : 1.0 - ps) * (m ? pm : 1.0 - pm);
                        int flip = cin ^ s ^ m;
                        int nn = n + flip > 2 ? 2 : n + flip;
                        int cout = (c == 71) ? 0 : m;
                        N[cout][nn] += pr;
                    }
            }
        for (int a = 0; a < 2; a++)
            for (int n = 0; n < 3; n++) P[a][n] = N[a][n];
    }
    if (p0) *p0 = P[0][0];
    if (p1) *p1 = P[0][1];
    return P[0][2];
}
