/*
 * Astron Built this AP
 * Project 023 | Error correction | CCSDS Reed-Solomon (255,223) codec with errors-and-erasures decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-08 IST
 *
 * Link-level view of the code:
 *  1. Monte Carlo codeword error rate against the bounded-distance formula
 *       P_W = sum_{i=17}^{255} C(255,i) ps^i (1-ps)^(255-i)
 *     for independent symbol errors (what an ideal interleaver delivers).
 *  2. Required Eb/N0 for a decoded bit error rate of 1e-5 on hard-decision
 *     BPSK, and the coding gain over uncoded BPSK.
 *  3. Interleaving depth vs. the longest burst a frame survives.
 * Exits nonzero if the simulation disagrees with the formula by > 4 sigma.
 */
#include "rs.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static uint64_t rng = 0xD1B54A32D192ED03ull;
static uint64_t rnd64(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return rng;
}
static double urand(void) { return (double)(rnd64() >> 11) * 0x1.0p-53; }

/* log C(n, i) via lgamma keeps the tail sum stable at small ps */
static double binom_term(int n, int i, double p)
{
    return exp(lgamma(n + 1.0) - lgamma(i + 1.0) - lgamma(n - i + 1.0)
               + i * log(p) + (n - i) * log1p(-p));
}

static double word_error(double ps)
{
    double s = 0;
    for (int i = RS_T + 1; i <= RS_N; i++)
        s += binom_term(RS_N, i, ps);
    return s;
}

/* Decoded bit error rate. A failed word passes its i symbol errors through
 * uncorrected; a wrong 8-bit symbol has 128/255 of its bits wrong on average. */
static double decoded_ber(double ps)
{
    double s = 0;
    for (int i = RS_T + 1; i <= RS_N; i++)
        s += (double)i / RS_N * binom_term(RS_N, i, ps);
    return s * 128.0 / 255.0;
}

static double qfunc(double x) { return 0.5 * erfc(x / sqrt(2.0)); }

/* hard-decision BPSK at a given information Eb/N0 (dB), code rate r */
static double symbol_error(double ebn0_db, double rate)
{
    double pb = qfunc(sqrt(2.0 * rate * pow(10.0, ebn0_db / 10.0)));
    return 1.0 - pow(1.0 - pb, 8);
}

static double solve_ebn0(double target, double rate, int coded)
{
    double lo = 0.0, hi = 15.0;
    for (int it = 0; it < 100; it++) {
        double mid = 0.5 * (lo + hi);
        double pe = coded ? decoded_ber(symbol_error(mid, rate))
                          : qfunc(sqrt(2.0 * pow(10.0, mid / 10.0)));
        if (pe > target)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

int main(void)
{
    rs_init();
    int bad = 0;

    printf("CCSDS RS(255,223), t = %d, rate = %.4f\n", RS_T, (double)RS_K / RS_N);
    printf("g(x) coefficients, x^32 first (conventional basis):\n ");
    for (int i = 0; i <= RS_NROOTS; i++)
        printf(" %02x", rs_genpoly[i]);
    printf("\n\n");

    /* 1. Monte Carlo vs. bounded-distance formula */
    printf("Independent symbol errors, 20000 words per point\n");
    printf("   ps     P_W formula   P_W sim     fails  miscorr   z\n");
    const double pss[] = {0.030, 0.035, 0.040, 0.045, 0.050};
    const int frames = 20000;
    uint8_t cw[RS_N], rx[RS_N];
    for (unsigned q = 0; q < sizeof pss / sizeof pss[0]; q++) {
        double ps = pss[q];
        int fail = 0, miscorr = 0;
        for (int f = 0; f < frames; f++) {
            for (int i = 0; i < RS_K; i++)
                cw[i] = (uint8_t)rnd64();
            rs_encode_conv(cw);
            memcpy(rx, cw, RS_N);
            for (int i = 0; i < RS_N; i++)
                if (urand() < ps)
                    rx[i] ^= (uint8_t)(1 + rnd64() % 255);
            int r = rs_decode_conv(rx, NULL, 0, 0);
            if (r < 0)
                fail++;
            else if (memcmp(rx, cw, RS_N) != 0)
                miscorr++;
        }
        double pw = word_error(ps);
        double sim = (double)(fail + miscorr) / frames;
        double z = (sim - pw) / sqrt(pw * (1 - pw) / frames);
        printf("  %.3f  %.4e    %.4e  %6d  %6d   %+5.2f\n", ps, pw, sim, fail, miscorr, z);
        if (fabs(z) > 4.0 || miscorr > 0)
            bad = 1;
    }

    /* 2. Coding gain on hard-decision BPSK */
    double rate = (double)RS_K / RS_N;
    double eb_unc = solve_ebn0(1e-5, 1.0, 0);
    double eb_rs = solve_ebn0(1e-5, rate, 1);
    printf("\nHard-decision BPSK, ideal symbol interleaving, target BER 1e-5\n");
    printf("  uncoded Eb/N0      %.2f dB\n", eb_unc);
    printf("  RS(255,223) Eb/N0  %.2f dB   (channel ps = %.3e)\n", eb_rs, symbol_error(eb_rs, rate));
    printf("  coding gain        %.2f dB\n", eb_unc - eb_rs);
    if (fabs(eb_unc - 9.59) > 0.01)   /* textbook BPSK value */
        bad = 1;

    /* 3. Burst tolerance vs. interleaving depth */
    printf("\nInterleaving depth vs. longest correctable channel burst\n");
    printf("  I   frame bytes  burst bytes  burst bits (guaranteed)\n");
    for (int depth = 1; depth <= 8; depth++) {
        if (depth == 6 || depth == 7)
            continue;   /* CCSDS allows I = 1..5 and 8 */
        /* a burst of b bits touches at most ceil((b-1)/8)+1 bytes */
        int bytes = RS_T * depth, bits = (bytes - 1) * 8 + 1;
        printf("  %d   %5d        %4d         %5d\n", depth, RS_N * depth, bytes, bits);
    }

    printf("\n%s\n", bad ? "MISMATCH" : "simulation agrees with theory");
    return bad;
}
