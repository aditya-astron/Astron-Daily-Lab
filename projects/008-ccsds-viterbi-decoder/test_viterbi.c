/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
#include "channel.h"
#include "viterbi.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;

#define CHECK(cond, ...) do { \
    if (cond) { printf("PASS  "); } else { printf("FAIL  "); failures++; } \
    printf(__VA_ARGS__); printf("\n"); } while (0)

static void test_impulse_response(void)
{
    /* A lone 1 walks through the register, so the output stream spells the
     * generator taps MSB first; the G2 stream is inverted per CCSDS. */
    uint8_t in[1] = { 1 }, out[2 * (1 + CC_MEM)];
    cc_encode(in, 1, out);
    int ok = 1;
    for (int k = 0; k < CC_K; k++) {
        ok &= out[2 * k] == ((CC_G1 >> (CC_MEM - k)) & 1);
        ok &= out[2 * k + 1] == (((CC_G2 >> (CC_MEM - k)) & 1) ^ 1);
    }
    CHECK(ok, "impulse response reproduces G1=171, G2=133 (inverted)");
}

static void test_spectrum(void)
{
    /* Reference: Proakis, Digital Communications, Table 8.2-x; Odenwalder. */
    const int    d[]  = { 10, 12, 14, 16, 18, 20 };
    const double ad[] = { 11, 38, 193, 1331, 7275, 40406 };
    const double bd[] = { 36, 211, 1404, 11633, 77433, 502690 };
    double a[21], b[21];
    int dfree = cc_distance_spectrum(20, a, b);
    CHECK(dfree == 10, "free distance = %d (expected 10)", dfree);
    int ok = 1;
    for (int i = 0; i < 6; i++) {
        ok &= a[d[i]] == ad[i] && b[d[i]] == bd[i];
        ok &= a[d[i] - 1] == 0.0;   /* only even weights for this code */
    }
    CHECK(ok, "weight spectrum a_d, B_d matches published values for d=10..20");
}

static void test_noiseless(void)
{
    enum { N = 5000 };
    static uint8_t info[N], code[2 * (N + CC_MEM)], dec[N];
    static double rx[2 * (N + CC_MEM)];
    rng_t r;
    rng_seed(&r, 1);
    for (int i = 0; i < N; i++)
        info[i] = (uint8_t)(rng_next(&r) >> 63);
    cc_encode(info, N, code);
    for (int i = 0; i < 2 * (N + CC_MEM); i++)
        rx[i] = 1.0 - 2.0 * code[i];
    cc_viterbi_decode(rx, N, dec);
    CHECK(memcmp(info, dec, N) == 0, "noiseless %d-bit block decodes exactly", N);
}

static void test_guaranteed_correction(void)
{
    /* With d_free = 10 and a terminated trellis, any pattern of up to 4 hard
     * symbol errors leaves the sent codeword strictly closest, so ML decoding
     * must recover it. Exhaustive over random placements. */
    enum { N = 64, M = 2 * (N + CC_MEM), TRIALS = 20000 };
    uint8_t info[N], code[M], dec[N];
    double rx[M];
    rng_t r;
    rng_seed(&r, 7);
    int bad = 0;
    for (int t = 0; t < TRIALS; t++) {
        for (int i = 0; i < N; i++)
            info[i] = (uint8_t)(rng_next(&r) >> 63);
        cc_encode(info, N, code);
        for (int i = 0; i < M; i++)
            rx[i] = 1.0 - 2.0 * code[i];
        int nerr = 1 + (int)(rng_next(&r) % 4);
        for (int e = 0; e < nerr; e++) {
            /* clustered half the time: bursts are the harder case */
            int base = (int)(rng_next(&r) % M);
            int pos = (t & 1) ? (base % (M - 12)) + (int)(rng_next(&r) % 12) : base;
            rx[pos] = 1.0 - 2.0 * code[pos] == rx[pos] ? -rx[pos] : rx[pos];
        }
        cc_viterbi_decode(rx, N, dec);
        bad += memcmp(info, dec, N) != 0;
    }
    CHECK(bad == 0, "all %d blocks with <=4 symbol errors corrected (%d failures)",
          TRIALS, bad);

    /* And 5 errors placed on the minimum-weight error event can fool it:
     * flipping 5 of the 10 differing symbols makes a tie, 6 makes a loss. */
    uint8_t zero[N] = { 0 }, one[N] = { 0 }, cz[M], co[M];
    one[10] = 1;
    cc_encode(zero, N, cz);
    cc_encode(one, N, co);
    int diff[M], nd = 0;
    for (int i = 0; i < M; i++)
        if (cz[i] != co[i])
            diff[nd++] = i;
    for (int i = 0; i < M; i++)
        rx[i] = 1.0 - 2.0 * cz[i];
    for (int k = 0; k < 6; k++)
        rx[diff[k]] = -rx[diff[k]];
    cc_viterbi_decode(rx, N, dec);
    CHECK(nd == 10 && memcmp(dec, one, N) == 0,
          "single-bit event has weight %d; 6 errors on it decode to the wrong word", nd);
}

static double sim_ber(enum rx_mode mode, double snr, long nbits, uint64_t seed)
{
    enum { N = 4096 };
    static uint8_t info[N], code[2 * (N + CC_MEM)], dec[N];
    static double rx[2 * (N + CC_MEM)];
    rng_t r;
    rng_seed(&r, seed);
    long errs = 0, bits = 0;
    while (bits < nbits) {
        for (int i = 0; i < N; i++)
            info[i] = (uint8_t)(rng_next(&r) >> 63);
        cc_encode(info, N, code);
        bpsk_awgn(code, 2 * (N + CC_MEM), snr, 0.5, mode, &r, rx);
        cc_viterbi_decode(rx, N, dec);
        for (int i = 0; i < N; i++)
            errs += info[i] != dec[i];
        bits += N;
    }
    return (double)errs / bits;
}

static void test_channel_and_ber(void)
{
    /* Uncoded hard decisions must hit Q(sqrt(2 Eb/N0)) at 4 dB: 1.25e-2. */
    enum { M = 1000000 };
    static uint8_t sym[M];
    static double y[M];
    rng_t r;
    rng_seed(&r, 3);
    bpsk_awgn(sym, M, 4.0, 1.0, RX_HARD, &r, y);
    long e = 0;
    for (int i = 0; i < M; i++)
        e += y[i] < 0.0;
    double p = (double)e / M, pt = q_func(sqrt(2.0 * pow(10.0, 0.4)));
    CHECK(fabs(p - pt) < 4.0 * sqrt(pt / M),
          "uncoded BPSK at 4 dB: %.4e vs theory %.4e", p, pt);

    /* Coded, 4 dB: soft decision should sit under the union bound and within
     * a factor of 3 of it (the bound is tight here); quantization costs. */
    double a[25], b[25];
    cc_distance_spectrum(24, a, b);
    double ub = 0.0, ebn0 = pow(10.0, 0.4);
    for (int d = 0; d <= 24; d++)
        ub += b[d] * q_func(sqrt(d * ebn0));
    long nb = 2000000;
    double soft = sim_ber(RX_SOFT, 4.0, nb, 11);
    double q3 = sim_ber(RX_Q3, 4.0, nb, 11);
    double hard = sim_ber(RX_HARD, 4.0, nb, 11);
    CHECK(soft < ub && soft > ub / 3.0,
          "soft BER at 4 dB %.2e vs union bound %.2e", soft, ub);
    CHECK(soft <= q3 && q3 < hard,
          "ordering soft %.2e <= 3-bit %.2e < hard %.2e", soft, q3, hard);
    /* Hard decisions: binary symmetric channel with crossover p, exact
     * pairwise error P_d (ties split evenly), bound Pb <= sum B_d P_d. */
    double pc = q_func(sqrt(ebn0)), hub = 0.0;
    for (int d = 10; d <= 24; d += 2) {
        double pd = 0.0;
        for (int k = d / 2; k <= d; k++) {
            double c = 1.0;
            for (int j = 0; j < k; j++)
                c = c * (d - j) / (j + 1);
            double term = c * pow(pc, k) * pow(1.0 - pc, d - k);
            pd += k == d / 2 ? 0.5 * term : term;
        }
        hub += b[d] * pd;
    }
    CHECK(hard < hub && hard < pt,
          "hard BER %.2e under BSC union bound %.2e and uncoded %.2e", hard, hub, pt);
}

int main(void)
{
    test_impulse_response();
    test_spectrum();
    test_noiseless();
    test_guaranteed_correction();
    test_channel_and_ber();
    printf(failures ? "\n%d test(s) failed\n" : "\nall tests passed\n", failures);
    return failures != 0;
}
