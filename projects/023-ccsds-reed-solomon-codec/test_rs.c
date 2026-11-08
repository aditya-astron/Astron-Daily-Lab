/*
 * Astron Built this AP
 * Project 023 | Error correction | CCSDS Reed-Solomon (255,223) codec with errors-and-erasures decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-08 IST
 *
 * Structural checks on the field, generator and dual basis, then randomized
 * decode checks right up to the 2e + s <= 32 boundary and one step past it.
 */
#include "rs.h"

#include <stdio.h>
#include <string.h>

static int fails = 0, checks = 0;

#define CHECK(cond, ...)                                       \
    do {                                                       \
        checks++;                                              \
        if (!(cond)) {                                         \
            fails++;                                           \
            printf("FAIL %s:%d: ", __FILE__, __LINE__);        \
            printf(__VA_ARGS__);                               \
            printf("\n");                                      \
        }                                                      \
    } while (0)

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static uint32_t rnd(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (uint32_t)(rng >> 16);
}
static uint8_t rnd_nonzero(void) { return (uint8_t)(1 + rnd() % 255); }

/* k distinct indices in [lo, hi) */
static void pick(int *idx, int k, int lo, int hi)
{
    for (int a = 0; a < k; a++) {
        int dup;
        do {
            idx[a] = lo + (int)(rnd() % (unsigned)(hi - lo));
            dup = 0;
            for (int b = 0; b < a; b++)
                dup |= idx[a] == idx[b];
        } while (dup);
    }
}

static void random_codeword(uint8_t cw[RS_N])
{
    for (int i = 0; i < RS_K; i++)
        cw[i] = (uint8_t)rnd();
    rs_encode_conv(cw);
}

static uint8_t gpoly_at(uint8_t x)
{
    uint8_t acc = 0;
    for (int i = 0; i <= RS_NROOTS; i++)
        acc = gf_mul(acc, x) ^ rs_genpoly[i];
    return acc;
}

static void test_field_and_generator(void)
{
    int seen[256] = {0};
    for (int i = 0; i < 255; i++)
        seen[gf_exp[i]]++;
    int ok = !seen[0];
    for (int v = 1; v < 256; v++)
        ok &= seen[v] == 1;
    CHECK(ok, "alpha is not primitive for x^8+x^7+x^2+x+1");

    /* beta = alpha^11 generates the whole group since gcd(11, 255) = 1 */
    uint8_t beta = gf_exp[RS_PRIM], y = beta;
    int order = 1;
    while (y != 1) {
        y = gf_mul(y, beta);
        order++;
    }
    CHECK(order == 255, "order of beta = %d", order);

    CHECK(rs_genpoly[0] == 1, "generator not monic");
    for (int j = 0; j < RS_NROOTS; j++)
        CHECK(gpoly_at(gf_exp[(RS_PRIM * (RS_FCR + j)) % 255]) == 0, "beta^%d not a root", RS_FCR + j);
    CHECK(gpoly_at(gf_exp[(RS_PRIM * (RS_FCR - 1)) % 255]) != 0, "beta^111 should not be a root");
    CHECK(gpoly_at(gf_exp[(RS_PRIM * (RS_FCR + 32)) % 255]) != 0, "beta^144 should not be a root");

    /* Roots come in reciprocal pairs beta^(112+j) * beta^(143-j) = beta^255,
     * so g(x) is self-reciprocal: the CCSDS generator is palindromic. */
    int pal = 1;
    for (int i = 0; i <= RS_NROOTS; i++)
        pal &= rs_genpoly[i] == rs_genpoly[RS_NROOTS - i];
    CHECK(pal, "generator is not palindromic");
}

static void test_dual_basis(void)
{
    int seen[256] = {0};
    for (int v = 0; v < 256; v++)
        seen[rs_to_dual[v]]++;
    int bij = 1;
    for (int v = 0; v < 256; v++)
        bij &= seen[v] == 1 && rs_to_conv[rs_to_dual[v]] == v;
    CHECK(bij, "dual-basis map is not a bijection");

    /* z_(7-m) = Tr(x alpha^(117 m)): the CCSDS symbol bits are the trace
     * coordinates against {1, alpha^117, ..., alpha^(117*7)}. */
    int ok = 1;
    for (int v = 0; v < 256; v++)
        for (int m = 0; m < 8; m++) {
            int bit = (rs_to_dual[v] >> (7 - m)) & 1;
            ok &= bit == gf_trace(gf_mul((uint8_t)v, gf_exp[(117 * m) % 255]));
        }
    CHECK(ok, "TAL matrix is not the trace dual of {alpha^(117 m)}");
}

static void test_encode_linear(void)
{
    uint8_t a[RS_N], b[RS_N], c[RS_N], s[RS_NROOTS];
    for (int trial = 0; trial < 200; trial++) {
        random_codeword(a);
        rs_syndromes(a, s);
        int z = 1;
        for (int j = 0; j < RS_NROOTS; j++)
            z &= s[j] == 0;
        CHECK(z, "codeword has nonzero syndrome");
        random_codeword(b);
        for (int i = 0; i < RS_K; i++)
            c[i] = a[i] ^ b[i];
        rs_encode_conv(c);
        int lin = 1;
        for (int i = RS_K; i < RS_N; i++)
            lin &= c[i] == (a[i] ^ b[i]);
        CHECK(lin, "encoder is not linear");
    }
}

static void test_errors_only(void)
{
    uint8_t cw[RS_N], rx[RS_N];
    int pos[RS_NROOTS];
    for (int e = 1; e <= RS_T; e++) {
        int good = 0, trials = 300;
        for (int trial = 0; trial < trials; trial++) {
            random_codeword(cw);
            memcpy(rx, cw, RS_N);
            pick(pos, e, 0, RS_N);
            for (int k = 0; k < e; k++)
                rx[pos[k]] ^= rnd_nonzero();
            int r = rs_decode_conv(rx, NULL, 0, 0);
            good += r == e && memcmp(rx, cw, RS_N) == 0;
        }
        CHECK(good == trials, "e=%d: %d/%d corrected", e, good, trials);
    }
}

static void test_errors_and_erasures(void)
{
    uint8_t cw[RS_N], rx[RS_N];
    int idx[RS_NROOTS];
    const int svals[] = {1, 2, 5, 10, 17, 24, 31, 32};
    for (unsigned q = 0; q < sizeof svals / sizeof svals[0]; q++) {
        int s = svals[q], e = (RS_NROOTS - s) / 2;
        int good = 0, trials = 200;
        for (int trial = 0; trial < trials; trial++) {
            random_codeword(cw);
            memcpy(rx, cw, RS_N);
            pick(idx, s + e, 0, RS_N);
            /* first s are erasures: about one in four left intact, since a
             * flagged symbol is often fine in practice */
            int wrong = 0;
            for (int k = 0; k < s + e; k++)
                if (k >= s || rnd() % 4) {
                    rx[idx[k]] ^= rnd_nonzero();
                    wrong++;
                }
            int r = rs_decode_conv(rx, idx, s, 0);
            good += r == wrong && memcmp(rx, cw, RS_N) == 0;
        }
        CHECK(good == trials, "s=%d e=%d: %d/%d corrected", s, e, good, trials);
    }
}

static void test_beyond_capacity(void)
{
    /* 17 errors: a bounded-distance decoder must refuse almost always. The
     * chance of landing inside another codeword's radius is about 1/16!. */
    uint8_t cw[RS_N], rx[RS_N];
    int pos[RS_T + 1];
    int refused = 0, trials = 2000;
    for (int trial = 0; trial < trials; trial++) {
        random_codeword(cw);
        memcpy(rx, cw, RS_N);
        pick(pos, RS_T + 1, 0, RS_N);
        for (int k = 0; k <= RS_T; k++)
            rx[pos[k]] ^= rnd_nonzero();
        uint8_t before[RS_N];
        memcpy(before, rx, RS_N);
        int r = rs_decode_conv(rx, NULL, 0, 0);
        refused += r < 0 && memcmp(rx, before, RS_N) == 0;
    }
    CHECK(refused == trials, "17 errors: refused %d/%d", refused, trials);

    /* s = 2, e = 16: 2e + s = 34 > 32 must also be refused */
    int idx[18];
    refused = 0;
    for (int trial = 0; trial < 500; trial++) {
        random_codeword(cw);
        memcpy(rx, cw, RS_N);
        pick(idx, 18, 0, RS_N);
        for (int k = 0; k < 18; k++)
            rx[idx[k]] ^= rnd_nonzero();
        refused += rs_decode_conv(rx, idx, 2, 0) < 0;
    }
    CHECK(refused == 500, "s=2 e=16: refused %d/500", refused);
}

static void test_ccsds_shortened(void)
{
    /* CCSDS virtual fill: a 255-byte code shortened by 100 bytes */
    const int pad = 100, k = RS_K - pad, n = RS_N - pad;
    uint8_t data[RS_K], blk[RS_N], tx[RS_N];
    int pos[RS_T], good = 0, trials = 300;
    for (int trial = 0; trial < trials; trial++) {
        for (int i = 0; i < k; i++)
            data[i] = (uint8_t)rnd();
        rs_ccsds_encode(data, tx, pad);
        memcpy(blk, tx, (size_t)n);
        pick(pos, RS_T, 0, n);
        for (int j = 0; j < RS_T; j++)
            blk[pos[j]] ^= rnd_nonzero();
        int r = rs_ccsds_decode(blk, NULL, 0, pad);
        good += r == RS_T && memcmp(blk, tx, (size_t)n) == 0 && memcmp(blk, data, (size_t)k) == 0;
    }
    CHECK(good == trials, "shortened (pad=100): %d/%d corrected", good, trials);

    /* The shortened codeword in conventional basis must equal the full
     * encoder's output with zero fill, i.e. a genuine (255,223) codeword. */
    uint8_t cw[RS_N] = {0}, s[RS_NROOTS];
    for (int j = 0; j < n; j++)
        cw[pad + j] = rs_to_conv[tx[j]];
    rs_syndromes(cw, s);
    int z = 1;
    for (int j = 0; j < RS_NROOTS; j++)
        z &= s[j] == 0;
    CHECK(z, "shortened codeword is not a (255,223) codeword");
}

static void test_interleaved_burst(void)
{
    enum { I = 5 };
    uint8_t data[I * RS_K], tx[I * RS_N], rx[I * RS_N];
    for (int i = 0; i < I * RS_K; i++)
        data[i] = (uint8_t)rnd();
    rs_ccsds_encode_interleaved(data, tx, I);

    /* any contiguous burst of 16*I bytes is correctable */
    int good = 0, trials = 200;
    for (int trial = 0; trial < trials; trial++) {
        int start = (int)(rnd() % (I * RS_N - RS_T * I + 1));
        memcpy(rx, tx, sizeof rx);
        for (int m = 0; m < RS_T * I; m++)
            rx[start + m] ^= rnd_nonzero();
        good += rs_ccsds_decode_interleaved(rx, I) == RS_T * I && memcmp(rx, tx, sizeof rx) == 0;
    }
    CHECK(good == trials, "I=5 burst of 80 bytes: %d/%d corrected", good, trials);

    /* 16*I + 1 bytes starting on a codeword-0 byte puts 17 in codeword 0 */
    memcpy(rx, tx, sizeof rx);
    for (int m = 0; m < RS_T * I + 1; m++)
        rx[10 * I + m] ^= rnd_nonzero();
    CHECK(rs_ccsds_decode_interleaved(rx, I) < 0, "I=5 burst of 81 bytes should fail");
}

int main(void)
{
    rs_init();
    test_field_and_generator();
    test_dual_basis();
    test_encode_linear();
    test_errors_only();
    test_errors_and_erasures();
    test_beyond_capacity();
    test_ccsds_shortened();
    test_interleaved_burst();
    printf("%d/%d checks passed\n", checks - fails, checks);
    return fails != 0;
}
