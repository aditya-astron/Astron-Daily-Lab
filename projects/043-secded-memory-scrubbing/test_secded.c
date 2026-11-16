/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * Exhaustive code-property checks and statistical checks of the closed-form
 * word failure model against a Monte Carlo of the real decoder.
 */
#include "scrub_sim.h"
#include "secded.h"

#include <math.h>
#include <stdio.h>

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static uint8_t syn(const int *pos, int k)
{
    uint8_t s = 0;
    for (int i = 0; i < k; i++) s ^= secded_col[pos[i]];
    return s;
}

static void test_matrix(void)
{
    for (int i = 0; i < 72; i++) {
        CHECK(__builtin_popcount(secded_col[i]) & 1, "column %d has even weight", i);
        for (int j = 0; j < i; j++)
            CHECK(secded_col[i] != secded_col[j], "columns %d and %d equal", i, j);
    }
    for (int i = 0; i < 8; i++)
        CHECK(__builtin_popcountll(secded_row[i]) == 26, "row %d weight %d != 26",
              i, __builtin_popcountll(secded_row[i]));
    printf("matrix: 72 distinct odd columns, every row 26 data ones + 1 check\n");
}

static void test_single_double(void)
{
    uint64_t s[4] = {1, 2, 3, 4};
    long singles = 0, doubles = 0;
    for (int trial = 0; trial < 64; trial++) {
        uint64_t d0 = rng_next(s);
        uint8_t c0 = secded_encode(d0);
        uint64_t d = d0; uint8_t c = c0; int pos;
        CHECK(secded_decode(&d, &c, &pos) == SECDED_CLEAN && d == d0, "clean word flagged");
        for (int i = 0; i < 72; i++) {
            d = d0; c = c0;
            if (i < 64) d ^= 1ULL << i; else c ^= (uint8_t)(1u << (i - 64));
            secded_status st = secded_decode(&d, &c, &pos);
            CHECK(st == SECDED_CORRECTED && pos == i && d == d0 && c == c0,
                  "single error at %d not corrected", i);
            singles++;
            for (int j = i + 1; j < 72; j++) {
                d = d0; c = c0;
                if (i < 64) d ^= 1ULL << i; else c ^= (uint8_t)(1u << (i - 64));
                if (j < 64) d ^= 1ULL << j; else c ^= (uint8_t)(1u << (j - 64));
                uint64_t dbad = d;
                st = secded_decode(&d, &c, &pos);
                CHECK(st == SECDED_DUE && d == dbad, "double error %d,%d not detected", i, j);
                doubles++;
            }
        }
    }
    printf("all %ld single errors corrected, all %ld double errors flagged DUE (64 random words)\n",
           singles, doubles);
}

/* Triple errors: odd syndrome, so the decoder either flags DUE (syndrome not
 * a column) or silently flips a fourth bit. Weight-4 errors whose syndrome is
 * zero are codewords and pass undetected. A4 is counted two independent ways:
 * direct enumeration, and from pairs sharing a syndrome (each weight-4
 * codeword splits into 3 pairs-of-pairs). */
static long A4_direct = 0, miscorrect3 = 0;
static void test_triples_quads(void)
{
    long due3 = 0;
    int p[4];
    for (p[0] = 0; p[0] < 72; p[0]++)
        for (p[1] = p[0] + 1; p[1] < 72; p[1]++)
            for (p[2] = p[1] + 1; p[2] < 72; p[2]++) {
                int k = secded_locate(syn(p, 3));
                if (k >= 0) miscorrect3++; else due3++;
                for (p[3] = p[2] + 1; p[3] < 72; p[3]++)
                    if (syn(p, 4) == 0) A4_direct++;
            }
    long npair[256] = {0}, A4_pairs = 0;
    for (int i = 0; i < 72; i++)
        for (int j = i + 1; j < 72; j++) npair[secded_col[i] ^ secded_col[j]]++;
    for (int v = 0; v < 256; v++) A4_pairs += npair[v] * (npair[v] - 1) / 2;
    CHECK(A4_pairs % 3 == 0 && A4_pairs / 3 == A4_direct,
          "A4 mismatch: direct %ld, pair count %ld/3", A4_direct, A4_pairs);
    CHECK(miscorrect3 + due3 == 59640, "triple count");
    /* Each weight-4 codeword holds 4 triples that miscorrect onto its 4th
     * bit, and each miscorrected triple belongs to exactly one codeword. */
    CHECK(miscorrect3 == 4 * A4_direct, "miscorrected triples %ld != 4*A4 %ld",
          miscorrect3, 4 * A4_direct);

    /* Decoder spot-check: a miscorrected triple must come back wrong but
     * unflagged, the failure the scrubber cannot see. */
    uint64_t d0 = 0x0123456789ABCDEFULL, d; uint8_t c0 = secded_encode(d0), c;
    int found = 0;
    for (p[0] = 0; p[0] < 64 && !found; p[0]++)
        for (p[1] = p[0] + 1; p[1] < 64 && !found; p[1]++)
            for (p[2] = p[1] + 1; p[2] < 64 && !found; p[2]++)
                if (secded_locate(syn(p, 3)) >= 0) {
                    int pos;
                    d = d0 ^ (1ULL << p[0]) ^ (1ULL << p[1]) ^ (1ULL << p[2]); c = c0;
                    CHECK(secded_decode(&d, &c, &pos) == SECDED_CORRECTED && d != d0,
                          "miscorrection not reproduced by decoder");
                    found = 1;
                }
    printf("triples: %ld miscorrect (%.2f %%), %ld DUE; weight-4 codewords A4 = %ld "
           "(direct = pair count / 3, triples = 4*A4)\n",
           miscorrect3, 100.0 * miscorrect3 / 59640.0, due3, A4_direct);
}

static void test_model_limits(void)
{
    /* f = 0: the correlated (I=1) recursion must equal the i.i.d. binomial. */
    double worst = 0.0;
    for (double lt = 1e-9; lt < 0.1; lt *= 3.7) {
        double a0, a1, b0, b1;
        double fa = word_fail_prob(1.0, 0.0, 1, lt, &a0, &a1);
        double fb = word_fail_prob(1.0, 0.0, 4, lt, &b0, &b1);
        double r = fabs(fa - fb) / fb;
        if (r > worst) worst = r;
        CHECK(fabs(a0 + a1 + fa - 1.0) < 1e-12, "probabilities do not sum to 1");
    }
    CHECK(worst < 1e-12, "I=1 recursion vs binomial at f=0: rel %.2e", worst);
    /* Small-lamT asymptote: P2 -> C(72,2) (lamT)^2 */
    double lt = 1e-7, pf = word_fail_prob(1.0, 0.0, 4, lt, NULL, NULL);
    CHECK(fabs(pf / (2556.0 * lt * lt) - 1.0) < 1e-4, "asymptote");
    /* I=1 with MCU: leading term is the MCU rate itself, 71 inner pairs plus
     * nothing else at first order (boundary spills only hit single bits). */
    lt = 1e-9;
    double pf1 = word_fail_prob(1.0, 0.2, 1, lt, NULL, NULL);
    CHECK(fabs(pf1 / (71.0 * 0.2 * lt) - 1.0) < 1e-5, "MCU first-order term %.6g", pf1 / (71.0 * 0.2 * lt));
    printf("model: I=1 recursion = binomial at f=0 (worst rel %.1e); C(72,2)(lamT)^2 and 71 f lamT asymptotes hold\n",
           worst);
}

static void test_monte_carlo(void)
{
    struct { int I; double lamT, f; } cases[] = {
        {1, 2e-3, 0.0}, {1, 2e-3, 0.3}, {4, 2e-3, 0.3}, {2, 8e-3, 0.5},
    };
    printf("Monte Carlo (8192 words x 400 intervals, real decoder) vs exact model:\n");
    printf("   I   lamT     f    P0 model  P0 MC     P1 model  P1 MC     Pfail model  Pfail MC    z\n");
    for (unsigned k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        scrub_tally t;
        scrub_mc(8192, cases[k].I, cases[k].lamT, cases[k].f, 400, 1000 + k, &t);
        double p0, p1, pf = word_fail_prob(1.0, cases[k].f, cases[k].I, cases[k].lamT, &p0, &p1);
        double n = (double)t.words, mf = (t.due + t.sdc) / n;
        double z = (mf - pf) / sqrt(pf * (1 - pf) / n);
        double z0 = (t.clean / n - p0) / sqrt(p0 * (1 - p0) / n);
        printf("  %2d  %.0e  %.1f   %.5f  %.5f   %.5f  %.5f   %.5f      %.5f  %+5.2f\n",
               cases[k].I, cases[k].lamT, cases[k].f, p0, t.clean / n, p1, t.corrected / n,
               pf, mf, z);
        /* Words are not fully independent (an I=1 MCU can straddle two
         * words) but the correlation is tiny; 4.5 sigma keeps false alarms
         * negligible across the cases. */
        CHECK(fabs(z) < 4.5 && fabs(z0) < 4.5, "MC disagrees with model, case %u", k);
    }
}

int main(void)
{
    secded_init();
    test_matrix();
    test_single_double();
    test_triples_quads();
    test_model_limits();
    test_monte_carlo();
    if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
    printf("all checks passed\n");
    return 0;
}
