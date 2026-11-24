/*
 * Astron Built this AP
 * Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-24 IST
 *
 * Firmware-level checks: bit-exactness against a wide reference, the tight
 * register width, Hogenauer's variance gains, and the pruning error model.
 */
#include "cic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures;

#define CHECK(cond, ...)                         \
    do {                                         \
        if (cond) {                              \
            printf("PASS  ");                    \
        } else {                                 \
            printf("FAIL  ");                    \
            failures++;                          \
        }                                        \
        printf(__VA_ARGS__);                     \
        printf("\n");                            \
    } while (0)

/* xorshift64*, deterministic */
static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static uint64_t rng(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 2685821657736338717ull;
}

/* Reference: direct form, N cascaded length-RM moving sums, then decimate,
 * in __int128 so nothing can overflow. Output k uses inputs up to index
 * (k+1)R - 1, matching the integrator-first structure. */
static size_t ref_cic(const int32_t *x, size_t n, int N, int R, int M, int64_t *out)
{
    const int RM = R * M;
    __int128 *a = malloc(n * sizeof *a), *b = malloc(n * sizeof *b);
    for (size_t i = 0; i < n; i++)
        a[i] = x[i];
    for (int s = 0; s < N; s++) {
        __int128 acc = 0;
        for (size_t i = 0; i < n; i++) {
            acc += a[i];
            if (i >= (size_t)RM)
                acc -= a[i - RM];
            b[i] = acc;
        }
        __int128 *t = a; a = b; b = t;
    }
    size_t m = 0;
    for (size_t i = (size_t)R - 1; i < n; i += (size_t)R)
        out[m++] = (int64_t)a[i];
    free(a);
    free(b);
    return m;
}

static void test_bit_exact(void)
{
    const size_t n = 200000;
    int32_t *x = malloc(n * sizeof *x);
    int64_t *y = malloc(n * sizeof *y), *r = malloc(n * sizeof *r);
    static const int cfg[][4] = {   /* N, R, M, Bin */
        {3, 64, 1, 2}, {4, 32, 2, 2}, {5, 16, 1, 12}, {2, 250, 1, 16},
    };
    for (size_t c = 0; c < sizeof cfg / sizeof cfg[0]; c++) {
        int N = cfg[c][0], R = cfg[c][1], M = cfg[c][2], bin = cfg[c][3];
        int32_t lo = -(1 << (bin - 1)), span = 1 << bin;
        for (size_t i = 0; i < n; i++)
            x[i] = lo + (int32_t)(rng() % (uint64_t)span);
        cic_t cic;
        int w = cic_full_width(N, R, M, bin);
        cic_init(&cic, N, R, M, w, NULL);
        size_t m = cic_run(&cic, x, n, y);
        size_t mr = ref_cic(x, n, N, R, M, r);
        size_t bad = (m == mr) ? 0 : 1;
        for (size_t i = 0; i < m && i < mr; i++)
            bad += (y[i] != r[i]);
        CHECK(bad == 0, "bit-exact vs 128-bit direct form: N=%d R=%d M=%d Bin=%d, width %d, %zu outputs",
              N, R, M, bin, w, m);
    }
    free(x); free(y); free(r);
}

/* Constant input drives the output to its extreme (RM)^N x. */
static int64_t const_output(int N, int R, int width, int32_t v)
{
    cic_t cic;
    cic_init(&cic, N, R, 1, width, NULL);
    int64_t y = 0;
    for (int i = 0; i < R * (N + 2); i++)
        cic_push(&cic, v, &y);
    return y;
}

static void test_register_width(void)
{
    const int N = 3, R = 96;
    const int64_t g = (int64_t)R * R * R;   /* 884736 */
    int w = cic_full_width(N, R, 1, 2);
    CHECK(w == 22, "Bmax for N=3, R=96, 2-bit input: %d bits (ceil(3 log2 96) + 2 = 22)", w);
    CHECK(const_output(N, R, w, -2) == -2 * g, "width %d, input -2: output %lld = -2 (RM)^N", w,
          (long long)const_output(N, R, w, -2));
    CHECK(const_output(N, R, w - 1, -2) != -2 * g,
          "width %d (one short), input -2: output %lld wraps", w - 1,
          (long long)const_output(N, R, w - 1, -2));
    /* A +/-1 bitstream never reaches -2, so one bit fewer is enough ... */
    CHECK(const_output(N, R, w - 1, 1) == g && const_output(N, R, w - 1, -1) == -g,
          "width %d, +/-1 stream: extremes +/-%lld represented", w - 1, (long long)g);
    /* ... but not two: +884736 > 2^19 - 1 */
    CHECK(const_output(N, R, w - 2, 1) != g, "width %d, +1 stream: wraps", w - 2);
}

static void test_hogenauer_gain(void)
{
    static const int cfg[][3] = {{3, 64, 1}, {4, 25, 2}, {5, 8, 1}};
    for (size_t c = 0; c < 3; c++) {
        int N = cfg[c][0], R = cfg[c][1], M = cfg[c][2], RM = R * M;
        double worst = 0.0;
        for (int j = 1; j <= 2 * N + 1; j++) {
            double f2;
            if (j <= N) {
                /* brute force: impulse through N-j+1 integrators then N combs */
                int len = (RM - 1) * N + j + 4 * RM;
                double *h = calloc((size_t)len, sizeof *h);
                h[0] = 1.0;
                for (int s = 0; s < N - j + 1; s++)
                    for (int i = 1; i < len; i++)
                        h[i] += h[i - 1];
                for (int s = 0; s < N; s++)
                    for (int i = len - 1; i >= RM; i--)
                        h[i] -= h[i - RM];
                f2 = 0.0;
                for (int i = 0; i < len; i++)
                    f2 += h[i] * h[i];
                free(h);
            } else {
                int m = 2 * N + 1 - j;
                double h[16] = {1.0};
                for (int s = 0; s < m; s++)
                    for (int i = m; i >= 1; i--)
                        h[i] -= h[i - 1];
                f2 = 0.0;
                for (int i = 0; i <= m; i++)
                    f2 += h[i] * h[i];
            }
            double rel = fabs(cic_hogenauer_F2(N, R, M, j) - f2) / f2;
            if (rel > worst)
                worst = rel;
        }
        CHECK(worst < 1e-12, "Hogenauer F_j^2 closed form vs brute-force impulse response, N=%d R=%d M=%d (max rel err %.1e)",
              N, R, M, worst);
    }
}

static void test_pruning_error(void)
{
    const int N = 4, R = 32, M = 1, bin = 2, bout = 16;
    const size_t n = (size_t)R * 400000;
    int w = cic_full_width(N, R, M, bin);
    int prune[2 * CIC_MAX_N + 2];
    cic_hogenauer_prune(N, R, M, w, bout, prune);

    /* predicted error: exact discrete truncation moments, propagated with
     * F_j^2 (variance) and the DC gain of each path (mean) */
    double var = 0.0, mean = 0.0;
    for (int j = 1; j <= 2 * N + 1; j++) {
        double a = ldexp(1.0, prune[j]), b = ldexp(1.0, prune[j - 1]);
        double vj = (a * a - b * b) / 12.0, mj = -(a - b) / 2.0;
        var += vj * cic_hogenauer_F2(N, R, M, j);
        if (j == 1)
            mean += mj * pow((double)(R * M), N);
        if (j == 2 * N + 1)
            mean += mj;
    }
    double sigma_t2 = ldexp(1.0, 2 * prune[2 * N + 1]) / 12.0;

    /* Uniform 2-bit input, as Hogenauer's model assumes. A +/-1 stream would
     * not do: a sum of k values of +/-1 has the parity of k, so the low
     * bits of the registers are not random. */
    int32_t *x = malloc(n * sizeof *x);
    for (size_t i = 0; i < n; i++)
        x[i] = (int32_t)(rng() & 3u) - 2;
    int64_t *yf = malloc(n / R * sizeof *yf), *yp = malloc(n / R * sizeof *yp);
    cic_t full, pr;
    cic_init(&full, N, R, M, w, NULL);
    cic_init(&pr, N, R, M, w, prune);
    size_t m = cic_run(&full, x, n, yf);
    cic_run(&pr, x, n, yp);

    double s = 0.0, ss = 0.0;
    size_t cnt = 0;
    for (size_t i = (size_t)N; i < m; i++) {
        double e = ldexp((double)yp[i], prune[2 * N + 1]) - (double)yf[i];
        s += e;
        ss += e * e;
        cnt++;
    }
    double em = s / cnt, ev = ss / cnt - em * em;
    printf("      N=%d R=%d width %d -> %d-bit output; discarded bits B_j:", N, R, w, bout);
    for (int j = 1; j <= 2 * N + 1; j++)
        printf(" %d", prune[j]);
    printf("\n");
    CHECK(var <= 2.0 * sigma_t2 * 1.0001,
          "Hogenauer budget: total error variance %.4g <= 2 sigma_T^2 = %.4g", var, 2.0 * sigma_t2);
    CHECK(fabs(ev / var - 1.0) < 0.03, "pruning error variance: measured %.5g vs predicted %.5g (%.2f %%)", ev,
          var, 100.0 * (ev / var - 1.0));
    double se_mean = sqrt(var / cnt);
    CHECK(fabs(em - mean) < 5.0 * se_mean, "pruning error mean: measured %.2f vs predicted %.2f (5 sigma = %.2f)",
          em, mean, 5.0 * se_mean);
    free(x); free(yf); free(yp);
}

static void test_modulator(void)
{
    const size_t n = 1u << 20;
    double *x = malloc(n * sizeof *x), *e = malloc(n * sizeof *e);
    int8_t *y = malloc(n);
    for (size_t i = 0; i < n; i++)
        x[i] = 0.5 * sin(2.0 * M_PI * 0.0013 * (double)i) + 0.1;
    sdm2_run(x, n, y, e);

    /* Y = z^-2 X + (1 - z^-1)^2 E must hold to round-off on every sample */
    double worst = 0.0, emax = 0.0;
    for (size_t i = 2; i < n; i++) {
        double r = (double)y[i] - x[i - 2] - (e[i] - 2.0 * e[i - 1] + e[i - 2]);
        if (fabs(r) > worst)
            worst = fabs(r);
        if (fabs(e[i]) > emax)
            emax = fabs(e[i]);
    }
    CHECK(worst < 1e-9, "modulator obeys Y = z^-2 X + (1-z^-1)^2 E exactly (max residual %.1e)", worst);
    CHECK(emax < 8.0, "quantizer error stays bounded for a 0.6 peak input (max |e| = %.2f, no overload)", emax);

    for (size_t i = 0; i < n; i++)
        x[i] = 0.3;
    sdm2_run(x, n, y, NULL);
    long sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += y[i];
    double dc = (double)sum / (double)n;
    CHECK(fabs(dc - 0.3) < 1e-4, "bitstream mean for DC 0.3 input: %.6f", dc);
    free(x); free(e); free(y);
}

int main(void)
{
    test_bit_exact();
    test_register_width();
    test_hogenauer_gain();
    test_pruning_error();
    test_modulator();
    printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
