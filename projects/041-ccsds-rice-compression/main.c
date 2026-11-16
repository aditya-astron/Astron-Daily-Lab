/*
 * Astron Built this AP
 * Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * Compresses a day of representative small-satellite telemetry products,
 * compares the coded rate with the zeroth-order entropy of the mapped
 * residuals, and checks the coder against Gaussian and geometric theory.
 */
#include "rice.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t rng_state = 0xA5A5F00DCAFE1234ull;
static uint64_t rng_u64(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545F4914F6CDD1Dull;
}
static double rng_unit(void) { return ((double)(rng_u64() >> 11) + 0.5) * 0x1.0p-53; }
static double rng_normal(void)
{
    return sqrt(-2.0 * log(rng_unit())) * cos(2.0 * M_PI * rng_unit());
}

static uint32_t quantize(double v, int n)
{
    double top = ldexp(1.0, n) - 1.0;
    v = floor(v + 0.5);
    return (uint32_t)(v < 0 ? 0 : v > top ? top : v);
}

/* Zeroth-order entropy (bits/sample) of the mapped residuals, the rate a
 * memoryless ideal coder would need after the same predictor. */
static double residual_entropy(const uint32_t *x, size_t N, const rice_params *p)
{
    size_t span = (size_t)1 << p->n;
    size_t *h = calloc(span, sizeof *h);
    size_t per_rsi = (size_t)p->J * (size_t)p->rsi, cnt = 0;
    for (size_t i = 0; i < N; i++) {
        if (i % per_rsi == 0)
            continue; /* reference samples are sent raw */
        h[rice_map(x[i], x[i - 1], p->n)]++;
        cnt++;
    }
    double H = 0;
    for (size_t v = 0; v < span; v++)
        if (h[v]) {
            double q = (double)h[v] / (double)cnt;
            H -= q * log2(q);
        }
    free(h);
    return H;
}

typedef struct {
    double rate;   /* coded bits per sample */
    double H0;     /* residual entropy      */
    rice_stats st;
    int ok;        /* lossless round trip   */
} result;

static result run(const uint32_t *x, size_t N, const rice_params *p)
{
    result r;
    uint8_t *buf;
    size_t len = rice_encode(x, N, p, &buf, &r.st);
    uint32_t *y = malloc(N * sizeof *y);
    r.ok = len && rice_decode(buf, len, p, y, N) == 0 && !memcmp(x, y, N * sizeof *x);
    r.rate = 8.0 * (double)len / (double)N;
    r.H0 = residual_entropy(x, N, p);
    free(y);
    free(buf);
    return r;
}

static void report(const char *name, const uint32_t *x, size_t N, const rice_params *p,
                   double *raw_tot, double *cmp_tot)
{
    result r = run(x, N, p);
    long split = 0;
    for (int k = 0; k < 32; k++)
        split += r.st.split[k];
    double blocks = (double)(N / (size_t)p->J);
    double raw_kB = (double)N * p->n / 8e3, cmp_kB = r.rate * (double)N / 8e3;
    *raw_tot += raw_kB;
    *cmp_tot += cmp_kB;
    printf("%-26s %2d %8zu %9.1f %9.1f %5.2f %6.3f %6.3f  %4.1f %4.1f %4.1f %4.1f  %s\n", name, p->n, N,
           raw_kB, cmp_kB, raw_kB / cmp_kB, r.rate, r.H0, 100.0 * r.st.zero_blocks / blocks,
           100.0 * r.st.second_ext / blocks, 100.0 * split / blocks,
           100.0 * r.st.uncompressed / blocks, r.ok ? "ok" : "MISMATCH");
}

int main(void)
{
    rice_params p = {12, 16, 128, -1};
    double raw_tot = 0, cmp_tot = 0;

    printf("CCSDS 121.0 adaptive Rice coding of one day of 3U CubeSat telemetry\n");
    printf("J = 16 samples/block, 128 blocks per reference sample interval\n\n");
    printf("%-26s %2s %8s %9s %9s %5s %6s %6s  %4s %4s %4s %4s\n", "product", "n", "samples",
           "raw kB", "coded kB", "CR", "b/smp", "H0", "%zb", "%se", "%spl", "%raw");

    /* 1. panel thermistor: 1 Hz, 0.02 degC/LSB, orbital swing plus ADC noise */
    {
        size_t N = 86400;
        uint32_t *x = malloc(N * sizeof *x);
        for (size_t i = 0; i < N; i++) {
            double t = (double)i;
            double degc = 15.0 + 9.0 * sin(2 * M_PI * t / 5640.0) + 1.5 * sin(2 * M_PI * t / 86400.0);
            x[i] = quantize((degc + 40.0) / 0.02 + 0.35 * rng_normal(), 12);
        }
        p.n = 12;
        report("thermistor 1 Hz", x, N, &p, &raw_tot, &cmp_tot);
        free(x);
    }

    /* 2. magnetometer axis: 10 Hz, 2 nT/LSB, dipole field sweep, 6 nT noise */
    {
        size_t N = 864000;
        uint32_t *x = malloc(N * sizeof *x);
        for (size_t i = 0; i < N; i++) {
            double t = (double)i / 10.0;
            double nT = 38000.0 * sin(2 * M_PI * t / 5640.0) + 9000.0 * sin(4 * M_PI * t / 5640.0 + 0.7);
            x[i] = quantize(32768.0 + nT / 2.0 + 3.0 * rng_normal(), 16);
        }
        p.n = 16;
        report("magnetometer 10 Hz", x, N, &p, &raw_tot, &cmp_tot);
        free(x);
    }

    /* 3. imager: 12-bit pushbroom, 1024 rows x 2048 px; cloud/land texture
     *    as a smoothed random walk plus shot noise at 4 e-/DN */
    {
        size_t W = 2048, H = 1024, N = W * H;
        uint32_t *x = malloc(N * sizeof *x);
        double *scene = malloc(W * sizeof *scene);
        double s = 1200.0;
        for (size_t c = 0; c < W; c++) {
            s += 6.0 * rng_normal() + (c % 400 == 0 ? 600.0 * (rng_unit() - 0.4) : 0.0);
            s = fmin(fmax(s, 200.0), 3800.0);
            scene[c] = s;
        }
        for (size_t r = 0; r < H; r++)
            for (size_t c = 0; c < W; c++) {
                double dn = scene[c] * (0.8 + 0.2 * sin(2 * M_PI * (double)r / 700.0));
                x[r * W + c] = quantize(dn + sqrt(dn / 4.0) * rng_normal(), 12);
            }
        p.n = 12;
        report("imager 2048 x 1024", x, N, &p, &raw_tot, &cmp_tot);
        free(scene);
        free(x);
    }

    /* 4. already-encrypted payload: incompressible, the bound is ID overhead */
    {
        size_t N = 262144;
        uint32_t *x = malloc(N * sizeof *x);
        for (size_t i = 0; i < N; i++)
            x[i] = (uint32_t)(rng_u64() & 0xFFFF);
        p.n = 16;
        report("encrypted payload", x, N, &p, &raw_tot, &cmp_tot);
        printf("%-26s worst-case expansion bound 1 + L/(nJ) = %.4f\n", "",
               1.0 + 4.0 / (16.0 * 16.0));
        free(x);
    }
    printf("\nday total: %.1f kB raw -> %.1f kB coded (%.2f:1)\n", raw_tot, cmp_tot,
           raw_tot / cmp_tot);

    /* Gaussian sweep: unit-delay residual of white N(0, s^2) noise is
     * N(0, 2 s^2); its differential entropy is 0.5 log2(2 pi e 2 s^2). */
    printf("\nWhite Gaussian noise about mid-scale (n = 16, J = 16): coded rate vs. entropy\n");
    printf("%8s %8s %8s %8s %8s\n", "sigma", "b/smp", "H0", "h(D)", "rate-H0");
    {
        size_t N = 16 * 128 * 64;
        uint32_t *x = malloc(N * sizeof *x);
        for (int e = 0; e <= 10; e += 2) {
            double sigma = ldexp(1.0, e);
            for (size_t i = 0; i < N; i++)
                x[i] = quantize(32768.0 + sigma * rng_normal(), 16);
            rice_params g = {16, 16, 128, -1};
            result r = run(x, N, &g);
            double h = 0.5 * log2(2 * M_PI * M_E * 2 * sigma * sigma);
            printf("%8.0f %8.3f %8.3f %8.3f %8.3f%s\n", sigma, r.rate, r.H0, h, r.rate - r.H0,
                   r.ok ? "" : "  MISMATCH");
        }
        printf("\nBlock size trade at sigma = 8:\n%8s %8s\n", "J", "b/smp");
        for (size_t i = 0; i < N; i++)
            x[i] = quantize(32768.0 + 8.0 * rng_normal(), 16);
        for (int J = 8; J <= 64; J *= 2) {
            rice_params g = {16, J, 4096 / J, -1};
            result r = run(x, N, &g);
            printf("%8d %8.3f\n", J, r.rate);
        }
        free(x);
    }

    /* Geometric sources: Kiely's closed-form optimal k and its redundancy */
    printf("\nGeometric residuals P(j) = (1-t) t^j: optimal Rice parameter\n");
    printf("%8s %8s %4s %10s %10s %10s\n", "theta", "mean", "k*", "L(k*)", "H", "L-H");
    {
        const double phi = (1.0 + sqrt(5.0)) / 2.0;
        const double th[] = {0.3, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999};
        for (int i = 0; i < 7; i++) {
            double t = th[i];
            int k = 1 + (int)floor(log2(log(phi - 1.0) / log(t)));
            if (k < 0)
                k = 0;
            double q = pow(t, ldexp(1.0, k));
            double L = k + 1.0 + q / (1.0 - q);
            double H = (-(1 - t) * log2(1 - t) - t * log2(t)) / (1 - t);
            printf("%8.3f %8.2f %4d %10.4f %10.4f %10.4f\n", t, t / (1 - t), k, L, H, L - H);
        }
    }
    return 0;
}
