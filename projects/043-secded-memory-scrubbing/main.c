/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * Scrub interval design for a 256 MiB EDAC-protected payload memory in GEO.
 * Expected uncorrectable words over an exposure of D days with scrub period T:
 *     E(T) = N_words * Pfail(T) * D / T
 * Pfail ~ C(72,2) (r T)^2 for independent bits, so E grows linearly in T and
 * scrubbing faster always helps -- unless multi-cell upsets land two flips in
 * one word, which puts a floor 71 f lam N D under E that no scrub rate removes.
 */
#include "secded.h"

#include <math.h>
#include <stdio.h>

/* Assumed environment and part data (order-of-magnitude, for illustration). */
static const double N_WORDS = 33554432.0;  /* 256 MiB of data = 2^25 x 64 bit */
static const double LAM_QUIET = 2e-10;     /* upset events / cell / day, quiet GEO */
static const double LAM_SPE = 1e-7;        /* worst-week solar particle event      */
static const double F_MCU = 0.15;          /* fraction of events hitting 2 cells   */
static const double D_SPE = 7.0;           /* days                                 */
static const double D_LIFE = 15.0 * 365.25;
static const double TARGET = 0.01;         /* allowed expected uncorrectable words */

static double expected(double lam, int I, double T, double D)
{
    return N_WORDS * word_fail_prob(lam, F_MCU, I, T, NULL, NULL) * D / T;
}

/* Largest T (days) with E(T) <= TARGET; E is increasing in T. Returns -1 if
 * even an arbitrarily fast scrub misses (the MCU floor), D if never binding. */
static double max_period(double lam, int I, double D)
{
    double lo = 1.0 / 86400.0, hi = D;
    if (expected(lam, I, lo, D) > TARGET) return -1.0;
    if (expected(lam, I, hi, D) <= TARGET) return hi;
    for (int k = 0; k < 200; k++) {
        double mid = sqrt(lo * hi);
        if (expected(lam, I, mid, D) <= TARGET) lo = mid; else hi = mid;
    }
    return lo;
}

static void fmt_period(double days, char *buf, size_t n)
{
    if (days < 0)              snprintf(buf, n, "infeasible");
    else if (days < 1.0 / 24)  snprintf(buf, n, "%.1f min", days * 1440);
    else if (days < 2.0)       snprintf(buf, n, "%.2f h", days * 24);
    else                       snprintf(buf, n, "%.1f d", days);
}

int main(void)
{
    secded_init();
    printf("Hsiao (72,64) SECDED scrub design: %.0f words (256 MiB data), f_MCU = %.2f\n",
           N_WORDS, F_MCU);
    printf("environment: quiet %.0e, SPE %.0e events/cell/day; requirement E <= %.2g uncorrectable words\n\n",
           LAM_QUIET, LAM_SPE, TARGET);

    const double periods[] = {1.0 / 1440, 10.0 / 1440, 1.0 / 24, 6.0 / 24, 1.0, 7.0};
    const char *names[] = {"1 min", "10 min", "1 h", "6 h", "1 d", "7 d"};
    printf("expected uncorrectable words during a %.0f-day SPE:\n", D_SPE);
    printf("  scrub T    no interleave (I=1)   4-way interleave (I=4)   no MCU (f=0 reference)\n");
    for (int k = 0; k < 6; k++) {
        double T = periods[k];
        double e0 = N_WORDS * word_fail_prob(LAM_SPE, 0.0, 4, T, NULL, NULL) * D_SPE / T;
        printf("  %-8s   %12.4g          %12.4g             %12.4g\n", names[k],
               expected(LAM_SPE, 1, T, D_SPE), expected(LAM_SPE, 4, T, D_SPE), e0);
    }
    printf("  MCU floor at I=1 (T -> 0): 71 f lam N D = %.4g\n\n",
           71.0 * F_MCU * LAM_SPE * N_WORDS * D_SPE);

    printf("longest scrub period meeting E <= %.2g:\n", TARGET);
    printf("  case                I=1            I=4            I=4 closed form\n");
    struct { const char *name; double lam, D; } env[] = {
        {"7-day SPE", LAM_SPE, D_SPE}, {"15-yr quiet", LAM_QUIET, D_LIFE}};
    for (int k = 0; k < 2; k++) {
        double t1 = max_period(env[k].lam, 1, env[k].D);
        double t4 = max_period(env[k].lam, 4, env[k].D);
        double r = env[k].lam * (1.0 + F_MCU);
        double tc = TARGET / (N_WORDS * 2556.0 * r * r * env[k].D);
        char a[32], b[32], c[32];
        fmt_period(t1, a, sizeof a); fmt_period(t4, b, sizeof b); fmt_period(tc, c, sizeof c);
        printf("  %-14s  %-13s  %-13s  %-13s\n", env[k].name, a, b, c);
        if (k == 0 && t4 > 0)
            printf("  SPE scrub load at I=4: %.0f words/s (%.1f kB/s of EDAC read bandwidth)\n",
                   N_WORDS / (t4 * 86400.0), N_WORDS * 9.0 / (t4 * 86400.0) / 1e3);
    }
    return 0;
}
