/*
 * Astron Built this AP
 * Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 */
#include "rice.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int checks, failures;

static void check(int ok, const char *what)
{
    checks++;
    if (!ok)
        failures++;
    printf("[%s] %s\n", ok ? "PASS" : "FAIL", what);
}

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static uint64_t rng_u64(void)
{
    /* xorshift64* */
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545F4914F6CDD1Dull;
}
static double rng_unit(void) { return (double)(rng_u64() >> 11) * 0x1.0p-53; }

/* Geometric(theta) on {0,1,...}: P(j) = (1-theta) theta^j. */
static uint32_t rng_geometric(double theta)
{
    return (uint32_t)floor(log(1.0 - rng_unit()) / log(theta));
}

static int roundtrip(const uint32_t *x, size_t N, const rice_params *p, rice_stats *st)
{
    uint8_t *buf;
    size_t len = rice_encode(x, N, p, &buf, st);
    if (!len)
        return 0;
    uint32_t *y = malloc(N * sizeof *y);
    int ok = rice_decode(buf, len, p, y, N) == 0 && memcmp(x, y, N * sizeof *x) == 0;
    free(y);
    free(buf);
    return ok;
}

/* Expected Rice codeword length for a geometric source, split parameter k:
 * k LSBs + unary quotient; P(q >= j) = theta^(j 2^k). */
static double rice_mean_len(double theta, int k)
{
    double t = pow(theta, ldexp(1.0, k));
    return k + 1.0 + t / (1.0 - t);
}

int main(void)
{
    char msg[160];

    /* 1. mapper is a bijection onto [0, 2^n - 1] for every prediction */
    {
        int ok = 1;
        for (int n = 1; n <= 12 && ok; n++) {
            uint32_t top = 1u << n;
            uint8_t *seen = malloc(top);
            for (uint32_t pred = 0; pred < top && ok; pred++) {
                memset(seen, 0, top);
                for (uint32_t x = 0; x < top; x++) {
                    uint32_t d = rice_map(x, pred, n);
                    if (d >= top || seen[d] || rice_unmap(d, pred, n) != x) {
                        ok = 0;
                        break;
                    }
                    seen[d] = 1;
                }
            }
            free(seen);
        }
        check(ok, "mapper: bijective with exact inverse for every (x, pred), n = 1..12");
        /* the interleave 0,-1,+1,-2,+2 near mid-scale */
        check(rice_map(100, 100, 16) == 0 && rice_map(99, 100, 16) == 1 &&
                  rice_map(101, 100, 16) == 2 && rice_map(98, 100, 16) == 3 &&
                  rice_map(0, 2, 16) == 3 && rice_map(65530, 65533, 16) == 5 &&
                  rice_map(65535, 65533, 16) == 4,
              "mapper: 2D / 2|D|-1 interleave and one-sided tail theta+|D| (spot values)");
    }

    /* 2. random round trips over every parameter combination */
    {
        static const int Js[] = {8, 16, 32, 64};
        static const int rsis[] = {1, 3, 64, 65, 200};
        int ok = 1, cases = 0;
        size_t N = 64 * 300;
        uint32_t *x = malloc(N * sizeof *x);
        for (int n = 1; n <= 16; n++)
            for (int ji = 0; ji < 4; ji++)
                for (int ri = 0; ri < 5; ri++)
                    for (int kind = 0; kind < 4; kind++) {
                        uint32_t top = (1u << n) - 1;
                        uint32_t v = (uint32_t)(rng_u64() % (top + 1u));
                        for (size_t i = 0; i < N; i++) {
                            switch (kind) {
                            case 0: x[i] = (uint32_t)(rng_u64() % (top + 1u)); break;
                            case 1: /* small steps, clipped at the rails */
                                v = (uint32_t)fmin(top, fmax(0, (double)v + (double)(rng_u64() % 5) - 2));
                                x[i] = v;
                                break;
                            case 2: /* plateaus with random lengths */
                                if (rng_u64() % 97 == 0)
                                    v = (uint32_t)(rng_u64() % (top + 1u));
                                x[i] = v;
                                break;
                            default: x[i] = (rng_u64() & 1) ? top : 0; break;
                            }
                        }
                        size_t NN = N - (size_t)Js[ji] * (size_t)(rng_u64() % 7);
                        rice_params p = {n, Js[ji], rsis[ri], -1};
                        ok &= roundtrip(x, NN, &p, NULL);
                        cases++;
                    }
        snprintf(msg, sizeof msg, "round trip: %d random cases, n = 1..16, J = 8..64, RSI 1..200", cases);
        check(ok, msg);

        int okk = 1;
        for (int k = 0; k <= 13; k++) {
            for (size_t i = 0; i < N; i++)
                x[i] = (uint32_t)(32768 + (int)(rng_u64() % 2001) - 1000);
            rice_params p = {16, 16, 64, k};
            okk &= roundtrip(x, N, &p, NULL);
        }
        check(okk, "round trip: forced split k = 0..13 at n = 16");
        free(x);
    }

    /* 3. exact sizes for white noise (all raw) and a constant (zero blocks) */
    {
        size_t N = 16 * 2557;
        uint32_t *x = malloc(N * sizeof *x);
        for (size_t i = 0; i < N; i++)
            x[i] = (uint32_t)(rng_u64() & 0xFFFF);
        rice_params p = {16, 16, 128, -1};
        rice_stats st;
        int ok = roundtrip(x, N, &p, &st);
        size_t expect = (size_t)st.blocks * (4 + 16 * 16);
        /* the mapped residual of white noise is not uniform (small values
         * are twice as likely), so a few blocks still beat raw with k = 13 */
        snprintf(msg, sizeof msg,
                 "white 16-bit noise: %ld/%ld blocks raw, size %zu <= blocks * (4 + nJ) = %zu bits",
                 st.uncompressed, st.blocks, st.bits, expect);
        check(ok && st.blocks == 2557 && st.uncompressed * 100 >= st.blocks * 99 &&
                  st.bits <= expect,
              msg);

        for (size_t i = 0; i < N; i++)
            x[i] = 1234;
        ok = roundtrip(x, N, &p, &st);
        /* full RSI of 128 blocks = two segments: [ID 0 ref ROS] + [ID 0 ROS];
         * the 125-block tail RSI splits 64 + 61, both still >= 5 -> ROS */
        size_t per_rsi = (4 + 1 + 16 + 5) + (4 + 1 + 5);
        size_t expect_c = 20 * per_rsi;
        snprintf(msg, sizeof msg, "constant: zero-block + ROS coding, %zu bits = 20 RSIs x %zu (libaec: 90 B)",
                 st.bits, per_rsi);
        check(ok && st.bits == expect_c && st.ros == 40 && st.zero_blocks == 2557 &&
                  (st.bits + 7) / 8 == 90,
              msg);
        free(x);
    }

    /* 4. zero-run length codes 1..4, ROS and >= 5 */
    {
        int ok = 1;
        for (int run = 1; run <= 12; run++) {
            /* block 0 busy, then `run` zero blocks, then busy blocks to RSI end */
            size_t J = 16, nb = 64, N = J * nb;
            uint32_t *x = malloc(N * sizeof *x);
            for (size_t i = 0; i < N; i++)
                x[i] = 500 + (uint32_t)(rng_u64() % 50);
            for (size_t i = J; i < J * (size_t)(1 + run) + 1 && i < N; i++)
                x[i] = x[J - 1];
            rice_params p = {12, 16, 64, -1};
            rice_stats st;
            ok &= roundtrip(x, N, &p, &st) && st.zero_blocks == run && st.ros == 0;
            free(x);
        }
        check(ok, "zero-block runs of 1..12 blocks mid-segment coded without ROS and recovered");
    }

    /* 5. second extension wins on sparse small residuals */
    {
        size_t N = 16 * 1000;
        uint32_t *x = malloc(N * sizeof *x);
        uint32_t v = 2000;
        for (size_t i = 0; i < N; i++) {
            if (rng_u64() % 20 == 0)
                v += (rng_u64() & 1) ? 1 : -1;
            x[i] = v;
        }
        rice_params p = {12, 16, 1000, -1};
        rice_stats st;
        int ok = roundtrip(x, N, &p, &st);
        snprintf(msg, sizeof msg, "second extension carries %ld of %ld coded blocks for p(step) = 5%% data",
                 st.second_ext, st.blocks);
        check(ok && st.second_ext > st.blocks / 2, msg);

        /* adaptive selection never loses to any fixed k on the same data */
        size_t best_forced = (size_t)-1;
        for (int k = 0; k <= 13; k++) {
            rice_params pf = {12, 16, 1000, k};
            rice_stats sf;
            roundtrip(x, N, &pf, &sf);
            if (sf.bits < best_forced)
                best_forced = sf.bits;
        }
        snprintf(msg, sizeof msg, "adaptive %zu bits <= best fixed-k %zu bits", st.bits, best_forced);
        check(st.bits <= best_forced, msg);
        free(x);
    }

    /* 6. fixed-k length vs. the closed-form geometric mean length */
    {
        int ok = 1;
        printf("    theta    k   measured  analytic   (bits/sample, payload only)\n");
        const double thetas[] = {0.5, 0.8, 0.95, 0.99};
        for (int t = 0; t < 4; t++) {
            double th = thetas[t];
            int J = 64;
            size_t nb = 4000, N = (size_t)J * nb;
            uint32_t *x = malloc(N * sizeof *x);
            x[0] = 32768;
            for (size_t i = 1; i < N; i++) {
                uint32_t d;
                do
                    d = rng_geometric(th);
                while (d > 8000);
                x[i] = rice_unmap(d, x[i - 1], 16);
            }
            for (int k = 0; k <= 7; k++) {
                rice_params p = {16, J, (int)nb, k};
                rice_stats st;
                roundtrip(x, N, &p, &st);
                /* strip IDs and the one reference sample */
                double payload = (double)st.bits - (double)nb * 4 - 16;
                double meas = payload / (double)(N - 1);
                double ana = rice_mean_len(th, k);
                /* std of the quotient: sqrt(t)/(1-t), t = theta^(2^k) */
                double tq = pow(th, ldexp(1.0, k));
                double sd = sqrt(tq) / (1 - tq) / sqrt((double)(N - 1));
                if (fabs(meas - ana) > 5 * sd + 1e-12)
                    ok = 0;
                if ((t == 1 && k <= 4) || (t == 3 && k >= 4 && k <= 7))
                    printf("    %.2f   %2d   %8.4f  %8.4f\n", th, k, meas, ana);
            }
            free(x);
        }
        check(ok, "fixed-k coded length = k + 1 + q/(1-q), q = theta^(2^k), within 5 sigma (32 cases)");
    }

    /* 7. Kiely's closed-form optimal k vs. brute-force argmin */
    {
        const double phi = (1.0 + sqrt(5.0)) / 2.0;
        int ok = 1, tested = 0;
        for (double lt = -6.0; lt <= -0.01; lt += 0.0137) {
            double th = 1.0 - pow(10.0, lt);
            int kf = (int)floor(log2(log(phi - 1.0) / log(th)));
            kf = kf + 1 > 0 ? kf + 1 : 0;
            int kb = 0;
            for (int k = 1; k < 30; k++)
                if (rice_mean_len(th, k) < rice_mean_len(th, kb))
                    kb = k;
            ok &= (kb == kf) ||
                  fabs(rice_mean_len(th, kb) - rice_mean_len(th, kf)) < 1e-9; /* exact ties */
            tested++;
        }
        snprintf(msg, sizeof msg, "k* = max(0, 1 + floor(log2(ln(phi-1)/ln theta))) is the argmin (%d thetas)", tested);
        check(ok, msg);
    }

    /* 8. malformed input is rejected, not over-read */
    {
        size_t N = 16 * 64;
        uint32_t *x = malloc(N * sizeof *x), *y = malloc(N * sizeof *y);
        for (size_t i = 0; i < N; i++)
            x[i] = (uint32_t)(rng_u64() & 0xFFF);
        rice_params p = {12, 16, 64, -1};
        uint8_t *buf;
        size_t len = rice_encode(x, N, &p, &buf, NULL);
        int ok = 1;
        for (size_t cut = 0; cut < len; cut += 7)
            ok &= rice_decode(buf, cut, &p, y, N) != 0;
        free(buf);
        free(x);
        free(y);
        check(ok, "truncated streams return an error at every cut point");
    }

    printf("\n%d/%d checks passed\n", checks - failures, checks);
    return failures ? 1 : 0;
}
