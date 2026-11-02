/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
#include "channel.h"

#include <math.h>

static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

void rng_seed(rng_t *r, uint64_t seed)
{
    /* splitmix64 expands the seed into the 256-bit state */
    for (int i = 0; i < 4; i++) {
        uint64_t z = (seed += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        r->s[i] = z ^ (z >> 31);
    }
}

uint64_t rng_next(rng_t *r)
{
    uint64_t *s = r->s;
    uint64_t res = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return res;
}

static double rng_unit(rng_t *r)   /* (0, 1) */
{
    return ((rng_next(r) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

double rng_gauss(rng_t *r)
{
    const double two_pi = 6.283185307179586;
    return sqrt(-2.0 * log(rng_unit(r))) * cos(two_pi * rng_unit(r));
}

void bpsk_awgn(const uint8_t *sym, size_t m, double ebn0_db, double rate,
               enum rx_mode mode, rng_t *r, double *out)
{
    /* Unit-energy symbols: Es = R * Eb, sigma^2 = N0 / 2 = 1 / (2 R Eb/N0). */
    double ebn0 = pow(10.0, ebn0_db / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * ebn0));
    for (size_t i = 0; i < m; i++) {
        double y = (1.0 - 2.0 * sym[i]) + sigma * rng_gauss(r);
        switch (mode) {
        case RX_HARD:
            out[i] = y >= 0.0 ? 1.0 : -1.0;
            break;
        case RX_Q3: {
            /* 8-level uniform quantizer, step 0.5 of the signal amplitude:
             * the classic 3-bit soft decision used in flight decoders. */
            double q = floor(y / 0.5);
            if (q < -4.0) q = -4.0;
            if (q > 3.0) q = 3.0;
            out[i] = 2.0 * q + 1.0;
            break;
        }
        case RX_SOFT:
            out[i] = y;
            break;
        }
    }
}

double q_func(double x) { return 0.5 * erfc(x / sqrt(2.0)); }
