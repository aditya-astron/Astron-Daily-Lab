/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 */
#include "scrub_sim.h"
#include "secded.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* xoshiro256** seeded by splitmix64; s points at 4 words of state. */
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

uint64_t rng_next(uint64_t *s)
{
    uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
}

double rng_uniform(uint64_t *s) { return (double)(rng_next(s) >> 11) * 0x1.0p-53; }

static void rng_seed(uint64_t *s, uint64_t seed)
{
    for (int i = 0; i < 4; i++) {
        uint64_t z = (seed += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        s[i] = z ^ (z >> 31);
    }
}

/* Physical cell p -> (word, bit). Cells are laid out in groups of 72*I;
 * inside a group, consecutive cells walk across the I words first, so two
 * neighbouring cells never share a word when I >= 2. */
static void cell_map(long p, int I, long *word, int *bit)
{
    long g = p / (72L * I);
    long r = p % (72L * I);
    *word = g * I + r % I;
    *bit = (int)(r / I);
}

static void flip(uint64_t *d, uint8_t *c, long w, int bit)
{
    if (bit < 64) d[w] ^= 1ULL << bit;
    else          c[w] ^= (uint8_t)(1u << (bit - 64));
}

void scrub_mc(int n_words, int interleave, double lamT, double f,
              int intervals, uint64_t seed, scrub_tally *t)
{
    uint64_t s[4];
    rng_seed(s, seed);
    long cells = 72L * n_words;
    uint64_t *d = malloc(sizeof *d * n_words), *d0 = malloc(sizeof *d0 * n_words);
    uint8_t  *c = malloc((size_t)n_words), *c0 = malloc((size_t)n_words);
    memset(t, 0, sizeof *t);

    for (int it = 0; it < intervals; it++) {
        for (int w = 0; w < n_words; w++) {
            d0[w] = d[w] = rng_next(s);
            c0[w] = c[w] = secded_encode(d0[w]);
        }
        /* Upsets over one interval: a Poisson process of total rate
         * cells*lamT on [0,1), drawn as exponential gaps (exact). */
        double M = (double)cells * lamT, tt = 0.0;
        for (;;) {
            tt += -log1p(-rng_uniform(s)) / M;
            if (tt >= 1.0) break;
            long p = (long)(rng_uniform(s) * (double)cells);
            long w; int b;
            cell_map(p, interleave, &w, &b);
            flip(d, c, w, b);
            if (rng_uniform(s) < f) {          /* 2-cell MCU, ring wrap */
                cell_map((p + 1) % cells, interleave, &w, &b);
                flip(d, c, w, b);
            }
        }
        for (int w = 0; w < n_words; w++) {
            int pos;
            secded_status st = secded_decode(&d[w], &c[w], &pos);
            int intact = d[w] == d0[w] && c[w] == c0[w];
            if (st == SECDED_DUE)              t->due++;
            else if (!intact)                  t->sdc++;
            else if (st == SECDED_CLEAN)       t->clean++;
            else                               t->corrected++;
        }
        t->words += n_words;
    }
    free(d); free(d0); free(c); free(c0);
}
