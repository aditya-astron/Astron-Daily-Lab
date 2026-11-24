/*
 * Astron Built this AP
 * Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-24 IST
 */
#include "cic.h"

#include <math.h>

static uint64_t width_mask(int w)
{
    return (w >= 64) ? ~(uint64_t)0 : (((uint64_t)1 << w) - 1u);
}

int cic_full_width(int N, int R, int M, int bin)
{
    /* |y| <= (RM)^N 2^(bin-1)  ->  ceil(N log2(RM)) + bin bits signed */
    double g = (double)N * log2((double)R * (double)M);
    return (int)ceil(g - 1e-12) + bin;
}

int cic_init(cic_t *c, int N, int R, int M, int width, const int *prune)
{
    if (N < 1 || N > CIC_MAX_N || R < 1 || M < 1 || M > CIC_MAX_M ||
        width < 2 || width > 63)
        return -1;
    c->N = N;
    c->R = R;
    c->M = M;
    c->width_full = width;
    c->prune[0] = 0;
    for (int j = 1; j <= 2 * N + 1; j++) {
        int b = prune ? prune[j] : 0;
        if (b < c->prune[j - 1] || b >= width)
            return -1;
        c->prune[j] = b;
    }
    for (int j = 0; j <= 2 * N + 1; j++)
        c->mask[j] = width_mask(width - c->prune[j]);
    for (int i = 0; i < N; i++) {
        c->integ[i] = 0;
        for (int k = 0; k < CIC_MAX_M; k++)
            c->comb_dly[i][k] = 0;
    }
    c->dly_pos = 0;
    c->phase = 0;
    return 0;
}

int cic_push(cic_t *c, int32_t x, int64_t *out)
{
    const int N = c->N;

    /* Input: floor(x / 2^B1), then reduce modulo the stage-1 width. The
     * arithmetic shift on a signed value is floor division. */
    int64_t xs = (int64_t)x >> c->prune[1];
    uint64_t v = (uint64_t)xs & c->mask[1];
    c->integ[0] = (c->integ[0] + v) & c->mask[1];
    for (int i = 1; i < N; i++) {
        /* Logical shift of the unsigned residue is floor division modulo the
         * narrower width, so pruning commutes with wraparound. */
        uint64_t in = (c->integ[i - 1] >> (c->prune[i + 1] - c->prune[i])) & c->mask[i + 1];
        c->integ[i] = (c->integ[i] + in) & c->mask[i + 1];
    }

    if (++c->phase < c->R)
        return 0;
    c->phase = 0;

    uint64_t s = c->integ[N - 1];
    int pos = c->dly_pos;
    for (int i = 0; i < N; i++) {
        int j = N + 1 + i;                         /* comb stage index */
        uint64_t in = (s >> (c->prune[j] - c->prune[j - 1])) & c->mask[j];
        uint64_t d = c->comb_dly[i][pos];
        c->comb_dly[i][pos] = in;
        s = (in - d) & c->mask[j];
    }
    c->dly_pos = (pos + 1) % c->M;

    int jo = 2 * N + 1;
    s = (s >> (c->prune[jo] - c->prune[jo - 1])) & c->mask[jo];
    int w = c->width_full - c->prune[jo];
    /* sign-extend from w bits */
    uint64_t sign = (uint64_t)1 << (w - 1);
    *out = (int64_t)((s ^ sign) - sign);
    return 1;
}

size_t cic_run(cic_t *c, const int32_t *x, size_t n, int64_t *out)
{
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        m += (size_t)cic_push(c, x[i], &out[m]);
    return m;
}

size_t cic_run_i8(cic_t *c, const int8_t *x, size_t n, int64_t *out)
{
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        m += (size_t)cic_push(c, x[i], &out[m]);
    return m;
}

static double binom(int n, int k)
{
    if (k < 0 || k > n)
        return 0.0;
    if (k > n - k)
        k = n - k;
    double r = 1.0;
    for (int i = 1; i <= k; i++)
        r = r * (double)(n - k + i) / (double)i;
    return r;
}

double cic_hogenauer_F2(int N, int R, int M, int j)
{
    const int RM = R * M;
    double f2 = 0.0;
    if (j == 2 * N + 1)
        return 1.0;
    if (j > N) {
        int m = 2 * N + 1 - j;
        for (int k = 0; k <= m; k++) {
            double h = binom(m, k);
            f2 += h * h;
        }
        return f2;
    }
    /* H_j(z) = (1 - z^-RM)^N / (1 - z^-1)^(N-j+1) */
    int kmax = (RM - 1) * N + j - 1;
    for (int k = 0; k <= kmax; k++) {
        double h = 0.0;
        for (int l = 0; l <= k / RM; l++) {
            double t = binom(N, l) * binom(N - j + k - RM * l, k - RM * l);
            h += (l & 1) ? -t : t;
        }
        f2 += h * h;
    }
    return f2;
}

void cic_hogenauer_prune(int N, int R, int M, int width, int bout, int *prune)
{
    int bt = width - bout;
    if (bt < 0)
        bt = 0;
    double log2_sigma_t = (double)bt - 0.5 * log2(12.0);
    prune[0] = 0;
    for (int j = 1; j <= 2 * N; j++) {
        double fj = sqrt(cic_hogenauer_F2(N, R, M, j));
        int b = (int)floor(-log2(fj) + log2_sigma_t + 0.5 * log2(6.0 / N));
        if (b < 0)
            b = 0;
        if (b > bt)
            b = bt;
        if (b < prune[j - 1])
            b = prune[j - 1];
        prune[j] = b;
    }
    prune[2 * N + 1] = bt;
}

void sdm2_run(const double *x, size_t n, int8_t *y, double *e)
{
    double s1 = 0.0, s2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        int8_t yi = (s2 >= 0.0) ? 1 : -1;
        y[i] = yi;
        if (e)
            e[i] = (double)yi - s2;
        double s1n = s1 + x[i] - yi;
        s2 = s2 + s1 - 2.0 * yi;
        s1 = s1n;
    }
}
