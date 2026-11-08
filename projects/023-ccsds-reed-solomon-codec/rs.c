/*
 * Astron Built this AP
 * Project 023 | Error correction | CCSDS Reed-Solomon (255,223) codec with errors-and-erasures decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-08 IST
 */
#include "rs.h"

#include <string.h>

#define GF_POLY 0x187   /* x^8 + x^7 + x^2 + x + 1 */

uint8_t gf_exp[2 * 255];
uint8_t gf_log[256];
uint8_t rs_genpoly[RS_NROOTS + 1];
uint8_t rs_to_dual[256];
uint8_t rs_to_conv[256];

/*
 * Berlekamp's transformation matrix from CCSDS 131.0-B (Annex F). Row 7-k is
 * the dual-basis image of the conventional basis element alpha^k. test_rs
 * checks that bit (7-m) of the image of x equals Tr(x * alpha^(117 m)), i.e.
 * that these rows really are the dual of the basis {alpha^(117 m)}.
 */
static const uint8_t TAL[8] = {0x8d, 0xef, 0xec, 0x86, 0xfa, 0x99, 0xaf, 0x7b};

static int mod255(int e)
{
    e %= 255;
    return e < 0 ? e + 255 : e;
}

/* alpha^e for any integer e */
static uint8_t apow(int e) { return gf_exp[mod255(e)]; }

uint8_t gf_mul(uint8_t a, uint8_t b)
{
    if (a == 0 || b == 0)
        return 0;
    return gf_exp[gf_log[a] + gf_log[b]];
}

uint8_t gf_div(uint8_t a, uint8_t b)
{
    if (a == 0)
        return 0;
    return gf_exp[gf_log[a] + 255 - gf_log[b]];
}

uint8_t gf_trace(uint8_t a)
{
    /* Tr(a) = a + a^2 + a^4 + ... + a^128, always 0 or 1 */
    uint8_t s = 0, y = a;
    for (int i = 0; i < 8; i++) {
        s ^= y;
        y = gf_mul(y, y);
    }
    return s;
}

void rs_init(void)
{
    unsigned x = 1;
    for (int i = 0; i < 255; i++) {
        gf_exp[i] = gf_exp[i + 255] = (uint8_t)x;
        gf_log[x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x100)
            x ^= GF_POLY;
    }

    /* g(x) = prod (x + beta^(FCR+j)), built up in ascending order then
     * reversed so rs_genpoly[0] is the x^32 coefficient. */
    uint8_t g[RS_NROOTS + 1] = {1};
    for (int j = 0; j < RS_NROOTS; j++) {
        uint8_t root = apow(RS_PRIM * (RS_FCR + j));
        for (int i = j + 1; i > 0; i--)
            g[i] = g[i - 1] ^ gf_mul(g[i], root);
        g[0] = gf_mul(g[0], root);
    }
    for (int i = 0; i <= RS_NROOTS; i++)
        rs_genpoly[i] = g[RS_NROOTS - i];

    for (int v = 0; v < 256; v++) {
        uint8_t d = 0;
        for (int k = 0; k < 8; k++)
            if (v & (1 << k))
                d ^= TAL[7 - k];
        rs_to_dual[v] = d;
        rs_to_conv[d] = (uint8_t)v;
    }
}

void rs_encode_conv(uint8_t cw[RS_N])
{
    /* Systematic encoding: parity = m(x) x^32 mod g(x), via the usual
     * division LFSR. g is monic so the feedback term needs no division. */
    uint8_t par[RS_NROOTS];
    memset(par, 0, sizeof par);
    for (int i = 0; i < RS_K; i++) {
        uint8_t fb = cw[i] ^ par[0];
        memmove(par, par + 1, RS_NROOTS - 1);
        par[RS_NROOTS - 1] = 0;
        if (fb)
            for (int j = 0; j < RS_NROOTS; j++)
                par[j] ^= gf_mul(fb, rs_genpoly[j + 1]);
    }
    memcpy(cw + RS_K, par, RS_NROOTS);
}

void rs_syndromes(const uint8_t cw[RS_N], uint8_t s[RS_NROOTS])
{
    for (int j = 0; j < RS_NROOTS; j++) {
        uint8_t root = apow(RS_PRIM * (RS_FCR + j));
        uint8_t acc = 0;
        for (int i = 0; i < RS_N; i++)          /* Horner, highest degree first */
            acc = gf_mul(acc, root) ^ cw[i];
        s[j] = acc;
    }
}

/* Locator of byte index i: X = beta^(254 - i), stored as a log to alpha. */
static int loc_log(int i) { return mod255(RS_PRIM * (RS_N - 1 - i)); }

static uint8_t poly_eval(const uint8_t *p, int deg, uint8_t x)
{
    uint8_t acc = 0;
    for (int i = deg; i >= 0; i--)
        acc = gf_mul(acc, x) ^ p[i];
    return acc;
}

int rs_decode_conv(uint8_t cw[RS_N], const int *eras, int neras, int pad)
{
    uint8_t s[RS_NROOTS];
    rs_syndromes(cw, s);

    int nonzero = 0;
    for (int j = 0; j < RS_NROOTS; j++)
        nonzero |= s[j];
    if (!nonzero)
        return 0;
    if (neras > RS_NROOTS)
        return -1;

    /* Erasure locator Gamma(x) = prod (1 + X_k x) seeds Berlekamp-Massey,
     * so the iteration only has to find the unknown errors. */
    uint8_t lambda[RS_NROOTS + 1] = {1};
    for (int k = 0; k < neras; k++) {
        if (eras[k] < pad || eras[k] >= RS_N)
            return -1;
        uint8_t xk = gf_exp[loc_log(eras[k])];
        for (int i = k + 1; i > 0; i--)
            lambda[i] ^= gf_mul(xk, lambda[i - 1]);
    }

    /* Berlekamp-Massey with erasures (Blahut's form). L counts errors plus
     * erasures; the update condition is shifted by the erasure count. */
    uint8_t b[RS_NROOTS + 1], t[RS_NROOTS + 1];
    memcpy(b, lambda, sizeof b);
    int L = neras;
    for (int r = neras + 1; r <= RS_NROOTS; r++) {
        uint8_t delta = 0;
        for (int i = 0; i < r && i <= RS_NROOTS; i++)
            delta ^= gf_mul(lambda[i], s[r - 1 - i]);

        /* b <- x b */
        memmove(b + 1, b, RS_NROOTS);
        b[0] = 0;
        if (delta == 0)
            continue;

        for (int i = 0; i <= RS_NROOTS; i++)
            t[i] = lambda[i] ^ gf_mul(delta, b[i]);
        if (2 * L <= r + neras - 1) {
            L = r + neras - L;
            /* length change: b <- old lambda / delta replaces x b */
            for (int i = 0; i <= RS_NROOTS; i++)
                b[i] = gf_div(lambda[i], delta);
        }
        memcpy(lambda, t, sizeof t);
    }

    int deg = 0;
    for (int i = 0; i <= RS_NROOTS; i++)
        if (lambda[i])
            deg = i;
    /* 2e + s <= 32 with e = deg - s */
    if (deg == 0 || 2 * deg - neras > RS_NROOTS)
        return -1;

    /* Chien search: byte i is in error iff Lambda(X_i^-1) = 0. Every nonzero
     * field element is some X_i, so deg distinct roots must all be found. */
    int pos[RS_NROOTS], npos = 0;
    for (int i = 0; i < RS_N; i++) {
        if (poly_eval(lambda, deg, gf_exp[mod255(-loc_log(i))]) == 0) {
            if (i < pad || npos == RS_NROOTS)
                return -1;
            pos[npos++] = i;
        }
    }
    if (npos != deg)
        return -1;

    /* Omega(x) = S(x) Lambda(x) mod x^32 */
    uint8_t omega[RS_NROOTS];
    for (int i = 0; i < RS_NROOTS; i++) {
        uint8_t acc = 0;
        for (int j = 0; j <= i && j <= deg; j++)
            acc ^= gf_mul(lambda[j], s[i - j]);
        omega[i] = acc;
    }

    /* Forney, for roots beta^(FCR+j):
     *   Y = X^(1-FCR) Omega(X^-1) / Lambda'(X^-1)
     * In characteristic 2 the formal derivative keeps only odd terms. */
    uint8_t fixed[RS_N];
    memcpy(fixed, cw, RS_N);
    int changed = 0;
    for (int k = 0; k < npos; k++) {
        int lx = loc_log(pos[k]);
        uint8_t xinv = gf_exp[mod255(-lx)];
        uint8_t num = poly_eval(omega, RS_NROOTS - 1, xinv);
        uint8_t den = 0;
        for (int i = 1; i <= deg; i += 2)
            den ^= gf_mul(lambda[i], gf_exp[mod255(-lx * (i - 1))]);
        if (den == 0)
            return -1;
        uint8_t y = gf_mul(gf_div(num, den), apow(lx * (1 - RS_FCR)));
        if (y) {
            fixed[pos[k]] ^= y;
            changed++;
        }
    }

    /* Cheap insurance against a decoder bug masquerading as a correction. */
    rs_syndromes(fixed, s);
    for (int j = 0; j < RS_NROOTS; j++)
        if (s[j])
            return -1;

    memcpy(cw, fixed, RS_N);
    return changed;
}

void rs_ccsds_encode(const uint8_t *data, uint8_t *out, int pad)
{
    uint8_t cw[RS_N] = {0};   /* virtual fill: leading pad bytes are zero */
    int k = RS_K - pad;
    for (int j = 0; j < k; j++)
        cw[pad + j] = rs_to_conv[data[j]];
    rs_encode_conv(cw);
    memcpy(out, data, (size_t)k);
    for (int j = 0; j < RS_NROOTS; j++)
        out[k + j] = rs_to_dual[cw[RS_K + j]];
}

int rs_ccsds_decode(uint8_t *blk, const int *eras, int neras, int pad)
{
    uint8_t cw[RS_N] = {0};
    int n = RS_N - pad;
    int e2[RS_NROOTS];
    if (neras > RS_NROOTS)
        return -1;
    for (int j = 0; j < n; j++)
        cw[pad + j] = rs_to_conv[blk[j]];
    for (int k = 0; k < neras; k++)
        e2[k] = eras[k] + pad;
    int r = rs_decode_conv(cw, e2, neras, pad);
    if (r > 0)
        for (int j = 0; j < n; j++)
            blk[j] = rs_to_dual[cw[pad + j]];
    return r;
}

void rs_ccsds_encode_interleaved(const uint8_t *data, uint8_t *blk, int depth)
{
    uint8_t d[RS_K], c[RS_N];
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < RS_K; j++)
            d[j] = data[j * depth + i];
        rs_ccsds_encode(d, c, 0);
        for (int j = 0; j < RS_N; j++)
            blk[j * depth + i] = c[j];
    }
}

int rs_ccsds_decode_interleaved(uint8_t *blk, int depth)
{
    uint8_t c[RS_N];
    int total = 0, failed = 0;
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < RS_N; j++)
            c[j] = blk[j * depth + i];
        int r = rs_ccsds_decode(c, NULL, 0, 0);
        if (r < 0) {
            failed = 1;
            continue;
        }
        for (int j = 0; j < RS_N; j++)
            blk[j * depth + i] = c[j];
        total += r;
    }
    return failed ? -1 : total;
}
