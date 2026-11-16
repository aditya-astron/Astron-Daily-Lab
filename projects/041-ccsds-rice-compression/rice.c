/*
 * Astron Built this AP
 * Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 */
#include "rice.h"

#include <stdlib.h>
#include <string.h>

#define SEGMENT_BLOCKS 64 /* zero-block runs never cross a 64-block segment */

/* ---------------------------------------------------------------- bits -- */

typedef struct {
    uint8_t *buf;
    size_t cap;
    size_t nbits;
} bitwriter;

static int bw_reserve(bitwriter *w, size_t extra_bits)
{
    size_t need = (w->nbits + extra_bits + 7) / 8;
    if (need <= w->cap)
        return 0;
    size_t cap = w->cap ? w->cap : 1024;
    while (cap < need)
        cap *= 2;
    uint8_t *b = realloc(w->buf, cap);
    if (!b)
        return -1;
    memset(b + w->cap, 0, cap - w->cap);
    w->buf = b;
    w->cap = cap;
    return 0;
}

/* MSB-first append of the low `len` bits of v (len <= 32). */
static void bw_put(bitwriter *w, uint32_t v, int len)
{
    for (int i = len - 1; i >= 0; i--) {
        if ((v >> i) & 1u)
            w->buf[w->nbits >> 3] |= (uint8_t)(0x80u >> (w->nbits & 7));
        w->nbits++;
    }
}

/* Fundamental sequence codeword: m zeros then a one. Buffer is pre-zeroed. */
static void bw_fs(bitwriter *w, uint64_t m)
{
    w->nbits += (size_t)m;
    bw_put(w, 1, 1);
}

typedef struct {
    const uint8_t *buf;
    size_t nbits;
    size_t pos;
} bitreader;

static int br_get(bitreader *r, int len, uint32_t *v)
{
    if (r->pos + (size_t)len > r->nbits)
        return -1;
    uint32_t x = 0;
    for (int i = 0; i < len; i++, r->pos++)
        x = (x << 1) | ((r->buf[r->pos >> 3] >> (7 - (r->pos & 7))) & 1u);
    *v = x;
    return 0;
}

static int br_fs(bitreader *r, uint64_t *m)
{
    uint64_t z = 0;
    for (;;) {
        if (r->pos >= r->nbits)
            return -1;
        int bit = (r->buf[r->pos >> 3] >> (7 - (r->pos & 7))) & 1;
        r->pos++;
        if (bit)
            break;
        z++;
    }
    *m = z;
    return 0;
}

/* -------------------------------------------------------- preprocessor -- */

int rice_id_len(int n)
{
    if (n > 16)
        return 5;
    if (n > 8)
        return 4;
    return 3;
}

int rice_kmax(int n)
{
    /* IDs 1 .. 2^L-2 carry k = ID-1; all-zeros is low entropy, all-ones raw */
    return (1 << rice_id_len(n)) - 3;
}

/*
 * theta is the distance from the prediction to the nearer end of the
 * sample range; residuals within +-theta interleave as 0,-1,+1,-2,...
 * (2D or 2|D|-1) and the one-sided remainder is appended linearly, so the
 * map is a bijection onto [0, 2^n - 1].
 */
uint32_t rice_map(uint32_t x, uint32_t pred, int n)
{
    int64_t xmax = ((int64_t)1 << n) - 1;
    int64_t theta = pred < (uint64_t)(xmax - pred) ? (int64_t)pred : xmax - pred;
    int64_t d = (int64_t)x - (int64_t)pred;
    if (d >= 0 && d <= theta)
        return (uint32_t)(2 * d);
    if (d < 0 && -d <= theta)
        return (uint32_t)(-2 * d - 1);
    return (uint32_t)(theta + (d < 0 ? -d : d));
}

uint32_t rice_unmap(uint32_t delta, uint32_t pred, int n)
{
    int64_t xmax = ((int64_t)1 << n) - 1;
    int64_t lo = pred, hi = xmax - (int64_t)pred;
    int64_t theta = lo < hi ? lo : hi;
    int64_t dl = delta;
    int64_t d;
    if (dl <= 2 * theta)
        d = (dl & 1) ? -(dl + 1) / 2 : dl / 2;
    else
        d = lo < hi ? dl - theta : -(dl - theta); /* only the roomy side */
    return (uint32_t)((int64_t)pred + d);
}

/* --------------------------------------------------------------- costs -- */

size_t rice_cost_split(const uint32_t *d, int J, int ref, int k)
{
    size_t bits = 0;
    for (int i = ref; i < J; i++)
        bits += (size_t)(d[i] >> k) + 1 + (size_t)k;
    return bits;
}

/* gamma = (a+b)(a+b+1)/2 + b enumerates pairs along anti-diagonals. */
static uint64_t se_gamma(uint32_t a, uint32_t b)
{
    uint64_t s = (uint64_t)a + b;
    return s * (s + 1) / 2 + b;
}

size_t rice_cost_se(const uint32_t *d, int J)
{
    size_t bits = 0;
    for (int i = 0; i < J; i += 2)
        bits += (size_t)se_gamma(d[i], d[i + 1]) + 1;
    return bits;
}

/* ------------------------------------------------------------- encoder -- */

static int valid_params(const rice_params *p)
{
    return p->n >= 1 && p->n <= 16 &&
           (p->J == 8 || p->J == 16 || p->J == 32 || p->J == 64) &&
           p->rsi >= 1 && p->rsi <= 4096 && p->force_k <= rice_kmax(p->n);
}

typedef struct {
    long count;      /* zero blocks pending                  */
    int has_ref;     /* run started on the RSI's first block */
    uint32_t ref;
} zero_run;

static int flush_run(bitwriter *w, const rice_params *p, zero_run *run,
                     int at_segment_end, rice_stats *st)
{
    int L = rice_id_len(p->n);
    if (bw_reserve(w, (size_t)(L + 1 + p->n + 70)))
        return -1;
    bw_put(w, 0, L);
    bw_put(w, 0, 1); /* low-entropy selector: 0 = zero-block */
    if (run->has_ref)
        bw_put(w, run->ref, p->n);
    uint64_t m;
    if (run->count >= 5 && at_segment_end) {
        m = 4; /* ROS: the rest of this segment is zero */
        st->ros++;
    } else if (run->count <= 4) {
        m = (uint64_t)run->count - 1;
    } else {
        m = (uint64_t)run->count;
    }
    bw_fs(w, m);
    st->zero_runs++;
    st->zero_blocks += run->count;
    run->count = 0;
    run->has_ref = 0;
    return 0;
}

static int emit_block(bitwriter *w, const rice_params *p, const uint32_t *d,
                      int ref, uint32_t ref_sample, rice_stats *st)
{
    const int J = p->J, n = p->n, L = rice_id_len(n);
    const int kmax = rice_kmax(n);

    int best_k = -1;
    size_t best = (size_t)(J - ref) * (size_t)n; /* no-compression */
    int use_se = 0;

    if (p->force_k >= 0) {
        best_k = p->force_k;
        best = rice_cost_split(d, J, ref, best_k);
    } else {
        for (int k = 0; k <= kmax && k < n; k++) {
            size_t c = rice_cost_split(d, J, ref, k);
            if (c < best) {
                best = c;
                best_k = k;
            }
        }
        size_t se = rice_cost_se(d, J) + 1; /* +1 selector bit */
        if (se < best) {
            best = se;
            use_se = 1;
        }
    }

    if (bw_reserve(w, (size_t)L + (size_t)n + best + 8))
        return -1;

    if (use_se) {
        bw_put(w, 0, L);
        bw_put(w, 1, 1);
        if (ref)
            bw_put(w, ref_sample, n);
        for (int i = 0; i < J; i += 2)
            bw_fs(w, se_gamma(d[i], d[i + 1]));
        st->second_ext++;
    } else if (best_k < 0) {
        bw_put(w, (1u << L) - 1, L);
        if (ref)
            bw_put(w, ref_sample, n);
        for (int i = ref; i < J; i++)
            bw_put(w, d[i], n);
        st->uncompressed++;
    } else {
        bw_put(w, (uint32_t)best_k + 1, L);
        if (ref)
            bw_put(w, ref_sample, n);
        /* all FS parts first, then all k-bit LSB fields */
        for (int i = ref; i < J; i++)
            bw_fs(w, d[i] >> best_k);
        if (best_k)
            for (int i = ref; i < J; i++)
                bw_put(w, d[i] & ((1u << best_k) - 1), best_k);
        st->split[best_k]++;
    }
    st->blocks++;
    return 0;
}

size_t rice_encode(const uint32_t *x, size_t N, const rice_params *p,
                   uint8_t **out, rice_stats *stats)
{
    rice_stats local;
    rice_stats *st = stats ? stats : &local;
    memset(st, 0, sizeof *st);
    *out = NULL;
    if (!valid_params(p) || N == 0 || N % (size_t)p->J)
        return 0;

    const int J = p->J;
    const uint32_t xmax = (1u << p->n) - 1;
    const size_t nblocks = N / (size_t)J;
    bitwriter w = {0};
    zero_run run = {0};
    uint32_t d[RICE_MAX_J];

    for (size_t b0 = 0; b0 < nblocks; b0 += (size_t)p->rsi) {
        size_t rsi_blocks = nblocks - b0 < (size_t)p->rsi ? nblocks - b0 : (size_t)p->rsi;
        for (size_t b = 0; b < rsi_blocks; b++) {
            const uint32_t *xb = x + (b0 + b) * (size_t)J;
            int ref = (b == 0);
            int all_zero = 1;
            for (int i = 0; i < J; i++) {
                if (xb[i] > xmax)
                    goto fail;
                if (ref && i == 0) {
                    d[0] = 0;
                    continue;
                }
                uint32_t pred = i ? xb[i - 1] : xb[-1];
                d[i] = rice_map(xb[i], pred, p->n);
                all_zero &= (d[i] == 0);
            }

            if (all_zero && p->force_k < 0) {
                if (run.count == 0) {
                    run.has_ref = ref;
                    run.ref = xb[0];
                }
                run.count++;
                int seg_end = ((b + 1) % SEGMENT_BLOCKS) == 0 || b + 1 == rsi_blocks;
                if (seg_end && flush_run(&w, p, &run, 1, st))
                    goto fail;
                continue;
            }
            if (run.count && flush_run(&w, p, &run, 0, st))
                goto fail;
            if (emit_block(&w, p, d, ref, xb[0], st))
                goto fail;
        }
    }

    st->bits = w.nbits;
    if (!w.buf && bw_reserve(&w, 8))
        goto fail;
    *out = w.buf;
    return (w.nbits + 7) / 8;

fail:
    free(w.buf);
    return 0;
}

/* ------------------------------------------------------------- decoder -- */

int rice_decode(const uint8_t *in, size_t len, const rice_params *p,
                uint32_t *x, size_t N)
{
    if (!valid_params(p) || N % (size_t)p->J)
        return -1;
    const int J = p->J, n = p->n, L = rice_id_len(n);
    const uint32_t raw_id = (1u << L) - 1;
    const size_t nblocks = N / (size_t)J;
    bitreader r = {in, len * 8, 0};
    uint32_t d[RICE_MAX_J];
    uint32_t v;
    uint64_t m;

    for (size_t b0 = 0; b0 < nblocks; b0 += (size_t)p->rsi) {
        size_t rsi_blocks = nblocks - b0 < (size_t)p->rsi ? nblocks - b0 : (size_t)p->rsi;
        size_t b = 0;
        while (b < rsi_blocks) {
            uint32_t *xb = x + (b0 + b) * (size_t)J;
            int ref = (b == 0);
            uint32_t ref_sample = 0;
            size_t count = 1;
            uint32_t id;
            if (br_get(&r, L, &id))
                return -1;

            if (id == 0) {
                uint32_t sel;
                if (br_get(&r, 1, &sel))
                    return -1;
                if (ref && br_get(&r, n, &ref_sample))
                    return -1;
                if (sel == 0) {
                    if (br_fs(&r, &m))
                        return -1;
                    if (m < 4) {
                        count = (size_t)m + 1;
                    } else if (m == 4) {
                        size_t seg_end = (b / SEGMENT_BLOCKS + 1) * SEGMENT_BLOCKS;
                        count = (seg_end < rsi_blocks ? seg_end : rsi_blocks) - b;
                    } else {
                        count = (size_t)m;
                    }
                    if (b + count > rsi_blocks)
                        return -1;
                    memset(d, 0, sizeof d);
                } else {
                    for (int i = 0; i < J; i += 2) {
                        if (br_fs(&r, &m))
                            return -1;
                        /* invert gamma: s = largest with s(s+1)/2 <= m */
                        uint64_t s = 0;
                        while ((s + 1) * (s + 2) / 2 <= m)
                            s++;
                        uint64_t bb = m - s * (s + 1) / 2;
                        if (s - bb > 0xFFFFFFFFu || bb > 0xFFFFFFFFu)
                            return -1;
                        d[i] = (uint32_t)(s - bb);
                        d[i + 1] = (uint32_t)bb;
                    }
                }
            } else if (id == raw_id) {
                if (ref && br_get(&r, n, &ref_sample))
                    return -1;
                for (int i = ref; i < J; i++)
                    if (br_get(&r, n, &d[i]))
                        return -1;
            } else {
                int k = (int)id - 1;
                if (ref && br_get(&r, n, &ref_sample))
                    return -1;
                for (int i = ref; i < J; i++) {
                    if (br_fs(&r, &m) || m > (0xFFFFFFFFu >> k))
                        return -1;
                    d[i] = (uint32_t)m << k;
                }
                if (k)
                    for (int i = ref; i < J; i++) {
                        if (br_get(&r, k, &v))
                            return -1;
                        d[i] |= v;
                    }
            }

            const uint32_t xmax = (1u << n) - 1;
            for (size_t c = 0; c < count; c++) {
                uint32_t *xc = xb + c * (size_t)J;
                for (int i = 0; i < J; i++) {
                    if (ref && c == 0 && i == 0) {
                        xc[0] = ref_sample;
                        continue;
                    }
                    if (d[i] > xmax)
                        return -1;
                    xc[i] = rice_unmap(d[i], i ? xc[i - 1] : xc[-1], n);
                }
            }
            b += count;
        }
    }
    return 0;
}
