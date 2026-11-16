/*
 * Astron Built this AP
 * Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 */
#ifndef RICE_H
#define RICE_H

#include <stddef.h>
#include <stdint.h>

/*
 * Lossless data compression per CCSDS 121.0-B (the "Rice" algorithm):
 * unit-delay predictor + prediction-error mapper, then per-block selection
 * among zero-block, second-extension, fundamental-sequence / split-sample
 * (k = 0..kmax) and no-compression options. Unsigned samples, 1 <= n <= 16.
 *
 * Stream layout matches the libaec convention with DATA_PREPROCESS and no
 * RSI padding: RSIs are concatenated bit-wise, the stream is zero-padded to
 * a byte at the very end.
 */

#define RICE_MAX_J 64

typedef struct {
    int n;        /* bits per sample                                   */
    int J;        /* block size: 8, 16, 32 or 64                        */
    int rsi;      /* blocks per reference sample interval (<= 4096)    */
    int force_k;  /* -1: adaptive selection; >=0: always split with k  */
} rice_params;

typedef struct {
    long blocks;         /* blocks coded                                  */
    long zero_blocks;    /* blocks absorbed in zero-block runs            */
    long zero_runs;      /* zero-block option codewords emitted           */
    long ros;            /* runs coded as "remainder of segment"          */
    long second_ext;     /* blocks coded with second extension            */
    long uncompressed;   /* blocks coded raw                              */
    long split[32];      /* blocks coded with split-sample k (k=0 is FS)  */
    size_t bits;         /* payload bits before final byte padding        */
} rice_stats;

/* ID length in bits for sample resolution n (non-restricted table). */
int rice_id_len(int n);
/* Largest split-sample k the ID field can express. */
int rice_kmax(int n);

/* Prediction-error mapper and its inverse (CCSDS 121.0 eq. 3-1). */
uint32_t rice_map(uint32_t x, uint32_t pred, int n);
uint32_t rice_unmap(uint32_t delta, uint32_t pred, int n);

/* Exact bit cost of one block under each option (ref = 1 if the block
 * carries the reference sample). d[] holds the J mapped residuals with
 * d[0] = 0 when ref is set. */
size_t rice_cost_split(const uint32_t *d, int J, int ref, int k);
size_t rice_cost_se(const uint32_t *d, int J);

/*
 * Encode N samples (N a multiple of J). Returns the byte count written to
 * *out (malloc'd, caller frees) or 0 on error. stats may be NULL.
 */
size_t rice_encode(const uint32_t *x, size_t N, const rice_params *p,
                   uint8_t **out, rice_stats *stats);

/* Decode exactly N samples. Returns 0 on success, -1 on malformed input. */
int rice_decode(const uint8_t *in, size_t len, const rice_params *p,
                uint32_t *x, size_t N);

#endif
