/*
 * Astron Built this AP
 * Project 043 | Embedded firmware / Radiation effects | Hsiao (72,64) SECDED EDAC and memory scrub interval design
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * Hsiao odd-weight-column (72,64) single-error-correcting, double-error-
 * detecting code, the EDAC scheme used in front of most flight SDRAM/SRAM,
 * plus the exact word-failure model used to size the scrub period.
 */
#ifndef SECDED_H
#define SECDED_H

#include <stdint.h>

#define SECDED_DATA_BITS  64
#define SECDED_CHECK_BITS 8
#define SECDED_N          72

typedef enum {
    SECDED_CLEAN = 0,     /* syndrome zero                        */
    SECDED_CORRECTED = 1, /* odd syndrome matching one column     */
    SECDED_DUE = 2        /* detected uncorrectable error         */
} secded_status;

/* Parity-check matrix. Column j (0..63) is the 8-bit syndrome of a flip in
 * data bit j; check bit i sits at codeword position 64+i with column 1<<i. */
extern uint8_t  secded_col[SECDED_N];
extern uint64_t secded_row[SECDED_CHECK_BITS];

void          secded_init(void);
uint8_t       secded_encode(uint64_t data);
/* Decode in place. *pos receives the corrected codeword position (0..71) or -1. */
secded_status secded_decode(uint64_t *data, uint8_t *check, int *pos);
/* Syndrome -> codeword position, -1 if not a column. */
int           secded_locate(uint8_t syndrome);

/* ---- word failure model between scrubs (see README) ----
 * Each physical cell sees upset events at rate lam (per cell per unit time).
 * A fraction f are multi-cell upsets that also flip the next physical cell.
 * Interleave I: physically adjacent cells belong to I different words.
 * Returns P(number of flipped bits in a word >= 2) after an interval T,
 * i.e. the probability the scrubber finds the word uncorrectable (DUE or
 * silent miscorrection). Also returns P(exactly 0) and P(exactly 1) when
 * the pointers are non-null. Every term is a sum of positive products,
 * so it stays accurate down to 1e-300 instead of cancelling against 1. */
double word_fail_prob(double lam, double f, int interleave, double T,
                      double *p0, double *p1);

#endif
