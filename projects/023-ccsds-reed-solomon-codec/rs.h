/*
 * Astron Built this AP
 * Project 023 | Error correction | CCSDS Reed-Solomon (255,223) codec with errors-and-erasures decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-08 IST
 *
 * CCSDS 131.0-B Reed-Solomon code, E = 16:
 *   field      GF(2^8), F(x) = x^8 + x^7 + x^2 + x + 1, alpha a root of F
 *   generator  g(x) = prod_{j=112}^{143} (x - beta^j),  beta = alpha^11
 *   symbols    transmitted in Berlekamp's dual basis
 *
 * Codewords are 255 bytes, byte 0 first on the wire and the coefficient of
 * x^254. The first 223 bytes are data, the last 32 are parity.
 */
#ifndef RS_H
#define RS_H

#include <stdint.h>

#define RS_N      255
#define RS_K      223
#define RS_NROOTS 32
#define RS_T      16
#define RS_FCR    112   /* first consecutive root exponent of beta */
#define RS_PRIM   11    /* beta = alpha^RS_PRIM */

/* Field tables, filled by rs_init(). exp is doubled so a sum of two logs
 * never needs a modulo. log[0] is unused. */
extern uint8_t gf_exp[2 * 255];
extern uint8_t gf_log[256];
extern uint8_t rs_genpoly[RS_NROOTS + 1];   /* conventional basis, g[0] = x^32 coeff */
extern uint8_t rs_to_dual[256];             /* conventional -> dual basis */
extern uint8_t rs_to_conv[256];             /* dual basis   -> conventional */

void rs_init(void);

uint8_t gf_mul(uint8_t a, uint8_t b);
uint8_t gf_div(uint8_t a, uint8_t b);       /* b != 0 */
uint8_t gf_trace(uint8_t a);

/* Conventional-basis codec on a full 255-byte block. */
void rs_encode_conv(uint8_t cw[RS_N]);      /* fills cw[223..254] from cw[0..222] */
void rs_syndromes(const uint8_t cw[RS_N], uint8_t s[RS_NROOTS]);

/*
 * Errors-and-erasures decoder. eras[] lists known-bad byte indices (0..254),
 * neras may be 0. Bytes with index < pad are a shortened code's virtual fill
 * and are treated as known zeros: an error located there is a failure.
 * Returns the number of symbols changed (errors + erasures that were wrong),
 * or -1 if the word is not within decoding distance; cw is untouched on -1.
 */
int rs_decode_conv(uint8_t cw[RS_N], const int *eras, int neras, int pad);

/*
 * CCSDS wire format: dual-basis symbols, shortened to 255 - pad bytes.
 * data has 223 - pad bytes, out has 255 - pad bytes.
 */
void rs_ccsds_encode(const uint8_t *data, uint8_t *out, int pad);
int  rs_ccsds_decode(uint8_t *blk, const int *eras, int neras, int pad);

/*
 * Interleaving depth I: byte m of the transmitted block belongs to codeword
 * m % I, so a channel burst of up to 16*I bytes lands at most 16 bytes in
 * any one codeword. data has I*223 bytes (codeword i's data at data[j*I+i]),
 * blk has I*255 bytes.
 */
void rs_ccsds_encode_interleaved(const uint8_t *data, uint8_t *blk, int depth);
int  rs_ccsds_decode_interleaved(uint8_t *blk, int depth);

#endif
