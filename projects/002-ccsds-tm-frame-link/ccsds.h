/*
 * Astron Built this AP
 * Project 002 | Telemetry framing | CCSDS TM link: packet mux, CADU sync, CRC-checked demux
 * Author: Aditya (aditya-astron)
 * Date: 2026-10-31 IST
 */
#ifndef CCSDS_H
#define CCSDS_H

#include <stddef.h>
#include <stdint.h>

/* Fixed-length TM Transfer Frame (CCSDS 132.0-B-3):
 *   6 B primary header | DATA_LEN B data field | 2 B FECF (CRC-16)        */
#define TM_HDR_LEN    6
#define TM_DATA_LEN   120
#define TM_FECF_LEN   2
#define TM_FRAME_LEN  (TM_HDR_LEN + TM_DATA_LEN + TM_FECF_LEN)   /* 128 */

#define ASM_WORD      0x1ACFFC1Du    /* attached sync marker, 131.0-B */
#define ASM_LEN       4
#define CADU_LEN      (ASM_LEN + TM_FRAME_LEN)

#define SP_HDR_LEN    6              /* space packet primary header, 133.0-B */
#define SP_MAX_LEN    (SP_HDR_LEN + 65536)
#define APID_IDLE     0x7FF

#define FHP_NO_START  0x7FF          /* no packet header starts in this frame */
#define FHP_IDLE      0x7FE          /* data field is idle fill only          */

uint16_t crc16_ccitt(const uint8_t *buf, size_t len);

/* Bernoulli pseudo-randomizer, h(x) = x^8 + x^7 + x^5 + x^3 + 1, seed all ones.
 * XOR is its own inverse, so the same call randomizes and derandomizes.      */
void pn_xor(uint8_t *buf, size_t len);

typedef struct {
    uint16_t apid;
    uint16_t seq_count;
    uint16_t data_len;               /* user data bytes, 1..65536 in the spec */
    uint8_t  data[SP_MAX_LEN - SP_HDR_LEN];
} space_packet_t;

/* Returns encoded length (SP_HDR_LEN + data_len). */
size_t sp_encode(const space_packet_t *p, uint8_t *out);

/* ---------- transmitter: packet stream -> fixed frames -> CADUs ---------- */

typedef struct {
    uint16_t scid;
    uint8_t  vcid;
    uint8_t  mc_count, vc_count;
    uint8_t  pending[SP_MAX_LEN];    /* bytes of a packet not yet framed */
    size_t   pend_len, pend_pos;
    uint16_t idle_seq;
} tm_tx_t;

void tm_tx_init(tm_tx_t *tx, uint16_t scid, uint8_t vcid);

/* Packet source: return 1 and fill *p, or 0 when nothing is queued. */
typedef int (*packet_source_fn)(void *ctx, space_packet_t *p);

/* Build one CADU (ASM + randomized frame). Packets are pulled from src as
 * room allows; a packet spills into following frames and the First Header
 * Pointer marks where the next new header begins. Unused space is filled
 * with an idle packet so the receiver never sees a dangling fragment.     */
void tm_tx_next_cadu(tm_tx_t *tx, packet_source_fn src, void *ctx, uint8_t cadu[CADU_LEN]);

/* ---------- receiver: bit stream -> frame sync -> frames -> packets ------ */

typedef struct {
    unsigned frames_ok, frames_crc_fail, frames_missed, packets_out, packets_dropped;
} tm_rx_stats_t;

typedef void (*packet_sink_fn)(void *ctx, const space_packet_t *p);

typedef struct {
    uint16_t scid;
    uint8_t  vcid;
    int      have_vc_count;
    uint8_t  next_vc_count;
    unsigned crc_fail_run;           /* CRC failures since last good frame */
    uint8_t  asm_tol;                /* allowed bit errors in the 32-bit ASM */
    uint8_t  buf[SP_MAX_LEN];        /* packet reassembly                    */
    size_t   buf_len;
    int      in_sync;                /* 0 until a header is located via FHP  */
    tm_rx_stats_t st;
    packet_sink_fn sink;
    void    *sink_ctx;
} tm_rx_t;

void tm_rx_init(tm_rx_t *rx, uint16_t scid, uint8_t vcid, uint8_t asm_tol,
                packet_sink_fn sink, void *ctx);

/* Process one frame exactly as it came off the channel (still randomized). */
void tm_rx_frame(tm_rx_t *rx, const uint8_t frame[TM_FRAME_LEN]);

/* Search a packed MSB-first bit stream for ASMs at any bit offset and feed
 * every frame found to tm_rx_frame. Returns number of ASMs locked.        */
unsigned tm_rx_bitstream(tm_rx_t *rx, const uint8_t *bits, size_t nbits);

#endif
