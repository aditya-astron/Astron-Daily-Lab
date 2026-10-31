/*
 * Astron Built this AP
 * Project 002 | Telemetry framing | CCSDS TM link: packet mux, CADU sync, CRC-checked demux
 * Author: Aditya (aditya-astron)
 * Date: 2026-10-31 IST
 */
#include "ccsds.h"
#include <string.h>

/* CRC-16-CCITT, poly 0x1021, init 0xFFFF, no reflection, no final XOR.
 * This is the FECF defined in 132.0-B section 4.1.6. Check("123456789") = 0x29B1. */
uint16_t crc16_ccitt(const uint8_t *buf, size_t len)
{
    uint16_t crc = 0xFFFF;
    for (size_t i = 0; i < len; i++) {
        crc ^= (uint16_t)buf[i] << 8;
        for (int b = 0; b < 8; b++)
            crc = (crc & 0x8000) ? (uint16_t)((crc << 1) ^ 0x1021) : (uint16_t)(crc << 1);
    }
    return crc;
}

/* Fibonacci LFSR for h(x) = x^8 + x^7 + x^5 + x^3 + 1, i.e. the recurrence
 * a[n+8] = a[n+7] ^ a[n+5] ^ a[n+3] ^ a[n]. The register holds a[n] in bit 7
 * (next output) down to a[n+7] in bit 0. Period is 255 bits.               */
void pn_xor(uint8_t *buf, size_t len)
{
    uint8_t s = 0xFF;
    for (size_t i = 0; i < len; i++) {
        uint8_t pn = 0;
        for (int b = 0; b < 8; b++) {
            uint8_t out = (s >> 7) & 1;
            uint8_t fb  = ((s >> 7) ^ (s >> 4) ^ (s >> 2) ^ s) & 1;
            pn = (uint8_t)((pn << 1) | out);
            s  = (uint8_t)((s << 1) | fb);
        }
        buf[i] ^= pn;
    }
}

size_t sp_encode(const space_packet_t *p, uint8_t *out)
{
    /* version 000, type 0 (TM), no secondary header, seq flags 11 (unsegmented) */
    out[0] = (uint8_t)((p->apid >> 8) & 0x07);
    out[1] = (uint8_t)(p->apid & 0xFF);
    out[2] = (uint8_t)(0xC0 | ((p->seq_count >> 8) & 0x3F));
    out[3] = (uint8_t)(p->seq_count & 0xFF);
    uint16_t lf = (uint16_t)(p->data_len - 1);    /* length field is octets - 1 */
    out[4] = (uint8_t)(lf >> 8);
    out[5] = (uint8_t)(lf & 0xFF);
    memcpy(out + SP_HDR_LEN, p->data, p->data_len);
    return SP_HDR_LEN + p->data_len;
}

/* ------------------------------ transmitter ------------------------------ */

void tm_tx_init(tm_tx_t *tx, uint16_t scid, uint8_t vcid)
{
    memset(tx, 0, sizeof *tx);
    tx->scid = scid & 0x3FF;
    tx->vcid = vcid & 0x7;
}

static void tx_load_idle(tm_tx_t *tx, size_t room)
{
    /* Smallest legal packet is 7 octets. If less room is left, the idle
     * packet simply spills a few bytes into the next frame.                */
    static space_packet_t idle;
    idle.apid = APID_IDLE;
    idle.seq_count = tx->idle_seq++ & 0x3FFF;
    idle.data_len = (uint16_t)(room >= SP_HDR_LEN + 1 ? room - SP_HDR_LEN : 1);
    memset(idle.data, 0x55, idle.data_len);
    tx->pend_len = sp_encode(&idle, tx->pending);
    tx->pend_pos = 0;
}

void tm_tx_next_cadu(tm_tx_t *tx, packet_source_fn src, void *ctx, uint8_t cadu[CADU_LEN])
{
    static space_packet_t pkt;
    uint8_t *f = cadu + ASM_LEN;
    uint8_t *data = f + TM_HDR_LEN;
    size_t pos = 0;
    uint16_t fhp = FHP_NO_START;
    int idle_only = 1;

    while (pos < TM_DATA_LEN) {
        if (tx->pend_pos < tx->pend_len) {
            size_t n = tx->pend_len - tx->pend_pos;
            if (n > TM_DATA_LEN - pos) n = TM_DATA_LEN - pos;
            memcpy(data + pos, tx->pending + tx->pend_pos, n);
            if (((tx->pending[0] & 0x07) << 8 | tx->pending[1]) != APID_IDLE) idle_only = 0;
            pos += n;
            tx->pend_pos += n;
            continue;
        }
        if (fhp == FHP_NO_START) fhp = (uint16_t)pos;
        if (src && src(ctx, &pkt)) {
            tx->pend_len = sp_encode(&pkt, tx->pending);
            tx->pend_pos = 0;
        } else {
            tx_load_idle(tx, TM_DATA_LEN - pos);
        }
    }
    if (idle_only && fhp == 0) fhp = FHP_IDLE;

    f[0] = (uint8_t)(((tx->scid >> 4) & 0x3F));               /* TFVN=00 | SCID[9:4] */
    f[1] = (uint8_t)(((tx->scid & 0x0F) << 4) | (tx->vcid << 1)); /* SCID[3:0]|VCID|OCF=0 */
    f[2] = tx->mc_count++;
    f[3] = tx->vc_count++;
    /* data field status: no sec hdr, sync=0, order=0, segment length id = 11 */
    f[4] = (uint8_t)(0x18 | ((fhp >> 8) & 0x07));
    f[5] = (uint8_t)(fhp & 0xFF);

    uint16_t crc = crc16_ccitt(f, TM_FRAME_LEN - TM_FECF_LEN);
    f[TM_FRAME_LEN - 2] = (uint8_t)(crc >> 8);
    f[TM_FRAME_LEN - 1] = (uint8_t)(crc & 0xFF);

    /* Randomizer covers the whole frame including FECF, never the ASM. */
    pn_xor(f, TM_FRAME_LEN);
    cadu[0] = (uint8_t)(ASM_WORD >> 24);
    cadu[1] = (uint8_t)(ASM_WORD >> 16);
    cadu[2] = (uint8_t)(ASM_WORD >> 8);
    cadu[3] = (uint8_t)(ASM_WORD);
}

/* -------------------------------- receiver ------------------------------- */

void tm_rx_init(tm_rx_t *rx, uint16_t scid, uint8_t vcid, uint8_t asm_tol,
                packet_sink_fn sink, void *ctx)
{
    memset(rx, 0, sizeof *rx);
    rx->scid = scid & 0x3FF;
    rx->vcid = vcid & 0x7;
    rx->asm_tol = asm_tol;
    rx->sink = sink;
    rx->sink_ctx = ctx;
}

static void rx_drop_partial(tm_rx_t *rx)
{
    if (rx->buf_len) rx->st.packets_dropped++;
    rx->buf_len = 0;
}

static void rx_lose_sync(tm_rx_t *rx)
{
    rx_drop_partial(rx);
    rx->in_sync = 0;
}

/* Append bytes to the reassembly buffer, emitting each packet as it closes. */
static void rx_feed(tm_rx_t *rx, const uint8_t *d, size_t n)
{
    static space_packet_t pkt;
    for (size_t i = 0; i < n; i++) {
        rx->buf[rx->buf_len++] = d[i];
        if (rx->buf_len < SP_HDR_LEN) continue;
        size_t total = SP_HDR_LEN + ((size_t)rx->buf[4] << 8 | rx->buf[5]) + 1;
        if (rx->buf_len < total) continue;

        pkt.apid = (uint16_t)((rx->buf[0] & 0x07) << 8 | rx->buf[1]);
        pkt.seq_count = (uint16_t)((rx->buf[2] & 0x3F) << 8 | rx->buf[3]);
        pkt.data_len = (uint16_t)(total - SP_HDR_LEN);
        if (pkt.apid != APID_IDLE) {
            memcpy(pkt.data, rx->buf + SP_HDR_LEN, pkt.data_len);
            rx->st.packets_out++;
            if (rx->sink) rx->sink(rx->sink_ctx, &pkt);
        }
        rx->buf_len = 0;
    }
}

void tm_rx_frame(tm_rx_t *rx, const uint8_t raw[TM_FRAME_LEN])
{
    uint8_t f[TM_FRAME_LEN];
    memcpy(f, raw, TM_FRAME_LEN);
    pn_xor(f, TM_FRAME_LEN);

    uint16_t crc = crc16_ccitt(f, TM_FRAME_LEN - TM_FECF_LEN);
    if (crc != (uint16_t)(f[TM_FRAME_LEN - 2] << 8 | f[TM_FRAME_LEN - 1])) {
        rx->st.frames_crc_fail++;
        rx->crc_fail_run++;
        rx_lose_sync(rx);
        return;
    }

    uint16_t scid = (uint16_t)((f[0] & 0x3F) << 4 | f[1] >> 4);
    uint8_t vcid = (f[1] >> 1) & 0x07;
    if (scid != rx->scid || vcid != rx->vcid) return;

    /* A gap in the VC frame counter means frames went missing. Part of the
     * gap is frames we already rejected on CRC; the rest never locked at
     * all (ASM hit too hard). Any packet straddling the gap is lost.       */
    uint8_t vcc = f[3];
    if (rx->have_vc_count && vcc != rx->next_vc_count) {
        unsigned gap = (uint8_t)(vcc - rx->next_vc_count);
        rx->st.frames_missed += gap > rx->crc_fail_run ? gap - rx->crc_fail_run : 0;
        rx_lose_sync(rx);
    }
    rx->crc_fail_run = 0;
    rx->have_vc_count = 1;
    rx->next_vc_count = (uint8_t)(vcc + 1);
    rx->st.frames_ok++;

    uint16_t fhp = (uint16_t)((f[4] & 0x07) << 8 | f[5]);
    const uint8_t *data = f + TM_HDR_LEN;

    if (fhp == FHP_IDLE) {
        rx_drop_partial(rx);
        rx->in_sync = 1;            /* next frame begins on a packet boundary */
        return;
    }
    if (fhp != FHP_NO_START && fhp >= TM_DATA_LEN) {   /* malformed, ignore */
        rx_lose_sync(rx);
        return;
    }

    if (!rx->in_sync) {
        if (fhp == FHP_NO_START) return;   /* still mid-packet, keep waiting */
        rx->in_sync = 1;
        rx_feed(rx, data + fhp, TM_DATA_LEN - fhp);
        return;
    }

    if (fhp == FHP_NO_START) {
        rx_feed(rx, data, TM_DATA_LEN);
        return;
    }
    rx_feed(rx, data, fhp);
    if (rx->buf_len) rx_drop_partial(rx);   /* continuation disagrees with FHP */
    rx_feed(rx, data + fhp, TM_DATA_LEN - fhp);
}

static inline unsigned get_bit(const uint8_t *b, size_t k)
{
    return (b[k >> 3] >> (7 - (k & 7))) & 1u;
}

unsigned tm_rx_bitstream(tm_rx_t *rx, const uint8_t *bits, size_t nbits)
{
    const size_t frame_bits = (size_t)TM_FRAME_LEN * 8;
    unsigned locks = 0;
    uint32_t win = 0;
    size_t k = 0, filled = 0;

    while (k < nbits) {
        win = (win << 1) | get_bit(bits, k++);
        if (++filled < 32) continue;
        if ((unsigned)__builtin_popcount(win ^ ASM_WORD) > rx->asm_tol) continue;
        if (k + frame_bits > nbits) break;           /* truncated tail */

        uint8_t frame[TM_FRAME_LEN] = {0};
        for (size_t i = 0; i < frame_bits; i++)
            frame[i >> 3] |= (uint8_t)(get_bit(bits, k + i) << (7 - (i & 7)));
        tm_rx_frame(rx, frame);
        locks++;
        k += frame_bits;
        filled = 0;
        win = 0;
    }
    return locks;
}
