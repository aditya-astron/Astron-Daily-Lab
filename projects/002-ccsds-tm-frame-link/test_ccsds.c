/*
 * Astron Built this AP
 * Project 002 | Telemetry framing | CCSDS TM link: packet mux, CADU sync, CRC-checked demux
 * Author: Aditya (aditya-astron)
 * Date: 2026-10-31 IST
 */
#include "ccsds.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails;
#define CHECK(c) do { if (!(c)) { printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #c); fails++; } } while (0)

/* ---- packet source / sink used by the link tests ---- */
#define NPKT 60
typedef struct { int next, n; } src_t;
typedef struct { int got[NPKT]; int bad; } sink_t;

static uint16_t pkt_len(int i) { return (uint16_t)(1 + (i * 37) % 290); } /* 1..290 B */

static void make_pkt(int i, space_packet_t *p)
{
    p->apid = (uint16_t)(0x100 + i % 3);
    p->seq_count = (uint16_t)i;
    p->data_len = pkt_len(i);
    for (int k = 0; k < p->data_len; k++) p->data[k] = (uint8_t)(i * 7 + k);
}

static int src_fn(void *ctx, space_packet_t *p)
{
    src_t *s = ctx;
    if (s->next >= s->n) return 0;
    make_pkt(s->next++, p);
    return 1;
}

static void sink_fn(void *ctx, const space_packet_t *p)
{
    static space_packet_t ref;
    sink_t *k = ctx;
    if (p->seq_count >= NPKT) { k->bad++; return; }
    make_pkt(p->seq_count, &ref);
    if (p->apid != ref.apid || p->data_len != ref.data_len ||
        memcmp(p->data, ref.data, ref.data_len) != 0) { k->bad++; return; }
    k->got[p->seq_count]++;
}

static int build_cadus(uint8_t *out, int nframes)
{
    static tm_tx_t tx;
    src_t s = { 0, NPKT };
    tm_tx_init(&tx, 0x2A5, 3);
    for (int i = 0; i < nframes; i++) tm_tx_next_cadu(&tx, src_fn, &s, out + i * CADU_LEN);
    return s.next;
}

#define NFRAMES 100

static void test_crc(void)
{
    CHECK(crc16_ccitt((const uint8_t *)"123456789", 9) == 0x29B1);
    CHECK(crc16_ccitt(NULL, 0) == 0xFFFF);
}

static void test_pn(void)
{
    /* First 16 octets of the CCSDS 131.0-B pseudo-randomizer sequence. */
    static const uint8_t ref[16] = { 0xFF, 0x48, 0x0E, 0xC0, 0x9A, 0x0D, 0x70, 0xBC,
                                     0x8E, 0x2C, 0x93, 0xAD, 0xA7, 0xB7, 0x46, 0xCE };
    uint8_t z[255 * 2] = {0};
    pn_xor(z, sizeof z);
    CHECK(memcmp(z, ref, sizeof ref) == 0);
    /* 255-bit period: 255 octets = 8 periods, so the byte pattern repeats at 255. */
    CHECK(memcmp(z, z + 255, 255) == 0);
    pn_xor(z, sizeof z);                      /* involution */
    int allzero = 1;
    for (size_t i = 0; i < sizeof z; i++) allzero &= z[i] == 0;
    CHECK(allzero);
}

static void test_sp_header(void)
{
    static space_packet_t p;
    uint8_t out[SP_HDR_LEN + 4];
    p.apid = 0x123; p.seq_count = 0x2ABC; p.data_len = 4;
    memcpy(p.data, "\x01\x02\x03\x04", 4);
    CHECK(sp_encode(&p, out) == 10);
    static const uint8_t hdr[6] = { 0x01, 0x23, 0xEA, 0xBC, 0x00, 0x03 };
    CHECK(memcmp(out, hdr, 6) == 0);
}

static void test_clean_roundtrip(void)
{
    static uint8_t cadus[NFRAMES * CADU_LEN];
    int sent = build_cadus(cadus, NFRAMES);
    CHECK(sent == NPKT);                       /* 100 frames is enough room */

    sink_t k = {0};
    tm_rx_t rx;
    tm_rx_init(&rx, 0x2A5, 3, 0, sink_fn, &k);
    unsigned locks = tm_rx_bitstream(&rx, cadus, sizeof cadus * 8);
    CHECK(locks == NFRAMES);
    CHECK(rx.st.frames_ok == NFRAMES && rx.st.frames_crc_fail == 0);
    CHECK(k.bad == 0);
    int all = 1;
    for (int i = 0; i < NPKT; i++) all &= k.got[i] == 1;
    CHECK(all);
}

static void test_bit_slip_and_asm_errors(void)
{
    /* Shift the stream by 5 bits behind random junk, flip 2 bits in one ASM. */
    static uint8_t cadus[NFRAMES * CADU_LEN];
    static uint8_t bits[NFRAMES * CADU_LEN + 16];
    build_cadus(cadus, NFRAMES);
    const int shift = 5 + 8 * 3;               /* 29-bit offset */
    memset(bits, 0, sizeof bits);
    for (int i = 0; i < 3; i++) bits[i] = (uint8_t)(0x93 ^ i * 0x31);
    for (size_t i = 0; i < sizeof cadus * 8; i++) {
        unsigned b = (cadus[i >> 3] >> (7 - (i & 7))) & 1;
        size_t j = i + shift;
        bits[j >> 3] |= (uint8_t)(b << (7 - (j & 7)));
    }
    size_t asm7 = shift + (size_t)7 * CADU_LEN * 8;   /* ASM of frame 7 */
    bits[(asm7 + 3) >> 3] ^= (uint8_t)(1 << (7 - ((asm7 + 3) & 7)));
    bits[(asm7 + 20) >> 3] ^= (uint8_t)(1 << (7 - ((asm7 + 20) & 7)));

    sink_t k = {0};
    tm_rx_t rx;
    tm_rx_init(&rx, 0x2A5, 3, 2, sink_fn, &k);
    unsigned locks = tm_rx_bitstream(&rx, bits, (sizeof cadus * 8) + shift);
    CHECK(locks == NFRAMES);
    CHECK(k.bad == 0 && rx.st.packets_out == NPKT);

    /* With zero tolerance that frame is missed and the VC counter flags it. */
    sink_t k0 = {0};
    tm_rx_init(&rx, 0x2A5, 3, 0, sink_fn, &k0);
    tm_rx_bitstream(&rx, bits, (sizeof cadus * 8) + shift);
    CHECK(rx.st.frames_missed == 1);
    CHECK(k0.bad == 0);
}

static void test_corrupted_frame(void)
{
    /* One payload bit flip: CRC must reject the frame, packets touching it are
     * lost, and the receiver must resync on the next First Header Pointer. */
    static uint8_t cadus[NFRAMES * CADU_LEN];
    build_cadus(cadus, NFRAMES);
    const int victim = 4;
    cadus[victim * CADU_LEN + ASM_LEN + 50] ^= 0x08;

    sink_t k = {0};
    tm_rx_t rx;
    tm_rx_init(&rx, 0x2A5, 3, 0, sink_fn, &k);
    tm_rx_bitstream(&rx, cadus, sizeof cadus * 8);
    CHECK(rx.st.frames_crc_fail == 1);
    CHECK(k.bad == 0);

    /* Work out which packets overlap the victim frame's data field. */
    size_t lo = (size_t)victim * TM_DATA_LEN, hi = lo + TM_DATA_LEN, off = 0;
    int expect_lost = 0, delivered = 0, lost_ok = 1;
    for (int i = 0; i < NPKT; i++) {
        size_t a = off, b = off + SP_HDR_LEN + pkt_len(i);
        int overlaps = a < hi && b > lo;
        /* a packet that started earlier but ends inside the victim frame is
         * also lost; one that starts after a bad frame is recovered via FHP */
        expect_lost += overlaps;
        if (overlaps && k.got[i]) lost_ok = 0;
        if (!overlaps && k.got[i] != 1) lost_ok = 0;
        delivered += k.got[i];
        off = b;
    }
    CHECK(lost_ok);
    CHECK(delivered == NPKT - expect_lost);
    CHECK(expect_lost >= 1);
}

static void test_idle_frames(void)
{
    /* No source at all: every frame is pure idle fill with FHP = 0x7FE. */
    static tm_tx_t tx;
    uint8_t cadu[CADU_LEN];
    tm_tx_init(&tx, 0x2A5, 3);
    tm_tx_next_cadu(&tx, NULL, NULL, cadu);
    uint8_t f[TM_FRAME_LEN];
    memcpy(f, cadu + ASM_LEN, TM_FRAME_LEN);
    pn_xor(f, TM_FRAME_LEN);
    CHECK(((f[4] & 0x07) << 8 | f[5]) == FHP_IDLE);

    tm_rx_t rx;
    tm_rx_init(&rx, 0x2A5, 3, 0, NULL, NULL);
    tm_rx_frame(&rx, cadu + ASM_LEN);
    CHECK(rx.st.frames_ok == 1 && rx.st.packets_out == 0);
}

int main(void)
{
    test_crc();
    test_pn();
    test_sp_header();
    test_clean_roundtrip();
    test_bit_slip_and_asm_errors();
    test_corrupted_frame();
    test_idle_frames();
    if (fails) { printf("%d check(s) failed\n", fails); return 1; }
    printf("all tests passed\n");
    return 0;
}
