/*
 * Astron Built this AP
 * Project 002 | Telemetry framing | CCSDS TM link: packet mux, CADU sync, CRC-checked demux
 * Author: Aditya (aditya-astron)
 * Date: 2026-10-31 IST
 *
 * End-to-end downlink: three housekeeping APIDs are packetized, multiplexed
 * into 128-octet TM frames, randomized and prefixed with the ASM, then pushed
 * through a binary symmetric channel with an arbitrary bit offset. The ground
 * side finds frame sync, checks the FECF and reassembles packets.
 *
 * Usage: ./tm_link [BER] [frames] [seed]     defaults: 2e-5 400 7
 */
#include "ccsds.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SCID 0x2A5
#define VCID 1

enum { APID_EPS = 0x101, APID_THERMAL = 0x102, APID_EVENT = 0x103 };

static uint64_t rng = 7;
static uint32_t rnd32(void)                 /* xorshift64*, reproducible runs */
{
    rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27;
    return (uint32_t)((rng * 0x2545F4914F6CDD1DULL) >> 32);
}
static double rnd01(void) { return (rnd32() + 0.5) / 4294967296.0; }

typedef struct {
    unsigned tick;
    uint16_t seq[3];
    unsigned sent;
    unsigned sent_bytes;
} hk_src_t;

static void put_u16(uint8_t *p, uint16_t v) { p[0] = (uint8_t)(v >> 8); p[1] = (uint8_t)v; }

/* Each call emits the next packet of a fixed schedule: EPS and thermal every
 * tick, an event message every fifth tick. Lengths vary so packets straddle
 * frame boundaries at many different offsets.                              */
static int hk_source(void *ctx, space_packet_t *p)
{
    hk_src_t *s = ctx;
    unsigned slot = s->tick % 3;
    if (slot == 2 && (s->tick / 3) % 5 != 0) { s->tick++; slot = 0; }

    double t = (s->tick / 3) * 1.0;          /* seconds */
    if (slot == 0) {
        /* bus voltage [mV], array current [mA], battery SoC [0.1 %] */
        p->apid = APID_EPS;
        p->data_len = 6;
        put_u16(p->data + 0, (uint16_t)(28000 + 400 * sin(t / 30.0)));
        put_u16(p->data + 2, (uint16_t)(3200 + 900 * cos(t / 45.0)));
        put_u16(p->data + 4, (uint16_t)(870 - t / 10.0));
    } else if (slot == 1) {
        /* 24 thermistors, centi-degC */
        p->apid = APID_THERMAL;
        p->data_len = 48;
        for (int i = 0; i < 24; i++)
            put_u16(p->data + 2 * i, (uint16_t)(int16_t)(1500 + 120 * i + 300 * sin(t / 60.0 + i)));
    } else {
        p->apid = APID_EVENT;
        int n = snprintf((char *)p->data, 300,
                         "T+%05.0f EVT heater zone %u cycled; duty %u%%; "
                         "reaction wheel %u speed %d rpm; star tracker quaternion valid",
                         t, (unsigned)(s->tick % 8), 30 + s->tick % 50,
                         (unsigned)(s->tick % 4), 2000 + (int)(s->tick * 13 % 900));
        p->data_len = (uint16_t)n;
    }
    p->seq_count = s->seq[p->apid - APID_EPS]++ & 0x3FFF;
    s->tick++;
    s->sent++;
    s->sent_bytes += SP_HDR_LEN + p->data_len;
    return 1;
}

typedef struct { unsigned per_apid[3]; } gnd_t;

static void gnd_sink(void *ctx, const space_packet_t *p)
{
    gnd_t *g = ctx;
    if (p->apid >= APID_EPS && p->apid <= APID_EVENT) g->per_apid[p->apid - APID_EPS]++;
}

int main(int argc, char **argv)
{
    double ber = argc > 1 ? atof(argv[1]) : 2e-5;
    int nframes = argc > 2 ? atoi(argv[2]) : 400;
    rng = argc > 3 ? strtoull(argv[3], NULL, 10) : 7;
    if (nframes <= 0 || ber < 0 || ber >= 0.5) {
        fprintf(stderr, "usage: %s [BER] [frames] [seed]\n", argv[0]);
        return 2;
    }

    /* Channel: random lead-in of 1..63 bits, then the CADU stream. */
    size_t lead = 1 + rnd32() % 63;
    size_t nbits = lead + (size_t)nframes * CADU_LEN * 8;
    uint8_t *bits = calloc((nbits + 7) / 8, 1);
    if (!bits) return 1;
    for (size_t i = 0; i < lead; i++)
        if (rnd32() & 1) bits[i >> 3] |= (uint8_t)(0x80 >> (i & 7));

    static tm_tx_t tx;
    hk_src_t src = {0};
    tm_tx_init(&tx, SCID, VCID);
    size_t k = lead;
    uint8_t cadu[CADU_LEN];
    for (int f = 0; f < nframes; f++) {
        tm_tx_next_cadu(&tx, hk_source, &src, cadu);
        for (int i = 0; i < CADU_LEN * 8; i++, k++)
            if ((cadu[i >> 3] >> (7 - (i & 7))) & 1) bits[k >> 3] |= (uint8_t)(0x80 >> (k & 7));
    }
    /* Packets still in the transmitter's spill buffer never left the craft. */
    unsigned in_flight = tx.pend_pos < tx.pend_len ? 1 : 0;

    /* BSC: geometric gaps between errors instead of a coin flip per bit. */
    unsigned flips = 0;
    if (ber > 0) {
        double lg = log(1.0 - ber);
        for (double pos = floor(log(rnd01()) / lg); pos < (double)nbits;
             pos += 1.0 + floor(log(rnd01()) / lg)) {
            size_t b = (size_t)pos;
            bits[b >> 3] ^= (uint8_t)(0x80 >> (b & 7));
            flips++;
        }
    }

    static tm_rx_t rx;
    gnd_t g = {0};
    tm_rx_init(&rx, SCID, VCID, 2, gnd_sink, &g);
    unsigned locks = tm_rx_bitstream(&rx, bits, nbits);
    free(bits);

    unsigned expected = src.sent - in_flight;
    /* A frame survives only if all 1024 of its bits do. ASM hits within the
     * 2-bit tolerance are harmless, so they are left out of the estimate. */
    double p_meas = (double)flips / nbits;
    double fer_pred = 1.0 - pow(1.0 - p_meas, TM_FRAME_LEN * 8);
    printf("CCSDS TM downlink  SCID 0x%03X VC %d  frame %d B  CADU %d B\n",
           SCID, VCID, TM_FRAME_LEN, CADU_LEN);
    printf("channel            BER %.1e  lead-in %zu bits  %u bit errors in %zu bits\n",
           ber, lead, flips, nbits);
    printf("frames             sent %d  locked %u  CRC ok %u  CRC fail %u  missed %u\n",
           nframes, locks, rx.st.frames_ok, rx.st.frames_crc_fail, rx.st.frames_missed);
    printf("frame error rate   measured %.4f  predicted 1-(1-p)^%d = %.4f (p = %.2e observed)\n",
           (double)(nframes - rx.st.frames_ok) / nframes, TM_FRAME_LEN * 8, fer_pred, p_meas);
    printf("packets            sent %u  delivered %u  dropped partial %u  lost %u\n",
           expected, rx.st.packets_out, rx.st.packets_dropped, expected - rx.st.packets_out);
    printf("  by APID          EPS %u  THERMAL %u  EVENT %u\n",
           g.per_apid[0], g.per_apid[1], g.per_apid[2]);
    printf("efficiency         %.1f%% of channel bits carry packet octets\n",
           100.0 * src.sent_bytes / ((double)nframes * CADU_LEN));
    return 0;
}
