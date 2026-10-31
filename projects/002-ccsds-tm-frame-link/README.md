# 002 - CCSDS TM link: packet mux, CADU sync, CRC-checked demux

## Problem

Carry variable-length housekeeping packets from a spacecraft to the ground
over a noisy bit stream, following the CCSDS telemetry stack:

- Space Packets (133.0-B) for three APIDs: EPS, thermal, event log
- fixed 128-octet TM Transfer Frames (132.0-B) with a CRC-16 FECF
- pseudo-randomization and a 0x1ACFFC1D attached sync marker (131.0-B)

The ground side receives an unaligned bit stream with random errors and has
to find frame sync, reject corrupted frames, and reassemble packets that
straddle frame boundaries without ever delivering a damaged one.

## Approach

**Transmitter.** Packets are written back to back into the 120-octet frame
data field. A packet that does not fit spills into the next frame. The
11-bit First Header Pointer (FHP) gives the offset of the first packet header
that starts in the frame (0x7FF if none, 0x7FE for idle-only frames). When
the queue is empty the remaining space is filled with an idle packet
(APID 0x7FF) so a frame never ends on a dangling fragment.

**FECF.** CRC-16-CCITT: `g(x) = x^16 + x^12 + x^5 + 1`, preset 0xFFFF, over
the header and data field. Check value for "123456789" is 0x29B1.

**Randomizer.** `h(x) = x^8 + x^7 + x^5 + x^3 + 1`, seeded with all ones,
XORed over the whole frame (not the ASM). The sequence starts
`FF 48 0E C0 9A 0D 70 BC ...` and repeats every 255 bits. It guarantees bit
transitions for the receiver's clock recovery when telemetry is all zeros.

**Receiver.**

1. Slide a 32-bit window one bit at a time; declare sync where
   `popcount(window XOR ASM) <= 2`. False-lock probability per position is
   `(1 + 32 + 496) / 2^32 = 1.2e-7`.
2. Derandomize, recompute CRC. Failure drops the frame and any partial packet.
3. A gap in the 8-bit VC frame counter means frames went missing; packets
   crossing the gap are discarded.
4. After any loss, wait for a frame with a valid FHP and restart parsing
   there. Bytes before the FHP must exactly complete the packet in progress,
   otherwise that packet is dropped.

**Sanity check.** On a binary symmetric channel a frame is good only if all
1024 bits are, so `FER = 1 - (1 - p)^1024`. The demo prints the measured FER
next to this prediction.

## Files

| File | Purpose |
|------|---------|
| `ccsds.h`, `ccsds.c` | CRC, randomizer, packet encoder, frame mux, sync and demux |
| `tm_link.c` | end-to-end downlink simulation over a noisy channel |
| `test_ccsds.c` | unit and link tests |

## Run

```
make test          # unit + link tests
make run           # default: BER 2e-5, 400 frames
./tm_link 1e-4 2000 3
```

Tests cover the CRC check value, the first 16 randomizer octets and its
255-bit period, the packet header bit layout, a clean round trip of 60
packets (1 to 290 octets), a 29-bit stream offset with 2 bit errors in one
ASM, a single payload bit flip (exactly the packets overlapping that frame
are lost and the rest recovered via FHP), and idle-only frames.

## Sample output

```
CCSDS TM downlink  SCID 0x2A5 VC 1  frame 128 B  CADU 132 B
channel            BER 2.0e-05  lead-in 13 bits  10 bit errors in 422413 bits
frames             sent 400  locked 400  CRC ok 391  CRC fail 9  missed 0
frame error rate   measured 0.0225  predicted 1-(1-p)^1024 = 0.0240 (p = 2.37e-05 observed)
packets            sent 1194  delivered 1155  dropped partial 8  lost 39
  by APID          EPS 527  THERMAL 523  EVENT 105
efficiency         91.0% of channel bits carry packet octets
```

Nine bad frames cost 39 packets: each 120-octet data field holds about three
packets (mean length ~40 octets), plus the ones straddling each boundary. Adding Reed-Solomon (255,223) under the CRC would recover nearly all
of them at this BER; that is the natural next step.
