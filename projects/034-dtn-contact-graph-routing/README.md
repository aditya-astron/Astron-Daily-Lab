# 034 - Contact Graph Routing for a Lunar Far-Side DTN Relay

**Domain:** Algorithmic problem solving

## Problem

A rover on the lunar far side produces a 1.2 Gbit image product every
20 minutes. It never sees Earth. Data reaches mission control only by
store-and-forward through two relay orbiters and the three DSN complexes,
over links that exist during scheduled windows:

| Link | Window | Rate |
|------|--------|------|
| ROVER -> ORB1 (12 h frozen orbit) | 4 h around each apolune | 2 Mbit/s |
| ROVER -> ORB2 (2 h low polar orbit) | 12 min every 2 h | 8 Mbit/s |
| ORB1 -> DSS14 / DSS34 (Ka) | 1.5 h each, once a day | 20 Mbit/s |
| ORB2 -> DSS63 (X) | 1 h once a day | 4 Mbit/s |
| DSS -> MOC | always | 50 Mbit/s |

Which relay should carry each product so it lands as early as possible,
without overbooking any pass? ORB2 picks data up quickly but its one X-band
pass holds only 12 products a day, so the obvious local choice is wrong.

## Approach

Each contact `c = (tx, rx, t_start, t_end, rate, owlt)` is a vertex. A
bundle of `B` bits held at `c.tx` from time `t` can use it if

    t_tx   = max(t, t_start, free_at)
    t_done = t_tx + B / rate  <=  t_end       (no fragmentation)
    arrive = t_done + owlt                    (store-and-forward)

`free_at` is when the transmitter finishes bundles already booked on that
contact (FIFO, no backfill). The map `t -> arrive` is non-decreasing, so
waiting never helps a later start overtake an earlier one, and Dijkstra over
the contact graph (edge `c -> c'` when `c.rx == c'.tx`, loop-free) returns
the true earliest arrival. Bundles are routed in creation order; each route
is booked by advancing `free_at` along it, so later bundles see the residual
capacity. This is the core of CCSDS SABR / ION CGR without route caching.

Baselines: a greedy next-hop forwarder (hand to whichever neighbour can be
reached first, prefer the destination) and CGR with ORB2 out of service.

**Checks** (`test_cgr.py`):

- Hand calculation of the first bundle on an empty network:
  `1 h + 600 s` to ORB1, wait for the 5.5 h Goldstone pass, `60 s + 1.28 s`,
  then `24 s + 0.05 s` to MOC = 19885.33 s (5.5237 h).
- Dijkstra equals an exhaustive DFS over every loop-free contact sequence on
  400 random plans (6 nodes, 22 contacts, partly pre-booked), 250 of them
  reachable with routes of 1 to 4 hops.
- FIFO (non-overtaking) property of the hop map on 2000 random cases.
- Full replay of every booking: inside its window, no overlap on a contact,
  never sent before it is held, booked volume within contact capacity, and
  hop times reproduce the reported arrival.
- The ORB2 X-band pass is filled to exactly `3600 s x 4 Mbit/s / 1.2 Gbit = 12`
  products; greedy delivers fewer bundles and later; losing ORB2 raises mean
  latency.

**Limits.** Contact plan is idealised and exact (no range-dependent rates,
no plan errors or link outages), no bundle fragmentation, priorities or
custody, and one bundle in flight per contact.

## How to run

    python3 cgr.py        # scenario
    python3 test_cgr.py   # verification

Standard library only, Python 3.9+.

## Sample output

```
Contact plan: 37 contacts over 48 h, 72 rover bundles of 1200 Mbit

First bundle, empty network:
  ROVER -> ORB1  tx  1.0000 h  arrive  1.1667 h
   ORB1 -> DSS14 tx  5.5000 h  arrive  5.5170 h
  DSS14 -> MOC   tx  5.5170 h  arrive  5.5237 h

CGR (volume-aware Dijkstra): delivered 72/72 bundles
  latency  mean   9.46 h   median   9.42 h   max  17.86 h
   52 x ROVER > ORB1 > DSS14 > MOC
   12 x ROVER > ORB2 > DSS63 > MOC
    8 x ROVER > ORB1 > DSS34 > MOC
  DSN pass utilisation, day 1:
    ORB1 -> DSS14 at   5.5 h:  16.7 %
    ORB1 -> DSS34 at  13.0 h:   8.9 %
    ORB2 -> DSS63 at  19.0 h: 100.0 %

Greedy next-hop: delivered 50/72 bundles
  latency  mean  14.55 h   median  13.67 h   max  35.42 h
   24 x ROVER > ORB2 > DSS63 > MOC
   18 x ROVER > ORB1 > DSS14 > MOC
    8 x ROVER > ORB1 > DSS34 > MOC

CGR, ORB2 out of service: delivered 71/72 bundles
  latency  mean  12.15 h   median  13.31 h   max  21.86 h
   55 x ROVER > ORB1 > DSS14 > MOC
   16 x ROVER > ORB1 > DSS34 > MOC
```
