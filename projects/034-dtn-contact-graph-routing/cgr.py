# Astron Built this AP
# Project 034 | Algorithmic problem solving | Contact graph routing for a lunar far-side DTN relay
# Author: Aditya (aditya-astron)
# Date: 2026-11-12 IST
"""
Contact graph routing (CGR) for a delay-tolerant lunar relay network.

A far-side rover never sees Earth. Its data reaches mission control by
store-and-forward through two relay orbiters and three DSN stations, over
links that exist only during scheduled contact windows. The contact plan is
known in advance, so routing is a shortest-path problem over time rather
than over space.

Model (one contact = one directed, scheduled link):

    c = (tx, rx, t_start, t_end, rate [bit/s], owlt [s])

A bundle of B bits that is fully held at node c.tx from time t can use c if

    t_tx  = max(t, c.t_start, c.free_at)     # wait for the window and queue
    t_end = t_tx + B / c.rate                 # serialisation
    t_end <= c.t_end                          # no fragmentation: fits whole
    arrival at c.rx = t_end + c.owlt          # store-and-forward

c.free_at is the time the transmitter finishes the bundles already booked
on that contact (FIFO, no backfill). The hop map t -> arrival is
non-decreasing in t (FIFO property), so Dijkstra on the contact graph,
with contacts as vertices and "c.rx == c'.tx" as edges, returns the true
earliest arrival. Booking a route advances free_at on each contact used;
routing bundles one after another in creation order therefore gives exact,
capacity-respecting delivery times.

All times are seconds from the plan epoch.
"""

from __future__ import annotations

import heapq
import math
from dataclasses import dataclass, field

INF = math.inf


@dataclass
class Contact:
    tx: str
    rx: str
    start: float
    end: float
    rate: float          # bit/s
    owlt: float          # one-way light time, s
    free_at: float = 0.0  # transmitter busy until (booked FIFO traffic)
    booked_bits: float = 0.0
    bookings: list = field(default_factory=list)  # (bundle id, t_tx, t_end)

    def hop(self, t_ready: float, bits: float) -> tuple[float, float] | None:
        """Return (t_tx, t_arrive) for a bundle ready at tx by t_ready, or None."""
        t_tx = max(t_ready, self.start, self.free_at)
        t_done = t_tx + bits / self.rate
        if t_done > self.end + 1e-9:
            return None
        return t_tx, t_done + self.owlt

    def residual_bits(self) -> float:
        return max(0.0, (self.end - max(self.start, self.free_at)) * self.rate)


@dataclass
class Route:
    contacts: list[int]      # indices into the plan
    arrival: float
    hops: list[tuple[float, float]]  # (t_tx, t_arrive) per contact

    def path(self, plan: list[Contact]) -> list[str]:
        if not self.contacts:
            return []
        nodes = [plan[self.contacts[0]].tx]
        nodes += [plan[i].rx for i in self.contacts]
        return nodes


def earliest_route(plan: list[Contact], src: str, dst: str, t0: float,
                   bits: float, exclude: set[str] | None = None) -> Route | None:
    """Dijkstra over contacts: earliest arrival of `bits` from src at t0 to dst.

    Loops are prevented by never re-entering a node already on the
    candidate path (the label at each contact carries its node set implicitly
    through the predecessor chain).
    """
    exclude = exclude or set()
    out: dict[str, list[int]] = {}
    for i, c in enumerate(plan):
        if c.rx not in exclude:
            out.setdefault(c.tx, []).append(i)

    best = [INF] * len(plan)
    pred = [-1] * len(plan)
    hopv: list[tuple[float, float] | None] = [None] * len(plan)
    done = [False] * len(plan)
    heap: list[tuple[float, int]] = []

    for i in out.get(src, []):
        h = plan[i].hop(t0, bits)
        if h and h[1] < best[i]:
            best[i], hopv[i] = h[1], h
            heapq.heappush(heap, (h[1], i))

    def on_path(i: int) -> set[str]:
        nodes = set()
        while i != -1:
            nodes.add(plan[i].tx)
            nodes.add(plan[i].rx)
            i = pred[i]
        return nodes

    while heap:
        t, i = heapq.heappop(heap)
        if done[i] or t > best[i]:
            continue
        done[i] = True
        c = plan[i]
        if c.rx == dst:
            # first final-hop contact popped is the earliest arrival at dst
            seq = []
            j = i
            while j != -1:
                seq.append(j)
                j = pred[j]
            seq.reverse()
            return Route(seq, t, [hopv[k] for k in seq])
        visited = on_path(i)
        for j in out.get(c.rx, []):
            if done[j] or plan[j].rx in visited:
                continue
            h = plan[j].hop(t, bits)
            if h and h[1] < best[j]:
                best[j], pred[j], hopv[j] = h[1], i, h
                heapq.heappush(heap, (h[1], j))
    return None


def book(plan: list[Contact], route: Route, bits: float, bid: int) -> None:
    for i, (t_tx, _) in zip(route.contacts, route.hops):
        c = plan[i]
        t_end = t_tx + bits / c.rate
        c.free_at = t_end
        c.booked_bits += bits
        c.bookings.append((bid, t_tx, t_end))


def route_and_book(plan, bundles, dst, exclude=None):
    """Route bundles (id, src, t_created, bits) in creation order; book each."""
    results = []
    for bid, src, t0, bits in sorted(bundles, key=lambda b: (b[2], b[0])):
        r = earliest_route(plan, src, dst, t0, bits, exclude)
        if r is not None:
            book(plan, r, bits, bid)
        results.append((bid, src, t0, bits, r))
    return results


def greedy_forward(plan, bundles, dst):
    """Baseline without a contact plan view: each node hands the bundle to
    whichever neighbour it can finish sending to first (never back to a
    visited node), preferring dst if it is reachable directly. Returns the
    same tuple format as route_and_book."""
    out: dict[str, list[int]] = {}
    for i, c in enumerate(plan):
        out.setdefault(c.tx, []).append(i)
    results = []
    for bid, src, t0, bits in sorted(bundles, key=lambda b: (b[2], b[0])):
        node, t, seen = src, t0, {src}
        seq, hops = [], []
        while node != dst:
            cands = []
            for i in out.get(node, []):
                if plan[i].rx in seen:
                    continue
                h = plan[i].hop(t, bits)
                if h:
                    cands.append((plan[i].rx != dst, h[1], i, h))
            if not cands:
                seq = None
                break
            _, _, i, h = min(cands)
            seq.append(i)
            hops.append(h)
            node, t = plan[i].rx, h[1]
            seen.add(node)
        r = Route(seq, t, hops) if seq is not None else None
        if r is not None:
            book(plan, r, bits, bid)
        results.append((bid, src, t0, bits, r))
    return results


# ---------------------------------------------------------------- scenario

H = 3600.0
OWLT_MOON = 1.28      # Earth-Moon one-way light time, s
OWLT_LOCAL = 0.02     # rover-orbiter, s (a few thousand km)
OWLT_GROUND = 0.05    # DSN complex to MOC over terrestrial network, s


def periodic(tx, rx, period, phase, dur, rate, owlt, horizon):
    out, k = [], 0
    while phase + k * period < horizon:
        s = phase + k * period
        out.append(Contact(tx, rx, s, min(s + dur, horizon), rate, owlt))
        k += 1
    return out


def lunar_plan(horizon: float = 48 * H) -> list[Contact]:
    """Idealised far-side relay plan.

    ORB1: 12 h frozen elliptical relay, apolune over the south far side. It
          sees the rover for 4 h around each apolune and holds Ka-band DSN
          passes of 1.5 h twice a day at 20 Mbit/s.
    ORB2: 2 h low polar orbiter, 12 min rover passes at 8 Mbit/s but only
          one 1 h X-band DSN pass per day at 4 Mbit/s.
    DSN:  each complex covers the Moon about 8 h a day; a pass is only
          booked inside that complex's visibility.
    Ground: DSN to MOC is always up at 50 Mbit/s.
    """
    p = []
    p += periodic("ROVER", "ORB1", 12 * H, 1 * H, 4 * H, 2e6, OWLT_LOCAL, horizon)
    p += periodic("ROVER", "ORB2", 2 * H, 0.25 * H, 0.2 * H, 8e6, OWLT_LOCAL, horizon)
    # ORB1 Ka passes, alternating Goldstone / Canberra / Madrid by time of day
    for day in range(int(math.ceil(horizon / (24 * H)))):
        d = day * 24 * H
        p.append(Contact("ORB1", "DSS14", d + 5.5 * H, d + 7.0 * H, 20e6, OWLT_MOON))
        p.append(Contact("ORB1", "DSS34", d + 13.0 * H, d + 14.5 * H, 20e6, OWLT_MOON))
        p.append(Contact("ORB2", "DSS63", d + 19.0 * H, d + 20.0 * H, 4e6, OWLT_MOON))
    for dss in ("DSS14", "DSS34", "DSS63"):
        p.append(Contact(dss, "MOC", 0.0, horizon, 50e6, OWLT_GROUND))
    p = [c for c in p if c.start < horizon]
    p.sort(key=lambda c: (c.start, c.tx, c.rx))
    return p


def rover_traffic(horizon: float, every: float = 20 * 60, bits: float = 1.2e9):
    """One 1.2 Gbit image product every 20 min for the first 24 h."""
    n = int(24 * H // every)
    return [(k, "ROVER", k * every, bits) for k in range(n) if k * every < horizon]


def summarise(name, results, plan):
    lat = [r.arrival - t0 for _, _, t0, _, r in results if r is not None]
    n, ok = len(results), len(lat)
    lines = [f"{name}: delivered {ok}/{n} bundles"]
    if lat:
        lat_sorted = sorted(lat)
        med = lat_sorted[len(lat) // 2]
        lines.append(f"  latency  mean {sum(lat)/ok/H:6.2f} h   median {med/H:6.2f} h"
                     f"   max {max(lat)/H:6.2f} h")
    via = {}
    for *_, r in results:
        if r is not None:
            key = " > ".join(r.path(plan))
            via[key] = via.get(key, 0) + 1
    for k, v in sorted(via.items(), key=lambda kv: -kv[1]):
        lines.append(f"  {v:3d} x {k}")
    return "\n".join(lines)


def main():
    horizon = 48 * H
    traffic = rover_traffic(horizon)

    plan = lunar_plan(horizon)
    print(f"Contact plan: {len(plan)} contacts over {horizon/H:.0f} h, "
          f"{len(traffic)} rover bundles of {traffic[0][3]/1e6:.0f} Mbit\n")

    # single-bundle example before any traffic is booked
    r = earliest_route(plan, "ROVER", "MOC", 0.0, traffic[0][3])
    print("First bundle, empty network:")
    for i, (t_tx, t_arr) in zip(r.contacts, r.hops):
        c = plan[i]
        print(f"  {c.tx:>5} -> {c.rx:<5} tx {t_tx/H:7.4f} h  arrive {t_arr/H:7.4f} h")
    print()

    cgr = route_and_book(plan, traffic, "MOC")
    print(summarise("CGR (volume-aware Dijkstra)", cgr, plan))
    util = [(c.tx, c.rx, c.start, c.booked_bits / ((c.end - c.start) * c.rate))
            for c in plan if c.rx != "MOC" and c.tx != "ROVER" and c.start < 24 * H]
    print("  DSN pass utilisation, day 1:")
    for tx, rx, s, u in util:
        print(f"    {tx} -> {rx} at {s/H:5.1f} h: {100*u:5.1f} %")
    print()

    plan_g = lunar_plan(horizon)
    print(summarise("Greedy next-hop", greedy_forward(plan_g, traffic, "MOC"), plan_g))
    print()

    plan_1 = lunar_plan(horizon)
    print(summarise("CGR, ORB2 out of service",
                    route_and_book(plan_1, traffic, "MOC", exclude={"ORB2"}), plan_1))


if __name__ == "__main__":
    main()
