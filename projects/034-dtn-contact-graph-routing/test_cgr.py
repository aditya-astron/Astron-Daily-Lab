# Astron Built this AP
# Project 034 | Algorithmic problem solving | Contact graph routing for a lunar far-side DTN relay
# Author: Aditya (aditya-astron)
# Date: 2026-11-12 IST
"""Verification for cgr.py. Run: python3 test_cgr.py"""

import math
import random

from cgr import (Contact, H, OWLT_GROUND, OWLT_LOCAL, OWLT_MOON, earliest_route,
                 greedy_forward, lunar_plan, route_and_book, rover_traffic)


def brute_force(plan, src, dst, t0, bits):
    """Reference earliest arrival: enumerate every loop-free contact sequence
    by DFS, with the hop arithmetic written out independently."""
    best = math.inf

    def dfs(node, t, visited):
        nonlocal best
        for c in plan:
            if c.tx != node or c.rx in visited:
                continue
            start = max(t, c.start, c.free_at)
            if start + bits / c.rate > c.end + 1e-9:
                continue
            arr = start + bits / c.rate + c.owlt
            if c.rx == dst:
                best = min(best, arr)
            else:
                dfs(c.rx, arr, visited | {c.rx})

    dfs(src, t0, {src})
    return best


def random_plan(rng, n_nodes=6, n_contacts=22, horizon=100.0):
    nodes = [f"N{i}" for i in range(n_nodes)]
    plan = []
    for _ in range(n_contacts):
        a, b = rng.sample(nodes, 2)
        s = rng.uniform(0, horizon)
        c = Contact(a, b, s, s + rng.uniform(1, 25), rng.choice([0.5, 1, 2, 5]),
                    rng.uniform(0, 3))
        c.free_at = s + rng.uniform(0, 10) if rng.random() < 0.3 else 0.0
        plan.append(c)
    return plan, nodes


def test_hand_example():
    plan = lunar_plan(48 * H)
    bits = 1.2e9
    r = earliest_route(plan, "ROVER", "MOC", 0.0, bits)
    assert r.path(plan) == ["ROVER", "ORB1", "DSS14", "MOC"]
    # ORB1 window opens at 1 h; 600 s on the 2 Mbit/s proximity link
    t1 = 1 * H + bits / 2e6 + OWLT_LOCAL
    # wait for the 5.5 h Goldstone Ka pass; 60 s at 20 Mbit/s
    t2 = 5.5 * H + bits / 20e6 + OWLT_MOON
    # terrestrial hop, 24 s at 50 Mbit/s
    t3 = t2 + bits / 50e6 + OWLT_GROUND
    assert abs(r.hops[0][1] - t1) < 1e-9
    assert abs(r.hops[1][1] - t2) < 1e-9
    assert abs(r.arrival - t3) < 1e-9
    assert abs(t3 - 19885.33) < 1e-6   # 5.5237 h


def test_dijkstra_matches_exhaustive():
    rng = random.Random(34)
    checked = reachable = 0
    for _ in range(400):
        plan, nodes = random_plan(rng)
        src, dst = rng.sample(nodes, 2)
        t0, bits = rng.uniform(0, 50), rng.choice([1.0, 4.0, 10.0])
        ref = brute_force(plan, src, dst, t0, bits)
        r = earliest_route(plan, src, dst, t0, bits)
        if math.isinf(ref):
            assert r is None
        else:
            reachable += 1
            assert r is not None and abs(r.arrival - ref) < 1e-9, (ref, r)
            p = r.path(plan)
            assert p[0] == src and p[-1] == dst and len(set(p)) == len(p)
        checked += 1
    assert reachable > 100   # make sure the test exercises real routes
    return checked, reachable


def test_fifo_property():
    rng = random.Random(7)
    for _ in range(2000):
        c = Contact("A", "B", rng.uniform(0, 10), rng.uniform(10, 30), 2.0, 1.0,
                    free_at=rng.uniform(0, 15))
        bits = rng.uniform(0.1, 10)
        ta, tb = sorted(rng.uniform(0, 30) for _ in range(2))
        ha, hb = c.hop(ta, bits), c.hop(tb, bits)
        if hb is not None:
            assert ha is not None and ha[1] <= hb[1] + 1e-12


def check_bookings(plan, results):
    """Replay every delivered bundle hop by hop against the contact windows."""
    for c in plan:
        last = -math.inf
        for _, t_tx, t_end in c.bookings:          # booking order = FIFO order
            assert t_tx >= c.start - 1e-9 and t_end <= c.end + 1e-9
            assert t_tx >= last - 1e-9               # no overlap on a contact
            last = t_end
        assert c.booked_bits <= (c.end - c.start) * c.rate + 1e-3
    for bid, src, t0, bits, r in results:
        if r is None:
            continue
        t = t0
        for i, (t_tx, t_arr) in zip(r.contacts, r.hops):
            c = plan[i]
            assert t_tx >= t - 1e-9                  # bundle is held before it is sent
            assert any(b == bid and abs(s - t_tx) < 1e-9 for b, s, _ in c.bookings)
            t = t_tx + bits / c.rate + c.owlt
            assert abs(t - t_arr) < 1e-9
        assert abs(t - r.arrival) < 1e-9


def test_scenario_capacity_and_routing():
    horizon = 48 * H
    traffic = rover_traffic(horizon)

    plan = lunar_plan(horizon)
    cgr = route_and_book(plan, traffic, "MOC")
    check_bookings(plan, cgr)
    assert all(r is not None for *_, r in cgr)

    plan_g = lunar_plan(horizon)
    greedy = greedy_forward(plan_g, traffic, "MOC")
    check_bookings(plan_g, greedy)

    plan_x = lunar_plan(horizon)
    no_orb2 = route_and_book(plan_x, traffic, "MOC", exclude={"ORB2"})
    check_bookings(plan_x, no_orb2)
    assert all("ORB2" not in r.path(plan_x) for *_, r in no_orb2 if r)

    def mean_lat(res):
        lat = [r.arrival - t0 for _, _, t0, _, r in res if r]
        return sum(lat) / len(lat), len(lat)

    m_cgr, n_cgr = mean_lat(cgr)
    m_g, n_g = mean_lat(greedy)
    m_x, n_x = mean_lat(no_orb2)
    # greedy floods ORB2, whose single X-band pass holds 12 products a day
    assert n_g < n_cgr and m_g > m_cgr
    # losing ORB2 must cost something
    assert m_x > m_cgr
    # the ORB2 downlink is a hard bottleneck: 1 h x 4 Mbit/s / 1.2 Gbit = 12
    orb2_dl = [c for c in plan if c.tx == "ORB2" and c.rx == "DSS63"]
    assert [len(c.bookings) for c in orb2_dl] == [12, 0]
    return (n_cgr, m_cgr / H), (n_g, m_g / H), (n_x, m_x / H)


if __name__ == "__main__":
    test_hand_example()
    print("hand example (1.2 Gbit, empty network): 5.5237 h  ok")
    n, k = test_dijkstra_matches_exhaustive()
    print(f"Dijkstra vs exhaustive DFS: {n} random plans, {k} reachable, all equal  ok")
    test_fifo_property()
    print("FIFO (non-overtaking) hop property: 2000 cases  ok")
    a, b, c = test_scenario_capacity_and_routing()
    print(f"scenario replay ok | CGR {a[0]}/72 mean {a[1]:.2f} h | "
          f"greedy {b[0]}/72 mean {b[1]:.2f} h | no ORB2 {c[0]}/72 mean {c[1]:.2f} h")
    print("all tests passed")
