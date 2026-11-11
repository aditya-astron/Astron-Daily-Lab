/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 *
 * Checks the analysis against hand calculations and the simulator against
 * the analysis: RTA is exact at the critical instant, tight under adversarial
 * blocking, and an upper bound under random phasing.
 */
#include <stdio.h>
#include <string.h>

#include "fsw_tasks.h"

static int fails;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static long hyper(const taskset_t *ts)
{
    long H = 1;
    for (int i = 0; i < ts->n; i++) H = lcm_l(H, ts->t[i].T);
    return H;
}

/* Classic three-task example (T,C) = (7,3), (12,3), (20,5): by hand
 * R3: 5 -> 11 -> 14 -> 17 -> 20 -> 20. */
static void test_textbook(void)
{
    taskset_t ts = {0};
    int T[3] = {7, 12, 20}, C[3] = {3, 3, 5};
    ts.n = 3;
    for (int i = 0; i < 3; i++)
        ts.t[i] = (task_t){"t", T[i], T[i], C[i], 0, 1, {{NO_RES, C[i]}}};
    CHECK(sched_setup(&ts) == 0, "textbook setup");
    long R[3] = {sched_rta(&ts, 0, 1), sched_rta(&ts, 1, 1), sched_rta(&ts, 2, 1)};
    CHECK(R[0] == 3 && R[1] == 6 && R[2] == 20, "textbook R = %ld %ld %ld, want 3 6 20", R[0], R[1], R[2]);
    /* one more tick on task 3 pushes R past 20 */
    ts.t[2].C = 6; ts.t[2].seg[0].len = 6;
    CHECK(sched_rta(&ts, 2, 1) < 0, "textbook C3=6 should be unschedulable");
}

static void test_setup(const taskset_t *ts)
{
    int want_prio[7] = {0, 1, 2, 3, 4, 5, 6};
    for (int i = 0; i < ts->n; i++) CHECK(ts->t[i].prio == want_prio[i], "DM prio %s", ts->t[i].name);
    CHECK(ts->ceil[R_STATE_DB] == 0 && ts->ceil[R_IMU_BUS] == 0 && ts->ceil[R_CAN_BUS] == 1, "ceilings");
    /* U = .18 + .125 + .18 + .06 + .15 + .095 + .10 */
    double U = sched_utilisation(ts);
    CHECK(U > 0.8899 && U < 0.8901, "U = %f", U);
    int who, res;
    CHECK(sched_blocking_icpp(ts, 0, &who, &res) == 20 && who == 6 && res == R_STATE_DB, "B aocs");
    CHECK(sched_blocking_icpp(ts, 6, &who, &res) == 0 && who == -1, "B lowest = 0");
    /* aocs: C + B = 18 + 20, no interference; rw: 25 + 20 + ceil(R/100)*18 = 63 */
    CHECK(sched_rta(ts, 0, 1) == 38, "R aocs %ld", sched_rta(ts, 0, 1));
    CHECK(sched_rta(ts, 1, 1) == 63, "R rw %ld", sched_rta(ts, 1, 1));
}

/* Synchronous release, nothing holding a lock: first job's response equals
 * RTA with B = 0 for every task (critical instant theorem, R <= D <= T). */
static void test_critical_instant(const taskset_t *ts)
{
    simcfg_t cfg;
    simres_t r;
    sim_default_cfg(ts, &cfg);
    for (int p = 0; p < 3; p++) {
        sim_run(ts, &cfg, (lockproto_t)p, hyper(ts), &r);
        for (int i = 0; i < ts->n; i++) {
            long R0 = sched_rta_b(ts, i, 1.0, 0);
            if (p == LOCK_ICPP)
                CHECK(r.first_resp[i] == R0, "sync %s sim %ld rta %ld", ts->t[i].name, r.first_resp[i], R0);
            CHECK(r.misses[i] == 0, "sync misses %s proto %d", ts->t[i].name, p);
        }
        CHECK(!r.overrun, "sync overrun");
    }
    CHECK(r.icpp_lock_conflicts == 0, "sync conflicts");
}

/* Adversarial phasing for each task i: the worst blocker enters its critical
 * section at t = 0, task i and everything above it is released at t = 1.
 * In integer ticks the blocker has L-1 ticks left, so R_sim = RTA(B = L-1). */
static void test_blocking_tight(const taskset_t *ts)
{
    for (int i = 0; i < ts->n; i++) {
        int who, res;
        int B = sched_blocking_icpp(ts, i, &who, &res);
        if (!B) continue;
        simcfg_t cfg;
        simres_t r;
        sim_default_cfg(ts, &cfg);
        for (int k = 0; k < ts->n; k++) cfg.offset[k] = ts->t[k].prio <= ts->t[i].prio ? 1 : -1;
        cfg.offset[who] = 0;
        cfg.only_first[who] = 1;
        cfg.first_nseg[who] = 1;
        cfg.first_seg[who][0] = (seg_t){res, B};
        sim_run(ts, &cfg, LOCK_ICPP, hyper(ts), &r);
        long want = sched_rta_b(ts, i, 1.0, B - 1);
        CHECK(r.first_resp[i] == want, "blocked %s sim %ld rta(B-1) %ld", ts->t[i].name, r.first_resp[i], want);
        CHECK(r.first_resp[i] <= sched_rta(ts, i, 1.0), "blocked %s above bound", ts->t[i].name);
        CHECK(r.icpp_lock_conflicts == 0, "blocked conflicts");
    }
}

/* Random phasing and execution times never exceed the ICPP bound. */
static void test_random(const taskset_t *ts)
{
    uint64_t s = 7;
    long worst[MAX_TASKS] = {0};
    for (int run = 0; run < 300; run++) {
        simcfg_t cfg;
        simres_t r;
        sim_default_cfg(ts, &cfg);
        for (int i = 0; i < ts->n; i++) {
            s = s * 6364136223846793005ull + 1442695040888963407ull;
            cfg.offset[i] = (int)((s >> 33) % ts->t[i].T);
        }
        cfg.exec_frac_min = (run % 3 == 0) ? 1.0 : 0.2;
        cfg.seed = s | 1;
        sim_run(ts, &cfg, LOCK_ICPP, 2 * hyper(ts), &r);
        CHECK(r.icpp_lock_conflicts == 0, "random conflicts run %d", run);
        for (int i = 0; i < ts->n; i++) {
            CHECK(r.misses[i] == 0, "random miss %s run %d", ts->t[i].name, run);
            if (r.max_resp[i] > worst[i]) worst[i] = r.max_resp[i];
        }
    }
    for (int i = 0; i < ts->n; i++)
        CHECK(worst[i] <= sched_rta(ts, i, 1.0), "random %s %ld > bound", ts->t[i].name, worst[i]);
}

/* Logger holds state_db, everything else released one tick later.
 * ICPP and inheritance: aocs = 19 left + 18 = 37 ticks.
 * Plain mutex: rw (25), nav up to its state_db section (80) and thermal up to
 * its section (170) run first, then logger's 19, then aocs 18:
 * 25 + 80 + 170 + 19 + 18 = 312 ticks. */
static void test_inversion(const taskset_t *ts)
{
    long want[3] = {37, 37, 312};
    for (int p = 0; p < 3; p++) {
        simcfg_t cfg;
        simres_t r;
        sim_default_cfg(ts, &cfg);
        for (int i = 0; i < ts->n; i++) { cfg.offset[i] = 1; cfg.only_first[i] = 1; }
        cfg.offset[6] = 0;
        cfg.first_nseg[6] = 1;
        cfg.first_seg[6][0] = (seg_t){R_STATE_DB, 20};
        sim_run(ts, &cfg, (lockproto_t)p, hyper(ts), &r);
        CHECK(r.first_resp[0] == want[p], "inversion proto %d aocs %ld want %ld", p, r.first_resp[0], want[p]);
    }
}

/* Breakdown scale is a true boundary of the exact test. */
static void test_breakdown(const taskset_t *ts)
{
    double s = sched_breakdown(ts);
    int ok_at = 1, ok_above = 1;
    for (int i = 0; i < ts->n; i++) {
        if (sched_rta(ts, i, s) < 0) ok_at = 0;
        if (sched_rta(ts, i, s * (1 + 1e-6) + 1e-9) < 0) ok_above = 0;
    }
    CHECK(s > 1.0 && s < 1.0 / sched_utilisation(ts) + 1e-9, "breakdown %f out of range", s);
    CHECK(ok_at && !ok_above, "breakdown not a boundary");
}

int main(void)
{
    taskset_t ts;
    fsw_taskset(&ts);
    CHECK(sched_setup(&ts) == 0, "setup");
    test_textbook();
    test_setup(&ts);
    test_critical_instant(&ts);
    test_blocking_tight(&ts);
    test_random(&ts);
    test_inversion(&ts);
    test_breakdown(&ts);
    printf(fails ? "%d check(s) failed\n" : "all checks passed\n", fails);
    return fails != 0;
}
