/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 */
/*
 * Schedulability report for the flight task set: utilisation bounds, ICPP
 * blocking terms, exact response times, breakdown scale, and a simulated
 * priority-inversion case comparing a plain mutex, priority inheritance
 * and the immediate ceiling protocol.
 */
#include <math.h>
#include <stdio.h>

#include "fsw_tasks.h"

#define TICK_MS 0.1

int main(void)
{
    taskset_t ts;
    fsw_taskset(&ts);
    if (sched_setup(&ts)) { fprintf(stderr, "bad task set\n"); return 1; }

    int n = ts.n;
    double U = sched_utilisation(&ts);
    double ll = n * (pow(2.0, 1.0 / n) - 1.0), hyp = 1.0;
    for (int i = 0; i < n; i++) hyp *= 1.0 + (double)ts.t[i].C / ts.t[i].T;
    long H = 1;
    for (int i = 0; i < n; i++) H = lcm_l(H, ts.t[i].T);

    printf("Tasks %d, utilisation U = %.4f, hyperperiod %.0f ms\n", n, U, H * TICK_MS);
    printf("Liu-Layland bound n(2^(1/n)-1) = %.4f -> %s\n", ll, U <= ll ? "pass" : "inconclusive");
    printf("Hyperbolic bound prod(1+Ui) = %.4f (<= 2 needed) -> %s\n", hyp, hyp <= 2 ? "pass" : "inconclusive");
    printf("Resource ceilings:");
    for (int r = 0; r < ts.nres; r++) printf(" %s=P%d", fsw_res_name[r], ts.ceil[r]);
    printf("\n\n%-13s %4s %6s %6s %6s %6s %8s %-18s %6s\n", "task", "prio", "T ms", "D ms",
           "C ms", "B ms", "R ms", "blocked by", "slack");

    simcfg_t cfg;
    sim_default_cfg(&ts, &cfg);
    simres_t sync;
    sim_run(&ts, &cfg, LOCK_ICPP, H, &sync);

    for (int p = 0; p < n; p++) {
        int i = 0;
        while (ts.t[i].prio != p) i++;
        const task_t *t = &ts.t[i];
        int who, res;
        int B = sched_blocking_icpp(&ts, i, &who, &res);
        long R = sched_rta(&ts, i, 1.0);
        char by[40] = "-";
        if (B) snprintf(by, sizeof by, "%s/%s", ts.t[who].name, fsw_res_name[res]);
        printf("%-13s %4d %6.1f %6.1f %6.1f %6.1f %8.1f %-18s %5.0f%%\n", t->name, t->prio,
               t->T * TICK_MS, t->D * TICK_MS, t->C * TICK_MS, B * TICK_MS, R * TICK_MS, by,
               R < 0 ? 0.0 : 100.0 * (t->D - R) / t->D);
    }

    double s = sched_breakdown(&ts);
    printf("\nBreakdown WCET scale: %.4f (every WCET can grow %.1f%% before a deadline is missed)\n",
           s, 100 * (s - 1));

    /* Monte Carlo: random phasing and execution times, ICPP. */
    long worst[MAX_TASKS] = {0}, misses = 0, conflicts = 0;
    uint64_t seed = 12345;
    for (int run = 0; run < 200; run++) {
        sim_default_cfg(&ts, &cfg);
        for (int i = 0; i < n; i++) {
            seed = seed * 6364136223846793005ull + 1442695040888963407ull;
            cfg.offset[i] = (int)((seed >> 33) % ts.t[i].T);
        }
        cfg.exec_frac_min = 0.3;
        cfg.seed = seed | 1;
        simres_t r;
        sim_run(&ts, &cfg, LOCK_ICPP, 2 * H, &r);
        for (int i = 0; i < n; i++) {
            if (r.max_resp[i] > worst[i]) worst[i] = r.max_resp[i];
            misses += r.misses[i];
        }
        conflicts += r.icpp_lock_conflicts;
    }
    printf("\nSimulated (200 random phasings x 2 hyperperiods, ICPP): misses %ld, lock conflicts %ld\n",
           misses, conflicts);
    for (int i = 0; i < n; i++)
        printf("  %-13s worst observed %6.1f ms   bound %6.1f ms\n", ts.t[i].name,
               worst[i] * TICK_MS, sched_rta(&ts, i, 1.0) * TICK_MS);

    /* Priority inversion: hk_logger takes state_db at t=0; aocs_ctrl and the
     * medium-priority tasks are released one tick later. */
    printf("\nInversion case: hk_logger holds state_db at t=0, everything else released once at t=1\n");
    const char *pname[] = {"ICPP", "inheritance", "plain mutex"};
    for (int pr = 0; pr < 3; pr++) {
        sim_default_cfg(&ts, &cfg);
        for (int i = 0; i < n; i++) { cfg.offset[i] = 1; cfg.only_first[i] = 1; }
        int lg = n - 1;
        cfg.offset[lg] = 0;
        cfg.first_nseg[lg] = 1;
        cfg.first_seg[lg][0] = (seg_t){R_STATE_DB, ts.t[lg].seg[0].len};
        simres_t r;
        sim_run(&ts, &cfg, (lockproto_t)pr, H, &r);
        printf("  %-12s aocs_ctrl response %6.1f ms (D = %.1f ms)%s\n", pname[pr],
               r.first_resp[0] * TICK_MS, ts.t[0].D * TICK_MS,
               r.first_resp[0] > ts.t[0].D ? "  DEADLINE MISSED" : "");
    }
    return 0;
}
