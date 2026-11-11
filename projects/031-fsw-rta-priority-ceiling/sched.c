/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 */
#include "sched.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

long lcm_l(long a, long b)
{
    long x = a, y = b;
    while (y) { long r = x % y; x = y; y = r; }
    return a / x * b;
}

int sched_setup(taskset_t *ts)
{
    int order[MAX_TASKS];
    for (int i = 0; i < ts->n; i++) {
        task_t *t = &ts->t[i];
        int sum = 0;
        for (int s = 0; s < t->nseg; s++) {
            if (t->seg[s].len <= 0) return -1;
            if (t->seg[s].res >= ts->nres) return -1;
            sum += t->seg[s].len;
        }
        if (sum != t->C || t->D > t->T || t->C > t->D) return -1;
        order[i] = i;
    }
    /* Deadline-monotonic: optimal fixed-priority order for D <= T. */
    for (int a = 1; a < ts->n; a++)
        for (int b = a; b > 0; b--) {
            const task_t *x = &ts->t[order[b - 1]], *y = &ts->t[order[b]];
            if (x->D > y->D || (x->D == y->D && x->T > y->T)) {
                int tmp = order[b]; order[b] = order[b - 1]; order[b - 1] = tmp;
            } else break;
        }
    for (int k = 0; k < ts->n; k++) ts->t[order[k]].prio = k;

    /* Ceiling of a resource = priority of its highest-priority user. */
    for (int r = 0; r < ts->nres; r++) ts->ceil[r] = MAX_TASKS;
    for (int i = 0; i < ts->n; i++)
        for (int s = 0; s < ts->t[i].nseg; s++) {
            int r = ts->t[i].seg[s].res;
            if (r != NO_RES && ts->t[i].prio < ts->ceil[r]) ts->ceil[r] = ts->t[i].prio;
        }
    return 0;
}

double sched_utilisation(const taskset_t *ts)
{
    double u = 0;
    for (int i = 0; i < ts->n; i++) u += (double)ts->t[i].C / ts->t[i].T;
    return u;
}

static int scaled(int x, double s) { return (int)ceil(x * s - 1e-9); }

static int blocking(const taskset_t *ts, int i, double s, int *who, int *res)
{
    int B = 0, pi = ts->t[i].prio;
    if (who) *who = -1;
    if (res) *res = NO_RES;
    for (int j = 0; j < ts->n; j++) {
        if (ts->t[j].prio <= pi) continue;
        for (int k = 0; k < ts->t[j].nseg; k++) {
            int r = ts->t[j].seg[k].res;
            if (r == NO_RES || ts->ceil[r] > pi) continue;
            int L = scaled(ts->t[j].seg[k].len, s);
            if (L > B) { B = L; if (who) *who = j; if (res) *res = r; }
        }
    }
    return B;
}

int sched_blocking_icpp(const taskset_t *ts, int i, int *who, int *res)
{
    return blocking(ts, i, 1.0, who, res);
}

long sched_rta_b(const taskset_t *ts, int i, double s, long B)
{
    const task_t *ti = &ts->t[i];
    long R = scaled(ti->C, s) + B, prev = -1;
    while (R != prev) {
        if (R > ti->D) return -1;
        prev = R;
        R = scaled(ti->C, s) + B;
        for (int j = 0; j < ts->n; j++)
            if (ts->t[j].prio < ti->prio)
                R += ((prev + ts->t[j].T - 1) / ts->t[j].T) * scaled(ts->t[j].C, s);
    }
    return R;
}

long sched_rta(const taskset_t *ts, int i, double s)
{
    return sched_rta_b(ts, i, s, blocking(ts, i, s, NULL, NULL));
}

static int all_ok(const taskset_t *ts, double s)
{
    for (int i = 0; i < ts->n; i++)
        if (sched_rta(ts, i, s) < 0) return 0;
    return 1;
}

double sched_breakdown(const taskset_t *ts)
{
    double lo = 0.0, hi = 1.0;
    if (!all_ok(ts, 1e-6)) return 0.0;
    while (all_ok(ts, hi)) { lo = hi; hi *= 2; }
    for (int it = 0; it < 60; it++) {
        double m = 0.5 * (lo + hi);
        if (all_ok(ts, m)) lo = m; else hi = m;
    }
    return lo;
}

/* ------------------------------ simulator ------------------------------ */

typedef struct {
    int task, nseg, seg, rem, holding;
    long release;
    seg_t s[MAX_SEGS];
} job_t;

#define POOL 512

static uint64_t xs(uint64_t *st)
{
    uint64_t x = *st;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return *st = x;
}

void sim_default_cfg(const taskset_t *ts, simcfg_t *cfg)
{
    memset(cfg, 0, sizeof *cfg);
    (void)ts;
    cfg->exec_frac_min = 1.0;
    cfg->seed = 0x9E3779B97F4A7C15ull;
}

static int need_res(const job_t *j)
{
    int r = j->s[j->seg].res;
    return (r != NO_RES && !j->holding) ? r : NO_RES;
}

void sim_run(const taskset_t *ts, const simcfg_t *cfg, lockproto_t proto,
             long horizon, simres_t *out)
{
    static job_t pool[POOL];
    int active[POOL], na = 0, nfree = POOL, freel[POOL];
    int owner[MAX_RES];
    uint64_t rng = cfg->seed ? cfg->seed : 1;

    memset(out, 0, sizeof *out);
    for (int i = 0; i < MAX_TASKS; i++) out->first_resp[i] = -1;
    for (int r = 0; r < MAX_RES; r++) owner[r] = -1;
    for (int k = 0; k < POOL; k++) freel[k] = POOL - 1 - k;

    for (long t = 0; t < horizon; t++) {
        /* releases */
        for (int i = 0; i < ts->n; i++) {
            const task_t *tk = &ts->t[i];
            int off = cfg->offset[i];
            if (off < 0 || t < off || (t - off) % tk->T) continue;
            int first = (t == off);
            if (!first && cfg->only_first[i]) continue;
            if (!nfree) abort();
            int id = freel[--nfree];
            job_t *j = &pool[id];
            memset(j, 0, sizeof *j);
            j->task = i;
            j->release = t;
            if (first && cfg->first_nseg[i]) {
                j->nseg = cfg->first_nseg[i];
                memcpy(j->s, cfg->first_seg[i], sizeof(seg_t) * j->nseg);
            } else {
                j->nseg = tk->nseg;
                memcpy(j->s, tk->seg, sizeof(seg_t) * j->nseg);
                for (int k = 0; k < j->nseg; k++) {
                    if (j->s[k].res != NO_RES || cfg->exec_frac_min >= 1.0) continue;
                    double u = (double)(xs(&rng) >> 11) / 9007199254740992.0;
                    double f = cfg->exec_frac_min + (1.0 - cfg->exec_frac_min) * u;
                    int L = (int)ceil(j->s[k].len * f);
                    j->s[k].len = L < 1 ? 1 : L;
                }
            }
            j->rem = j->s[0].len;
            active[na++] = id;
        }
        if (!na) continue;

        /* eligibility and effective priority */
        int best = -1, bestp = 1 << 30;
        for (int a = 0; a < na; a++) {
            job_t *j = &pool[active[a]];
            int base = ts->t[j->task].prio, eff = base;
            int r = need_res(j);
            /* Under ICPP a job never waits on a mutex; whether the chosen job
             * finds its lock free is checked after selection instead. */
            if (proto != LOCK_ICPP && r != NO_RES && owner[r] >= 0 && owner[r] != active[a])
                continue; /* blocked on the mutex */
            if (j->holding) {
                int hr = j->s[j->seg].res;
                if (proto == LOCK_ICPP) {
                    if (ts->ceil[hr] < eff) eff = ts->ceil[hr];
                } else if (proto == LOCK_PIP) {
                    /* inherit from every job waiting on the held lock */
                    for (int b = 0; b < na; b++) {
                        job_t *w = &pool[active[b]];
                        if (need_res(w) == hr && ts->t[w->task].prio < eff)
                            eff = ts->t[w->task].prio;
                    }
                }
            }
            if (eff < bestp || (eff == bestp && j->release < pool[active[best]].release)) {
                bestp = eff;
                best = a;
            }
        }
        if (best < 0) continue; /* every job blocked: cannot happen without nesting */

        /* run the chosen job for one tick */
        int id = active[best];
        job_t *j = &pool[id];
        int r = need_res(j);
        if (r != NO_RES && owner[r] >= 0 && owner[r] != id) {
            out->icpp_lock_conflicts++; /* ICPP guarantee violated: lose the tick */
            continue;
        }
        if (r != NO_RES) { owner[r] = id; j->holding = 1; }
        if (--j->rem == 0) {
            if (j->holding) { owner[j->s[j->seg].res] = -1; j->holding = 0; }
            if (++j->seg < j->nseg) {
                j->rem = j->s[j->seg].len;
            } else {
                const task_t *tk = &ts->t[j->task];
                long resp = t + 1 - j->release;
                if (out->first_resp[j->task] < 0) out->first_resp[j->task] = resp;
                if (resp > out->max_resp[j->task]) out->max_resp[j->task] = resp;
                if (resp > tk->D) out->misses[j->task]++;
                out->jobs[j->task]++;
                active[best] = active[--na];
                freel[nfree++] = id;
            }
        }
    }
    for (int a = 0; a < na; a++) {
        const job_t *j = &pool[active[a]];
        out->overrun = 1;
        if (horizon - j->release > ts->t[j->task].D) out->misses[j->task]++;
    }
}
