/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 */
#ifndef SCHED_H
#define SCHED_H

#include <stdint.h>

#define MAX_TASKS 16
#define MAX_SEGS 8
#define MAX_RES 8
#define NO_RES (-1)

/* One stretch of a job: either plain code (res == NO_RES) or a critical
 * section on shared resource `res`. Sections are not nested. */
typedef struct {
    int res;
    int len;
} seg_t;

/* Periodic task. Times are integer ticks. prio: 0 is highest, all distinct. */
typedef struct {
    const char *name;
    int T, D, C;
    int prio;
    int nseg;
    seg_t seg[MAX_SEGS];
} task_t;

typedef struct {
    int n, nres;
    task_t t[MAX_TASKS];
    int ceil[MAX_RES]; /* filled by sched_setup */
} taskset_t;

typedef enum { LOCK_ICPP, LOCK_PIP, LOCK_PLAIN } lockproto_t;

/* Validates segment sums, assigns deadline-monotonic priorities and
 * resource ceilings. Returns 0 on success. */
int sched_setup(taskset_t *ts);

double sched_utilisation(const taskset_t *ts);

/* Worst-case blocking under the immediate ceiling protocol: the longest
 * critical section of a lower-priority task on a resource whose ceiling is
 * at least as high as task i's priority. Also returns which (task, res). */
int sched_blocking_icpp(const taskset_t *ts, int i, int *who, int *res);

/* Exact response-time analysis R = C + B + sum_hp ceil(R/Tj) Cj.
 * cscale scales every WCET (for breakdown search); B from ICPP. Returns R,
 * or -1 if the iteration passes D. */
long sched_rta(const taskset_t *ts, int i, double cscale);

/* Same with an explicit blocking term B (unscaled ticks). */
long sched_rta_b(const taskset_t *ts, int i, double cscale, long B);

/* Largest uniform WCET scale for which every task still meets D (bisection). */
double sched_breakdown(const taskset_t *ts);

/* ---- tick-level simulator of a preemptive fixed-priority kernel ---- */

typedef struct {
    int offset[MAX_TASKS];   /* first release; <0 means the task never runs */
    int only_first[MAX_TASKS]; /* release a single job only */
    /* optional override of the first job's segments (adversarial phasing) */
    int first_nseg[MAX_TASKS];
    seg_t first_seg[MAX_TASKS][MAX_SEGS];
    double exec_frac_min;    /* plain segments run U(min,1)*len; 1 = WCET */
    uint64_t seed;
} simcfg_t;

typedef struct {
    long max_resp[MAX_TASKS];
    long first_resp[MAX_TASKS];
    long jobs[MAX_TASKS];
    long misses[MAX_TASKS];
    long icpp_lock_conflicts; /* must stay 0: ICPP never meets a held lock */
    int overrun;              /* job still pending at horizon end */
} simres_t;

void sim_default_cfg(const taskset_t *ts, simcfg_t *cfg);
void sim_run(const taskset_t *ts, const simcfg_t *cfg, lockproto_t proto,
             long horizon, simres_t *out);

long lcm_l(long a, long b);

#endif
