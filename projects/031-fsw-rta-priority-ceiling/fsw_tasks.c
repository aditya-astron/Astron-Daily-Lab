/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 */
#include "fsw_tasks.h"

const char *const fsw_res_name[N_RES] = {"state_db", "imu_bus", "can_bus"};

#define P(n) {NO_RES, n}
#define S(r, n) {r, n}

/* A small-satellite flight software build: 10 ms attitude loop down to a
 * 500 ms housekeeping logger. WCETs are measured worst cases in ticks;
 * critical sections are the mutex-protected stretches of each job. */
void fsw_taskset(taskset_t *ts)
{
    static const task_t tasks[] = {
        /* name             T     D     C   prio nseg  segments */
        {"aocs_ctrl",     100,  100,   18, 0, 3, {S(R_STATE_DB, 3), S(R_IMU_BUS, 4), P(11)}},
        {"rw_cmd",        200,  150,   25, 0, 2, {S(R_CAN_BUS, 4), P(21)}},
        {"nav_filter",    500,  500,   90, 0, 3, {S(R_IMU_BUS, 6), P(74), S(R_STATE_DB, 10)}},
        {"fdir",         1000,  800,   60, 0, 2, {S(R_STATE_DB, 5), P(55)}},
        {"tm_packetize", 1000, 1000,  150, 0, 3, {S(R_STATE_DB, 12), P(130), S(R_CAN_BUS, 8)}},
        {"thermal_ctrl", 2000, 2000,  190, 0, 3, {P(170), S(R_STATE_DB, 15), S(R_CAN_BUS, 5)}},
        {"hk_logger",    5000, 5000,  500, 0, 2, {S(R_STATE_DB, 20), P(480)}},
    };
    ts->n = (int)(sizeof tasks / sizeof tasks[0]);
    ts->nres = N_RES;
    for (int i = 0; i < ts->n; i++) ts->t[i] = tasks[i];
}
