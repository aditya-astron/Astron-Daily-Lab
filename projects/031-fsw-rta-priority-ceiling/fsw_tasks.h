/*
 * Astron Built this AP
 * Project 031 | Embedded firmware | Flight software schedulability: response-time analysis with priority ceiling
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-11 IST
 */
#ifndef FSW_TASKS_H
#define FSW_TASKS_H
#include "sched.h"

/* 1 tick = 100 us. Resources shared between tasks: */
enum { R_STATE_DB, R_IMU_BUS, R_CAN_BUS, N_RES };

void fsw_taskset(taskset_t *ts);
extern const char *const fsw_res_name[N_RES];

#endif
