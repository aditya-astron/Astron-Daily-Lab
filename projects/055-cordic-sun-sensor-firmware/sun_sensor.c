/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 */
#include "sun_sensor.h"
#include "cordic.h"

/*
 * Two chained vectoring passes. The first gives azimuth and |D|; |D| feeds the
 * second pass directly, so there is no square root, division or libm call.
 * Counts are lifted to Q12 so |D| keeps 12 fractional bits between passes.
 */
void sun_sensor_solve(const sun_sensor_cfg *cfg, const uint16_t counts[4],
                      uint32_t *azimuth, uint32_t *off_axis)
{
    /* Multiply rather than shift: left-shifting a negative int is UB in C. */
    int32_t dx = ((int32_t)counts[0] - (int32_t)counts[1]) * 4096;
    int32_t dy = ((int32_t)counts[2] - (int32_t)counts[3]) * 4096;
    int32_t s = ((int32_t)counts[0] + counts[1] + counts[2] + counts[3]) * 4096;

    uint32_t d_mag, unused;
    cordic_atan2_mag(dx, dy, cfg->n_iter, azimuth, &d_mag);

    /* |D| <= sqrt(2)*4095*2^12 < 2^25 and c < 16, so c|D| < 2^29 fits int32. */
    int32_t cd = (int32_t)(((uint64_t)d_mag * cfg->c_q28 + (1u << 27)) >> 28);
    cordic_atan2_mag(s, cd, cfg->n_iter, off_axis, &unused);
}
