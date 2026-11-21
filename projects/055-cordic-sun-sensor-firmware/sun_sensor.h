/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * Four-cell coarse sun sensor. Cell normals are tilted by t0 from the sensor
 * boresight (+z) toward +x, -x, +y, -y:
 *     n_{+-x} = (+-sin t0, 0, cos t0),  n_{+-y} = (0, +-sin t0, cos t0)
 * For cosine-law cells all lit, I_k = G * n_k . s, so
 *     Dx = I+x - I-x = 2G sin t0 s_x,   Dy = 2G sin t0 s_y
 *     S  = sum I_k   = 4G cos t0 s_z
 * azimuth   = atan2(Dy, Dx)
 * off-axis  = atan2(c*|D|, S),  c = 2 cot t0
 * The gain G (sun intensity, cell ageing, distance) cancels in both ratios.
 */
#ifndef SUN_SENSOR_H
#define SUN_SENSOR_H

#include <stdint.h>

typedef struct {
    uint32_t c_q28;   /* 2 cot(t0) in Q4.28 */
    int n_iter;       /* CORDIC iterations */
} sun_sensor_cfg;

/* counts: {+x, -x, +y, -y} ADC readings. Outputs in BAM. */
void sun_sensor_solve(const sun_sensor_cfg *cfg, const uint16_t counts[4],
                      uint32_t *azimuth, uint32_t *off_axis);

#endif
