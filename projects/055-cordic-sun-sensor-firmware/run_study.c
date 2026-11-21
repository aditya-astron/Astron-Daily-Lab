/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * 1. CORDIC accuracy vs iteration count, measured against libm and the bounds.
 * 2. Sun-sensor error budget: how many iterations the firmware needs before
 *    the CORDIC stops adding to the 12-bit ADC floor.
 */
#include <stdio.h>
#include "common.h"
#include "sun_sensor.h"

#define DEG (180.0 / M_PI)

static void accuracy_table(void)
{
    printf("CORDIC accuracy vs iterations (worst case over 65536 angles / 100000 vectors)\n");
    printf("  N   sincos err   bound      atan2 err    bound      |r| rel err  bound\n");
    const int ns[] = {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    for (size_t k = 0; k < sizeof ns / sizeof ns[0]; k++) {
        int n = ns[k];
        double ws = 0, wa = 0, wm = 0;
        for (uint32_t j = 0; j < 65536; j++) {
            uint32_t a = j * 65536u + (uint32_t)(rng_u64() & 0xFFFF);
            int32_t c, s;
            cordic_sincos(a, n, &c, &s);
            double e = hypot(c * 0x1.0p-30 - cos(bam_to_rad(a)), s * 0x1.0p-30 - sin(bam_to_rad(a)));
            if (e > ws) ws = e;
        }
        for (int j = 0; j < 100000; j++) {
            int32_t x = (int32_t)(rng_u64() >> 33) - (1 << 30);
            int32_t y = (int32_t)(rng_u64() >> 33) - (1 << 30);
            uint32_t a, m;
            cordic_atan2_mag(x, y, n, &a, &m);
            double ea = fabs(ang_err(bam_to_rad(a), atan2((double)y, (double)x)));
            double r = hypot((double)x, (double)y);
            double em = fabs(m - r) - 0.5;      /* integer output rounding excluded */
            em = em > 0 ? em / r : 0.0;
            if (ea > wa) wa = ea;
            if (em > wm) wm = em;
        }
        printf("  %2d  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n", n, ws, bound_rotation(n), wa,
               bound_vector_angle(n), wm, bound_vector_mag_rel(n));
    }
}

typedef struct { double sum2, worst; long count; } stat;

static void acc(stat *s, double e)
{
    s->sum2 += e * e;
    if (e > s->worst) s->worst = e;
    s->count++;
}

static double angle_between(double ax, double ay, double az, double bx, double by, double bz)
{
    double cx = ay * bz - az * by, cy = az * bx - ax * bz, cz = ax * by - ay * bx;
    return atan2(sqrt(cx * cx + cy * cy + cz * cz), ax * bx + ay * by + az * bz);
}

static uint16_t adc(double i, double sigma)
{
    double v = i + sigma * rng_gauss();
    if (v < 0) v = 0;
    if (v > 4095) v = 4095;
    return (uint16_t)lround(v);
}

static void sun_budget(double sigma)
{
    const double t0 = 30.0 / DEG, g = 4000.0, fov = 55.0 / DEG;
    const int ns[] = {6, 8, 10, 11, 12, 13, 14, 16, 20, 24};
    enum { NN = sizeof ns / sizeof ns[0], TRIALS = 200000 };
    stat floor_s = {0}, fw[NN] = {{0}};

    for (int k = 0; k < TRIALS; k++) {
        /* Uniform over the spherical cap out to the 55 deg field-of-view edge. */
        double off = acos(1.0 - rng_unit() * (1.0 - cos(fov)));
        double azt = 2.0 * M_PI * rng_unit();
        double sx = sin(off) * cos(azt), sy = sin(off) * sin(azt), sz = cos(off);
        double st = sin(t0), ct = cos(t0);
        uint16_t c[4] = {adc(g * (st * sx + ct * sz), sigma), adc(g * (-st * sx + ct * sz), sigma),
                         adc(g * (st * sy + ct * sz), sigma), adc(g * (-st * sy + ct * sz), sigma)};

        /* ADC floor: the same counts through double-precision math. */
        double dx = c[0] - (double)c[1], dy = c[2] - (double)c[3];
        double ssum = (double)c[0] + c[1] + c[2] + c[3];
        double ex = dx / (2 * st), ey = dy / (2 * st), ez = ssum / (4 * ct);
        acc(&floor_s, angle_between(ex, ey, ez, sx, sy, sz));

        for (int j = 0; j < NN; j++) {
            sun_sensor_cfg cfg = {(uint32_t)llround(2.0 / tan(t0) * 0x1.0p28), ns[j]};
            uint32_t a, o;
            sun_sensor_solve(&cfg, c, &a, &o);
            double oa = bam_to_rad(o), aa = bam_to_rad(a);
            acc(&fw[j], angle_between(sin(oa) * cos(aa), sin(oa) * sin(aa), cos(oa), sx, sy, sz));
        }
    }

    printf("\nSun sensor pointing error, t0 = 30 deg, G = 4000 counts, 12-bit ADC, sigma = %.1f count, FOV 55 deg\n", sigma);
    printf("  solver            RMS [deg]   max [deg]   RMS over floor\n");
    double fr = sqrt(floor_s.sum2 / (double)floor_s.count);
    printf("  double precision  %.5f     %.5f     -\n", fr * DEG, floor_s.worst * DEG);
    for (int j = 0; j < NN; j++) {
        double r = sqrt(fw[j].sum2 / (double)fw[j].count);
        printf("  CORDIC N=%-2d       %.5f     %.5f     %+6.2f %%\n", ns[j], r * DEG, fw[j].worst * DEG,
               100.0 * (r / fr - 1.0));
    }
}

int main(void)
{
    accuracy_table();
    sun_budget(0.0);
    sun_budget(1.0);
    return 0;
}
