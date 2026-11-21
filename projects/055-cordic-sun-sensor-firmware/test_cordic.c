/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * Checks the integer CORDIC against libm and against the analytic error
 * bounds in common.h, then checks the sun-sensor solver on hand cases.
 */
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "sun_sensor.h"

static int failures = 0;

#define CHECK(cond, ...) do { if (!(cond)) { failures++; \
    printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static void test_tables(void)
{
    double k = 1.0;
    for (int i = 0; i < CORDIC_MAX_ITER; i++) {
        long a = lround(atan(ldexp(1.0, -i)) / (2.0 * M_PI) * BAM_PER_TURN);
        CHECK(a == cordic_atan_bam[i], "atan table[%d] %ld vs %d", i, a, cordic_atan_bam[i]);
        k *= sqrt(1.0 + ldexp(1.0, -2 * i));
        long g = lround(ldexp(1.0, 31) / k);
        CHECK(g == cordic_inv_gain_q31[i], "gain table[%d] %ld vs %d", i, g, cordic_inv_gain_q31[i]);
    }
    CHECK(fabs(k - 1.6467602581210656) < 1e-15, "K_inf %.16f", k);
    CHECK((-7 >> 1) == -4, "arithmetic right shift assumed");
}

static void test_rotation(void)
{
    double prev = 0.0;
    for (int n = 4; n <= CORDIC_MAX_ITER; n++) {
        double worst = 0.0;
        for (uint32_t j = 0; j < 65536; j++) {
            uint32_t a = j * 65536u + (uint32_t)(rng_u64() & 0xFFFF);
            int32_t c, s;
            cordic_sincos(a, n, &c, &s);
            double th = bam_to_rad(a);
            double e = hypot(c * 0x1.0p-30 - cos(th), s * 0x1.0p-30 - sin(th));
            if (e > worst) worst = e;
        }
        double b = bound_rotation(n);
        CHECK(worst <= b, "sincos N=%d err %.3e > bound %.3e", n, worst, b);
        /* While the residual-angle term dominates, each step halves the error. */
        if (n >= 6 && n <= 20) {
            double r = prev / worst;
            CHECK(r > 1.6 && r < 2.5, "sincos N=%d halving ratio %.2f", n, r);
        }
        prev = worst;
    }
    /* Exact quadrant points at full precision. */
    int32_t c, s;
    cordic_sincos(0x40000000u, 30, &c, &s);
    CHECK(abs(c) <= 40 && abs(s - (1 << 30)) <= 40, "90 deg: %d %d", c, s);
    cordic_sincos(0x80000000u, 30, &c, &s);
    CHECK(abs(c + (1 << 30)) <= 40 && abs(s) <= 40, "180 deg: %d %d", c, s);
}

static void check_vector(int32_t x, int32_t y, int n, double *worst_a, double *worst_m)
{
    uint32_t a, m;
    cordic_atan2_mag(x, y, n, &a, &m);
    double ref_a = atan2((double)y, (double)x);
    double ref_m = hypot((double)x, (double)y);
    double ea = fabs(ang_err(bam_to_rad(a), ref_a));
    double em = fabs(m - ref_m) - 0.5;          /* output is rounded to integer units */
    em = em > 0 ? em / ref_m : 0.0;
    if (ea > *worst_a) *worst_a = ea;
    if (em > *worst_m) *worst_m = em;
}

static void test_vectoring(void)
{
    const int32_t edge[][2] = {
        {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {-1, -1}, {1, -1}, {-5, 0}, {-5, 1},
        {-5, -1}, {INT32_MAX, INT32_MAX}, {INT32_MIN, INT32_MIN}, {INT32_MIN, 0},
        {0, INT32_MIN}, {INT32_MAX, INT32_MIN}, {INT32_MIN, 1}, {3, 4},
        {1 << 28, 1 << 28}, {(1 << 29) - 1, 1 << 29}, {-4095, 4095}};

    double prev = 0.0;
    for (int n = 4; n <= CORDIC_MAX_ITER; n++) {
        double wa = 0.0, wm = 0.0;
        for (size_t k = 0; k < sizeof edge / sizeof edge[0]; k++)
            check_vector(edge[k][0], edge[k][1], n, &wa, &wm);
        for (int k = 0; k < 200000; k++) {
            /* Random scale from 1 to 2^31 so both normalisation paths run. */
            int bits = 1 + (int)(rng_u64() % 31);
            int64_t lim = (bits == 31) ? INT32_MAX : ((int64_t)1 << bits);
            int32_t x = (int32_t)((int64_t)(rng_u64() % (uint64_t)(2 * lim + 1)) - lim);
            int32_t y = (int32_t)((int64_t)(rng_u64() % (uint64_t)(2 * lim + 1)) - lim);
            if (x == 0 && y == 0) continue;
            /* Small integer inputs have their own quantisation; skip |v| < 2^10. */
            if (fabs((double)x) + fabs((double)y) < 1024 && n > 10) continue;
            check_vector(x, y, n, &wa, &wm);
        }
        double ba = bound_vector_angle(n), bm = bound_vector_mag_rel(n);
        CHECK(wa <= ba, "atan2 N=%d err %.3e > bound %.3e", n, wa, ba);
        CHECK(wm <= bm, "mag N=%d rel err %.3e > bound %.3e", n, wm, bm);
        if (n >= 6 && n <= 20) {
            double r = prev / wa;
            CHECK(r > 1.6 && r < 2.5, "atan2 N=%d halving ratio %.2f", n, r);
        }
        prev = wa;
    }

    uint32_t a, m;
    cordic_atan2_mag(0, 0, 16, &a, &m);
    CHECK(a == 0 && m == 0, "zero vector -> %u %u", a, m);
    cordic_atan2_mag(3, 4, 30, &a, &m);
    CHECK(m == 5, "3-4-5 magnitude %u", m);
    cordic_atan2_mag(INT32_MIN, INT32_MIN, 30, &a, &m);
    CHECK(fabs(m - 3037000499.97605) <= 3037000499.97605 * bound_vector_mag_rel(30) + 0.5,
          "INT32_MIN diagonal magnitude %u", m);
    CHECK(fabs(ang_err(bam_to_rad(a), -0.75 * M_PI)) < 1e-8, "INT32_MIN diagonal angle");
}

static uint16_t cell(double g, double nx, double ny, double nz, double sx, double sy, double sz)
{
    double i = g * (nx * sx + ny * sy + nz * sz);
    return (uint16_t)lround(i < 0 ? 0 : i);
}

static void test_sun_sensor(void)
{
    const double t0 = 30.0 * M_PI / 180.0;
    sun_sensor_cfg cfg = {(uint32_t)llround(2.0 / tan(t0) * 0x1.0p28), 24};
    CHECK(cfg.c_q28 == 929887697u, "c_q28 for 30 deg = 2*sqrt(3)*2^28 -> %u", cfg.c_q28);

    /* Boresight: all four cells read cos(30) * 4000 = 3464.1 -> 3464. */
    uint16_t c0[4] = {3464, 3464, 3464, 3464};
    uint32_t az, off;
    sun_sensor_solve(&cfg, c0, &az, &off);
    /* D = 0 takes the zero-vector path (az = 0); the off-axis pass sees y = 0 and
       leaves only the residual of the final micro-rotation. */
    CHECK(az == 0 && fabs(bam_to_rad(off)) <= bound_vector_angle(24), "boresight off=%u az=%u", off, az);

    /*
     * Hand case: sun 40 deg off axis at azimuth 0, G = 4000 counts.
     * I+x = 4000 cos(30-40) = 3939.23, I-x = 4000 cos(70) = 1368.08,
     * I+-y = 4000 cos30 cos40 = 2653.68.  Rounded: 3939, 1368, 2654, 2654.
     * Dx = 2571, S = 10615, off = atan(2 sqrt3 * 2571 / 10615) = 39.9994 deg.
     */
    uint16_t c1[4] = {3939, 1368, 2654, 2654};
    sun_sensor_solve(&cfg, c1, &az, &off);
    double off_deg = bam_to_rad(off) * 180.0 / M_PI;
    double hand = atan(2.0 * sqrt(3.0) * 2571.0 / 10615.0) * 180.0 / M_PI;
    CHECK(fabs(off_deg - hand) < 1e-5, "40 deg hand case: %.6f vs %.6f", off_deg, hand);
    CHECK(fabs(off_deg - 40.0) < 0.01, "40 deg truth: %.6f", off_deg);
    CHECK(fabs(bam_to_rad(az)) <= bound_vector_angle(24), "40 deg azimuth %u", az);

    /* Gain invariance and sweep: noiseless, finely quantised (G = 60000). */
    double worst = 0.0;
    for (int k = 0; k < 20000; k++) {
        double off_t = 55.0 * M_PI / 180.0 * sqrt(rng_unit());
        double az_t = 2.0 * M_PI * rng_unit();
        double sx = sin(off_t) * cos(az_t), sy = sin(off_t) * sin(az_t), sz = cos(off_t);
        double g = 60000.0;
        uint16_t c[4] = {cell(g, sin(t0), 0, cos(t0), sx, sy, sz), cell(g, -sin(t0), 0, cos(t0), sx, sy, sz),
                         cell(g, 0, sin(t0), cos(t0), sx, sy, sz), cell(g, 0, -sin(t0), cos(t0), sx, sy, sz)};
        /* Firmware path versus double-precision math on the same counts. */
        sun_sensor_solve(&cfg, c, &az, &off);
        double dx = c[0] - (double)c[1], dy = c[2] - (double)c[3], s = c[0] + c[1] + (double)c[2] + c[3];
        double off_ref = atan2(2.0 / tan(t0) * hypot(dx, dy), s);
        double e = fabs(ang_err(bam_to_rad(off), off_ref));
        if (e > worst) worst = e;
    }
    /* Two chained N=24 passes; Q12 counts keep the c|D| rounding below 1e-7 rad. */
    CHECK(worst < 2.0 * bound_vector_angle(24) + 1e-7, "solver vs double: %.3e rad", worst);
}

int main(void)
{
    test_tables();
    test_rotation();
    test_vectoring();
    test_sun_sensor();
    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("all CORDIC and sun-sensor checks passed\n");
    return 0;
}
