/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * Micro-rotation i:  x' = x - d*(y >> i),  y' = y + d*(x >> i),  z' = z - d*atan(2^-i)
 * Each step stretches the vector by k_i = sqrt(1 + 2^-2i); after N steps the
 * gain is K_N = prod k_i (-> 1.6467602581), removed with one Q31 multiply.
 */
#include "cordic.h"

/* round(atan(2^-i) / 2pi * 2^32); regenerated and checked in test_cordic.c */
const int32_t cordic_atan_bam[CORDIC_MAX_ITER] = {
    536870912, 316933406, 167458907, 85004756, 42667331, 21354465, 10679838,
    5340245, 2670163, 1335087, 667544, 333772, 166886, 83443, 41722, 20861,
    10430, 5215, 2608, 1304, 652, 326, 163, 81, 41, 20, 10, 5, 3, 1};

/* round(2^31 / K_N) for N = 1..30 */
const int32_t cordic_inv_gain_q31[CORDIC_MAX_ITER] = {
    1518500250, 1358187913, 1317635818, 1307460871, 1304914694, 1304277995,
    1304118810, 1304079014, 1304069065, 1304066577, 1304065955, 1304065800,
    1304065761, 1304065751, 1304065749, 1304065748, 1304065748, 1304065748,
    1304065748, 1304065748, 1304065748, 1304065748, 1304065748, 1304065748,
    1304065748, 1304065748, 1304065748, 1304065748, 1304065748, 1304065748};

static int clamp_iter(int n)
{
    if (n < 1) return 1;
    if (n > CORDIC_MAX_ITER) return CORDIC_MAX_ITER;
    return n;
}

void cordic_sincos(uint32_t angle, int n_iter, int32_t *cos_q30, int32_t *sin_q30)
{
    int n = clamp_iter(n_iter);
    int32_t z = (int32_t)angle;
    int flip = 0;

    /* Convergence range is about +-99.9 deg, so fold |angle| > 90 deg by pi. */
    if (z > (1 << 30) || z < -(1 << 30)) {
        z = (int32_t)(angle + 0x80000000u);
        flip = 1;
    }

    /* Pre-scale by 1/K_N so the output lands on the unit circle. */
    int32_t x = (cordic_inv_gain_q31[n - 1] + 1) >> 1;
    int32_t y = 0;

    for (int i = 0; i < n; i++) {
        int32_t xs = x >> i, ys = y >> i;
        if (z >= 0) {
            x -= ys; y += xs; z -= cordic_atan_bam[i];
        } else {
            x += ys; y -= xs; z += cordic_atan_bam[i];
        }
    }
    *cos_q30 = flip ? -x : x;
    *sin_q30 = flip ? -y : y;
}

void cordic_atan2_mag(int32_t x, int32_t y, int n_iter, uint32_t *angle, uint32_t *mag)
{
    int n = clamp_iter(n_iter);
    uint32_t ax = x < 0 ? 0u - (uint32_t)x : (uint32_t)x;
    uint32_t ay = y < 0 ? 0u - (uint32_t)y : (uint32_t)y;
    uint32_t m = ax > ay ? ax : ay;

    if (m == 0) {
        *angle = 0;
        *mag = 0;
        return;
    }

    /*
     * Block-normalise: shift so max(|x|,|y|) is in [2^28, 2^29). The headroom
     * covers the sqrt(2) diagonal times K_N < 1.647: 2^29 * 2.33 < 2^31.
     * A right shift (large inputs) truncates at most 3 LSBs per component.
     */
    int lz = __builtin_clz(m);           /* m in [2^(31-lz), 2^(32-lz)) */
    int shift = lz - 3;                  /* >0: shift left, <0: shift right */
    if (shift >= 0) {
        x = (int32_t)((uint32_t)x << shift);
        y = (int32_t)((uint32_t)y << shift);
    } else {
        x >>= -shift;
        y >>= -shift;
    }

    /* Fold the left half-plane onto the right with a pi rotation. */
    uint32_t z = 0;
    if (x < 0) {
        x = -x;
        y = -y;
        z = 0x80000000u;
    }

    /* Drive y to zero; the accumulated rotation is the angle. */
    for (int i = 0; i < n; i++) {
        int32_t xs = x >> i, ys = y >> i;
        if (y > 0) {
            x += ys; y -= xs; z += (uint32_t)cordic_atan_bam[i];
        } else {
            x -= ys; y += xs; z -= (uint32_t)cordic_atan_bam[i];
        }
    }
    *angle = z;

    /* x_N = K_N * r; remove the gain (32x32->64, one SMULL on Cortex-M) and undo the normalisation. */
    uint64_t r = (uint64_t)(uint32_t)x * (uint32_t)cordic_inv_gain_q31[n - 1];
    int out_shift = 31 + shift;          /* >= 28 since shift >= -3 */
    *mag = (uint32_t)((r + (1ull << (out_shift - 1))) >> out_shift);
}
