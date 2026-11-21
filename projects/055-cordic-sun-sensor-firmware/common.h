/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * Host-side helpers shared by the test and the study: a reproducible PRNG,
 * BAM conversions and the analytic CORDIC error bounds.
 */
#ifndef COMMON_H
#define COMMON_H

#include <math.h>
#include <stdint.h>
#include "cordic.h"

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;

static inline uint64_t rng_u64(void)
{
    /* xorshift64* */
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545F4914F6CDD1Dull;
}

static inline double rng_unit(void) { return (double)(rng_u64() >> 11) * 0x1.0p-53; }

static inline double rng_gauss(void)
{
    double u1 = rng_unit(), u2 = rng_unit();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static inline double bam_to_rad(uint32_t a) { return (int32_t)a * (2.0 * M_PI / BAM_PER_TURN); }

static inline uint32_t rad_to_bam(double r)
{
    double t = r / (2.0 * M_PI);
    t -= floor(t);
    return (uint32_t)(uint64_t)llround(t * BAM_PER_TURN);
}

/* Wrapped angle difference in radians. */
static inline double ang_err(double a, double b)
{
    double d = fmod(a - b, 2.0 * M_PI);
    if (d > M_PI) d -= 2.0 * M_PI;
    if (d < -M_PI) d += 2.0 * M_PI;
    return d;
}

static inline double cordic_gain(int n)
{
    double k = 1.0;
    for (int i = 0; i < n; i++) k *= sqrt(1.0 + ldexp(1.0, -2 * i));
    return k;
}

/*
 * Rotation mode, worst-case vector error on the unit circle:
 *   atan(2^-(N-1))           residual angle after N steps (convergence bound)
 * + N/2 BAM                  atan table rounding
 * + (sqrt2*(N-1) + 1/2)*K_N  truncation of the N-1 nonzero shifts, each <1 LSB
 *                            per component and amplified at most by K_N, plus
 *                            the Q30 rounding of 1/K_N, in units of 2^-30.
 */
static inline double bound_rotation(int n)
{
    double k = cordic_gain(n);
    return atan(ldexp(1.0, -(n - 1))) + 0.5 * n * (2.0 * M_PI / BAM_PER_TURN)
         + (M_SQRT2 * (n - 1) + 0.5) * k * ldexp(1.0, -30);
}

/*
 * Vectoring mode, worst-case angle error. After normalisation r >= 2^28, so a
 * datapath error vector E (<= sqrt2*(N-1)*K_N LSB plus up to sqrt2*8 LSB from
 * a right-shift normalisation) tilts the final vector by at most E/(K_N r).
 * The final computed vector may itself sit up to atan(2^-(N-1)) + E/(K_N r)
 * off the axis, hence the factor 2.
 */
static inline double vec_datapath_rel(int n)
{
    double k = cordic_gain(n);
    return (M_SQRT2 * (n - 1) * k + M_SQRT2 * 8.0) / (k * ldexp(1.0, 28));
}

static inline double bound_vector_angle(int n)
{
    return atan(ldexp(1.0, -(n - 1))) + 0.5 * n * (2.0 * M_PI / BAM_PER_TURN)
         + 2.0 * vec_datapath_rel(n);
}

/* Relative magnitude error: cos of the residual angle, datapath, 1/K rounding. */
static inline double bound_vector_mag_rel(int n)
{
    double phi = atan(ldexp(1.0, -(n - 1))) + vec_datapath_rel(n);
    return (1.0 - cos(phi)) + vec_datapath_rel(n) + ldexp(1.0, -31) * cordic_gain(n);
}

#endif
