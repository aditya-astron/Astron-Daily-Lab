/*
 * Astron Built this AP
 * Project 039 | Rocket propulsion | BATES solid motor internal ballistics
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-15 IST
 */
#ifndef BALLISTICS_H
#define BALLISTICS_H

#include <stddef.h>

/* Propellant and combustion-gas properties (SI). */
typedef struct {
    double rho;    /* propellant density, kg/m^3                      */
    double a;      /* burn-rate coefficient, r = a * Pc^n, m/s/Pa^n  */
    double n;      /* pressure exponent                               */
    double gamma;  /* ratio of specific heats of the products         */
    double R;      /* specific gas constant of the products, J/kg/K   */
    double Tc;     /* chamber temperature, K                          */
    double eta_c;  /* c* efficiency (combustion losses)               */
} Propellant;

/* N identical cylindrical BATES segments in a case of inner diameter D.
 * Outer surface bonded/inhibited; core and both end faces burn. */
typedef struct {
    int    N;
    double D;      /* grain outer diameter, m       */
    double d;      /* initial core diameter, m      */
    double L;      /* initial segment length, m     */
    double Lcase;  /* free length of the case, m    */
} Grain;

typedef struct {
    double At;     /* throat area, m^2              */
    double eps;    /* expansion ratio Ae/At         */
} Nozzle;

typedef struct {
    double t_burnout;   /* web burnout time, s                       */
    double t_end;       /* end of simulation (Pc < cutoff), s        */
    double P_max;       /* peak chamber pressure, Pa                 */
    double F_max;       /* peak thrust, N                            */
    double I_total;     /* total impulse at ambient pa, N s          */
    double I_vac;       /* total impulse in vacuum, N s              */
    double m_prop;      /* initial propellant mass, kg               */
    double m_nozzle;    /* mass expelled through the throat, kg      */
    double m_gas_end;   /* gas left in the chamber at t_end, kg      */
    double P_burnout;   /* chamber pressure at web burnout, Pa       */
    double V_free_burnout; /* chamber free volume at burnout, m^3   */
    size_t n_samples;   /* trace length                              */
} MotorResult;

/* One trace sample, filled when the caller passes a buffer. */
typedef struct {
    double t, x, Pc, F, Ab;
} TraceSample;

double cstar_ideal(const Propellant *p);
double cstar(const Propellant *p);

/* Isentropic area ratio A/A* at Mach M. */
double area_ratio(double M, double gamma);
/* Supersonic Mach number for area ratio eps (eps >= 1). */
double mach_supersonic(double eps, double gamma);
/* Exit-to-chamber pressure ratio for a nozzle of expansion ratio eps. */
double pe_over_pc(double eps, double gamma);
/* Thrust coefficient at chamber pressure pc and ambient pa. */
double thrust_coeff(double eps, double gamma, double pc, double pa);

double grain_web(const Grain *g);
double grain_burn_area(const Grain *g, double x);   /* x = web burned */
double grain_prop_volume(const Grain *g, double x); /* remaining      */
double grain_free_volume(const Grain *g, double x); /* chamber gas    */
/* Segment length that makes Ab(0) == Ab(web): L = (3D + d) / 2. */
double grain_neutral_length(double D, double d);

/* Equilibrium (quasi-steady, rho_gas << rho_p) chamber pressure for Kn. */
double pc_equilibrium(const Propellant *p, double Kn);
/* Same with the rho_gas correction: Pc = (a c* Kn (rho_p - rho_g))^(1/(1-n)),
 * solved by fixed-point iteration. */
double pc_equilibrium_exact(const Propellant *p, double Kn);

/* Throat area that puts the peak equilibrium pressure at Pmax. */
double size_throat(const Propellant *p, const Grain *g, double Pmax);

/* Time-march the coupled chamber mass balance and regression through
 * ignition fill, burn and tail-off. trace may be NULL. */
MotorResult simulate_motor(const Propellant *p, const Grain *g,
                           const Nozzle *nz, double pa, double dt,
                           double P_ign, TraceSample *trace,
                           size_t trace_cap, size_t trace_every);

#endif
