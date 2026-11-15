/*
 * Astron Built this AP
 * Project 039 | Rocket propulsion | BATES solid motor internal ballistics
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-15 IST
 */
#ifndef MOTOR_DESIGN_H
#define MOTOR_DESIGN_H

#include "ballistics.h"

#define PA_SEA_LEVEL 101325.0
#define P_DESIGN_MAX 5.0e6     /* MEOP target for the case, Pa          */
#define P_IGNITER    0.3e6     /* chamber pressure the igniter leaves, Pa */

/* Generic AP/HTPB composite, textbook-level values: r = 7 mm/s at 6.9 MPa
 * with n = 0.35, products gamma 1.21, M = 24 g/mol, Tc = 3000 K. */
static inline Propellant design_propellant(void)
{
    Propellant p;
    p.rho   = 1750.0;
    p.n     = 0.35;
    p.a     = 0.007 / pow(6.9e6, p.n);
    p.gamma = 1.21;
    p.R     = 8314.462618 / 24.0;
    p.Tc    = 3000.0;
    p.eta_c = 0.95;
    return p;
}

/* 98 mm case, four segments of neutral length, 2.5 mm spacer per joint. */
static inline Grain design_grain(void)
{
    Grain g;
    g.N = 4;
    g.D = 0.098;
    g.d = 0.035;
    g.L = grain_neutral_length(g.D, g.d);
    g.Lcase = g.N * g.L + 0.010;
    return g;
}

/* Expansion ratio that makes pe == pa at chamber pressure pc. */
static inline double eps_for_optimum(double pc, double pa, double gamma)
{
    double M = sqrt(2.0 / (gamma - 1.0)
                    * (pow(pc / pa, (gamma - 1.0) / gamma) - 1.0));
    return area_ratio(M, gamma);
}

#endif
