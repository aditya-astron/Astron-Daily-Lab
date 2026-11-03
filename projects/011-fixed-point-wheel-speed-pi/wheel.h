/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
/* Reaction wheel speed dynamics with an ideal inner current loop:
 *   J dw/dt = Kt*i - b*w - tau_load
 * Between controller samples i and tau_load are held constant, so the ODE
 * is linear with constant input and is advanced with its exact solution. */
#ifndef WHEEL_H
#define WHEEL_H

typedef struct {
    double J;   /* rotor inertia, kg m^2 */
    double kt;  /* torque constant, N m / A */
    double b;   /* viscous drag, N m s / rad */
    double w;   /* speed, rad/s */
} wheel_t;

void wheel_step(wheel_t *wh, double i, double tau_load, double dt);

#endif
